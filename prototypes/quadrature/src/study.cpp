// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Compares quadrature rule generators on cut cells against exact per-cell values.
//
// Test problem: sphere |x - centre| = radius in [-1, 1]^3, meshed with n^3 hexahedra
// or 6 n^3 Kuhn tetrahedra. Parts are selected with the same expressions as
// HOCutResult, e.g. "phi < 0", "phi > 0", "phi = 0".
//
// Usage:
//   quadrature_study [--mesh hex|tet] [--n 16,32] [--q 3,5] [--centre x,y,z]
//                    [--radius r] [--gen algoim-auto,alpha-split,...]
//                    [--part "phi < 0"]... [--csv file] [--plane] [--vtk prefix]
//
// --plane uses phi = x + 0.3 y - 0.2 z instead of the sphere. Every rule must then
// be exact: the volume below is 4 and the area 4 sqrt(1.13). Only global errors
// are reported for it.
//
// Generators: algoim-auto, algoim-gl, gl-cellmask, alpha, split, alpha-split, quadgen (hex only),
// certify (certify-and-bisect, margin 0.25) or certify:<margin>.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "certified.h"
#include "clipped_box.h"
#include "exact_reference.h"
#include "generators.h"
#include "selection_expr.h"
#include "vtk_output.h"

using namespace cutcells;
using namespace cutcells::proto;

namespace
{
struct StudyConfig
{
    std::string mesh = "tet";
    std::vector<int> ns = {16, 32};
    std::vector<int> qs = {3};
    Vec3 centre = {0.0123, -0.0371, 0.0217};
    double radius = 0.7;
    std::vector<std::string> generators = {"algoim-auto", "algoim-gl", "gl-cellmask", "alpha", "alpha-split"};
    std::vector<std::string> parts;
    std::string csv;
    bool plane = false; ///< planar level set x + 0.3 y - 0.2 z instead of the sphere
    std::string vtk;    ///< prefix for point-cloud .vtu files, one per run
};

struct Metrics
{
    long cut_cells = 0;
    long points = 0;
    long splits = 0;
    long uncertified = 0;
    double sum_abs = 0;    ///< sum over cut cells of |generated - exact|
    double sum_gen = 0;    ///< generated part measure over cut cells
    double sum_uncut = 0;  ///< exact part measure over uncut cells
    double max_rel = 0;    ///< worst cell, relative to the cell's exact value
    double seconds = 0;
};

std::vector<std::string> split_list(const std::string& s)
{
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ','))
        if (!item.empty())
            out.push_back(item);
    return out;
}

StudyConfig parse_args(int argc, char** argv)
{
    StudyConfig cfg;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        auto next = [&]() -> std::string
        {
            if (i + 1 >= argc)
                throw std::runtime_error("missing value after " + a);
            return argv[++i];
        };
        if (a == "--mesh")
            cfg.mesh = next();
        else if (a == "--n")
        {
            cfg.ns.clear();
            for (const auto& s : split_list(next()))
                cfg.ns.push_back(std::stoi(s));
        }
        else if (a == "--q")
        {
            cfg.qs.clear();
            for (const auto& s : split_list(next()))
                cfg.qs.push_back(std::stoi(s));
        }
        else if (a == "--centre")
        {
            const auto c = split_list(next());
            if (c.size() != 3)
                throw std::runtime_error("--centre needs x,y,z");
            for (int d = 0; d < 3; ++d)
                cfg.centre[d] = std::stod(c[d]);
        }
        else if (a == "--radius")
            cfg.radius = std::stod(next());
        else if (a == "--gen")
            cfg.generators = split_list(next());
        else if (a == "--part")
            cfg.parts.push_back(next());
        else if (a == "--csv")
            cfg.csv = next();
        else if (a == "--plane")
            cfg.plane = true;
        else if (a == "--vtk")
            cfg.vtk = next();
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    if (cfg.parts.empty())
        cfg.parts = {"phi < 0", "phi = 0"};
    return cfg;
}

/// A cell of the test mesh: the clipped box handed to the generators and its
/// physical polytope (relative to the sphere centre) for the exact reference.
struct TestCell
{
    ClippedBox box;
    std::vector<exact::Face> faces;
    Vec3 centroid;
    double radius = 0; ///< circumradius about the centroid
    double volume = 0;
};

std::vector<TestCell> grid_cells(const std::string& mesh, const Vec3& lo, double h, const Vec3& centre)
{
    auto rel = [&](const Vec3& x) { return exact::V3{x[0] - centre[0], x[1] - centre[1], x[2] - centre[2]}; };
    std::vector<TestCell> cells;
    if (mesh == "hex")
    {
        TestCell c;
        c.box = hex_cell(lo, h);
        c.faces = exact::box_faces(rel(lo), h);
        c.centroid = {lo[0] + 0.5 * h, lo[1] + 0.5 * h, lo[2] + 0.5 * h};
        c.radius = 0.5 * std::sqrt(3.0) * h;
        c.volume = h * h * h;
        cells.push_back(c);
        return cells;
    }
    if (mesh != "tet")
        throw std::runtime_error("unknown mesh: " + mesh);
    std::array<int, 3> p = {0, 1, 2};
    do
    {
        // Kuhn tetrahedron: path from lo along the axes p[0], p[1], p[2]
        std::array<Vec3, 4> X;
        X[0] = lo;
        for (int k = 0; k < 3; ++k)
        {
            X[k + 1] = X[k];
            X[k + 1][p[k]] += h;
        }
        TestCell c;
        c.box = tet_cell(X);
        c.faces = exact::tet_faces({rel(X[0]), rel(X[1]), rel(X[2]), rel(X[3])});
        c.centroid = {0, 0, 0};
        for (const Vec3& x : X)
            for (int d = 0; d < 3; ++d)
                c.centroid[d] += 0.25 * x[d];
        for (const Vec3& x : X)
        {
            double r2 = 0;
            for (int d = 0; d < 3; ++d)
                r2 += (x[d] - c.centroid[d]) * (x[d] - c.centroid[d]);
            c.radius = std::max(c.radius, std::sqrt(r2));
        }
        c.volume = h * h * h / 6.0;
        cells.push_back(c);
    } while (std::next_permutation(p.begin(), p.end()));
    return cells;
}
/// Planar level set: every generator must reproduce the totals to rounding.
int plane_study(const StudyConfig& cfg)
{
    LevelSet ls;
    ls.degree = 1;
    ls.value = [](const Vec3& x) { return x[0] + 0.3 * x[1] - 0.2 * x[2]; };
    std::printf("Plane x + 0.3 y - 0.2 z: global relative errors of the totals.\n");
    for (int n : cfg.ns)
        for (int q : cfg.qs)
            for (const std::string& gen : cfg.generators)
            {
                const bool certify = gen.rfind("certify", 0) == 0;
                CertifyOptions copt;
                if (certify && gen.find(':') != std::string::npos)
                    copt.margin = std::stod(gen.substr(gen.find(':') + 1));
                if (gen == "quadgen")
                    continue;
                const GeneratorOptions opt = certify ? GeneratorOptions{} : generator_preset(gen);
                for (const std::string& part : cfg.parts)
                {
                    SelectionExpr expr = parse_selection_expr(part);
                    compile_selection_expr(expr, {"phi"});
                    const SelectionTerm& term = expr.terms.front();
                    const PartKind kind = part_kind(term);
                    const double exact = kind == PartKind::interface ? 4.0 * std::sqrt(1.13) : 4.0;
                    const double h = 2.0 / n;
                    double total = 0;
                    long cut = 0;
                    for (int i0 = 0; i0 < n; ++i0)
                        for (int i1 = 0; i1 < n; ++i1)
                            for (int i2 = 0; i2 < n; ++i2)
                            {
                                const Vec3 lo = {-1 + h * i0, -1 + h * i1, -1 + h * i2};
                                for (const TestCell& cell : grid_cells(cfg.mesh, lo, h, {0, 0, 0}))
                                {
                                    // vertex signs decide cut / inside / outside for a plane
                                    bool pos = false, neg = false;
                                    for (int v = 0; v < 8; ++v)
                                    {
                                        const Vec3 u = {double(v & 1), double((v >> 1) & 1), double((v >> 2) & 1)};
                                        if (!inside_clips(cell.box, u) && !cell.box.clips.empty())
                                            continue;
                                        const double val = ls.value(physical_point(cell.box, u));
                                        pos |= val > 0;
                                        neg |= val < 0;
                                    }
                                    if (!(pos && neg))
                                    {
                                        if (kind == PartKind::negative && neg)
                                            total += cell.volume;
                                        if (kind == PartKind::positive && pos)
                                            total += cell.volume;
                                        continue;
                                    }
                                    ++cut;
                                    Rule rule;
                                    GeneratorStats stats;
                                    if (certify)
                                    {
                                        CertifyStats cs;
                                        certified_bisection(cell.box, ls, term, q, copt, rule, cs);
                                    }
                                    else
                                        algoim_clipped_box(cell.box, ls, term, q, opt, rule, stats);
                                    for (double w : rule.weights)
                                        total += w;
                                }
                            }
                    const double exact_part = kind == PartKind::positive ? 4.0 : exact;
                    std::printf("%-4s %4d %2d %-12s %-8s | cut %7ld | global %.1e\n", cfg.mesh.c_str(), n, q, gen.c_str(),
                                part.c_str(), cut, std::abs(total - exact_part) / exact_part);
                }
            }
    return 0;
}
} // namespace

int main(int argc, char** argv)
{
    StudyConfig cfg;
    try
    {
        cfg = parse_args(argc, argv);
    }
    catch (const std::exception& e)
    {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }
    if (cfg.plane)
        return plane_study(cfg);
    const double r = cfg.radius;
    const Vec3 c = cfg.centre;
    const double ball = 4.0 / 3.0 * M_PI * r * r * r, sphere = 4.0 * M_PI * r * r;

    LevelSet ls;
    ls.degree = 2;
    ls.value = [&](const Vec3& x)
    { return (x[0] - c[0]) * (x[0] - c[0]) + (x[1] - c[1]) * (x[1] - c[1]) + (x[2] - c[2]) * (x[2] - c[2]) - r * r; };

    std::ofstream csv;
    if (!cfg.csv.empty())
    {
        csv.open(cfg.csv);
        csv << "mesh,n,q,generator,part,cut_cells,points_per_cut_cell,global_rel,l1_rel,max_cell_rel,splits,uncertified,us_per_cut_cell\n";
    }
    std::printf("Sphere centre (%g, %g, %g), radius %g. Errors relative to the ball volume (volume parts)\n"
                "or the sphere area (interface); per-cell L1 = sum over cut cells of |cell error|.\n",
                c[0], c[1], c[2], r);
    std::printf("%-4s %4s %2s %-12s %-8s | %7s | %8s | %9s %9s %9s | %6s %6s | %8s\n", "mesh", "n", "q", "generator",
                "part", "cut", "pts/cut", "global", "L1", "worst", "splits", "uncert", "us/cut");

    for (int n : cfg.ns)
        for (int q : cfg.qs)
            for (const std::string& gen : cfg.generators)
            {
                if (gen == "quadgen" && cfg.mesh != "hex")
                {
                    std::printf("(quadgen skipped: hexahedra only)\n");
                    continue;
                }
                const bool certify = gen.rfind("certify", 0) == 0;
                CertifyOptions copt;
                if (certify && gen.find(':') != std::string::npos)
                    copt.margin = std::stod(gen.substr(gen.find(':') + 1));
                const GeneratorOptions opt = gen == "quadgen" || certify ? GeneratorOptions{} : generator_preset(gen);
                for (const std::string& part : cfg.parts)
                {
                    SelectionExpr expr = parse_selection_expr(part);
                    compile_selection_expr(expr, {"phi"});
                    if (expr.terms.size() != 1)
                        throw std::runtime_error("prototype: one selection term per part");
                    const SelectionTerm& term = expr.terms.front();
                    const PartKind kind = part_kind(term);
                    const bool surface = kind == PartKind::interface;
                    const double scale = surface ? sphere : ball;
                    const double total = kind == PartKind::negative   ? ball
                                         : kind == PartKind::positive ? 8.0 - ball
                                         : kind == PartKind::interface ? sphere
                                                                       : 8.0;

                    const double h = 2.0 / n;
                    Metrics m;
                    std::vector<double> vtk_points, vtk_weights;
                    std::vector<std::int32_t> vtk_cells;
                    std::int32_t cell_index = -1;
                    for (int i0 = 0; i0 < n; ++i0)
                        for (int i1 = 0; i1 < n; ++i1)
                            for (int i2 = 0; i2 < n; ++i2)
                            {
                                const Vec3 lo = {-1 + h * i0, -1 + h * i1, -1 + h * i2};
                                for (const TestCell& cell : grid_cells(cfg.mesh, lo, h, c))
                                {
                                    ++cell_index;
                                    double dc = 0;
                                    for (int d = 0; d < 3; ++d)
                                        dc += (cell.centroid[d] - c[d]) * (cell.centroid[d] - c[d]);
                                    dc = std::sqrt(dc);
                                    double area = 0, inside = 0;
                                    if (dc - cell.radius >= r)
                                        inside = 0; // outside the ball
                                    else if (dc + cell.radius <= r)
                                        inside = cell.volume;
                                    else
                                    {
                                        area = exact::sphere_area(cell.faces, r);
                                        inside = exact::ball_volume(cell.faces, r, area);
                                    }
                                    const double exact_part = kind == PartKind::negative   ? inside
                                                              : kind == PartKind::positive ? cell.volume - inside
                                                              : kind == PartKind::interface ? area
                                                                                            : cell.volume;
                                    if (area <= 0.0)
                                    {
                                        m.sum_uncut += exact_part;
                                        continue;
                                    }
                                    Rule rule;
                                    GeneratorStats stats;
                                    const auto t0 = std::chrono::steady_clock::now();
                                    if (gen == "quadgen")
                                        algoim_quadgen_sphere(cell.box, c, r, term, q, rule);
                                    else if (certify)
                                    {
                                        CertifyStats cs;
                                        certified_bisection(cell.box, ls, term, q, copt, rule, cs);
                                        stats.splits += cs.bisections;
                                        m.uncertified += cs.uncertified;
                                    }
                                    else
                                        algoim_clipped_box(cell.box, ls, term, q, opt, rule, stats);
                                    m.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                                    double value = 0;
                                    for (double w : rule.weights)
                                        value += w;
                                    if (!cfg.vtk.empty())
                                        for (int i = 0; i < rule.n_points(); ++i)
                                        {
                                            // reference coordinates of these test cells are box coordinates
                                            const Vec3 x = physical_point(
                                                cell.box, {rule.points[3 * i], rule.points[3 * i + 1], rule.points[3 * i + 2]});
                                            vtk_points.insert(vtk_points.end(), x.begin(), x.end());
                                            vtk_weights.push_back(rule.weights[i]);
                                            vtk_cells.push_back(cell_index);
                                        }
                                    ++m.cut_cells;
                                    m.points += rule.n_points();
                                    m.splits += stats.splits;
                                    m.sum_gen += value;
                                    m.sum_abs += std::abs(value - exact_part);
                                    const double floor = surface ? 1e-3 * h * h : 1e-3 * h * h * h;
                                    if (exact_part > floor)
                                        m.max_rel = std::max(m.max_rel, std::abs(value - exact_part) / exact_part);
                                }
                            }
                    const double global = std::abs(m.sum_gen + m.sum_uncut - total) / scale;
                    const double l1 = m.sum_abs / scale;
                    const double pts = m.cut_cells ? static_cast<double>(m.points) / m.cut_cells : 0.0;
                    const double us = m.cut_cells ? m.seconds / m.cut_cells * 1e6 : 0.0;
                    std::printf("%-4s %4d %2d %-12s %-8s | %7ld | %8.1f | %9.1e %9.1e %9.1e | %6ld %6ld | %8.1f\n",
                                cfg.mesh.c_str(), n, q, gen.c_str(), part.c_str(), m.cut_cells, pts, global, l1, m.max_rel,
                                m.splits, m.uncertified, us);
                    std::fflush(stdout);
                    if (!cfg.vtk.empty())
                    {
                        const char* slug = kind == PartKind::negative   ? "negative"
                                           : kind == PartKind::positive ? "positive"
                                           : kind == PartKind::interface ? "interface"
                                                                         : "whole";
                        std::string name = gen;
                        std::replace(name.begin(), name.end(), ':', '_');
                        write_point_cloud(cfg.vtk + "_" + cfg.mesh + "_" + name + "_" + slug + "_n" + std::to_string(n) + "_q"
                                              + std::to_string(q) + ".vtu",
                                          vtk_points, vtk_weights, vtk_cells);
                    }
                    if (csv)
                        csv << cfg.mesh << ',' << n << ',' << q << ',' << gen << ",\"" << part << "\"," << m.cut_cells << ','
                            << pts << ',' << global << ',' << l1 << ',' << m.max_rel << ',' << m.splits << ','
                            << m.uncertified << ',' << us << '\n';
                }
            }
    return 0;
}
