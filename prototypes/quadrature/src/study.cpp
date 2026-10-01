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
//                    [--leaves prefix] [--leaf-degree p]
//                    [--diagnose] [--only cell] [--masks M] [--no-diagonal]
//
// --diagnose (certify) prints why bisections happen and the worst cell; --only
// restricts the run to one cell index; --masks and --no-diagonal set
// CertifyOptions::mask_subdivisions and diagonal_frames.
//
// --leaves writes the certify engine's leaf cells (Lagrange cells of degree p, 3 by
// default) for every cut cell, plus the uncut cells of volume parts as linear cells.
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
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "certified.h"
#include "clipped_box.h"
#include "exact_reference.h"
#include "generators.h"
#include "leaf_mesh.h"
#include "selection_expr.h"
#include "test_mesh.h"
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
    std::string leaves; ///< prefix for leaf-cell .vtu files (certify generators)
    int leaf_degree = 3;
    bool diagnose = false; ///< certify: print why bisections happen
    int only = -1;         ///< restrict the study to this cell index
    int masks = 1;         ///< certify: CertifyOptions::mask_subdivisions
    bool diagonal = true;  ///< certify: CertifyOptions::diagonal_frames
};

struct Metrics
{
    long cut_cells = 0;
    long points = 0;
    long splits = 0;
    long uncertified = 0;
    long rotations = 0;
    double sum_abs = 0;    ///< sum over cut cells of |generated - exact|
    double sum_gen = 0;    ///< generated part measure over cut cells
    double sum_uncut = 0;  ///< exact part measure over uncut cells
    double max_rel = 0;    ///< worst cell, relative to the cell's exact value
    double seconds = 0;
    // the worst cell, for --diagnose
    std::int32_t worst_cell = -1;
    double worst_exact = 0, worst_value = 0;
    long worst_points = 0, worst_splits = 0, worst_uncertified = 0;
    std::vector<std::int32_t> uncertified_cells; ///< cells with uncertified boxes
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
        else if (a == "--leaves")
            cfg.leaves = next();
        else if (a == "--leaf-degree")
            cfg.leaf_degree = std::stoi(next());
        else if (a == "--diagnose")
            cfg.diagnose = true;
        else if (a == "--only")
            cfg.only = std::stoi(next());
        else if (a == "--masks")
            cfg.masks = std::stoi(next());
        else if (a == "--no-diagonal")
            cfg.diagonal = false;
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    if (cfg.parts.empty())
        cfg.parts = {"phi < 0", "phi = 0"};
    return cfg;
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
                copt.mask_subdivisions = cfg.masks;
                copt.diagonal_frames = cfg.diagonal;
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
                copt.diagnose = cfg.diagnose;
                copt.mask_subdivisions = cfg.masks;
                copt.diagonal_frames = cfg.diagonal;
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
                    std::map<std::string, long> causes;
                    std::vector<double> vtk_points, vtk_weights;
                    std::vector<std::int32_t> vtk_cells;
                    std::int32_t cell_index = -1;
                    const bool want_leaves = !cfg.leaves.empty() && certify;
                    LeafMesh leaf_mesh, cut_cells; // leaves, and the cut background cells for reference
                    long incomplete = 0;
                    auto add_linear_cell = [&](const TestCell& cell, LeafMesh& leaf_mesh)
                    {
                        // an uncut cell of a volume part, as a plain VTK cell
                        const bool hex = cell.box.clips.empty();
                        const std::vector<Vec3> corners
                            = hex ? std::vector<Vec3>{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}
                                  : std::vector<Vec3>{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
                        const std::int32_t first = leaf_mesh.n_points();
                        for (const Vec3& u : corners)
                        {
                            const Vec3 x = physical_point(cell.box, u);
                            leaf_mesh.connectivity.push_back(leaf_mesh.n_points());
                            leaf_mesh.points.insert(leaf_mesh.points.end(), x.begin(), x.end());
                        }
                        if (!hex && jacobian_determinant(cell.box) < 0)
                            std::swap(leaf_mesh.connectivity[first + 1], leaf_mesh.connectivity[first + 2]);
                        leaf_mesh.offsets.push_back(static_cast<std::int32_t>(leaf_mesh.connectivity.size()));
                        leaf_mesh.types.push_back(hex ? vtk_hexahedron : vtk_tetra);
                        leaf_mesh.parent.push_back(cell_index);
                        leaf_mesh.degree.push_back(1);
                    };
                    for (int i0 = 0; i0 < n; ++i0)
                        for (int i1 = 0; i1 < n; ++i1)
                            for (int i2 = 0; i2 < n; ++i2)
                            {
                                const Vec3 lo = {-1 + h * i0, -1 + h * i1, -1 + h * i2};
                                for (const TestCell& cell : grid_cells(cfg.mesh, lo, h, c))
                                {
                                    ++cell_index;
                                    if (cfg.only >= 0 && cell_index != cfg.only)
                                        continue;
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
                                        if (want_leaves && !surface && exact_part > 0.5 * cell.volume)
                                            add_linear_cell(cell, leaf_mesh);
                                        continue;
                                    }
                                    if (want_leaves)
                                    {
                                        CertifyStats ls_stats;
                                        certified_leaves(cell.box, ls, term, cfg.leaf_degree, copt, cell_index, leaf_mesh, ls_stats);
                                        incomplete += ls_stats.incomplete_leaves;
                                        add_linear_cell(cell, cut_cells);
                                    }
                                    Rule rule;
                                    GeneratorStats stats;
                                    const long uncertified_before = m.uncertified;
                                    const auto t0 = std::chrono::steady_clock::now();
                                    if (gen == "quadgen")
                                        algoim_quadgen_sphere(cell.box, c, r, term, q, rule);
                                    else if (certify)
                                    {
                                        CertifyStats cs;
                                        certified_bisection(cell.box, ls, term, q, copt, rule, cs);
                                        stats.splits += cs.bisections;
                                        m.uncertified += cs.uncertified;
                                        m.rotations += cs.rotations;
                                        for (const auto& [key, count] : cs.causes)
                                            causes[key] += count;
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
                                    if (m.uncertified > uncertified_before)
                                        m.uncertified_cells.push_back(cell_index);
                                    ++m.cut_cells;
                                    m.points += rule.n_points();
                                    m.splits += stats.splits;
                                    m.sum_gen += value;
                                    m.sum_abs += std::abs(value - exact_part);
                                    const double floor = surface ? 1e-3 * h * h : 1e-3 * h * h * h;
                                    if (exact_part > floor && std::abs(value - exact_part) / exact_part > m.max_rel)
                                    {
                                        m.max_rel = std::abs(value - exact_part) / exact_part;
                                        m.worst_cell = cell_index;
                                        m.worst_exact = exact_part;
                                        m.worst_value = value;
                                        m.worst_points = rule.n_points();
                                        m.worst_splits = stats.splits;
                                        m.worst_uncertified = m.uncertified - uncertified_before;
                                    }
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
                    if (cfg.diagnose)
                        std::printf("     %ld level-2 boxes in the diagonal frame\n", m.rotations);
                    if (cfg.diagnose)
                        std::printf("     worst cell %d: exact %.6e, rule %.6e, %ld points, %ld bisections, %ld uncertified\n",
                                    m.worst_cell, m.worst_exact, m.worst_value, m.worst_points, m.worst_splits,
                                    m.worst_uncertified);
                    if (cfg.diagnose && !m.uncertified_cells.empty())
                    {
                        std::printf("     %zu cells with uncertified boxes:", m.uncertified_cells.size());
                        for (std::size_t i = 0; i < std::min<std::size_t>(12, m.uncertified_cells.size()); ++i)
                            std::printf(" %d", m.uncertified_cells[i]);
                        std::printf("\n");
                    }
                    if (!causes.empty())
                    {
                        std::printf("     bisections by cause:\n");
                        for (const auto& [key, count] : causes)
                            std::printf("       %7ld  %s\n", count, key.c_str());
                    }
                    if (want_leaves)
                    {
                        const char* slug = kind == PartKind::negative   ? "negative"
                                           : kind == PartKind::positive ? "positive"
                                           : kind == PartKind::interface ? "interface"
                                                                         : "whole";
                        std::string name = gen;
                        std::replace(name.begin(), name.end(), ':', '_');
                        const std::string path = cfg.leaves + "_" + cfg.mesh + "_" + name + "_" + slug + "_n" + std::to_string(n)
                                                 + "_p" + std::to_string(cfg.leaf_degree) + ".vtu";
                        write_leaf_mesh(path, leaf_mesh);
                        write_leaf_mesh(cfg.leaves + "_" + cfg.mesh + "_cut_cells_n" + std::to_string(n) + ".vtu", cut_cells);
                        std::printf("     leaves: %d cells, %d nodes, %ld incomplete -> %s\n", leaf_mesh.n_cells(),
                                    leaf_mesh.n_points(), incomplete, path.c_str());
                    }
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
