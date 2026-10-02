// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Compares quadrature rule generators on cut cells against exact per-cell values.
//
// Test problem: sphere |x - centre| = radius in [-1, 1]^3, meshed with n^3 hexahedra
// or 6 n^3 Kuhn tetrahedra. Parts are selected with the same expressions as
// HOCutResult, e.g. "phi < 0", "phi > 0", "phi = 0". The level set reaches
// quadrays as Bernstein coefficients on each cell, as from a finite-element
// level set; algoim's generators interpolate the sphere on their boxes.
//
// Usage:
//   quadrays_study [--mesh hex|tet] [--n 16,32] [--q 3,5] [--centre x,y,z]
//                  [--radius r] [--gen quadrays,quadrays:0.1,...]
//                  [--part "phi < 0"]... [--csv file] [--plane] [--vtk prefix]
//                  [--leaves prefix] [--leaf-degree p]
//                  [--diagnose] [--only cell] [--masks M] [--no-diagonal]
//
// Generators: quadrays (margin 0.25) or quadrays:<margin>; with
// CUTCELLS_WITH_ALGOIM also algoim-auto, algoim-gl, gl-cellmask, alpha, split,
// alpha-split and quadgen (hex only).
//
// Metrics per run: per-cell L1 (sum over cut cells of |cell error|, relative to
// the ball volume or the sphere area), worst cell (relative to the cell's exact
// value), the global error, points and microseconds per cut cell. The time of a
// quadrays cell covers the conversion of its Bernstein coefficients to the box
// and the engine; the coefficients themselves come from the front end.
//
// --diagnose prints why bisections happen and the worst cell; --only restricts
// the run to one cell index; --masks and --no-diagonal set
// Options::mask_subdivisions and diagonal_frames. --leaves writes the leaf cells
// (Lagrange cells of degree p, 3 by default) of every cut cell, plus the uncut
// cells of volume parts as linear cells. --plane uses phi = x + 0.3 y - 0.2 z:
// every rule must then be exact (volume below 4, area 4 sqrt(1.13)).

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

#include <cutcells/quadrays/leaves.h>
#include <cutcells/quadrays/rules.h>
#include <cutcells/selection_expr.h>

#include "exact_reference.h"
#include "test_mesh.h"

#ifdef CUTCELLS_WITH_ALGOIM
#include "generators.h"
#endif

using namespace cutcells;
using namespace cutcells::quadrays;
using namespace cutcells::quadrays::support;

namespace
{
struct StudyConfig
{
    std::string mesh = "tet";
    std::vector<int> ns = {16, 32};
    std::vector<int> qs = {3};
    V3 centre = {0.0123, -0.0371, 0.0217};
    double radius = 0.7;
    std::vector<std::string> generators = {"quadrays"};
    std::vector<std::string> parts;
    std::string csv;
    bool plane = false;
    std::string vtk;
    std::string leaves;
    int leaf_degree = 3;
    bool diagnose = false;
    int only = -1;
    int masks = 1;
    bool diagonal = true;
};

struct Metrics
{
    long cut_cells = 0;
    long points = 0;
    long bisections = 0;
    long uncertified = 0;
    long rotations = 0;
    double sum_abs = 0;   ///< sum over cut cells of |generated - exact|
    double sum_gen = 0;   ///< generated part measure over cut cells
    double sum_uncut = 0; ///< exact part measure over uncut cells
    double max_rel = 0;   ///< worst cell, relative to the cell's exact value
    double seconds = 0;
    std::int32_t worst_cell = -1;
    double worst_exact = 0, worst_value = 0;
    long worst_points = 0, worst_bisections = 0, worst_uncertified = 0;
    std::vector<std::int32_t> uncertified_cells;
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

bool is_quadrays(const std::string& gen) { return gen.rfind("quadrays", 0) == 0; }

Options quadrays_options(const std::string& gen, const StudyConfig& cfg)
{
    Options opt;
    if (gen.find(':') != std::string::npos)
        opt.margin = std::stod(gen.substr(gen.find(':') + 1));
    opt.diagnose = cfg.diagnose;
    opt.mask_subdivisions = cfg.masks;
    opt.diagonal_frames = cfg.diagonal;
    return opt;
}

const char* part_slug(Part part)
{
    switch (part)
    {
    case Part::negative:
        return "negative";
    case Part::positive:
        return "positive";
    case Part::interface:
        return "interface";
    default:
        return "whole";
    }
}

/// Quadrature points as a VTK point cloud (ASCII .vtu) with weights and cells.
void write_point_cloud(const std::string& path, const std::vector<double>& points,
                       const std::vector<double>& weights, const std::vector<std::int32_t>& cells)
{
    const std::size_t n = weights.size();
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("write_point_cloud: cannot open " + path);
    out.precision(17);
    out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
           "byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n<Piece NumberOfPoints=\""
        << n << "\" NumberOfCells=\"" << n << "\">\n<PointData Scalars=\"weight\">\n"
        << "<DataArray type=\"Float64\" Name=\"weight\" format=\"ascii\">\n";
    for (double w : weights)
        out << w << '\n';
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"cell\" format=\"ascii\">\n";
    for (std::int32_t c : cells)
        out << c << '\n';
    out << "</DataArray>\n</PointData>\n<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" "
           "format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i)
        out << points[3 * i] << ' ' << points[3 * i + 1] << ' ' << points[3 * i + 2] << '\n';
    out << "</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i)
        out << i << '\n';
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i)
        out << i + 1 << '\n';
    out << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i)
        out << "1\n";
    out << "</DataArray>\n</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
}

/// One cell through quadrays: rule in reference coordinates (box coordinates of
/// these test cells) with physical weights.
void quadrays_rule(const TestCell& cell, int degree, const std::vector<double>& coeffs,
                   const SelectionTerm& term, int q, const Options& opt,
                   quadrature::QuadratureRules<double>& rule, Stats& stats)
{
    append_cell_rules<double>(cell.type, cell.vertices, degree, coeffs, term, 0, q, opt, 0, rule, stats);
}

/// Planar level set: every generator must reproduce the totals to rounding.
int plane_study(const StudyConfig& cfg)
{
    auto phi = [](const V3& x) { return x[0] + 0.3 * x[1] - 0.2 * x[2]; };
    std::printf("Plane x + 0.3 y - 0.2 z: global relative errors of the totals.\n");
    std::vector<double> coeffs;
    for (int n : cfg.ns)
        for (int q : cfg.qs)
            for (const std::string& gen : cfg.generators)
            {
                if (gen == "quadgen")
                    continue;
                for (const std::string& part_text : cfg.parts)
                {
                    SelectionExpr expr = parse_selection_expr(part_text);
                    compile_selection_expr(expr, {"phi"});
                    const SelectionTerm& term = expr.terms.front();
                    const Part part = part_of(term);
                    const double exact = part == Part::interface ? 4.0 * std::sqrt(1.13) : 4.0;
                    const double h = 2.0 / n;
                    double total = 0;
                    long cut = 0;
                    for (int i0 = 0; i0 < n; ++i0)
                        for (int i1 = 0; i1 < n; ++i1)
                            for (int i2 = 0; i2 < n; ++i2)
                            {
                                const V3 lo = {-1 + h * i0, -1 + h * i1, -1 + h * i2};
                                for (const TestCell& cell : grid_cells(cfg.mesh, lo, h, {0, 0, 0}))
                                {
                                    // vertex signs decide cut / inside / outside for a plane
                                    bool pos = false, neg = false;
                                    for (std::size_t v = 0; v < cell.vertices.size(); v += 3)
                                    {
                                        const double val = phi({cell.vertices[v], cell.vertices[v + 1], cell.vertices[v + 2]});
                                        pos |= val > 0;
                                        neg |= val < 0;
                                    }
                                    if (!(pos && neg))
                                    {
                                        if ((part == Part::negative && neg) || (part == Part::positive && pos))
                                            total += cell.volume;
                                        continue;
                                    }
                                    ++cut;
                                    quadrature::QuadratureRules<double> rule;
                                    if (is_quadrays(gen))
                                    {
                                        Stats stats;
                                        cell_coefficients(cell, 1, phi, coeffs);
                                        quadrays_rule(cell, 1, coeffs, term, q, quadrays_options(gen, cfg), rule, stats);
                                    }
#ifdef CUTCELLS_WITH_ALGOIM
                                    else
                                    {
                                        benchmarks::LevelSet ls;
                                        ls.degree = 1;
                                        ls.value = [&](const Vec3<double>& x) { return phi({x[0], x[1], x[2]}); };
                                        ClippedBox<double> box;
                                        make_clipped_box<double>(cell.type, cell.vertices, 3, box);
                                        benchmarks::GeneratorStats gstats;
                                        benchmarks::algoim_clipped_box(box, ls, term, q, benchmarks::generator_preset(gen),
                                                                       rule, gstats);
                                    }
#else
                                    else
                                        throw std::runtime_error("generator " + gen + " needs CUTCELLS_WITH_ALGOIM");
#endif
                                    for (double w : rule._weights)
                                        total += w;
                                }
                            }
                    const double exact_part = part == Part::positive ? 4.0 : exact;
                    std::printf("%-4s %4d %2d %-14s %-8s | cut %7ld | global %.1e\n", cfg.mesh.c_str(), n, q,
                                gen.c_str(), part_text.c_str(), cut, std::abs(total - exact_part) / exact_part);
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
    const V3 c = cfg.centre;
    const double ball = 4.0 / 3.0 * M_PI * r * r * r, sphere = 4.0 * M_PI * r * r;
    auto phi = [&](const V3& x)
    { return (x[0] - c[0]) * (x[0] - c[0]) + (x[1] - c[1]) * (x[1] - c[1]) + (x[2] - c[2]) * (x[2] - c[2]) - r * r; };

    std::ofstream csv;
    if (!cfg.csv.empty())
    {
        csv.open(cfg.csv);
        csv << "mesh,n,q,generator,part,cut_cells,points_per_cut_cell,global_rel,l1_rel,max_cell_rel,bisections,"
               "uncertified,us_per_cut_cell\n";
    }
    std::printf("Sphere centre (%g, %g, %g), radius %g. Errors relative to the ball volume (volume parts)\n"
                "or the sphere area (interface); per-cell L1 = sum over cut cells of |cell error|.\n",
                c[0], c[1], c[2], r);
    std::printf("%-4s %4s %2s %-14s %-8s | %7s | %8s | %9s %9s %9s | %6s %6s | %8s\n", "mesh", "n", "q", "generator",
                "part", "cut", "pts/cut", "global", "L1", "worst", "bisect", "uncert", "us/cut");

    std::vector<double> coeffs;
    for (int n : cfg.ns)
        for (int q : cfg.qs)
            for (const std::string& gen : cfg.generators)
            {
                if (gen == "quadgen" && cfg.mesh != "hex")
                {
                    std::printf("(quadgen skipped: hexahedra only)\n");
                    continue;
                }
                const bool quad = is_quadrays(gen);
#ifndef CUTCELLS_WITH_ALGOIM
                if (!quad)
                {
                    std::fprintf(stderr, "generator %s needs CUTCELLS_WITH_ALGOIM\n", gen.c_str());
                    return 1;
                }
#endif
                const Options opt = quadrays_options(gen, cfg);
                for (const std::string& part_text : cfg.parts)
                {
                    SelectionExpr expr = parse_selection_expr(part_text);
                    compile_selection_expr(expr, {"phi"});
                    if (expr.terms.size() != 1)
                        throw std::runtime_error("one selection term per part");
                    const SelectionTerm& term = expr.terms.front();
                    const Part part = part_of(term);
                    const bool surface = part == Part::interface;
                    const double scale = surface ? sphere : ball;
                    const double total = part == Part::negative    ? ball
                                         : part == Part::positive  ? 8.0 - ball
                                         : part == Part::interface ? sphere
                                                                   : 8.0;

                    const double h = 2.0 / n;
                    Metrics m;
                    Stats stats_all;
                    std::vector<double> vtk_points, vtk_weights;
                    std::vector<std::int32_t> vtk_cells;
                    std::int32_t cell_index = -1;
                    const bool want_leaves = !cfg.leaves.empty() && quad;
                    LeafMesh<double> leaf_mesh, cut_cells;
                    for (int i0 = 0; i0 < n; ++i0)
                        for (int i1 = 0; i1 < n; ++i1)
                            for (int i2 = 0; i2 < n; ++i2)
                            {
                                const V3 lo = {-1 + h * i0, -1 + h * i1, -1 + h * i2};
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
                                    const double exact_part = part == Part::negative    ? inside
                                                              : part == Part::positive  ? cell.volume - inside
                                                              : part == Part::interface ? area
                                                                                        : cell.volume;
                                    if (area <= 0.0)
                                    {
                                        m.sum_uncut += exact_part;
                                        if (want_leaves && !surface && exact_part > 0.5 * cell.volume)
                                            append_linear_cell<double>(cell.type, cell.vertices, cell_index, leaf_mesh);
                                        continue;
                                    }
                                    if (quad)
                                        cell_coefficients(cell, 2, phi, coeffs);
                                    if (want_leaves)
                                    {
                                        ClippedBox<double> box;
                                        BoxBernstein<double> form;
                                        make_clipped_box<double>(cell.type, cell.vertices, 3, box);
                                        cell_bernstein_on_box<double>(cell.type, 2, coeffs, form);
                                        Stats leaf_stats;
                                        append_leaves<double>(box, form, part, cfg.leaf_degree, opt, cell_index, leaf_mesh,
                                                              leaf_stats);
                                        stats_all.incomplete_leaves += leaf_stats.incomplete_leaves;
                                        append_linear_cell<double>(cell.type, cell.vertices, cell_index, cut_cells);
                                    }
                                    quadrature::QuadratureRules<double> rule;
                                    Stats stats;
                                    const auto t0 = std::chrono::steady_clock::now();
                                    if (quad)
                                        quadrays_rule(cell, 2, coeffs, term, q, opt, rule, stats);
#ifdef CUTCELLS_WITH_ALGOIM
                                    else
                                    {
                                        ClippedBox<double> box;
                                        make_clipped_box<double>(cell.type, cell.vertices, 3, box);
                                        if (gen == "quadgen")
                                            benchmarks::algoim_quadgen_sphere(box, {c[0], c[1], c[2]}, r, term, q, rule);
                                        else
                                        {
                                            benchmarks::LevelSet ls;
                                            ls.degree = 2;
                                            ls.value = [&](const Vec3<double>& x) { return phi({x[0], x[1], x[2]}); };
                                            benchmarks::GeneratorStats gstats;
                                            benchmarks::algoim_clipped_box(box, ls, term, q,
                                                                           benchmarks::generator_preset(gen), rule, gstats);
                                            stats.bisections += gstats.splits;
                                        }
                                    }
#endif
                                    m.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                                    double value = 0;
                                    for (double w : rule._weights)
                                        value += w;
                                    const int n_points = static_cast<int>(rule._weights.size());
                                    if (!cfg.vtk.empty())
                                        for (int i = 0; i < n_points; ++i)
                                        {
                                            const V3 x = physical(cell, {rule._points[3 * i], rule._points[3 * i + 1],
                                                                         rule._points[3 * i + 2]});
                                            vtk_points.insert(vtk_points.end(), x.begin(), x.end());
                                            vtk_weights.push_back(rule._weights[i]);
                                            vtk_cells.push_back(cell_index);
                                        }
                                    if (stats.uncertified > 0)
                                        m.uncertified_cells.push_back(cell_index);
                                    ++m.cut_cells;
                                    m.points += n_points;
                                    m.bisections += stats.bisections;
                                    m.uncertified += stats.uncertified;
                                    m.rotations += stats.rotations;
                                    for (const auto& [key, count] : stats.causes)
                                        stats_all.causes[key] += count;
                                    m.sum_gen += value;
                                    m.sum_abs += std::abs(value - exact_part);
                                    const double floor = surface ? 1e-3 * h * h : 1e-3 * h * h * h;
                                    if (exact_part > floor && std::abs(value - exact_part) / exact_part > m.max_rel)
                                    {
                                        m.max_rel = std::abs(value - exact_part) / exact_part;
                                        m.worst_cell = cell_index;
                                        m.worst_exact = exact_part;
                                        m.worst_value = value;
                                        m.worst_points = n_points;
                                        m.worst_bisections = stats.bisections;
                                        m.worst_uncertified = stats.uncertified;
                                    }
                                }
                            }
                    const double global = std::abs(m.sum_gen + m.sum_uncut - total) / scale;
                    const double l1 = m.sum_abs / scale;
                    const double pts = m.cut_cells ? static_cast<double>(m.points) / m.cut_cells : 0.0;
                    const double us = m.cut_cells ? m.seconds / m.cut_cells * 1e6 : 0.0;
                    std::printf("%-4s %4d %2d %-14s %-8s | %7ld | %8.1f | %9.1e %9.1e %9.1e | %6ld %6ld | %8.1f\n",
                                cfg.mesh.c_str(), n, q, gen.c_str(), part_text.c_str(), m.cut_cells, pts, global, l1,
                                m.max_rel, m.bisections, m.uncertified, us);
                    std::fflush(stdout);
                    if (cfg.diagnose)
                    {
                        std::printf("     %ld level-2 boxes in the diagonal frame\n", m.rotations);
                        std::printf("     worst cell %d: exact %.6e, rule %.6e, %ld points, %ld bisections, %ld "
                                    "uncertified\n",
                                    m.worst_cell, m.worst_exact, m.worst_value, m.worst_points, m.worst_bisections,
                                    m.worst_uncertified);
                        if (!m.uncertified_cells.empty())
                        {
                            std::printf("     %zu cells with uncertified boxes:", m.uncertified_cells.size());
                            for (std::size_t i = 0; i < std::min<std::size_t>(12, m.uncertified_cells.size()); ++i)
                                std::printf(" %d", m.uncertified_cells[i]);
                            std::printf("\n");
                        }
                        if (!stats_all.causes.empty())
                        {
                            std::printf("     bisections by cause:\n");
                            for (const auto& [key, count] : stats_all.causes)
                                std::printf("       %7lld  %s\n", static_cast<long long>(count), key.c_str());
                        }
                    }
                    std::string name = gen;
                    std::replace(name.begin(), name.end(), ':', '_');
                    if (want_leaves)
                    {
                        const std::string path = cfg.leaves + "_" + cfg.mesh + "_" + name + "_" + part_slug(part) + "_n"
                                                 + std::to_string(n) + "_p" + std::to_string(cfg.leaf_degree) + ".vtu";
                        write_leaves(path, leaf_mesh);
                        write_leaves(cfg.leaves + "_" + cfg.mesh + "_cut_cells_n" + std::to_string(n) + ".vtu", cut_cells);
                        std::printf("     leaves: %d cells, %d nodes, %d incomplete -> %s\n", leaf_mesh.n_cells(),
                                    leaf_mesh.n_points(), stats_all.incomplete_leaves, path.c_str());
                    }
                    if (!cfg.vtk.empty())
                        write_point_cloud(cfg.vtk + "_" + cfg.mesh + "_" + name + "_" + part_slug(part) + "_n"
                                              + std::to_string(n) + "_q" + std::to_string(q) + ".vtu",
                                          vtk_points, vtk_weights, vtk_cells);
                    if (csv)
                        csv << cfg.mesh << ',' << n << ',' << q << ',' << gen << ",\"" << part_text << "\"," << m.cut_cells
                            << ',' << pts << ',' << global << ',' << l1 << ',' << m.max_rel << ',' << m.bisections << ','
                            << m.uncertified << ',' << us << '\n';
                }
            }
    return 0;
}
