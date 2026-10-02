// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Robustness report (docs/quadrays/RESULTS.md, v1.2). Every cell of [-1, 1]^3,
// cut or not, goes to every generator for every case of
// cpp/tests/quadrays/support/robustness_cases.h and is checked against exact
// per-cell values; for cones, the double root and the torus only the totals are
// known, so L1, max and bad stay empty.
//
// Per run (case, mesh, q, generator, part):
//   fail   cells where the generator threw, or returned non-finite points or weights
//   neg    negative weights
//   out    points outside the cell (beyond 1e-12 in reference coordinates)
//   side   volume points on the wrong side of the level set, or interface points off
//          it, by more than 1e-9 h
//   L1     sum over cells of |error|, relative to the exact total
//   max    largest cell error, relative to h^3 (volume) or h^2 (interface)
//   bad    cells with an error above 1e-4 h^3 (volume) or 1e-4 h^2 (interface)
//   total  error of the total over the mesh, relative to the exact total
//   us     mean and largest time per cell, over all cells
// The interface of a plane lying in faces has no well-defined owner: those cells
// are left out of L1, max and bad, and the exact total counts each such face once.
//
// Usage: quadrays_robustness_report [--case a,b] [--mesh tet,hex] [--n 16] [--q 3]
//                                   [--gen quadrays,quadrays:0.1,...]
//                                   [--part "phi < 0"]... [--list]
// Generators: quadrays, quadrays:<margin>; with CUTCELLS_WITH_ALGOIM also
// algoim-auto, algoim-gl and quadgen (hexahedra, unscaled spheres).

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <cutcells/quadrays/rules.h>
#include <cutcells/selection_expr.h>

#include "robustness_cases.h"
#include "test_mesh.h"

#ifdef CUTCELLS_WITH_ALGOIM
#include "generators.h"
#endif

using namespace cutcells;
using namespace cutcells::quadrays;
using namespace cutcells::quadrays::support;

namespace
{
struct Tally
{
    long cut = 0, fail = 0, negative = 0, outside = 0, side = 0, bad = 0, points = 0;
    long bisections = 0, uncertified = 0;
    double l1 = 0, max_error = 0, total = 0, exact_total = 0, seconds = 0, max_seconds = 0;
    long cells = 0;
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
} // namespace

int main(int argc, char** argv)
{
    std::vector<std::string> case_names, meshes = {"tet", "hex"};
    std::vector<std::string> generators = {"quadrays", "quadrays:0.1"};
#ifdef CUTCELLS_WITH_ALGOIM
    generators.insert(generators.end(), {"algoim-auto", "algoim-gl", "quadgen"});
#endif
    std::vector<std::string> parts;
    std::vector<int> qs = {3};
    int n_cells = 16;
    const std::vector<Case> cases = all_cases();
    try
    {
        for (int i = 1; i < argc; ++i)
        {
            const std::string a = argv[i];
            auto next = [&]() -> std::string
            {
                if (i + 1 >= argc)
                    throw std::runtime_error("missing value after " + a);
                return argv[++i];
            };
            if (a == "--case")
                case_names = split_list(next());
            else if (a == "--mesh")
                meshes = split_list(next());
            else if (a == "--gen")
                generators = split_list(next());
            else if (a == "--part")
                parts.push_back(next());
            else if (a == "--n")
                n_cells = std::stoi(next());
            else if (a == "--q")
            {
                qs.clear();
                for (const std::string& s : split_list(next()))
                    qs.push_back(std::stoi(s));
            }
            else if (a == "--list")
            {
                for (const Case& c : cases)
                    std::printf("%-22s %s\n", c.name.c_str(), c.description.c_str());
                return 0;
            }
            else
                throw std::runtime_error("unknown argument: " + a);
        }
    }
    catch (const std::exception& e)
    {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }
    if (parts.empty())
        parts = {"phi < 0", "phi = 0"};
    const double h = 2.0 / n_cells;

    std::printf("%-22s %-4s %2s %-14s %-8s | %5s | %4s %4s %4s %5s | %8s %8s %5s | %8s | %7s | %7s %8s | %s\n", "case",
                "mesh", "q", "generator", "part", "cut", "fail", "neg", "out", "side", "L1", "max", "bad", "total",
                "pts/cut", "us/cell", "max us", "bisect/uncert");
    std::vector<double> coeffs;
    for (const Case& c : cases)
    {
        if (!case_names.empty() && std::find(case_names.begin(), case_names.end(), c.name) == case_names.end())
            continue;
        const std::array<double, 2> totals = exact_totals(c);
        for (const std::string& mesh : meshes)
            for (int q : qs)
                for (const std::string& gen : generators)
                {
                    const bool quad = gen.rfind("quadrays", 0) == 0;
#ifndef CUTCELLS_WITH_ALGOIM
                    if (!quad)
                    {
                        std::fprintf(stderr, "generator %s needs CUTCELLS_WITH_ALGOIM\n", gen.c_str());
                        return 1;
                    }
#endif
                    if (gen == "quadgen" && (mesh != "hex" || c.shape != Shape::sphere || c.scale != 1.0))
                        continue; // algoim's 2015 engine here takes the sphere itself, on boxes
                    Options opt;
                    if (quad && gen.find(':') != std::string::npos)
                        opt.margin = std::stod(gen.substr(gen.find(':') + 1));
                    for (const std::string& part_text : parts)
                    {
                        SelectionExpr expr = parse_selection_expr(part_text);
                        compile_selection_expr(expr, {"phi"});
                        const SelectionTerm& term = expr.terms.front();
                        const Part part = part_of(term);
                        const bool surface = part == Part::interface;
                        if (part != Part::negative && part != Part::interface)
                            throw std::runtime_error("robustness report: parts phi < 0 and phi = 0 only");
                        const double cell_scale = surface ? h * h : h * h * h;
                        auto generate = [&](const TestCell& cell, quadrature::QuadratureRules<double>& rule, Stats& stats)
                        {
                            if (quad)
                            {
                                cell_coefficients(cell, degree(c), [&c](const V3& x) { return phi(c, x); }, coeffs);
                                append_cell_rules<double>(cell.type, cell.vertices, degree(c), coeffs, term, 0, q, opt, 0,
                                                          rule, stats);
                                return;
                            }
#ifdef CUTCELLS_WITH_ALGOIM
                            ClippedBox<double> box;
                            make_clipped_box<double>(cell.type, cell.vertices, 3, box);
                            if (gen == "quadgen")
                                benchmarks::algoim_quadgen_sphere(box, {c.centre[0], c.centre[1], c.centre[2]}, c.radius,
                                                                  term, q, rule);
                            else
                            {
                                benchmarks::LevelSet ls;
                                ls.degree = degree(c);
                                ls.value = [&c](const Vec3<double>& x) { return phi(c, {x[0], x[1], x[2]}); };
                                benchmarks::GeneratorStats gstats;
                                benchmarks::algoim_clipped_box(box, ls, term, q, benchmarks::generator_preset(gen), rule,
                                                               gstats);
                            }
#endif
                        };
                        Tally t;
                        for (int i0 = 0; i0 < n_cells; ++i0)
                            for (int i1 = 0; i1 < n_cells; ++i1)
                                for (int i2 = 0; i2 < n_cells; ++i2)
                                {
                                    const V3 lo = {-1 + h * i0, -1 + h * i1, -1 + h * i2};
                                    for (const TestCell& cell : grid_cells(mesh, lo, h, {0, 0, 0}))
                                    {
                                        ++t.cells;
                                        const Reference ref = reference(c, cell);
                                        const bool on_face = surface && ref.face_area > 0;
                                        const double exact = surface ? ref.area : ref.volume;
                                        bool cut = false;
                                        if (ref.known)
                                        {
                                            t.exact_total += surface ? ref.area + 0.5 * ref.face_area : ref.volume;
                                            cut = surface ? (ref.area > 0 || on_face)
                                                          : (ref.volume > 0 && ref.volume < cell.volume);
                                        }
                                        else
                                        {
                                            // a sign change at 5^3 points of the cell's box, inside the cell
                                            bool pos = false, neg = false;
                                            for (int a = 0; a < 125; ++a)
                                            {
                                                const V3 xi = {(a % 5) / 4.0, (a / 5 % 5) / 4.0, (a / 25) / 4.0};
                                                if (cell.type == cell::type::tetrahedron && xi[0] + xi[1] + xi[2] > 1)
                                                    continue;
                                                const double f = phi(c, physical(cell, xi));
                                                pos |= f > 0;
                                                neg |= f < 0;
                                            }
                                            cut = pos && neg;
                                        }
                                        t.cut += cut;

                                        quadrature::QuadratureRules<double> rule;
                                        Stats stats;
                                        bool failed = false;
                                        const auto t0 = std::chrono::steady_clock::now();
                                        try
                                        {
                                            generate(cell, rule, stats);
                                        }
                                        catch (const std::exception&)
                                        {
                                            failed = true;
                                        }
                                        const double dt
                                            = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                                        t.seconds += dt;
                                        t.max_seconds = std::max(t.max_seconds, dt);
                                        t.bisections += stats.bisections;
                                        t.uncertified += stats.uncertified;
                                        if (cut)
                                            t.points += static_cast<long>(rule._weights.size());
                                        const RuleCheck check = failed ? RuleCheck{true} : check_rule(c, cell, rule, surface, h);
                                        if (check.fail)
                                        {
                                            ++t.fail;
                                            continue;
                                        }
                                        t.negative += check.negative;
                                        t.outside += check.outside;
                                        t.side += check.side;
                                        t.total += check.value;
                                        if (on_face || !ref.known)
                                            continue;
                                        const double error = std::abs(check.value - exact);
                                        t.l1 += error;
                                        t.max_error = std::max(t.max_error, error / cell_scale);
                                        t.bad += error > 1e-4 * cell_scale;
                                    }
                                }
                        const double exact_total
                            = c.shape == Shape::cone || c.shape == Shape::double_root || c.shape == Shape::torus
                                  ? totals[surface ? 1 : 0]
                                  : t.exact_total;
                        std::printf("%-22s %-4s %2d %-14s %-8s | %5ld | %4ld %4ld %4ld %5ld | %8.1e %8.1e %5ld | %8.1e | "
                                    "%7.1f | %7.1f %8.0f | %ld/%ld\n",
                                    c.name.c_str(), mesh.c_str(), q, gen.c_str(), part_text.c_str(), t.cut, t.fail,
                                    t.negative, t.outside, t.side, exact_total > 0 ? t.l1 / exact_total : t.l1,
                                    t.max_error, t.bad,
                                    exact_total > 0 ? std::abs(t.total - exact_total) / exact_total : std::abs(t.total),
                                    t.cut ? double(t.points) / t.cut : 0.0, t.seconds / t.cells * 1e6,
                                    t.max_seconds * 1e6, t.bisections, t.uncertified);
                        std::fflush(stdout);
                    }
                }
    }
    return 0;
}
