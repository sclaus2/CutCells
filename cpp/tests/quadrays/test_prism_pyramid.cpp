// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// quadrays on prisms (two per grid cube) and pyramids (six per grid cube, apex
// at its centre) of [-1, 1]^3, for Bernstein coefficients (a prism's basis is
// the triangle's times the interval's; a pyramid's spans its rational Lagrange
// space, and the engine reads it on the cube of collapsed coordinates) and for
// analytic level sets.
//  - Planes in general position are exact on every cell, cut or not, as
//    degree-1 and degree-2 coefficients and analytically: per cell relative to
//    h^3 or h^2, and in total, to 1e-13.
//  - The sphere of radius 0.7 off the grid's symmetry converges with q: on
//    n = 8, the per-cell L1 of the ball's volume and the sphere's area,
//    relative to their totals, stays below the bounds given for q = 3, 5, 8.
//  - The robustness cases of support/robustness_cases.h, as test_robustness
//    runs them on hexahedra and tetrahedra (q = 3, n = 8 for batch 1, n = 4 for
//    batch 2, the same checks; the thin shells and the double root left out;
//    Bernstein forms on pyramids n = 4 throughout: they bisect towards the apex
//    where a zero set passes near it), with a per-cell L1 below 2e-4 for the
//    placements and planes: the sphere tangent to grid planes at vertices has
//    1.3e-4 on prisms (8.4e-5 on hexahedra).
// The argument "bernstein" or "analytic" runs only those sources.
// Exits non-zero on failure.

#include <cutcells/quadrays/analytic.h>
#include <cutcells/quadrays/rules.h>
#include <cutcells/selection_expr.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "support/exact_reference.h"
#include "support/robustness_cases.h"
#include "support/test_mesh.h"

using namespace cutcells;
using namespace cutcells::quadrays;
using namespace cutcells::quadrays::support;

namespace
{
struct Plane
{
    V3 a;
    double b;

    template <typename V>
    V operator()(const std::array<V, 3>& x) const
    {
        return x[0] * a[0] + x[1] * a[1] + x[2] * a[2] - b;
    }
};

SelectionTerm term_of(const char* text)
{
    SelectionExpr expr = parse_selection_expr(text);
    compile_selection_expr(expr, {"phi"});
    return expr.terms.front();
}
} // namespace

int main(int argc, char** argv)
{
    int failures = 0;
    const Options opt;
    const std::string only = argc > 1 ? argv[1] : "";
    std::vector<bool> sources; // analytic or not
    if (only != "analytic")
        sources.push_back(false);
    if (only != "bernstein")
        sources.push_back(true);

    // planes: exact on every cell
    {
        double worst_overall = 0;
        const std::vector<Plane> planes = {{{1.0, 0.3, -0.2}, 0.0}, {{-0.45, 1.0, 0.7}, 0.1}};
        const std::vector<std::pair<const char*, int>> meshes = {{"prism", 5}, {"pyramid", 4}};
        std::vector<double> coeffs;
        for (const Plane& plane : planes)
            for (const auto& [mesh, n] : meshes)
                for (const int source : {1, 2, 0}) // Bernstein of degree 1 and 2, analytic
                {
                    if (!only.empty() && (source == 0) != (only == "analytic"))
                        continue;
                    for (const int q : {3, 5})
                        for (const char* part : {"phi < 0", "phi > 0", "phi = 0"})
                        {
                            const SelectionTerm term = term_of(part);
                            const Part p = part_of(term);
                            const AnalyticLevelSet ls = analytic_level_set(plane);
                            const double h = 2.0 / n;
                            double total = 0, exact_total = 0, worst = 0;
                            for (int i0 = 0; i0 < n; ++i0)
                                for (int i1 = 0; i1 < n; ++i1)
                                    for (int i2 = 0; i2 < n; ++i2)
                                        for (const TestCell& cell :
                                             grid_cells(mesh, {-1 + h * i0, -1 + h * i1, -1 + h * i2}, h, {0, 0, 0}))
                                        {
                                            const exact::PlaneCut cut = exact::plane_cut(cell.faces, plane.a, plane.b);
                                            const double exact_value = p == Part::negative ? cut.volume_below
                                                                       : p == Part::positive
                                                                           ? cell.volume - cut.volume_below
                                                                           : cut.cut_area;
                                            quadrature::QuadratureRules<double> rule;
                                            Stats stats;
                                            if (source == 0)
                                                append_cell_rules<double>(cell.type, cell.vertices, ls, term, 0, q, opt, 0,
                                                                          rule, stats);
                                            else
                                            {
                                                cell_coefficients(
                                                    cell, source,
                                                    [&plane](const V3& x) { return plane(std::array<double, 3>{x[0], x[1], x[2]}); },
                                                    coeffs);
                                                append_cell_rules<double>(cell.type, cell.vertices, source, coeffs, term, 0,
                                                                          q, opt, 0, rule, stats);
                                            }
                                            double value = 0;
                                            for (double w : rule._weights)
                                                value += w;
                                            total += value;
                                            exact_total += exact_value;
                                            worst = std::max(worst, std::abs(value - exact_value)
                                                                        / (p == Part::interface ? h * h : h * h * h));
                                        }
                            const double rel_total = std::abs(total - exact_total) / exact_total;
                            worst_overall = std::max({worst_overall, worst, rel_total});
                            if (worst > 1e-13 || rel_total > 1e-13)
                            {
                                std::printf("FAILED: plane (%g, %g, %g) . x = %g, %s, %s n = %d, q = %d, %s: worst cell "
                                            "%.1e, total %.1e\n",
                                            plane.a[0], plane.a[1], plane.a[2], plane.b,
                                            source == 0 ? "analytic" : source == 1 ? "P1" : "P2", mesh, n, q, part, worst,
                                            rel_total);
                                ++failures;
                            }
                        }
                }
        std::printf("planes: largest error %.1e (per cell relative to h^3 or h^2, totals relative)\n", worst_overall);
    }

    // the sphere: per-cell L1 against the exact ball and sphere, by q
    {
        // about three times the larger L1 of the two meshes
        struct Bound
        {
            int q;
            const char* part;
            double l1;
        };
        const Bound bounds[] = {{3, "phi < 0", 2e-5}, {3, "phi = 0", 1.5e-4}, {5, "phi < 0", 5e-8},
                                {5, "phi = 0", 1e-6}, {8, "phi < 0", 4e-11},  {8, "phi = 0", 2.5e-9}};
        Case sphere;
        sphere.centre = {0.0123, -0.0371, 0.0217};
        sphere.radius = 0.7;
        for (const bool analytic : sources)
            for (const std::string mesh : {"prism", "pyramid"})
                for (const Bound& b : bounds)
                {
                    const CaseRun r = run_case(sphere, mesh, 8, term_of(b.part), b.q, opt, analytic);
                    const double l1 = r.l1 / r.exact_total;
                    const bool ok = l1 <= b.l1 && r.fail + r.negative + r.outside + r.side == 0;
                    std::printf("sphere %-9s %-7s q = %d %-8s L1 %.1e (bound %.1e), max bisections per cell %d%s\n",
                                analytic ? "analytic" : "P2", mesh.c_str(), b.q, b.part, l1, b.l1, r.max_bisections,
                                ok ? "" : " FAILED");
                    failures += !ok;
                }
    }

    // robustness cases
    std::map<std::string, double> baseline;
    for (const bool analytic : sources)
    for (const Case& c : all_cases())
    {
        if (c.shape == Shape::double_root || c.shape == Shape::shell)
            continue; // see test_two_roots and the report
        const bool batch2 = c.shape != Shape::sphere && c.shape != Shape::plane;
        for (const std::string mesh : {"prism", "pyramid"})
            for (const char* part : {"phi < 0", "phi = 0"})
            {
                // a pyramid's Bernstein form bisects towards the apex where a
                // zero set passes near it: batch 1 on n = 4 there too
                const int n = batch2 || (mesh == "pyramid" && !analytic) ? 4 : 8;
                const SelectionTerm term = term_of(part);
                const bool surface = part_of(term) == Part::interface;
                const CaseRun r = run_case(c, mesh, n, term, 3, opt, analytic);
                std::vector<std::string> problems;
                if (r.fail + r.negative + r.outside + r.side > 0)
                    problems.push_back("fail/neg/out/side " + std::to_string(r.fail) + "/" + std::to_string(r.negative)
                                       + "/" + std::to_string(r.outside) + "/" + std::to_string(r.side));
                if (r.max_bisections > opt.max_bisections)
                    problems.push_back(std::to_string(r.max_bisections) + " bisections in one cell");
                const bool known = c.shape == Shape::sphere || c.shape == Shape::plane || c.shape == Shape::two_spheres;
                const double l1 = r.exact_total > 0 ? r.l1 / r.exact_total : r.l1;
                if (known && !batch2 && l1 > 2e-4)
                    problems.push_back("per-cell L1 " + std::to_string(l1));
                const std::string key = (analytic ? "analytic " : "") + mesh + part;
                if (c.name == "sphere")
                    baseline[key] = r.total;
                if (c.scale != 1.0 && std::abs(r.total - baseline[key]) > 1e-12 * baseline[key])
                    problems.push_back("total differs from the unscaled sphere");
                if (c.shape == Shape::cone || c.shape == Shape::torus)
                {
                    const double exact = exact_totals(c)[surface ? 1 : 0];
                    const double tol = c.shape == Shape::cone && surface ? 2e-2 : 1e-3;
                    if (std::abs(r.total - exact) > tol * exact)
                        problems.push_back("total off by " + std::to_string(std::abs(r.total - exact) / exact));
                }
                std::string message;
                for (const std::string& p : problems)
                    message += " " + p + ";";
                std::printf("%-9s %-22s %-7s %-8s n = %d: L1 %.1e, max bisections per cell %d%s%s\n",
                            analytic ? "analytic" : "Bernstein", c.name.c_str(), mesh.c_str(), part, n,
                            known ? l1 : 0.0, r.max_bisections, problems.empty() ? "" : " FAILED:", message.c_str());
                failures += !problems.empty();
            }
    }
    return failures == 0 ? 0 : 1;
}
