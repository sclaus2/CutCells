// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// The robustness cases of support/robustness_cases.h as pass/fail, the cases
// the quadrays plan lists: placements on the grid, phi scaled by 1e+-150 and
// 1e+-200, planes (in faces, through vertices, slivers), touching spheres,
// cones and the torus. Every cell of [-1, 1]^3, cut or not, goes through
// quadrays (q = 3) on n = 8, and on n = 4 for the degree-4 products, cones and
// the torus. Every run must have no exception, no non-finite value, no negative
// weight, no point outside its cell or on the wrong side, and at most
// Options::max_bisections bisections per cell. Scaled level sets must give the
// unscaled totals, cones and the torus their exact totals (1e-3; 2e-2 for the
// area near a cone's apex on n = 4), and the placements and planes a per-cell
// L1 below 1e-4 of the exact total. On n = 4 the touching spheres are smaller
// than the cells, which then hold two sheets, the engine's known weak spot:
// they take the checks only, as do the double root and the thin shells, which
// are left to benchmarks/quadrays/robustness_report.
// Exits non-zero on failure.

#include <cutcells/quadrays/rules.h>
#include <cutcells/selection_expr.h>

#include <cmath>
#include <cstdio>
#include <exception>
#include <map>
#include <string>
#include <vector>

#include "support/robustness_cases.h"
#include "support/test_mesh.h"

using namespace cutcells;
using namespace cutcells::quadrays;
using namespace cutcells::quadrays::support;

namespace
{
struct Run
{
    long fail = 0, negative = 0, outside = 0, side = 0;
    int max_bisections = 0;
    double total = 0, exact_total = 0, l1 = 0;
};

Run run_case(const Case& c, const std::string& mesh, int n, const SelectionTerm& term, const Options& opt)
{
    const double h = 2.0 / n;
    const bool surface = part_of(term) == Part::interface;
    std::vector<double> coeffs;
    Run r;
    for (int i0 = 0; i0 < n; ++i0)
        for (int i1 = 0; i1 < n; ++i1)
            for (int i2 = 0; i2 < n; ++i2)
                for (const TestCell& cell : grid_cells(mesh, {-1 + h * i0, -1 + h * i1, -1 + h * i2}, h, {0, 0, 0}))
                {
                    const Reference ref = reference(c, cell);
                    const bool on_face = surface && ref.face_area > 0;
                    if (ref.known)
                        r.exact_total += surface ? ref.area + 0.5 * ref.face_area : ref.volume;
                    quadrature::QuadratureRules<double> rule;
                    Stats stats;
                    try
                    {
                        cell_coefficients(cell, degree(c), [&c](const V3& x) { return phi(c, x); }, coeffs);
                        append_cell_rules<double>(cell.type, cell.vertices, degree(c), coeffs, term, 0, 3, opt, 0, rule,
                                                  stats);
                    }
                    catch (const std::exception&)
                    {
                        ++r.fail;
                        continue;
                    }
                    r.max_bisections = std::max(r.max_bisections, stats.bisections);
                    const RuleCheck check = check_rule(c, cell, rule, surface, h);
                    if (check.fail)
                    {
                        ++r.fail;
                        continue;
                    }
                    r.negative += check.negative;
                    r.outside += check.outside;
                    r.side += check.side;
                    r.total += check.value;
                    if (ref.known && !on_face)
                        r.l1 += std::abs(check.value - (surface ? ref.area : ref.volume));
                }
    return r;
}
} // namespace

int main()
{
    const Options opt;
    int failures = 0;
    std::map<std::string, double> baseline; // unscaled sphere totals by mesh and part
    for (const Case& c : all_cases())
    {
        if (c.shape == Shape::double_root || c.shape == Shape::shell)
            continue; // weak spots, see the report
        const bool batch2 = c.shape != Shape::sphere && c.shape != Shape::plane;
        const int n = batch2 ? 4 : 8;
        for (const std::string mesh : {"tet", "hex"})
            for (const std::string part_text : {"phi < 0", "phi = 0"})
            {
                SelectionExpr expr = parse_selection_expr(part_text);
                compile_selection_expr(expr, {"phi"});
                const SelectionTerm& term = expr.terms.front();
                const bool surface = part_of(term) == Part::interface;
                const Run r = run_case(c, mesh, n, term, opt);

                std::vector<std::string> problems;
                if (r.fail + r.negative + r.outside + r.side > 0)
                    problems.push_back("fail/neg/out/side " + std::to_string(r.fail) + "/" + std::to_string(r.negative)
                                       + "/" + std::to_string(r.outside) + "/" + std::to_string(r.side));
                if (r.max_bisections > opt.max_bisections)
                    problems.push_back(std::to_string(r.max_bisections) + " bisections in one cell");
                const bool known = c.shape == Shape::sphere || c.shape == Shape::plane || c.shape == Shape::two_spheres;
                const double l1 = r.exact_total > 0 ? r.l1 / r.exact_total : r.l1;
                if (known && !batch2 && l1 > 1e-4)
                    problems.push_back("per-cell L1 " + std::to_string(l1));
                const std::string key = mesh + part_text;
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
                std::printf("%-22s %-4s %-8s n = %d: L1 %.1e, max bisections per cell %d%s%s\n", c.name.c_str(),
                            mesh.c_str(), part_text.c_str(), n, known ? l1 : 0.0, r.max_bisections,
                            problems.empty() ? "" : " FAILED:", message.c_str());
                failures += !problems.empty();
            }
    }
    return failures == 0 ? 0 : 1;
}
