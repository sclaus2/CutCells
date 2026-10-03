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
// than the cells, which then hold two sheets: they take the checks only. The
// thin shells are in test_two_roots; the double root is left to
// benchmarks/quadrays/robustness_report.
// The cases run twice: as Bernstein coefficients and as analytic level sets
// through quadrays/analytic.h (Taylor-model bounds).
// Exits non-zero on failure.

#include <cutcells/selection_expr.h>

#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "support/robustness_cases.h"

using namespace cutcells;
using namespace cutcells::quadrays;
using namespace cutcells::quadrays::support;

int main()
{
    const Options opt;
    int failures = 0;
    std::map<std::string, double> baseline; // unscaled sphere totals by source, mesh and part
    for (const bool analytic : {false, true})
    {
        for (const Case& c : all_cases())
        {
            if (c.shape == Shape::double_root || c.shape == Shape::shell)
                continue; // see the report and test_two_roots
            const bool batch2 = c.shape != Shape::sphere && c.shape != Shape::plane;
            const int n = batch2 ? 4 : 8;
            for (const std::string mesh : {"tet", "hex"})
                for (const std::string part_text : {"phi < 0", "phi = 0"})
                {
                    SelectionExpr expr = parse_selection_expr(part_text);
                    compile_selection_expr(expr, {"phi"});
                    const SelectionTerm& term = expr.terms.front();
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
                    if (known && !batch2 && l1 > 1e-4)
                        problems.push_back("per-cell L1 " + std::to_string(l1));
                    const std::string key = (analytic ? "analytic " : "") + mesh + part_text;
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
                    std::printf("%-9s %-22s %-4s %-8s n = %d: L1 %.1e, max bisections per cell %d%s%s\n",
                                analytic ? "analytic" : "Bernstein", c.name.c_str(), mesh.c_str(), part_text.c_str(), n,
                                known ? l1 : 0.0, r.max_bisections, problems.empty() ? "" : " FAILED:", message.c_str());
                    failures += !problems.empty();
                }
        }
    }
    return failures == 0 ? 0 : 1;
}
