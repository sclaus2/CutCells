// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Two sheets of one level set in a cell: the spherical shells of the
// robustness cases (support/robustness_cases.h, products of degree 4) on the
// cells of the octant x, y, z >= 0 of the n = 8 meshes, q = 3, as Bernstein
// coefficients and as analytic level sets. Every rule must pass the checks of
// the robustness cases. The shell 1e-3 wide on hexahedra and tetrahedra: boxes
// must be certified with two roots per line, and the per-cell L1 must stay
// below the bounds given; on hexahedra with Bernstein coefficients it must be
// ten times smaller than with one root per line only (Options::two_roots_depth
// past max_depth), and with analytic level sets not larger. The shell 1e-6
// wide takes the checks only, on hexahedra.
// Exits non-zero on failure.

#include <cutcells/selection_expr.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "support/robustness_cases.h"

using namespace cutcells;
using namespace cutcells::quadrays;
using namespace cutcells::quadrays::support;

int main()
{
    const int n = 8, first = 4, q = 3;
    struct Bound
    {
        const char* mesh;
        bool analytic;
        const char* part;
        double l1;
    };
    // about three times the measured L1
    const Bound bounds[] = {
        {"hex", false, "phi < 0", 3e-5}, {"hex", false, "phi = 0", 2.5e-3}, {"tet", false, "phi < 0", 1.5e-5},
        {"tet", false, "phi = 0", 1e-3}, {"hex", true, "phi < 0", 3e-2},    {"hex", true, "phi = 0", 6e-2},
        {"tet", true, "phi < 0", 1.5e-5}, {"tet", true, "phi = 0", 1e-3},
    };
    Case shell, thinner;
    for (const Case& c : all_cases())
    {
        if (c.name == "shell-1e-3")
            shell = c;
        if (c.name == "shell-1e-6")
            thinner = c;
    }
    Options one_root;
    one_root.two_roots_depth = one_root.max_depth + 1;

    int failures = 0;
    auto checks = [](const CaseRun& r, std::vector<std::string>& problems)
    {
        if (r.fail + r.negative + r.outside + r.side > 0)
            problems.push_back("fail/neg/out/side " + std::to_string(r.fail) + "/" + std::to_string(r.negative) + "/"
                               + std::to_string(r.outside) + "/" + std::to_string(r.side));
    };
    auto report = [&](const char* what, const std::vector<std::string>& problems)
    {
        std::string message;
        for (const std::string& p : problems)
            message += " " + p + ";";
        std::printf("%s%s%s\n", what, problems.empty() ? "" : " FAILED:", message.c_str());
        failures += !problems.empty();
    };
    for (const Bound& b : bounds)
    {
        SelectionExpr expr = parse_selection_expr(b.part);
        compile_selection_expr(expr, {"phi"});
        const SelectionTerm& term = expr.terms.front();
        const CaseRun r = run_case(shell, b.mesh, n, term, q, Options{}, b.analytic, first);
        const double l1 = r.l1 / r.exact_total;
        std::vector<std::string> problems;
        checks(r, problems);
        if (r.two_roots == 0)
            problems.push_back("no box certified with two roots");
        if (l1 > b.l1)
            problems.push_back("per-cell L1 above the bound");
        char what[256];
        int len = std::snprintf(what, sizeof what, "%-9s shell-1e-3 %s %-8s L1 %.1e (bound %.1e), two-root boxes %d",
                                b.analytic ? "analytic" : "Bernstein", b.mesh, b.part, l1, b.l1, r.two_roots);
        if (std::string(b.mesh) == "hex")
        {
            const CaseRun r1 = run_case(shell, b.mesh, n, term, q, one_root, b.analytic, first);
            const double l1_one = r1.l1 / r1.exact_total;
            checks(r1, problems);
            if (b.analytic ? l1 > l1_one : 10 * l1 > l1_one)
                problems.push_back(b.analytic ? "larger than with one root per line"
                                              : "not ten times smaller than with one root per line");
            std::snprintf(what + len, sizeof what - len, ", with one root per line %.1e", l1_one);
        }
        report(what, problems);
    }
    for (const bool analytic : {false, true})
        for (const char* part : {"phi < 0", "phi = 0"})
        {
            SelectionExpr expr = parse_selection_expr(part);
            compile_selection_expr(expr, {"phi"});
            const CaseRun r = run_case(thinner, "hex", n, expr.terms.front(), q, Options{}, analytic, first);
            std::vector<std::string> problems;
            checks(r, problems);
            char what[256];
            std::snprintf(what, sizeof what, "%-9s shell-1e-6 hex %-8s L1 %.1e (checks only)",
                          analytic ? "analytic" : "Bernstein", part, r.l1 / r.exact_total);
            report(what, problems);
        }
    return failures == 0 ? 0 : 1;
}
