// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Plane exactness: quadrays integrates a planar level set exactly. Every cell of
// hexahedral and Kuhn-tetrahedral meshes of [-1, 1]^3, cut or not, for three
// parts and planes in general position, given as degree-1 and degree-2 level
// sets; per cell against exact polytope cuts and in total, to 1e-13.
// Exits non-zero on failure.

#include <cutcells/quadrays/rules.h>
#include <cutcells/selection_expr.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "support/exact_reference.h"
#include "support/test_mesh.h"

using namespace cutcells;
using namespace cutcells::quadrays;
using namespace cutcells::quadrays::support;

int main()
{
    struct Plane
    {
        V3 a;
        double b;
    };
    const std::vector<Plane> planes = {{{1.0, 0.3, -0.2}, 0.0}, {{-0.45, 1.0, 0.7}, 0.1}};
    struct Mesh
    {
        const char* name;
        int n;
    };
    const std::vector<Mesh> meshes = {{"hex", 5}, {"tet", 4}, {"tet", 7}};
    const std::vector<std::string> parts = {"phi < 0", "phi > 0", "phi = 0"};

    int failures = 0;
    std::vector<double> coeffs;
    double worst_overall = 0;
    for (const Plane& plane : planes)
        for (int degree : {1, 2})
            for (const Mesh& mesh : meshes)
                for (int q : {3, 5})
                    for (const std::string& part_text : parts)
                    {
                        SelectionExpr expr = parse_selection_expr(part_text);
                        compile_selection_expr(expr, {"phi"});
                        const SelectionTerm& term = expr.terms.front();
                        const Part part = part_of(term);
                        auto phi = [&](const V3& x)
                        { return plane.a[0] * x[0] + plane.a[1] * x[1] + plane.a[2] * x[2] - plane.b; };

                        const double h = 2.0 / mesh.n;
                        const double cell_scale = part == Part::interface ? h * h : h * h * h;
                        double total = 0, exact_total = 0, worst = 0;
                        for (int i0 = 0; i0 < mesh.n; ++i0)
                            for (int i1 = 0; i1 < mesh.n; ++i1)
                                for (int i2 = 0; i2 < mesh.n; ++i2)
                                {
                                    const V3 lo = {-1 + h * i0, -1 + h * i1, -1 + h * i2};
                                    for (const TestCell& cell : grid_cells(mesh.name, lo, h, {0, 0, 0}))
                                    {
                                        const exact::PlaneCut cut = exact::plane_cut(cell.faces, plane.a, plane.b);
                                        const double exact_value = part == Part::negative ? cut.volume_below
                                                                   : part == Part::positive
                                                                       ? cell.volume - cut.volume_below
                                                                       : cut.cut_area;
                                        cell_coefficients(cell, degree, phi, coeffs);
                                        quadrature::QuadratureRules<double> rules;
                                        Stats stats;
                                        append_cell_rules<double>(cell.type, cell.vertices, degree, coeffs, term, 0, q,
                                                                  Options{}, 0, rules, stats);
                                        double value = 0;
                                        for (double w : rules._weights)
                                            value += w;
                                        total += value;
                                        exact_total += exact_value;
                                        worst = std::max(worst, std::abs(value - exact_value) / cell_scale);
                                    }
                                }
                        const double rel_total = std::abs(total - exact_total) / exact_total;
                        worst_overall = std::max({worst_overall, worst, rel_total});
                        if (worst > 1e-13 || rel_total > 1e-13)
                        {
                            std::printf("FAILED: plane (%g, %g, %g) . x = %g, degree %d, %s n = %d, q = %d, %s: "
                                        "worst cell %.1e, total %.1e\n",
                                        plane.a[0], plane.a[1], plane.a[2], plane.b, degree, mesh.name, mesh.n, q,
                                        part_text.c_str(), worst, rel_total);
                            ++failures;
                        }
                    }
    std::printf("test_plane_exactness: largest error %.1e (per cell relative to h^3 or h^2, totals relative)\n",
                worst_overall);
    return failures == 0 ? 0 : 1;
}
