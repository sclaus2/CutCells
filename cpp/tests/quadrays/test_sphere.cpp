// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// The sphere of radius 0.7 off the grid's symmetry against exact per-cell values,
// on hexahedra and Kuhn tetrahedra of [-1, 1]^3 with n = 8. Per-cell L1 (relative
// to the ball volume or the sphere area) and the worst cell (relative to the
// cell's exact value) must stay within 25% of the prototype's numbers
// (docs/quadrays/RESULTS.md), and so must the bisections.
// Exits non-zero on failure.

#include <cutcells/quadrays/rules.h>
#include <cutcells/selection_expr.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "support/exact_reference.h"
#include "support/test_mesh.h"

using namespace cutcells;
using namespace cutcells::quadrays;
using namespace cutcells::quadrays::support;

namespace
{
struct Expected
{
    const char* mesh;
    int q;
    const char* part;
    double l1, worst;
    long bisections;
};

// The prototype (branch quadrature-1d, commit 8b4d248) at n = 8, margin 0.25.
const Expected expected[] = {
    {"tet", 3, "phi < 0", 3.4e-7, 5.2e-5, 3156}, {"tet", 3, "phi = 0", 7.5e-6, 3.0e-4, 3156},
    {"tet", 5, "phi < 0", 1.1e-9, 3.9e-7, 3156}, {"tet", 5, "phi = 0", 5.5e-8, 3.9e-6, 3156},
    {"hex", 3, "phi < 0", 6.5e-6, 2.1e-4, 0},    {"hex", 3, "phi = 0", 4.8e-5, 7.2e-4, 0},
    {"hex", 5, "phi < 0", 2.2e-8, 9.4e-7, 0},    {"hex", 5, "phi = 0", 3.7e-7, 8.6e-6, 0},
};
} // namespace

int main()
{
    const int n = 8;
    const double h = 2.0 / n, r = 0.7;
    const V3 c = {0.0123, -0.0371, 0.0217};
    const double ball = 4.0 / 3.0 * M_PI * r * r * r, sphere = 4.0 * M_PI * r * r;
    auto phi = [&](const V3& x)
    { return (x[0] - c[0]) * (x[0] - c[0]) + (x[1] - c[1]) * (x[1] - c[1]) + (x[2] - c[2]) * (x[2] - c[2]) - r * r; };

    int failures = 0;
    std::vector<double> coeffs;
    for (const Expected& e : expected)
    {
        SelectionExpr expr = parse_selection_expr(e.part);
        compile_selection_expr(expr, {"phi"});
        const SelectionTerm& term = expr.terms.front();
        const bool surface = part_of(term) == Part::interface;
        double sum_abs = 0, worst = 0;
        Stats stats;
        for (int i0 = 0; i0 < n; ++i0)
            for (int i1 = 0; i1 < n; ++i1)
                for (int i2 = 0; i2 < n; ++i2)
                {
                    const V3 lo = {-1 + h * i0, -1 + h * i1, -1 + h * i2};
                    for (const TestCell& cell : grid_cells(e.mesh, lo, h, c))
                    {
                        double dc = 0;
                        for (int d = 0; d < 3; ++d)
                            dc += (cell.centroid[d] - c[d]) * (cell.centroid[d] - c[d]);
                        dc = std::sqrt(dc);
                        if (dc - cell.radius >= r || dc + cell.radius <= r)
                            continue;
                        const double area = exact::sphere_area(cell.faces, r);
                        if (area <= 0.0)
                            continue;
                        const double exact_value = surface ? area : exact::ball_volume(cell.faces, r, area);
                        cell_coefficients(cell, 2, phi, coeffs);
                        quadrature::QuadratureRules<double> rules;
                        append_cell_rules<double>(cell.type, cell.vertices, 2, coeffs, term, 0, e.q, Options{}, 0,
                                                  rules, stats);
                        double value = 0;
                        for (double w : rules._weights)
                            value += w;
                        sum_abs += std::abs(value - exact_value);
                        const double floor = surface ? 1e-3 * h * h : 1e-3 * h * h * h;
                        if (exact_value > floor)
                            worst = std::max(worst, std::abs(value - exact_value) / exact_value);
                    }
                }
        const double l1 = sum_abs / (surface ? sphere : ball);
        const bool ok = l1 <= 1.25 * e.l1 && worst <= 1.25 * e.worst && stats.bisections <= 1.25 * e.bisections;
        std::printf("%s q = %d %-8s L1 %.1e (prototype %.1e), worst %.1e (%.1e), bisections %d (%ld) %s\n", e.mesh, e.q,
                    e.part, l1, e.l1, worst, e.worst, stats.bisections, e.bisections, ok ? "" : "FAILED");
        failures += !ok;
    }
    return failures == 0 ? 0 : 1;
}
