// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// A zero set that passes by a pyramid's apex. The sphere of radius 0.7 about
// c = h (1/2, 0.309, 0.207) on the grid of [-1, 1]^3 with n = 16 (the sweep of
// benchmarks/quadrays at s = 1/2) passes 6e-4 from the centres of two grid
// cubes, the apexes of their six pyramids. A pyramid's form (1 - z)^n phi has
// a degenerate zero at the apex, so that no direction certifies the boxes
// there, which are integrated uncertified at the depth limit. As c is centred
// in x, the level set does not vary along the first box axis of four of the
// pyramids: the height direction must not be that axis, whose lines find no
// roots (the whole interface of a pyramid was lost at depths 6, 9, 12, ...).
//  - The P1 interpolant, affine on these pyramids, as degree-1 and degree-2
//    forms: the area of its zero set exact to 1e-13 h^2 for every depth limit
//    from 4 to 16.
//  - The P2 interpolant, the sphere itself: its area against the exact
//    reference at q = 8, to 1.5e-10 h^2.
// Exits non-zero on failure.

#include <cutcells/quadrays/rules.h>
#include <cutcells/selection_expr.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include "support/exact_reference.h"
#include "support/test_mesh.h"

using namespace cutcells;
using namespace cutcells::quadrays;
using namespace cutcells::quadrays::support;

namespace
{
/// Area of the zero set of the level set with these coefficients in the cell.
double interface_area(const TestCell& cell, int degree, const std::vector<double>& coeffs, int q, const Options& opt)
{
    SelectionExpr expr = parse_selection_expr("phi = 0");
    compile_selection_expr(expr, {"phi"});
    quadrature::QuadratureRules<double> rule;
    Stats stats;
    append_cell_rules<double>(cell.type, cell.vertices, degree, coeffs, expr.terms.front(), 0, q, opt, 0, rule,
                              stats);
    double sum = 0;
    for (const double w : rule._weights)
        sum += w;
    return sum;
}
} // namespace

int main()
{
    const double h = 0.125, r = 0.7;
    const V3 c = {0.5 * h, 0.25 * h * (std::sqrt(5.0) - 1.0), 0.5 * h * (std::sqrt(2.0) - 1.0)};
    const auto sphere = [&](const V3& x)
    { return (x[0] - c[0]) * (x[0] - c[0]) + (x[1] - c[1]) * (x[1] - c[1]) + (x[2] - c[2]) * (x[2] - c[2]) - r * r; };
    // the two grid cubes whose centres lie 6e-4 inside the sphere
    const std::vector<V3> cubes = {{0, 4 * h, -4 * h}, {0, -2 * h, 5 * h}};

    int failures = 0;
    double worst_p1 = 0, worst_p2 = 0;
    std::vector<double> coeffs;
    for (const V3& lo : cubes)
    {
        const std::vector<TestCell> cells = grid_cells("pyramid", lo, h, {0, 0, 0});
        const std::vector<TestCell> centred = grid_cells("pyramid", lo, h, c); // faces relative to c
        for (std::size_t i = 0; i < cells.size(); ++i)
        {
            const TestCell& cell = cells[i];
            auto vertex = [&](int v)
            { return V3{cell.vertices[3 * v], cell.vertices[3 * v + 1], cell.vertices[3 * v + 2]}; };

            // P1: the affine a . x - b through vertices 0, 1, 2 and the apex 4,
            // which takes the sphere's value at vertex 3 too
            const V3 x0 = vertex(0);
            const std::array<V3, 3> e = {exact::sub(vertex(1), x0), exact::sub(vertex(2), x0), exact::sub(vertex(4), x0)};
            const std::array<V3, 3> dual = {exact::cross(e[1], e[2]), exact::cross(e[2], e[0]), exact::cross(e[0], e[1])};
            const double det = exact::dot(e[0], dual[0]);
            V3 a = {0, 0, 0};
            for (int j = 0; j < 3; ++j)
                a = exact::add(a, exact::mul((sphere(exact::add(x0, e[j])) - sphere(x0)) / det, dual[j]));
            const double b = exact::dot(a, x0) - sphere(x0);
            const auto plane = [&](const V3& x) { return exact::dot(a, x) - b; };
            if (std::abs(plane(vertex(3)) - sphere(vertex(3))) > 1e-14)
            {
                std::printf("FAILED: the P1 interpolant is not affine on pyramid %zu of the cube at (%g, %g, %g)\n", i,
                            lo[0], lo[1], lo[2]);
                ++failures;
            }
            const double exact_p1 = exact::plane_cut(cell.faces, a, b).cut_area;
            for (const int degree : {1, 2})
            {
                cell_coefficients(cell, degree, plane, coeffs);
                for (int max_depth = 4; max_depth <= 16; ++max_depth)
                {
                    Options opt;
                    opt.max_depth = max_depth;
                    const double err = std::abs(interface_area(cell, degree, coeffs, 3, opt) - exact_p1) / (h * h);
                    worst_p1 = std::max(worst_p1, err);
                    if (err > 1e-13)
                    {
                        std::printf("FAILED: P1 as degree %d, pyramid %zu of the cube at (%g, %g, %g), max_depth %d: "
                                    "error %.1e h^2\n",
                                    degree, i, lo[0], lo[1], lo[2], max_depth, err);
                        ++failures;
                    }
                }
            }

            // P2: the sphere
            cell_coefficients(cell, 2, sphere, coeffs);
            const double exact_p2 = exact::sphere_area(centred[i].faces, r);
            const double err = std::abs(interface_area(cell, 2, coeffs, 8, Options{}) - exact_p2) / (h * h);
            worst_p2 = std::max(worst_p2, err);
            if (err > 1.5e-10)
            {
                std::printf("FAILED: P2, pyramid %zu of the cube at (%g, %g, %g): error %.1e h^2\n", i, lo[0], lo[1],
                            lo[2], err);
                ++failures;
            }
        }
    }
    std::printf("pyramid apex: largest error %.1e h^2 (P1, all depth limits), %.1e h^2 (P2, q = 8)\n", worst_p1,
                worst_p2);
    return failures == 0 ? 0 : 1;
}
