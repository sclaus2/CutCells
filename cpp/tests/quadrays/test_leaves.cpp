// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Leaf cells of the decomposition: VTK's Lagrange node order (bijective, corners
// first in VTK's order), and the leaves of a sphere on hexahedra and Kuhn
// tetrahedra: complete, hexahedra with a positive Jacobian, interface
// quadrilaterals with their normal along grad phi, interface nodes on the
// sphere and volume nodes inside the ball. VTK itself checks the node order in
// python/tests/test_quadrays_leaves.py. Exits non-zero on failure.

#include <cutcells/quadrays/leaves.h>
#include <cutcells/selection_expr.h>
#include <cutcells/write_vtk.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "support/test_mesh.h"

using namespace cutcells;
using namespace cutcells::quadrays;
using namespace cutcells::quadrays::support;

namespace
{
int failures = 0;

void check(bool ok, const std::string& what)
{
    if (!ok)
    {
        std::printf("FAILED: %s\n", what.c_str());
        ++failures;
    }
}

void test_node_order()
{
    // VTK's corner order of a hexahedron: bottom face counter-clockwise, then top
    const int corners[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    for (int p = 1; p <= 6; ++p)
    {
        const int p1 = p + 1;
        std::vector<int> seen(p1 * p1 * p1, 0);
        for (int k = 0; k <= p; ++k)
            for (int j = 0; j <= p; ++j)
                for (int i = 0; i <= p; ++i)
                {
                    const int index = io::vtk_lagrange_hexahedron_index(i, j, k, p);
                    if (index >= 0 && index < p1 * p1 * p1)
                        ++seen[index];
                }
        bool bijective = true;
        for (int s : seen)
            bijective &= s == 1;
        check(bijective, "hexahedron node order is a bijection, degree " + std::to_string(p));
        for (int c = 0; c < 8; ++c)
            check(io::vtk_lagrange_hexahedron_index(p * corners[c][0], p * corners[c][1], p * corners[c][2], p) == c,
                  "hexahedron corner " + std::to_string(c));

        std::vector<int> seen2(p1 * p1, 0);
        for (int j = 0; j <= p; ++j)
            for (int i = 0; i <= p; ++i)
            {
                const int index = io::vtk_lagrange_quadrilateral_index(i, j, p);
                if (index >= 0 && index < p1 * p1)
                    ++seen2[index];
            }
        bijective = true;
        for (int s : seen2)
            bijective &= s == 1;
        check(bijective, "quadrilateral node order is a bijection, degree " + std::to_string(p));
        for (int c = 0; c < 4; ++c)
            check(io::vtk_lagrange_quadrilateral_index(p * corners[c][0], p * corners[c][1], p) == c,
                  "quadrilateral corner " + std::to_string(c));
    }
}

V3 sub(const V3& a, const V3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
V3 cross(const V3& a, const V3& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
double dot(const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

void test_sphere_leaves(const char* mesh, int degree)
{
    const int n = 4;
    const double h = 2.0 / n, r = 0.7;
    const V3 c = {0.0123, -0.0371, 0.0217};
    auto phi = [&](const V3& x) { return dot(sub(x, c), sub(x, c)) - r * r; };
    std::vector<double> coeffs;
    for (const char* part_text : {"phi < 0", "phi = 0"})
    {
        SelectionExpr expr = parse_selection_expr(part_text);
        compile_selection_expr(expr, {"phi"});
        const Part part = part_of(expr.terms.front());
        LeafMesh<double> leaves;
        Stats stats;
        std::int32_t index = 0;
        for (int i0 = 0; i0 < n; ++i0)
            for (int i1 = 0; i1 < n; ++i1)
                for (int i2 = 0; i2 < n; ++i2)
                    for (const TestCell& cell : grid_cells(mesh, {-1 + h * i0, -1 + h * i1, -1 + h * i2}, h, c))
                    {
                        cell_coefficients(cell, 2, phi, coeffs);
                        ClippedBox<double> box;
                        BoxBernstein<double> form;
                        make_clipped_box<double>(cell.type, cell.vertices, 3, box);
                        cell_bernstein_on_box<double>(cell.type, 2, coeffs, form);
                        append_leaves(box, form, part, degree, Options{}, index++, leaves, stats);
                    }
        const std::string what = std::string(mesh) + " " + part_text + ", degree " + std::to_string(degree);
        check(leaves.n_cells() > 0 && stats.incomplete_leaves == 0, what + ": leaves complete");
        int wrong_side = 0, wrong_orientation = 0;
        double off_sphere = 0;
        for (int cell = 0; cell < leaves.n_cells(); ++cell)
        {
            const std::int32_t* nodes = leaves.connectivity.data() + leaves.offsets[cell];
            auto x = [&](int local) -> V3
            {
                const double* p = leaves.points.data() + 3 * nodes[local];
                return {p[0], p[1], p[2]};
            };
            const int count = leaves.offsets[cell + 1] - leaves.offsets[cell];
            for (int a = 0; a < count; ++a)
            {
                const double f = phi(x(a));
                if (part == Part::interface)
                    off_sphere = std::max(off_sphere, std::abs(std::sqrt(f + r * r) - r));
                else
                    wrong_side += f > 1e-12;
            }
            if (part == Part::interface)
            {
                // corners 0, 1, 3 in VTK order span the quadrilateral: normal along grad phi
                const V3 normal = cross(sub(x(1), x(0)), sub(x(3), x(0)));
                const V3 centre = {0.25 * (x(0)[0] + x(1)[0] + x(2)[0] + x(3)[0]),
                                   0.25 * (x(0)[1] + x(1)[1] + x(2)[1] + x(3)[1]),
                                   0.25 * (x(0)[2] + x(1)[2] + x(2)[2] + x(3)[2])};
                wrong_orientation += dot(normal, sub(centre, c)) < 0;
            }
            else
            {
                // corners 0, 1, 3, 4 in VTK order span the hexahedron
                wrong_orientation += dot(cross(sub(x(1), x(0)), sub(x(3), x(0))), sub(x(4), x(0))) < 0;
            }
        }
        check(wrong_side == 0, what + ": volume nodes inside the ball (" + std::to_string(wrong_side) + " outside)");
        check(off_sphere < 1e-14, what + ": interface nodes on the sphere");
        check(wrong_orientation == 0, what + ": orientation (" + std::to_string(wrong_orientation) + " flipped)");
    }
}
} // namespace

int main()
{
    test_node_order();
    for (int degree : {1, 2, 3})
    {
        test_sphere_leaves("hex", degree);
        test_sphere_leaves("tet", degree);
    }
    if (failures == 0)
        std::printf("test_leaves: ok\n");
    return failures == 0 ? 0 : 1;
}
