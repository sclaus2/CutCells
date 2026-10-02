// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Test meshes of [-1, 1]^3 (hexahedra and Kuhn tetrahedra) and the Bernstein
// coefficients of polynomial level sets on their cells. Shared by the tests in
// cpp/tests/quadrays and the drivers in benchmarks/quadrays.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <cutcells/bernstein.h>
#include <cutcells/cell_types.h>

#include "exact_reference.h"

namespace cutcells::quadrays::support
{

using V3 = std::array<double, 3>;

/// A cell of a test mesh: type and vertices for the engine, and its polytope
/// for the exact references.
struct TestCell
{
    cell::type type = cell::type::hexahedron;
    std::vector<double> vertices;   ///< Basix order, 3 per vertex
    std::vector<exact::Face> faces; ///< relative to the origin given to grid_cells
    V3 centroid = {0, 0, 0};
    double radius = 0; ///< circumradius about the centroid
    double volume = 0;
};

/// The cells of the grid cube [lo, lo + h]^3: the cube itself ("hex") or its
/// six Kuhn tetrahedra ("tet"). Faces are given relative to @p origin.
inline std::vector<TestCell> grid_cells(const std::string& mesh, const V3& lo, double h, const V3& origin)
{
    auto rel = [&](const V3& x) { return exact::V3{x[0] - origin[0], x[1] - origin[1], x[2] - origin[2]}; };
    std::vector<TestCell> cells;
    if (mesh == "hex")
    {
        TestCell c;
        c.type = cell::type::hexahedron;
        for (int v = 0; v < 8; ++v) // Basix order: vertex v at (v & 1, (v >> 1) & 1, v >> 2)
            for (int d = 0; d < 3; ++d)
                c.vertices.push_back(lo[d] + h * ((v >> d) & 1));
        c.faces = exact::box_faces(rel(lo), h);
        c.centroid = {lo[0] + 0.5 * h, lo[1] + 0.5 * h, lo[2] + 0.5 * h};
        c.radius = 0.5 * std::sqrt(3.0) * h;
        c.volume = h * h * h;
        cells.push_back(c);
        return cells;
    }
    if (mesh != "tet")
        throw std::runtime_error("unknown mesh: " + mesh);
    std::array<int, 3> p = {0, 1, 2};
    do
    {
        // Kuhn tetrahedron: the path from lo along the axes p[0], p[1], p[2]
        std::array<V3, 4> X;
        X[0] = lo;
        for (int k = 0; k < 3; ++k)
        {
            X[k + 1] = X[k];
            X[k + 1][p[k]] += h;
        }
        TestCell c;
        c.type = cell::type::tetrahedron;
        for (const V3& x : X)
            c.vertices.insert(c.vertices.end(), x.begin(), x.end());
        c.faces = exact::tet_faces({rel(X[0]), rel(X[1]), rel(X[2]), rel(X[3])});
        for (const V3& x : X)
            for (int d = 0; d < 3; ++d)
                c.centroid[d] += 0.25 * x[d];
        for (const V3& x : X)
        {
            double r2 = 0;
            for (int d = 0; d < 3; ++d)
                r2 += (x[d] - c.centroid[d]) * (x[d] - c.centroid[d]);
            c.radius = std::max(c.radius, std::sqrt(r2));
        }
        c.volume = h * h * h / 6.0;
        cells.push_back(c);
    } while (std::next_permutation(p.begin(), p.end()));
    return cells;
}

/// Physical point of the reference coordinates @p xi of a test cell (affine).
inline V3 physical(const TestCell& cell, const V3& xi)
{
    const std::array<int, 3> axes = cell.type == cell::type::hexahedron ? std::array<int, 3>{1, 2, 4}
                                                                         : std::array<int, 3>{1, 2, 3};
    V3 x;
    for (int d = 0; d < 3; ++d)
    {
        x[d] = cell.vertices[d];
        for (int k = 0; k < 3; ++k)
            x[d] += xi[k] * (cell.vertices[3 * axes[k] + d] - cell.vertices[d]);
    }
    return x;
}

/// Bernstein coefficients, in CutCells' order, of the polynomial @p phi of the
/// given degree on a cell: from its values on the equispaced lattice, as
/// make_cell_level_set computes them for a finite-element level set.
inline void cell_coefficients(const TestCell& cell, int degree, const std::function<double(const V3&)>& phi,
                              std::vector<double>& coeffs)
{
    if (degree == 0)
    {
        coeffs.assign(1, phi(cell.centroid));
        return;
    }
    std::vector<double> points, values;
    const bool hex = cell.type == cell::type::hexahedron;
    for (int k = 0; k <= degree; ++k)
        for (int j = 0; j <= degree; ++j)
            for (int i = 0; i <= degree; ++i)
            {
                if (!hex && i + j + k > degree)
                    continue;
                const V3 xi = {double(i) / degree, double(j) / degree, double(k) / degree};
                points.insert(points.end(), xi.begin(), xi.end());
                values.push_back(phi(physical(cell, xi)));
            }
    bernstein::lagrange_to_bernstein<double>(cell.type, degree, points, values, coeffs);
}

} // namespace cutcells::quadrays::support
