// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Test meshes of [-1, 1]^3 (hexahedra, Kuhn tetrahedra, prisms, pyramids) and
// [-1, 1]^2 (quadrilaterals, triangles), and the Bernstein coefficients of
// polynomial level sets on their cells. Shared by the tests in
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

/// The cells of the grid cube [lo, lo + h]^3: the cube itself ("hex"), its
/// six Kuhn tetrahedra ("tet"), two prisms ("prism": the triangles of the
/// bottom face extruded along z) or six pyramids ("pyramid": the faces with the
/// cube's centre as apex). Faces are given relative to @p origin.
inline std::vector<TestCell> grid_cells(const std::string& mesh, const V3& lo, double h, const V3& origin)
{
    auto rel = [&](const V3& x) { return exact::V3{x[0] - origin[0], x[1] - origin[1], x[2] - origin[2]}; };
    std::vector<TestCell> cells;
    // a cell from its vertices (Basix order) and the loops of its faces
    auto polytope = [&](cell::type type, const std::vector<V3>& vertices, const std::vector<std::vector<V3>>& loops,
                        double volume)
    {
        TestCell c;
        c.type = type;
        std::vector<std::vector<exact::V3>> rel_loops;
        for (const auto& loop : loops)
        {
            rel_loops.emplace_back();
            for (const V3& x : loop)
                rel_loops.back().push_back(rel(x));
        }
        c.faces = exact::make_faces(rel_loops);
        for (const V3& x : vertices)
        {
            c.vertices.insert(c.vertices.end(), x.begin(), x.end());
            for (int d = 0; d < 3; ++d)
                c.centroid[d] += x[d] / double(vertices.size());
        }
        for (const V3& x : vertices)
        {
            double r2 = 0;
            for (int d = 0; d < 3; ++d)
                r2 += (x[d] - c.centroid[d]) * (x[d] - c.centroid[d]);
            c.radius = std::max(c.radius, std::sqrt(r2));
        }
        c.volume = volume;
        cells.push_back(c);
    };
    auto P = [&](int i, int j, int k) { return V3{lo[0] + i * h, lo[1] + j * h, lo[2] + k * h}; };
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
    if (mesh == "prism")
    {
        for (const auto& b : {std::array{P(0, 0, 0), P(1, 0, 0), P(0, 1, 0)}, std::array{P(1, 1, 0), P(0, 1, 0), P(1, 0, 0)}})
        {
            std::array<V3, 3> t = b;
            for (V3& x : t)
                x[2] += h;
            polytope(cell::type::prism, {b[0], b[1], b[2], t[0], t[1], t[2]},
                     {{b[0], b[1], b[2]},
                      {t[0], t[1], t[2]},
                      {b[0], b[1], t[1], t[0]},
                      {b[1], b[2], t[2], t[1]},
                      {b[2], b[0], t[0], t[2]}},
                     0.5 * h * h * h);
        }
        return cells;
    }
    if (mesh == "pyramid")
    {
        // the cube's faces as bases v0, v1, v2, v3 = v1 + v2 - v0 (Basix order)
        const V3 apex = {lo[0] + 0.5 * h, lo[1] + 0.5 * h, lo[2] + 0.5 * h};
        const std::array<std::array<V3, 4>, 6> bases = {{
            {P(0, 0, 0), P(1, 0, 0), P(0, 1, 0), P(1, 1, 0)},
            {P(0, 0, 1), P(1, 0, 1), P(0, 1, 1), P(1, 1, 1)},
            {P(0, 0, 0), P(1, 0, 0), P(0, 0, 1), P(1, 0, 1)},
            {P(0, 1, 0), P(1, 1, 0), P(0, 1, 1), P(1, 1, 1)},
            {P(0, 0, 0), P(0, 1, 0), P(0, 0, 1), P(0, 1, 1)},
            {P(1, 0, 0), P(1, 1, 0), P(1, 0, 1), P(1, 1, 1)},
        }};
        for (const auto& b : bases)
            polytope(cell::type::pyramid, {b[0], b[1], b[2], b[3], apex},
                     {{b[0], b[1], b[3], b[2]}, {b[0], b[1], apex}, {b[1], b[3], apex}, {b[3], b[2], apex}, {b[2], b[0], apex}},
                     h * h * h / 6.0);
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

/// The cells of the grid square [lo, lo + h]^2 (z ignored): the square
/// ("quad") or its two triangles ("tri"), vertices in Basix order with 2
/// coordinates each; loop is the polygon in the plane z = 0.
struct TestCell2D
{
    cell::type type = cell::type::quadrilateral;
    std::vector<double> vertices; ///< Basix order, 2 per vertex
    std::vector<exact::V3> loop;  ///< the polygon, z = 0
    double area = 0;
};

inline std::vector<TestCell2D> grid_cells_2d(const std::string& mesh, const V3& lo, double h)
{
    std::vector<TestCell2D> cells;
    auto P = [&](int i, int j) { return exact::V3{lo[0] + i * h, lo[1] + j * h, 0}; };
    if (mesh == "quad")
    {
        TestCell2D c;
        c.type = cell::type::quadrilateral;
        for (const auto& [i, j] : {std::pair{0, 0}, std::pair{1, 0}, std::pair{0, 1}, std::pair{1, 1}})
            c.vertices.insert(c.vertices.end(), {lo[0] + i * h, lo[1] + j * h});
        c.loop = {P(0, 0), P(1, 0), P(1, 1), P(0, 1)};
        c.area = h * h;
        cells.push_back(c);
        return cells;
    }
    if (mesh != "tri")
        throw std::runtime_error("unknown 2D mesh: " + mesh);
    for (const auto& tri : {std::array{P(0, 0), P(1, 0), P(1, 1)}, std::array{P(0, 0), P(1, 1), P(0, 1)}})
    {
        TestCell2D c;
        c.type = cell::type::triangle;
        for (const exact::V3& v : tri)
            c.vertices.insert(c.vertices.end(), {v[0], v[1]});
        c.loop = {tri[0], tri[1], tri[2]};
        c.area = 0.5 * h * h;
        cells.push_back(c);
    }
    return cells;
}

/// Physical point (z = 0) of the reference coordinates @p xi of a 2D test cell (affine).
inline V3 physical(const TestCell2D& cell, const std::array<double, 2>& xi)
{
    // vertices 1 and 2 are next to vertex 0 along xi0 and xi1 on both cell types
    V3 x = {cell.vertices[0], cell.vertices[1], 0};
    for (int d = 0; d < 2; ++d)
        x[d] += xi[0] * (cell.vertices[2 + d] - cell.vertices[d]) + xi[1] * (cell.vertices[4 + d] - cell.vertices[d]);
    return x;
}

/// True if the reference point @p xi (tdim coordinates) lies in the reference
/// cell of @p type, up to @p tol.
inline bool in_reference_cell(cell::type type, std::span<const double> xi, double tol)
{
    for (double x : xi)
        if (x < -tol)
            return false;
    switch (type)
    {
    case cell::type::triangle:
        return xi[0] + xi[1] <= 1 + tol;
    case cell::type::tetrahedron:
        return xi[0] + xi[1] + xi[2] <= 1 + tol;
    case cell::type::prism:
        return xi[0] + xi[1] <= 1 + tol && xi[2] <= 1 + tol;
    case cell::type::pyramid:
        return xi[0] + xi[2] <= 1 + tol && xi[1] + xi[2] <= 1 + tol;
    default:
        for (double x : xi)
            if (x > 1 + tol)
                return false;
        return true;
    }
}

/// Physical point of the reference coordinates @p xi of a test cell (affine).
inline V3 physical(const TestCell& cell, const V3& xi)
{
    // the vertices next to vertex 0 along the reference axes
    const bool hex = cell.type == cell::type::hexahedron || cell.type == cell::type::pyramid;
    const std::array<int, 3> axes = hex ? std::array<int, 3>{1, 2, 4} : std::array<int, 3>{1, 2, 3};
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
    // the equispaced lattice of the cell's Lagrange space
    auto in_lattice = [&](int i, int j, int k)
    {
        switch (cell.type)
        {
        case cell::type::tetrahedron:
            return i + j + k <= degree;
        case cell::type::prism:
            return i + j <= degree;
        case cell::type::pyramid:
            return i + k <= degree && j + k <= degree;
        default:
            return true;
        }
    };
    for (int k = 0; k <= degree; ++k)
        for (int j = 0; j <= degree; ++j)
            for (int i = 0; i <= degree; ++i)
            {
                if (!in_lattice(i, j, k))
                    continue;
                const V3 xi = {double(i) / degree, double(j) / degree, double(k) / degree};
                points.insert(points.end(), xi.begin(), xi.end());
                values.push_back(phi(physical(cell, xi)));
            }
    bernstein::lagrange_to_bernstein<double>(cell.type, degree, points, values, coeffs);
}

/// cell_coefficients on a 2D test cell; @p phi takes points with z = 0.
inline void cell_coefficients(const TestCell2D& cell, int degree, const std::function<double(const V3&)>& phi,
                              std::vector<double>& coeffs)
{
    std::vector<double> points, values;
    const bool quad = cell.type == cell::type::quadrilateral;
    for (int j = 0; j <= degree; ++j)
        for (int i = 0; i <= degree; ++i)
        {
            if (!quad && i + j > degree)
                continue;
            const std::array<double, 2> xi = {double(i) / degree, double(j) / degree};
            points.insert(points.end(), xi.begin(), xi.end());
            values.push_back(phi(physical(cell, xi)));
        }
    bernstein::lagrange_to_bernstein<double>(cell.type, degree, points, values, coeffs);
}

} // namespace cutcells::quadrays::support
