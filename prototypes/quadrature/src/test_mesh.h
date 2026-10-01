// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "clipped_box.h"
#include "exact_reference.h"

namespace cutcells::proto
{

/// A cell of the test mesh: the clipped box handed to the generators and its
/// physical polytope (relative to a given origin) for the exact references.
struct TestCell
{
    ClippedBox box;
    std::vector<exact::Face> faces;
    Vec3 centroid;
    double radius = 0; ///< circumradius about the centroid
    double volume = 0;
};

/// The cells of the grid cube [lo, lo + h]^3: the cube itself ("hex") or its six Kuhn
/// tetrahedra ("tet"). Faces are given relative to @p origin.
inline std::vector<TestCell> grid_cells(const std::string& mesh, const Vec3& lo, double h, const Vec3& origin)
{
    auto rel = [&](const Vec3& x) { return exact::V3{x[0] - origin[0], x[1] - origin[1], x[2] - origin[2]}; };
    std::vector<TestCell> cells;
    if (mesh == "hex")
    {
        TestCell c;
        c.box = hex_cell(lo, h);
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
        // Kuhn tetrahedron: path from lo along the axes p[0], p[1], p[2]
        std::array<Vec3, 4> X;
        X[0] = lo;
        for (int k = 0; k < 3; ++k)
        {
            X[k + 1] = X[k];
            X[k + 1][p[k]] += h;
        }
        TestCell c;
        c.box = tet_cell(X);
        c.faces = exact::tet_faces({rel(X[0]), rel(X[1]), rel(X[2]), rel(X[3])});
        c.centroid = {0, 0, 0};
        for (const Vec3& x : X)
            for (int d = 0; d < 3; ++d)
                c.centroid[d] += 0.25 * x[d];
        for (const Vec3& x : X)
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

} // namespace cutcells::proto
