// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// [-1, 1]^3 as a MeshView of n^3 hexahedra or 6 n^3 Kuhn tetrahedra, with its
// storage, cell by cell in the order of quadrays' test cells (grid_cells), so
// that cell c of the mesh is the c-th test cell with its exact references.

#pragma once

#include <cmath>
#include <string>
#include <vector>

#include <cutcells/mesh_view.h>

#include "../../quadrays/support/test_mesh.h"

namespace cutcells::part::support
{

using quadrays::support::TestCell;
using quadrays::support::V3;

/// A mesh and the storage its view points to: build in place, do not copy.
struct BoxMesh
{
    std::vector<double> coordinates;
    std::vector<int> connectivity, offsets = {0};
    std::vector<cell::type> types;
    std::vector<TestCell> cells; ///< relative to origin
    MeshView<double, int> view;
};

/// @brief Fill @p mesh; the test cells' faces are relative to @p origin.
inline void make_box_mesh(const std::string& kind, int n, const V3& origin, BoxMesh& mesh)
{
    const double h = 2.0 / n;
    mesh.coordinates.clear();
    for (int k = 0; k <= n; ++k)
        for (int j = 0; j <= n; ++j)
            for (int i = 0; i <= n; ++i)
                mesh.coordinates.insert(mesh.coordinates.end(), {-1 + h * i, -1 + h * j, -1 + h * k});
    auto node = [&](const double* x)
    {
        int idx[3];
        for (int d = 0; d < 3; ++d)
            idx[d] = static_cast<int>(std::lround((x[d] + 1) / h));
        return idx[0] + (n + 1) * (idx[1] + (n + 1) * idx[2]);
    };
    for (int i0 = 0; i0 < n; ++i0)
        for (int i1 = 0; i1 < n; ++i1)
            for (int i2 = 0; i2 < n; ++i2)
                for (TestCell& cell : quadrays::support::grid_cells(kind, {-1 + h * i0, -1 + h * i1, -1 + h * i2}, h,
                                                                    origin))
                {
                    for (std::size_t v = 0; v < cell.vertices.size(); v += 3)
                        mesh.connectivity.push_back(node(cell.vertices.data() + v));
                    mesh.offsets.push_back(static_cast<int>(mesh.connectivity.size()));
                    mesh.types.push_back(cell.type);
                    mesh.cells.push_back(std::move(cell));
                }
    mesh.view.gdim = 3;
    mesh.view.tdim = 3;
    mesh.view.coordinates = mesh.coordinates;
    mesh.view.connectivity = mesh.connectivity;
    mesh.view.offsets = mesh.offsets;
    mesh.view.cell_types = mesh.types;
}

} // namespace cutcells::part::support
