// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <vector>

#include "../level_set.h"
#include "../level_set_cell.h"
#include "../mesh_view.h"
#include "../quadrays/box_bernstein.h"
#include "../quadrays/clipped_box.h"
#include "../quadrays/source.h"

namespace cutcells::part
{

/// One cell and one level set as quadrays reads them: the cell's clipped box
/// and the level set on it, a Bernstein form on the box or the analytic level
/// set. Reused from cell to cell.
template <std::floating_point T, std::integral I = int>
struct CellSource
{
    quadrays::ClippedBox<T> box;
    quadrays::BoxBernstein<T> form;
    quadrays::Source<T> source;
    LevelSetCell<T, I> ls_cell;
    std::vector<T> vertices; ///< Basix order, 3 per vertex
    std::vector<I> nodes;    ///< scratch
};

/// @brief Box coordinates of vertex @p v (Basix numbering) of a tetrahedron's
/// or hexahedron's box; they equal the cell's reference coordinates.
template <std::floating_point T>
quadrays::Vec3<T> box_vertex(cell::type type, int v)
{
    if (type == cell::type::tetrahedron)
    {
        quadrays::Vec3<T> u = {0, 0, 0};
        if (v > 0)
            u[v - 1] = T(1);
        return u;
    }
    return {T(v & 1), T((v >> 1) & 1), T((v >> 2) & 1)};
}

/// @brief Fill @p out for level set @p ls on cell @p cell_id.
/// @return false if quadrays does not take the cell (only tetrahedra and
///         hexahedra in 3D)
/// @throws std::invalid_argument if @p ls has neither dof values nor an
///         analytic level set
template <std::floating_point T, std::integral I>
bool cell_source(const MeshView<T, I>& mesh, const LevelSetFunction<T, I>& ls, I cell_id,
                 CellSource<T, I>& out);

} // namespace cutcells::part
