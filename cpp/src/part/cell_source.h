// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <span>
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
    std::vector<T> vertices; ///< Basix order, gdim per vertex
    std::vector<I> nodes;    ///< scratch
};

/// One cell and several level sets as quadrays reads them, in the order given.
template <std::floating_point T, std::integral I = int>
struct CellSources
{
    quadrays::ClippedBox<T> box;
    std::vector<quadrays::BoxBernstein<T>> forms; ///< per level set; unused for analytic ones
    std::vector<quadrays::Source<T>> sources;
    LevelSetCell<T, I> ls_cell;
    std::vector<T> vertices; ///< Basix order, gdim per vertex
    std::vector<I> nodes;    ///< scratch
};

/// @brief Box coordinates of vertex @p v (Basix numbering) of a cell that
/// quadrays takes; they equal the cell's reference coordinates (u2 = 0 on 2D
/// cells).
template <std::floating_point T>
quadrays::Vec3<T> box_vertex(cell::type type, int v)
{
    switch (type)
    {
    case cell::type::triangle:
    case cell::type::tetrahedron:
    {
        quadrays::Vec3<T> u = {0, 0, 0};
        if (v > 0)
            u[v - 1] = T(1);
        return u;
    }
    case cell::type::prism:
        return {T(v % 3 == 1), T(v % 3 == 2), T(v >= 3)};
    case cell::type::pyramid:
        return v == 4 ? quadrays::Vec3<T>{0, 0, 1} : quadrays::Vec3<T>{T(v & 1), T((v >> 1) & 1), 0};
    default: // quadrilaterals and hexahedra
        return {T(v & 1), T((v >> 1) & 1), T((v >> 2) & 1)};
    }
}

/// @brief True if quadrays takes level sets on cells of type @p type in a mesh
/// of geometric dimension @p gdim: triangles and quadrilaterals in 2D,
/// tetrahedra, hexahedra, prisms and pyramids in 3D.
template <std::floating_point T, std::integral I>
bool quadrays_takes(cell::type type, int gdim, const LevelSetFunction<T, I>&)
{
    return quadrays::supported_cell(type) && gdim == cell::get_tdim(type);
}

/// @brief Fill @p out for level set @p ls on cell @p cell_id.
/// @return false if quadrays does not take the cell or the level set on it
///         (quadrays_takes)
/// @throws std::invalid_argument if @p ls has neither dof values nor an
///         analytic level set
template <std::floating_point T, std::integral I>
bool cell_source(const MeshView<T, I>& mesh, const LevelSetFunction<T, I>& ls, I cell_id,
                 CellSource<T, I>& out);

/// @brief Fill @p out for the level sets @p level_sets on cell @p cell_id.
/// @return false if quadrays does not take the cell or one of the level sets
///         on it
template <std::floating_point T, std::integral I>
bool cell_sources(const MeshView<T, I>& mesh, std::span<const LevelSetFunction<T, I>* const> level_sets,
                  I cell_id, CellSources<T, I>& out);

} // namespace cutcells::part
