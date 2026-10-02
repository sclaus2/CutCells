// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <string>

#include "../cut_mesh.h"
#include "../lut/cell_pieces.h"
#include "../quadrature.h"
#include "../quadrays/engine.h"
#include "../quadrays/leaves.h"
#include "mesh_part.h"

namespace cutcells::part
{

/// @brief Quadrature rules of a part, one rule per cell, ordered by cell: the
/// backend's rules on the cut cells, rules on the zero faces the cells own,
/// and with @p include_uncut_cells the rules of the whole cells. The type of
/// the options picks the backend: quadrays (the default) or the lookup tables.
///
/// Points are in the cells' reference coordinates, weights are physical.
///
/// @param order  quadrays: Gauss-Legendre points per segment of each height
///               line; the lookup tables' straight pieces, whole cells and zero
///               faces get the reference rules exact for degree 2 order - 1
///               (at most 10)
/// @throws std::runtime_error (quadrays) where a cell would need two level
///         sets at once (several level sets per cell come in phase 6)
template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> quadrature_rules(const MeshPart<T, I>& part, int order, bool include_uncut_cells,
                                                const quadrays::Options& options = {});

/// @brief Quadrature rules of a part from the lookup tables: the selected
/// straight pieces of each cut cell (lut::cut_cell, with the level sets that
/// cut the cell and that the expression names).
/// @throws std::invalid_argument for cells other than triangles,
///         quadrilaterals, tetrahedra and hexahedra, and for parts where two
///         level sets vanish
template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> quadrature_rules(const MeshPart<T, I>& part, int order, bool include_uncut_cells,
                                                const lut::Options& options);

/// @brief A part as cells for visualisation (quadrays): on cut cells the
/// leaves of the quadrays decomposition as Lagrange cells of @p degree, zero
/// faces as linear faces, and with @p include_uncut_cells the whole cells as
/// linear cells.
template <std::floating_point T, std::integral I>
quadrays::LeafMesh<T> visualization_mesh(const MeshPart<T, I>& part, int degree, bool include_uncut_cells,
                                         const quadrays::Options& options = {});

/// @brief A part as straight cells for visualisation (lookup tables): the
/// selected pieces of cut cells, zero faces, and with @p include_uncut_cells
/// the whole cells; _parent_map holds each cell's background cell.
template <std::floating_point T, std::integral I>
mesh::CutMesh<T> visualization_mesh(const MeshPart<T, I>& part, bool include_uncut_cells,
                                    const lut::Options& options);

/// @brief Write visualization_mesh (quadrays) to a .vtu file.
template <std::floating_point T, std::integral I>
void write_vtu(const std::string& filename, const MeshPart<T, I>& part, int degree, bool include_uncut_cells,
               const quadrays::Options& options = {});

/// @brief Write visualization_mesh (lookup tables) to a .vtu file.
template <std::floating_point T, std::integral I>
void write_vtu(const std::string& filename, const MeshPart<T, I>& part, bool include_uncut_cells,
               const lut::Options& options);

} // namespace cutcells::part
