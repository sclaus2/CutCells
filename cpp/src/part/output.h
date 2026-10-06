// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <string>
#include <vector>

#include "../lut/cell_pieces.h"
#include "../lut/cut_mesh.h"
#include "../quadrature.h"
#include "../quadrays/engine.h"
#include "../quadrays/leaves.h"
#include "entity_rules.h"
#include "mesh_part.h"

namespace cutcells::part
{

/// @brief Quadrature rules of a part, one rule per cell, ordered by cell: the
/// backend's rules on the cut cells, with @p include_uncut_cells the rules of
/// the whole cells, and the rules on the facets, edges and vertices the part
/// holds (entity_rules), taken by the cells that share them as @p entities
/// says. The type of the options picks the backend: quadrays (the default) or
/// the lookup tables.
///
/// Points are in the cells' reference coordinates, weights are physical.
/// quadrays takes every level set that cuts a cell and that a term holding on
/// a piece of the cell names; an interface part is integrated once per level
/// set whose zero set it asks for, a curve part once per pair. The points
/// where two curves cross on 2D cells weigh 1.
///
/// Cells keep to their interior: a piece on a facet, edge or vertex is
/// integrated once, on that entity, and comes into the cells through
/// @p entities. With split (the default) every sharing cell takes the entity's
/// rule with weight 1/n, so each piece counts once in total.
///
/// @param order  the polynomial degree integrated exactly on flat pieces, 1 to
///               10: quadrays takes ceil((order + part.dim) / 2) Gauss-Legendre
///               points per segment of each height line (each nested level's
///               linear bounds raise the degree by one), or
///               options.points_per_segment; whole cells and entities get the
///               reference rules of degree order
/// @throws std::invalid_argument (quadrays) for parts of dimension below
///         tdim - 2, and for cells quadrays does not take
template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> quadrature_rules(const MeshPart<T, I>& part, int order, bool include_uncut_cells,
                                                const quadrays::Options& options = {},
                                                EntityPolicy entities = EntityPolicy::split);

/// @brief Quadrature rules of a part from the lookup tables: the selected
/// straight pieces of each cut cell (lut::cut_cell, with the level sets that
/// cut the cell and that the expression names), with the reference rules of
/// degree @p order (1 to 10) on the pieces and whole cells, and the rules on
/// entities as @p entities says. A piece lying in a facet or an edge of its
/// cell on which its zero sets all vanish is the entity's, and dropped.
/// @throws std::invalid_argument for cells other than intervals, triangles,
///         quadrilaterals, tetrahedra and hexahedra, and for parts where
///         three level sets vanish
template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> quadrature_rules(const MeshPart<T, I>& part, int order, bool include_uncut_cells,
                                                const lut::Options& options,
                                                EntityPolicy entities = EntityPolicy::split);

/// @brief The rules of a part on facets, edges and vertices, one EntityRules
/// per entity dimension (index d for dimension d < tdim): the entities it holds
/// whole (MeshPart::entities), and on 2D cells the points where two curves
/// cross that the backend's cells find on an edge or a vertex (within
/// 2e-10 in reference coordinates), one point per position, weight 1.
template <std::floating_point T, std::integral I>
std::vector<EntityRules<T>> entity_rules(const MeshPart<T, I>& part, int order, const quadrays::Options& options = {});

/// @brief entity_rules with the lookup tables finding the points on 2D cells.
template <std::floating_point T, std::integral I>
std::vector<EntityRules<T>> entity_rules(const MeshPart<T, I>& part, int order, const lut::Options& options);

/// @brief A part as cells for visualisation (quadrays): on cut cells the
/// leaves of the quadrays decomposition as Lagrange cells of @p degree, the
/// entities the part holds as linear cells, once each (their lowest sharing
/// cell as parent), and with @p include_uncut_cells the whole cells as linear
/// cells.
template <std::floating_point T, std::integral I>
quadrays::LeafMesh<T> visualization_mesh(const MeshPart<T, I>& part, int degree, bool include_uncut_cells,
                                         const quadrays::Options& options = {});

/// @brief A part as straight cells for visualisation (lookup tables): the
/// selected pieces of cut cells, the entities the part holds, once each, and
/// with @p include_uncut_cells the whole cells; _parent_map holds each cell's
/// background cell.
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
