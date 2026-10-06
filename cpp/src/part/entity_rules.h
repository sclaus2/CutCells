// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "../cell_types.h"
#include "mesh_part.h"

namespace cutcells::part
{

/// How the cells that share a facet, edge or vertex take a rule on it, in the
/// cell rules of quadrature_rules.
enum class EntityPolicy : std::uint8_t
{
    split,    ///< every sharing cell, with weight 1/n for n sharing cells
    negative, ///< the cells on the negative side of the zero level sets of the piece, 1/k each; all if none is
    positive, ///< the cells on their positive side, 1/k each; all if none is
    owner,    ///< the lowest sharing cell, with weight 1
    exclude   ///< none: cells keep to their interior, entity_rules has the rules
};

/// @brief "split", "negative", "positive", "owner" or "exclude".
std::string entity_policy_to_str(EntityPolicy policy);

/// @throws std::invalid_argument for other strings
EntityPolicy string_to_entity_policy(const std::string& policy);

/// The rules of a part on the entities of one dimension (facets, edges or
/// vertices), one rule per entity: its points in the entity's reference cell
/// (the vertices in the order of the entity's chart, EntityTable), weights
/// physical, and the cells that share the entity with their local numbers,
/// their sides and the points in their reference coordinates.
///
/// Sharing cell s of rule r is cells[s] for s in [cell_offsets[r],
/// cell_offsets[r + 1]); its points are cell_points[cell_point_offsets[s] * tdim,
/// cell_point_offsets[s + 1] * tdim), in the order of the rule's points.
template <std::floating_point T>
struct EntityRules
{
    int tdim = 0;                                       ///< of the background cells
    int entity_dim = 0;                                 ///< tdim - 1 for facets, ..., 0 for vertices
    std::vector<cell::type> entity_types;               ///< per rule: point, interval, triangle, quadrilateral
    std::vector<std::int32_t> entity_nodes;             ///< per rule, the entity's mesh nodes in chart order
    std::vector<std::int32_t> node_offsets = {0};       ///< CSR over entity_nodes
    std::vector<std::uint64_t> zero;                    ///< per rule: the level sets whose zero sets hold its points
    std::vector<T> points;                              ///< entity reference coordinates, entity_dim each
    std::vector<T> weights;                             ///< physical (1 for a point)
    std::vector<std::int32_t> offsets = {0};            ///< per rule, CSR over points
    std::vector<std::int32_t> cell_offsets = {0};       ///< per rule, CSR over its sharing cells
    std::vector<std::int32_t> cells;                    ///< the sharing cells, ascending per rule
    std::vector<std::int8_t> local_index;               ///< the entity in each of them (Basix numbering)
    std::vector<std::uint64_t> negative;                ///< per sharing cell: bit i, level set i < 0 there
    std::vector<std::uint64_t> positive;                ///< per sharing cell: bit i, level set i > 0 there
    std::vector<T> cell_points;                         ///< per sharing cell, the rule's points in its
                                                        ///< reference coordinates (tdim each)
    std::vector<std::int32_t> cell_point_offsets = {0}; ///< per sharing cell, CSR over cell_points

    int n_rules() const { return static_cast<int>(entity_types.size()); }
};

/// @brief The rules of @p part on the entities of dimension part.dim that it
/// holds whole (MeshPart::entities): each entity's reference rule of degree
/// @p order (the lowest at least as high with positive weights), mapped to
/// physical space for the weights (a Jacobian per point on bilinear
/// quadrilaterals), and to each sharing cell's reference coordinates; a vertex
/// has one point of weight 1.
template <std::floating_point T, std::integral I>
EntityRules<T> whole_entity_rules(const MeshPart<T, I>& part, int order);

/// @brief Appends rule @p r's points in sharing cell @p s (an index into
/// rules.cells) to @p points and its weights times @p factor to @p weights.
template <std::floating_point T>
void append_cell_contribution(const EntityRules<T>& rules, int r, std::int32_t s, T factor, std::vector<T>& points,
                              std::vector<T>& weights);

/// @brief The sharing cells of rule @p r that take it under @p policy, as
/// indices into rules.cells, with the factor of the weights each takes.
template <std::floating_point T>
void policy_cells(const EntityRules<T>& rules, int r, EntityPolicy policy, std::vector<std::int32_t>& slots,
                  T& factor);

/// @brief Appends the sharing cells of the last rule of @p rules, whose points
/// are in: those of entity @p e of @p table, their local numbers and sides,
/// and the rule's points in each of them, mapped from the entity's chart.
template <std::floating_point T, std::integral I>
void append_sharing_cells(const MeshView<T, I>& mesh, const EntityTable<I>& table, std::int32_t e,
                          EntityRules<T>& rules);

} // namespace cutcells::part
