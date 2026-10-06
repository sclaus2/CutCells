// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "../cell_types.h"
#include "../mesh_view.h"

namespace cutcells::part
{

/// Facets, edges or vertices of a mesh, of one dimension: in a CutResult those
/// on which a level set vanishes. Each is keyed by its sorted mesh nodes and
/// stored with every cell that shares it; entities are sorted by key.
///
/// The chart of an entity is the order of its nodes in nodes: sorted, except
/// that a quadrilateral starts at its lowest node, then that node's lower and
/// higher neighbours, then the opposite node (Basix's vertex order of a
/// quadrilateral). Rules on an entity have their points in the reference cell
/// of its type with the vertices in that order.
template <std::integral I = int>
struct EntityTable
{
    int dim = 0;                                  ///< of the entities
    std::vector<cell::type> types;                ///< per entity: point, interval, triangle or quadrilateral
    std::vector<I> nodes;                         ///< per entity, its mesh nodes in chart order
    std::vector<std::int32_t> node_offsets = {0}; ///< CSR over nodes
    std::vector<std::int32_t> cell_offsets = {0}; ///< CSR over cells
    std::vector<I> cells;                         ///< the cells sharing each entity, ascending
    std::vector<std::int8_t> local_index;         ///< the entity's number in each of them (Basix numbering)
    std::vector<std::uint64_t> vanishing;         ///< per entity: bit l, level set l vanishes on it
    /// Per entity: the vanishing level sets that count there, those that change
    /// sign across it (some sharing cell lies on each side), or all of them on
    /// the mesh boundary.
    std::vector<std::uint64_t> zero;
    std::vector<std::uint64_t> negative; ///< per sharing cell: bit l, the cell lies where vanishing l < 0
    std::vector<std::uint64_t> positive; ///< per sharing cell: bit l, where vanishing l > 0
    std::vector<std::uint8_t> boundary;  ///< per entity: 1 on the mesh boundary
    std::vector<std::array<I, 4>> keys;  ///< per entity, its nodes sorted, padded with the largest I

    int n_entities() const { return static_cast<int>(types.size()); }

    std::span<const I> entity_nodes(int e) const
    {
        return std::span<const I>(nodes).subspan(static_cast<std::size_t>(node_offsets[e]),
                                                 static_cast<std::size_t>(node_offsets[e + 1] - node_offsets[e]));
    }

    std::span<const I> entity_cells(int e) const
    {
        return std::span<const I>(cells).subspan(static_cast<std::size_t>(cell_offsets[e]),
                                                 static_cast<std::size_t>(cell_offsets[e + 1] - cell_offsets[e]));
    }
};

/// The cells around some vertices (nodes) of a mesh, with the vertex nodes of
/// those cells: what entity tables need to find the cells sharing an entity.
template <std::integral I = int>
struct NodeStars
{
    std::vector<I> nodes;                         ///< sorted
    std::vector<std::int32_t> offsets = {0};      ///< CSR over cells, per node
    std::vector<I> cells;                         ///< the cells around each node, ascending
    std::vector<I> star_cells;                    ///< all of them, ascending, each once
    std::vector<std::int32_t> cell_offsets = {0}; ///< CSR over cell_nodes, per star cell
    std::vector<I> cell_nodes;                    ///< their vertex nodes, Basix order
    /// Per star cell, bit f: its facet f (Basix) holds one of the nodes given
    /// for boundaries (node_stars) and is on the mesh boundary: no other cell
    /// holds it.
    std::vector<std::uint16_t> boundary_facets;

    /// The cells around @p node (one of nodes).
    std::span<const I> cells_around(I node) const
    {
        const std::size_t k = static_cast<std::size_t>(std::lower_bound(nodes.begin(), nodes.end(), node) - nodes.begin());
        return std::span<const I>(cells).subspan(static_cast<std::size_t>(offsets[k]),
                                                 static_cast<std::size_t>(offsets[k + 1] - offsets[k]));
    }

    /// The position of @p cell (one of star_cells) in star_cells.
    std::size_t star(I cell) const
    {
        return static_cast<std::size_t>(std::lower_bound(star_cells.begin(), star_cells.end(), cell)
                                        - star_cells.begin());
    }

    /// The vertex nodes of @p cell (one of star_cells), in Basix order.
    std::span<const I> vertex_nodes(I cell) const
    {
        const std::size_t k = star(cell);
        return std::span<const I>(cell_nodes).subspan(static_cast<std::size_t>(cell_offsets[k]),
                                                      static_cast<std::size_t>(cell_offsets[k + 1] - cell_offsets[k]));
    }
};

/// @brief The mesh nodes of a cell's vertices, in Basix order.
/// @param scratch  node indices, for meshes whose rows are not contiguous
template <std::floating_point T, std::integral I>
void cell_vertex_nodes(const MeshView<T, I>& mesh, I cell_id, std::vector<I>& nodes, std::vector<I>& scratch);

/// @brief The key of an entity: its nodes, sorted, padded with the largest I.
template <std::integral I>
std::array<I, 4> entity_key(std::span<const I> nodes)
{
    std::array<I, 4> key;
    key.fill(std::numeric_limits<I>::max());
    for (std::size_t k = 0; k < nodes.size(); ++k)
        key[k] = nodes[k];
    for (std::size_t a = 1; a < nodes.size(); ++a)
        for (std::size_t b = a; b > 0 && key[b] < key[b - 1]; --b)
            std::swap(key[b], key[b - 1]);
    return key;
}

/// @brief The chart order of an entity's nodes, given in the Basix vertex
/// order of a cell's sub-entity: sorted; a quadrilateral's lowest node, its
/// lower and higher neighbours, the opposite node.
template <std::integral I>
std::array<I, 4> chart_order(cell::type type, std::span<const I> local_nodes);

/// @brief The entity of @p table with key @p key, or -1.
template <std::integral I>
std::int32_t find_entity(const EntityTable<I>& table, const std::array<I, 4>& key);

/// @brief The cells of @p mesh around the vertices @p nodes (sorted, unique)
/// and their vertex nodes, in one pass over the cells, and which of their
/// facets holding one of @p boundary_nodes (sorted, among nodes) lie on the
/// mesh boundary (all the cells holding such a facet are star cells).
template <std::floating_point T, std::integral I>
void node_stars(const MeshView<T, I>& mesh, std::vector<I> nodes, std::span<const I> boundary_nodes,
                NodeStars<I>& out);

/// @brief The table of the entities of dimension @p dim with the keys @p keys:
/// their types, charts, sharing cells (every cell holding all the entity's
/// nodes), local numbers and boundary flags. vanishing, zero, negative and
/// positive are zero.
/// @param keys   sorted and unique
/// @param stars  around every node of the keys; below the facets' dimension,
///               given for boundaries every node of the keys
template <std::floating_point T, std::integral I>
EntityTable<I> make_entity_table(const MeshView<T, I>& mesh, int dim, std::vector<std::array<I, 4>> keys,
                                 const NodeStars<I>& stars);

} // namespace cutcells::part
