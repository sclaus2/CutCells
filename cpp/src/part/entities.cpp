// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "entities.h"

#include <algorithm>
#include <array>
#include <bit>
#include <stdexcept>
#include <utility>

#include "../cell_topology.h"
#include "../reference_cell.h"

namespace cutcells::part
{

template <std::floating_point T, std::integral I>
void cell_vertex_nodes(const MeshView<T, I>& mesh, I cell_id, std::vector<I>& nodes, std::vector<I>& scratch)
{
    const cell::type type = mesh.cell_type(cell_id);
    const int nv = cell::get_num_vertices(type);
    const std::span<const I> row = mesh.cell_nodes(cell_id, scratch);
    nodes.resize(static_cast<std::size_t>(nv));
    for (int v = 0; v < nv; ++v)
        nodes[static_cast<std::size_t>(v)]
            = row[static_cast<std::size_t>(mesh.vtk_vertex_order ? cell::basix_to_vtk_vertex(type, v) : v)];
}

template <std::integral I>
std::array<I, 4> chart_order(cell::type type, std::span<const I> local_nodes)
{
    if (type != cell::type::quadrilateral)
        return entity_key(local_nodes);
    // Basix's quadrilateral (0, 0), (1, 0), (0, 1), (1, 1): vertex p's
    // neighbours and the vertex opposite it
    constexpr int neighbours[4][2] = {{1, 2}, {0, 3}, {0, 3}, {1, 2}};
    const int p = static_cast<int>(std::min_element(local_nodes.begin(), local_nodes.end()) - local_nodes.begin());
    const I a = local_nodes[static_cast<std::size_t>(neighbours[p][0])],
            b = local_nodes[static_cast<std::size_t>(neighbours[p][1])];
    return {local_nodes[static_cast<std::size_t>(p)], std::min(a, b), std::max(a, b),
            local_nodes[static_cast<std::size_t>(3 - p)]};
}

template <std::integral I>
std::int32_t find_entity(const EntityTable<I>& table, const std::array<I, 4>& key)
{
    const auto it = std::lower_bound(table.keys.begin(), table.keys.end(), key);
    if (it == table.keys.end() || *it != key)
        return -1;
    return static_cast<std::int32_t>(it - table.keys.begin());
}

template <std::floating_point T, std::integral I>
void node_stars(const MeshView<T, I>& mesh, std::vector<I> nodes, std::span<const I> boundary_nodes,
                NodeStars<I>& out)
{
    out = NodeStars<I>{};
    out.nodes = std::move(nodes);
    std::vector<char> marked(static_cast<std::size_t>(mesh.num_nodes()), 0);
    for (const I n : out.nodes)
        marked[static_cast<std::size_t>(n)] = 1;
    for (const I n : boundary_nodes)
        marked[static_cast<std::size_t>(n)] = 2;
    std::vector<std::pair<I, I>> pairs; // (node, cell)
    std::vector<I> vertices, scratch;
    for (I c = 0; c < mesh.num_cells(); ++c)
    {
        cell_vertex_nodes(mesh, c, vertices, scratch);
        bool star = false;
        for (const I n : vertices)
            if (marked[static_cast<std::size_t>(n)])
            {
                pairs.emplace_back(n, c);
                star = true;
            }
        if (!star)
            continue;
        out.star_cells.push_back(c);
        out.cell_nodes.insert(out.cell_nodes.end(), vertices.begin(), vertices.end());
        out.cell_offsets.push_back(static_cast<std::int32_t>(out.cell_nodes.size()));
    }
    std::sort(pairs.begin(), pairs.end());
    std::size_t k = 0;
    for (const I n : out.nodes)
    {
        for (; k < pairs.size() && pairs[k].first == n; ++k)
            out.cells.push_back(pairs[k].second);
        out.offsets.push_back(static_cast<std::int32_t>(out.cells.size()));
    }

    // the facets holding one of the boundary nodes, by key: those of one cell
    // are on the boundary
    std::vector<std::pair<std::array<I, 4>, std::pair<std::int32_t, int>>> facets; // key, (star, facet)
    std::vector<I> local;
    for (std::size_t s = 0; s < out.star_cells.size(); ++s)
    {
        const cell::type type = mesh.cell_type(out.star_cells[s]);
        const int tdim = cell::get_tdim(type);
        for (int f = 0; f < cell::num_sub_entities(type, tdim - 1); ++f)
        {
            local.clear();
            bool held = false;
            for (const int v : cell::sub_entity_vertices(type, tdim - 1, f))
            {
                local.push_back(out.cell_nodes[static_cast<std::size_t>(out.cell_offsets[s] + v)]);
                held |= marked[static_cast<std::size_t>(local.back())] == 2;
            }
            if (held)
                facets.push_back({entity_key(std::span<const I>(local)), {static_cast<std::int32_t>(s), f}});
        }
    }
    std::sort(facets.begin(), facets.end());
    out.boundary_facets.assign(out.star_cells.size(), 0);
    for (std::size_t i = 0; i < facets.size();)
    {
        std::size_t j = i;
        while (j < facets.size() && facets[j].first == facets[i].first)
            ++j;
        if (j - i == 1)
            out.boundary_facets[static_cast<std::size_t>(facets[i].second.first)]
                |= static_cast<std::uint16_t>(1u << facets[i].second.second);
        i = j;
    }
}

namespace
{
/// The vertices (bits) of sub-entity @p i of dimension @p dim of a cell, made
/// once per cell type.
unsigned vertex_mask(cell::type type, int dim, int i)
{
    static const std::array<std::array<std::vector<unsigned>, 3>, 8> masks = []
    {
        std::array<std::array<std::vector<unsigned>, 3>, 8> m;
        for (const cell::type t : {cell::type::interval, cell::type::triangle, cell::type::quadrilateral,
                                   cell::type::tetrahedron, cell::type::hexahedron, cell::type::prism,
                                   cell::type::pyramid})
            for (int d = 0; d < cell::get_tdim(t); ++d)
                for (int e = 0; e < cell::num_sub_entities(t, d); ++e)
                {
                    unsigned mask = 0;
                    for (const int v : cell::sub_entity_vertices(t, d, e))
                        mask |= 1u << v;
                    m[static_cast<std::size_t>(t)][static_cast<std::size_t>(d)].push_back(mask);
                }
        return m;
    }();
    return masks[static_cast<std::size_t>(type)][static_cast<std::size_t>(dim)][static_cast<std::size_t>(i)];
}
} // namespace

template <std::floating_point T, std::integral I>
EntityTable<I> make_entity_table(const MeshView<T, I>& mesh, int dim, std::vector<std::array<I, 4>> keys,
                                 const NodeStars<I>& stars)
{
    EntityTable<I> table;
    table.dim = dim;
    const I pad = std::numeric_limits<I>::max();
    table.types.reserve(keys.size());
    table.nodes.reserve(keys.size() * static_cast<std::size_t>(dim + 1));
    table.node_offsets.reserve(keys.size() + 1);
    table.cell_offsets.reserve(keys.size() + 1);
    table.boundary.reserve(keys.size());
    std::vector<I> local;
    for (const std::array<I, 4>& key : keys)
    {
        const int size = static_cast<int>(std::count_if(key.begin(), key.end(), [pad](I n) { return n != pad; }));
        const std::size_t first = table.cells.size();
        cell::type type = cell::type::point;
        std::array<I, 4> chart{};
        bool boundary = false;
        int tdim = 0;
        for (const I c : stars.cells_around(key[0]))
        {
            // the cell's vertices at the entity's nodes: a sub-entity if the cell holds them all
            const std::size_t star = stars.star(c);
            const std::span<const I> nodes
                = std::span<const I>(stars.cell_nodes)
                      .subspan(static_cast<std::size_t>(stars.cell_offsets[star]),
                               static_cast<std::size_t>(stars.cell_offsets[star + 1] - stars.cell_offsets[star]));
            unsigned mask = 0;
            for (std::size_t v = 0; v < nodes.size(); ++v)
                if (std::find(key.begin(), key.begin() + size, nodes[v]) != key.begin() + size)
                    mask |= 1u << v;
            if (std::popcount(mask) != size)
                continue;
            const cell::type ctype = mesh.cell_type(c);
            tdim = cell::get_tdim(ctype);
            int index = -1;
            for (int i = 0; i < cell::num_sub_entities(ctype, dim) && index < 0; ++i)
                if (vertex_mask(ctype, dim, i) == mask)
                    index = i;
            if (index < 0)
                continue; // the nodes are not one sub-entity of the cell
            if (table.cells.size() == first)
            {
                type = cell::sub_entity_type(ctype, dim, index);
                local.clear();
                for (const int v : cell::sub_entity_vertices(ctype, dim, index))
                    local.push_back(nodes[static_cast<std::size_t>(v)]);
                chart = chart_order(type, std::span<const I>(local));
            }
            table.cells.push_back(c);
            table.local_index.push_back(static_cast<std::int8_t>(index));
            // on the mesh boundary if one of the cell's boundary facets holds it
            for (unsigned bits = stars.boundary_facets[star]; bits != 0 && !boundary; bits &= bits - 1)
                boundary = (vertex_mask(ctype, tdim - 1, std::countr_zero(bits)) & mask) == mask;
        }
        if (table.cells.size() == first)
            throw std::runtime_error("part: no cell holds an entity of the mesh");
        if (dim == tdim - 1)
            boundary = table.cells.size() - first == 1;
        table.types.push_back(type);
        table.nodes.insert(table.nodes.end(), chart.begin(), chart.begin() + size);
        table.node_offsets.push_back(static_cast<std::int32_t>(table.nodes.size()));
        table.cell_offsets.push_back(static_cast<std::int32_t>(table.cells.size()));
        table.boundary.push_back(boundary ? 1 : 0);
    }
    table.keys = std::move(keys);
    table.vanishing.assign(table.keys.size(), 0);
    table.zero.assign(table.keys.size(), 0);
    table.negative.assign(table.cells.size(), 0);
    table.positive.assign(table.cells.size(), 0);
    return table;
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template void cell_vertex_nodes<float, int>(const MeshView<float, int>&, int, std::vector<int>&, std::vector<int>&);
template void cell_vertex_nodes<double, int>(const MeshView<double, int>&, int, std::vector<int>&,
                                             std::vector<int>&);
template std::array<int, 4> chart_order<int>(cell::type, std::span<const int>);
template std::int32_t find_entity<int>(const EntityTable<int>&, const std::array<int, 4>&);
template void node_stars<float, int>(const MeshView<float, int>&, std::vector<int>, std::span<const int>,
                                     NodeStars<int>&);
template void node_stars<double, int>(const MeshView<double, int>&, std::vector<int>, std::span<const int>,
                                      NodeStars<int>&);
template EntityTable<int> make_entity_table<float, int>(const MeshView<float, int>&, int,
                                                        std::vector<std::array<int, 4>>, const NodeStars<int>&);
template EntityTable<int> make_entity_table<double, int>(const MeshView<double, int>&, int,
                                                         std::vector<std::array<int, 4>>, const NodeStars<int>&);

} // namespace cutcells::part
