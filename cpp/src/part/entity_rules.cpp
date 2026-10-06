// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "entity_rules.h"

#include <algorithm>
#include <stdexcept>

#include "../cell_topology.h"
#include "../lut/piece_rules.h"
#include "../reference_cell.h"

namespace cutcells::part
{

std::string entity_policy_to_str(EntityPolicy policy)
{
    switch (policy)
    {
    case EntityPolicy::split:
        return "split";
    case EntityPolicy::negative:
        return "negative";
    case EntityPolicy::positive:
        return "positive";
    case EntityPolicy::owner:
        return "owner";
    case EntityPolicy::exclude:
        return "exclude";
    }
    throw std::invalid_argument("part: unknown entity policy");
}

EntityPolicy string_to_entity_policy(const std::string& policy)
{
    for (const EntityPolicy p : {EntityPolicy::split, EntityPolicy::negative, EntityPolicy::positive,
                                 EntityPolicy::owner, EntityPolicy::exclude})
        if (policy == entity_policy_to_str(p))
            return p;
    throw std::invalid_argument("part: entities takes 'split', 'negative', 'positive', 'owner' or 'exclude'; got '"
                                + policy + "'");
}

template <std::floating_point T, std::integral I>
void append_sharing_cells(const MeshView<T, I>& mesh, const EntityTable<I>& table, std::int32_t e,
                          EntityRules<T>& rules)
{
    const int r = rules.n_rules() - 1;
    const int m = rules.entity_dim, tdim = rules.tdim;
    const std::span<const I> chart = table.entity_nodes(e);
    const std::size_t first = static_cast<std::size_t>(rules.offsets[r]);
    const std::size_t np = static_cast<std::size_t>(rules.offsets[r + 1] - rules.offsets[r]);
    const std::span<const T> xi = std::span<const T>(rules.points).subspan(first * static_cast<std::size_t>(m),
                                                                           np * static_cast<std::size_t>(m));
    // the entity's chart into the cell: its vertices' reference coordinates there
    lut::CellMap<T> map;
    map.type = table.types[static_cast<std::size_t>(e)];
    map.gdim = tdim;
    map.vertices.resize(chart.size() * static_cast<std::size_t>(tdim));
    std::vector<I> nodes, scratch;
    std::vector<T> x;
    for (std::int32_t s = table.cell_offsets[e]; s < table.cell_offsets[e + 1]; ++s)
    {
        const I c = table.cells[static_cast<std::size_t>(s)];
        const cell::type type = mesh.cell_type(c);
        const std::vector<T> ref = cell::reference_vertices<T>(type);
        cell_vertex_nodes(mesh, c, nodes, scratch);
        const std::span<const int> sv
            = cell::sub_entity_vertices(type, table.dim, table.local_index[static_cast<std::size_t>(s)]);
        for (std::size_t j = 0; j < chart.size(); ++j)
        {
            const auto k = std::find_if(sv.begin(), sv.end(),
                                        [&](int v) { return nodes[static_cast<std::size_t>(v)] == chart[j]; });
            std::copy_n(ref.begin() + *k * tdim, tdim, map.vertices.begin() + static_cast<std::ptrdiff_t>(j * tdim));
        }
        if (m == 0)
        {
            x.clear();
            for (std::size_t q = 0; q < np; ++q)
                x.insert(x.end(), map.vertices.begin(), map.vertices.begin() + tdim);
        }
        else
            lut::push_forward(map, xi, x);
        rules.cells.push_back(static_cast<std::int32_t>(c));
        rules.local_index.push_back(table.local_index[static_cast<std::size_t>(s)]);
        rules.negative.push_back(table.negative[static_cast<std::size_t>(s)]);
        rules.positive.push_back(table.positive[static_cast<std::size_t>(s)]);
        rules.cell_points.insert(rules.cell_points.end(), x.begin(), x.end());
        rules.cell_point_offsets.push_back(rules.cell_point_offsets.back() + static_cast<std::int32_t>(np));
    }
    rules.cell_offsets.push_back(static_cast<std::int32_t>(rules.cells.size()));
}

template <std::floating_point T, std::integral I>
EntityRules<T> whole_entity_rules(const MeshPart<T, I>& part, int order)
{
    const CutResult<T, I>& r = *part.result;
    const MeshView<T, I>& mesh = *r.mesh;
    EntityRules<T> rules;
    rules.tdim = r.tdim();
    rules.entity_dim = part.dim;
    if (part.dim < 0 || part.dim >= rules.tdim)
        return rules;
    const std::size_t d = static_cast<std::size_t>(part.dim);
    const EntityTable<I>& table = r.entities[d];
    lut::CellMap<T> map;
    map.gdim = mesh.gdim;
    for (std::size_t k = 0; k < part.entities[d].size(); ++k)
    {
        const std::int32_t e = part.entities[d][k];
        const cell::type type = table.types[static_cast<std::size_t>(e)];
        const std::span<const I> chart = table.entity_nodes(e);
        rules.entity_types.push_back(type);
        for (const I n : chart)
            rules.entity_nodes.push_back(static_cast<std::int32_t>(n));
        rules.node_offsets.push_back(static_cast<std::int32_t>(rules.entity_nodes.size()));
        rules.zero.push_back(part.entity_zero[d][k]);
        if (part.dim == 0)
            rules.weights.push_back(T(1));
        else
        {
            // the reference rule through the entity's map: physical weights,
            // with a Jacobian per point on bilinear quadrilaterals
            map.type = type;
            map.vertices.clear();
            for (const I n : chart)
                map.vertices.insert(map.vertices.end(), mesh.node(n), mesh.node(n) + mesh.gdim);
            const std::vector<T> ref = cell::reference_vertices<T>(type);
            lut::append_piece_rule(map, type, std::span<const T>(ref), order, rules.points, rules.weights);
        }
        rules.offsets.push_back(static_cast<std::int32_t>(rules.weights.size()));
        append_sharing_cells(mesh, table, e, rules);
    }
    return rules;
}

template <std::floating_point T>
void append_cell_contribution(const EntityRules<T>& rules, int r, std::int32_t s, T factor, std::vector<T>& points,
                              std::vector<T>& weights)
{
    const std::size_t tdim = static_cast<std::size_t>(rules.tdim);
    points.insert(points.end(),
                  rules.cell_points.begin() + static_cast<std::ptrdiff_t>(rules.cell_point_offsets[s] * tdim),
                  rules.cell_points.begin() + static_cast<std::ptrdiff_t>(rules.cell_point_offsets[s + 1] * tdim));
    for (std::int32_t q = rules.offsets[r]; q < rules.offsets[r + 1]; ++q)
        weights.push_back(rules.weights[static_cast<std::size_t>(q)] * factor);
}

template <std::floating_point T>
void policy_cells(const EntityRules<T>& rules, int r, EntityPolicy policy, std::vector<std::int32_t>& slots,
                  T& factor)
{
    slots.clear();
    factor = T(0);
    const std::int32_t first = rules.cell_offsets[r], last = rules.cell_offsets[r + 1];
    if (policy == EntityPolicy::exclude || last == first)
        return;
    if (policy == EntityPolicy::owner)
    {
        // the cells are ascending: the lowest
        slots.push_back(first);
        factor = T(1);
        return;
    }
    if (policy == EntityPolicy::negative || policy == EntityPolicy::positive)
    {
        const std::uint64_t zero = rules.zero[static_cast<std::size_t>(r)];
        const std::vector<std::uint64_t>& side = policy == EntityPolicy::negative ? rules.negative : rules.positive;
        for (std::int32_t s = first; s < last && zero != 0; ++s)
            if ((side[static_cast<std::size_t>(s)] & zero) == zero)
                slots.push_back(s);
    }
    // split, and a side no sharing cell lies on (the mesh boundary): all
    if (slots.empty())
        for (std::int32_t s = first; s < last; ++s)
            slots.push_back(s);
    factor = T(1) / T(slots.size());
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template void append_sharing_cells<float, int>(const MeshView<float, int>&, const EntityTable<int>&, std::int32_t,
                                               EntityRules<float>&);
template void append_sharing_cells<double, int>(const MeshView<double, int>&, const EntityTable<int>&, std::int32_t,
                                                EntityRules<double>&);
template EntityRules<float> whole_entity_rules<float, int>(const MeshPart<float, int>&, int);
template EntityRules<double> whole_entity_rules<double, int>(const MeshPart<double, int>&, int);
template void append_cell_contribution<float>(const EntityRules<float>&, int, std::int32_t, float,
                                              std::vector<float>&, std::vector<float>&);
template void append_cell_contribution<double>(const EntityRules<double>&, int, std::int32_t, double,
                                               std::vector<double>&, std::vector<double>&);
template void policy_cells<float>(const EntityRules<float>&, int, EntityPolicy, std::vector<std::int32_t>&, float&);
template void policy_cells<double>(const EntityRules<double>&, int, EntityPolicy, std::vector<std::int32_t>&,
                                   double&);

} // namespace cutcells::part
