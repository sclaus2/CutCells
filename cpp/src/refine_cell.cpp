// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "refine_cell.h"
#include "cell_subdivision.h"
#include "cell_topology.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <map>
#include <stdexcept>

namespace cutcells
{
namespace
{

template <std::floating_point T>
struct EdgeStateData
{
    std::array<std::int32_t, 2> verts = {0, 0};
    std::vector<EdgeRootTag> tags;
    std::vector<T> split_params;
    std::vector<std::uint8_t> split_has_value;
    std::vector<T> root_params;
    std::vector<std::int32_t> root_vertex_ids;
    std::vector<std::uint8_t> root_has_value;
};

template <std::floating_point T>
struct CapturedEdgeState
{
    int num_level_sets = 0;
    std::vector<std::pair<std::pair<int, int>, EdgeStateData<T>>> edges_in_old_order;
    std::map<std::pair<int, int>, const EdgeStateData<T>*> by_key;
};

template <std::floating_point T>
CapturedEdgeState<T> capture_edge_state(const AdaptCell<T>& adapt_cell)
{
    CapturedEdgeState<T> state;
    const int n_edges = adapt_cell.n_entities(1);
    state.num_level_sets = adapt_cell.edge_root_tag_num_level_sets;
    state.edges_in_old_order.reserve(static_cast<std::size_t>(n_edges));

    for (int e = 0; e < n_edges; ++e)
    {
        auto ev = adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(e)];
        const std::pair<int, int> key = {
            std::min(static_cast<int>(ev[0]), static_cast<int>(ev[1])),
            std::max(static_cast<int>(ev[0]), static_cast<int>(ev[1]))};

        EdgeStateData<T> data;
        data.verts = {ev[0], ev[1]};
        data.tags.resize(static_cast<std::size_t>(state.num_level_sets),
                         EdgeRootTag::not_classified);
        data.split_params.resize(static_cast<std::size_t>(state.num_level_sets), T(0));
        data.split_has_value.resize(static_cast<std::size_t>(state.num_level_sets),
                                    std::uint8_t(0));
        data.root_params.resize(static_cast<std::size_t>(state.num_level_sets), T(0));
        data.root_vertex_ids.resize(static_cast<std::size_t>(state.num_level_sets),
                                    std::int32_t(-1));
        data.root_has_value.resize(static_cast<std::size_t>(state.num_level_sets),
                                   std::uint8_t(0));

        for (int ls = 0; ls < state.num_level_sets; ++ls)
        {
            data.tags[static_cast<std::size_t>(ls)] =
                adapt_cell.get_edge_root_tag(ls, e);

            const auto idx = static_cast<std::size_t>(ls * n_edges + e);
            if (idx < adapt_cell.edge_green_split_param.size())
                data.split_params[static_cast<std::size_t>(ls)] =
                    adapt_cell.edge_green_split_param[idx];
            if (idx < adapt_cell.edge_green_split_has_value.size())
                data.split_has_value[static_cast<std::size_t>(ls)] =
                    adapt_cell.edge_green_split_has_value[idx];
            if (idx < adapt_cell.edge_one_root_param.size())
                data.root_params[static_cast<std::size_t>(ls)] =
                    adapt_cell.edge_one_root_param[idx];
            if (idx < adapt_cell.edge_one_root_vertex_id.size())
                data.root_vertex_ids[static_cast<std::size_t>(ls)] =
                    adapt_cell.edge_one_root_vertex_id[idx];
            if (idx < adapt_cell.edge_one_root_has_value.size())
                data.root_has_value[static_cast<std::size_t>(ls)] =
                    adapt_cell.edge_one_root_has_value[idx];
        }

        state.edges_in_old_order.emplace_back(key, std::move(data));
    }

    for (auto& [key, data] : state.edges_in_old_order)
        state.by_key[key] = &data;

    return state;
}

template <std::floating_point T>
void clear_topology_caches(AdaptCell<T>& adapt_cell)
{
    for (auto& row : adapt_cell.connectivity)
    {
        for (auto& conn : row)
        {
            conn.offsets.clear();
            conn.indices.clear();
        }
    }
    adapt_cell.has_connectivity = {};

    adapt_cell.zero_entity_dim.clear();
    adapt_cell.zero_entity_id.clear();
    adapt_cell.zero_entity_zero_mask.clear();
    adapt_cell.zero_entity_is_owned.clear();
    adapt_cell.zero_entity_parent_dim.clear();
    adapt_cell.zero_entity_parent_id.clear();
    adapt_cell.zero_entity_host_cell_id.clear();
    adapt_cell.zero_entity_host_cell_type.clear();
    adapt_cell.zero_entity_host_face_id.clear();
    adapt_cell.zero_entity_source_level_set.clear();
    adapt_cell.zero_entity_host_cell_vertices.offsets.clear();
    adapt_cell.zero_entity_host_cell_vertices.indices.clear();
    ++adapt_cell.zero_entity_version;
}

template <std::floating_point T>
int append_vertex_with_parent_info(AdaptCell<T>& adapt_cell,
                                   std::span<const T> coords,
                                   int parent_dim,
                                   int parent_id,
                                   std::span<const T> parent_param,
                                   int source_edge_id)
{
    const int new_v = adapt_cell.n_vertices();
    adapt_cell.vertex_coords.insert(adapt_cell.vertex_coords.end(),
                                    coords.begin(), coords.end());
    adapt_cell.vertex_parent_dim.push_back(static_cast<std::int8_t>(parent_dim));
    adapt_cell.vertex_parent_id.push_back(static_cast<std::int32_t>(parent_id));
    const std::int32_t param_offset =
        static_cast<std::int32_t>(adapt_cell.vertex_parent_param.size());
    adapt_cell.vertex_parent_param.insert(adapt_cell.vertex_parent_param.end(),
                                          parent_param.begin(),
                                          parent_param.end());
    adapt_cell.vertex_parent_param_offset.push_back(
        param_offset + static_cast<std::int32_t>(parent_param.size()));
    adapt_cell.zero_mask_per_vertex.push_back(0);
    adapt_cell.negative_mask_per_vertex.push_back(0);
    adapt_cell.vertex_source_edge_id.push_back(static_cast<std::int32_t>(source_edge_id));
    return new_v;
}

template <std::floating_point T>
bool vertex_parameter_on_parent_edge(const AdaptCell<T>& adapt_cell,
                                     int vertex_id,
                                     int parent_edge_id,
                                     T& t)
{
    const int dim = adapt_cell.vertex_parent_dim[static_cast<std::size_t>(vertex_id)];
    if (dim == 0)
    {
        const auto parent_edge =
            cell::edges(adapt_cell.parent_cell_type)[static_cast<std::size_t>(parent_edge_id)];
        if (parent_edge[0] == vertex_id)
        {
            t = T(0);
            return true;
        }
        if (parent_edge[1] == vertex_id)
        {
            t = T(1);
            return true;
        }
        return false;
    }

    if (dim != 1)
        return false;
    if (adapt_cell.vertex_parent_id[static_cast<std::size_t>(vertex_id)] != parent_edge_id)
        return false;

    const auto begin =
        static_cast<std::size_t>(adapt_cell.vertex_parent_param_offset[static_cast<std::size_t>(vertex_id)]);
    const auto end =
        static_cast<std::size_t>(adapt_cell.vertex_parent_param_offset[static_cast<std::size_t>(vertex_id + 1)]);
    if (end - begin != 1)
        return false;

    t = adapt_cell.vertex_parent_param[begin];
    return true;
}

template <std::floating_point T>
bool edge_is_on_single_parent_edge(const AdaptCell<T>& adapt_cell,
                                   int edge_id,
                                   int& parent_edge_id)
{
    auto edge_vertices = adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(edge_id)];
    if (edge_vertices.size() != 2)
        return false;

    const auto parent_edges = cell::edges(adapt_cell.parent_cell_type);
    for (int pe = 0; pe < static_cast<int>(parent_edges.size()); ++pe)
    {
        T unused = T(0);
        if (vertex_parameter_on_parent_edge(adapt_cell, edge_vertices[0], pe, unused)
            && vertex_parameter_on_parent_edge(adapt_cell, edge_vertices[1], pe, unused))
        {
            parent_edge_id = pe;
            return true;
        }
    }

    return false;
}

template <std::floating_point T>
int append_interpolated_vertex_on_edge(AdaptCell<T>& adapt_cell, int edge_id, T t)
{
    auto edge_vertices = adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(edge_id)];
    if (edge_vertices.size() != 2)
        throw std::runtime_error("append_interpolated_vertex_on_edge: invalid edge");
    const int v0 = edge_vertices[0];
    const int v1 = edge_vertices[1];

    std::vector<T> x(static_cast<std::size_t>(adapt_cell.tdim), T(0));
    for (int d = 0; d < adapt_cell.tdim; ++d)
    {
        const T x0 = adapt_cell.vertex_coords[static_cast<std::size_t>(v0 * adapt_cell.tdim + d)];
        const T x1 = adapt_cell.vertex_coords[static_cast<std::size_t>(v1 * adapt_cell.tdim + d)];
        x[static_cast<std::size_t>(d)] = (T(1) - t) * x0 + t * x1;
    }

    int parent_dim = adapt_cell.tdim;
    int parent_id = adapt_cell.parent_cell_id;
    std::vector<T> parent_param(x.begin(), x.end());

    int parent_edge_id = -1;
    T parent_t0 = T(0);
    T parent_t1 = T(1);
    if (edge_is_on_single_parent_edge(adapt_cell, edge_id, parent_edge_id)
        && vertex_parameter_on_parent_edge(adapt_cell, v0, parent_edge_id, parent_t0)
        && vertex_parameter_on_parent_edge(adapt_cell, v1, parent_edge_id, parent_t1))
    {
        parent_dim = 1;
        parent_id = parent_edge_id;
        parent_param.assign(1, (T(1) - t) * parent_t0 + t * parent_t1);
    }

    return append_vertex_with_parent_info(
        adapt_cell,
        std::span<const T>(x),
        parent_dim,
        parent_id,
        std::span<const T>(parent_param),
        edge_id);
}

template <std::floating_point T>
int append_cell_center_vertex(AdaptCell<T>& adapt_cell, std::span<const std::int32_t> verts)
{
    std::vector<T> x(static_cast<std::size_t>(adapt_cell.tdim), T(0));
    const T inv = T(1) / T(static_cast<int>(verts.size()));
    for (auto v : verts)
    {
        for (int d = 0; d < adapt_cell.tdim; ++d)
            x[static_cast<std::size_t>(d)] +=
                adapt_cell.vertex_coords[static_cast<std::size_t>(v * adapt_cell.tdim + d)] * inv;
    }
    return append_vertex_with_parent_info(
        adapt_cell,
        std::span<const T>(x),
        adapt_cell.tdim,
        adapt_cell.parent_cell_id,
        std::span<const T>(x),
        -1);
}

template <std::floating_point T>
std::map<std::pair<int, int>, int> build_edge_lookup(const AdaptCell<T>& adapt_cell)
{
    std::map<std::pair<int, int>, int> edge_lookup;
    const int n_edges = adapt_cell.n_entities(1);
    for (int e = 0; e < n_edges; ++e)
    {
        auto ev = adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(e)];
        const int a = std::min(static_cast<int>(ev[0]), static_cast<int>(ev[1]));
        const int b = std::max(static_cast<int>(ev[0]), static_cast<int>(ev[1]));
        edge_lookup[{a, b}] = e;
    }
    return edge_lookup;
}

template <std::floating_point T>
void rebuild_leaf_edges_preserve_certification(
    AdaptCell<T>& adapt_cell,
    const CapturedEdgeState<T>& old_edge_state)
{
    const int tdim = adapt_cell.tdim;
    std::map<std::pair<int, int>, std::array<std::int32_t, 2>> new_leaf_edges;
    std::vector<std::pair<int, int>> new_edge_keys_in_encounter_order;

    const int n_cells = adapt_cell.n_entities(tdim);
    for (int c = 0; c < n_cells; ++c)
    {
        const cell::type ctype = adapt_cell.entity_types[tdim][static_cast<std::size_t>(c)];
        auto cell_verts = adapt_cell.entity_to_vertex[tdim][static_cast<std::int32_t>(c)];
        auto cell_edges = cell::edges(ctype);

        for (const auto& ce : cell_edges)
        {
            const std::int32_t lv0 = cell_verts[static_cast<std::size_t>(ce[0])];
            const std::int32_t lv1 = cell_verts[static_cast<std::size_t>(ce[1])];
            const std::pair<int, int> key = {
                std::min(static_cast<int>(lv0), static_cast<int>(lv1)),
                std::max(static_cast<int>(lv0), static_cast<int>(lv1))};

            if (!new_leaf_edges.contains(key))
            {
                new_leaf_edges[key] = {lv0, lv1};
                new_edge_keys_in_encounter_order.push_back(key);
            }
        }
    }

    std::vector<std::array<std::int32_t, 2>> final_edges;
    std::vector<const EdgeStateData<T>*> final_old_state;
    final_edges.reserve(new_leaf_edges.size());
    final_old_state.reserve(new_leaf_edges.size());

    for (const auto& [key, old_data] : old_edge_state.edges_in_old_order)
    {
        auto it = new_leaf_edges.find(key);
        if (it == new_leaf_edges.end())
            continue;

        final_edges.push_back(old_data.verts);
        final_old_state.push_back(&old_data);
        new_leaf_edges.erase(it);
    }

    for (const auto& key : new_edge_keys_in_encounter_order)
    {
        auto it = new_leaf_edges.find(key);
        if (it == new_leaf_edges.end())
            continue;

        final_edges.push_back(it->second);
        final_old_state.push_back(nullptr);
    }

    adapt_cell.entity_types[1].assign(final_edges.size(), cell::type::interval);
    adapt_cell.entity_to_vertex[1].indices.clear();
    adapt_cell.entity_to_vertex[1].offsets.clear();
    adapt_cell.entity_to_vertex[1].offsets.push_back(0);
    for (const auto& edge : final_edges)
    {
        adapt_cell.entity_to_vertex[1].indices.push_back(edge[0]);
        adapt_cell.entity_to_vertex[1].indices.push_back(edge[1]);
        adapt_cell.entity_to_vertex[1].offsets.push_back(
            static_cast<std::int32_t>(adapt_cell.entity_to_vertex[1].indices.size()));
    }

    const int n_edges = static_cast<int>(final_edges.size());
    const int nls = old_edge_state.num_level_sets;
    adapt_cell.edge_root_tag_num_level_sets = nls;
    adapt_cell.edge_root_tag.assign(static_cast<std::size_t>(nls * n_edges),
                                    EdgeRootTag::not_classified);
    adapt_cell.edge_green_split_param.assign(static_cast<std::size_t>(nls * n_edges), T(0));
    adapt_cell.edge_green_split_has_value.assign(static_cast<std::size_t>(nls * n_edges),
                                                 std::uint8_t(0));
    adapt_cell.edge_one_root_param.assign(static_cast<std::size_t>(nls * n_edges), T(0));
    adapt_cell.edge_one_root_vertex_id.assign(static_cast<std::size_t>(nls * n_edges),
                                              std::int32_t(-1));
    adapt_cell.edge_one_root_has_value.assign(static_cast<std::size_t>(nls * n_edges),
                                              std::uint8_t(0));

    for (int e = 0; e < n_edges; ++e)
    {
        const auto* old_data = final_old_state[static_cast<std::size_t>(e)];
        if (!old_data)
            continue;

        for (int ls = 0; ls < nls; ++ls)
        {
            const auto idx = static_cast<std::size_t>(ls * n_edges + e);
            adapt_cell.edge_root_tag[idx] = old_data->tags[static_cast<std::size_t>(ls)];
            adapt_cell.edge_green_split_param[idx] =
                old_data->split_params[static_cast<std::size_t>(ls)];
            adapt_cell.edge_green_split_has_value[idx] =
                old_data->split_has_value[static_cast<std::size_t>(ls)];
            adapt_cell.edge_one_root_param[idx] =
                old_data->root_params[static_cast<std::size_t>(ls)];
            adapt_cell.edge_one_root_vertex_id[idx] =
                old_data->root_vertex_ids[static_cast<std::size_t>(ls)];
            adapt_cell.edge_one_root_has_value[idx] =
                old_data->root_has_value[static_cast<std::size_t>(ls)];
        }
    }
}

template <std::floating_point T>
void rebuild_leaf_cell_certification(
    AdaptCell<T>& adapt_cell,
    std::span<const CellCertTag> old_cell_tags,
    int old_num_level_sets,
    int old_num_cells,
    std::span<const int> old_cell_ids_for_new_cells)
{
    const int new_num_cells = adapt_cell.n_entities(adapt_cell.tdim);
    adapt_cell.cell_cert_tag_num_level_sets = old_num_level_sets;
    adapt_cell.cell_cert_tag.assign(
        static_cast<std::size_t>(old_num_level_sets * new_num_cells),
        CellCertTag::not_classified);

    for (int c = 0; c < new_num_cells; ++c)
    {
        const int old_c = old_cell_ids_for_new_cells[static_cast<std::size_t>(c)];
        if (old_c < 0 || old_c >= old_num_cells)
            continue;

        for (int ls = 0; ls < old_num_level_sets; ++ls)
        {
            adapt_cell.cell_cert_tag[static_cast<std::size_t>(ls * new_num_cells + c)] =
                old_cell_tags[static_cast<std::size_t>(ls * old_num_cells + old_c)];
        }
    }
}

template <std::floating_point T>
void rebuild_leaf_cell_provenance(
    AdaptCell<T>& adapt_cell,
    std::span<const cell::type> old_cell_types,
    const EntityAdjacency& old_cells,
    std::span<const std::int32_t> old_cell_source_cell_id,
    std::span<const std::int32_t> old_cell_refinement_generation,
    std::span<const CellRefinementReason> old_cell_refinement_reason,
    std::span<const std::int32_t> old_cell_host_parent_cell_id,
    const EntityAdjacency& old_host_cell_vertices,
    std::span<const cell::type> old_host_cell_types,
    std::span<const std::int32_t> old_source_level_sets,
    std::span<const int> old_cell_ids_for_new_cells,
    std::span<const int> source_cell_ids_for_new_cells,
    std::span<const CellRefinementReason> refinement_reasons_for_new_cells)
{
    const int tdim = adapt_cell.tdim;
    const int n_cells = adapt_cell.n_entities(tdim);

    adapt_cell.cell_source_cell_id.assign(
        static_cast<std::size_t>(n_cells), std::int32_t(-1));
    adapt_cell.cell_refinement_generation.assign(
        static_cast<std::size_t>(n_cells), std::int32_t(0));
    adapt_cell.cell_refinement_reason.assign(
        static_cast<std::size_t>(n_cells), CellRefinementReason::none);
    adapt_cell.cell_host_parent_cell_id.assign(
        static_cast<std::size_t>(n_cells), adapt_cell.parent_cell_id);

    adapt_cell.entity_host_cell_id[tdim].clear();
    adapt_cell.entity_host_cell_type[tdim].clear();
    adapt_cell.entity_host_face_id[tdim].clear();
    adapt_cell.entity_source_level_set[tdim].clear();
    adapt_cell.entity_host_cell_vertices[tdim].offsets.clear();
    adapt_cell.entity_host_cell_vertices[tdim].indices.clear();
    adapt_cell.entity_host_cell_vertices[tdim].offsets.push_back(std::int32_t(0));

    auto append_host_vertices = [&](std::span<const std::int32_t> vertices)
    {
        for (const auto v : vertices)
            adapt_cell.entity_host_cell_vertices[tdim].indices.push_back(v);
        adapt_cell.entity_host_cell_vertices[tdim].offsets.push_back(
            static_cast<std::int32_t>(
                adapt_cell.entity_host_cell_vertices[tdim].indices.size()));
    };

    const int old_num_cells = static_cast<int>(old_cell_types.size());
    for (int c = 0; c < n_cells; ++c)
    {
        const int old_c =
            (c < static_cast<int>(old_cell_ids_for_new_cells.size()))
                ? old_cell_ids_for_new_cells[static_cast<std::size_t>(c)]
                : -1;
        const int source_c =
            (c < static_cast<int>(source_cell_ids_for_new_cells.size()))
                ? source_cell_ids_for_new_cells[static_cast<std::size_t>(c)]
                : old_c;
        const CellRefinementReason reason =
            (c < static_cast<int>(refinement_reasons_for_new_cells.size()))
                ? refinement_reasons_for_new_cells[static_cast<std::size_t>(c)]
                : CellRefinementReason::none;
        const CellRefinementReason stored_reason =
            (reason == CellRefinementReason::none
             && old_c >= 0
             && old_c < static_cast<int>(old_cell_refinement_reason.size()))
                ? old_cell_refinement_reason[static_cast<std::size_t>(old_c)]
                : reason;
        const int meta_c =
            (source_c >= 0 && source_c < old_num_cells)
                ? source_c
                : ((old_c >= 0 && old_c < old_num_cells) ? old_c : -1);

        adapt_cell.cell_source_cell_id[static_cast<std::size_t>(c)] =
            static_cast<std::int32_t>(source_c);
        if (meta_c >= 0
            && meta_c < static_cast<int>(old_cell_refinement_generation.size()))
        {
            const bool unchanged_copy =
                old_c == source_c && old_c >= 0
                && reason == CellRefinementReason::none;
            adapt_cell.cell_refinement_generation[static_cast<std::size_t>(c)] =
                old_cell_refinement_generation[static_cast<std::size_t>(meta_c)]
                + (unchanged_copy ? 0 : 1);
        }
        adapt_cell.cell_refinement_reason[static_cast<std::size_t>(c)] = stored_reason;
        if (meta_c >= 0
            && meta_c < static_cast<int>(old_cell_host_parent_cell_id.size()))
        {
            adapt_cell.cell_host_parent_cell_id[static_cast<std::size_t>(c)] =
                old_cell_host_parent_cell_id[static_cast<std::size_t>(meta_c)];
        }

        adapt_cell.entity_host_cell_id[tdim].push_back(
            static_cast<std::int32_t>(source_c));
        const cell::type host_type =
            (meta_c >= 0 && meta_c < static_cast<int>(old_host_cell_types.size()))
                ? old_host_cell_types[static_cast<std::size_t>(meta_c)]
                : ((meta_c >= 0 && meta_c < static_cast<int>(old_cell_types.size()))
                ? old_cell_types[static_cast<std::size_t>(meta_c)]
                : adapt_cell.entity_types[tdim][static_cast<std::size_t>(c)]);
        adapt_cell.entity_host_cell_type[tdim].push_back(host_type);
        adapt_cell.entity_host_face_id[tdim].push_back(std::int32_t(-1));
        const int source_level_set =
            (meta_c >= 0 && meta_c < static_cast<int>(old_source_level_sets.size()))
                ? old_source_level_sets[static_cast<std::size_t>(meta_c)]
                : -1;
        adapt_cell.entity_source_level_set[tdim].push_back(
            reason == CellRefinementReason::cut_level_set
                ? source_level_set
                : std::int32_t(-1));

        std::span<const std::int32_t> host_vertices;
        if (reason == CellRefinementReason::cut_level_set
            && meta_c >= 0 && meta_c < old_cells.size())
        {
            host_vertices = old_cells[static_cast<std::int32_t>(meta_c)];
        }
        else if (old_c >= 0 && old_c < old_host_cell_vertices.size()
                 && reason == CellRefinementReason::none)
        {
            host_vertices = old_host_cell_vertices[static_cast<std::int32_t>(old_c)];
        }
        if (host_vertices.empty())
            host_vertices = adapt_cell.entity_to_vertex[tdim][static_cast<std::int32_t>(c)];
        append_host_vertices(host_vertices);
    }
}

void append_top_cell(std::vector<cell::type>& types,
                     EntityAdjacency& adj,
                     cell::type ctype,
                     std::span<const int> verts)
{
    types.push_back(ctype);
    for (int v : verts)
        adj.indices.push_back(static_cast<std::int32_t>(v));
    adj.offsets.push_back(static_cast<std::int32_t>(adj.indices.size()));
}

} // anonymous namespace

template <std::floating_point T>
void apply_topology_update_preserve_certification(
    AdaptCell<T>& adapt_cell,
    std::vector<cell::type>&& new_types,
    EntityAdjacency&& new_cells,
    std::span<const int> old_cell_ids_for_new_cells,
    std::span<const int> source_cell_ids_for_new_cells,
    std::span<const CellRefinementReason> refinement_reasons_for_new_cells)
{
    const int tdim = adapt_cell.tdim;
    const CapturedEdgeState<T> old_edge_state = capture_edge_state(adapt_cell);
    const int old_num_level_sets = adapt_cell.cell_cert_tag_num_level_sets;
    const int old_num_cells = adapt_cell.n_entities(tdim);
    const std::vector<CellCertTag> old_cell_tags = adapt_cell.cell_cert_tag;
    const std::vector<cell::type> old_cell_types = adapt_cell.entity_types[tdim];
    const EntityAdjacency old_cells = adapt_cell.entity_to_vertex[tdim];
    const std::vector<std::int32_t> old_cell_source_cell_id =
        adapt_cell.cell_source_cell_id;
    const std::vector<std::int32_t> old_cell_refinement_generation =
        adapt_cell.cell_refinement_generation;
    const std::vector<CellRefinementReason> old_cell_refinement_reason =
        adapt_cell.cell_refinement_reason;
    const std::vector<std::int32_t> old_cell_host_parent_cell_id =
        adapt_cell.cell_host_parent_cell_id;
    const EntityAdjacency old_host_cell_vertices =
        adapt_cell.entity_host_cell_vertices[tdim];
    const std::vector<cell::type> old_host_cell_types =
        adapt_cell.entity_host_cell_type[tdim];
    const std::vector<std::int32_t> old_source_level_sets =
        adapt_cell.entity_source_level_set[tdim];

    adapt_cell.entity_types[tdim] = std::move(new_types);
    adapt_cell.entity_to_vertex[tdim] = std::move(new_cells);

    rebuild_leaf_cell_provenance<T>(
        adapt_cell,
        std::span<const cell::type>(old_cell_types),
        old_cells,
        std::span<const std::int32_t>(old_cell_source_cell_id),
        std::span<const std::int32_t>(old_cell_refinement_generation),
        std::span<const CellRefinementReason>(old_cell_refinement_reason),
        std::span<const std::int32_t>(old_cell_host_parent_cell_id),
        old_host_cell_vertices,
        std::span<const cell::type>(old_host_cell_types),
        std::span<const std::int32_t>(old_source_level_sets),
        old_cell_ids_for_new_cells,
        source_cell_ids_for_new_cells,
        refinement_reasons_for_new_cells);
    rebuild_leaf_edges_preserve_certification(adapt_cell, old_edge_state);
    rebuild_leaf_cell_certification(adapt_cell,
                                    std::span<const CellCertTag>(old_cell_tags),
                                    old_num_level_sets,
                                    old_num_cells,
                                    old_cell_ids_for_new_cells);
    clear_topology_caches(adapt_cell);

    const int nls = std::max(adapt_cell.cell_cert_tag_num_level_sets,
                             adapt_cell.edge_root_tag_num_level_sets);
    recompute_active_level_set_masks(adapt_cell, nls);
    rebuild_zero_entity_inventory(adapt_cell);
}


// =====================================================================
// Invalidation helpers
// =====================================================================

template <std::floating_point T>
void invalidate_edge_tags_for_new_edges(AdaptCell<T>& adapt_cell,
                                        std::span<const int> new_edge_ids)
{
    const int nls = adapt_cell.edge_root_tag_num_level_sets;
    const int n_edges = adapt_cell.n_entities(1);
    for (int ls = 0; ls < nls; ++ls)
    {
        for (int e : new_edge_ids)
        {
            if (e >= 0 && e < n_edges)
                adapt_cell.set_edge_root_tag(ls, e, EdgeRootTag::not_classified);
        }
    }
}

template <std::floating_point T>
void invalidate_cell_tags_for_new_cells(AdaptCell<T>& adapt_cell,
                                        std::span<const int> new_cell_ids)
{
    const int nls = adapt_cell.cell_cert_tag_num_level_sets;
    const int tdim = adapt_cell.tdim;
    const int n_cells = adapt_cell.n_entities(tdim);
    for (int ls = 0; ls < nls; ++ls)
    {
        for (int c : new_cell_ids)
        {
            if (c >= 0 && c < n_cells)
                adapt_cell.set_cell_cert_tag(ls, c, CellCertTag::not_classified);
        }
    }
}

// =====================================================================
// Green refinement
// =====================================================================

template <std::floating_point T>
bool refine_green_on_multiple_root_edges(AdaptCell<T>& adapt_cell,
                                         int level_set_id)
{
    if (adapt_cell.tdim != 2 && adapt_cell.tdim != 3)
        return false;

    const int n_edges = adapt_cell.n_entities(1);
    int split_edge = -1;
    T split_t = T(0);
    for (int e = 0; e < n_edges; ++e)
    {
        if (adapt_cell.get_edge_root_tag(level_set_id, e) != EdgeRootTag::multiple_roots)
            continue;
        const auto idx = static_cast<std::size_t>(level_set_id * n_edges + e);
        if (!adapt_cell.edge_green_split_has_value[idx])
            continue;
        split_edge = e;
        split_t = adapt_cell.edge_green_split_param[idx];
        break;
    }

    if (split_edge < 0)
        return false;

    const int tdim = adapt_cell.tdim;
    const std::vector<cell::type> old_types = adapt_cell.entity_types[tdim];
    const EntityAdjacency old_cells = adapt_cell.entity_to_vertex[tdim];

    if (tdim == 2
        && std::any_of(old_types.begin(), old_types.end(),
                       [](cell::type ct) { return ct != cell::type::triangle; }))
    {
        return false;
    }
    if (tdim == 3
        && std::any_of(old_types.begin(), old_types.end(),
                       [](cell::type ct) { return ct != cell::type::tetrahedron; }))
    {
        return false;
    }

    auto ev = adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(split_edge)];
    const int v0 = ev[0];
    const int v1 = ev[1];
    const int new_v = append_interpolated_vertex_on_edge(adapt_cell, split_edge, split_t);

    EntityAdjacency new_cells;
    new_cells.offsets.push_back(0);
    std::vector<cell::type> new_types;
    std::vector<int> old_cell_ids_for_new_cells;
    std::vector<int> source_cell_ids_for_new_cells;
    std::vector<CellRefinementReason> refinement_reasons_for_new_cells;

    for (int c = 0; c < static_cast<int>(old_types.size()); ++c)
    {
        auto verts = old_cells[static_cast<std::int32_t>(c)];

        if (tdim == 2)
        {
            const int a = verts[0];
            const int b = verts[1];
            const int cvert = verts[2];
            bool split = false;

            if ((a == v0 && b == v1) || (a == v1 && b == v0))
            {
                const std::array<int, 3> t0 = {a, new_v, cvert};
                const std::array<int, 3> t1 = {new_v, b, cvert};
                append_top_cell(new_types, new_cells, cell::type::triangle, std::span<const int>(t0));
                old_cell_ids_for_new_cells.push_back(-1);
                source_cell_ids_for_new_cells.push_back(c);
                refinement_reasons_for_new_cells.push_back(CellRefinementReason::green_edge);
                append_top_cell(new_types, new_cells, cell::type::triangle, std::span<const int>(t1));
                old_cell_ids_for_new_cells.push_back(-1);
                source_cell_ids_for_new_cells.push_back(c);
                refinement_reasons_for_new_cells.push_back(CellRefinementReason::green_edge);
                split = true;
            }
            else if ((b == v0 && cvert == v1) || (b == v1 && cvert == v0))
            {
                const std::array<int, 3> t0 = {b, new_v, a};
                const std::array<int, 3> t1 = {new_v, cvert, a};
                append_top_cell(new_types, new_cells, cell::type::triangle, std::span<const int>(t0));
                old_cell_ids_for_new_cells.push_back(-1);
                source_cell_ids_for_new_cells.push_back(c);
                refinement_reasons_for_new_cells.push_back(CellRefinementReason::green_edge);
                append_top_cell(new_types, new_cells, cell::type::triangle, std::span<const int>(t1));
                old_cell_ids_for_new_cells.push_back(-1);
                source_cell_ids_for_new_cells.push_back(c);
                refinement_reasons_for_new_cells.push_back(CellRefinementReason::green_edge);
                split = true;
            }
            else if ((cvert == v0 && a == v1) || (cvert == v1 && a == v0))
            {
                const std::array<int, 3> t0 = {cvert, new_v, b};
                const std::array<int, 3> t1 = {new_v, a, b};
                append_top_cell(new_types, new_cells, cell::type::triangle, std::span<const int>(t0));
                old_cell_ids_for_new_cells.push_back(-1);
                source_cell_ids_for_new_cells.push_back(c);
                refinement_reasons_for_new_cells.push_back(CellRefinementReason::green_edge);
                append_top_cell(new_types, new_cells, cell::type::triangle, std::span<const int>(t1));
                old_cell_ids_for_new_cells.push_back(-1);
                source_cell_ids_for_new_cells.push_back(c);
                refinement_reasons_for_new_cells.push_back(CellRefinementReason::green_edge);
                split = true;
            }

            if (!split)
            {
                std::vector<int> copy(verts.begin(), verts.end());
                append_top_cell(new_types, new_cells, cell::type::triangle, std::span<const int>(copy));
                old_cell_ids_for_new_cells.push_back(c);
                source_cell_ids_for_new_cells.push_back(c);
                refinement_reasons_for_new_cells.push_back(CellRefinementReason::none);
            }
        }
        else
        {
            const int a = verts[0];
            const int b = verts[1];
            const int cvert = verts[2];
            const int dvert = verts[3];
            bool split = false;

            const std::array<std::array<int, 4>, 6> edge_positions = {{
                {{0, 1, 2, 3}},
                {{0, 2, 1, 3}},
                {{0, 3, 1, 2}},
                {{1, 2, 0, 3}},
                {{1, 3, 0, 2}},
                {{2, 3, 0, 1}},
            }};

            const std::array<int, 4> tet = {a, b, cvert, dvert};
            for (const auto& pos : edge_positions)
            {
                const int e0 = tet[static_cast<std::size_t>(pos[0])];
                const int e1 = tet[static_cast<std::size_t>(pos[1])];
                if (!((e0 == v0 && e1 == v1) || (e0 == v1 && e1 == v0)))
                    continue;

                // Preserve orientation: replace one endpoint at a time
                // in the parent's vertex ordering.
                std::array<int, 4> t0 = tet;
                t0[static_cast<std::size_t>(pos[1])] = new_v;
                std::array<int, 4> t1 = tet;
                t1[static_cast<std::size_t>(pos[0])] = new_v;
                append_top_cell(new_types, new_cells, cell::type::tetrahedron, std::span<const int>(t0));
                old_cell_ids_for_new_cells.push_back(-1);
                source_cell_ids_for_new_cells.push_back(c);
                refinement_reasons_for_new_cells.push_back(CellRefinementReason::green_edge);
                append_top_cell(new_types, new_cells, cell::type::tetrahedron, std::span<const int>(t1));
                old_cell_ids_for_new_cells.push_back(-1);
                source_cell_ids_for_new_cells.push_back(c);
                refinement_reasons_for_new_cells.push_back(CellRefinementReason::green_edge);
                split = true;
                break;
            }

            if (!split)
            {
                std::vector<int> copy(verts.begin(), verts.end());
                append_top_cell(new_types, new_cells, cell::type::tetrahedron, std::span<const int>(copy));
                old_cell_ids_for_new_cells.push_back(c);
                source_cell_ids_for_new_cells.push_back(c);
                refinement_reasons_for_new_cells.push_back(CellRefinementReason::none);
            }
        }
    }

    apply_topology_update_preserve_certification(
        adapt_cell, std::move(new_types), std::move(new_cells),
        std::span<const int>(old_cell_ids_for_new_cells),
        std::span<const int>(source_cell_ids_for_new_cells),
        std::span<const CellRefinementReason>(refinement_reasons_for_new_cells));
    return true;
}

// =====================================================================
// Red refinement
// =====================================================================

template <std::floating_point T>
bool refine_red_on_ambiguous_cells(AdaptCell<T>& adapt_cell,
                                   int level_set_id)
{
    const int tdim = adapt_cell.tdim;
    const int n_cells = adapt_cell.n_entities(tdim);

    std::vector<int> ambiguous_cells;
    for (int c = 0; c < n_cells; ++c)
    {
        if (adapt_cell.get_cell_cert_tag(level_set_id, c)
            == CellCertTag::ambiguous)
            ambiguous_cells.push_back(c);
    }

    if (ambiguous_cells.empty())
        return false;

    const std::vector<cell::type> old_types = adapt_cell.entity_types[tdim];
    const EntityAdjacency old_cells = adapt_cell.entity_to_vertex[tdim];
    const auto edge_lookup = build_edge_lookup(adapt_cell);

    std::map<int, int> midpoint_vertex_by_edge;
    auto get_midpoint_vertex = [&](int edge_id) -> int
    {
        auto it = midpoint_vertex_by_edge.find(edge_id);
        if (it != midpoint_vertex_by_edge.end())
            return it->second;

        const int mid = append_interpolated_vertex_on_edge(adapt_cell, edge_id, T(0.5));
        midpoint_vertex_by_edge[edge_id] = mid;
        return mid;
    };

    EntityAdjacency new_cells;
    new_cells.offsets.push_back(0);
    std::vector<cell::type> new_types;
    std::vector<int> old_cell_ids_for_new_cells;
    std::vector<int> source_cell_ids_for_new_cells;
    std::vector<CellRefinementReason> refinement_reasons_for_new_cells;

    for (int c = 0; c < n_cells; ++c)
    {
        auto verts = old_cells[static_cast<std::int32_t>(c)];
        const cell::type ctype = old_types[static_cast<std::size_t>(c)];
        const bool is_ambiguous =
            (adapt_cell.get_cell_cert_tag(level_set_id, c) == CellCertTag::ambiguous);

        if (!is_ambiguous)
        {
            std::vector<int> copy(verts.begin(), verts.end());
            append_top_cell(new_types, new_cells, ctype, std::span<const int>(copy));
            old_cell_ids_for_new_cells.push_back(c);
            source_cell_ids_for_new_cells.push_back(c);
            refinement_reasons_for_new_cells.push_back(CellRefinementReason::none);
            continue;
        }

        switch (ctype)
        {
        case cell::type::interval:
        {
            const int a = verts[0];
            const int b = verts[1];
            const std::pair<int, int> key = {std::min(a, b), std::max(a, b)};
            auto it = edge_lookup.find(key);
            if (it == edge_lookup.end())
            {
                throw std::runtime_error(
                    "refine_red_on_ambiguous_cells: missing interval edge key=("
                    + std::to_string(key.first) + "," + std::to_string(key.second)
                    + ") for cell " + std::to_string(c));
            }
            const int e = it->second;
            const int m = get_midpoint_vertex(e);
            const std::array<std::array<int, 2>, 2> children = {{{a, m}, {m, b}}};
            for (const auto& child : children)
            {
                append_top_cell(new_types, new_cells, cell::type::interval, std::span<const int>(child));
                old_cell_ids_for_new_cells.push_back(-1);
                source_cell_ids_for_new_cells.push_back(c);
                refinement_reasons_for_new_cells.push_back(CellRefinementReason::red_cell);
            }
            break;
        }
        case cell::type::triangle:
        {
            std::array<int, 6> local_to_global = {
                static_cast<int>(verts[0]), static_cast<int>(verts[1]), static_cast<int>(verts[2]),
                -1, -1, -1};
            const auto ledges = cell::edges(cell::type::triangle);
            for (int le = 0; le < 3; ++le)
            {
                const int a = verts[static_cast<std::size_t>(ledges[static_cast<std::size_t>(le)][0])];
                const int b = verts[static_cast<std::size_t>(ledges[static_cast<std::size_t>(le)][1])];
                const std::pair<int, int> key = {std::min(a, b), std::max(a, b)};
                auto it = edge_lookup.find(key);
                if (it == edge_lookup.end())
                {
                    throw std::runtime_error(
                        "refine_red_on_ambiguous_cells: missing triangle edge key=("
                        + std::to_string(key.first) + "," + std::to_string(key.second)
                        + ") for cell " + std::to_string(c) + " local_edge="
                        + std::to_string(le));
                }
                local_to_global[static_cast<std::size_t>(3 + le)] =
                    get_midpoint_vertex(it->second);
            }

            for (const auto& child : cell::triangle_subdivision_table)
            {
                const std::array<int, 3> g = {
                    local_to_global[static_cast<std::size_t>(child[0])],
                    local_to_global[static_cast<std::size_t>(child[1])],
                    local_to_global[static_cast<std::size_t>(child[2])]};
                append_top_cell(new_types, new_cells, cell::type::triangle, std::span<const int>(g));
                old_cell_ids_for_new_cells.push_back(-1);
                source_cell_ids_for_new_cells.push_back(c);
                refinement_reasons_for_new_cells.push_back(CellRefinementReason::red_cell);
            }
            break;
        }
        case cell::type::quadrilateral:
        {
            std::array<int, 9> local_to_global = {
                static_cast<int>(verts[0]), static_cast<int>(verts[1]),
                static_cast<int>(verts[2]), static_cast<int>(verts[3]),
                -1, -1, -1, -1, -1};
            const auto ledges = cell::edges(cell::type::quadrilateral);
            for (int le = 0; le < 4; ++le)
            {
                const int a = verts[static_cast<std::size_t>(ledges[static_cast<std::size_t>(le)][0])];
                const int b = verts[static_cast<std::size_t>(ledges[static_cast<std::size_t>(le)][1])];
                const std::pair<int, int> key = {std::min(a, b), std::max(a, b)};
                auto it = edge_lookup.find(key);
                if (it == edge_lookup.end())
                {
                    throw std::runtime_error(
                        "refine_red_on_ambiguous_cells: missing quadrilateral edge key=("
                        + std::to_string(key.first) + "," + std::to_string(key.second)
                        + ") for cell " + std::to_string(c) + " local_edge="
                        + std::to_string(le));
                }
                local_to_global[static_cast<std::size_t>(4 + le)] =
                    get_midpoint_vertex(it->second);
            }
            local_to_global[8] = append_cell_center_vertex(adapt_cell, verts);

            for (const auto& child : cell::quadrilateral_subdivision_table)
            {
                const std::array<int, 4> g = {
                    local_to_global[static_cast<std::size_t>(child[0])],
                    local_to_global[static_cast<std::size_t>(child[1])],
                    local_to_global[static_cast<std::size_t>(child[2])],
                    local_to_global[static_cast<std::size_t>(child[3])]};
                append_top_cell(new_types, new_cells, cell::type::quadrilateral, std::span<const int>(g));
                old_cell_ids_for_new_cells.push_back(-1);
                source_cell_ids_for_new_cells.push_back(c);
                refinement_reasons_for_new_cells.push_back(CellRefinementReason::red_cell);
            }
            break;
        }
        case cell::type::tetrahedron:
        {
            std::array<int, 10> local_to_global = {
                static_cast<int>(verts[0]), static_cast<int>(verts[1]),
                static_cast<int>(verts[2]), static_cast<int>(verts[3]),
                -1, -1, -1, -1, -1, -1};
            const auto ledges = cell::edges(cell::type::tetrahedron);
            for (int le = 0; le < 6; ++le)
            {
                const int a = verts[static_cast<std::size_t>(ledges[static_cast<std::size_t>(le)][0])];
                const int b = verts[static_cast<std::size_t>(ledges[static_cast<std::size_t>(le)][1])];
                const std::pair<int, int> key = {std::min(a, b), std::max(a, b)};
                auto it = edge_lookup.find(key);
                if (it == edge_lookup.end())
                {
                    throw std::runtime_error(
                        "refine_red_on_ambiguous_cells: missing tetrahedron edge key=("
                        + std::to_string(key.first) + "," + std::to_string(key.second)
                        + ") for cell " + std::to_string(c) + " local_edge="
                        + std::to_string(le));
                }
                local_to_global[static_cast<std::size_t>(4 + le)] =
                    get_midpoint_vertex(it->second);
            }

            for (const auto& child : cell::tetrahedron_subdivision_table)
            {
                const std::array<int, 4> g = {
                    local_to_global[static_cast<std::size_t>(child[0])],
                    local_to_global[static_cast<std::size_t>(child[1])],
                    local_to_global[static_cast<std::size_t>(child[2])],
                    local_to_global[static_cast<std::size_t>(child[3])]};
                append_top_cell(new_types, new_cells, cell::type::tetrahedron, std::span<const int>(g));
                old_cell_ids_for_new_cells.push_back(-1);
                source_cell_ids_for_new_cells.push_back(c);
                refinement_reasons_for_new_cells.push_back(CellRefinementReason::red_cell);
            }
            break;
        }
        default:
            throw std::runtime_error(
                "refine_red_on_ambiguous_cells: unsupported cell type");
        }
    }

    apply_topology_update_preserve_certification(
        adapt_cell, std::move(new_types), std::move(new_cells),
        std::span<const int>(old_cell_ids_for_new_cells),
        std::span<const int>(source_cell_ids_for_new_cells),
        std::span<const CellRefinementReason>(refinement_reasons_for_new_cells));
    return true;
}

// =====================================================================
// Explicit template instantiations
// =====================================================================

template void invalidate_edge_tags_for_new_edges(AdaptCell<double>&,
                                                  std::span<const int>);
template void invalidate_edge_tags_for_new_edges(AdaptCell<float>&,
                                                  std::span<const int>);

template void invalidate_cell_tags_for_new_cells(AdaptCell<double>&,
                                                  std::span<const int>);
template void invalidate_cell_tags_for_new_cells(AdaptCell<float>&,
                                                  std::span<const int>);

template <std::floating_point T>
void invalidate_face_tags_for_new_faces(AdaptCell<T>& adapt_cell,
                                        std::span<const int> new_face_ids)
{
    const int nls = adapt_cell.face_cert_tag_num_level_sets;
    const int n_faces = adapt_cell.n_entities(2);
    for (int ls = 0; ls < nls; ++ls)
    {
        for (int f : new_face_ids)
        {
            if (f >= 0 && f < n_faces)
                adapt_cell.set_face_cert_tag(ls, f, FaceCertTag::not_classified);
        }
    }
}

template void invalidate_face_tags_for_new_faces(AdaptCell<double>&,
                                                  std::span<const int>);
template void invalidate_face_tags_for_new_faces(AdaptCell<float>&,
                                                  std::span<const int>);

template bool refine_green_on_multiple_root_edges(AdaptCell<double>&, int);
template bool refine_green_on_multiple_root_edges(AdaptCell<float>&, int);

template bool refine_red_on_ambiguous_cells(AdaptCell<double>&, int);
template bool refine_red_on_ambiguous_cells(AdaptCell<float>&, int);

template void apply_topology_update_preserve_certification(
    AdaptCell<double>&,
    std::vector<cell::type>&&,
    EntityAdjacency&&,
    std::span<const int>,
    std::span<const int>,
    std::span<const CellRefinementReason>);
template void apply_topology_update_preserve_certification(
    AdaptCell<float>&,
    std::vector<cell::type>&&,
    EntityAdjacency&&,
    std::span<const int>,
    std::span<const int>,
    std::span<const CellRefinementReason>);

} // namespace cutcells
