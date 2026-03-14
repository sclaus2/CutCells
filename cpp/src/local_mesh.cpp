// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "local_mesh.h"
#include "cell_topology.h"
#include "mapping.h"
#include "cell_flags.h"
#include "cut_cell.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace cutcells
{

namespace
{
struct EdgeKey
{
    int32_t a = -1;
    int32_t b = -1;

    bool operator==(const EdgeKey& other) const noexcept
    {
        return a == other.a && b == other.b;
    }
};

struct EdgeKeyHash
{
    std::size_t operator()(const EdgeKey& k) const noexcept
    {
        const auto h0 = std::hash<int32_t>{}(k.a);
        const auto h1 = std::hash<int32_t>{}(k.b);
        return h0 ^ (h1 + 0x9e3779b97f4a7c15ULL + (h0 << 6) + (h0 >> 2));
    }
};

template <std::floating_point T>
inline uint64_t bits_key(T v)
{
    if constexpr (sizeof(T) == sizeof(uint64_t))
        return std::bit_cast<uint64_t>(v);
    else
        return static_cast<uint64_t>(std::bit_cast<uint32_t>(v));
}

template <std::floating_point T>
struct VertexDedupKey
{
    int32_t parent_dim = -1;
    int32_t parent_id = -1;
    std::array<uint64_t, 3> x = {0, 0, 0};
    std::array<uint64_t, 3> xref = {0, 0, 0};
    uint8_t gdim = 0;
    uint8_t tdim = 0;
    uint8_t has_ref = 0;

    bool operator==(const VertexDedupKey& other) const noexcept
    {
        return parent_dim == other.parent_dim
               && parent_id == other.parent_id
               && x == other.x
               && xref == other.xref
               && gdim == other.gdim
               && tdim == other.tdim
               && has_ref == other.has_ref;
    }
};

template <std::floating_point T>
struct VertexDedupKeyHash
{
    std::size_t operator()(const VertexDedupKey<T>& k) const noexcept
    {
        std::size_t h = std::hash<int32_t>{}(k.parent_dim);
        auto mix = [&](std::size_t v) {
            h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        };
        mix(std::hash<int32_t>{}(k.parent_id));
        mix(std::hash<uint64_t>{}(k.x[0]));
        mix(std::hash<uint64_t>{}(k.x[1]));
        mix(std::hash<uint64_t>{}(k.x[2]));
        mix(std::hash<uint64_t>{}(k.xref[0]));
        mix(std::hash<uint64_t>{}(k.xref[1]));
        mix(std::hash<uint64_t>{}(k.xref[2]));
        mix(std::hash<uint8_t>{}(k.gdim));
        mix(std::hash<uint8_t>{}(k.tdim));
        mix(std::hash<uint8_t>{}(k.has_ref));
        return h;
    }
};

template <std::floating_point T>
void build_local_edges(LocalMesh<T>& mesh)
{
    const int nc = mesh.n_cells();
    std::vector<std::pair<int32_t, int32_t>> unique_edges;
    std::unordered_map<EdgeKey, int32_t, EdgeKeyHash> edge_to_id;
    edge_to_id.reserve(static_cast<std::size_t>(nc * 8));
    mesh.cell_edge_offsets.resize(nc + 1);
    mesh.cell_edges_flat.clear();

    for (int i = 0; i < nc; ++i)
    {
        mesh.cell_edge_offsets[i] = static_cast<int32_t>(mesh.cell_edges_flat.size());
        const int v_start = mesh.cell_offsets[i];
        const auto cell_typ = mesh.cell_types[i];
        const auto edge_patterns = cell::edges(cell_typ);
        for (const auto& ep : edge_patterns)
        {
            int32_t v0 = mesh.cell_vertices[v_start + ep[0]];
            int32_t v1 = mesh.cell_vertices[v_start + ep[1]];
            if (v0 > v1)
                std::swap(v0, v1);

            const EdgeKey key{v0, v1};
            const auto it = edge_to_id.find(key);
            int32_t edge_idx = 0;
            if (it == edge_to_id.end())
            {
                edge_idx = static_cast<int32_t>(unique_edges.size());
                unique_edges.push_back({v0, v1});
                edge_to_id.emplace(key, edge_idx);
            }
            else
            {
                edge_idx = it->second;
            }
            mesh.cell_edges_flat.push_back(edge_idx);
        }
    }
    mesh.cell_edge_offsets[nc] = static_cast<int32_t>(mesh.cell_edges_flat.size());

    const int n_edges = static_cast<int>(unique_edges.size());
    mesh.edge_vertices.resize(n_edges * 2);
    for (int i = 0; i < n_edges; ++i)
    {
        mesh.edge_vertices[i * 2] = unique_edges[i].first;
        mesh.edge_vertices[i * 2 + 1] = unique_edges[i].second;
    }
    mesh.edge_parent_dim.assign(n_edges, 1);
    mesh.edge_parent_id.resize(n_edges);
    for (int i = 0; i < n_edges; ++i)
        mesh.edge_parent_id[i] = i;
    mesh.edge_state.assign(n_edges, static_cast<uint8_t>(EdgeState::uncertain));
    mesh.edge_root_vertex.assign(n_edges, -1);
    mesh.edge_root_parameter.assign(n_edges, T(-1));
    mesh.edge_root_iterations.assign(n_edges, 0);
    mesh.edge_root_evaluations.assign(n_edges, 0);
    mesh.edge_root_converged.assign(n_edges, 0);
    mesh.edge_root_residual.assign(n_edges, T(0));
}

inline std::span<const std::array<int, 2>> cut_edge_order(cell::type ct)
{
    // Edge numbering follows the clipping/cut tables (VTK-style), not Basix.
    static constexpr std::array<std::array<int, 2>, 1> interval = {{{0, 1}}};
    static constexpr std::array<std::array<int, 2>, 3> triangle = {{{0, 1}, {1, 2}, {2, 0}}};
    static constexpr std::array<std::array<int, 2>, 4> quadrilateral = {{{0, 1}, {1, 2}, {2, 3}, {3, 0}}};
    static constexpr std::array<std::array<int, 2>, 6> tetrahedron = {{{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}}};
    static constexpr std::array<std::array<int, 2>, 12> hexahedron = {{
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4},
        {0, 4}, {1, 5}, {3, 7}, {2, 6}
    }};
    static constexpr std::array<std::array<int, 2>, 9> prism = {{
        {0, 1}, {1, 2}, {2, 0},
        {3, 4}, {4, 5}, {5, 3},
        {0, 3}, {1, 4}, {2, 5}
    }};
    static constexpr std::array<std::array<int, 2>, 8> pyramid = {{
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {0, 4}, {1, 4}, {2, 4}, {3, 4}
    }};

    switch (ct)
    {
    case cell::type::interval: return std::span(interval);
    case cell::type::triangle: return std::span(triangle);
    case cell::type::quadrilateral: return std::span(quadrilateral);
    case cell::type::tetrahedron: return std::span(tetrahedron);
    case cell::type::hexahedron: return std::span(hexahedron);
    case cell::type::prism: return std::span(prism);
    case cell::type::pyramid: return std::span(pyramid);
    default:
        throw std::invalid_argument("cut_edge_order: unsupported cell type");
    }
}

template <std::floating_point T>
int local_edge_for_cut_edge_token(const LocalMesh<T>& mesh,
                                  const std::span<const int32_t> parent_cell_vertices,
                                  const std::span<const int32_t> parent_cell_edges,
                                  const cell::type ct,
                                  const int token_edge)
{
    const auto cut_edges = cut_edge_order(ct);
    if (token_edge < 0 || token_edge >= static_cast<int>(cut_edges.size()))
        return -1;

    const int lv0 = cut_edges[static_cast<std::size_t>(token_edge)][0];
    const int lv1 = cut_edges[static_cast<std::size_t>(token_edge)][1];
    if (lv0 < 0 || lv1 < 0
        || lv0 >= static_cast<int>(parent_cell_vertices.size())
        || lv1 >= static_cast<int>(parent_cell_vertices.size()))
        return -1;

    const int32_t gv0 = parent_cell_vertices[static_cast<std::size_t>(lv0)];
    const int32_t gv1 = parent_cell_vertices[static_cast<std::size_t>(lv1)];

    for (const int32_t local_mesh_e : parent_cell_edges)
    {
        if (local_mesh_e < 0 || local_mesh_e >= mesh.n_edges())
            continue;
        const int32_t e0 = mesh.edge_vertices[static_cast<std::size_t>(2 * local_mesh_e)];
        const int32_t e1 = mesh.edge_vertices[static_cast<std::size_t>(2 * local_mesh_e + 1)];
        if ((e0 == gv0 && e1 == gv1) || (e0 == gv1 && e1 == gv0))
            return local_mesh_e;
    }

    return -1;
}

template <std::floating_point T>
void deduplicate_vertices_exact(LocalMesh<T>& mesh)
{
    const int old_nv = mesh.n_vertices();
    if (old_nv <= 1)
        return;

    const bool has_ref = !mesh.vertex_ref_x.empty();
    const int tdim = mesh.tdim;
    const int gdim = mesh.gdim;

    std::vector<int32_t> old_to_new(static_cast<std::size_t>(old_nv), -1);
    std::vector<T> new_x;
    std::vector<T> new_ref;
    std::vector<int32_t> new_pdim;
    std::vector<int32_t> new_pid;
    std::vector<int32_t> new_root_edge_id;
    new_x.reserve(mesh.vertex_x.size());
    if (has_ref)
        new_ref.reserve(mesh.vertex_ref_x.size());
    new_pdim.reserve(mesh.vertex_parent_dim.size());
    new_pid.reserve(mesh.vertex_parent_id.size());
    new_root_edge_id.reserve(mesh.vertex_root_edge_id.size());

    std::unordered_map<VertexDedupKey<T>, int32_t, VertexDedupKeyHash<T>> key_to_new;
    key_to_new.reserve(static_cast<std::size_t>(old_nv * 2));

    auto make_key = [&](int old_v) -> VertexDedupKey<T>
    {
        VertexDedupKey<T> key;
        key.parent_dim = mesh.vertex_parent_dim[static_cast<std::size_t>(old_v)];
        key.parent_id = mesh.vertex_parent_id[static_cast<std::size_t>(old_v)];
        key.gdim = static_cast<uint8_t>(gdim);
        key.tdim = static_cast<uint8_t>(tdim);
        key.has_ref = static_cast<uint8_t>(has_ref ? 1 : 0);
        for (int d = 0; d < gdim; ++d)
            key.x[static_cast<std::size_t>(d)] = bits_key(
                mesh.vertex_x[static_cast<std::size_t>(old_v * gdim + d)]);
        if (has_ref)
        {
            for (int d = 0; d < tdim; ++d)
                key.xref[static_cast<std::size_t>(d)] = bits_key(
                    mesh.vertex_ref_x[static_cast<std::size_t>(old_v * tdim + d)]);
        }
        return key;
    };

    int new_nv = 0;
    for (int old_v = 0; old_v < old_nv; ++old_v)
    {
        const auto key = make_key(old_v);
        const auto it = key_to_new.find(key);
        int found = -1;
        if (it == key_to_new.end())
        {
            found = new_nv++;
            key_to_new.emplace(key, found);
            for (int d = 0; d < gdim; ++d)
                new_x.push_back(mesh.vertex_x[static_cast<std::size_t>(old_v * gdim + d)]);
            if (has_ref)
            {
                for (int d = 0; d < tdim; ++d)
                    new_ref.push_back(mesh.vertex_ref_x[static_cast<std::size_t>(old_v * tdim + d)]);
            }
            new_pdim.push_back(mesh.vertex_parent_dim[static_cast<std::size_t>(old_v)]);
            new_pid.push_back(mesh.vertex_parent_id[static_cast<std::size_t>(old_v)]);
            new_root_edge_id.push_back(mesh.vertex_root_edge_id[static_cast<std::size_t>(old_v)]);
        }
        else
            found = it->second;
        old_to_new[static_cast<std::size_t>(old_v)] = found;
    }

    if (new_nv == old_nv)
        return;

    for (auto& cv : mesh.cell_vertices)
        cv = old_to_new[static_cast<std::size_t>(cv)];

    mesh.vertex_x.swap(new_x);
    mesh.vertex_ref_x.swap(new_ref);
    mesh.vertex_parent_dim.swap(new_pdim);
    mesh.vertex_parent_id.swap(new_pid);
    mesh.vertex_root_edge_id.swap(new_root_edge_id);
}

template <std::floating_point T>
void rebuild_parent_entity_maps(LocalMesh<T>& mesh)
{
    const int n_parent_vertices = cell::get_num_vertices(mesh.parent_cell_type);
    const int n_parent_edges = cell::num_edges(mesh.parent_cell_type);

    mesh.parent_vertex_to_local_vertex.assign(static_cast<std::size_t>(n_parent_vertices), -1);
    mesh.parent_edge_to_local_edge.assign(static_cast<std::size_t>(n_parent_edges), -1);

    for (int lv = 0; lv < mesh.n_vertices(); ++lv)
    {
        if (mesh.vertex_parent_dim[static_cast<std::size_t>(lv)] != 0)
            continue;
        const int pv = mesh.vertex_parent_id[static_cast<std::size_t>(lv)];
        if (pv >= 0 && pv < n_parent_vertices && mesh.parent_vertex_to_local_vertex[static_cast<std::size_t>(pv)] < 0)
            mesh.parent_vertex_to_local_vertex[static_cast<std::size_t>(pv)] = lv;
    }

    const auto parent_edges = cell::edges(mesh.parent_cell_type);

    // First attempt: match local edges that connect two parent corner vertices.
    for (int le = 0; le < mesh.n_edges(); ++le)
    {
        const int lv0 = mesh.edge_vertices[static_cast<std::size_t>(2 * le)];
        const int lv1 = mesh.edge_vertices[static_cast<std::size_t>(2 * le + 1)];

        if (mesh.vertex_parent_dim[static_cast<std::size_t>(lv0)] != 0
            || mesh.vertex_parent_dim[static_cast<std::size_t>(lv1)] != 0)
            continue;

        const int pv0 = mesh.vertex_parent_id[static_cast<std::size_t>(lv0)];
        const int pv1 = mesh.vertex_parent_id[static_cast<std::size_t>(lv1)];
        for (int pe = 0; pe < static_cast<int>(parent_edges.size()); ++pe)
        {
            const int e0 = parent_edges[static_cast<std::size_t>(pe)][0];
            const int e1 = parent_edges[static_cast<std::size_t>(pe)][1];
            if ((pv0 == e0 && pv1 == e1) || (pv0 == e1 && pv1 == e0))
            {
                if (mesh.parent_edge_to_local_edge[static_cast<std::size_t>(pe)] < 0)
                    mesh.parent_edge_to_local_edge[static_cast<std::size_t>(pe)] = le;
                break;
            }
        }
    }

    // Rewrite vertex parent ids:
    // - vertices on original edges should point to local edge id (when known)
    // - cell-interior vertices point to parent cell id
    for (int lv = 0; lv < mesh.n_vertices(); ++lv)
    {
        const int pdim = mesh.vertex_parent_dim[static_cast<std::size_t>(lv)];
        int& pid = mesh.vertex_parent_id[static_cast<std::size_t>(lv)];
        if (pdim == 1)
        {
            if (pid >= 0 && pid < n_parent_edges)
            {
                const int local_e = mesh.parent_edge_to_local_edge[static_cast<std::size_t>(pid)];
                if (local_e >= 0)
                    pid = local_e;
            }
        }
        else if (pdim == mesh.tdim)
        {
            pid = mesh.parent_cell_id;
        }
    }
}

int infer_lagrange_order(cell::type ct, int n_vertices)
{
    if (n_vertices == cell::get_num_vertices(ct))
        return 1;

    for (int p = 2; p <= 4; ++p)
    {
        if (iso_p1_template(ct, p).n_vertices == n_vertices)
            return p;
    }
    throw std::invalid_argument("init_local_mesh_from_cell: could not infer interpolation order from node count");
}

std::pair<int, int> infer_gdim_and_order(cell::type ct, std::size_t flat_coord_size)
{
    std::vector<std::pair<int, int>> matches;
    for (int gdim = 1; gdim <= 3; ++gdim)
    {
        if (flat_coord_size % static_cast<std::size_t>(gdim) != 0)
            continue;
        const int n_vertices = static_cast<int>(flat_coord_size / static_cast<std::size_t>(gdim));
        try
        {
            const int order = infer_lagrange_order(ct, n_vertices);
            matches.push_back({gdim, order});
        }
        catch (const std::invalid_argument&)
        {
            // Try next candidate.
        }
    }

    if (matches.empty())
        throw std::invalid_argument("init_local_mesh_from_cell: could not infer geometric dimension and interpolation order");
    if (matches.size() > 1)
        throw std::invalid_argument("init_local_mesh_from_cell: ambiguous geometric dimension/order from coordinate array");
    return matches[0];
}
} // namespace

// ============================================================================
// LocalMesh Implementation
// ============================================================================

template <std::floating_point T>
void init_local_mesh_from_template(
    LocalMesh<T>&             mesh,
    const RefinementTemplate& tpl,
    std::span<const T>        parent_cell_coords,
    cell::type                parent_cell_type,
    int                       parent_cell_id,
    int                       n_level_sets)
{
    if (parent_cell_type != tpl.bg_cell_type)
        throw std::invalid_argument("init_local_mesh_from_template: parent_cell_type does not match template background cell type");

    const int nv = tpl.n_vertices;
    const int nc = tpl.n_cells;
    const int tdim = tpl.tdim;
    if (nv <= 0)
        throw std::invalid_argument("init_local_mesh_from_template: template has no vertices");
    if (parent_cell_coords.size() % static_cast<std::size_t>(nv) != 0)
        throw std::invalid_argument("init_local_mesh_from_template: parent_cell_coords size is not divisible by n_vertices");

    const int gdim = static_cast<int>(parent_cell_coords.size() / static_cast<std::size_t>(nv));
    if (gdim <= 0)
        throw std::invalid_argument("init_local_mesh_from_template: invalid geometric dimension");

    mesh.gdim = gdim;
    mesh.tdim = tdim;
    mesh.parent_cell_id = parent_cell_id;
    mesh.parent_cell_type = parent_cell_type;
    mesh.n_level_sets = n_level_sets;
    const int p1_nv = cell::get_num_vertices(parent_cell_type);
    mesh.parent_cell_coords_p1.resize(static_cast<std::size_t>(p1_nv * gdim));
    for (int i = 0; i < p1_nv; ++i)
    {
        for (int d = 0; d < gdim; ++d)
            mesh.parent_cell_coords_p1[static_cast<std::size_t>(i * gdim + d)] = parent_cell_coords[static_cast<std::size_t>(i * gdim + d)];
    }

    // Reset and resize vertices
    mesh.vertex_ref_x.resize(static_cast<std::size_t>(nv * tdim));
    mesh.vertex_x.resize(nv * gdim);
    mesh.vertex_parent_dim.assign(tpl.vertex_parent_dim.begin(), tpl.vertex_parent_dim.end());
    mesh.vertex_parent_id.assign(tpl.vertex_parent_id.begin(), tpl.vertex_parent_id.end());
    mesh.vertex_root_edge_id.assign(static_cast<std::size_t>(nv), -1);

    for (int i = 0; i < nv; ++i)
    {
        for (int d = 0; d < tdim; ++d)
            mesh.vertex_ref_x[static_cast<std::size_t>(i * tdim + d)] = static_cast<T>(
                tpl.ref_vertex_coords[static_cast<std::size_t>(i * tdim + d)]);
        for (int d = 0; d < gdim; ++d)
            mesh.vertex_x[i * gdim + d] = parent_cell_coords[i * gdim + d];
    }

    // Build cells
    mesh.cell_offsets.resize(nc + 1);
    mesh.cell_vertices.assign(tpl.cell_connectivity.begin(), tpl.cell_connectivity.end());
    for (int i = 0; i <= nc; ++i) {
        mesh.cell_offsets[i] = i * tpl.vertices_per_cell;
    }
    mesh.cell_types.assign(nc, tpl.child_cell_type);
    mesh.cell_domain.assign(nc, static_cast<uint8_t>(cell::domain::unset));

    deduplicate_vertices_exact(mesh);

    // Classification masks
    const int nvd = mesh.n_vertices();
    mesh.vertex_zero_mask.assign(nvd, 0);
    mesh.vertex_inside_mask.assign(nvd, 0);
    mesh.vertex_phi.assign(static_cast<std::size_t>(nvd * n_level_sets), T(0));

    build_local_edges(mesh);
    rebuild_parent_entity_maps(mesh);
}

template <std::floating_point T>
void init_local_mesh_from_cell(
    LocalMesh<T>&      mesh,
    std::span<const T> parent_cell_coords,
    cell::type         parent_cell_type,
    int                parent_cell_id,
    int                n_level_sets)
{
    const auto [gdim, order] = infer_gdim_and_order(parent_cell_type, parent_cell_coords.size());

    const RefinementTemplate& meta_tpl = (order == 1) ? p1_template(parent_cell_type) : iso_p1_template(parent_cell_type, order);
    const RefinementTemplate& p1_tpl = p1_template(parent_cell_type);

    mesh.gdim = gdim;
    mesh.tdim = meta_tpl.tdim;
    mesh.parent_cell_id = parent_cell_id;
    mesh.parent_cell_type = parent_cell_type;
    mesh.n_level_sets = n_level_sets;

    mesh.vertex_x.assign(parent_cell_coords.begin(), parent_cell_coords.end());
    const int p1_nv = cell::get_num_vertices(parent_cell_type);
    mesh.parent_cell_coords_p1.resize(static_cast<std::size_t>(p1_nv * gdim));
    for (int i = 0; i < p1_nv; ++i)
    {
        for (int d = 0; d < gdim; ++d)
            mesh.parent_cell_coords_p1[static_cast<std::size_t>(i * gdim + d)] = parent_cell_coords[static_cast<std::size_t>(i * gdim + d)];
    }
    mesh.vertex_ref_x.resize(static_cast<std::size_t>(meta_tpl.n_vertices * meta_tpl.tdim));
    if (order == 1)
    {
        const auto xref = p1_ref_coords(parent_cell_type);
        for (std::size_t i = 0; i < xref.size(); ++i)
            mesh.vertex_ref_x[i] = static_cast<T>(xref[i]);
    }
    else
    {
        const auto xref = iso_p1_ref_coords(parent_cell_type, order);
        for (std::size_t i = 0; i < xref.size(); ++i)
            mesh.vertex_ref_x[i] = static_cast<T>(xref[i]);
    }

    mesh.vertex_parent_dim.assign(meta_tpl.vertex_parent_dim.begin(), meta_tpl.vertex_parent_dim.end());
    mesh.vertex_parent_id.assign(meta_tpl.vertex_parent_id.begin(), meta_tpl.vertex_parent_id.end());
    mesh.vertex_root_edge_id.assign(static_cast<std::size_t>(meta_tpl.n_vertices), -1);
    mesh.cell_offsets = {0, p1_tpl.vertices_per_cell};
    mesh.cell_vertices.assign(p1_tpl.cell_connectivity.begin(), p1_tpl.cell_connectivity.end());
    mesh.cell_types = {parent_cell_type};
    mesh.cell_domain = {static_cast<uint8_t>(cell::domain::unset)};

    deduplicate_vertices_exact(mesh);

    const int nvd = mesh.n_vertices();
    mesh.vertex_zero_mask.assign(nvd, 0);
    mesh.vertex_inside_mask.assign(nvd, 0);
    mesh.vertex_phi.assign(static_cast<std::size_t>(nvd * n_level_sets), T(0));

    build_local_edges(mesh);
    rebuild_parent_entity_maps(mesh);
}

template <std::floating_point T>
void refine_local_mesh_from_template(
    LocalMesh<T>&             mesh,
    const RefinementTemplate& tpl)
{
    const cell::type parent_cell_type = mesh.parent_cell_type;
    const int parent_cell_id = mesh.parent_cell_id;
    const int n_level_sets = mesh.n_level_sets;
    const int p1_nv = cell::get_num_vertices(parent_cell_type);
    if (p1_nv <= 0)
        throw std::invalid_argument("refine_local_mesh_from_template: invalid parent cell type");

    // If current local-mesh vertices already match template-node count, reuse them
    // directly as parent-cell interpolation-node coordinates.
    if (mesh.gdim > 0
        && static_cast<int>(mesh.vertex_x.size()) == tpl.n_vertices * mesh.gdim)
    {
        const std::vector<T> current_coords = mesh.vertex_x;
        init_local_mesh_from_template(
            mesh,
            tpl,
            std::span<const T>(current_coords.data(), current_coords.size()),
            parent_cell_type,
            parent_cell_id,
            n_level_sets);
        return;
    }

    const std::span<const T> parent_cell_coords(mesh.parent_cell_coords_p1.data(), mesh.parent_cell_coords_p1.size());
    int gdim = -1;
    if (parent_cell_coords.size() % static_cast<std::size_t>(p1_nv) == 0)
        gdim = static_cast<int>(parent_cell_coords.size() / static_cast<std::size_t>(p1_nv));
    if (gdim <= 0)
        throw std::invalid_argument("refine_local_mesh_from_template: missing valid stored parent-cell P1 coordinates");

    std::vector<T> mapped_template_coords(static_cast<std::size_t>(tpl.n_vertices * gdim), T(0));
    std::vector<T> parent_p1(parent_cell_coords.begin(), parent_cell_coords.end());
    const std::vector<T> ref_coords_t(
        tpl.ref_vertex_coords.begin(), tpl.ref_vertex_coords.end());

    cell::push_forward_affine(
        parent_cell_type, parent_p1, gdim,
        std::span<const T>(ref_coords_t.data(), ref_coords_t.size()),
        std::span<T>(mapped_template_coords.data(), mapped_template_coords.size()));

    init_local_mesh_from_template(
        mesh,
        tpl,
        std::span<const T>(mapped_template_coords.data(), mapped_template_coords.size()),
        parent_cell_type,
        parent_cell_id,
        n_level_sets);
    rebuild_parent_entity_maps(mesh);
}

template <std::floating_point T>
void evaluate_levelsets_on_vertices(
    LocalMesh<T>&                           mesh,
    const std::vector<LevelSetFunction<T>>& phi,
    int                                     n_ls,
    T                                       tol)
{
    int nv = mesh.n_vertices();
    int gdim = mesh.gdim;
    mesh.vertex_phi.assign(nv * n_ls, 0.0);
    mesh.vertex_zero_mask.assign(nv, 0);
    mesh.vertex_inside_mask.assign(nv, 0);

    for (int i = 0; i < nv; ++i) {
        const T* x = &mesh.vertex_x[i * gdim];
        for (int ls = 0; ls < n_ls; ++ls) {
            T val = 0.0;
            val = phi[ls].value(x, mesh.parent_cell_id);
            
            mesh.vertex_phi[i * n_ls + ls] = val;
            if (std::abs(val) < tol) {
                mesh.vertex_zero_mask[i] |= (1ULL << ls);
            } else if (val < 0) {
                mesh.vertex_inside_mask[i] |= (1ULL << ls);
            }
        }
    }
}

template <std::floating_point T>
void classify_local_edges(LocalMesh<T>& mesh, int level_set_id)
{
    int ne = mesh.n_edges();
    uint64_t mask = (1ULL << level_set_id);

    for (int i = 0; i < ne; ++i) {
        int v0 = mesh.edge_vertices[i * 2];
        int v1 = mesh.edge_vertices[i * 2 + 1];

        bool zero0 = (mesh.vertex_zero_mask[v0] & mask);
        bool zero1 = (mesh.vertex_zero_mask[v1] & mask);
        bool inside0 = (mesh.vertex_inside_mask[v0] & mask);
        bool inside1 = (mesh.vertex_inside_mask[v1] & mask);

        if (zero0 || zero1) {
            mesh.edge_state[i] = static_cast<uint8_t>(EdgeState::near_tangency);
        } else if (inside0 == inside1) {
            mesh.edge_state[i] = static_cast<uint8_t>(EdgeState::no_root);
        } else {
            mesh.edge_state[i] = static_cast<uint8_t>(EdgeState::one_root);
        }
    }
}

template <std::floating_point T>
int compute_edge_root_linear(LocalMesh<T>& mesh, int edge_id, int level_set_id)
{
    if (edge_id < 0 || edge_id >= mesh.n_edges())
        throw std::invalid_argument("compute_edge_root_linear: invalid edge id");
    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("compute_edge_root_linear: invalid level_set_id");

    if (mesh.edge_root_vertex[static_cast<std::size_t>(edge_id)] >= 0)
        return mesh.edge_root_vertex[static_cast<std::size_t>(edge_id)];

    int v0_idx = mesh.edge_vertices[edge_id * 2];
    int v1_idx = mesh.edge_vertices[edge_id * 2 + 1];

    int n_ls = mesh.n_level_sets;
    T val0 = mesh.vertex_phi[v0_idx * n_ls + level_set_id];
    T val1 = mesh.vertex_phi[v1_idx * n_ls + level_set_id];

    // Linear interpolation factor s in [0, 1]
    const T denom = (val1 - val0);
    T s = T(0.5);
    if (std::abs(denom) > std::numeric_limits<T>::epsilon())
        s = -val0 / denom;
    s = std::clamp(s, T(0), T(1));

    // Create new vertex
    int new_v_idx = mesh.n_vertices();
    int gdim = mesh.gdim;
    int tdim = mesh.tdim;

    for (int d = 0; d < gdim; ++d) {
        mesh.vertex_x.push_back(mesh.vertex_x[v0_idx * gdim + d] + s * (mesh.vertex_x[v1_idx * gdim + d] - mesh.vertex_x[v0_idx * gdim + d]));
    }
    for (int d = 0; d < tdim; ++d) {
        mesh.vertex_ref_x.push_back(mesh.vertex_ref_x[v0_idx * tdim + d] + s * (mesh.vertex_ref_x[v1_idx * tdim + d] - mesh.vertex_ref_x[v0_idx * tdim + d]));
    }

    // Parent info for the new root
    mesh.vertex_parent_dim.push_back(1); // Edge
    mesh.vertex_parent_id.push_back(edge_id);
    mesh.vertex_root_edge_id.push_back(edge_id);

    // Level set values at root are approx 0
    mesh.vertex_phi.resize(mesh.vertex_phi.size() + n_ls, 0.0);
    mesh.vertex_phi[(new_v_idx * n_ls) + level_set_id] = 0.0;
    mesh.vertex_zero_mask.push_back(1ULL << level_set_id);
    mesh.vertex_inside_mask.push_back(0);

    mesh.edge_root_vertex[edge_id] = new_v_idx;
    mesh.edge_root_parameter[static_cast<std::size_t>(edge_id)] = s;
    mesh.edge_root_iterations[static_cast<std::size_t>(edge_id)] = 0;
    mesh.edge_root_evaluations[static_cast<std::size_t>(edge_id)] = 0;
    mesh.edge_root_converged[static_cast<std::size_t>(edge_id)] = 1;
    mesh.edge_root_residual[static_cast<std::size_t>(edge_id)] = std::abs(val0 + s * (val1 - val0));
    return new_v_idx;
}

template <std::floating_point T, std::integral I>
int compute_edge_root(LocalMesh<T>& mesh,
                      const LevelSetFunction<T, I>& level_set,
                      int edge_id,
                      int level_set_id,
                      cell::edge_root::method root_method)
{
    if (root_method == cell::edge_root::method::linear || !level_set.has_value())
        return compute_edge_root_linear(mesh, edge_id, level_set_id);

    if (edge_id < 0 || edge_id >= mesh.n_edges())
        throw std::invalid_argument("compute_edge_root: invalid edge id");
    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("compute_edge_root: invalid level_set_id");
    if (mesh.edge_root_vertex[static_cast<std::size_t>(edge_id)] >= 0)
        return mesh.edge_root_vertex[static_cast<std::size_t>(edge_id)];

    const int v0_idx = mesh.edge_vertices[static_cast<std::size_t>(edge_id * 2)];
    const int v1_idx = mesh.edge_vertices[static_cast<std::size_t>(edge_id * 2 + 1)];
    const int gdim = mesh.gdim;
    const int tdim = mesh.tdim;
    const int n_ls = mesh.n_level_sets;

    const std::span<const T> p0(mesh.vertex_x.data() + static_cast<std::size_t>(v0_idx * gdim),
                                static_cast<std::size_t>(gdim));
    const std::span<const T> p1(mesh.vertex_x.data() + static_cast<std::size_t>(v1_idx * gdim),
                                static_cast<std::size_t>(gdim));
    auto phi = [&](std::span<const T> x) -> T
    {
        return level_set.value(x.data(), static_cast<I>(mesh.parent_cell_id));
    };

    const auto root_info = cell::edge_root::find_root_parameter_info<T>(
        p0, p1, phi, root_method, T(0), 64, static_cast<T>(1e-12), static_cast<T>(1e-12));
    const T t = root_info.t;

    const int new_v_idx = mesh.n_vertices();
    for (int d = 0; d < gdim; ++d)
    {
        mesh.vertex_x.push_back(mesh.vertex_x[static_cast<std::size_t>(v0_idx * gdim + d)]
                                + t * (mesh.vertex_x[static_cast<std::size_t>(v1_idx * gdim + d)]
                                       - mesh.vertex_x[static_cast<std::size_t>(v0_idx * gdim + d)]));
    }
    for (int d = 0; d < tdim; ++d)
    {
        mesh.vertex_ref_x.push_back(mesh.vertex_ref_x[static_cast<std::size_t>(v0_idx * tdim + d)]
                                    + t * (mesh.vertex_ref_x[static_cast<std::size_t>(v1_idx * tdim + d)]
                                           - mesh.vertex_ref_x[static_cast<std::size_t>(v0_idx * tdim + d)]));
    }

    mesh.vertex_parent_dim.push_back(1);
    mesh.vertex_parent_id.push_back(edge_id);
    mesh.vertex_root_edge_id.push_back(edge_id);

    mesh.vertex_phi.resize(mesh.vertex_phi.size() + static_cast<std::size_t>(n_ls), T(0));
    mesh.vertex_zero_mask.push_back(0);
    mesh.vertex_inside_mask.push_back(0);
    const T val = phi(std::span<const T>(
        mesh.vertex_x.data() + static_cast<std::size_t>(new_v_idx * gdim),
        static_cast<std::size_t>(gdim)));
    mesh.vertex_phi[static_cast<std::size_t>(new_v_idx * n_ls + level_set_id)] = val;
    const uint64_t mask = (uint64_t(1) << level_set_id);
    const T tol = static_cast<T>(1e-12);
    if (std::abs(val) <= tol)
        mesh.vertex_zero_mask.back() |= mask;
    else if (val < T(0))
        mesh.vertex_inside_mask.back() |= mask;

    mesh.edge_root_vertex[static_cast<std::size_t>(edge_id)] = new_v_idx;
    mesh.edge_root_parameter[static_cast<std::size_t>(edge_id)] = t;
    mesh.edge_root_iterations[static_cast<std::size_t>(edge_id)] = root_info.iterations;
    mesh.edge_root_evaluations[static_cast<std::size_t>(edge_id)] = root_info.evaluations;
    mesh.edge_root_converged[static_cast<std::size_t>(edge_id)] = root_info.converged ? 1 : 0;
    mesh.edge_root_residual[static_cast<std::size_t>(edge_id)] = root_info.residual;
    return new_v_idx;
}

template <std::floating_point T>
void compute_all_roots_linear(LocalMesh<T>& mesh, int level_set_id)
{
    int ne = mesh.n_edges();
    for (int i = 0; i < ne; ++i) {
        if (mesh.edge_state[i] == static_cast<uint8_t>(EdgeState::one_root) && mesh.edge_root_vertex[i] == -1) {
            compute_edge_root_linear(mesh, i, level_set_id);
        }
    }
}

template <std::floating_point T, std::integral I>
void compute_all_roots(LocalMesh<T>& mesh,
                       const LevelSetFunction<T, I>& level_set,
                       int level_set_id,
                       cell::edge_root::method root_method)
{
    if (root_method == cell::edge_root::method::linear)
    {
        compute_all_roots_linear(mesh, level_set_id);
        return;
    }

    const int ne = mesh.n_edges();
    for (int i = 0; i < ne; ++i)
    {
        if (mesh.edge_state[static_cast<std::size_t>(i)] != static_cast<uint8_t>(EdgeState::one_root))
            continue;
        if (mesh.edge_root_vertex[static_cast<std::size_t>(i)] >= 0)
            continue;
        compute_edge_root<T, I>(mesh, level_set, i, level_set_id, root_method);
    }
}

namespace
{
template <std::floating_point T>
int find_vertex_by_coords(const LocalMesh<T>& mesh, const std::span<const T> x, T tol = static_cast<T>(1e-14))
{
    const int nv = mesh.n_vertices();
    const int gdim = mesh.gdim;
    for (int v = 0; v < nv; ++v)
    {
        bool same = true;
        for (int d = 0; d < gdim; ++d)
        {
            const T dv = mesh.vertex_x[static_cast<std::size_t>(v * gdim + d)] - x[static_cast<std::size_t>(d)];
            if (std::abs(dv) > tol)
            {
                same = false;
                break;
            }
        }
        if (same)
            return v;
    }
    return -1;
}

template <std::floating_point T>
int map_cut_vertex_to_local_mesh(
    LocalMesh<T>&               mesh,
    const cell::CutCell<T>&     cut_cell,
    const int                   lv,
    const cell::type            parent_cell_type,
    const std::span<const int32_t> parent_cell_vertices,
    const std::span<const int32_t> parent_cell_edges,
    const int                   level_set_id)
{
    const int gdim = mesh.gdim;
    const int tdim = mesh.tdim;
    const int nls = mesh.n_level_sets;
    const int ncutv = static_cast<int>(cut_cell._vertex_coords.size()) / gdim;

    if (lv < 0 || lv >= ncutv)
        throw std::invalid_argument("map_cut_vertex_to_local_mesh: invalid local cut-vertex id");

    const bool has_tokens = static_cast<int>(cut_cell._vertex_parent_entity.size()) == ncutv;
    const int32_t token = has_tokens ? cut_cell._vertex_parent_entity[static_cast<std::size_t>(lv)] : -1;

    if (token >= 100 && token < 200)
    {
        const int local_parent_v = static_cast<int>(token - 100);
        if (local_parent_v >= 0 && local_parent_v < static_cast<int>(parent_cell_vertices.size()))
            return parent_cell_vertices[static_cast<std::size_t>(local_parent_v)];
    }
    else if (token >= 0 && token < 100)
    {
        const int local_parent_e = static_cast<int>(token);
        const int local_mesh_e = local_edge_for_cut_edge_token(
            mesh, parent_cell_vertices, parent_cell_edges, parent_cell_type, local_parent_e);
        if (local_mesh_e >= 0 && local_mesh_e < mesh.n_edges())
        {
            int rv = mesh.edge_root_vertex[static_cast<std::size_t>(local_mesh_e)];
            if (rv < 0
                && mesh.edge_state[static_cast<std::size_t>(local_mesh_e)] == static_cast<uint8_t>(EdgeState::one_root))
            {
                rv = compute_edge_root_linear(mesh, local_mesh_e, level_set_id);
            }
            if (rv >= 0)
                return rv;
        }
    }

    // Fallback (rare): lookup by coordinates and add if not found.
    const std::span<const T> xv(
        cut_cell._vertex_coords.data() + static_cast<std::size_t>(lv * gdim),
        static_cast<std::size_t>(gdim));
    const int existing = find_vertex_by_coords(mesh, xv);
    if (existing >= 0)
        return existing;

    const int new_v = mesh.n_vertices();
    for (int d = 0; d < gdim; ++d)
        mesh.vertex_x.push_back(xv[static_cast<std::size_t>(d)]);

    // Keep reference coordinates sized and coherent. If unavailable, append zeros.
    if (!mesh.vertex_ref_x.empty())
    {
        for (int d = 0; d < tdim; ++d)
            mesh.vertex_ref_x.push_back(T(0));
    }

    int32_t parent_dim = -1;
    int32_t parent_id = -1;
    int32_t root_edge_id = -1;
    if (token >= 100 && token < 200)
    {
        parent_dim = 0;
        const int local_parent_v = static_cast<int>(token - 100);
        if (local_parent_v >= 0 && local_parent_v < static_cast<int>(parent_cell_vertices.size()))
            parent_id = parent_cell_vertices[static_cast<std::size_t>(local_parent_v)];
    }
    else if (token >= 0 && token < 100)
    {
        parent_dim = 1;
        const int local_parent_e = static_cast<int>(token);
        parent_id = local_edge_for_cut_edge_token(
            mesh, parent_cell_vertices, parent_cell_edges, parent_cell_type, local_parent_e);
        root_edge_id = parent_id;
    }
    mesh.vertex_parent_dim.push_back(parent_dim);
    mesh.vertex_parent_id.push_back(parent_id);
    mesh.vertex_root_edge_id.push_back(root_edge_id);

    mesh.vertex_phi.resize(mesh.vertex_phi.size() + static_cast<std::size_t>(nls), T(0));
    mesh.vertex_zero_mask.push_back(0);
    mesh.vertex_inside_mask.push_back(0);
    if (token >= 0 && token < 100 && level_set_id >= 0 && level_set_id < nls)
    {
        mesh.vertex_zero_mask.back() = (uint64_t(1) << level_set_id);
        mesh.vertex_phi[static_cast<std::size_t>(new_v * nls + level_set_id)] = T(0);
    }
    return new_v;
}

template <std::floating_point T>
void append_cut_fragments(
    LocalMesh<T>&                     mesh,
    const cell::CutCell<T>&           cut_cell,
    const cell::type                  parent_cell_type,
    const std::span<const int32_t>    parent_cell_vertices,
    const std::span<const int32_t>    parent_cell_edges,
    const cell::domain                dom,
    const int                         level_set_id,
    std::vector<int32_t>&             out_cell_vertices,
    std::vector<int32_t>&             out_cell_offsets,
    std::vector<cell::type>&          out_cell_types,
    std::vector<uint8_t>&             out_cell_domain)
{
    const int ncutv = static_cast<int>(cut_cell._vertex_coords.size()) / mesh.gdim;
    if (ncutv == 0 || cut_cell._offset.empty())
        return;

    std::vector<int32_t> cut_to_local(static_cast<std::size_t>(ncutv), -1);
    for (int lv = 0; lv < ncutv; ++lv)
    {
        cut_to_local[static_cast<std::size_t>(lv)] = map_cut_vertex_to_local_mesh(
            mesh, cut_cell, lv, parent_cell_type, parent_cell_vertices, parent_cell_edges, level_set_id);
    }

    const int nsub = cell::num_cells(cut_cell);
    for (int i = 0; i < nsub; ++i)
    {
        const auto verts = cell::cell_vertices(cut_cell, i);
        out_cell_types.push_back(cut_cell._types[static_cast<std::size_t>(i)]);
        out_cell_domain.push_back(static_cast<uint8_t>(dom));
        for (const int lv : verts)
            out_cell_vertices.push_back(cut_to_local[static_cast<std::size_t>(lv)]);
        out_cell_offsets.push_back(static_cast<int32_t>(out_cell_vertices.size()));
    }
}
} // namespace

template <std::floating_point T>
void decompose_local_mesh_linear(LocalMesh<T>& mesh, int level_set_id, bool triangulate)
{
    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("decompose_local_mesh_linear: invalid level_set_id");
    if (static_cast<int>(mesh.vertex_root_edge_id.size()) != mesh.n_vertices())
        mesh.vertex_root_edge_id.assign(static_cast<std::size_t>(mesh.n_vertices()), -1);

    compute_all_roots_linear(mesh, level_set_id);

    const std::vector<int32_t> old_cell_vertices = mesh.cell_vertices;
    const std::vector<int32_t> old_cell_offsets = mesh.cell_offsets;
    const std::vector<cell::type> old_cell_types = mesh.cell_types;
    const std::vector<int32_t> old_cell_edge_offsets = mesh.cell_edge_offsets;
    const std::vector<int32_t> old_cell_edges_flat = mesh.cell_edges_flat;

    std::vector<int32_t> new_cell_vertices;
    std::vector<int32_t> new_cell_offsets(1, 0);
    std::vector<cell::type> new_cell_types;
    std::vector<uint8_t> new_cell_domain;

    const int old_nc = static_cast<int>(old_cell_types.size());
    for (int c = 0; c < old_nc; ++c)
    {
        const int c0 = old_cell_offsets[static_cast<std::size_t>(c)];
        const int c1 = old_cell_offsets[static_cast<std::size_t>(c + 1)];
        const int nv = c1 - c0;
        const cell::type ct = old_cell_types[static_cast<std::size_t>(c)];

        const std::span<const int32_t> parent_cell_vertices(
            old_cell_vertices.data() + static_cast<std::size_t>(c0),
            static_cast<std::size_t>(nv));

        const int ce0 = old_cell_edge_offsets[static_cast<std::size_t>(c)];
        const int ce1 = old_cell_edge_offsets[static_cast<std::size_t>(c + 1)];
        const std::span<const int32_t> parent_cell_edges(
            old_cell_edges_flat.data() + static_cast<std::size_t>(ce0),
            static_cast<std::size_t>(ce1 - ce0));

        std::vector<T> ls_values(static_cast<std::size_t>(nv), T(0));
        std::vector<T> cell_vertex_coords(static_cast<std::size_t>(nv * mesh.gdim), T(0));
        for (int j = 0; j < nv; ++j)
        {
            const int gv = parent_cell_vertices[static_cast<std::size_t>(j)];
            ls_values[static_cast<std::size_t>(j)] = mesh.vertex_phi[static_cast<std::size_t>(gv * mesh.n_level_sets + level_set_id)];
            for (int d = 0; d < mesh.gdim; ++d)
            {
                cell_vertex_coords[static_cast<std::size_t>(j * mesh.gdim + d)]
                    = mesh.vertex_x[static_cast<std::size_t>(gv * mesh.gdim + d)];
            }
        }

        const cell::domain dom = cell::classify_cell_domain<T>(
            std::span<const T>(ls_values.data(), ls_values.size()));
        if (dom == cell::domain::inside || dom == cell::domain::outside)
        {
            new_cell_types.push_back(ct);
            new_cell_domain.push_back(static_cast<uint8_t>(dom));
            for (int j = 0; j < nv; ++j)
                new_cell_vertices.push_back(parent_cell_vertices[static_cast<std::size_t>(j)]);
            new_cell_offsets.push_back(static_cast<int32_t>(new_cell_vertices.size()));
            continue;
        }

        const int n_cut_edges = cell::num_edges(ct);
        std::vector<T> edge_root_coords(static_cast<std::size_t>(n_cut_edges * mesh.gdim), T(0));
        std::vector<uint8_t> edge_has_root(static_cast<std::size_t>(n_cut_edges), 0);
        for (int te = 0; te < n_cut_edges; ++te)
        {
            const int local_mesh_e = local_edge_for_cut_edge_token(
                mesh, parent_cell_vertices, parent_cell_edges, ct, te);
            if (local_mesh_e < 0 || local_mesh_e >= mesh.n_edges())
                continue;
            const int rv = mesh.edge_root_vertex[static_cast<std::size_t>(local_mesh_e)];
            if (rv < 0 || rv >= mesh.n_vertices())
                continue;
            edge_has_root[static_cast<std::size_t>(te)] = 1;
            for (int d = 0; d < mesh.gdim; ++d)
            {
                edge_root_coords[static_cast<std::size_t>(te * mesh.gdim + d)]
                    = mesh.vertex_x[static_cast<std::size_t>(rv * mesh.gdim + d)];
            }
        }

        cell::CutCell<T> cut_inside;
        cell::CutCell<T> cut_outside;
        cell::cut_from_cached_roots<T>(
            ct,
            std::span<const T>(cell_vertex_coords.data(), cell_vertex_coords.size()),
            mesh.gdim,
            std::span<const T>(ls_values.data(), ls_values.size()),
            "phi<0",
            std::span<const T>(edge_root_coords.data(), edge_root_coords.size()),
            std::span<const uint8_t>(edge_has_root.data(), edge_has_root.size()),
            cut_inside,
            triangulate);
        cell::cut_from_cached_roots<T>(
            ct,
            std::span<const T>(cell_vertex_coords.data(), cell_vertex_coords.size()),
            mesh.gdim,
            std::span<const T>(ls_values.data(), ls_values.size()),
            "phi>0",
            std::span<const T>(edge_root_coords.data(), edge_root_coords.size()),
            std::span<const uint8_t>(edge_has_root.data(), edge_has_root.size()),
            cut_outside,
            triangulate);

        append_cut_fragments(
            mesh, cut_inside, ct, parent_cell_vertices, parent_cell_edges,
            cell::domain::inside, level_set_id,
            new_cell_vertices, new_cell_offsets, new_cell_types, new_cell_domain);
        append_cut_fragments(
            mesh, cut_outside, ct, parent_cell_vertices, parent_cell_edges,
            cell::domain::outside, level_set_id,
            new_cell_vertices, new_cell_offsets, new_cell_types, new_cell_domain);
    }

    mesh.cell_vertices.swap(new_cell_vertices);
    mesh.cell_offsets.swap(new_cell_offsets);
    mesh.cell_types.swap(new_cell_types);
    mesh.cell_domain.swap(new_cell_domain);

    build_local_edges(mesh);
    rebuild_parent_entity_maps(mesh);
}

template <std::floating_point T, std::integral I>
void decompose_local_mesh(LocalMesh<T>& mesh,
                          const LevelSetFunction<T, I>& level_set,
                          int level_set_id,
                          cell::edge_root::method root_method,
                          bool triangulate)
{
    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("decompose_local_mesh: invalid level_set_id");
    if (static_cast<int>(mesh.vertex_root_edge_id.size()) != mesh.n_vertices())
        mesh.vertex_root_edge_id.assign(static_cast<std::size_t>(mesh.n_vertices()), -1);

    compute_all_roots<T, I>(mesh, level_set, level_set_id, root_method);
    decompose_local_mesh_linear(mesh, level_set_id, triangulate);
}

// Explicit instantiations for common types
template void init_local_mesh_from_template<double>(LocalMesh<double>&, const RefinementTemplate&, std::span<const double>, cell::type, int, int);
template void init_local_mesh_from_template<float>(LocalMesh<float>&, const RefinementTemplate&, std::span<const float>, cell::type, int, int);
template void init_local_mesh_from_cell<double>(LocalMesh<double>&, std::span<const double>, cell::type, int, int);
template void init_local_mesh_from_cell<float>(LocalMesh<float>&, std::span<const float>, cell::type, int, int);
template void refine_local_mesh_from_template<double>(LocalMesh<double>&, const RefinementTemplate&);
template void refine_local_mesh_from_template<float>(LocalMesh<float>&, const RefinementTemplate&);
template void evaluate_levelsets_on_vertices<double>(LocalMesh<double>&, const std::vector<LevelSetFunction<double>>&, int, double);
template void evaluate_levelsets_on_vertices<float>(LocalMesh<float>&, const std::vector<LevelSetFunction<float>>&, int, float);
template void classify_local_edges<double>(LocalMesh<double>&, int);
template void classify_local_edges<float>(LocalMesh<float>&, int);
template int compute_edge_root_linear<double>(LocalMesh<double>&, int, int);
template int compute_edge_root_linear<float>(LocalMesh<float>&, int, int);
template int compute_edge_root<double, int>(LocalMesh<double>&, const LevelSetFunction<double, int>&, int, int, cell::edge_root::method);
template int compute_edge_root<float, int>(LocalMesh<float>&, const LevelSetFunction<float, int>&, int, int, cell::edge_root::method);
template void compute_all_roots_linear<double>(LocalMesh<double>&, int);
template void compute_all_roots_linear<float>(LocalMesh<float>&, int);
template void compute_all_roots<double, int>(LocalMesh<double>&, const LevelSetFunction<double, int>&, int, cell::edge_root::method);
template void compute_all_roots<float, int>(LocalMesh<float>&, const LevelSetFunction<float, int>&, int, cell::edge_root::method);
template void decompose_local_mesh_linear<double>(LocalMesh<double>&, int, bool);
template void decompose_local_mesh_linear<float>(LocalMesh<float>&, int, bool);
template void decompose_local_mesh<double, int>(LocalMesh<double>&, const LevelSetFunction<double, int>&, int, cell::edge_root::method, bool);
template void decompose_local_mesh<float, int>(LocalMesh<float>&, const LevelSetFunction<float, int>&, int, cell::edge_root::method, bool);

} // namespace cutcells
