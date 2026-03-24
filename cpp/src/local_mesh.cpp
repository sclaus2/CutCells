// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "local_mesh.h"
#include "cell_subdivision.h"
#include "cell_topology.h"
#include "edge_classification.h"
#include "mapping.h"
#include "cell_flags.h"
#include "cut_cell.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <functional>
#include <iostream>
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

inline int background_edge_id_from_corners(cell::type ct, int pv0, int pv1)
{
    const auto parent_edges = cell::edges(ct);
    for (int pe = 0; pe < static_cast<int>(parent_edges.size()); ++pe)
    {
        const int e0 = parent_edges[static_cast<std::size_t>(pe)][0];
        const int e1 = parent_edges[static_cast<std::size_t>(pe)][1];
        if ((pv0 == e0 && pv1 == e1) || (pv0 == e1 && pv1 == e0))
            return pe;
    }
    return -1;
}

inline bool background_edge_on_hexahedron_face(int edge_id, int face_id)
{
    if (face_id < 0 || face_id >= static_cast<int>(cell::hexahedron_faces.size()))
        return false;
    const auto edge = cell::hexahedron_edges[static_cast<std::size_t>(edge_id)];
    const auto face = cell::hexahedron_faces[static_cast<std::size_t>(face_id)];
    const auto on_face = [&](int v)
    {
        for (const int fv : face)
        {
            if (fv == v)
                return true;
        }
        return false;
    };
    return on_face(edge[0]) && on_face(edge[1]);
}

template <std::floating_point T>
bool vertex_lies_on_background_edge(const LocalMesh<T>& mesh, int lv, int bg_edge_id)
{
    if (bg_edge_id < 0 || bg_edge_id >= cell::num_edges(mesh.parent_cell_type))
        return false;

    const int pdim = mesh.vertex_parent_dim[static_cast<std::size_t>(lv)];
    const int pid = mesh.vertex_parent_id[static_cast<std::size_t>(lv)];
    if (pdim == 1 && pid == bg_edge_id)
        return true;
    if (pdim == 0)
    {
        const auto edge = cell::edges(mesh.parent_cell_type)[static_cast<std::size_t>(bg_edge_id)];
        return pid == edge[0] || pid == edge[1];
    }
    return false;
}

template <std::floating_point T>
std::pair<int32_t, int32_t> infer_background_edge_parent(
    const LocalMesh<T>& mesh,
    int                 lv0,
    int                 lv1)
{
    const int n_parent_edges = cell::num_edges(mesh.parent_cell_type);
    for (int pe = 0; pe < n_parent_edges; ++pe)
    {
        if (vertex_lies_on_background_edge(mesh, lv0, pe)
            && vertex_lies_on_background_edge(mesh, lv1, pe))
        {
            return {1, pe};
        }
    }
    return {-1, -1};
}

template <std::floating_point T>
bool vertex_lies_on_background_face(const LocalMesh<T>& mesh, int lv, int bg_face_id)
{
    if (mesh.parent_cell_type != cell::type::hexahedron)
        return false;
    if (bg_face_id < 0 || bg_face_id >= static_cast<int>(cell::hexahedron_faces.size()))
        return false;

    const int pdim = mesh.vertex_parent_dim[static_cast<std::size_t>(lv)];
    const int pid = mesh.vertex_parent_id[static_cast<std::size_t>(lv)];
    if (pdim == 2 && pid == bg_face_id)
        return true;
    if (pdim == 1 && pid >= 0 && pid < cell::num_edges(mesh.parent_cell_type))
        return background_edge_on_hexahedron_face(pid, bg_face_id);
    if (pdim == 0)
    {
        const auto face = cell::hexahedron_faces[static_cast<std::size_t>(bg_face_id)];
        for (const int fv : face)
        {
            if (pid == fv)
                return true;
        }
    }
    return false;
}

template <std::floating_point T>
std::pair<int32_t, int32_t> infer_background_face_parent(
    const LocalMesh<T>&             mesh,
    std::span<const int32_t>        face_vertices)
{
    if (mesh.parent_cell_type != cell::type::hexahedron)
        return {-1, -1};

    for (int f = 0; f < static_cast<int>(cell::hexahedron_faces.size()); ++f)
    {
        bool all_on_face = true;
        for (const int32_t lv : face_vertices)
        {
            if (!vertex_lies_on_background_face(mesh, lv, f))
            {
                all_on_face = false;
                break;
            }
        }
        if (all_on_face)
            return {2, f};
    }
    return {-1, -1};
}

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
    mesh.edge_parent_dim.assign(n_edges, -1);
    mesh.edge_parent_id.assign(n_edges, -1);
    for (int i = 0; i < n_edges; ++i)
    {
        const auto [pdim, pid] = infer_background_edge_parent(
            mesh,
            mesh.edge_vertices[static_cast<std::size_t>(2 * i)],
            mesh.edge_vertices[static_cast<std::size_t>(2 * i + 1)]);
        mesh.edge_parent_dim[static_cast<std::size_t>(i)] = pdim;
        mesh.edge_parent_id[static_cast<std::size_t>(i)] = pid;
    }
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

    // Rewrite only parent-cell interior ids so they remain stable across later
    // local refinement. Edge ancestry stays in background-edge numbering.
    for (int lv = 0; lv < mesh.n_vertices(); ++lv)
    {
        const int pdim = mesh.vertex_parent_dim[static_cast<std::size_t>(lv)];
        int& pid = mesh.vertex_parent_id[static_cast<std::size_t>(lv)];
        if (pdim == mesh.tdim)
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

template <std::floating_point T>
int32_t append_generated_vertex(
    std::vector<T>&        vertex_x,
    std::vector<T>&        vertex_ref_x,
    std::vector<int32_t>&  vertex_parent_dim,
    std::vector<int32_t>&  vertex_parent_id,
    std::vector<int32_t>&  vertex_root_edge_id,
    std::span<const T>     x_ref,
    std::span<const T>     x_phys,
    int32_t                parent_dim,
    int32_t                parent_id)
{
    const int32_t new_id = static_cast<int32_t>(vertex_parent_dim.size());
    vertex_x.insert(vertex_x.end(), x_phys.begin(), x_phys.end());
    vertex_ref_x.insert(vertex_ref_x.end(), x_ref.begin(), x_ref.end());
    vertex_parent_dim.push_back(parent_dim);
    vertex_parent_id.push_back(parent_id);
    vertex_root_edge_id.push_back(-1);
    return new_id;
}

template <std::floating_point T>
int32_t append_segment_midpoint(
    const LocalMesh<T>&    old,
    int32_t                va,
    int32_t                vb,
    std::vector<T>&        vertex_x,
    std::vector<T>&        vertex_ref_x,
    std::vector<int32_t>&  vertex_parent_dim,
    std::vector<int32_t>&  vertex_parent_id,
    std::vector<int32_t>&  vertex_root_edge_id)
{
    std::array<T, 3> x_ref = {T(0), T(0), T(0)};
    std::array<T, 3> x_phys = {T(0), T(0), T(0)};
    for (int d = 0; d < old.tdim; ++d)
    {
        x_ref[static_cast<std::size_t>(d)]
            = T(0.5) * (old.vertex_ref_x[static_cast<std::size_t>(va * old.tdim + d)]
                        + old.vertex_ref_x[static_cast<std::size_t>(vb * old.tdim + d)]);
    }
    for (int d = 0; d < old.gdim; ++d)
    {
        x_phys[static_cast<std::size_t>(d)]
            = T(0.5) * (old.vertex_x[static_cast<std::size_t>(va * old.gdim + d)]
                        + old.vertex_x[static_cast<std::size_t>(vb * old.gdim + d)]);
    }

    const auto [pdim, pid] = infer_background_edge_parent(old, va, vb);
    return append_generated_vertex<T>(
        vertex_x,
        vertex_ref_x,
        vertex_parent_dim,
        vertex_parent_id,
        vertex_root_edge_id,
        std::span<const T>(x_ref.data(), static_cast<std::size_t>(old.tdim)),
        std::span<const T>(x_phys.data(), static_cast<std::size_t>(old.gdim)),
        pdim,
        pid);
}

template <std::floating_point T>
int32_t append_average_vertex(
    const LocalMesh<T>&    old,
    std::span<const int32_t> verts,
    int32_t                parent_dim,
    int32_t                parent_id,
    std::vector<T>&        vertex_x,
    std::vector<T>&        vertex_ref_x,
    std::vector<int32_t>&  vertex_parent_dim,
    std::vector<int32_t>&  vertex_parent_id,
    std::vector<int32_t>&  vertex_root_edge_id)
{
    std::array<T, 3> x_ref = {T(0), T(0), T(0)};
    std::array<T, 3> x_phys = {T(0), T(0), T(0)};
    const T inv = T(1) / static_cast<T>(verts.size());
    for (const int32_t v : verts)
    {
        for (int d = 0; d < old.tdim; ++d)
            x_ref[static_cast<std::size_t>(d)]
                += inv * old.vertex_ref_x[static_cast<std::size_t>(v * old.tdim + d)];
        for (int d = 0; d < old.gdim; ++d)
            x_phys[static_cast<std::size_t>(d)]
                += inv * old.vertex_x[static_cast<std::size_t>(v * old.gdim + d)];
    }

    return append_generated_vertex<T>(
        vertex_x,
        vertex_ref_x,
        vertex_parent_dim,
        vertex_parent_id,
        vertex_root_edge_id,
        std::span<const T>(x_ref.data(), static_cast<std::size_t>(old.tdim)),
        std::span<const T>(x_phys.data(), static_cast<std::size_t>(old.gdim)),
        parent_dim,
        parent_id);
}

template <typename Table>
void append_subdivided_cells(
    const Table&               table,
    std::span<const int32_t>   local_to_global,
    cell::type                 child_cell_type,
    std::vector<int32_t>&      cell_vertices,
    std::vector<int32_t>&      cell_offsets,
    std::vector<cell::type>&   cell_types)
{
    for (const auto& child : table)
    {
        for (const int lid : child)
            cell_vertices.push_back(local_to_global[static_cast<std::size_t>(lid)]);
        cell_offsets.push_back(static_cast<int32_t>(cell_vertices.size()));
        cell_types.push_back(child_cell_type);
    }
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
void red_refine_marked_cells(
    LocalMesh<T>&            mesh,
    std::span<const uint8_t> marked_cells)
{
    if (static_cast<int>(marked_cells.size()) != mesh.n_cells())
        throw std::invalid_argument("red_refine_marked_cells: marked_cells size mismatch");

    const LocalMesh<T> old = mesh;

    std::vector<int32_t> old_to_kept(static_cast<std::size_t>(old.n_vertices()), -1);
    std::vector<T> new_vertex_x;
    std::vector<T> new_vertex_ref_x;
    std::vector<int32_t> new_vertex_parent_dim;
    std::vector<int32_t> new_vertex_parent_id;
    std::vector<int32_t> new_vertex_root_edge_id;
    new_vertex_x.reserve(old.vertex_x.size());
    new_vertex_ref_x.reserve(old.vertex_ref_x.size());
    new_vertex_parent_dim.reserve(old.vertex_parent_dim.size());
    new_vertex_parent_id.reserve(old.vertex_parent_id.size());
    new_vertex_root_edge_id.reserve(old.vertex_root_edge_id.size());

    for (int32_t old_v = 0; old_v < old.n_vertices(); ++old_v)
    {
        if (old.vertex_root_edge_id[static_cast<std::size_t>(old_v)] >= 0)
            continue;
        const int32_t new_v = static_cast<int32_t>(new_vertex_parent_dim.size());
        old_to_kept[static_cast<std::size_t>(old_v)] = new_v;
        for (int d = 0; d < old.gdim; ++d)
            new_vertex_x.push_back(old.vertex_x[static_cast<std::size_t>(old_v * old.gdim + d)]);
        for (int d = 0; d < old.tdim; ++d)
            new_vertex_ref_x.push_back(old.vertex_ref_x[static_cast<std::size_t>(old_v * old.tdim + d)]);
        new_vertex_parent_dim.push_back(old.vertex_parent_dim[static_cast<std::size_t>(old_v)]);
        new_vertex_parent_id.push_back(old.vertex_parent_id[static_cast<std::size_t>(old_v)]);
        // Drop all cached roots from rejected candidate cuts.
        new_vertex_root_edge_id.push_back(-1);
    }

    std::vector<int32_t> new_cell_vertices;
    std::vector<int32_t> new_cell_offsets = {0};
    std::vector<cell::type> new_cell_types;
    new_cell_vertices.reserve(old.cell_vertices.size() * 4);
    new_cell_types.reserve(old.n_cells() * 4);

    for (int c = 0; c < old.n_cells(); ++c)
    {
        const int c0 = old.cell_offsets[static_cast<std::size_t>(c)];
        const int c1 = old.cell_offsets[static_cast<std::size_t>(c + 1)];
        const std::span<const int32_t> old_cell_verts(
            old.cell_vertices.data() + static_cast<std::size_t>(c0),
            static_cast<std::size_t>(c1 - c0));
        const cell::type ct = old.cell_types[static_cast<std::size_t>(c)];
        std::vector<int32_t> cell_verts(static_cast<std::size_t>(old_cell_verts.size()), -1);
        for (std::size_t i = 0; i < old_cell_verts.size(); ++i)
        {
            const int32_t old_v = old_cell_verts[i];
            if (old_v < 0 || old_v >= old.n_vertices())
                throw std::invalid_argument("red_refine_marked_cells: invalid cell vertex id");
            cell_verts[i] = old_to_kept[static_cast<std::size_t>(old_v)];
            if (cell_verts[i] < 0)
                throw std::invalid_argument("red_refine_marked_cells: cell references dropped root vertex");
        }

        if (!marked_cells[static_cast<std::size_t>(c)])
        {
            new_cell_vertices.insert(new_cell_vertices.end(), cell_verts.begin(), cell_verts.end());
            new_cell_offsets.push_back(static_cast<int32_t>(new_cell_vertices.size()));
            new_cell_types.push_back(ct);
            continue;
        }

        switch (ct)
        {
        case cell::type::interval:
        {
            std::array<int32_t, 3> local = {
                cell_verts[0],
                cell_verts[1],
                append_segment_midpoint(
                    old, old_cell_verts[0], old_cell_verts[1],
                    new_vertex_x, new_vertex_ref_x,
                    new_vertex_parent_dim, new_vertex_parent_id, new_vertex_root_edge_id)};
            append_subdivided_cells(
                cell::interval_subdivision_table,
                std::span<const int32_t>(local.data(), local.size()),
                cell::type::interval,
                new_cell_vertices,
                new_cell_offsets,
                new_cell_types);
            break;
        }
        case cell::type::triangle:
        {
            std::array<int32_t, 6> local = {
                cell_verts[0], cell_verts[1], cell_verts[2], -1, -1, -1};
            const auto edges = cell::edges(cell::type::triangle);
            for (int e = 0; e < 3; ++e)
            {
                const int a = edges[static_cast<std::size_t>(e)][0];
                const int b = edges[static_cast<std::size_t>(e)][1];
                local[static_cast<std::size_t>(3 + e)] = append_segment_midpoint(
                    old, old_cell_verts[static_cast<std::size_t>(a)], old_cell_verts[static_cast<std::size_t>(b)],
                    new_vertex_x, new_vertex_ref_x,
                    new_vertex_parent_dim, new_vertex_parent_id, new_vertex_root_edge_id);
            }
            append_subdivided_cells(
                cell::triangle_subdivision_table,
                std::span<const int32_t>(local.data(), local.size()),
                cell::type::triangle,
                new_cell_vertices,
                new_cell_offsets,
                new_cell_types);
            break;
        }
        case cell::type::quadrilateral:
        {
            std::array<int32_t, 9> local = {
                cell_verts[0], cell_verts[1], cell_verts[2], cell_verts[3],
                -1, -1, -1, -1, -1};
            const auto edges = cell::edges(cell::type::quadrilateral);
            for (int e = 0; e < 4; ++e)
            {
                const int a = edges[static_cast<std::size_t>(e)][0];
                const int b = edges[static_cast<std::size_t>(e)][1];
                local[static_cast<std::size_t>(4 + e)] = append_segment_midpoint(
                    old, old_cell_verts[static_cast<std::size_t>(a)], old_cell_verts[static_cast<std::size_t>(b)],
                    new_vertex_x, new_vertex_ref_x,
                    new_vertex_parent_dim, new_vertex_parent_id, new_vertex_root_edge_id);
            }
            local[8] = append_average_vertex(
                old,
                old_cell_verts,
                old.tdim,
                old.parent_cell_id,
                new_vertex_x, new_vertex_ref_x,
                new_vertex_parent_dim, new_vertex_parent_id, new_vertex_root_edge_id);
            append_subdivided_cells(
                cell::quadrilateral_subdivision_table,
                std::span<const int32_t>(local.data(), local.size()),
                cell::type::quadrilateral,
                new_cell_vertices,
                new_cell_offsets,
                new_cell_types);
            break;
        }
        case cell::type::tetrahedron:
        {
            std::array<int32_t, 10> local = {
                cell_verts[0], cell_verts[1], cell_verts[2], cell_verts[3],
                -1, -1, -1, -1, -1, -1};
            const auto edges = cell::edges(cell::type::tetrahedron);
            for (int e = 0; e < 6; ++e)
            {
                const int a = edges[static_cast<std::size_t>(e)][0];
                const int b = edges[static_cast<std::size_t>(e)][1];
                local[static_cast<std::size_t>(4 + e)] = append_segment_midpoint(
                    old, old_cell_verts[static_cast<std::size_t>(a)], old_cell_verts[static_cast<std::size_t>(b)],
                    new_vertex_x, new_vertex_ref_x,
                    new_vertex_parent_dim, new_vertex_parent_id, new_vertex_root_edge_id);
            }
            append_subdivided_cells(
                cell::tetrahedron_subdivision_table,
                std::span<const int32_t>(local.data(), local.size()),
                cell::type::tetrahedron,
                new_cell_vertices,
                new_cell_offsets,
                new_cell_types);
            break;
        }
        case cell::type::hexahedron:
        {
            std::array<int32_t, 27> local = {
                cell_verts[0], cell_verts[1], cell_verts[2], cell_verts[3],
                cell_verts[4], cell_verts[5], cell_verts[6], cell_verts[7],
                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1, -1};
            const auto edges = cell::edges(cell::type::hexahedron);
            for (int e = 0; e < 12; ++e)
            {
                const int a = edges[static_cast<std::size_t>(e)][0];
                const int b = edges[static_cast<std::size_t>(e)][1];
                local[static_cast<std::size_t>(8 + e)] = append_segment_midpoint(
                    old, old_cell_verts[static_cast<std::size_t>(a)], old_cell_verts[static_cast<std::size_t>(b)],
                    new_vertex_x, new_vertex_ref_x,
                    new_vertex_parent_dim, new_vertex_parent_id, new_vertex_root_edge_id);
            }
            for (int f = 0; f < 6; ++f)
            {
                const auto face = cell::hexahedron_faces[static_cast<std::size_t>(f)];
                std::array<int32_t, 4> face_verts = {
                    old_cell_verts[static_cast<std::size_t>(face[0])],
                    old_cell_verts[static_cast<std::size_t>(face[1])],
                    old_cell_verts[static_cast<std::size_t>(face[2])],
                    old_cell_verts[static_cast<std::size_t>(face[3])]};
                const auto [pdim, pid] = infer_background_face_parent(
                    old, std::span<const int32_t>(face_verts.data(), face_verts.size()));
                local[static_cast<std::size_t>(20 + f)] = append_average_vertex(
                    old,
                    std::span<const int32_t>(face_verts.data(), face_verts.size()),
                    pdim,
                    pid,
                    new_vertex_x, new_vertex_ref_x,
                    new_vertex_parent_dim, new_vertex_parent_id, new_vertex_root_edge_id);
            }
            local[26] = append_average_vertex(
                old,
                old_cell_verts,
                old.tdim,
                old.parent_cell_id,
                new_vertex_x, new_vertex_ref_x,
                new_vertex_parent_dim, new_vertex_parent_id, new_vertex_root_edge_id);
            append_subdivided_cells(
                cell::hexahedron_subdivision_table,
                std::span<const int32_t>(local.data(), local.size()),
                cell::type::hexahedron,
                new_cell_vertices,
                new_cell_offsets,
                new_cell_types);
            break;
        }
        default:
            throw std::invalid_argument("red_refine_marked_cells: unsupported cell type");
        }
    }

    mesh.vertex_x = std::move(new_vertex_x);
    mesh.vertex_ref_x = std::move(new_vertex_ref_x);
    mesh.vertex_parent_dim = std::move(new_vertex_parent_dim);
    mesh.vertex_parent_id = std::move(new_vertex_parent_id);
    mesh.vertex_root_edge_id = std::move(new_vertex_root_edge_id);
    mesh.cell_vertices = std::move(new_cell_vertices);
    mesh.cell_offsets = std::move(new_cell_offsets);
    mesh.cell_types = std::move(new_cell_types);
    mesh.cell_domain.assign(static_cast<std::size_t>(mesh.cell_types.size()),
                            static_cast<uint8_t>(cell::domain::unset));

    deduplicate_vertices_exact(mesh);

    const int nv = mesh.n_vertices();
    mesh.vertex_zero_mask.assign(static_cast<std::size_t>(nv), 0);
    mesh.vertex_inside_mask.assign(static_cast<std::size_t>(nv), 0);
    mesh.vertex_phi.assign(static_cast<std::size_t>(nv * mesh.n_level_sets), T(0));

    build_local_edges(mesh);
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
void clear_edge_root_cache(LocalMesh<T>& mesh)
{
    std::fill(mesh.edge_root_vertex.begin(), mesh.edge_root_vertex.end(), -1);
    std::fill(mesh.edge_root_parameter.begin(), mesh.edge_root_parameter.end(), T(-1));
    std::fill(mesh.edge_root_iterations.begin(), mesh.edge_root_iterations.end(), 0);
    std::fill(mesh.edge_root_evaluations.begin(), mesh.edge_root_evaluations.end(), 0);
    std::fill(mesh.edge_root_converged.begin(), mesh.edge_root_converged.end(), uint8_t(0));
    std::fill(mesh.edge_root_residual.begin(), mesh.edge_root_residual.end(), T(0));
}

template <std::floating_point T>
void drop_unreferenced_root_vertices(LocalMesh<T>& mesh)
{
    const int old_nv = mesh.n_vertices();
    if (old_nv == 0 || static_cast<int>(mesh.vertex_root_edge_id.size()) != old_nv)
        return;

    std::vector<uint8_t> referenced(static_cast<std::size_t>(old_nv), uint8_t(0));
    for (const int32_t v : mesh.cell_vertices)
    {
        if (v >= 0 && v < old_nv)
            referenced[static_cast<std::size_t>(v)] = 1;
    }

    std::vector<int32_t> old_to_new(static_cast<std::size_t>(old_nv), -1);
    std::vector<T> new_x;
    std::vector<T> new_ref_x;
    std::vector<int32_t> new_parent_dim;
    std::vector<int32_t> new_parent_id;
    std::vector<int32_t> new_root_edge_id;
    std::vector<T> new_phi;
    std::vector<uint64_t> new_zero_mask;
    std::vector<uint64_t> new_inside_mask;

    new_x.reserve(mesh.vertex_x.size());
    new_ref_x.reserve(mesh.vertex_ref_x.size());
    new_parent_dim.reserve(mesh.vertex_parent_dim.size());
    new_parent_id.reserve(mesh.vertex_parent_id.size());
    new_root_edge_id.reserve(mesh.vertex_root_edge_id.size());
    new_phi.reserve(mesh.vertex_phi.size());
    new_zero_mask.reserve(mesh.vertex_zero_mask.size());
    new_inside_mask.reserve(mesh.vertex_inside_mask.size());

    int new_nv = 0;
    for (int old_v = 0; old_v < old_nv; ++old_v)
    {
        const bool is_root = mesh.vertex_root_edge_id[static_cast<std::size_t>(old_v)] >= 0;
        if (is_root && !referenced[static_cast<std::size_t>(old_v)])
            continue;

        old_to_new[static_cast<std::size_t>(old_v)] = new_nv++;
        for (int d = 0; d < mesh.gdim; ++d)
        {
            new_x.push_back(
                mesh.vertex_x[static_cast<std::size_t>(old_v * mesh.gdim + d)]);
        }
        if (!mesh.vertex_ref_x.empty())
        {
            for (int d = 0; d < mesh.tdim; ++d)
            {
                new_ref_x.push_back(
                    mesh.vertex_ref_x[static_cast<std::size_t>(old_v * mesh.tdim + d)]);
            }
        }
        new_parent_dim.push_back(mesh.vertex_parent_dim[static_cast<std::size_t>(old_v)]);
        new_parent_id.push_back(mesh.vertex_parent_id[static_cast<std::size_t>(old_v)]);
        new_root_edge_id.push_back(mesh.vertex_root_edge_id[static_cast<std::size_t>(old_v)]);
        if (!mesh.vertex_phi.empty())
        {
            for (int ls = 0; ls < mesh.n_level_sets; ++ls)
            {
                new_phi.push_back(
                    mesh.vertex_phi[static_cast<std::size_t>(old_v * mesh.n_level_sets + ls)]);
            }
        }
        if (!mesh.vertex_zero_mask.empty())
            new_zero_mask.push_back(mesh.vertex_zero_mask[static_cast<std::size_t>(old_v)]);
        if (!mesh.vertex_inside_mask.empty())
            new_inside_mask.push_back(mesh.vertex_inside_mask[static_cast<std::size_t>(old_v)]);
    }

    if (new_nv == old_nv)
        return;

    for (int32_t& v : mesh.cell_vertices)
        v = old_to_new[static_cast<std::size_t>(v)];

    mesh.vertex_x.swap(new_x);
    mesh.vertex_ref_x.swap(new_ref_x);
    mesh.vertex_parent_dim.swap(new_parent_dim);
    mesh.vertex_parent_id.swap(new_parent_id);
    mesh.vertex_root_edge_id.swap(new_root_edge_id);
    mesh.vertex_phi.swap(new_phi);
    mesh.vertex_zero_mask.swap(new_zero_mask);
    mesh.vertex_inside_mask.swap(new_inside_mask);
}

template <std::floating_point T>
void deduplicate_zero_vertices(LocalMesh<T>& mesh,
                               int           level_set_id,
                               T             coord_tol,
                               T)
{
    const int old_nv = mesh.n_vertices();
    if (old_nv <= 1
        || static_cast<int>(mesh.vertex_phi.size()) != old_nv * mesh.n_level_sets
        || level_set_id < 0
        || level_set_id >= mesh.n_level_sets)
    {
        return;
    }

    std::vector<int32_t> old_to_new(static_cast<std::size_t>(old_nv), -1);
    std::vector<T> new_x;
    std::vector<T> new_ref_x;
    std::vector<int32_t> new_parent_dim;
    std::vector<int32_t> new_parent_id;
    std::vector<int32_t> new_root_edge_id;
    std::vector<T> new_phi;
    std::vector<uint64_t> new_zero_mask;
    std::vector<uint64_t> new_inside_mask;

    new_x.reserve(mesh.vertex_x.size());
    new_ref_x.reserve(mesh.vertex_ref_x.size());
    new_parent_dim.reserve(mesh.vertex_parent_dim.size());
    new_parent_id.reserve(mesh.vertex_parent_id.size());
    new_root_edge_id.reserve(mesh.vertex_root_edge_id.size());
    new_phi.reserve(mesh.vertex_phi.size());
    new_zero_mask.reserve(mesh.vertex_zero_mask.size());
    new_inside_mask.reserve(mesh.vertex_inside_mask.size());

    const auto append_vertex = [&](int old_v)
    {
        for (int d = 0; d < mesh.gdim; ++d)
        {
            new_x.push_back(mesh.vertex_x[static_cast<std::size_t>(old_v * mesh.gdim + d)]);
        }
        if (!mesh.vertex_ref_x.empty())
        {
            for (int d = 0; d < mesh.tdim; ++d)
            {
                new_ref_x.push_back(
                    mesh.vertex_ref_x[static_cast<std::size_t>(old_v * mesh.tdim + d)]);
            }
        }
        new_parent_dim.push_back(mesh.vertex_parent_dim[static_cast<std::size_t>(old_v)]);
        new_parent_id.push_back(mesh.vertex_parent_id[static_cast<std::size_t>(old_v)]);
        new_root_edge_id.push_back(mesh.vertex_root_edge_id[static_cast<std::size_t>(old_v)]);
        for (int ls = 0; ls < mesh.n_level_sets; ++ls)
        {
            new_phi.push_back(mesh.vertex_phi[static_cast<std::size_t>(old_v * mesh.n_level_sets + ls)]);
        }
        new_zero_mask.push_back(mesh.vertex_zero_mask[static_cast<std::size_t>(old_v)]);
        new_inside_mask.push_back(mesh.vertex_inside_mask[static_cast<std::size_t>(old_v)]);
    };

    auto overwrite_vertex = [&](int new_v, int old_v)
    {
        for (int d = 0; d < mesh.gdim; ++d)
        {
            new_x[static_cast<std::size_t>(new_v * mesh.gdim + d)]
                = mesh.vertex_x[static_cast<std::size_t>(old_v * mesh.gdim + d)];
        }
        if (!mesh.vertex_ref_x.empty())
        {
            for (int d = 0; d < mesh.tdim; ++d)
            {
                new_ref_x[static_cast<std::size_t>(new_v * mesh.tdim + d)]
                    = mesh.vertex_ref_x[static_cast<std::size_t>(old_v * mesh.tdim + d)];
            }
        }
        new_parent_dim[static_cast<std::size_t>(new_v)]
            = mesh.vertex_parent_dim[static_cast<std::size_t>(old_v)];
        new_parent_id[static_cast<std::size_t>(new_v)]
            = mesh.vertex_parent_id[static_cast<std::size_t>(old_v)];
        new_root_edge_id[static_cast<std::size_t>(new_v)]
            = mesh.vertex_root_edge_id[static_cast<std::size_t>(old_v)];
        for (int ls = 0; ls < mesh.n_level_sets; ++ls)
        {
            new_phi[static_cast<std::size_t>(new_v * mesh.n_level_sets + ls)]
                = mesh.vertex_phi[static_cast<std::size_t>(old_v * mesh.n_level_sets + ls)];
        }
        new_zero_mask[static_cast<std::size_t>(new_v)]
            = mesh.vertex_zero_mask[static_cast<std::size_t>(old_v)];
        new_inside_mask[static_cast<std::size_t>(new_v)]
            = mesh.vertex_inside_mask[static_cast<std::size_t>(old_v)];
    };

    int new_nv = 0;
    for (int old_v = 0; old_v < old_nv; ++old_v)
    {
        int found = -1;
        for (int new_v = 0; new_v < new_nv; ++new_v)
        {
            bool same = true;
            for (int d = 0; d < mesh.gdim; ++d)
            {
                const T dv
                    = new_x[static_cast<std::size_t>(new_v * mesh.gdim + d)]
                      - mesh.vertex_x[static_cast<std::size_t>(old_v * mesh.gdim + d)];
                if (std::abs(dv) > coord_tol)
                {
                    same = false;
                    break;
                }
            }
            if (same)
            {
                found = new_v;
                break;
            }
        }

        if (found < 0)
        {
            old_to_new[static_cast<std::size_t>(old_v)] = new_nv++;
            append_vertex(old_v);
            continue;
        }

        old_to_new[static_cast<std::size_t>(old_v)] = found;
        const T old_abs_phi = std::abs(
            mesh.vertex_phi[static_cast<std::size_t>(old_v * mesh.n_level_sets + level_set_id)]);
        const T new_abs_phi = std::abs(
            new_phi[static_cast<std::size_t>(found * mesh.n_level_sets + level_set_id)]);
        if (old_abs_phi < new_abs_phi)
            overwrite_vertex(found, old_v);

        new_zero_mask[static_cast<std::size_t>(found)]
            |= mesh.vertex_zero_mask[static_cast<std::size_t>(old_v)];
        new_inside_mask[static_cast<std::size_t>(found)]
            |= mesh.vertex_inside_mask[static_cast<std::size_t>(old_v)];
    }

    if (new_nv == old_nv)
        return;

    for (int32_t& v : mesh.cell_vertices)
        v = old_to_new[static_cast<std::size_t>(v)];

    mesh.vertex_x.swap(new_x);
    mesh.vertex_ref_x.swap(new_ref_x);
    mesh.vertex_parent_dim.swap(new_parent_dim);
    mesh.vertex_parent_id.swap(new_parent_id);
    mesh.vertex_root_edge_id.swap(new_root_edge_id);
    mesh.vertex_phi.swap(new_phi);
    mesh.vertex_zero_mask.swap(new_zero_mask);
    mesh.vertex_inside_mask.swap(new_inside_mask);
}

template <std::floating_point T>
constexpr T root_vertex_merge_tol()
{
    if constexpr (std::is_same_v<T, float>)
        return static_cast<T>(1e-5);
    else
        return static_cast<T>(1e-11);
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
    mesh.vertex_parent_dim.push_back(mesh.edge_parent_dim[static_cast<std::size_t>(edge_id)]);
    mesh.vertex_parent_id.push_back(mesh.edge_parent_id[static_cast<std::size_t>(edge_id)]);
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
int compute_edge_root_bernstein(
    LocalMesh<T>&                      mesh,
    const LocalLevelSetFunction<T, I>& level_set,
    int                                edge_id,
    int                                level_set_id,
    T                                  tol)
{
    if (edge_id < 0 || edge_id >= mesh.n_edges())
        throw std::invalid_argument("compute_edge_root_bernstein: invalid edge id");
    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("compute_edge_root_bernstein: invalid level_set_id");
    if (mesh.edge_root_vertex[static_cast<std::size_t>(edge_id)] >= 0)
        return mesh.edge_root_vertex[static_cast<std::size_t>(edge_id)];
    if (mesh.vertex_ref_x.size() != static_cast<std::size_t>(mesh.n_vertices() * mesh.tdim))
        throw std::invalid_argument("compute_edge_root_bernstein: missing reference coordinates");

    const int v0_idx = mesh.edge_vertices[static_cast<std::size_t>(2 * edge_id)];
    const int v1_idx = mesh.edge_vertices[static_cast<std::size_t>(2 * edge_id + 1)];

    const std::span<const T> x0_ref(
        mesh.vertex_ref_x.data() + static_cast<std::size_t>(v0_idx * mesh.tdim),
        static_cast<std::size_t>(mesh.tdim));
    const std::span<const T> x1_ref(
        mesh.vertex_ref_x.data() + static_cast<std::size_t>(v1_idx * mesh.tdim),
        static_cast<std::size_t>(mesh.tdim));

    std::vector<T> edge_coeffs;
    if (level_set.has_segment_restriction())
    {
        level_set.segment_restriction(x0_ref, x1_ref, edge_coeffs);
    }
    else if (level_set.has_edge_restriction()
             && mesh.edge_parent_dim[static_cast<std::size_t>(edge_id)] == 1
             && mesh.edge_parent_id[static_cast<std::size_t>(edge_id)] >= 0)
    {
        edge_coeffs = level_set.edge_restriction(
            mesh.edge_parent_id[static_cast<std::size_t>(edge_id)]).coeffs;
    }
    else
    {
        throw std::invalid_argument(
            "compute_edge_root_bernstein: no Bernstein restriction available for local edge");
    }
    const auto root_info = bernstein_edge_root_info<T>(
        std::span<const T>(edge_coeffs.data(), edge_coeffs.size()), 64, tol, tol);
    const T t = std::clamp(root_info.t, T(0), T(1));

    const int new_v_idx = mesh.n_vertices();
    const int gdim = mesh.gdim;
    const int tdim = mesh.tdim;
    const int n_ls = mesh.n_level_sets;

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

    mesh.vertex_parent_dim.push_back(mesh.edge_parent_dim[static_cast<std::size_t>(edge_id)]);
    mesh.vertex_parent_id.push_back(mesh.edge_parent_id[static_cast<std::size_t>(edge_id)]);
    mesh.vertex_root_edge_id.push_back(edge_id);

    mesh.vertex_phi.resize(mesh.vertex_phi.size() + static_cast<std::size_t>(n_ls), T(0));
    mesh.vertex_zero_mask.push_back(0);
    mesh.vertex_inside_mask.push_back(0);

    const T* x_root_ref = mesh.vertex_ref_x.data() + static_cast<std::size_t>(new_v_idx * tdim);
    const T val = level_set.value(x_root_ref);
    mesh.vertex_phi[static_cast<std::size_t>(new_v_idx * n_ls + level_set_id)] = val;
    const uint64_t mask = (uint64_t(1) << level_set_id);
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

    mesh.vertex_parent_dim.push_back(mesh.edge_parent_dim[static_cast<std::size_t>(edge_id)]);
    mesh.vertex_parent_id.push_back(mesh.edge_parent_id[static_cast<std::size_t>(edge_id)]);
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
    drop_unreferenced_root_vertices(mesh);
    clear_edge_root_cache(mesh);
    int ne = mesh.n_edges();
    for (int i = 0; i < ne; ++i) {
        if (mesh.edge_state[i] == static_cast<uint8_t>(EdgeState::one_root)) {
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

    drop_unreferenced_root_vertices(mesh);
    clear_edge_root_cache(mesh);

    const int ne = mesh.n_edges();
    for (int i = 0; i < ne; ++i)
    {
        if (mesh.edge_state[static_cast<std::size_t>(i)] != static_cast<uint8_t>(EdgeState::one_root))
            continue;
        compute_edge_root<T, I>(mesh, level_set, i, level_set_id, root_method);
    }
}

template <std::floating_point T, std::integral I>
void compute_all_roots_with_backend(LocalMesh<T>& mesh,
                                    const LevelSetFunction<T, I>& level_set,
                                    LocalLevelSetBackend backend,
                                    int level_set_id,
                                    cell::edge_root::method root_method,
                                    T tol)
{
    if (backend != LocalLevelSetBackend::bernstein)
    {
        compute_all_roots<T, I>(mesh, level_set, level_set_id, root_method);
        return;
    }

    if (!level_set.has_nodal_values())
        throw std::invalid_argument("compute_all_roots_with_backend: Bernstein backend requires nodal values");
    const auto local_level_set
        = level_set.mesh != nullptr
              ? make_local_level_set_function_bernstein<T, I>(
                    *level_set.mesh, level_set, static_cast<I>(mesh.parent_cell_id))
              : make_local_level_set_function_bernstein<T, I>(
                    mesh.parent_cell_type, mesh.gdim, level_set,
                    static_cast<I>(mesh.parent_cell_id));

    drop_unreferenced_root_vertices(mesh);
    clear_edge_root_cache(mesh);

    const int ne = mesh.n_edges();
    for (int i = 0; i < ne; ++i)
    {
        if (mesh.edge_state[static_cast<std::size_t>(i)] != static_cast<uint8_t>(EdgeState::one_root))
            continue;
        compute_edge_root_bernstein<T, I>(mesh, local_level_set, i, level_set_id, tol);
    }
}

namespace
{
template <std::floating_point T>
int find_vertex_by_coords(
    const LocalMesh<T>& mesh,
    const std::span<const T> x,
    T tol = root_vertex_merge_tol<T>())
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

inline int sign_with_cut_tol(double v, double tol)
{
    if (v > tol)
        return 1;
    if (v < -tol)
        return -1;
    return 0;
}

template <std::floating_point T>
inline int sign_with_cut_tol(T v, T tol)
{
    if (v > tol)
        return 1;
    if (v < -tol)
        return -1;
    return 0;
}

inline bool is_cut_original_vertex_token(cell::type parent_cell_type, int32_t token)
{
    return token >= 100
           && token < 100 + cell::get_num_vertices(parent_cell_type);
}

inline bool is_cut_root_token(cell::type parent_cell_type, int32_t token)
{
    return token >= 0 && token < cell::num_edges(parent_cell_type);
}

template <std::floating_point T>
EdgeOrigin classify_cut_edge_origin(
    cell::type                parent_cell_type,
    const cell::CutCell<T>&   cut_cell,
    int                       lv0,
    int                       lv1)
{
    if (static_cast<int>(cut_cell._vertex_parent_entity.size())
        != static_cast<int>(cut_cell._vertex_coords.size()) / cut_cell._gdim)
    {
        return EdgeOrigin::lut_new;
    }

    const int32_t token0 = cut_cell._vertex_parent_entity[static_cast<std::size_t>(lv0)];
    const int32_t token1 = cut_cell._vertex_parent_entity[static_cast<std::size_t>(lv1)];

    if (is_cut_original_vertex_token(parent_cell_type, token0)
        && is_cut_original_vertex_token(parent_cell_type, token1))
    {
        const int pv0 = token0 - 100;
        const int pv1 = token1 - 100;
        return background_edge_id_from_corners(parent_cell_type, pv0, pv1) >= 0
                   ? EdgeOrigin::original
                   : EdgeOrigin::lut_new;
    }

    if (is_cut_root_token(parent_cell_type, token0)
        && is_cut_original_vertex_token(parent_cell_type, token1))
    {
        const auto edge = cut_edge_order(parent_cell_type)[static_cast<std::size_t>(token0)];
        const int pv = token1 - 100;
        return (pv == edge[0] || pv == edge[1])
                   ? EdgeOrigin::root_split
                   : EdgeOrigin::lut_new;
    }

    if (is_cut_root_token(parent_cell_type, token1)
        && is_cut_original_vertex_token(parent_cell_type, token0))
    {
        const auto edge = cut_edge_order(parent_cell_type)[static_cast<std::size_t>(token1)];
        const int pv = token0 - 100;
        return (pv == edge[0] || pv == edge[1])
                   ? EdgeOrigin::root_split
                   : EdgeOrigin::lut_new;
    }

    return EdgeOrigin::lut_new;
}

template <std::floating_point T>
void build_subcell_edges(
    const cell::CutCell<T>&   cut_cell,
    cell::type                parent_cell_type,
    int                       cell_id,
    std::vector<int32_t>&     edges,
    std::vector<uint8_t>&     edge_origin)
{
    const auto subcell_type = cut_cell._types[static_cast<std::size_t>(cell_id)];
    const auto subcell_vertices = cell::cell_vertices(cut_cell, cell_id);
    const auto edge_patterns = cell::edges(subcell_type);

    edges.clear();
    edge_origin.clear();
    edges.reserve(edge_patterns.size() * 2);
    edge_origin.reserve(edge_patterns.size());

    for (const auto& ep : edge_patterns)
    {
        const int lv0 = subcell_vertices[static_cast<std::size_t>(ep[0])];
        const int lv1 = subcell_vertices[static_cast<std::size_t>(ep[1])];
        edges.push_back(lv0);
        edges.push_back(lv1);
        edge_origin.push_back(static_cast<uint8_t>(
            classify_cut_edge_origin(parent_cell_type, cut_cell, lv0, lv1)));
    }
}

template <std::floating_point T, std::integral I>
bool cut_vertex_level_set_value(
    const cell::CutCell<T>&    cut_cell,
    std::span<const T>         parent_ls_values,
    const LevelSetFunction<T, I>& level_set,
    int                        parent_cell_id,
    int                        local_vertex_id,
    T&                         value)
{
    if (static_cast<int>(cut_cell._vertex_parent_entity.size())
        == static_cast<int>(cut_cell._vertex_coords.size()) / cut_cell._gdim)
    {
        const int32_t token = cut_cell._vertex_parent_entity[static_cast<std::size_t>(local_vertex_id)];
        if (is_cut_root_token(cut_cell._parent_cell_type, token))
        {
            value = T(0);
            return true;
        }
        if (token >= 100
            && token < 100 + static_cast<int32_t>(parent_ls_values.size()))
        {
            value = parent_ls_values[static_cast<std::size_t>(token - 100)];
            return true;
        }
    }

    if (level_set.has_value())
    {
        const T* x = cut_cell._vertex_coords.data()
                     + static_cast<std::size_t>(local_vertex_id * cut_cell._gdim);
        value = level_set.value(x, static_cast<I>(parent_cell_id));
        return true;
    }

    return false;
}

template <std::floating_point T, std::integral I>
bool edge_contains_root(
    const LevelSetFunction<T, I>&  level_set,
    const cell::CutCell<T>&        cut_cell,
    std::span<const T>             parent_ls_values,
    int                            parent_cell_id,
    int                            lv0,
    int                            lv1,
    T                              tol)
{
    T phi0 = T(0);
    T phi1 = T(0);
    if (!cut_vertex_level_set_value(
            cut_cell, parent_ls_values, level_set, parent_cell_id, lv0, phi0)
        || !cut_vertex_level_set_value(
            cut_cell, parent_ls_values, level_set, parent_cell_id, lv1, phi1))
    {
        return true;
    }

    const int s0 = sign_with_cut_tol(phi0, tol);
    const int s1 = sign_with_cut_tol(phi1, tol);
    if (s0 * s1 < 0)
        return true;

    if (!level_set.has_value())
        return false;

    std::vector<T> xmid(static_cast<std::size_t>(cut_cell._gdim), T(0));
    const T* x0 = cut_cell._vertex_coords.data()
                  + static_cast<std::size_t>(lv0 * cut_cell._gdim);
    const T* x1 = cut_cell._vertex_coords.data()
                  + static_cast<std::size_t>(lv1 * cut_cell._gdim);
    for (int d = 0; d < cut_cell._gdim; ++d)
        xmid[static_cast<std::size_t>(d)] = T(0.5) * (x0[d] + x1[d]);

    const T phim = level_set.value(xmid.data(), static_cast<I>(parent_cell_id));
    const int sm = sign_with_cut_tol(phim, tol);
    if (sm == 0)
        return true;
    if (s0 != 0 && s0 * sm < 0)
        return true;
    if (s1 != 0 && sm * s1 < 0)
        return true;
    return false;
}

template <std::floating_point T, std::integral I>
bool certify_subcell(
    const LevelSetFunction<T, I>&  level_set,
    const cell::CutCell<T>&        cut_cell,
    cell::type                     parent_cell_type,
    std::span<const T>             parent_ls_values,
    int                            parent_cell_id,
    int                            cell_id,
    T                              tol,
    bool                           debug)
{
    std::vector<int32_t> edges;
    std::vector<uint8_t> edge_origin;
    build_subcell_edges(cut_cell, parent_cell_type, cell_id, edges, edge_origin);

    for (std::size_t e = 0; e < edge_origin.size(); ++e)
    {
        const auto origin = static_cast<EdgeOrigin>(edge_origin[e]);
        if (origin == EdgeOrigin::original)
            continue;

        const int lv0 = edges[2 * e];
        const int lv1 = edges[2 * e + 1];
        if (edge_contains_root(
                level_set, cut_cell, parent_ls_values, parent_cell_id, lv0, lv1, tol))
        {
            if (debug)
                std::cout << "invalid edge detected\n";
            return false;
        }
    }

    return true;
}

template <std::floating_point T, std::integral I>
bool certify_cut(
    const LevelSetFunction<T, I>&  level_set,
    const cell::CutCell<T>&        cut_cell,
    cell::type                     parent_cell_type,
    std::span<const T>             parent_ls_values,
    int                            parent_cell_id,
    T                              tol,
    bool                           debug)
{
    const int nsub = cell::num_cells(cut_cell);
    for (int i = 0; i < nsub; ++i)
    {
        if (!certify_subcell(
                level_set, cut_cell, parent_cell_type, parent_ls_values,
                parent_cell_id, i, tol, debug))
        {
            return false;
        }
    }
    return true;
}

template <std::floating_point T, std::integral I>
bool certify_cut_on_local_mesh(
    LocalMesh<T>&                 mesh,
    const LevelSetFunction<T, I>& level_set,
    int                           level_set_id,
    bool                          triangulate,
    T                             tol,
    std::vector<uint8_t>&         marked_cells,
    int&                          invalid_cells,
    bool                          debug)
{
    marked_cells.assign(static_cast<std::size_t>(mesh.n_cells()), uint8_t(0));
    invalid_cells = 0;

    const int nc = mesh.n_cells();
    for (int c = 0; c < nc; ++c)
    {
        const int c0 = mesh.cell_offsets[static_cast<std::size_t>(c)];
        const int c1 = mesh.cell_offsets[static_cast<std::size_t>(c + 1)];
        const int nv = c1 - c0;
        const cell::type ct = mesh.cell_types[static_cast<std::size_t>(c)];

        const std::span<const int32_t> parent_cell_vertices(
            mesh.cell_vertices.data() + static_cast<std::size_t>(c0),
            static_cast<std::size_t>(nv));

        const int ce0 = mesh.cell_edge_offsets[static_cast<std::size_t>(c)];
        const int ce1 = mesh.cell_edge_offsets[static_cast<std::size_t>(c + 1)];
        const std::span<const int32_t> parent_cell_edges(
            mesh.cell_edges_flat.data() + static_cast<std::size_t>(ce0),
            static_cast<std::size_t>(ce1 - ce0));

        std::vector<T> ls_values(static_cast<std::size_t>(nv), T(0));
        std::vector<T> cell_vertex_coords(static_cast<std::size_t>(nv * mesh.gdim), T(0));
        for (int j = 0; j < nv; ++j)
        {
            const int gv = parent_cell_vertices[static_cast<std::size_t>(j)];
            ls_values[static_cast<std::size_t>(j)]
                = mesh.vertex_phi[static_cast<std::size_t>(gv * mesh.n_level_sets + level_set_id)];
            for (int d = 0; d < mesh.gdim; ++d)
            {
                cell_vertex_coords[static_cast<std::size_t>(j * mesh.gdim + d)]
                    = mesh.vertex_x[static_cast<std::size_t>(gv * mesh.gdim + d)];
            }
        }

        const cell::domain dom = cell::classify_cell_domain<T>(
            std::span<const T>(ls_values.data(), ls_values.size()));
        if (dom != cell::domain::intersected)
            continue;

        const int n_cut_edges = cell::num_edges(ct);
        std::vector<T> edge_root_coords(static_cast<std::size_t>(n_cut_edges * mesh.gdim), T(0));
        std::vector<uint8_t> edge_has_root(static_cast<std::size_t>(n_cut_edges), uint8_t(0));
        bool missing_root = false;
        for (int te = 0; te < n_cut_edges; ++te)
        {
            const int local_mesh_e = local_edge_for_cut_edge_token(
                mesh, parent_cell_vertices, parent_cell_edges, ct, te);
            if (local_mesh_e < 0 || local_mesh_e >= mesh.n_edges())
                continue;

            if (mesh.edge_state[static_cast<std::size_t>(local_mesh_e)]
                == static_cast<uint8_t>(EdgeState::one_root))
            {
                const int rv = mesh.edge_root_vertex[static_cast<std::size_t>(local_mesh_e)];
                if (rv < 0 || rv >= mesh.n_vertices())
                {
                    missing_root = true;
                    continue;
                }
                edge_has_root[static_cast<std::size_t>(te)] = 1;
                for (int d = 0; d < mesh.gdim; ++d)
                {
                    edge_root_coords[static_cast<std::size_t>(te * mesh.gdim + d)]
                        = mesh.vertex_x[static_cast<std::size_t>(rv * mesh.gdim + d)];
                }
            }
        }

        if (missing_root)
        {
            marked_cells[static_cast<std::size_t>(c)] = 1;
            ++invalid_cells;
            continue;
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

        const bool inside_ok = certify_cut(
            level_set, cut_inside, ct,
            std::span<const T>(ls_values.data(), ls_values.size()),
            mesh.parent_cell_id, tol, debug);
        const bool outside_ok = certify_cut(
            level_set, cut_outside, ct,
            std::span<const T>(ls_values.data(), ls_values.size()),
            mesh.parent_cell_id, tol, debug);
        if (!inside_ok || !outside_ok)
        {
            marked_cells[static_cast<std::size_t>(c)] = 1;
            ++invalid_cells;
        }
    }

    return invalid_cells > 0;
}

} // namespace

template <std::floating_point T, std::integral I>
bool certify_cut_subcells(
    const LevelSetFunction<T, I>& level_set,
    const cell::CutCell<T>&       cut_cell,
    cell::type                    parent_cell_type,
    std::span<const T>            parent_ls_values,
    int                           parent_cell_id,
    T                             tol,
    bool                          debug)
{
    return certify_cut(
        level_set, cut_cell, parent_cell_type,
        parent_ls_values, parent_cell_id, tol, debug);
}

template <std::floating_point T>
void decompose_local_mesh_from_cached_roots(
    LocalMesh<T>& mesh,
    int           level_set_id,
    bool          triangulate,
    bool          fill_missing_linear)
{
    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("decompose_local_mesh_linear: invalid level_set_id");
    if (static_cast<int>(mesh.vertex_root_edge_id.size()) != mesh.n_vertices())
        mesh.vertex_root_edge_id.assign(static_cast<std::size_t>(mesh.n_vertices()), -1);

    if (fill_missing_linear)
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
            if (mesh.edge_state[static_cast<std::size_t>(local_mesh_e)]
                != static_cast<uint8_t>(EdgeState::one_root))
            {
                continue;
            }
            int rv = mesh.edge_root_vertex[static_cast<std::size_t>(local_mesh_e)];
            if (rv < 0 && fill_missing_linear)
                rv = compute_edge_root_linear(mesh, local_mesh_e, level_set_id);
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

    deduplicate_zero_vertices(
        mesh, level_set_id, root_vertex_merge_tol<T>(), root_vertex_merge_tol<T>());
    build_local_edges(mesh);
    rebuild_parent_entity_maps(mesh);
}

template <std::floating_point T>
void decompose_local_mesh_linear(LocalMesh<T>& mesh, int level_set_id, bool triangulate)
{
    decompose_local_mesh_from_cached_roots(mesh, level_set_id, triangulate, true);
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

template <std::floating_point T, std::integral I>
LUTCertificationResult<T, I> decompose_local_mesh_with_backend(
    LocalMesh<T>&                 mesh,
    const LevelSetFunction<T, I>& level_set,
    LocalLevelSetBackend          backend,
    int                           level_set_id,
    cell::edge_root::method       root_method,
    bool                          triangulate,
    int                           max_refine_depth,
    T                             tol,
    bool                          debug)
{
    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("decompose_local_mesh_with_backend: invalid level_set_id");
    if (max_refine_depth < 0)
        throw std::invalid_argument("decompose_local_mesh_with_backend: max_refine_depth must be >= 0");

    LUTCertificationResult<T, I> result;
    if (static_cast<int>(mesh.vertex_root_edge_id.size()) != mesh.n_vertices())
        mesh.vertex_root_edge_id.assign(static_cast<std::size_t>(mesh.n_vertices()), -1);

    for (int depth = 0; ; ++depth)
    {
        if (backend == LocalLevelSetBackend::bernstein)
        {
            certify_local_mesh_polynomial<T, I>(
                mesh, level_set, level_set_id, 8, tol, true);
        }
        else
        {
            classify_edges_with_backend<T, I>(
                mesh, level_set, backend, level_set_id, tol, false);
        }

        compute_all_roots_with_backend<T, I>(
            mesh, level_set, backend, level_set_id, root_method, tol);

        std::vector<uint8_t> marked_cells;
        int invalid_cells = 0;
        const bool need_refine = certify_cut_on_local_mesh<T, I>(
            mesh, level_set, level_set_id, triangulate, tol,
            marked_cells, invalid_cells, debug);

        if (!need_refine)
        {
            decompose_local_mesh_from_cached_roots(mesh, level_set_id, triangulate, false);
            result.certified = true;
            result.invalid_cells = 0;
            result.refine_iterations = depth;
            return result;
        }

        result.invalid_cells = invalid_cells;
        result.refine_iterations = depth + 1;
        if (debug)
        {
            std::cerr << "rejecting LUT cut and refining local cells: depth=" << depth
                      << " invalid_cells=" << invalid_cells << "\n";
        }

        if (depth >= max_refine_depth)
        {
            decompose_local_mesh_linear(mesh, level_set_id, triangulate);
            result.certified = false;
            result.hit_max_depth = true;
            return result;
        }

        red_refine_marked_cells(
            mesh, std::span<const uint8_t>(marked_cells.data(), marked_cells.size()));
    }
}

// Explicit instantiations for common types
template void init_local_mesh_from_template<double>(LocalMesh<double>&, const RefinementTemplate&, std::span<const double>, cell::type, int, int);
template void init_local_mesh_from_template<float>(LocalMesh<float>&, const RefinementTemplate&, std::span<const float>, cell::type, int, int);
template void init_local_mesh_from_cell<double>(LocalMesh<double>&, std::span<const double>, cell::type, int, int);
template void init_local_mesh_from_cell<float>(LocalMesh<float>&, std::span<const float>, cell::type, int, int);
template void refine_local_mesh_from_template<double>(LocalMesh<double>&, const RefinementTemplate&);
template void refine_local_mesh_from_template<float>(LocalMesh<float>&, const RefinementTemplate&);
template void red_refine_marked_cells<double>(LocalMesh<double>&, std::span<const uint8_t>);
template void red_refine_marked_cells<float>(LocalMesh<float>&, std::span<const uint8_t>);
template void evaluate_levelsets_on_vertices<double>(LocalMesh<double>&, const std::vector<LevelSetFunction<double>>&, int, double);
template void evaluate_levelsets_on_vertices<float>(LocalMesh<float>&, const std::vector<LevelSetFunction<float>>&, int, float);
template void classify_local_edges<double>(LocalMesh<double>&, int);
template void classify_local_edges<float>(LocalMesh<float>&, int);
template int compute_edge_root_linear<double>(LocalMesh<double>&, int, int);
template int compute_edge_root_linear<float>(LocalMesh<float>&, int, int);
template int compute_edge_root_bernstein<double, int>(LocalMesh<double>&, const LocalLevelSetFunction<double, int>&, int, int, double);
template int compute_edge_root_bernstein<float, int>(LocalMesh<float>&, const LocalLevelSetFunction<float, int>&, int, int, float);
template int compute_edge_root<double, int>(LocalMesh<double>&, const LevelSetFunction<double, int>&, int, int, cell::edge_root::method);
template int compute_edge_root<float, int>(LocalMesh<float>&, const LevelSetFunction<float, int>&, int, int, cell::edge_root::method);
template void compute_all_roots_with_backend<double, int>(LocalMesh<double>&, const LevelSetFunction<double, int>&, LocalLevelSetBackend, int, cell::edge_root::method, double);
template void compute_all_roots_with_backend<float, int>(LocalMesh<float>&, const LevelSetFunction<float, int>&, LocalLevelSetBackend, int, cell::edge_root::method, float);
template bool certify_cut_subcells<double, int>(const LevelSetFunction<double, int>&, const cell::CutCell<double>&, cell::type, std::span<const double>, int, double, bool);
template bool certify_cut_subcells<float, int>(const LevelSetFunction<float, int>&, const cell::CutCell<float>&, cell::type, std::span<const float>, int, float, bool);
template void compute_all_roots_linear<double>(LocalMesh<double>&, int);
template void compute_all_roots_linear<float>(LocalMesh<float>&, int);
template void compute_all_roots<double, int>(LocalMesh<double>&, const LevelSetFunction<double, int>&, int, cell::edge_root::method);
template void compute_all_roots<float, int>(LocalMesh<float>&, const LevelSetFunction<float, int>&, int, cell::edge_root::method);
template void decompose_local_mesh_linear<double>(LocalMesh<double>&, int, bool);
template void decompose_local_mesh_linear<float>(LocalMesh<float>&, int, bool);
template void decompose_local_mesh<double, int>(LocalMesh<double>&, const LevelSetFunction<double, int>&, int, cell::edge_root::method, bool);
template void decompose_local_mesh<float, int>(LocalMesh<float>&, const LevelSetFunction<float, int>&, int, cell::edge_root::method, bool);
template LUTCertificationResult<double, int> decompose_local_mesh_with_backend<double, int>(LocalMesh<double>&, const LevelSetFunction<double, int>&, LocalLevelSetBackend, int, cell::edge_root::method, bool, int, double, bool);
template LUTCertificationResult<float, int> decompose_local_mesh_with_backend<float, int>(LocalMesh<float>&, const LevelSetFunction<float, int>&, LocalLevelSetBackend, int, cell::edge_root::method, bool, int, float, bool);

} // namespace cutcells
