// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "cut_mesh_curved.h"

#include "bernstein_backend.h"
#include "cell_flags.h"
#include "local_mesh.h"
#include "mapping.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{

struct SurfacePointKey
{
    int64_t x = 0;
    int64_t y = 0;
    int64_t z = 0;

    bool operator==(const SurfacePointKey& other) const noexcept
    {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct SurfacePointKeyHash
{
    std::size_t operator()(const SurfacePointKey& k) const noexcept
    {
        std::size_t h = std::hash<int64_t>{}(k.x);
        h ^= std::hash<int64_t>{}(k.y) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<int64_t>{}(k.z) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

template <std::floating_point T>
SurfacePointKey make_surface_point_key(const T* x_phys, int gdim, T tol)
{
    const T inv = T(1) / tol;
    SurfacePointKey key;
    key.x = static_cast<int64_t>(std::llround(x_phys[0] * inv));
    key.y = (gdim >= 2) ? static_cast<int64_t>(std::llround(x_phys[1] * inv)) : 0;
    key.z = (gdim >= 3) ? static_cast<int64_t>(std::llround(x_phys[2] * inv)) : 0;
    return key;
}

template <std::floating_point T>
int add_or_get_surface_point(cutcells::mesh::CutMesh<T>& out,
                             std::unordered_map<SurfacePointKey, int, SurfacePointKeyHash>& point_map,
                             const T* x_phys,
                             T tol)
{
    const auto key = make_surface_point_key(x_phys, out._gdim, tol);
    const auto it = point_map.find(key);
    if (it != point_map.end())
        return it->second;

    const int gv = static_cast<int>(
        out._vertex_coords.size() / static_cast<std::size_t>(out._gdim));
    for (int d = 0; d < out._gdim; ++d)
        out._vertex_coords.push_back(x_phys[static_cast<std::size_t>(d)]);
    point_map.emplace(key, gv);
    return gv;
}

struct BackgroundEdgeKey
{
    int32_t a = -1;
    int32_t b = -1;

    bool operator==(const BackgroundEdgeKey& other) const noexcept
    {
        return a == other.a && b == other.b;
    }
};

struct BackgroundEdgeKeyHash
{
    std::size_t operator()(const BackgroundEdgeKey& k) const noexcept
    {
        std::size_t h = std::hash<int32_t>{}(k.a);
        h ^= std::hash<int32_t>{}(k.b) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

struct BackgroundFaceKey
{
    std::array<int32_t, 4> verts = {-1, -1, -1, -1};
    uint8_t n = 0;

    bool operator==(const BackgroundFaceKey& other) const noexcept
    {
        return n == other.n && verts == other.verts;
    }
};

struct BackgroundFaceKeyHash
{
    std::size_t operator()(const BackgroundFaceKey& k) const noexcept
    {
        std::size_t h = std::hash<uint8_t>{}(k.n);
        for (int32_t v : k.verts)
            h ^= std::hash<int32_t>{}(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

struct BackgroundEntityMaps
{
    std::vector<int32_t> cell_vertex_offsets;
    std::vector<int32_t> cell_vertex_ids;
    std::vector<int32_t> cell_edge_offsets;
    std::vector<int32_t> cell_edge_ids;
    std::vector<int32_t> cell_face_offsets;
    std::vector<int32_t> cell_face_ids;
};

struct SurfaceParentKey
{
    std::array<int64_t, 8> data = {
        -1, -1, -1, -1, -1, -1, -1, -1};

    bool operator==(const SurfaceParentKey& other) const noexcept
    {
        return data == other.data;
    }
};

struct SurfaceParentKeyHash
{
    std::size_t operator()(const SurfaceParentKey& k) const noexcept
    {
        std::size_t h = 0;
        for (const int64_t v : k.data)
            h ^= std::hash<int64_t>{}(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

struct SurfacePointOrigin
{
    int parent_cell_id = -1;
    int zero_entity_id = -1;
    int tpl_vertex_id = -1;
    int tpl_parent_dim = -1;
    int tpl_parent_id = -1;
};

inline void log_surface_parent_mismatch(
    const char*               label,
    const SurfaceParentKey&   key,
    const SurfacePointOrigin& old_origin,
    const SurfacePointOrigin& new_origin,
    double                    max_diff)
{
    static int printed = 0;
    if (printed >= 32)
        return;
    ++printed;

    std::cerr << "Warning: " << label << "; max_diff=" << max_diff
              << ", key_kind=" << key.data[0]
              << ", key_data=(" << key.data[1]
              << "," << key.data[2]
              << "," << key.data[3]
              << "," << key.data[4]
              << "," << key.data[5]
              << "," << key.data[6]
              << "," << key.data[7]
              << ")"
              << ", first_origin=(cell=" << old_origin.parent_cell_id
              << ", ze=" << old_origin.zero_entity_id
              << ", tpl_vi=" << old_origin.tpl_vertex_id
              << ", tpl_parent_dim=" << old_origin.tpl_parent_dim
              << ", tpl_parent_id=" << old_origin.tpl_parent_id
              << ")"
              << ", new_origin=(cell=" << new_origin.parent_cell_id
              << ", ze=" << new_origin.zero_entity_id
              << ", tpl_vi=" << new_origin.tpl_vertex_id
              << ", tpl_parent_dim=" << new_origin.tpl_parent_dim
              << ", tpl_parent_id=" << new_origin.tpl_parent_id
              << ")\n";
}

inline int64_t pack_parent_ref(int32_t dim, int32_t id)
{
    const uint64_t hi = static_cast<uint32_t>(dim + 2);
    const uint64_t lo = static_cast<uint32_t>(id + 1);
    return static_cast<int64_t>((hi << 32) | lo);
}

template <std::floating_point T, std::integral I>
BackgroundEntityMaps build_background_entity_maps(const cutcells::MeshView<T, I>& mesh)
{
    BackgroundEntityMaps maps;
    const int nc = static_cast<int>(mesh.num_cells());
    maps.cell_vertex_offsets.reserve(static_cast<std::size_t>(nc + 1));
    maps.cell_edge_offsets.reserve(static_cast<std::size_t>(nc + 1));
    maps.cell_face_offsets.reserve(static_cast<std::size_t>(nc + 1));
    maps.cell_vertex_offsets.push_back(0);
    maps.cell_edge_offsets.push_back(0);
    maps.cell_face_offsets.push_back(0);

    std::unordered_map<BackgroundEdgeKey, int32_t, BackgroundEdgeKeyHash> edge_map;
    std::unordered_map<BackgroundFaceKey, int32_t, BackgroundFaceKeyHash> face_map;

    for (int c = 0; c < nc; ++c)
    {
        const auto ct = cutcells::cell::map_vtk_type_to_cell_type(
            static_cast<cutcells::cell::vtk_types>(mesh.cell_type(c)));
        const auto cell_nodes = mesh.cell_nodes(c);
        const int p1_nv = cutcells::cell::get_num_vertices(ct);
        if (static_cast<int>(cell_nodes.size()) < p1_nv)
            throw std::runtime_error(
                "build_background_entity_maps: cell has fewer than P1 vertices");

        for (int i = 0; i < p1_nv; ++i)
            maps.cell_vertex_ids.push_back(static_cast<int32_t>(cell_nodes[static_cast<std::size_t>(i)]));
        maps.cell_vertex_offsets.push_back(
            static_cast<int32_t>(maps.cell_vertex_ids.size()));

        const auto edges = cutcells::cell::edges(ct);
        for (const auto& edge : edges)
        {
            int32_t a = static_cast<int32_t>(cell_nodes[static_cast<std::size_t>(edge[0])]);
            int32_t b = static_cast<int32_t>(cell_nodes[static_cast<std::size_t>(edge[1])]);
            if (b < a)
                std::swap(a, b);
            const BackgroundEdgeKey key{a, b};
            auto [it, inserted] = edge_map.emplace(
                key, static_cast<int32_t>(edge_map.size()));
            maps.cell_edge_ids.push_back(it->second);
        }
        maps.cell_edge_offsets.push_back(
            static_cast<int32_t>(maps.cell_edge_ids.size()));

        const auto faces = cutcells::cell::faces(ct);
        const auto face_sizes = cutcells::cell::face_vertex_counts(ct);
        for (std::size_t fi = 0; fi < faces.size(); ++fi)
        {
            const auto& face = faces[fi];
            BackgroundFaceKey key;
            key.n = static_cast<uint8_t>(face_sizes[fi]);
            for (std::size_t i = 0; i < static_cast<std::size_t>(key.n); ++i)
                key.verts[i] = static_cast<int32_t>(
                    cell_nodes[static_cast<std::size_t>(face[i])]);
            std::sort(
                key.verts.begin(),
                key.verts.begin() + static_cast<std::ptrdiff_t>(key.n));
            auto [it, inserted] = face_map.emplace(
                key, static_cast<int32_t>(face_map.size()));
            maps.cell_face_ids.push_back(it->second);
        }
        maps.cell_face_offsets.push_back(
            static_cast<int32_t>(maps.cell_face_ids.size()));
    }

    return maps;
}

inline int32_t background_cell_vertex_id(
    const BackgroundEntityMaps& maps,
    int                         cell_id,
    int                         local_vertex_id)
{
    const int32_t off = maps.cell_vertex_offsets[static_cast<std::size_t>(cell_id)];
    return maps.cell_vertex_ids[static_cast<std::size_t>(off + local_vertex_id)];
}

inline int32_t background_cell_edge_id(
    const BackgroundEntityMaps& maps,
    int                         cell_id,
    int                         local_edge_id)
{
    const int32_t off = maps.cell_edge_offsets[static_cast<std::size_t>(cell_id)];
    return maps.cell_edge_ids[static_cast<std::size_t>(off + local_edge_id)];
}

inline int32_t background_cell_face_id(
    const BackgroundEntityMaps& maps,
    int                         cell_id,
    int                         local_face_id)
{
    const int32_t off = maps.cell_face_offsets[static_cast<std::size_t>(cell_id)];
    return maps.cell_face_ids[static_cast<std::size_t>(off + local_face_id)];
}

template <std::floating_point T>
std::pair<int32_t, int32_t> globalize_parent_ref(
    const BackgroundEntityMaps& maps,
    const cutcells::LocalMesh<T>& mesh,
    int32_t                     parent_dim,
    int32_t                     parent_id)
{
    if (parent_dim == 0 && parent_id >= 0)
        return {0, background_cell_vertex_id(maps, mesh.parent_cell_id, parent_id)};
    if (parent_dim == 1 && parent_id >= 0)
        return {1, background_cell_edge_id(maps, mesh.parent_cell_id, parent_id)};
    if (parent_dim == 2 && parent_id >= 0)
        return {2, background_cell_face_id(maps, mesh.parent_cell_id, parent_id)};
    if (parent_dim == 3)
        return {3, mesh.parent_cell_id};
    return {-1, -1};
}

template <std::floating_point T>
std::pair<int32_t, int32_t> globalize_vertex_parent_ref(
    const BackgroundEntityMaps& maps,
    const cutcells::LocalMesh<T>& mesh,
    int32_t                     lv)
{
    return globalize_parent_ref(
        maps,
        mesh,
        mesh.vertex_parent_dim[static_cast<std::size_t>(lv)],
        mesh.vertex_parent_id[static_cast<std::size_t>(lv)]);
}

template <std::floating_point T>
SurfaceParentKey make_zero_face_template_node_key(
    const BackgroundEntityMaps& maps,
    const cutcells::LocalMesh<T>& mesh,
    int                          zero_face_id,
    const cutcells::RefinementTemplate& tpl,
    int                          tpl_vertex_id,
    const std::array<int, 4>&    edge_node_totals,
    std::array<int, 4>&          edge_node_cursor,
    int&                         face_interior_id)
{
    SurfaceParentKey key;
    const int parent_dim = tpl.vertex_parent_dim[static_cast<std::size_t>(tpl_vertex_id)];
    const int parent_id = tpl.vertex_parent_id[static_cast<std::size_t>(tpl_vertex_id)];
    const int z0 = mesh.zero_entity_offsets[static_cast<std::size_t>(zero_face_id)];

    if (parent_dim == 0 && parent_id >= 0)
    {
        const int32_t lv = mesh.zero_entity_vertices[static_cast<std::size_t>(z0 + parent_id)];
        const auto [gdim, gid] = globalize_vertex_parent_ref(maps, mesh, lv);
        if (gdim == 2)
        {
            const int z1 = mesh.zero_entity_offsets[static_cast<std::size_t>(zero_face_id + 1)];
            const int n_face_verts = z1 - z0;
            const int prev_local = (parent_id + n_face_verts - 1) % n_face_verts;
            const int next_local = (parent_id + 1) % n_face_verts;
            const int32_t prev_lv = mesh.zero_entity_vertices[static_cast<std::size_t>(z0 + prev_local)];
            const int32_t next_lv = mesh.zero_entity_vertices[static_cast<std::size_t>(z0 + next_local)];
            auto [pd0, pid0] = globalize_vertex_parent_ref(maps, mesh, prev_lv);
            auto [pd1, pid1] = globalize_vertex_parent_ref(maps, mesh, next_lv);
            int64_t p0 = pack_parent_ref(pd0, pid0);
            int64_t p1 = pack_parent_ref(pd1, pid1);
            if (p1 < p0)
                std::swap(p0, p1);

            key.data[0] = 3;
            key.data[1] = gdim;
            key.data[2] = gid;
            key.data[3] = p0;
            key.data[4] = p1;
        }
        else
        {
            key.data[0] = 0;
            key.data[1] = gdim;
            key.data[2] = gid;
        }
        return key;
    }

    if (parent_dim == 1 && parent_id >= 0)
    {
        int edge_ze = -1;
        bool reversed = false;
        if (!cutcells::locate_zero_face_boundary_edge_entity(
                mesh, zero_face_id, parent_id, &edge_ze, &reversed)
            || edge_ze < 0)
        {
            key.data[0] = 9;
            key.data[1] = mesh.parent_cell_id;
            key.data[2] = zero_face_id;
            key.data[3] = tpl_vertex_id;
            return key;
        }

        int local_idx = edge_node_cursor[static_cast<std::size_t>(parent_id)]++;
        const int n_edge_nodes = edge_node_totals[static_cast<std::size_t>(parent_id)];
        if (reversed)
            local_idx = n_edge_nodes - 1 - local_idx;

        const int32_t ev0 = mesh.zero_entity_endpoint_v0[static_cast<std::size_t>(edge_ze)];
        const int32_t ev1 = mesh.zero_entity_endpoint_v1[static_cast<std::size_t>(edge_ze)];
        const auto [d0, id0] = globalize_vertex_parent_ref(maps, mesh, ev0);
        const auto [d1, id1] = globalize_vertex_parent_ref(maps, mesh, ev1);
        int64_t p0 = pack_parent_ref(d0, id0);
        int64_t p1 = pack_parent_ref(d1, id1);
        if (p1 < p0)
        {
            std::swap(p0, p1);
            local_idx = n_edge_nodes - 1 - local_idx;
        }

        const auto [carrier_dim, carrier_id] = globalize_parent_ref(
            maps,
            mesh,
            mesh.zero_entity_parent_dim[static_cast<std::size_t>(edge_ze)],
            mesh.zero_entity_parent_id[static_cast<std::size_t>(edge_ze)]);
        key.data[0] = 1;
        key.data[1] = carrier_dim;
        key.data[2] = carrier_id;
        key.data[3] = p0;
        key.data[4] = p1;
        key.data[5] = local_idx;
        // Temporary workaround: do not merge curved edge nodes across parent
        // cells until shared-face edge curving is made consistent.
        key.data[6] = mesh.parent_cell_id;
        return key;
    }

    if (parent_dim == 2)
    {
        const auto [carrier_dim, carrier_id] = globalize_parent_ref(
            maps,
            mesh,
            mesh.zero_entity_parent_dim[static_cast<std::size_t>(zero_face_id)],
            mesh.zero_entity_parent_id[static_cast<std::size_t>(zero_face_id)]);
        const int z1 = mesh.zero_entity_offsets[static_cast<std::size_t>(zero_face_id + 1)];
        std::array<int64_t, 4> corner_keys = {-1, -1, -1, -1};
        for (int i = 0; i < z1 - z0 && i < 4; ++i)
        {
            const int32_t lv = mesh.zero_entity_vertices[static_cast<std::size_t>(z0 + i)];
            const auto [cdim, cid] = globalize_vertex_parent_ref(maps, mesh, lv);
            corner_keys[static_cast<std::size_t>(i)] = pack_parent_ref(cdim, cid);
        }
        std::sort(corner_keys.begin(), corner_keys.end());
        key.data[0] = 2;
        key.data[1] = carrier_dim;
        key.data[2] = carrier_id;
        key.data[3] = corner_keys[0];
        key.data[4] = corner_keys[1];
        key.data[5] = corner_keys[2];
        key.data[6] = corner_keys[3];
        key.data[7] = face_interior_id++;
        return key;
    }

    key.data[0] = 9;
    key.data[1] = mesh.parent_cell_id;
    key.data[2] = zero_face_id;
    key.data[3] = tpl_vertex_id;
    return key;
}

template <std::floating_point T>
int32_t add_or_get_surface_point(
    cutcells::mesh::CurvedGlobalMesh<T>& out,
    std::unordered_map<SurfaceParentKey, int32_t, SurfaceParentKeyHash>& point_map,
    std::vector<SurfacePointOrigin>& origins,
    const SurfaceParentKey& key,
    const SurfacePointOrigin& origin,
    const T* x_phys,
    T tol)
{
    const auto it = point_map.find(key);
    if (it != point_map.end())
    {
        const int32_t gv = it->second;
        const T* x_old = &out.vertex_coords[static_cast<std::size_t>(gv * out.gdim)];
        T max_diff = T(0);
        for (int d = 0; d < out.gdim; ++d)
            max_diff = std::max(max_diff, std::abs(x_old[d] - x_phys[d]));
        if (max_diff > tol)
        {
            log_surface_parent_mismatch(
                "parent-mapped interface node mismatch on shared carrier",
                key,
                origins[static_cast<std::size_t>(gv)],
                origin,
                static_cast<double>(max_diff));
        }
        return gv;
    }

    const int32_t gv = static_cast<int32_t>(out.vertex_coords.size() / static_cast<std::size_t>(out.gdim));
    for (int d = 0; d < out.gdim; ++d)
        out.vertex_coords.push_back(x_phys[static_cast<std::size_t>(d)]);
    origins.push_back(origin);
    point_map.emplace(key, gv);
    return gv;
}

template <std::floating_point T>
int add_or_get_surface_point(
    cutcells::mesh::CutMesh<T>& out,
    std::unordered_map<SurfaceParentKey, int, SurfaceParentKeyHash>& point_map,
    std::vector<SurfacePointOrigin>& origins,
    const SurfaceParentKey& key,
    const SurfacePointOrigin& origin,
    const T* x_phys,
    T tol)
{
    const auto it = point_map.find(key);
    if (it != point_map.end())
    {
        const int gv = it->second;
        const T* x_old = &out._vertex_coords[static_cast<std::size_t>(gv * out._gdim)];
        T max_diff = T(0);
        for (int d = 0; d < out._gdim; ++d)
            max_diff = std::max(max_diff, std::abs(x_old[d] - x_phys[d]));
        if (max_diff > tol)
        {
            log_surface_parent_mismatch(
                "parent-mapped sampled interface node mismatch on shared carrier",
                key,
                origins[static_cast<std::size_t>(gv)],
                origin,
                static_cast<double>(max_diff));
        }
        return gv;
    }

    const int gv = static_cast<int>(
        out._vertex_coords.size() / static_cast<std::size_t>(out._gdim));
    for (int d = 0; d < out._gdim; ++d)
        out._vertex_coords.push_back(x_phys[static_cast<std::size_t>(d)]);
    origins.push_back(origin);
    point_map.emplace(key, gv);
    return gv;
}

template <std::floating_point T>
void init_curved_mesh(cutcells::mesh::CurvedGlobalMesh<T>& mesh,
                      int gdim,
                      int geom_order,
                      int level_set_id)
{
    mesh.vertex_coords.clear();
    mesh.connectivity.clear();
    mesh.offsets.clear();
    mesh.cell_types.clear();
    mesh.gdim = gdim;
    mesh.geom_order = geom_order;
    mesh.level_set_id = level_set_id;
    mesh.n_fallback_nodes = 0;
    mesh.offsets.push_back(0);
}

inline std::string to_lower_ascii(std::string_view text)
{
    std::string out(text);
    std::transform(out.begin(), out.end(), out.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return out;
}

inline bool bernstein_supported_cell(cutcells::cell::type ct)
{
    using cutcells::cell::type;
    return ct == type::interval || ct == type::triangle
           || ct == type::tetrahedron || ct == type::hexahedron;
}

template <std::floating_point T>
void ref_to_phys_affine_embedded(cutcells::cell::type ct,
                                 std::span<const T> parent_coords,
                                 int gdim,
                                 std::span<const T> x_ref,
                                 T* x_phys)
{
    if (gdim <= 0 || gdim > 3)
        throw std::invalid_argument("ref_to_phys_affine_embedded: gdim must be in [1, 3]");

    const int tdim = cutcells::cell::get_tdim(ct);
    if (static_cast<int>(x_ref.size()) != tdim)
        throw std::invalid_argument("ref_to_phys_affine_embedded: reference point has invalid dimension");

    const int p1_vertices = cutcells::cell::get_num_vertices(ct);
    if (static_cast<int>(parent_coords.size()) != p1_vertices * gdim)
        throw std::invalid_argument("ref_to_phys_affine_embedded: invalid parent_coords size");

    const T* x0 = parent_coords.data();
    for (int d = 0; d < gdim; ++d)
        x_phys[d] = x0[d];

    if (tdim == 0)
        return;

    const auto cols = cutcells::cell::jacobian_col_indices(ct);
    for (int k = 0; k < tdim; ++k)
    {
        const int vk = cols[static_cast<std::size_t>(k)];
        if (vk < 0 || vk >= p1_vertices)
            throw std::runtime_error("ref_to_phys_affine_embedded: invalid Jacobian column index");

        const T* xk = parent_coords.data() + static_cast<std::size_t>(vk * gdim);
        const T xi = x_ref[static_cast<std::size_t>(k)];
        for (int d = 0; d < gdim; ++d)
            x_phys[d] += xi * (xk[d] - x0[d]);
    }
}

template <std::floating_point T>
void ref_to_phys_single_local(const cutcells::LocalMesh<T>& mesh,
                              const T*                      x_ref,
                              T*                            x_phys)
{
    ref_to_phys_affine_embedded(
        mesh.parent_cell_type,
        std::span<const T>(
            mesh.parent_cell_coords_p1.data(),
            mesh.parent_cell_coords_p1.size()),
        mesh.gdim,
        std::span<const T>(x_ref, static_cast<std::size_t>(mesh.tdim)),
        x_phys);
}

template <std::floating_point T>
inline cutcells::cell::type zero_face_type_local(
    const cutcells::LocalMesh<T>& mesh,
    int                           zero_entity_id)
{
    const int z0 = mesh.zero_entity_offsets[static_cast<std::size_t>(zero_entity_id)];
    const int z1 = mesh.zero_entity_offsets[static_cast<std::size_t>(zero_entity_id + 1)];
    const int n_face_verts = z1 - z0;
    if (n_face_verts == 3)
        return cutcells::cell::type::triangle;
    if (n_face_verts == 4)
        return cutcells::cell::type::quadrilateral;
    return cutcells::cell::type::point;
}

template <std::floating_point T>
inline void build_zero_face_coordinate_polys_local(
    const cutcells::LocalMesh<T>&            mesh,
    int                                      zero_entity_id,
    int                                      geom_order,
    std::vector<cutcells::BernsteinCell<T>>& coord_polys)
{
    const auto face_type = zero_face_type_local(mesh, zero_entity_id);
    if (face_type != cutcells::cell::type::triangle
        && face_type != cutcells::cell::type::quadrilateral)
        throw std::invalid_argument(
            "build_zero_face_coordinate_polys_local: unsupported zero-face type");

    const int tdim = mesh.tdim;
    std::vector<T> face_nodes_ref;
    const int degree = cutcells::build_zero_face_nodes_ref(
        mesh, zero_entity_id, geom_order, face_nodes_ref);
    const auto& tpl = (degree == 1)
        ? cutcells::p1_template(face_type)
        : cutcells::iso_p1_template(face_type, degree);

    coord_polys.assign(static_cast<std::size_t>(tdim), cutcells::BernsteinCell<T>{});
    std::vector<T> lagrange_values(static_cast<std::size_t>(tpl.n_vertices), T(0));
    for (int d = 0; d < tdim; ++d)
    {
        for (int i = 0; i < tpl.n_vertices; ++i)
            lagrange_values[static_cast<std::size_t>(i)]
                = face_nodes_ref[static_cast<std::size_t>(i * tdim + d)];
        coord_polys[static_cast<std::size_t>(d)] = cutcells::make_bernstein_cell<T>(
            face_type, degree,
            std::span<const T>(lagrange_values.data(), lagrange_values.size()));
    }
}

template <std::floating_point T>
cutcells::mesh::CurvedGlobalMesh<T> assemble_curved_interface_mesh_parent_mapped(
    std::span<const cutcells::LocalMesh<T>*> local_meshes,
    const BackgroundEntityMaps&              bg_maps,
    int                                      level_set_id,
    int                                      geom_order)
{
    cutcells::mesh::CurvedGlobalMesh<T> out;
    out.geom_order = geom_order;
    out.level_set_id = level_set_id;
    out.offsets.push_back(0);

    if (local_meshes.empty())
        return out;

    out.gdim = local_meshes.front()->gdim;
    constexpr T surface_mismatch_tol = static_cast<T>(1e-10);
    std::unordered_map<SurfaceParentKey, int32_t, SurfaceParentKeyHash> point_map;
    std::vector<SurfacePointOrigin> point_origins;
    const uint64_t zero_mask = uint64_t(1) << level_set_id;

    for (const cutcells::LocalMesh<T>* mesh_ptr : local_meshes)
    {
        if (mesh_ptr == nullptr)
            continue;
        const auto& mesh = *mesh_ptr;
        if (mesh.tdim != 3 || mesh.n_zero_patches() == 0)
            continue;

        for (int patch = 0; patch < mesh.n_zero_patches(); ++patch)
        {
            if (mesh.zero_patch_zero_mask[static_cast<std::size_t>(patch)] != zero_mask)
                continue;

            const int p0 = mesh.zero_patch_offsets[static_cast<std::size_t>(patch)];
            const int p1 = mesh.zero_patch_offsets[static_cast<std::size_t>(patch + 1)];
            for (int pi = p0; pi < p1; ++pi)
            {
                const int ze = mesh.zero_patch_entity_ids[static_cast<std::size_t>(pi)];
                const auto face_type = zero_face_type_local(mesh, ze);
                std::vector<T> face_nodes_ref;
                const int degree = cutcells::build_zero_face_nodes_ref(
                    mesh, ze, geom_order, face_nodes_ref, &out.n_fallback_nodes);
                const auto& tpl = (degree == 1)
                    ? cutcells::p1_template(face_type)
                    : cutcells::iso_p1_template(face_type, degree);
                const auto& vtk_to_basix
                    = cutcells::vtk_lagrange_to_basix_permutation(face_type, degree);

                std::array<int, 4> edge_node_totals = {0, 0, 0, 0};
                for (int vi = 0; vi < tpl.n_vertices; ++vi)
                {
                    const int pdim = tpl.vertex_parent_dim[static_cast<std::size_t>(vi)];
                    const int pid = tpl.vertex_parent_id[static_cast<std::size_t>(vi)];
                    if (pdim == 1 && pid >= 0
                        && pid < static_cast<int>(edge_node_totals.size()))
                        ++edge_node_totals[static_cast<std::size_t>(pid)];
                }

                std::array<int, 4> edge_node_cursor = {0, 0, 0, 0};
                int face_interior_id = 0;
                std::array<T, 3> x_phys = {T(0), T(0), T(0)};
                std::vector<int32_t> point_ids(
                    static_cast<std::size_t>(tpl.n_vertices), int32_t(-1));

                for (int vi = 0; vi < tpl.n_vertices; ++vi)
                {
                    const SurfaceParentKey key = make_zero_face_template_node_key(
                        bg_maps, mesh, ze, tpl, vi, edge_node_totals,
                        edge_node_cursor, face_interior_id);
                    SurfacePointOrigin origin;
                    origin.parent_cell_id = mesh.parent_cell_id;
                    origin.zero_entity_id = ze;
                    origin.tpl_vertex_id = vi;
                    origin.tpl_parent_dim =
                        tpl.vertex_parent_dim[static_cast<std::size_t>(vi)];
                    origin.tpl_parent_id =
                        tpl.vertex_parent_id[static_cast<std::size_t>(vi)];
                    const T* x_ref = &face_nodes_ref[static_cast<std::size_t>(vi * mesh.tdim)];
                    ref_to_phys_single_local(mesh, x_ref, x_phys.data());
                    point_ids[static_cast<std::size_t>(vi)] = add_or_get_surface_point(
                        out,
                        point_map,
                        point_origins,
                        key,
                        origin,
                        x_phys.data(),
                        surface_mismatch_tol);
                }

                for (int vtk_vi = 0; vtk_vi < tpl.n_vertices; ++vtk_vi)
                {
                    const int vi = vtk_to_basix[static_cast<std::size_t>(vtk_vi)];
                    out.connectivity.push_back(point_ids[static_cast<std::size_t>(vi)]);
                }

                out.offsets.push_back(static_cast<int32_t>(out.connectivity.size()));
                out.cell_types.push_back(face_type);
            }
        }
    }

    return out;
}

template <std::floating_point T>
cutcells::mesh::CutMesh<T> sample_curved_local_interfaces_parent_mapped(
    std::span<const cutcells::LocalMesh<T>*> local_meshes,
    const BackgroundEntityMaps&              bg_maps,
    int                                      level_set_id,
    int                                      geom_order,
    int                                      subdivision)
{
    cutcells::mesh::CutMesh<T> out;
    out._offset.push_back(0);

    if (local_meshes.empty())
    {
        out._gdim = 0;
        out._tdim = 0;
        out._num_cells = 0;
        out._num_vertices = 0;
        return out;
    }

    out._gdim = local_meshes.front()->gdim;
    out._tdim = std::max(0, local_meshes.front()->tdim - 1);

    constexpr T surface_mismatch_tol = static_cast<T>(1e-10);
    std::unordered_map<SurfaceParentKey, int, SurfaceParentKeyHash> point_map;
    std::vector<SurfacePointOrigin> point_origins;
    const uint64_t zero_mask = uint64_t(1) << level_set_id;
    const int n_sub = std::max(1, subdivision);

    for (const cutcells::LocalMesh<T>* mesh_ptr : local_meshes)
    {
        if (mesh_ptr == nullptr)
            continue;
        const auto& mesh = *mesh_ptr;
        if (mesh.tdim != 3 || mesh.n_zero_entities() == 0)
            continue;

        const int n_zero = mesh.n_zero_entities();
        for (int ze = 0; ze < n_zero; ++ze)
        {
            if (mesh.zero_entity_dim[static_cast<std::size_t>(ze)] != 2)
                continue;
            if ((mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] & zero_mask) == 0)
                continue;

            const auto face_type = zero_face_type_local(mesh, ze);
            if (face_type != cutcells::cell::type::triangle
                && face_type != cutcells::cell::type::quadrilateral)
                continue;

            std::vector<cutcells::BernsteinCell<T>> face_coord_polys;
            build_zero_face_coordinate_polys_local(
                mesh, ze, geom_order, face_coord_polys);

            const auto& tpl = (n_sub <= 1)
                ? cutcells::p1_template(face_type)
                : cutcells::iso_p1_template(face_type, std::min(n_sub, 4));

            std::array<int, 4> edge_node_totals = {0, 0, 0, 0};
            for (int vi = 0; vi < tpl.n_vertices; ++vi)
            {
                const int pdim = tpl.vertex_parent_dim[static_cast<std::size_t>(vi)];
                const int pid = tpl.vertex_parent_id[static_cast<std::size_t>(vi)];
                if (pdim == 1 && pid >= 0
                    && pid < static_cast<int>(edge_node_totals.size()))
                    ++edge_node_totals[static_cast<std::size_t>(pid)];
            }

            std::array<int, 4> edge_node_cursor = {0, 0, 0, 0};
            int face_interior_id = 0;
            std::vector<int32_t> sample_node_ids(
                static_cast<std::size_t>(tpl.n_vertices), int32_t(-1));
            std::array<T, 3> x_ref = {T(0), T(0), T(0)};
            std::array<T, 3> x_phys = {T(0), T(0), T(0)};

            for (int vi = 0; vi < tpl.n_vertices; ++vi)
            {
                const SurfaceParentKey key = make_zero_face_template_node_key(
                    bg_maps, mesh, ze, tpl, vi, edge_node_totals,
                    edge_node_cursor, face_interior_id);
                SurfacePointOrigin origin;
                origin.parent_cell_id = mesh.parent_cell_id;
                origin.zero_entity_id = ze;
                origin.tpl_vertex_id = vi;
                origin.tpl_parent_dim =
                    tpl.vertex_parent_dim[static_cast<std::size_t>(vi)];
                origin.tpl_parent_id =
                    tpl.vertex_parent_id[static_cast<std::size_t>(vi)];
                const T u = static_cast<T>(
                    tpl.ref_vertex_coords[static_cast<std::size_t>(vi * tpl.tdim + 0)]);
                const T v = static_cast<T>(
                    tpl.ref_vertex_coords[static_cast<std::size_t>(vi * tpl.tdim + 1)]);
                const std::array<T, 2> uv = {u, v};
                for (int d = 0; d < mesh.tdim; ++d)
                {
                    x_ref[static_cast<std::size_t>(d)] =
                        cutcells::evaluate_bernstein_cell<T>(
                            face_coord_polys[static_cast<std::size_t>(d)],
                            std::span<const T>(uv.data(), uv.size()));
                }

                ref_to_phys_single_local(mesh, x_ref.data(), x_phys.data());
                sample_node_ids[static_cast<std::size_t>(vi)] = add_or_get_surface_point(
                    out,
                    point_map,
                    point_origins,
                    key,
                    origin,
                    x_phys.data(),
                    surface_mismatch_tol);
            }

            for (int child = 0; child < tpl.n_cells; ++child)
            {
                const int base = child * tpl.vertices_per_cell;
                for (int j = 0; j < tpl.vertices_per_cell; ++j)
                {
                    const int lv = tpl.cell_connectivity[
                        static_cast<std::size_t>(base + j)];
                    out._connectivity.push_back(
                        sample_node_ids[static_cast<std::size_t>(lv)]);
                }
                out._offset.push_back(static_cast<int>(out._connectivity.size()));
                out._types.push_back(tpl.child_cell_type);
                out._parent_map.push_back(-1);
            }
        }
    }

    out._num_cells = static_cast<int>(out._types.size());
    out._num_vertices = static_cast<int>(
        out._vertex_coords.size() / static_cast<std::size_t>(out._gdim));
    return out;
}

enum class ParentCellDomain : uint8_t
{
    inside = 0,
    intersected = 1,
    outside = 2
};

template <std::floating_point T>
ParentCellDomain classify_parent_cell_bernstein(cutcells::cell::type ct,
                                                int degree,
                                                std::span<const T> nodal_values,
                                                T tol)
{
    const auto poly = cutcells::make_bernstein_cell<T>(ct, degree, nodal_values);

    T scale = T(0);
    for (const T c : poly.coeffs)
        scale = std::max(scale, std::abs(c));
    const T tol_eff = std::max(
        tol,
        std::numeric_limits<T>::epsilon() * std::max(T(1), scale) * T(32));

    bool all_neg = true;
    bool all_pos = true;
    for (const T c : poly.coeffs)
    {
        all_neg = all_neg && (c < -tol_eff);
        all_pos = all_pos && (c > tol_eff);
    }

    if (all_neg)
        return ParentCellDomain::inside;
    if (all_pos)
        return ParentCellDomain::outside;
    return ParentCellDomain::intersected;
}

template <std::floating_point T>
void sample_callable_on_parent_nodes(
    cutcells::cell::type                 ct,
    std::span<const T>                   parent_coords_p1,
    int                                  gdim,
    int                                  degree,
    const cutcells::LevelSetFunction<T>& level_set,
    int                                  cell_id,
    std::vector<T>&                      sampled_node_coords,
    std::vector<T>&                      sampled_nodal_values)
{
    const int tdim = cutcells::cell::get_tdim(ct);
    const int n_nodes = cutcells::lagrange_node_count(ct, degree);
    const auto ref_nodes = cutcells::bernstein_reference_nodes(ct, degree);
    if (static_cast<int>(ref_nodes.size()) != n_nodes * tdim)
        throw std::runtime_error("sample_callable_on_parent_nodes: unexpected reference-node count");

    sampled_node_coords.resize(static_cast<std::size_t>(n_nodes * gdim));
    sampled_nodal_values.resize(static_cast<std::size_t>(n_nodes));

    std::array<T, 3> x_phys = {T(0), T(0), T(0)};
    std::array<T, 3> x_query = {T(0), T(0), T(0)};

    for (int i = 0; i < n_nodes; ++i)
    {
        std::array<T, 3> x_ref = {T(0), T(0), T(0)};
        for (int d = 0; d < tdim; ++d)
        {
            x_ref[static_cast<std::size_t>(d)] = static_cast<T>(
                ref_nodes[static_cast<std::size_t>(i * tdim + d)]);
        }

        ref_to_phys_affine_embedded<T>(
            ct,
            parent_coords_p1,
            gdim,
            std::span<const T>(x_ref.data(), static_cast<std::size_t>(tdim)),
            x_phys.data());

        for (int d = 0; d < gdim; ++d)
        {
            const T value = x_phys[static_cast<std::size_t>(d)];
            sampled_node_coords[static_cast<std::size_t>(i * gdim + d)] = value;
            x_query[static_cast<std::size_t>(d)] = value;
        }
        for (int d = gdim; d < 3; ++d)
            x_query[static_cast<std::size_t>(d)] = T(0);

        sampled_nodal_values[static_cast<std::size_t>(i)] = level_set.value(
            x_query.data(), static_cast<int>(cell_id));
    }
}

template <std::floating_point T>
void append_cell_to_curved_mesh(cutcells::mesh::CurvedGlobalMesh<T>& mesh,
                                cutcells::cell::type                 ct,
                                std::span<const T>                   coords,
                                int                                  gdim)
{
    if (mesh.offsets.empty())
        mesh.offsets.push_back(0);

    const int n_vertices = static_cast<int>(coords.size()) / gdim;
    const int32_t v0 = static_cast<int32_t>(mesh.vertex_coords.size() / static_cast<std::size_t>(gdim));

    mesh.vertex_coords.insert(mesh.vertex_coords.end(), coords.begin(), coords.end());
    for (int i = 0; i < n_vertices; ++i)
        mesh.connectivity.push_back(v0 + static_cast<int32_t>(i));

    mesh.offsets.push_back(static_cast<int32_t>(mesh.connectivity.size()));
    mesh.cell_types.push_back(ct);
}

template <std::floating_point T>
void append_local_cells_by_domain(const cutcells::LocalMesh<T>&      local_mesh,
                                  cutcells::cell::domain             domain,
                                  cutcells::mesh::CurvedGlobalMesh<T>& out,
                                  int                                  geom_order)
{
    using cutcells::mapping::eval_gh_tet_one_curved_tri_face;
    using cutcells::mapping::eval_gh_triangle_one_curved_edge;

    const int gdim = local_mesh.gdim;
    const int tdim = local_mesh.tdim;
    const int nc = local_mesh.n_cells();
    const bool has_curved = local_mesh.curved_geometry_order >= 2
        && !local_mesh.curved_zero_offsets.empty();
    const uint64_t zero_mask = uint64_t(1) << 0;
    const int n_zero = local_mesh.n_zero_entities();
    const int codim1_dim = std::max(0, tdim - 1);

    for (int c = 0; c < nc; ++c)
    {
        if (local_mesh.cell_domain[static_cast<std::size_t>(c)] != static_cast<uint8_t>(domain))
            continue;

        const int c0 = local_mesh.cell_offsets[static_cast<std::size_t>(c)];
        const int c1 = local_mesh.cell_offsets[static_cast<std::size_t>(c + 1)];
        const int n_vertices = c1 - c0;
        const auto ctype = local_mesh.cell_types[static_cast<std::size_t>(c)];

        int touching_ze = -1;
        int cell_if_a = -1;
        int cell_if_b = -1;
        bool iface_reversed = false;
        std::vector<int32_t> touching_face_verts;

        if (has_curved && geom_order >= 2
            && tdim == 2 && ctype == cutcells::cell::type::triangle && n_vertices == 3)
        {
            for (int ze = 0; ze < n_zero && touching_ze < 0; ++ze)
            {
                if (local_mesh.zero_entity_dim[static_cast<std::size_t>(ze)] != codim1_dim)
                    continue;
                if ((local_mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] & zero_mask) == 0)
                    continue;

                const int if_a = local_mesh.zero_entity_endpoint_v0[static_cast<std::size_t>(ze)];
                const int if_b = local_mesh.zero_entity_endpoint_v1[static_cast<std::size_t>(ze)];
                if (if_a < 0 || if_b < 0)
                    continue;

                for (int e = 0; e < n_vertices; ++e)
                {
                    const int ea = local_mesh.cell_vertices[static_cast<std::size_t>(c0 + e)];
                    const int eb = local_mesh.cell_vertices[static_cast<std::size_t>(c0 + (e + 1) % n_vertices)];
                    if ((ea == if_a && eb == if_b) || (ea == if_b && eb == if_a))
                    {
                        touching_ze = ze;
                        cell_if_a = ea;
                        cell_if_b = eb;
                        iface_reversed = (ea == if_b && eb == if_a);
                        break;
                    }
                }
            }
        }
        else if (has_curved && geom_order >= 2
                 && tdim == 3 && ctype == cutcells::cell::type::tetrahedron && n_vertices == 4)
        {
            std::vector<int32_t> cell_verts(static_cast<std::size_t>(n_vertices));
            for (int j = 0; j < n_vertices; ++j)
                cell_verts[static_cast<std::size_t>(j)]
                    = local_mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];

            for (int ze = 0; ze < n_zero && touching_ze < 0; ++ze)
            {
                if (local_mesh.zero_entity_dim[static_cast<std::size_t>(ze)] != codim1_dim)
                    continue;
                if ((local_mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] & zero_mask) == 0)
                    continue;

                const int z0 = local_mesh.zero_entity_offsets[static_cast<std::size_t>(ze)];
                const int z1 = local_mesh.zero_entity_offsets[static_cast<std::size_t>(ze + 1)];
                const int n_fv = z1 - z0;
                if (n_fv != 3)
                    continue;

                bool all_in_cell = true;
                std::vector<int32_t> fv(static_cast<std::size_t>(n_fv));
                for (int i = 0; i < n_fv; ++i)
                {
                    fv[static_cast<std::size_t>(i)]
                        = local_mesh.zero_entity_vertices[static_cast<std::size_t>(z0 + i)];
                    if (std::find(cell_verts.begin(), cell_verts.end(), fv[static_cast<std::size_t>(i)])
                        == cell_verts.end())
                    {
                        all_in_cell = false;
                        break;
                    }
                }
                if (all_in_cell)
                {
                    touching_ze = ze;
                    touching_face_verts = std::move(fv);
                }
            }
        }

        if (touching_ze >= 0
            && tdim == 2 && ctype == cutcells::cell::type::triangle && geom_order >= 2)
        {
            int opp_idx = -1;
            for (int j = 0; j < n_vertices; ++j)
            {
                const int lv = local_mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];
                if (lv != cell_if_a && lv != cell_if_b)
                {
                    opp_idx = j;
                    break;
                }
            }
            if (opp_idx >= 0)
            {
                const int lv_opp = local_mesh.cell_vertices[static_cast<std::size_t>(c0 + opp_idx)];
                const T* r_va = &local_mesh.vertex_ref_x[static_cast<std::size_t>(cell_if_a * tdim)];
                const T* r_vb = &local_mesh.vertex_ref_x[static_cast<std::size_t>(cell_if_b * tdim)];
                const T* r_opp = &local_mesh.vertex_ref_x[static_cast<std::size_t>(lv_opp * tdim)];

                const int node_start = local_mesh.curved_zero_offsets[static_cast<std::size_t>(touching_ze)];
                const int node_end = local_mesh.curved_zero_offsets[static_cast<std::size_t>(touching_ze + 1)];
                const int n_int = node_end - node_start;
                std::vector<T> x_int_ref(static_cast<std::size_t>(n_int * tdim), T(0));
                for (int ni = 0; ni < n_int; ++ni)
                {
                    const int src_ni = iface_reversed ? (n_int - 1 - ni) : ni;
                    const int ci_node = node_start + src_ni;
                    const T* xr = &local_mesh.curved_zero_ref_nodes[
                        static_cast<std::size_t>(ci_node * tdim)];
                    T* rp = &x_int_ref[static_cast<std::size_t>(ni * tdim)];
                    for (int d = 0; d < tdim; ++d)
                        rp[d] = xr[d];
                }

                const auto& tpl = cutcells::iso_p1_template(cutcells::cell::type::triangle, geom_order);
                const auto gl_pts_vec = ::cutcells::detail::gauss_lobatto_interior_points_1d<T>(geom_order);
                const T* gl_pts = gl_pts_vec.data();
                const int32_t v0 = static_cast<int32_t>(out.vertex_coords.size() / static_cast<std::size_t>(gdim));
                for (int vi = 0; vi < tpl.n_vertices; ++vi)
                {
                    const T xi = static_cast<T>(
                        tpl.ref_vertex_coords[static_cast<std::size_t>(vi * tpl.tdim + 0)]);
                    const T eta = static_cast<T>(
                        tpl.ref_vertex_coords[static_cast<std::size_t>(vi * tpl.tdim + 1)]);

                    T F_ref[3] = {};
                    T dF_dxi[3] = {};
                    T dF_deta[3] = {};
                    eval_gh_triangle_one_curved_edge(
                        xi, eta,
                        r_va, r_vb,
                        x_int_ref.data(), gl_pts, local_mesh.curved_geometry_order,
                        r_opp, tdim,
                        F_ref, dF_dxi, dF_deta);

                    T x_phys[3] = {};
                    ref_to_phys_single_local(local_mesh, F_ref, x_phys);
                    for (int d = 0; d < gdim; ++d)
                        out.vertex_coords.push_back(x_phys[d]);
                    out.connectivity.push_back(v0 + static_cast<int32_t>(vi));
                }
                out.offsets.push_back(static_cast<int32_t>(out.connectivity.size()));
                out.cell_types.push_back(ctype);
                continue;
            }
        }

        if (touching_ze >= 0
            && tdim == 3 && ctype == cutcells::cell::type::tetrahedron && geom_order >= 2)
        {
            int apex_vid = -1;
            for (int j = 0; j < n_vertices; ++j)
            {
                const int lv = local_mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];
                if (std::find(touching_face_verts.begin(), touching_face_verts.end(), lv)
                    == touching_face_verts.end())
                {
                    apex_vid = lv;
                    break;
                }
            }
            if (apex_vid >= 0)
            {
                std::vector<cutcells::BernsteinCell<T>> face_coord_polys;
                build_zero_face_coordinate_polys_local(
                    local_mesh, touching_ze, local_mesh.curved_geometry_order, face_coord_polys);
                const auto& tpl = cutcells::iso_p1_template(cutcells::cell::type::tetrahedron, geom_order);
                const T* r_apex = &local_mesh.vertex_ref_x[
                    static_cast<std::size_t>(apex_vid * tdim)];
                const int32_t v0 = static_cast<int32_t>(out.vertex_coords.size() / static_cast<std::size_t>(gdim));
                for (int vi = 0; vi < tpl.n_vertices; ++vi)
                {
                    const T xi = static_cast<T>(
                        tpl.ref_vertex_coords[static_cast<std::size_t>(vi * tpl.tdim + 0)]);
                    const T eta = static_cast<T>(
                        tpl.ref_vertex_coords[static_cast<std::size_t>(vi * tpl.tdim + 1)]);
                    const T zeta = static_cast<T>(
                        tpl.ref_vertex_coords[static_cast<std::size_t>(vi * tpl.tdim + 2)]);
                    T F_ref[3] = {};
                    T dF_dxi[3] = {};
                    T dF_deta[3] = {};
                    T dF_dzeta[3] = {};
                    eval_gh_tet_one_curved_tri_face<T>(
                        xi, eta, zeta,
                        std::span<const cutcells::BernsteinCell<T>>(
                            face_coord_polys.data(), face_coord_polys.size()),
                        r_apex,
                        F_ref, dF_dxi, dF_deta, dF_dzeta);
                    T x_phys[3] = {};
                    ref_to_phys_single_local(local_mesh, F_ref, x_phys);
                    for (int d = 0; d < gdim; ++d)
                        out.vertex_coords.push_back(x_phys[d]);
                    out.connectivity.push_back(v0 + static_cast<int32_t>(vi));
                }
                out.offsets.push_back(static_cast<int32_t>(out.connectivity.size()));
                out.cell_types.push_back(ctype);
                continue;
            }
        }

        const int32_t v0 = static_cast<int32_t>(out.vertex_coords.size() / static_cast<std::size_t>(gdim));

        for (int i = c0; i < c1; ++i)
        {
            const int lv = local_mesh.cell_vertices[static_cast<std::size_t>(i)];
            if (lv < 0 || lv >= local_mesh.n_vertices())
                throw std::runtime_error("append_local_cells_by_domain: local cell references invalid vertex id");
            const std::size_t base = static_cast<std::size_t>(lv * gdim);
            out.vertex_coords.insert(
                out.vertex_coords.end(),
                local_mesh.vertex_x.begin() + static_cast<std::ptrdiff_t>(base),
                local_mesh.vertex_x.begin() + static_cast<std::ptrdiff_t>(base + gdim));
        }

        for (int i = 0; i < n_vertices; ++i)
            out.connectivity.push_back(v0 + static_cast<int32_t>(i));
        out.offsets.push_back(static_cast<int32_t>(out.connectivity.size()));
        out.cell_types.push_back(ctype);
    }
}

template <std::floating_point T>
cutcells::mesh::CutMesh<T> curved_to_cut_mesh_direct(
    const cutcells::mesh::CurvedGlobalMesh<T>& curved)
{
    cutcells::mesh::CutMesh<T> out;
    out._gdim = curved.gdim;
    out._tdim = curved.cell_types.empty() ? 0 : cutcells::cell::get_tdim(curved.cell_types.front());
    out._num_cells = curved.n_cells();
    out._num_vertices = curved.n_vertices();
    out._vertex_coords = curved.vertex_coords;

    out._connectivity.resize(curved.connectivity.size());
    std::transform(curved.connectivity.begin(), curved.connectivity.end(),
                   out._connectivity.begin(),
                   [](int32_t v) { return static_cast<int>(v); });

    out._offset.resize(curved.offsets.size());
    std::transform(curved.offsets.begin(), curved.offsets.end(),
                   out._offset.begin(),
                   [](int32_t v) { return static_cast<int>(v); });

    out._types = curved.cell_types;
    out._parent_map.assign(static_cast<std::size_t>(out._num_cells), -1);
    return out;
}

template <std::floating_point T>
std::vector<T> gauss_lobatto_nodes_01(int order)
{
    std::vector<T> nodes;
    if (order <= 1)
    {
        nodes = {T(0), T(1)};
        return nodes;
    }

    nodes.reserve(static_cast<std::size_t>(order + 1));
    nodes.push_back(T(0));
    const auto interior = cutcells::detail::gauss_lobatto_interior_points_1d<T>(order);
    nodes.insert(nodes.end(), interior.begin(), interior.end());
    nodes.push_back(T(1));
    return nodes;
}

template <std::floating_point T>
void lagrange_eval_nd(std::span<const T> t_nodes,
                      std::span<const T> values,
                      int                gdim,
                      T                  t,
                      T*                 out)
{
    const int n = static_cast<int>(t_nodes.size());
    for (int d = 0; d < gdim; ++d)
        out[d] = T(0);

    for (int i = 0; i < n; ++i)
    {
        T li = T(1);
        const T ti = t_nodes[static_cast<std::size_t>(i)];
        for (int j = 0; j < n; ++j)
        {
            if (i == j)
                continue;
            const T tj = t_nodes[static_cast<std::size_t>(j)];
            li *= (t - tj) / (ti - tj);
        }

        for (int d = 0; d < gdim; ++d)
            out[d] += li * values[static_cast<std::size_t>(i * gdim + d)];
    }
}

template <std::floating_point T>
void append_line_segments_from_polyline(std::span<const T> polyline,
                                        int                gdim,
                                        cutcells::mesh::CutMesh<T>& out)
{
    const int n = static_cast<int>(polyline.size()) / gdim;
    if (n < 2)
        return;

    if (out._offset.empty())
        out._offset.push_back(0);

    for (int i = 0; i < n - 1; ++i)
    {
        const int v0 = static_cast<int>(out._vertex_coords.size() / static_cast<std::size_t>(gdim));
        const std::size_t base0 = static_cast<std::size_t>(i * gdim);
        const std::size_t base1 = static_cast<std::size_t>((i + 1) * gdim);

        out._vertex_coords.insert(
            out._vertex_coords.end(),
            polyline.begin() + static_cast<std::ptrdiff_t>(base0),
            polyline.begin() + static_cast<std::ptrdiff_t>(base0 + gdim));
        out._vertex_coords.insert(
            out._vertex_coords.end(),
            polyline.begin() + static_cast<std::ptrdiff_t>(base1),
            polyline.begin() + static_cast<std::ptrdiff_t>(base1 + gdim));

        out._connectivity.push_back(v0);
        out._connectivity.push_back(v0 + 1);
        out._offset.push_back(static_cast<int>(out._connectivity.size()));
        out._types.push_back(cutcells::cell::type::interval);
        out._parent_map.push_back(-1);
    }
}

template <std::floating_point T>
cutcells::mesh::CutMesh<T> linearize_curved_interface_mesh(
    const cutcells::mesh::CurvedGlobalMesh<T>& curved,
    int                                        subdivision)
{
    cutcells::mesh::CutMesh<T> out;
    out._gdim = curved.gdim;
    out._tdim = curved.cell_types.empty() ? 0 : cutcells::cell::get_tdim(curved.cell_types.front());
    out._offset.push_back(0);

    const int gdim = curved.gdim;

    for (int c = 0; c < curved.n_cells(); ++c)
    {
        const int32_t c0 = curved.offsets[static_cast<std::size_t>(c)];
        const int32_t c1 = curved.offsets[static_cast<std::size_t>(c + 1)];
        const int n_nodes = static_cast<int>(c1 - c0);
        const auto ct = curved.cell_types[static_cast<std::size_t>(c)];

        if (ct == cutcells::cell::type::interval)
        {
            if (n_nodes < 2)
                continue;

            std::vector<T> cell_points(static_cast<std::size_t>(n_nodes * gdim), T(0));
            for (int i = 0; i < n_nodes; ++i)
            {
                const int32_t gv = curved.connectivity[static_cast<std::size_t>(c0 + i)];
                for (int d = 0; d < gdim; ++d)
                {
                    cell_points[static_cast<std::size_t>(i * gdim + d)] =
                        curved.vertex_coords[static_cast<std::size_t>(gv * gdim + d)];
                }
            }

            if (n_nodes == 2)
            {
                append_line_segments_from_polyline<T>(
                    std::span<const T>(cell_points.data(), cell_points.size()), gdim, out);
                continue;
            }

            const int order = n_nodes - 1;
            const int n_segments = std::max(1, subdivision * order);
            std::vector<T> sampled(static_cast<std::size_t>((n_segments + 1) * gdim), T(0));
            auto t_nodes = gauss_lobatto_nodes_01<T>(order);
            if (static_cast<int>(t_nodes.size()) != n_nodes)
            {
                t_nodes.resize(static_cast<std::size_t>(n_nodes));
                for (int i = 0; i < n_nodes; ++i)
                    t_nodes[static_cast<std::size_t>(i)] = static_cast<T>(i) / static_cast<T>(n_nodes - 1);
            }

            for (int i = 0; i <= n_segments; ++i)
            {
                const T t = static_cast<T>(i) / static_cast<T>(n_segments);
                lagrange_eval_nd<T>(
                    std::span<const T>(t_nodes.data(), t_nodes.size()),
                    std::span<const T>(cell_points.data(), cell_points.size()),
                    gdim,
                    t,
                    sampled.data() + static_cast<std::size_t>(i * gdim));
            }

            append_line_segments_from_polyline<T>(
                std::span<const T>(sampled.data(), sampled.size()), gdim, out);
            continue;
        }

        const int v0 = static_cast<int>(out._vertex_coords.size() / static_cast<std::size_t>(gdim));
        for (int i = 0; i < n_nodes; ++i)
        {
            const int32_t gv = curved.connectivity[static_cast<std::size_t>(c0 + i)];
            const std::size_t base = static_cast<std::size_t>(gv * gdim);
            out._vertex_coords.insert(
                out._vertex_coords.end(),
                curved.vertex_coords.begin() + static_cast<std::ptrdiff_t>(base),
                curved.vertex_coords.begin() + static_cast<std::ptrdiff_t>(base + gdim));
            out._connectivity.push_back(v0 + i);
        }
        out._offset.push_back(static_cast<int>(out._connectivity.size()));
        out._types.push_back(ct);
        out._parent_map.push_back(-1);
    }

    out._num_cells = static_cast<int>(out._types.size());
    out._num_vertices = static_cast<int>(out._vertex_coords.size() / static_cast<std::size_t>(out._gdim));
    if (out._types.empty())
        out._tdim = 0;
    return out;
}

template <std::floating_point T>
cutcells::mesh::CutMesh<T> linearize_local_interfaces(
    std::span<const cutcells::LocalMesh<T>*> local_meshes,
    int                                      level_set_id,
    int                                      subdivision)
{
    cutcells::mesh::CutMesh<T> out;
    out._offset.push_back(0);
    (void)subdivision;

    if (local_meshes.empty())
    {
        out._gdim = 0;
        out._tdim = 0;
        out._num_cells = 0;
        out._num_vertices = 0;
        return out;
    }

    out._gdim = local_meshes.front()->gdim;
    out._tdim = std::max(0, local_meshes.front()->tdim - 1);

    const uint64_t zero_mask = uint64_t(1) << level_set_id;

    for (const cutcells::LocalMesh<T>* mesh_ptr : local_meshes)
    {
        if (mesh_ptr == nullptr)
            continue;

        const auto& mesh = *mesh_ptr;
        if (mesh.n_zero_entities() == 0)
            continue;

        if (mesh.tdim == 2)
        {
            if (mesh.n_zero_chains() == 0)
            {
                auto& nonconst_mesh = const_cast<cutcells::LocalMesh<T>&>(mesh);
                build_zero_chains(nonconst_mesh, zero_mask);
            }

            for (int chain = 0; chain < mesh.n_zero_chains(); ++chain)
            {
                if (mesh.zero_chain_zero_mask[static_cast<std::size_t>(chain)] != zero_mask)
                    continue;

                const int c0 = mesh.zero_chain_offsets[static_cast<std::size_t>(chain)];
                const int c1 = mesh.zero_chain_offsets[static_cast<std::size_t>(chain + 1)];

                for (int ci = c0; ci < c1; ++ci)
                {
                    const int ze = mesh.zero_chain_entity_ids[static_cast<std::size_t>(ci)];
                    const bool reversed = mesh.zero_chain_entity_reversed[static_cast<std::size_t>(ci)] != 0;

                    const int32_t lv0 = reversed
                        ? mesh.zero_entity_endpoint_v1[static_cast<std::size_t>(ze)]
                        : mesh.zero_entity_endpoint_v0[static_cast<std::size_t>(ze)];
                    const int32_t lv1 = reversed
                        ? mesh.zero_entity_endpoint_v0[static_cast<std::size_t>(ze)]
                        : mesh.zero_entity_endpoint_v1[static_cast<std::size_t>(ze)];
                    if (lv0 < 0 || lv1 < 0)
                        continue;

                    std::vector<T> polyline(static_cast<std::size_t>(2 * mesh.gdim), T(0));
                    for (int d = 0; d < mesh.gdim; ++d)
                    {
                        polyline[static_cast<std::size_t>(d)] =
                            mesh.vertex_x[static_cast<std::size_t>(lv0 * mesh.gdim + d)];
                        polyline[static_cast<std::size_t>(mesh.gdim + d)] =
                            mesh.vertex_x[static_cast<std::size_t>(lv1 * mesh.gdim + d)];
                    }

                    append_line_segments_from_polyline<T>(
                        std::span<const T>(polyline.data(), polyline.size()), mesh.gdim, out);
                }
            }
        }
        else if (mesh.tdim == 3)
        {
            const int n_zero = mesh.n_zero_entities();
            for (int ze = 0; ze < n_zero; ++ze)
            {
                if (mesh.zero_entity_dim[static_cast<std::size_t>(ze)] != 2)
                    continue;
                if ((mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] & zero_mask) == 0)
                    continue;

                const int z0 = mesh.zero_entity_offsets[static_cast<std::size_t>(ze)];
                const int z1 = mesh.zero_entity_offsets[static_cast<std::size_t>(ze + 1)];
                const int n_face_verts = z1 - z0;
                if (n_face_verts != 3 && n_face_verts != 4)
                    continue;

                const int32_t base
                    = static_cast<int32_t>(out._vertex_coords.size()
                                           / static_cast<std::size_t>(mesh.gdim));
                for (int i = z0; i < z1; ++i)
                {
                    const int lv = mesh.zero_entity_vertices[static_cast<std::size_t>(i)];
                    if (lv < 0 || lv >= mesh.n_vertices())
                        throw std::runtime_error(
                            "linearize_local_interfaces: zero face references invalid vertex id");
                    const std::size_t xb = static_cast<std::size_t>(lv * mesh.gdim);
                    for (int d = 0; d < mesh.gdim; ++d)
                        out._vertex_coords.push_back(
                            mesh.vertex_x[xb + static_cast<std::size_t>(d)]);
                    out._connectivity.push_back(base + static_cast<int32_t>(i - z0));
                }
                out._offset.push_back(static_cast<int>(out._connectivity.size()));
                out._types.push_back(
                    n_face_verts == 3 ? cutcells::cell::type::triangle
                                      : cutcells::cell::type::quadrilateral);
                out._parent_map.push_back(mesh.parent_cell_id);
            }
        }
    }

    out._num_cells = static_cast<int>(out._types.size());
    out._num_vertices = static_cast<int>(out._vertex_coords.size() / static_cast<std::size_t>(out._gdim));
    return out;
}

template <std::floating_point T>
cutcells::mesh::CutMesh<T> sample_curved_local_interfaces(
    std::span<const cutcells::LocalMesh<T>*> local_meshes,
    int                                      level_set_id,
    int                                      geom_order,
    int                                      subdivision)
{
    cutcells::mesh::CutMesh<T> out;
    out._offset.push_back(0);

    if (local_meshes.empty())
    {
        out._gdim = 0;
        out._tdim = 0;
        out._num_cells = 0;
        out._num_vertices = 0;
        return out;
    }

    out._gdim = local_meshes.front()->gdim;
    out._tdim = std::max(0, local_meshes.front()->tdim - 1);

    const uint64_t zero_mask = uint64_t(1) << level_set_id;
    const int n_sub = std::max(1, subdivision);
    constexpr T surface_merge_tol = static_cast<T>(1e-10);
    std::unordered_map<SurfacePointKey, int, SurfacePointKeyHash> point_map;

    for (const cutcells::LocalMesh<T>* mesh_ptr : local_meshes)
    {
        if (mesh_ptr == nullptr)
            continue;

        const auto& mesh = *mesh_ptr;
        if (mesh.n_zero_entities() == 0)
            continue;

        if (mesh.tdim == 2)
        {
            // 2D curved interface line sampling is not changed here.
            continue;
        }
        if (mesh.tdim != 3)
            continue;

        const int n_zero = mesh.n_zero_entities();
        for (int ze = 0; ze < n_zero; ++ze)
        {
            if (mesh.zero_entity_dim[static_cast<std::size_t>(ze)] != 2)
                continue;
            if ((mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] & zero_mask) == 0)
                continue;

            const auto face_type = zero_face_type_local(mesh, ze);
            if (face_type != cutcells::cell::type::triangle
                && face_type != cutcells::cell::type::quadrilateral)
                continue;

            std::vector<cutcells::BernsteinCell<T>> face_coord_polys;
            build_zero_face_coordinate_polys_local(
                mesh, ze, geom_order, face_coord_polys);

            const auto& tpl = (n_sub <= 1)
                ? cutcells::p1_template(face_type)
                : cutcells::iso_p1_template(face_type, std::min(n_sub, 4));

            std::vector<int32_t> sample_node_ids(
                static_cast<std::size_t>(tpl.n_vertices), int32_t(-1));
            std::array<T, 3> x_ref = {T(0), T(0), T(0)};
            std::array<T, 3> x_phys = {T(0), T(0), T(0)};

            for (int vi = 0; vi < tpl.n_vertices; ++vi)
            {
                const T u = static_cast<T>(
                    tpl.ref_vertex_coords[static_cast<std::size_t>(vi * tpl.tdim + 0)]);
                const T v = static_cast<T>(
                    tpl.ref_vertex_coords[static_cast<std::size_t>(vi * tpl.tdim + 1)]);
                const std::array<T, 2> uv = {u, v};
                for (int d = 0; d < mesh.tdim; ++d)
                {
                    x_ref[static_cast<std::size_t>(d)] =
                        cutcells::evaluate_bernstein_cell<T>(
                            face_coord_polys[static_cast<std::size_t>(d)],
                            std::span<const T>(uv.data(), uv.size()));
                }

                ref_to_phys_single_local(mesh, x_ref.data(), x_phys.data());
                sample_node_ids[static_cast<std::size_t>(vi)] = add_or_get_surface_point(
                    out, point_map, x_phys.data(), surface_merge_tol);
            }

            for (int child = 0; child < tpl.n_cells; ++child)
            {
                const int base = child * tpl.vertices_per_cell;
                for (int j = 0; j < tpl.vertices_per_cell; ++j)
                {
                    const int lv = tpl.cell_connectivity[
                        static_cast<std::size_t>(base + j)];
                    out._connectivity.push_back(
                        sample_node_ids[static_cast<std::size_t>(lv)]);
                }
                out._offset.push_back(static_cast<int>(out._connectivity.size()));
                out._types.push_back(tpl.child_cell_type);
                out._parent_map.push_back(-1);
            }
        }
    }

    out._num_cells = static_cast<int>(out._types.size());
    out._num_vertices = static_cast<int>(
        out._vertex_coords.size() / static_cast<std::size_t>(out._gdim));
    return out;
}

} // namespace

namespace cutcells::mesh
{

template <std::floating_point T>
CurvedCutMeshResult<T> cut_mesh_view_curved(
    const MeshView<T, int>&         mesh,
    const LevelSetFunction<T, int>& level_set,
    int                             geom_order,
    std::string_view                backend,
    int                             vis_subdivision,
    T                               tol,
    bool                            repair_diagonals,
    int                             max_repair_depth)
{
    if (!mesh.has_cell_types())
        throw std::invalid_argument("cut_mesh_view_curved: mesh cell_types are required");
    if (!level_set.has_value() && !level_set.has_nodal_values())
    {
        throw std::invalid_argument(
            "cut_mesh_view_curved: level_set requires either callable value(x) or nodal_values");
    }
    if (geom_order < 1)
        throw std::invalid_argument("cut_mesh_view_curved: geom_order must be >= 1");
    if (vis_subdivision < 1)
        throw std::invalid_argument("cut_mesh_view_curved: vis_subdivision must be >= 1");

    const std::string backend_name = to_lower_ascii(backend);
    if (backend_name != "bernstein")
        throw std::invalid_argument("cut_mesh_view_curved: only backend='bernstein' is supported");

    if (mesh.gdim <= 0 || mesh.gdim > 3)
        throw std::invalid_argument("cut_mesh_view_curved: mesh.gdim must be in [1, 3]");

    CurvedCutMeshResult<T> result;
    init_curved_mesh(result.inside, mesh.gdim, geom_order, 0);
    init_curved_mesh(result.interface, mesh.gdim, geom_order, 0);
    init_curved_mesh(result.outside, mesh.gdim, geom_order, 0);
    const BackgroundEntityMaps bg_maps = build_background_entity_maps(mesh);
    const bool use_fem_nodal_level_set
        = level_set.has_nodal_values()
          && level_set.kind == LevelSetFunction<T, int>::Kind::fem_nodal;
    if (use_fem_nodal_level_set)
    {
        if (level_set.degree < 1)
            throw std::invalid_argument(
                "cut_mesh_view_curved: fem_nodal level set requires degree >= 1");
        if (level_set.nodal_values.size()
            != static_cast<std::size_t>(mesh.num_nodes()))
        {
            throw std::invalid_argument(
                "cut_mesh_view_curved: fem_nodal nodal_values size must match mesh num_nodes");
        }
    }

    const int n_cells = mesh.num_cells();
    result.local_meshes.reserve(static_cast<std::size_t>(n_cells));

    const bool use_repair_callbacks
        = repair_diagonals && level_set.has_value() && !use_fem_nodal_level_set;
    const bool use_repair_bernstein = repair_diagonals && !use_repair_callbacks;
    if (repair_diagonals && !use_repair_callbacks && !use_fem_nodal_level_set)
    {
        std::cerr
            << "Warning: cut_mesh_view_curved: repair_diagonals requested but the "
            << "level set has no callable value(x) or FEM nodal values; proceeding without interface-split repair\n";
    }

    for (int cell_id = 0; cell_id < n_cells; ++cell_id)
    {
        const auto ct = cell::map_vtk_type_to_cell_type(
            static_cast<cell::vtk_types>(mesh.cell_type(cell_id)));
        if (!bernstein_supported_cell(ct))
        {
            throw std::invalid_argument(
                "cut_mesh_view_curved: bernstein backend currently supports "
                "interval/triangle/tetrahedron/hexahedron");
        }

        const int p1_vertices = cell::get_num_vertices(ct);
        if (mesh.cell_num_nodes(cell_id) < p1_vertices)
            throw std::runtime_error("cut_mesh_view_curved: parent cell has too few nodes");

        std::vector<T> parent_coords_p1(static_cast<std::size_t>(p1_vertices * mesh.gdim), T(0));
        for (int i = 0; i < p1_vertices; ++i)
        {
            const int node_id = mesh.cell_node(cell_id, i);
            if (node_id < 0 || node_id >= mesh.num_nodes())
                throw std::runtime_error("cut_mesh_view_curved: invalid node id in connectivity");

            const T* x = mesh.node(node_id);
            for (int d = 0; d < mesh.gdim; ++d)
                parent_coords_p1[static_cast<std::size_t>(i * mesh.gdim + d)] = x[d];
        }

        std::vector<T> sampled_node_coords;
        std::vector<T> cell_nodal_values;
        int cell_level_set_degree = std::max(geom_order, 2);
        if (use_fem_nodal_level_set)
        {
            cell_level_set_degree = level_set.degree;
            cell_nodal_values = cutcells::extract_cell_nodal_values(
                mesh, level_set.nodal_values, cell_id);
            const int expected_nodes
                = cutcells::lagrange_node_count(ct, cell_level_set_degree);
            if (static_cast<int>(cell_nodal_values.size()) != expected_nodes)
            {
                throw std::invalid_argument(
                    "cut_mesh_view_curved: fem_nodal cell node count does not match level_set degree");
            }
        }
        else
        {
            sample_callable_on_parent_nodes<T>(
                ct,
                std::span<const T>(parent_coords_p1.data(), parent_coords_p1.size()),
                mesh.gdim,
                cell_level_set_degree,
                level_set,
                cell_id,
                sampled_node_coords,
                cell_nodal_values);
        }

        const auto domain = classify_parent_cell_bernstein<T>(
            ct,
            cell_level_set_degree,
            std::span<const T>(cell_nodal_values.data(), cell_nodal_values.size()),
            tol);

        if (domain == ParentCellDomain::inside)
        {
            append_cell_to_curved_mesh<T>(
                result.inside,
                ct,
                std::span<const T>(parent_coords_p1.data(), parent_coords_p1.size()),
                mesh.gdim);
            result.n_parent_inside += 1;
            continue;
        }

        if (domain == ParentCellDomain::outside)
        {
            append_cell_to_curved_mesh<T>(
                result.outside,
                ct,
                std::span<const T>(parent_coords_p1.data(), parent_coords_p1.size()),
                mesh.gdim);
            result.n_parent_outside += 1;
            continue;
        }

        result.n_parent_intersected += 1;

        LevelSetFunction<T, int> local_level_set;
        local_level_set.nodal_values
            = std::span<const T>(cell_nodal_values.data(), cell_nodal_values.size());
        local_level_set.gdim = mesh.gdim;
        local_level_set.degree = cell_level_set_degree;
        local_level_set.kind = LevelSetFunction<T, int>::Kind::fem_nodal;
        local_level_set.owner = level_set.owner;

        LocalMesh<T> local_mesh;
        init_local_mesh_from_cell<T>(
            local_mesh,
            std::span<const T>(parent_coords_p1.data(), parent_coords_p1.size()),
            ct,
            cell_id,
            1);

        if (use_repair_callbacks)
        {
            decompose_local_mesh_with_backend<T, int>(
                local_mesh,
                level_set,
                LocalLevelSetBackend::analytical_callbacks,
                0,
                cell::edge_root::method::itp,
                true,
                6,
                tol,
                false,
                true,
                max_repair_depth);
        }
        else
        {
            decompose_local_mesh_with_backend<T, int>(
                local_mesh,
                local_level_set,
                LocalLevelSetBackend::bernstein,
                0,
                cell::edge_root::method::itp,
                true,
                6,
                tol,
                false,
                use_repair_bernstein,
                max_repair_depth);
        }

        build_zero_entities<T>(local_mesh, 0);
        if (level_set.has_value() && level_set.has_gradient())
        {
            curve_zero_entities_with_backend<T, int>(
                local_mesh,
                level_set,
                LocalLevelSetBackend::analytical_callbacks,
                0,
                geom_order,
                tol);
        }
        else
        {
            curve_zero_entities_with_backend<T, int>(
                local_mesh,
                local_level_set,
                LocalLevelSetBackend::bernstein,
                0,
                geom_order,
                tol);
        }

        {
            const uint64_t zero_mask = uint64_t(1) << 0;
            if (local_mesh.tdim == 2)
                build_zero_chains(local_mesh, zero_mask);
            else if (local_mesh.tdim == 3)
                build_zero_patches(local_mesh, zero_mask);
        }

        append_local_cells_by_domain<T>(
            local_mesh, cell::domain::inside, result.inside, geom_order);
        append_local_cells_by_domain<T>(
            local_mesh, cell::domain::outside, result.outside, geom_order);

        result.local_meshes.push_back(std::move(local_mesh));
    }

    if (!result.local_meshes.empty())
    {
        std::vector<const LocalMesh<T>*> local_ptrs;
        local_ptrs.reserve(result.local_meshes.size());
        for (const LocalMesh<T>& local_mesh : result.local_meshes)
            local_ptrs.push_back(&local_mesh);
        const bool interface_is_3d = !result.local_meshes.empty()
            && result.local_meshes.front().tdim == 3;

        if (interface_is_3d)
        {
            result.interface = assemble_curved_interface_mesh_parent_mapped<T>(
                std::span<const LocalMesh<T>*>(local_ptrs.data(), local_ptrs.size()),
                bg_maps,
                0,
                geom_order);
        }
        else
        {
            result.interface = assemble_curved_interface_mesh<T>(
                std::span<const LocalMesh<T>*>(local_ptrs.data(), local_ptrs.size()),
                0,
                geom_order);
        }
        result.interface_vis = linearize_local_interfaces<T>(
            std::span<const LocalMesh<T>*>(local_ptrs.data(), local_ptrs.size()),
            0,
            vis_subdivision);
        if (interface_is_3d)
        {
            result.interface_curved_vis = sample_curved_local_interfaces_parent_mapped<T>(
                std::span<const LocalMesh<T>*>(local_ptrs.data(), local_ptrs.size()),
                bg_maps,
                0,
                geom_order,
                vis_subdivision);
        }
        else
        {
            result.interface_curved_vis = sample_curved_local_interfaces<T>(
                std::span<const LocalMesh<T>*>(local_ptrs.data(), local_ptrs.size()),
                0,
                geom_order,
                vis_subdivision);
        }
    }

    result.inside_vis = curved_to_cut_mesh_direct<T>(result.inside);
    result.outside_vis = curved_to_cut_mesh_direct<T>(result.outside);
    if (result.interface_vis._types.empty())
        result.interface_vis = linearize_curved_interface_mesh<T>(result.interface, vis_subdivision);
    if (result.interface_curved_vis._types.empty())
        result.interface_curved_vis = result.interface_vis;

    return result;
}

template <std::floating_point T>
CurvedCutMeshResult<T> cut_vtk_mesh_curved(
    std::span<const T>              points,
    std::span<const int>            connectivity,
    std::span<const int>            offset,
    std::span<const int>            vtk_type,
    const LevelSetFunction<T, int>& level_set,
    int                             geom_order,
    std::string_view                backend,
    int                             vis_subdivision,
    T                               tol,
    bool                            repair_diagonals,
    int                             max_repair_depth)
{
    if (offset.empty())
        throw std::invalid_argument("cut_vtk_mesh_curved: offset must not be empty");
    if (vtk_type.size() + 1 != offset.size())
        throw std::invalid_argument("cut_vtk_mesh_curved: offset size must be vtk_type size + 1");

    int max_node = -1;
    for (const int node : connectivity)
    {
        if (node < 0)
            throw std::invalid_argument("cut_vtk_mesh_curved: connectivity contains negative node index");
        max_node = std::max(max_node, node);
    }

    const int n_nodes = max_node + 1;
    if (n_nodes <= 0)
        throw std::invalid_argument("cut_vtk_mesh_curved: could not infer node count from connectivity");
    if (points.size() % static_cast<std::size_t>(n_nodes) != 0)
        throw std::invalid_argument("cut_vtk_mesh_curved: points size is not divisible by inferred node count");

    MeshView<T, int> mesh;
    mesh.gdim = static_cast<int>(points.size() / static_cast<std::size_t>(n_nodes));
    mesh.tdim = 0;
    mesh.coordinates = points;
    mesh.connectivity = connectivity;
    mesh.offsets = offset;
    mesh.cell_types = vtk_type;

    return cut_mesh_view_curved<T>(
        mesh,
        level_set,
        geom_order,
        backend,
        vis_subdivision,
        tol,
        repair_diagonals,
        max_repair_depth);
}

// explicit instantiations

template CurvedCutMeshResult<float> cut_mesh_view_curved<float>(
    const MeshView<float, int>&,
    const LevelSetFunction<float, int>&,
    int,
    std::string_view,
    int,
    float,
    bool,
    int);

template CurvedCutMeshResult<double> cut_mesh_view_curved<double>(
    const MeshView<double, int>&,
    const LevelSetFunction<double, int>&,
    int,
    std::string_view,
    int,
    double,
    bool,
    int);

template CurvedCutMeshResult<float> cut_vtk_mesh_curved<float>(
    std::span<const float>,
    std::span<const int>,
    std::span<const int>,
    std::span<const int>,
    const LevelSetFunction<float, int>&,
    int,
    std::string_view,
    int,
    float,
    bool,
    int);

template CurvedCutMeshResult<double> cut_vtk_mesh_curved<double>(
    std::span<const double>,
    std::span<const int>,
    std::span<const int>,
    std::span<const int>,
    const LevelSetFunction<double, int>&,
    int,
    std::string_view,
    int,
    double,
    bool,
    int);

} // namespace cutcells::mesh
