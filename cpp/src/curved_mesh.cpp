// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "curved_mesh.h"
#include "bernstein_backend.h"
#include "cell_flags.h"
#include "edge_root.h"
#include "local_level_set.h"
#include "mapping.h"
#include "mapping_curved.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
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
int32_t add_or_get_surface_point(
    cutcells::mesh::CurvedGlobalMesh<T>& mesh,
    std::unordered_map<SurfacePointKey, int32_t, SurfacePointKeyHash>& point_map,
    int32_t& gv_count,
    const T* x_phys,
    int gdim,
    T tol)
{
    const auto key = make_surface_point_key(x_phys, gdim, tol);
    const auto it = point_map.find(key);
    if (it != point_map.end())
        return it->second;

    const int32_t gv = gv_count++;
    for (int d = 0; d < gdim; ++d)
        mesh.vertex_coords.push_back(x_phys[d]);
    point_map.emplace(key, gv);
    return gv;
}

/// Push-forward a single reference point to physical space using parent P1
/// affine map for the given LocalMesh.
template <std::floating_point T>
void ref_to_phys_single(const cutcells::LocalMesh<T>& mesh,
                        const T* x_ref,
                        T* x_phys)
{
    const int gdim = mesh.gdim;
    const int tdim = mesh.tdim;
    if (gdim <= 0 || gdim > 3)
        throw std::invalid_argument("ref_to_phys_single: gdim must be in [1, 3]");
    if (tdim < 0 || tdim > 3)
        throw std::invalid_argument("ref_to_phys_single: tdim must be in [0, 3]");
    if (tdim > gdim)
        throw std::invalid_argument("ref_to_phys_single: tdim must be <= gdim");

    cutcells::cell::ref_to_phys_affine<T>(
        mesh.parent_cell_type,
        std::span<const T>(mesh.parent_cell_coords_p1.data(),
                           mesh.parent_cell_coords_p1.size()),
        gdim,
        tdim,
        x_ref,
        x_phys);
}

template <std::floating_point T>
inline cutcells::cell::type zero_face_type(const cutcells::LocalMesh<T>& mesh,
                                           int zero_entity_id)
{
    return cutcells::zero_entity_face_type(mesh, zero_entity_id);
}

template <std::floating_point T>
inline void build_zero_face_coordinate_polys(
    const cutcells::LocalMesh<T>&            mesh,
    int                                      zero_entity_id,
    int                                      geom_order,
    std::vector<cutcells::BernsteinCell<T>>& coord_polys,
    int*                                     fallback_nodes = nullptr)
{
    const auto face_type = zero_face_type(mesh, zero_entity_id);
    if (face_type != cutcells::cell::type::triangle
        && face_type != cutcells::cell::type::quadrilateral)
        throw std::invalid_argument(
            "build_zero_face_coordinate_polys: unsupported zero-face type");

    const int tdim = mesh.tdim;
    std::vector<T> face_nodes_ref;
    const int degree = cutcells::build_zero_face_nodes_ref(
        mesh, zero_entity_id, geom_order, face_nodes_ref, fallback_nodes);
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
} // anonymous namespace

namespace cutcells::mesh
{

// ============================================================================
// assemble_curved_interface_mesh — cache-reading implementation
// ============================================================================

template <std::floating_point T>
CurvedGlobalMesh<T> assemble_curved_interface_mesh(
    std::span<const LocalMesh<T>*> local_meshes,
    int                            level_set_id,
    int                            geom_order)
{
    CurvedGlobalMesh<T> result;
    result.geom_order = geom_order;
    result.level_set_id = level_set_id;

    if (local_meshes.empty())
    {
        result.offsets.push_back(0);
        return result;
    }

    result.gdim = local_meshes[0]->gdim;
    const int gdim = result.gdim;

    int32_t gv_count = 0;
    constexpr T surface_merge_tol = static_cast<T>(1e-10);
    std::unordered_map<SurfacePointKey, int32_t, SurfacePointKeyHash> point_map;
    result.offsets.push_back(0);

    for (std::size_t mi = 0; mi < local_meshes.size(); ++mi)
    {
        const auto& mesh = *local_meshes[mi];

        if (mesh.gdim != gdim)
            throw std::invalid_argument(
                "assemble_curved_interface_mesh: inconsistent gdim");

        const int n_zero = mesh.n_zero_entities();
        if (n_zero == 0)
            continue;

        const int tdim = mesh.tdim;
        const bool has_curved = !mesh.curved_zero_offsets.empty()
            && static_cast<int>(mesh.curved_zero_offsets.size()) == n_zero + 1;

        if (tdim == 2)
        {
            const uint64_t zero_mask = uint64_t(1) << level_set_id;
            if (mesh.n_zero_chains() == 0)
                throw std::invalid_argument(
                    "assemble_curved_interface_mesh: zero chains not built; "
                    "call build_zero_chains before cache-reading assembly");

            std::unordered_map<int32_t, int32_t> endpoint_vertex_map;
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

                    auto get_or_add_endpoint = [&](int32_t lv) -> int32_t
                    {
                        const auto it = endpoint_vertex_map.find(lv);
                        if (it != endpoint_vertex_map.end())
                            return it->second;
                        const int32_t gv = gv_count++;
                        for (int d = 0; d < gdim; ++d)
                            result.vertex_coords.push_back(
                                mesh.vertex_x[static_cast<std::size_t>(lv * gdim + d)]);
                        endpoint_vertex_map.emplace(lv, gv);
                        return gv;
                    };

                    result.connectivity.push_back(get_or_add_endpoint(lv0));
                    result.connectivity.push_back(get_or_add_endpoint(lv1));
                    const int32_t node_start = has_curved ? mesh.curved_zero_offsets[static_cast<std::size_t>(ze)] : 0;
                    const int32_t node_end = has_curved ? mesh.curved_zero_offsets[static_cast<std::size_t>(ze + 1)] : 0;
                    auto append_curved_node = [&](int32_t ni)
                    {
                        const bool converged = !mesh.curved_zero_converged.empty()
                            && mesh.curved_zero_converged[static_cast<std::size_t>(ni)];
                        if (!converged)
                            ++result.n_fallback_nodes;
                        const T* x_ref = &mesh.curved_zero_ref_nodes[static_cast<std::size_t>(ni * tdim)];
                        std::vector<T> x_phys(static_cast<std::size_t>(gdim));
                        ref_to_phys_single(mesh, x_ref, x_phys.data());
                        for (int d = 0; d < gdim; ++d)
                            result.vertex_coords.push_back(x_phys[static_cast<std::size_t>(d)]);
                        result.connectivity.push_back(gv_count++);
                    };

                    if (reversed)
                    {
                        for (int32_t ni = node_end; ni > node_start; --ni)
                            append_curved_node(ni - 1);
                    }
                    else
                    {
                        for (int32_t ni = node_start; ni < node_end; ++ni)
                            append_curved_node(ni);
                    }
                    result.offsets.push_back(static_cast<int32_t>(result.connectivity.size()));
                    result.cell_types.push_back(cell::type::interval);
                }
            }
        }
        else if (tdim == 3)
        {
            const uint64_t zero_mask = uint64_t(1) << level_set_id;
            if (mesh.n_zero_patches() == 0)
                throw std::invalid_argument(
                    "assemble_curved_interface_mesh: zero patches not built; "
                    "call build_zero_patches before cache-reading assembly");

            for (int patch = 0; patch < mesh.n_zero_patches(); ++patch)
            {
                if (mesh.zero_patch_zero_mask[static_cast<std::size_t>(patch)] != zero_mask)
                    continue;
                const int p0 = mesh.zero_patch_offsets[static_cast<std::size_t>(patch)];
                const int p1 = mesh.zero_patch_offsets[static_cast<std::size_t>(patch + 1)];
                for (int pi = p0; pi < p1; ++pi)
                {
                    const int ze = mesh.zero_patch_entity_ids[static_cast<std::size_t>(pi)];
                    const auto face_type = zero_face_type(mesh, ze);
                    std::vector<T> face_nodes_ref;
                    const int degree = cutcells::build_zero_face_nodes_ref(
                        mesh, ze, geom_order, face_nodes_ref, &result.n_fallback_nodes);
                    const auto& tpl = (degree == 1)
                        ? cutcells::p1_template(face_type)
                        : cutcells::iso_p1_template(face_type, degree);
                    const auto& vtk_to_basix
                        = cutcells::vtk_lagrange_to_basix_permutation(face_type, degree);

                    for (int vtk_vi = 0; vtk_vi < tpl.n_vertices; ++vtk_vi)
                    {
                        const int vi = vtk_to_basix[static_cast<std::size_t>(vtk_vi)];
                        const T* x_ref = &face_nodes_ref[static_cast<std::size_t>(vi * tdim)];
                        std::vector<T> x_phys(static_cast<std::size_t>(gdim));
                        ref_to_phys_single(mesh, x_ref, x_phys.data());
                        result.connectivity.push_back(add_or_get_surface_point(
                            result, point_map, gv_count, x_phys.data(), gdim,
                            surface_merge_tol));
                    }

                    result.offsets.push_back(static_cast<int32_t>(result.connectivity.size()));
                    result.cell_types.push_back(face_type);
                }
            }
        }
    }

    return result;
}

// ============================================================================
// assemble_curved_interface_mesh — convenience wrapper (builds cache if needed)
// ============================================================================

template <std::floating_point T, std::integral I>
CurvedGlobalMesh<T> assemble_curved_interface_mesh(
    std::span<LocalMesh<T>*>          local_meshes,
    const LevelSetFunction<T, I>&     level_set,
    LocalLevelSetBackend              backend,
    int                               level_set_id,
    int                               geom_order,
    T                                 tol)
{
    // Ensure each mesh has interface entities and curved cache built
    for (auto* pm : local_meshes)
    {
        if (pm->n_zero_entities() == 0)
            build_zero_entities(*pm, level_set_id);
        if (pm->curved_cache_version == 0
            || pm->curved_cache_zero_version != pm->zero_entity_version
            || pm->curved_geometry_order != geom_order)
            curve_zero_entities_with_backend<T, I>(
                *pm, level_set, backend, level_set_id, geom_order, tol);
        const uint64_t zero_mask = uint64_t(1) << level_set_id;
        if (pm->tdim == 2 && pm->n_zero_chains() == 0)
            build_zero_chains(*pm, zero_mask);
        if (pm->tdim == 3 && pm->n_zero_patches() == 0)
            build_zero_patches(*pm, zero_mask);
    }

    // Delegate to the cache-reading assembler
    std::vector<const LocalMesh<T>*> const_ptrs(local_meshes.size());
    for (std::size_t i = 0; i < local_meshes.size(); ++i)
        const_ptrs[i] = local_meshes[i];
    return assemble_curved_interface_mesh<T>(
        std::span<const LocalMesh<T>*>(const_ptrs.data(), const_ptrs.size()),
        level_set_id, geom_order);
}

// ============================================================================
// assemble_curved_volume_mesh
// ============================================================================

template <std::floating_point T>
CurvedGlobalMesh<T> assemble_curved_volume_mesh(
    std::span<const LocalMesh<T>*> local_meshes,
    int level_set_id,
    int geom_order)
{
    CurvedGlobalMesh<T> result;
    result.geom_order = geom_order;
    result.level_set_id = level_set_id;

    if (local_meshes.empty())
    {
        result.offsets.push_back(0);
        return result;
    }

    result.gdim = local_meshes[0]->gdim;
    const int gdim = result.gdim;

    int32_t gv_count = 0;
    result.offsets.push_back(0);

    for (std::size_t mi = 0; mi < local_meshes.size(); ++mi)
    {
        const auto& mesh = *local_meshes[mi];

        if (mesh.gdim != gdim)
            throw std::invalid_argument(
                "assemble_curved_volume_mesh: inconsistent gdim");

        const int nc = mesh.n_cells();

        for (int c = 0; c < nc; ++c)
        {
            if (mesh.cell_domain[static_cast<std::size_t>(c)]
                != static_cast<uint8_t>(cell::domain::inside))
                continue;

            const int c0 = mesh.cell_offsets[static_cast<std::size_t>(c)];
            const int c1 = mesh.cell_offsets[static_cast<std::size_t>(c + 1)];

            for (int ci = c0; ci < c1; ++ci)
            {
                const int lv = mesh.cell_vertices[static_cast<std::size_t>(ci)];
                const int32_t gv = gv_count++;
                for (int d = 0; d < gdim; ++d)
                    result.vertex_coords.push_back(
                        mesh.vertex_x[static_cast<std::size_t>(lv * gdim + d)]);
                result.connectivity.push_back(gv);
            }

            result.offsets.push_back(
                static_cast<int32_t>(result.connectivity.size()));
            result.cell_types.push_back(
                mesh.cell_types[static_cast<std::size_t>(c)]);
        }
    }

    return result;
}

// ============================================================================
// assemble_curved_volume_mesh — with mapping backend selection
// ============================================================================

template <std::floating_point T>
CurvedGlobalMesh<T> assemble_curved_volume_mesh(
    std::span<const LocalMesh<T>*> local_meshes,
    int level_set_id,
    int geom_order,
    mapping::CurvedMappingBackend volume_backend,
    int vis_subdivision)
{
    using mapping::eval_pk_map;
    using mapping::eval_pk_deriv;
    using mapping::eval_gh_triangle_one_curved_edge;
    using mapping::eval_gh_quad_one_curved_edge;
    using mapping::eval_gh_tet_one_curved_tri_face;

    CurvedGlobalMesh<T> result;
    result.geom_order = geom_order;
    result.level_set_id = level_set_id;

    if (local_meshes.empty())
    {
        result.offsets.push_back(0);
        return result;
    }

    result.gdim = local_meshes[0]->gdim;
    const int gdim = result.gdim;

    int32_t gv_count = 0;
    result.offsets.push_back(0);

    const int n_sub = std::max(1, vis_subdivision);

    for (std::size_t mi = 0; mi < local_meshes.size(); ++mi)
    {
        const auto& mesh = *local_meshes[mi];

        if (mesh.gdim != gdim)
            throw std::invalid_argument(
                "assemble_curved_volume_mesh: inconsistent gdim");

        const int nc = mesh.n_cells();
        const int tdim = mesh.tdim;
        const bool has_curved = mesh.curved_geometry_order >= 2
            && !mesh.curved_zero_offsets.empty();
        const uint64_t zero_mask = uint64_t(1) << level_set_id;
        const int n_zero = mesh.n_zero_entities();
        const int codim1_dim = std::max(0, tdim - 1);

        for (int c = 0; c < nc; ++c)
        {
            if (mesh.cell_domain[static_cast<std::size_t>(c)]
                != static_cast<uint8_t>(cell::domain::inside))
                continue;

            const int c0 = mesh.cell_offsets[static_cast<std::size_t>(c)];
            const int c1 = mesh.cell_offsets[static_cast<std::size_t>(c + 1)];
            const int nv = c1 - c0;
            const cell::type ctype = mesh.cell_types[static_cast<std::size_t>(c)];

            // Check if this cell touches a zero entity
            int touching_ze = -1;
            int cell_if_a = -1;
            int cell_if_b = -1;
            bool iface_reversed = false;
            std::vector<int32_t> touching_face_verts;

            if (has_curved && tdim == 2 && ctype == cell::type::triangle && nv == 3)
            {
                for (int ze = 0; ze < n_zero && touching_ze < 0; ++ze)
                {
                    if (mesh.zero_entity_dim[static_cast<std::size_t>(ze)] != codim1_dim)
                        continue;
                    if ((mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] & zero_mask) == 0)
                        continue;

                    const int if_a = mesh.zero_entity_endpoint_v0[static_cast<std::size_t>(ze)];
                    const int if_b = mesh.zero_entity_endpoint_v1[static_cast<std::size_t>(ze)];
                    if (if_a < 0 || if_b < 0) continue;

                    for (int e = 0; e < nv; ++e)
                    {
                        const int ea = mesh.cell_vertices[static_cast<std::size_t>(c0 + e)];
                        const int eb = mesh.cell_vertices[static_cast<std::size_t>(c0 + (e + 1) % nv)];
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
            else if (has_curved && tdim == 3 && ctype == cell::type::tetrahedron)
            {
                std::vector<int32_t> cell_verts(static_cast<std::size_t>(nv));
                for (int j = 0; j < nv; ++j)
                    cell_verts[static_cast<std::size_t>(j)]
                        = mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];

                for (int ze = 0; ze < n_zero && touching_ze < 0; ++ze)
                {
                    if (mesh.zero_entity_dim[static_cast<std::size_t>(ze)] != codim1_dim)
                        continue;
                    if ((mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] & zero_mask) == 0)
                        continue;

                    const int z0 = mesh.zero_entity_offsets[static_cast<std::size_t>(ze)];
                    const int z1 = mesh.zero_entity_offsets[static_cast<std::size_t>(ze + 1)];
                    const int n_fv = z1 - z0;
                    if (n_fv != 3)
                        continue;

                    bool all_in_cell = true;
                    std::vector<int32_t> fv(static_cast<std::size_t>(n_fv));
                    for (int i = 0; i < n_fv; ++i)
                    {
                        fv[static_cast<std::size_t>(i)]
                            = mesh.zero_entity_vertices[static_cast<std::size_t>(z0 + i)];
                        bool found = false;
                        for (int j = 0; j < nv; ++j)
                        {
                            if (fv[static_cast<std::size_t>(i)] == cell_verts[static_cast<std::size_t>(j)])
                            {
                                found = true;
                                break;
                            }
                        }
                        if (!found)
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

            // Non-curved or unsupported cell type: emit P1 geometry
            if (touching_ze < 0
                || volume_backend == mapping::CurvedMappingBackend::collapsed)
            {
                for (int ci = c0; ci < c1; ++ci)
                {
                    const int lv = mesh.cell_vertices[static_cast<std::size_t>(ci)];
                    const int32_t gv = gv_count++;
                    for (int d = 0; d < gdim; ++d)
                        result.vertex_coords.push_back(
                            mesh.vertex_x[static_cast<std::size_t>(lv * gdim + d)]);
                    result.connectivity.push_back(gv);
                }
                result.offsets.push_back(
                    static_cast<int32_t>(result.connectivity.size()));
                result.cell_types.push_back(ctype);
                continue;
            }

            if (tdim == 3 && ctype == cell::type::tetrahedron)
            {
                int apex_vid = -1;
                for (int j = 0; j < nv; ++j)
                {
                    const int lv = mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];
                    bool on_face = false;
                    for (int i = 0; i < static_cast<int>(touching_face_verts.size()); ++i)
                    {
                        if (lv == touching_face_verts[static_cast<std::size_t>(i)])
                        {
                            on_face = true;
                            break;
                        }
                    }
                    if (!on_face)
                    {
                        apex_vid = lv;
                        break;
                    }
                }
                if (apex_vid < 0)
                {
                    for (int ci = c0; ci < c1; ++ci)
                    {
                        const int lv = mesh.cell_vertices[static_cast<std::size_t>(ci)];
                        const int32_t gv = gv_count++;
                        for (int d = 0; d < gdim; ++d)
                            result.vertex_coords.push_back(
                                mesh.vertex_x[static_cast<std::size_t>(lv * gdim + d)]);
                        result.connectivity.push_back(gv);
                    }
                    result.offsets.push_back(
                        static_cast<int32_t>(result.connectivity.size()));
                    result.cell_types.push_back(ctype);
                    continue;
                }

                std::vector<cutcells::BernsteinCell<T>> face_coord_polys;
                build_zero_face_coordinate_polys(
                    mesh, touching_ze, mesh.curved_geometry_order,
                    face_coord_polys, &result.n_fallback_nodes);

                const auto& tpl = (n_sub <= 1)
                    ? cutcells::p1_template(cell::type::tetrahedron)
                    : cutcells::iso_p1_template(cell::type::tetrahedron, std::min(n_sub, 4));

                std::vector<int32_t> sample_node_ids(
                    static_cast<std::size_t>(tpl.n_vertices), int32_t(-1));
                const T* r_apex = &mesh.vertex_ref_x[
                    static_cast<std::size_t>(apex_vid * tdim)];

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
                    ref_to_phys_single(mesh, F_ref, x_phys);

                    const int32_t gv = gv_count++;
                    sample_node_ids[static_cast<std::size_t>(vi)] = gv;
                    for (int d = 0; d < gdim; ++d)
                        result.vertex_coords.push_back(x_phys[d]);
                }

                for (int child = 0; child < tpl.n_cells; ++child)
                {
                    const int base = child * tpl.vertices_per_cell;
                    for (int j = 0; j < tpl.vertices_per_cell; ++j)
                    {
                        const int lv = tpl.cell_connectivity[
                            static_cast<std::size_t>(base + j)];
                        result.connectivity.push_back(
                            sample_node_ids[static_cast<std::size_t>(lv)]);
                    }
                    result.offsets.push_back(
                        static_cast<int32_t>(result.connectivity.size()));
                    result.cell_types.push_back(tpl.child_cell_type);
                }
                continue;
            }

            // --- GH-mapped curved triangle with vis_subdivision ---
            // Identify vertices: interface endpoints and opposite vertex
            const int if_a = cell_if_a;
            const int if_b = cell_if_b;

            int opp_idx = -1;
            for (int j = 0; j < nv; ++j)
            {
                const int lv = mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];
                if (lv != if_a && lv != if_b)
                {
                    opp_idx = j;
                    break;
                }
            }
            if (opp_idx < 0)
            {
                // Fallback: emit P1
                for (int ci = c0; ci < c1; ++ci)
                {
                    const int lv = mesh.cell_vertices[static_cast<std::size_t>(ci)];
                    const int32_t gv = gv_count++;
                    for (int d = 0; d < gdim; ++d)
                        result.vertex_coords.push_back(
                            mesh.vertex_x[static_cast<std::size_t>(lv * gdim + d)]);
                    result.connectivity.push_back(gv);
                }
                result.offsets.push_back(
                    static_cast<int32_t>(result.connectivity.size()));
                result.cell_types.push_back(ctype);
                continue;
            }

            const int lv_opp = mesh.cell_vertices[static_cast<std::size_t>(c0 + opp_idx)];
            const T* r_va = &mesh.vertex_ref_x[static_cast<std::size_t>(if_a * tdim)];
            const T* r_vb = &mesh.vertex_ref_x[static_cast<std::size_t>(if_b * tdim)];
            const T* r_opp = &mesh.vertex_ref_x[static_cast<std::size_t>(lv_opp * tdim)];

            // Interior curved nodes
            const int node_start = mesh.curved_zero_offsets[static_cast<std::size_t>(touching_ze)];
            const int node_end   = mesh.curved_zero_offsets[static_cast<std::size_t>(touching_ze + 1)];
            const int n_int = node_end - node_start;
            const int go = mesh.curved_geometry_order;

            std::vector<T> x_int_ref(static_cast<std::size_t>(n_int * tdim));
            for (int ni = 0; ni < n_int; ++ni)
            {
                const int src_ni = iface_reversed ? (n_int - 1 - ni) : ni;
                const int ci_node = node_start + src_ni;
                const T* xr = &mesh.curved_zero_ref_nodes[
                    static_cast<std::size_t>(ci_node * tdim)];
                T* rp = &x_int_ref[static_cast<std::size_t>(ni * tdim)];
                for (int d = 0; d < tdim; ++d)
                    rp[d] = xr[d];
            }

            const auto gl_pts_vec = ::cutcells::detail::gauss_lobatto_interior_points_1d<T>(go);
            const T* gl_pts = gl_pts_vec.data();

            // Sample a grid of (n_sub+1) x (n_sub+1) points on the (xi, eta)
            // triangle, emit sub-triangles.
            // Grid nodes: (i, j) with i + j <= n_sub
            // xi = i / n_sub, eta = j / n_sub (in collapsed coords)
            const int n1 = n_sub + 1;

            // Emit all grid vertices via GH map
            auto grid_vertex = [&](int i, int j) -> int32_t {
                // Compute (xi, eta) = (s, t) for the GH mapping
                const T xi  = static_cast<T>(i) / static_cast<T>(n_sub);
                const T eta = static_cast<T>(j) / static_cast<T>(n_sub);

                T F_ref[3] = {};
                T dF_dxi[3] = {};
                T dF_deta[3] = {};
                eval_gh_triangle_one_curved_edge(
                    xi, eta,
                    r_va, r_vb,
                    x_int_ref.data(), gl_pts, go,
                    r_opp, tdim,
                    F_ref, dF_dxi, dF_deta);

                // Push forward to physical space
                T x_phys[3] = {};
                ref_to_phys_single(mesh, F_ref, x_phys);

                const int32_t gv = gv_count++;
                for (int d = 0; d < gdim; ++d)
                    result.vertex_coords.push_back(x_phys[d]);
                return gv;
            };

            // Generate all grid vertices
            // Node index: row j, col i, with i + j <= n_sub
            // Total nodes: (n_sub+1)*(n_sub+2)/2
            std::vector<int32_t> node_ids;
            node_ids.reserve(static_cast<std::size_t>(n1 * (n_sub + 2) / 2));
            for (int j = 0; j < n1; ++j)
                for (int i = 0; i + j < n1; ++i)
                    node_ids.push_back(grid_vertex(i, j));

            // Index into node_ids: (i, j) → row_offset[j] + i
            auto node_id = [&](int i, int j) -> int32_t {
                // Row j starts at index j*(n_sub+1) - j*(j-1)/2
                int row_start = j * n1 - j * (j - 1) / 2;
                return node_ids[static_cast<std::size_t>(row_start + i)];
            };

            // Emit sub-triangles
            for (int j = 0; j < n_sub; ++j)
            {
                for (int i = 0; i + j < n_sub; ++i)
                {
                    // Lower triangle: (i,j), (i+1,j), (i,j+1)
                    result.connectivity.push_back(node_id(i, j));
                    result.connectivity.push_back(node_id(i + 1, j));
                    result.connectivity.push_back(node_id(i, j + 1));
                    result.offsets.push_back(
                        static_cast<int32_t>(result.connectivity.size()));
                    result.cell_types.push_back(cell::type::triangle);

                    // Upper triangle: (i+1,j), (i+1,j+1), (i,j+1)
                    if (i + 1 + j < n_sub)
                    {
                        result.connectivity.push_back(node_id(i + 1, j));
                        result.connectivity.push_back(node_id(i + 1, j + 1));
                        result.connectivity.push_back(node_id(i, j + 1));
                        result.offsets.push_back(
                            static_cast<int32_t>(result.connectivity.size()));
                        result.cell_types.push_back(cell::type::triangle);
                    }
                }
            }
        }
    }

    return result;
}

// ============================================================================
// Explicit instantiations
// ============================================================================

// Cache-reading assembler
template CurvedGlobalMesh<float> assemble_curved_interface_mesh<float>(
    std::span<const LocalMesh<float>*>, int, int);
template CurvedGlobalMesh<double> assemble_curved_interface_mesh<double>(
    std::span<const LocalMesh<double>*>, int, int);

// Convenience wrapper
template CurvedGlobalMesh<float> assemble_curved_interface_mesh<float, int>(
    std::span<LocalMesh<float>*>, const LevelSetFunction<float, int>&,
    LocalLevelSetBackend, int, int, float);
template CurvedGlobalMesh<double> assemble_curved_interface_mesh<double, int>(
    std::span<LocalMesh<double>*>, const LevelSetFunction<double, int>&,
    LocalLevelSetBackend, int, int, double);

template CurvedGlobalMesh<float> assemble_curved_volume_mesh<float>(
    std::span<const LocalMesh<float>*>, int, int);
template CurvedGlobalMesh<double> assemble_curved_volume_mesh<double>(
    std::span<const LocalMesh<double>*>, int, int);

template CurvedGlobalMesh<float> assemble_curved_volume_mesh<float>(
    std::span<const LocalMesh<float>*>, int, int,
    mapping::CurvedMappingBackend, int);
template CurvedGlobalMesh<double> assemble_curved_volume_mesh<double>(
    std::span<const LocalMesh<double>*>, int, int,
    mapping::CurvedMappingBackend, int);

} // namespace cutcells::mesh
