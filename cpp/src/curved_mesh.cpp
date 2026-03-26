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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace
{

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
                    result.connectivity.push_back(get_or_add_endpoint(lv1));
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
                    const int z0 = mesh.zero_entity_offsets[static_cast<std::size_t>(ze)];
                    const int z1 = mesh.zero_entity_offsets[static_cast<std::size_t>(ze + 1)];
                    for (int i = z0; i < z1; ++i)
                    {
                        const int lv = mesh.zero_entity_vertices[static_cast<std::size_t>(i)];
                        for (int d = 0; d < gdim; ++d)
                            result.vertex_coords.push_back(mesh.vertex_x[static_cast<std::size_t>(lv * gdim + d)]);
                        result.connectivity.push_back(gv_count++);
                    }
                    result.offsets.push_back(static_cast<int32_t>(result.connectivity.size()));
                    result.cell_types.push_back((z1 - z0) == 3 ? cell::type::triangle : cell::type::quadrilateral);
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
// Explicit instantiations
// ============================================================================

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

} // namespace cutcells::mesh
