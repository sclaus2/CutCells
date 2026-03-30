// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <span>
#include <vector>

#include "cell_types.h"
#include "level_set.h"
#include "local_mesh.h"
#include "mapping_curved.h"

namespace cutcells::mesh
{

/// A global mesh with Pk (curved) geometry nodes.
///
/// Stores physical coordinates of all geometry nodes (corners + high-order).
/// Connectivity uses VTK Lagrange node ordering.
template <std::floating_point T>
struct CurvedGlobalMesh
{
    std::vector<T>          vertex_coords;  // physical, flat, n_vertices * gdim
    std::vector<int32_t>    connectivity;   // CSR flat (VTK Lagrange ordering)
    std::vector<int32_t>    offsets;        // size = n_cells + 1
    std::vector<cell::type> cell_types;
    int gdim = 0;
    int geom_order = 1;
    int level_set_id = -1;
    int n_fallback_nodes = 0;

    int n_vertices() const
    {
        return gdim > 0 ? static_cast<int>(vertex_coords.size()) / gdim : 0;
    }
    int n_cells() const
    {
        return offsets.empty() ? 0 : static_cast<int>(offsets.size()) - 1;
    }
};

/// Assemble the curved zero-interface mesh (phi = 0 surface) from all
/// LocalMeshes. Reads curved GL nodes from the per-entity cache built by
/// build_zero_entities + curve_zero_entities_with_backend.
/// Does NOT re-project — only reads from the cache.
///
/// In 2D the interface consists of curved line segments (VTK Lagrange curves).
/// In 3D it consists of curved triangles or quadrilaterals.
///
/// Prerequisite: each LocalMesh must have had
///   build_zero_entities(mesh, level_set_id),
///   curve_zero_entities_with_backend(mesh, ..., level_set_id, geom_order, ...),
///   build_zero_chains(mesh, zero_mask)   [2D], or
///   build_zero_patches(mesh, zero_mask)  [3D]
///   called before this function.
///
/// @param local_meshes  pointers to decomposed + curved LocalMeshes
/// @param level_set_id  which level set to extract the interface for
/// @param geom_order    Pk geometry order (>= 1)
template <std::floating_point T>
CurvedGlobalMesh<T> assemble_curved_interface_mesh(
    std::span<const LocalMesh<T>*> local_meshes,
    int                            level_set_id,
    int                            geom_order);

/// Convenience overload: builds and curves interface entities on each
/// LocalMesh if not already done, then reads from the cache.
///
/// @param local_meshes  pointers to LocalMeshes (may be mutated to build cache)
/// @param level_set     level-set function for projection
/// @param backend       evaluation backend
/// @param level_set_id  which level set to extract the interface for
/// @param geom_order    Pk geometry order (>= 2)
/// @param tol           convergence tolerance for root projection
template <std::floating_point T, std::integral I = int>
CurvedGlobalMesh<T> assemble_curved_interface_mesh(
    std::span<LocalMesh<T>*>           local_meshes,
    const LevelSetFunction<T, I>&      level_set,
    LocalLevelSetBackend               backend,
    int                                level_set_id,
    int                                geom_order,
    T                                  tol = static_cast<T>(1e-12));

/// Assemble the curved volume mesh (phi < 0 sub-cells) from all LocalMeshes.
///
/// Sub-cells that touch the interface incorporate curved edge/face nodes
/// via a linear blend mapping. Sub-cells away from the interface remain
/// straight (P1).
///
/// @param local_meshes  pointers to decomposed LocalMeshes with curved cache
/// @param level_set_id  which level set
/// @param geom_order    Pk order for curved elements
template <std::floating_point T>
CurvedGlobalMesh<T> assemble_curved_volume_mesh(
    std::span<const LocalMesh<T>*> local_meshes,
    int level_set_id,
    int geom_order);

/// Overload with explicit mapping backend selection.
///
/// When backend is gordon_hall, sub-cells touching the interface are
/// output as high-order elements with geometry nodes sampled from the
/// Gordon-Hall mapping. Otherwise, P1 geometry is emitted for all cells.
///
/// @param local_meshes     pointers to decomposed LocalMeshes with curved cache
/// @param level_set_id     which level set
/// @param geom_order       Pk order for curved elements
/// @param volume_backend   mapping backend for curved cells
/// @param vis_subdivision  number of visualisation subdivisions per edge (>= 1)
template <std::floating_point T>
CurvedGlobalMesh<T> assemble_curved_volume_mesh(
    std::span<const LocalMesh<T>*> local_meshes,
    int level_set_id,
    int geom_order,
    mapping::CurvedMappingBackend volume_backend,
    int vis_subdivision = 1);

} // namespace cutcells::mesh
