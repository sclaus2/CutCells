// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <span>
#include <string_view>
#include <vector>

#include "curved_mesh.h"
#include "cut_mesh.h"
#include "level_set.h"
#include "local_mesh.h"
#include "mesh_view.h"

namespace cutcells::mesh
{

template <std::floating_point T>
struct CurvedCutMeshResult
{
    CurvedGlobalMesh<T> inside;
    CurvedGlobalMesh<T> interface;
    CurvedGlobalMesh<T> outside;

    CutMesh<T> inside_vis;
    CutMesh<T> interface_vis;
    CutMesh<T> interface_curved_vis;
    CutMesh<T> outside_vis;

    std::vector<LocalMesh<T>> local_meshes;

    int n_parent_inside = 0;
    int n_parent_outside = 0;
    int n_parent_intersected = 0;
};

template <std::floating_point T>
CurvedCutMeshResult<T> cut_mesh_view_curved(
    const MeshView<T, int>&         mesh,
    const LevelSetFunction<T, int>& level_set,
    int                             geom_order = 4,
    std::string_view                backend = "bernstein",
    int                             vis_subdivision = 3,
    T                               tol = static_cast<T>(1e-12),
    bool                            repair_diagonals = false,
    int                             max_repair_depth = 3);

template <std::floating_point T>
CurvedCutMeshResult<T> cut_vtk_mesh_curved(
    std::span<const T>              points,
    std::span<const int>            connectivity,
    std::span<const int>            offset,
    std::span<const int>            vtk_type,
    const LevelSetFunction<T, int>& level_set,
    int                             geom_order = 4,
    std::string_view                backend = "bernstein",
    int                             vis_subdivision = 3,
    T                               tol = static_cast<T>(1e-12),
    bool                            repair_diagonals = false,
    int                             max_repair_depth = 3);

} // namespace cutcells::mesh
