// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <span>
#include <vector>

#include "local_mesh.h"
#include "level_set.h"
#include "cell_flags.h"
#include "cell_types.h"
#include "edge_root.h"

namespace cutcells
{

// ============================================================================
// Ray-based interior-point split for tetrahedra
// ============================================================================
//
// For an intersected tetrahedron, pick the "odd" vertex (the one whose
// sign differs from the majority in the 1-vs-3 case, or the vertex
// with the largest |phi| in the 2-vs-2 case).  Shoot a ray from that
// vertex to the centroid of the opposite face.  Find the root of phi
// along the ray, creating a new interior vertex on the interface.
// Decompose the original tet into 4 child tets using this interior vertex.
//
// The 1-to-4 stellar subdivision of tet [v0,v1,v2,v3] with interior point p:
//   child 0: [p, v1, v2, v3]   (replace v0 with p)
//   child 1: [v0, p, v2, v3]   (replace v1 with p)
//   child 2: [v0, v1, p, v3]   (replace v2 with p)
//   child 3: [v0, v1, v2, p]   (replace v3 with p)

/// @brief Pick the apex vertex index for ray refinement of a tet.
///
/// 1-vs-3: the lone vertex (same as green refinement choice).
/// 2-vs-2: the vertex with the largest |phi| from among the pair whose
///         sign is the same as the vertex with globally largest |phi|.
///
/// @param phi   Level-set values at the 4 tet vertices (ordered as in cell_vertices).
/// @return      Local index (0..3) of the chosen apex vertex.
template <std::floating_point T>
int pick_ray_apex(std::span<const T> phi);

/// @brief Split one tet into 4 children using an interior point p.
///
/// The interior point p replaces each of the 4 vertices in turn,
/// producing 4 child tets.  Appends children to the output arrays.
///
/// @param cv     Cell vertex indices of the parent tet (size 4).
/// @param p      Global vertex index of the interior point.
/// @param out_cv Output cell_vertices (appended).
/// @param out_off Output cell_offsets (appended).
/// @param out_ct  Output cell_types (appended).
/// @return Number of children created (always 4).
int stellar_split_tetrahedron(
    std::span<const int32_t> cv,
    int32_t p,
    std::vector<int32_t>& out_cv,
    std::vector<int32_t>& out_off,
    std::vector<cell::type>& out_ct);

/// @brief Perform ray-based interior-point refinement on a single
///        intersected tet in the local mesh.
///
/// Finds the root of the level-set function along the ray from the
/// chosen apex vertex to the centroid of the opposite face, adds a new
/// interior vertex at the root location, and replaces the original tet
/// with 4 child tets.
///
/// @param mesh        Local mesh (modified in place).
/// @param cell_id     Index of the tet to split.
/// @param ls          Level-set function (callable).
/// @param level_set_id Which level set to use.
/// @param tol         Root-finding tolerance.
/// @return True on success, false if no root is found on the ray.
template <std::floating_point T>
bool ray_split_one_tet(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol = static_cast<T>(1e-13));

/// @brief Iteratively apply ray-based refinement to all intersected tets
///        in the local mesh, up to a maximum depth.
///
/// At each level, every intersected tet (according to vertex_phi signs)
/// is split into 4 children.  Edges and faces are rebuilt after each level.
///
/// @param mesh         Local mesh (modified in place).
/// @param ls           Level-set function (callable).
/// @param level_set_id Which level set to use.
/// @param max_depth    Maximum number of refinement levels.
/// @param tol          Root-finding tolerance.
/// @return Number of ray splits actually performed across all levels.
template <std::floating_point T>
int ray_refine_local_mesh(
    LocalMesh<T>& mesh,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    int max_depth,
    T tol = static_cast<T>(1e-13));

} // namespace cutcells
