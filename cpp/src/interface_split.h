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

namespace cutcells
{

/// Interface-adapted macro split for a 1-vs-3 tetrahedron.
///
/// For a tet with apex v0 (lone sign) and opposite face (v1,v2,v3):
///   - compute roots r01, r02, r03 on the 3 crossing edges (on interface)
///   - compute interior root c on the ray v0 → centroid(v1,v2,v3) (on interface)
///   - green-refine each intersected face using its 2 edge roots
///   - cone from c to all boundary sub-triangles → 10 child tets
///
/// The interface surface is represented by 3 triangles:
///   (c, r01, r02), (c, r02, r03), (c, r01, r03)
/// giving a refined piecewise-linear approximation of a curved interface.
///
/// Apex-side tets (3): (c, v0, r0i, r0j) — contain v0
/// Base-side tets (7): all vertices phi ≤ 0 — immediately clean
///
/// @param mesh          local mesh (modified in place)
/// @param cell_id       index of the tet to split
/// @param ls            level-set function (callable)
/// @param level_set_id  index of the active level set
/// @param tol           tolerance for zero classification and root finding
/// @return true if the split was performed
template <std::floating_point T>
bool interface_split_topology1_tet(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id = 0,
    T tol = static_cast<T>(1e-13));

/// Iteratively apply the interface-adapted split to all 1-vs-3 tets.
///
/// At each depth level, every tet with a 1-vs-3 sign pattern is split
/// into 10 children using interface-adapted root vertices.  Edge root
/// vertices shared between adjacent cells are deduplicated internally.
///
/// @param mesh          local mesh (modified in place)
/// @param ls            level-set function (callable)
/// @param level_set_id  index of the active level set
/// @param max_depth     maximum number of refinement passes
/// @param tol           tolerance for zero classification and root finding
/// @return total number of cells split across all depths
template <std::floating_point T>
int interface_refine_topology1(
    LocalMesh<T>& mesh,
    const LevelSetFunction<T>& ls,
    int level_set_id = 0,
    int max_depth = 5,
    T tol = static_cast<T>(1e-13));

/// Interface-adapted macro split for a 2-vs-2 tetrahedron.
///
/// For a tet with 2 positive vertices (p0,p1) and 2 negative (n0,n1):
///   - compute roots on the 4 crossing edges: (p0,n0), (p0,n1), (p1,n0), (p1,n1)
///   - compute interior root c on a ray through the centroid of the 4 edge roots
///   - each face (all 4 are intersected, 1-vs-2 pattern) is green-refined
///   - cone from c to all 12 boundary sub-triangles → 12 child tets
///
/// The interface is represented by 4 triangles (fan from c through the quad):
///   (c, r_p0n0, r_p1n0), (c, r_p0n1, r_p1n1),
///   (c, r_p0n0, r_p0n1), (c, r_p1n0, r_p1n1)
///
/// Positive-side tets (6): all vertices phi >= 0 — clean
/// Negative-side tets (6): all vertices phi <= 0 — clean
///
/// @param mesh          local mesh (modified in place)
/// @param cell_id       index of the tet to split
/// @param ls            level-set function (callable)
/// @param level_set_id  index of the active level set
/// @param tol           tolerance for zero classification and root finding
/// @return true if the split was performed
template <std::floating_point T>
bool interface_split_topology2_tet(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id = 0,
    T tol = static_cast<T>(1e-13));

/// Iteratively apply the interface-adapted split to all intersected tets.
///
/// Handles both 1-vs-3 (→10 children) and 2-vs-2 (→12 children)
/// topologies in each pass.  Edge root vertices are shared across cells
/// via an internal cache.
///
/// @param mesh          local mesh (modified in place)
/// @param ls            level-set function (callable)
/// @param level_set_id  index of the active level set
/// @param max_depth     maximum number of refinement passes
/// @param tol           tolerance for zero classification and root finding
/// @return total number of cells split across all depths
template <std::floating_point T>
int interface_refine(
    LocalMesh<T>& mesh,
    const LevelSetFunction<T>& ls,
    int level_set_id = 0,
    int max_depth = 5,
    T tol = static_cast<T>(1e-13));

// ============================================================================
// Face-enriched variants
// ============================================================================

/// Interface-adapted split for a 1-vs-3 tetrahedron with face-interior vertices.
///
/// Extends interface_split_topology1_tet by adding one interface vertex inside
/// each of the 3 intersected lateral faces (f01, f12, f20).  Each intersected
/// face is subdivided into 5 sub-triangles instead of 3, giving better
/// anchoring of the interface for non-linear (curved) level sets.
///
/// Total: 3 × 5 + 1 base = 16 child tetrahedra.
/// New vertices beyond the 10-tet scheme: f01, f12, f20.
///
/// @param mesh          local mesh (modified in place)
/// @param cell_id       index of the tet to split
/// @param ls            level-set function (callable)
/// @param level_set_id  index of the active level set
/// @param tol           tolerance for root finding
/// @return true if the split was performed
template <std::floating_point T>
bool interface_split_topology1_tet_face(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id = 0,
    T tol = static_cast<T>(1e-13));

/// Interface-adapted split for a 2-vs-2 tetrahedron with face-interior vertices.
///
/// Extends interface_split_topology2_tet by adding one interface vertex inside
/// each of the 4 intersected faces (fpn0, fpn1, fn0p0, fn0p1).  Each face is
/// subdivided into 5 sub-triangles, providing better anchoring of the interface
/// for non-linear (curved) level sets.
///
/// Total: 4 × 5 = 20 child tetrahedra.
/// New vertices beyond the 12-tet scheme: fpn0, fpn1, fn0p0, fn0p1.
///
/// @param mesh          local mesh (modified in place)
/// @param cell_id       index of the tet to split
/// @param ls            level-set function (callable)
/// @param level_set_id  index of the active level set
/// @param tol           tolerance for root finding
/// @return true if the split was performed
template <std::floating_point T>
bool interface_split_topology2_tet_face(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id = 0,
    T tol = static_cast<T>(1e-13));

/// Iteratively apply the face-enriched interface split to all intersected tets.
///
/// Handles both 1-vs-3 (→16 children) and 2-vs-2 (→20 children) topologies.
/// Equivalent to interface_refine but each intersected face also receives an
/// interior interface vertex, preventing non-linear level sets from
/// re-intersecting child tet faces.
///
/// @param mesh          local mesh (modified in place)
/// @param ls            level-set function (callable)
/// @param level_set_id  index of the active level set
/// @param max_depth     maximum number of refinement passes
/// @param tol           tolerance for root finding
/// @return total number of cells split across all depths
template <std::floating_point T>
int interface_refine_face(
    LocalMesh<T>& mesh,
    const LevelSetFunction<T>& ls,
    int level_set_id = 0,
    int max_depth = 5,
    T tol = static_cast<T>(1e-13));

/// Median-based face-enriched split for a 1-vs-3 tetrahedron.
///
/// Extends interface_split_topology1_tet by inserting:
///   - one midpoint on the opposite edge of each intersected lateral face,
///   - one interface root on the median from the lone-sign vertex to that midpoint,
///   - one interior interface root as in the existing 10-child split.
///
/// Each intersected lateral face contributes 6 boundary triangles and the base
/// face contributes 1, so coning from the interior root yields 19 children.
///
/// The split is rejected if post-split edge reclassification leaves any child
/// edge in {single_cross, multi_cross, uncertain}.
template <std::floating_point T>
bool interface_split_topology1_tet_faces(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id = 0,
    T tol = static_cast<T>(1e-13));

/// Median-based face-enriched split for a 2-vs-2 tetrahedron.
///
/// Each of the 4 intersected faces receives one midpoint and one face-interface
/// root on the median from the lone-sign face vertex to the opposite-edge
/// midpoint. Each face contributes 6 boundary triangles, so coning from the
/// interior root yields 24 children.
///
/// The split is rejected if post-split edge reclassification leaves any child
/// edge in {single_cross, multi_cross, uncertain}.
template <std::floating_point T>
bool interface_split_topology2_tet_faces(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id = 0,
    T tol = static_cast<T>(1e-13));

/// Iteratively apply the median-based face-enriched split to all intersected
/// tetrahedra with 1-vs-3 or 2-vs-2 sign topology.
template <std::floating_point T>
int interface_refine_faces(
    LocalMesh<T>& mesh,
    const LevelSetFunction<T>& ls,
    int level_set_id = 0,
    int max_depth = 5,
    T tol = static_cast<T>(1e-13));

} // namespace cutcells
