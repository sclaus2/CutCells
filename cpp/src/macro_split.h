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

/// Pick the apex vertex of a topology-1 tet (1-vs-3 sign pattern).
///
/// @param phi  level-set values at the 4 tet vertices
/// @return local index (0..3) of the lone-sign vertex, or -1 if the
///         tet does not have a 1-vs-3 pattern
template <std::floating_point T>
int pick_topology1_apex(std::span<const T> phi);

/// Split one topology-1 intersected tet into 6 subtets.
///
/// The opposite face of the apex vertex is subdivided by its edge
/// midpoints and centroid into 6 triangles, which are coned from the
/// apex to produce 6 child tetrahedra.  Each child inherits the
/// apex as the lone opposite-sign vertex.
///
/// @param mesh          local mesh (modified in place)
/// @param cell_id       index of the tet to split
/// @param ls            level-set function (callable)
/// @param level_set_id  index of the active level set
/// @param tol           tolerance for zero classification
/// @return true if the split was performed
template <std::floating_point T>
bool macro_split_topology1_tet(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id = 0,
    T tol = static_cast<T>(1e-13));

/// Iteratively apply the topology-1 macro split to all 1-vs-3 tets.
///
/// At each depth level, every tet with a 1-vs-3 sign pattern is split
/// into 6 children.  Midpoint vertices shared between adjacent cells
/// are deduplicated internally.
///
/// @param mesh          local mesh (modified in place)
/// @param ls            level-set function (callable)
/// @param level_set_id  index of the active level set
/// @param max_depth     maximum number of refinement passes
/// @param tol           tolerance for zero classification
/// @return total number of tets that were split
template <std::floating_point T>
int macro_refine_topology1(
    LocalMesh<T>& mesh,
    const LevelSetFunction<T>& ls,
    int level_set_id = 0,
    int max_depth = 5,
    T tol = static_cast<T>(1e-13));

} // namespace cutcells
