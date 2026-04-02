// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <functional>
#include <span>
#include <vector>

#include "cell_topology.h"
#include "cell_types.h"
#include "cut_cell.h"

namespace cutcells::cell
{

// ============================================================================
// Result struct
// ============================================================================

/// Diagnostic output of a triangulation diagonal repair pass on a CutCell.
struct TriangulationRepairInfo
{
    int checked_diagonals = 0;   ///< diagonals inspected
    int invalid_diagonals = 0;   ///< diagonals that cross the interface
    int swapped_diagonals = 0;   ///< fixed by swapping the diagonal
    int unresolved = 0;          ///< could not be fixed (needs green refinement)
};

// ============================================================================
// Free functions — diagonal classification
// ============================================================================

/// Test whether an edge in a CutCell is a triangulation diagonal
/// (a new volume edge introduced by decomposing a quad or prism into
/// simplices), as opposed to an original parent-cell edge, a root-split
/// edge, or an interface edge.
///
/// Uses the vertex_parent_entity tokens of the CutCell.
///
/// @param parent_cell_type   type of the cell that was cut
/// @param token_a            vertex_parent_entity of the first endpoint
/// @param token_b            vertex_parent_entity of the second endpoint
/// @return true if the edge is a triangulation diagonal
bool is_triangulation_diagonal(type parent_cell_type,
                               int32_t token_a,
                               int32_t token_b);

// ============================================================================
// Free functions — midpoint check
// ============================================================================

/// Check whether the midpoint of a CutCell edge lies on the wrong side
/// of the level-set interface.
///
/// @param cut_cell        the CutCell with vertex coordinates
/// @param va, vb          local vertex indices within the CutCell
/// @param eval_phi        level-set evaluator: T(const T* x_phys)
/// @param expected_sign   expected sign of phi in this domain
///                        (-1 for phi<0, +1 for phi>0)
/// @param tol             tolerance for zero classification
/// @return true if the midpoint contradicts the expected domain
template <std::floating_point T>
bool diagonal_crosses_interface(
    const CutCell<T>& cut_cell,
    int va, int vb,
    std::function<T(const T*)> eval_phi,
    int expected_sign,
    T tol = static_cast<T>(1e-14));

// ============================================================================
// Free functions — 2D diagonal swap
// ============================================================================

/// Swap the shared diagonal between two triangles in a CutCell.
///
/// Given sub-cells cell_i and cell_j that are both triangles sharing
/// edge (a, b):
///   T_i = {a, b, c}   T_j = {a, b, d}
/// Replaces with:
///   T_i' = {a, c, d}  T_j' = {b, d, c}
///
/// @return true if the swap was performed
template <std::floating_point T>
bool swap_diagonal_2d(CutCell<T>& cut_cell, int cell_i, int cell_j);

// ============================================================================
// Free functions — 3D alternative prism triangulation
// ============================================================================

/// Replace the triangulation of 3 tets from a prism with an alternative
/// decomposition that uses different quad-face diagonals.
///
/// Current decomposition (A):
///   {v0,v2,v1,v3}, {v1,v3,v5,v4}, {v1,v2,v5,v3}
///
/// Alternative decomposition (B):
///   {v0,v1,v2,v4}, {v0,v4,v2,v5}, {v0,v4,v5,v3}
///
/// The 3 tets must be consecutive starting at cell_start (cell_start,
/// cell_start+1, cell_start+2).
///
/// @param cut_cell    the CutCell to modify
/// @param cell_start  index of the first of the 3 consecutive tets
/// @param prism_verts ordered prism vertices [v0..v5] as they appear
///                    in the CutCell local numbering
/// @return true if replacement was performed
template <std::floating_point T>
bool swap_prism_triangulation_3d(CutCell<T>& cut_cell,
                                 int cell_start,
                                 std::span<const int> prism_verts);

// ============================================================================
// Free functions — combined repair
// ============================================================================

/// Check all triangulation diagonals in a CutCell and try to repair
/// by swapping.
///
/// For 2D (quad → 2 triangles): tries the alternative diagonal.
/// For 3D (prism → 3 tets): tries the alternative prism decomposition.
/// If neither alternative works, marks the diagonal as unresolved
/// (requiring green refinement or child cutting).
///
/// @param cut_cell          modified in place if swaps are made
/// @param parent_cell_type  the parent cell type
/// @param eval_phi          level-set evaluator: T(const T* x_phys)
/// @param expected_sign     -1 for phi<0 domain, +1 for phi>0 domain
/// @param tol               tolerance for zero classification
/// @param debug             if true, print diagnostic messages to stderr
/// @return repair info
template <std::floating_point T>
TriangulationRepairInfo repair_cut_cell_diagonals(
    CutCell<T>& cut_cell,
    type parent_cell_type,
    std::function<T(const T*)> eval_phi,
    int expected_sign,
    T tol = static_cast<T>(1e-14),
    bool debug = false);

} // namespace cutcells::cell
