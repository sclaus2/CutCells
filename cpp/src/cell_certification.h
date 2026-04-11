// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <span>
#include <vector>

#include "adapt_cell.h"
#include "cell_types.h"
#include "level_set_cell.h"

namespace cutcells
{

// =====================================================================
// Exact subcell Bernstein restriction
// =====================================================================

/// Compute the Bernstein coefficients of a parent polynomial restricted
/// to a leaf subcell.
///
/// The subcell is defined by its vertex reference coordinates (in the
/// parent reference frame). The restriction is exact: sample at
/// (degree+1)^tdim control points, then invert the Bernstein matrix.
///
/// @param parent_cell_type        Type of the parent cell.
/// @param degree                  Polynomial degree.
/// @param parent_coeffs           Bernstein coefficients on the parent cell.
/// @param subcell_type            Type of the subcell (triangle, tet, quad, ...).
/// @param subcell_vertex_coords   Reference coordinates of subcell vertices (flat: nverts * tdim).
/// @param[out] subcell_coeffs     Bernstein coefficients on the subcell.
template <std::floating_point T>
void restrict_subcell_bernstein_exact(cell::type parent_cell_type,
                                      int degree,
                                      std::span<const T> parent_coeffs,
                                      cell::type subcell_type,
                                      std::span<const T> subcell_vertex_coords,
                                      std::vector<T>& subcell_coeffs);

// =====================================================================
// Cell classifier
// =====================================================================

/// Classify a single leaf cell for one level set.
///
/// Logic:
///   A. If the incident edge pattern is a directly cuttable simplex case:
///        - triangle with exactly 2 one_root edges and no multiple_roots/zero
///        - tetrahedron with exactly 3 or 4 one_root edges and no multiple_roots/zero
///      → ready_to_cut.
///   B. Else if any incident edge has tag one_root, multiple_roots, or zero → cut.
///   B. Otherwise, restrict the parent Bernstein to the subcell and check
///      the sign hull:
///        - all positive → positive
///        - all negative → negative
///        - all zero → zero
///        - mixed → ambiguous
///
/// @param adapt_cell     The AdaptCell.
/// @param ls_cell        LevelSetCell providing Bernstein coefficients.
/// @param level_set_id   Which level set.
/// @param cell_id        Index of the leaf cell in entity_to_vertex[tdim].
/// @param zero_tol       Tolerance for all-zero.
/// @param sign_tol       Tolerance for all-positive / all-negative.
/// @return CellCertTag.
template <std::floating_point T, std::integral I>
CellCertTag classify_leaf_cell(const AdaptCell<T>& adapt_cell,
                               const LevelSetCell<T, I>& ls_cell,
                               int level_set_id,
                               int cell_id,
                               T zero_tol, T sign_tol);

/// Classify all not-yet-classified leaf cells for one level set.
///
/// @param adapt_cell     The AdaptCell (modified in place).
/// @param ls_cell        LevelSetCell providing Bernstein coefficients.
/// @param level_set_id   Which level set.
/// @param zero_tol       Tolerance for all-zero.
/// @param sign_tol       Tolerance for all-positive / all-negative.
template <std::floating_point T, std::integral I>
void classify_leaf_cells(AdaptCell<T>& adapt_cell,
                         const LevelSetCell<T, I>& ls_cell,
                         int level_set_id,
                         T zero_tol, T sign_tol);

/// Evaluate the level set on every current AdaptCell vertex and update the
/// sign/zero masks for that level set id.
template <std::floating_point T, std::integral I>
void fill_all_vertex_signs_from_level_set(AdaptCell<T>& adapt_cell,
                                          const LevelSetCell<T, I>& ls_cell,
                                          int level_set_id,
                                          T zero_tol);

/// Replace leaf cells marked ready_to_cut by the non-triangulated LUT cut
/// decomposition on the positive and negative side.
template <std::floating_point T, std::integral I>
void process_ready_to_cut_cells(AdaptCell<T>& adapt_cell,
                                const LevelSetCell<T, I>& ls_cell,
                                int level_set_id,
                                T zero_tol,
                                T sign_tol,
                                int edge_max_depth);

// =====================================================================
// Top-level certification + refinement driver
// =====================================================================

/// Iterative certification and refinement loop.
///
/// Per iteration:
///   1. classify_new_edges(...)
///   2. classify_leaf_cells(...)
///   3. If any multiple_roots edges → green refine (priority)
///   4. Else if any ambiguous cells → red refine
///   5. Stop if no refinement occurred.
///
/// @param adapt_cell       The AdaptCell (modified in place).
/// @param ls_cell          LevelSetCell providing Bernstein coefficients.
/// @param level_set_id     Which level set.
/// @param max_iterations   Maximum number of refine-reclassify iterations.
/// @param zero_tol         Tolerance for all-zero.
/// @param sign_tol         Tolerance for all-positive / all-negative.
/// @param edge_max_depth   Maximum subdivision depth for edge root search.
template <std::floating_point T, std::integral I>
void certify_and_refine(AdaptCell<T>& adapt_cell,
                        const LevelSetCell<T, I>& ls_cell,
                        int level_set_id,
                        int max_iterations,
                        T zero_tol, T sign_tol,
                        int edge_max_depth);

/// Full local single-level-set pipeline:
///   1. stamp vertex signs
///   2. certify + green/red refine until stable
///   3. replace ready_to_cut cells by the LUT decomposition
///   4. restamp vertex signs on the final leaf mesh
template <std::floating_point T, std::integral I>
void certify_refine_and_process_ready_cells(AdaptCell<T>& adapt_cell,
                                            const LevelSetCell<T, I>& ls_cell,
                                            int level_set_id,
                                            int max_iterations,
                                            T zero_tol, T sign_tol,
                                            int edge_max_depth);

} // namespace cutcells
