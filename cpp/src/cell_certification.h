// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "adapt_cell.h"
#include "cell_types.h"
#include "graph_criteria.h"
#include "level_set_cell.h"

namespace cutcells
{

enum class GraphProjectionMode : std::uint8_t
{
    straight_zero_entity_normal = 0,
    level_set_gradient = 1
};

enum class GraphRefinementMode : std::uint8_t
{
    green_edge = 0,
    red_failed_cell = 1
};

template <std::floating_point T>
struct ReadyCellGraphOptions
{
    bool enabled = true;
    int geometry_order = 2;
    int max_refinements = 5;
    GraphProjectionMode projection_mode =
        GraphProjectionMode::level_set_gradient;
    GraphRefinementMode refinement_mode =
        GraphRefinementMode::green_edge;
    graph_criteria::Options<T> criteria = {};
    T min_level_set_gradient_host_alignment = T(0.9);
};

enum class GraphNodeKind : std::uint8_t
{
    edge_endpoint = 0,
    edge_interior = 1,
    face_boundary_edge = 2,
    face_interior = 3
};

template <std::floating_point T>
struct GraphNodeDiagnostics
{
    int node_index = -1;
    GraphNodeKind node_kind = GraphNodeKind::edge_interior;
    bool accepted = true;
    bool fallback_used = false;
    graph_criteria::DirectionKind selected_direction_kind =
        graph_criteria::DirectionKind::projected_level_set_gradient;
    graph_criteria::FailureReason failure_reason =
        graph_criteria::FailureReason::none;
    int parent_entity_dim = -1;
    int parent_entity_id = -1;
    int requested_refinement_entity_dim = -1;
    int requested_refinement_entity_id = -1;

    T level_set_gradient_host_alignment = std::numeric_limits<T>::quiet_NaN();
    T level_set_gradient_angle_to_tangent_deg =
        std::numeric_limits<T>::quiet_NaN();
    T selected_host_alignment = std::numeric_limits<T>::quiet_NaN();
    T drift_amplification = std::numeric_limits<T>::quiet_NaN();
    T relative_correction_distance = std::numeric_limits<T>::quiet_NaN();
    T relative_tangential_shift = std::numeric_limits<T>::quiet_NaN();
    T true_transversality = std::numeric_limits<T>::quiet_NaN();

    std::vector<T> seed;
    std::vector<T> corrected;
    std::vector<T> selected_direction;
    std::vector<T> level_set_gradient_direction;
    std::vector<T> straight_helper_normal;
};

template <std::floating_point T>
struct ZeroEntityGraphDiagnostics
{
    int local_zero_entity_id = -1;
    int level_set_id = -1;
    int dimension = -1;
    std::uint64_t zero_mask = 0;
    bool accepted = true;
    int checked_edges = 0;
    int checked_faces = 0;
    int failed_checks = 0;
    graph_criteria::FailureReason failure_reason =
        graph_criteria::FailureReason::none;

    T min_true_transversality = std::numeric_limits<T>::infinity();
    T min_host_normal_alignment = std::numeric_limits<T>::infinity();
    T max_drift_amplification = T(0);
    T max_relative_correction_distance = T(0);
    T max_relative_tangential_shift = T(0);
    T min_edge_gap_ratio = std::numeric_limits<T>::infinity();
    T min_face_area_ratio = std::numeric_limits<T>::infinity();
    T min_level_set_gradient_host_alignment =
        std::numeric_limits<T>::infinity();
    int failed_face_triangle_index = -1;
    T failed_face_area_ratio = std::numeric_limits<T>::quiet_NaN();

    std::vector<T> failed_projection_seed;
    std::vector<T> failed_projection_direction;
    T failed_projection_clip_lo = std::numeric_limits<T>::quiet_NaN();
    T failed_projection_clip_hi = std::numeric_limits<T>::quiet_NaN();
    T failed_projection_root_t = std::numeric_limits<T>::quiet_NaN();

    int requested_refinement_entity_dim = -1;
    int requested_refinement_entity_id = -1;
    std::vector<GraphNodeDiagnostics<T>> nodes;
};

template <std::floating_point T>
struct ReadyCellGraphDiagnostics
{
    bool accepted = true;
    int checked_cells = 0;
    int checked_edges = 0;
    int checked_faces = 0;
    int failed_checks = 0;
    int graph_refinements = 0;
    int first_failed_cell = -1;
    graph_criteria::FailureReason first_failure_reason =
        graph_criteria::FailureReason::none;

    T min_true_transversality = std::numeric_limits<T>::infinity();
    T min_host_normal_alignment = std::numeric_limits<T>::infinity();
    T max_drift_amplification = T(0);
    T max_relative_correction_distance = T(0);
    T max_relative_tangential_shift = T(0);
    T min_edge_gap_ratio = std::numeric_limits<T>::infinity();
    T min_face_area_ratio = std::numeric_limits<T>::infinity();
    T min_level_set_gradient_host_alignment =
        std::numeric_limits<T>::infinity();
    int first_failed_face_triangle_index = -1;
    T first_failed_face_area_ratio = std::numeric_limits<T>::quiet_NaN();

    std::vector<T> first_failed_projection_seed;
    std::vector<T> first_failed_projection_direction;
    T first_failed_projection_clip_lo = std::numeric_limits<T>::quiet_NaN();
    T first_failed_projection_clip_hi = std::numeric_limits<T>::quiet_NaN();
    T first_failed_projection_root_t = std::numeric_limits<T>::quiet_NaN();
    int first_requested_refinement_entity_dim = -1;
    int first_requested_refinement_entity_id = -1;

    /// Diagnostics attached to the committed AdaptCell zero entities. These
    /// are populated after the linear cut is committed and before any lazy
    /// high-order curving can run, so visualization can display the graph
    /// criterion on the zero edge/face where it is evaluated.
    std::vector<ZeroEntityGraphDiagnostics<T>> zero_entities;
    std::vector<GraphNodeDiagnostics<T>> nodes;
};

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
///        - triangle with exactly 2 one_root edges and no multiple_roots
///        - tetrahedron with exactly 3 or 4 one_root edges and no multiple_roots
///      → ready_to_cut.
///   B. Else if any incident edge has tag one_root or multiple_roots → cut.
///      Zero edges alone fall through so later level sets can still classify
///      or refine cells whose interface is inherited from an earlier cut.
///   C. Otherwise, restrict the parent Bernstein to the subcell and check
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

// =====================================================================
// Face (facet) classifier (3D only)
// =====================================================================

/// Classify a single leaf face for one level set.
///
/// Logic mirrors classify_leaf_cell but operates on a 2D face entity
/// (dimension 2) of a 3D AdaptCell:
///   1. Check incident-edge root topology for a valid 2D cut pattern.
///   2. Restrict the parent Bernstein to the face and test the sign hull.
///   3. Apply monotonicity filter.
///   4. If the centroid is needed (all-zero case), evaluate there.
///   5. Otherwise → ambiguous.
///
/// @param adapt_cell     The AdaptCell.
/// @param ls_cell        LevelSetCell providing Bernstein coefficients.
/// @param level_set_id   Which level set.
/// @param face_id        Index of the leaf face in entity_to_vertex[2].
/// @param zero_tol       Tolerance for all-zero.
/// @param sign_tol       Tolerance for all-positive / all-negative.
/// @return FaceCertTag.
template <std::floating_point T, std::integral I>
FaceCertTag classify_leaf_face(const AdaptCell<T>& adapt_cell,
                               const LevelSetCell<T, I>& ls_cell,
                               int level_set_id,
                               int face_id,
                               T zero_tol, T sign_tol);

/// Classify all not-yet-classified leaf faces for one level set.
///
/// Face entities (entity_to_vertex[2]) must already be populated via
/// build_faces(). Only operates when tdim == 3.
///
/// @param adapt_cell     The AdaptCell (modified in place).
/// @param ls_cell        LevelSetCell providing Bernstein coefficients.
/// @param level_set_id   Which level set.
/// @param zero_tol       Tolerance for all-zero.
/// @param sign_tol       Tolerance for all-positive / all-negative.
template <std::floating_point T, std::integral I>
void classify_leaf_faces(AdaptCell<T>& adapt_cell,
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

/// Replace leaf cells marked ready_to_cut by the LUT cut decomposition on the
/// positive and negative side.
template <std::floating_point T, std::integral I>
void process_ready_to_cut_cells(AdaptCell<T>& adapt_cell,
                                const LevelSetCell<T, I>& ls_cell,
                                int level_set_id,
                                T zero_tol,
                                T sign_tol,
                                int edge_max_depth,
                                bool triangulate_cut_parts = false);

/// Run graph diagnostics on all current ready_to_cut cells without mutating
/// AdaptCell topology, masks, or provenance. The diagnostic cut is a temporary
/// P1/linear interface built from the current uncut ready leaf cells.
template <std::floating_point T, std::integral I>
ReadyCellGraphDiagnostics<T> check_ready_to_cut_cell_graphs(
    const AdaptCell<T>& adapt_cell,
    const LevelSetCell<T, I>& ls_cell,
    int level_set_id,
    const ReadyCellGraphOptions<T>& graph_options = {});

/// Populate per-committed-zero-entity graph diagnostics on an existing
/// diagnostic object. This is used after final zero-entity inventory rebuilds
/// so visualization data is keyed to the AdaptCell zero entities that users
/// inspect.
template <std::floating_point T, std::integral I>
void populate_committed_zero_entity_graph_diagnostics(
    const AdaptCell<T>& adapt_cell,
    const LevelSetCell<T, I>& ls_cell,
    int level_set_id,
    const ReadyCellGraphOptions<T>& graph_options,
    ReadyCellGraphDiagnostics<T>& diagnostics);

/// Green-refine a ready leaf cell through the midpoint of the cell edge whose
/// midpoint has the largest level-set value.
template <std::floating_point T, std::integral I>
bool refine_ready_cell_on_largest_midpoint_value(
    AdaptCell<T>& adapt_cell,
    const LevelSetCell<T, I>& ls_cell,
    int level_set_id,
    int cell_id);

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
                                            int edge_max_depth,
                                            bool triangulate_cut_parts = false);

/// Full local pipeline with graph preflight:
///   1. stamp vertex signs
///   2. certify + green/red refine until stable
///   3. build temporary P1 cuts for graph checks on ready cells
///   4. on graph failure, discard the temporary cut and green-refine the
///      original uncut ready cell, then recurse
///   5. replace ready_to_cut cells by the committed cut decomposition
template <std::floating_point T, std::integral I>
ReadyCellGraphDiagnostics<T> certify_refine_graph_check_and_process_ready_cells(
    AdaptCell<T>& adapt_cell,
    const LevelSetCell<T, I>& ls_cell,
    int level_set_id,
    int max_iterations,
    T zero_tol, T sign_tol,
    int edge_max_depth,
    bool triangulate_cut_parts = false,
    const ReadyCellGraphOptions<T>& graph_options = {});

} // namespace cutcells
