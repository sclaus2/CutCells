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
// Root-interval helper
// =====================================================================

/// A parameter interval [t0, t1] containing a root of a 1D Bernstein polynomial.
template <std::floating_point T>
struct EdgeRootInterval
{
    T t0 = 0;
    T t1 = 0;
};

// =====================================================================
// 1D Bernstein subdivision
// =====================================================================

/// Exact de Casteljau subdivision of a degree-p 1D Bernstein polynomial at t_split.
///
/// Given coefficients c_0, ..., c_p, computes left[0..p] and right[0..p]
/// such that:
///   left  is the Bernstein representation on [0, t_split]
///   right is the Bernstein representation on [t_split, 1]
///
/// @param coeffs   Bernstein coefficients of degree p (size p+1).
/// @param t_split  Parameter in (0,1).
/// @param[out] left   Left child coefficients (resized to p+1).
/// @param[out] right  Right child coefficients (resized to p+1).
template <std::floating_point T>
void subdivide_bernstein_1d(std::span<const T> coeffs,
                            T t_split,
                            std::vector<T>& left,
                            std::vector<T>& right);

// =====================================================================
// 1D sign-hull helpers
// =====================================================================

/// True if all coefficients are > tol.
template <std::floating_point T>
bool bernstein_all_positive(std::span<const T> coeffs, T tol);

/// True if all coefficients are < -tol.
template <std::floating_point T>
bool bernstein_all_negative(std::span<const T> coeffs, T tol);

/// True if |coeff| <= tol for every coefficient.
template <std::floating_point T>
bool bernstein_all_zero(std::span<const T> coeffs, T tol);

// =====================================================================
// Recursive root-interval finder
// =====================================================================

/// Recursively locate parameter intervals containing roots of a 1D
/// Bernstein polynomial, using convex-hull exclusion and de Casteljau
/// subdivision.
///
/// @param coeffs    Bernstein coefficients on the sub-interval [t0, t1].
/// @param t0, t1    Current parameter bounds (start with 0, 1).
/// @param zero_tol  Tolerance for the all-zero test.
/// @param sign_tol  Tolerance for the all-positive / all-negative test.
/// @param depth     Current recursion depth.
/// @param max_depth Maximum recursion depth.
/// @param[out] intervals     Collected root intervals.
/// @param[out] has_zero_segment  Set to true if an identically-zero segment is found.
template <std::floating_point T>
void find_root_intervals_1d(std::span<const T> coeffs,
                            T t0, T t1,
                            T zero_tol, T sign_tol,
                            int depth, int max_depth,
                            std::vector<EdgeRootInterval<T>>& intervals,
                            bool& has_zero_segment);

// =====================================================================
// Exact edge Bernstein extraction
// =====================================================================

/// Extract the Bernstein coefficients of a parent polynomial restricted to
/// a parent edge (local edge id, canonical numbering).
///
/// For simplex cells this is multi-index slicing.
/// For tensor-product cells this is row/column extraction.
///
/// @param parent_cell_type  Type of the parent cell.
/// @param degree            Polynomial degree.
/// @param parent_coeffs     Bernstein coefficients on the parent cell.
/// @param parent_local_edge_id  Local edge index (cell_topology.h ordering).
/// @param[out] edge_coeffs  Bernstein coefficients on the edge (degree+1 entries).
template <std::floating_point T>
void extract_parent_edge_bernstein(cell::type parent_cell_type,
                                   int degree,
                                   std::span<const T> parent_coeffs,
                                   int parent_local_edge_id,
                                   std::vector<T>& edge_coeffs);

/// Exact Bernstein restriction to an arbitrary edge [xi_a, xi_b] inside
/// the parent reference cell.
///
/// Computes the Bernstein coefficients of q(t) = phi((1-t)*xi_a + t*xi_b).
///
/// @param parent_cell_type  Type of the parent cell.
/// @param degree            Polynomial degree.
/// @param parent_coeffs     Bernstein coefficients on the parent cell.
/// @param xi_a, xi_b        Endpoints in parent reference coordinates (size = tdim).
/// @param[out] edge_coeffs  Bernstein coefficients on the edge (degree+1 entries).
template <std::floating_point T>
void restrict_edge_bernstein_exact(cell::type parent_cell_type,
                                   int degree,
                                   std::span<const T> parent_coeffs,
                                   std::span<const T> xi_a,
                                   std::span<const T> xi_b,
                                   std::vector<T>& edge_coeffs);

/// Test whether an adaptive edge lies entirely on a single parent edge.
///
/// Uses vertex provenance: both endpoints must have parent_dim == 1 (or 0
/// with the vertex on the same parent edge) and the same parent_edge_id.
///
/// @param adapt_cell       The AdaptCell.
/// @param edge_id          Index of the edge in entity_to_vertex[1].
/// @param[out] parent_edge_id  If true, the parent edge id.
/// @return true if the edge lies on a single parent edge.
template <std::floating_point T>
bool edge_is_on_single_parent_edge(const AdaptCell<T>& adapt_cell,
                                   int edge_id,
                                   int& parent_edge_id);

// =====================================================================
// Edge classifier
// =====================================================================

/// Classify a 1D Bernstein polynomial for root structure.
///
/// @param edge_coeffs      Bernstein coefficients (degree+1 entries).
/// @param zero_tol         Tolerance for all-zero.
/// @param sign_tol         Tolerance for all-positive / all-negative.
/// @param max_depth        Maximum subdivision depth for root search.
/// @param[out] green_split_t      Parameter between two distinct root intervals
///                                (valid only if tag == multiple_roots).
/// @param[out] has_green_split_t  True if green_split_t was computed.
/// @return EdgeRootTag.
template <std::floating_point T>
EdgeRootTag classify_edge_roots(std::span<const T> edge_coeffs,
                                T zero_tol, T sign_tol,
                                int max_depth,
                                T& green_split_t,
                                bool& has_green_split_t);

/// Localize the unique root on a 1D Bernstein polynomial already known to
/// have exactly one isolated root.
///
/// Returns false if the polynomial does not represent a unique isolated root.
template <std::floating_point T>
bool locate_one_root_parameter(std::span<const T> edge_coeffs,
                               T zero_tol, T sign_tol,
                               int max_depth,
                               T& root_t);

/// Classify all not-yet-classified edges of an AdaptCell for one level set.
///
/// For each edge with tag == not_classified:
///   1. extract exact edge Bernstein from the LevelSetCell
///   2. classify root structure
///   3. store tag and optional green-split parameter
///
/// @param adapt_cell   The AdaptCell (modified in place).
/// @param ls_cell      LevelSetCell providing Bernstein coefficients.
/// @param level_set_id Which level set (bit position / name index).
/// @param zero_tol     Tolerance for all-zero.
/// @param sign_tol     Tolerance for all-positive / all-negative.
/// @param max_depth    Maximum subdivision depth for root search.
template <std::floating_point T, std::integral I>
void classify_new_edges(AdaptCell<T>& adapt_cell,
                        const LevelSetCell<T, I>& ls_cell,
                        int level_set_id,
                        T zero_tol, T sign_tol,
                        int max_depth);

} // namespace cutcells
