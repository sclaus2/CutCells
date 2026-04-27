// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

#include "adapt_cell.h"
#include "level_set_cell.h"

namespace cutcells::curving
{

enum class NodeFamily : std::uint8_t
{
    gll = 0,
    equispaced = 1,
    lagrange = 2
};

enum class CurvingStatus : std::uint8_t
{
    not_built = 0,
    in_progress = 1,
    curved = 2,
    failed = 3
};

enum class CurvingFailureCode : std::uint8_t
{
    none = 0,
    exact_vertex = 1,
    boundary_from_edge = 2,
    invalid_constraint_count = 3,
    missing_level_set_cell = 4,
    empty_zero_mask = 5,
    unsupported_entity = 6,
    no_host_interval = 7,
    no_sign_changing_bracket = 8,
    brent_failed = 9,
    outside_host_domain = 10,
    singular_gradient_system = 11,
    line_search_failed = 12,
    max_iterations = 13,
    missing_boundary_edge = 14,
    boundary_edge_failed = 15,
    projection_failed = 16,
    closest_face_retry_failed = 17,
    constrained_newton_failed = 18
};

enum class CurvingProjectionMode : std::uint8_t
{
    none = 0,
    safe_line = 1,
    closest_face_retry = 2,
    constrained_newton = 3,
    vector_newton = 4
};

template <std::floating_point T>
struct CurvingOptions
{
    int geometry_order = 2;
    NodeFamily node_family = NodeFamily::gll;
    int max_iter = 32;
    T xtol = T(1e-12);
    T ftol = T(1e-12);
    T domain_tol = T(1e-10);
    T active_face_tol = T(1e-9);
    int max_subdivision_depth = 3;
};

template <std::floating_point T>
struct CurvedZeroEntityState
{
    CurvingStatus status = CurvingStatus::not_built;
    int geometry_order = -1;
    NodeFamily node_family = NodeFamily::gll;
    std::uint32_t zero_entity_version = 0;
    std::uint64_t zero_mask = 0;
    std::string failure_reason;

    // Full curved interpolation nodes, including endpoints/boundary nodes.
    // Coordinates are in the parent background cell reference frame.
    std::vector<T> ref_nodes;

    // Per intended interpolation node, including nodes that failed before the
    // entity as a whole was accepted. Successful vertices/boundary nodes have
    // zero iterations and a non-error failure code documenting their origin.
    std::vector<std::int32_t> node_iterations;
    std::vector<std::uint8_t> node_status;
    std::vector<std::uint8_t> node_failure_code;
    std::vector<T> node_residual;
    std::vector<std::uint32_t> node_active_face_mask;
    std::vector<std::int32_t> node_closest_face_id;
    std::vector<std::int32_t> node_safe_subspace_dim;
    std::vector<std::uint8_t> node_projection_mode;
    std::vector<std::int32_t> node_retry_count;
};

struct CurvingIdentity
{
    int cut_cell_id = -1;
    int local_zero_entity_id = -1;
    int dim = -1;
    std::int8_t parent_dim = -1;
    std::int32_t parent_id = -1;
    std::uint64_t zero_mask = 0;
};

template <std::floating_point T, std::integral I = int>
struct CurvingData
{
    std::vector<CurvingIdentity> identities;
    std::vector<CurvedZeroEntityState<T>> states;

    // CSR-like local lookup: local_to_canonical_offsets[k]..[k+1] belongs to
    // cut cell k and stores canonical ids for local zero entity ids.
    std::vector<int> local_to_canonical_offsets;
    std::vector<int> local_to_canonical;

    int num_cut_cells = 0;
    bool identity_valid = false;

    void clear()
    {
        identities.clear();
        states.clear();
        local_to_canonical_offsets.clear();
        local_to_canonical.clear();
        num_cut_cells = 0;
        identity_valid = false;
    }
};

NodeFamily node_family_from_string(std::string_view name);
std::string_view node_family_name(NodeFamily family);

template <std::floating_point T, std::integral I>
void rebuild_identity(CurvingData<T, I>& curving,
                      std::span<const I> parent_cell_ids,
                      std::span<const AdaptCell<T>> adapt_cells);

template <std::floating_point T, std::integral I>
const CurvedZeroEntityState<T>& ensure_curved(
    CurvingData<T, I>& curving,
    std::span<const I> parent_cell_ids,
    std::span<const AdaptCell<T>> adapt_cells,
    std::span<const LevelSetCell<T, I>> level_set_cells,
    std::span<const int> ls_offsets,
    int cut_cell_id,
    int local_zero_entity_id,
    const CurvingOptions<T>& options);

template <std::floating_point T, std::integral I>
void ensure_all_curved(CurvingData<T, I>& curving,
                       std::span<const I> parent_cell_ids,
                       std::span<const AdaptCell<T>> adapt_cells,
                       std::span<const LevelSetCell<T, I>> level_set_cells,
                       std::span<const int> ls_offsets,
                       const CurvingOptions<T>& options);

} // namespace cutcells::curving
