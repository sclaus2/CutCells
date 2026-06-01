// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <array>
#include <concepts>
#include <span>
#include <stdexcept>
#include <vector>

namespace cutcells::cell::prism_midpoint
{

/// Prism family induced by cutting a tetrahedron.
///
/// roots3:
///   The prism carries three root vertices and three original tetra vertices.
///   The quasi-normal midpoint strategy inserts one midpoint on each prism edge
///   connecting two original vertices.
///
/// roots4:
///   The prism carries four root vertices and two original tetra vertices.
///   The quasi-normal midpoint strategy inserts one midpoint on the single prism
///   edge connecting the two original vertices.
enum class family
{
    roots3,
    roots4,
};

/// Local prism edge with local vertex ids in the prism ordering used by the
/// tetrahedron LUT subcell definition.
struct LocalEdge
{
    int a = -1;
    int b = -1;
};

/// Static classification of a tetra-derived prism.
struct PrismRoleAnalysis
{
    family split_family = family::roots3;
    std::array<int, 6> root_local_ids = {-1, -1, -1, -1, -1, -1};
    int num_roots = 0;
    std::array<int, 6> non_root_local_ids = {-1, -1, -1, -1, -1, -1};
    int num_non_roots = 0;
    std::vector<LocalEdge> midpoint_edges;
};

/// True for tetra-cut root tokens. In the current tetra LUT these are the
/// intersection-point edge ids 0..5. Original tetra vertices use 100+vid.
bool is_root_token(int token);

/// Prism edge list in the local ordering expected by the tetrahedron cutter and
/// by cutcells::cell::triangulation(type::prism, ...).
std::span<const LocalEdge> prism_edges();

/// Classify a tetra-derived prism purely from its 6 tokens.
family classify_family(std::span<const int> prism_tokens);

/// Return root/non-root local ids and the default uncut-edge midpoint placement
/// suggested by the tetra-prism analysis.
PrismRoleAnalysis analyze_tetra_derived_prism(std::span<const int> prism_tokens);

inline int parent_tetrahedron_edge_token(int token_a, int token_b)
{
    if (is_root_token(token_a) || is_root_token(token_b))
    {
        throw std::runtime_error(
            "parent_tetrahedron_edge_token: midpoint endpoints must be original vertices");
    }

    const int a = token_a - 100;
    const int b = token_b - 100;
    if (a < 0 || a >= 4 || b < 0 || b >= 4)
    {
        throw std::runtime_error(
            "parent_tetrahedron_edge_token: invalid original-vertex token");
    }

    constexpr std::array<LocalEdge, 6> parent_edges = {{
        {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3},
    }};
    for (std::size_t e = 0; e < parent_edges.size(); ++e)
    {
        const LocalEdge edge = parent_edges[e];
        if ((edge.a == a && edge.b == b) || (edge.a == b && edge.b == a))
            return static_cast<int>(e);
    }

    throw std::runtime_error("parent_tetrahedron_edge_token: parent edge not found");
}

template <std::floating_point T>
struct MidpointInsertionResult
{
    PrismRoleAnalysis analysis;
    std::vector<T> added_vertex_coords;   // flat coord_dim storage
    std::vector<int> added_vertex_tokens; // synthetic tokens chosen by caller
    std::vector<LocalEdge> added_vertex_edges;
    std::vector<std::array<int, 4>> tets; // token-based tetrahedra
};

/// Insert the default midpoint set for a tetra-derived prism.
///
/// The returned midpoints are not yet wired into any tetra connectivity. This
/// helper exists so the new prism strategy can be implemented in a standalone
/// module first, then integrated at the prism-to-tet call sites.
///
/// @param prism_vertex_coords flat 6*coord_dim local coordinate storage
/// @param coord_dim prism coordinate dimension, expected 3 for tetra cuts
/// @param prism_tokens 6 tokens from the tetra LUT prism subcell
/// @param next_token_base caller-owned synthetic-token counter; incremented for
///        each created midpoint
template <std::floating_point T>
MidpointInsertionResult<T> create_default_midpoints(
    std::span<const T> prism_vertex_coords,
    int coord_dim,
    std::span<const int> prism_tokens,
    int& next_token_base)
{
    if (coord_dim <= 0)
        throw std::invalid_argument("create_default_midpoints: coord_dim must be positive");
    if (prism_vertex_coords.size() != static_cast<std::size_t>(6 * coord_dim))
    {
        throw std::invalid_argument(
            "create_default_midpoints: expected 6 prism vertices in flat storage");
    }
    if (prism_tokens.size() != 6)
        throw std::invalid_argument("create_default_midpoints: expected 6 prism tokens");

    MidpointInsertionResult<T> result;
    result.analysis = analyze_tetra_derived_prism(prism_tokens);
    result.added_vertex_coords.reserve(
        static_cast<std::size_t>(result.analysis.midpoint_edges.size() * coord_dim));
    result.added_vertex_tokens.reserve(result.analysis.midpoint_edges.size());
    result.added_vertex_edges.reserve(result.analysis.midpoint_edges.size());

    for (const LocalEdge edge : result.analysis.midpoint_edges)
    {
        result.added_vertex_tokens.push_back(parent_tetrahedron_edge_token(
            prism_tokens[static_cast<std::size_t>(edge.a)],
            prism_tokens[static_cast<std::size_t>(edge.b)]));
        result.added_vertex_edges.push_back(edge);
        for (int d = 0; d < coord_dim; ++d)
        {
            const T xa = prism_vertex_coords[static_cast<std::size_t>(edge.a * coord_dim + d)];
            const T xb = prism_vertex_coords[static_cast<std::size_t>(edge.b * coord_dim + d)];
            result.added_vertex_coords.push_back(T(0.5) * (xa + xb));
        }
    }
    (void) next_token_base;

    return result;
}

template <std::floating_point T>
T signed_tet_jacobian(std::span<const T> vertex_coords,
                      int coord_dim,
                      const std::array<int, 4>& tet_local_ids);

template <std::floating_point T>
void orient_tet_positive(std::span<const T> vertex_coords,
                         int coord_dim,
                         std::array<int, 4>& tet_local_ids);

/// Build the standalone tetrahedral split for a tetra-derived prism using the
/// default midpoint insertion pattern for the detected family.
///
/// The returned tetrahedra are expressed in tokens: original prism tokens are
/// reused and new midpoint vertices receive synthetic tokens chosen by the
/// caller. This lets the tetrahedron cutter integrate the result by only
/// appending the new coordinates and extending its token->local-index map.
template <std::floating_point T>
MidpointInsertionResult<T> split_tetra_derived_prism(
    std::span<const T> prism_vertex_coords,
    int coord_dim,
    std::span<const int> prism_tokens,
    int& next_token_base)
{
    constexpr std::array<std::array<int, 2>, 6> tetra_edges = {{
        {{0, 1}},
        {{1, 2}},
        {{2, 0}},
        {{0, 3}},
        {{1, 3}},
        {{2, 3}},
    }};

    MidpointInsertionResult<T> result = create_default_midpoints(
        prism_vertex_coords, coord_dim, prism_tokens, next_token_base);

    std::array<int, 6> canonical_to_actual = {0, 1, 2, 3, 4, 5};

    auto root_local_for_original_local = [&](int non_root_local_id) -> int
    {
        const int original_token = prism_tokens[static_cast<std::size_t>(non_root_local_id)];
        if (is_root_token(original_token))
        {
            throw std::runtime_error(
                "split_tetra_derived_prism: expected original tetra vertex token for non-root prism vertex");
        }

        const int original_vertex = original_token - 100;
        for (int i = 0; i < result.analysis.num_roots; ++i)
        {
            const int root_local_id = result.analysis.root_local_ids[static_cast<std::size_t>(i)];
            const int root_token = prism_tokens[static_cast<std::size_t>(root_local_id)];
            const auto& edge = tetra_edges[static_cast<std::size_t>(root_token)];
            if (edge[0] == original_vertex || edge[1] == original_vertex)
                return root_local_id;
        }

        throw std::runtime_error(
            "split_tetra_derived_prism: failed to match original prism vertex to a root edge");
    };

    if (result.analysis.split_family == family::roots3)
    {
        const bool bottom_non_roots
            = result.analysis.num_non_roots == 3
           && result.analysis.non_root_local_ids[0] == 0
           && result.analysis.non_root_local_ids[1] == 1
           && result.analysis.non_root_local_ids[2] == 2;
        const bool top_non_roots
            = result.analysis.num_non_roots == 3
           && result.analysis.non_root_local_ids[0] == 3
           && result.analysis.non_root_local_ids[1] == 4
           && result.analysis.non_root_local_ids[2] == 5;

        if (bottom_non_roots)
        {
            canonical_to_actual = {
                0,
                1,
                2,
                root_local_for_original_local(0),
                root_local_for_original_local(1),
                root_local_for_original_local(2),
            };
        }
        else if (top_non_roots)
        {
            canonical_to_actual = {
                3,
                4,
                5,
                root_local_for_original_local(3),
                root_local_for_original_local(4),
                root_local_for_original_local(5),
            };
        }
        else
        {
            throw std::runtime_error(
                "split_tetra_derived_prism: roots3 family expected one full prism triangle of non-root vertices");
        }
    }
    else
    {
        if (result.analysis.num_non_roots != 2)
        {
            throw std::runtime_error(
                "split_tetra_derived_prism: roots4 family expected exactly two non-root vertices");
        }

        const int a = result.analysis.non_root_local_ids[0];
        const int b = result.analysis.non_root_local_ids[1];

        int k = -1;
        if ((a == 0 && b == 3) || (a == 3 && b == 0))
            k = 0;
        else if ((a == 1 && b == 4) || (a == 4 && b == 1))
            k = 1;
        else if ((a == 2 && b == 5) || (a == 5 && b == 2))
            k = 2;

        if (k < 0)
        {
            throw std::runtime_error(
                "split_tetra_derived_prism: roots4 family expected non-root vertices on one prism vertical edge");
        }

        canonical_to_actual = {
            k,
            (k + 1) % 3,
            (k + 2) % 3,
            k + 3,
            ((k + 1) % 3) + 3,
            ((k + 2) % 3) + 3,
        };
    }

    std::vector<int> all_tokens(6);
    for (int i = 0; i < 6; ++i)
        all_tokens[static_cast<std::size_t>(i)] = prism_tokens[static_cast<std::size_t>(i)];
    for (int token : result.added_vertex_tokens)
        all_tokens.push_back(token);

    std::vector<T> all_coords(prism_vertex_coords.begin(), prism_vertex_coords.end());
    all_coords.insert(
        all_coords.end(),
        result.added_vertex_coords.begin(),
        result.added_vertex_coords.end());

    auto midpoint_local_id_for_actual_edge = [&](int actual_a, int actual_b) -> int
    {
        for (std::size_t i = 0; i < result.added_vertex_edges.size(); ++i)
        {
            const LocalEdge edge = result.added_vertex_edges[i];
            const bool matches
                = (edge.a == actual_a && edge.b == actual_b)
               || (edge.a == actual_b && edge.b == actual_a);
            if (matches)
                return 6 + static_cast<int>(i);
        }

        throw std::runtime_error(
            "split_tetra_derived_prism: failed to match canonical midpoint edge");
    };

    auto map_canonical_local_id = [&](int canonical_id) -> int
    {
        if (canonical_id < 6)
            return canonical_to_actual[static_cast<std::size_t>(canonical_id)];

        if (result.analysis.split_family == family::roots3)
        {
            switch (canonical_id)
            {
                case 6:
                    return midpoint_local_id_for_actual_edge(
                        canonical_to_actual[0], canonical_to_actual[1]);
                case 7:
                    return midpoint_local_id_for_actual_edge(
                        canonical_to_actual[1], canonical_to_actual[2]);
                case 8:
                    return midpoint_local_id_for_actual_edge(
                        canonical_to_actual[2], canonical_to_actual[0]);
                default:
                    break;
            }
        }
        else
        {
            if (canonical_id == 6)
            {
                return midpoint_local_id_for_actual_edge(
                    canonical_to_actual[0], canonical_to_actual[3]);
            }
        }

        throw std::runtime_error("split_tetra_derived_prism: unsupported canonical local id");
    };

    std::vector<std::array<int, 4>> canonical_tets;
    if (result.analysis.split_family == family::roots3)
    {
        // Build the three original-vertex corner tetrahedra first.
        //
        // The remaining central cell has vertices {3,4,5,6,7,8}:
        // - {3,4,5} are the root vertices (one on each tetra edge)
        // - {6,7,8} are midpoints on the non-root triangle edges
        //
        // This remaining polyhedron is a triangular antiprism with boundary
        // triangles:
        // - top:    {3,4,5}
        // - bottom: {6,7,8}
        // - sides:  {3,4,6}, {4,5,7}, {5,3,8}, {3,6,8}, {4,7,6}, {5,8,7}
        //
        // A 4-tet decomposition (total 7 tets incl. corners) is required here:
        // a 3-tet "prism" split does not cover the antiprism volume.
        canonical_tets = {
            std::array<int, 4>{0, 6, 8, 3},
            std::array<int, 4>{1, 7, 6, 4},
            std::array<int, 4>{2, 8, 7, 5},
            // Central triangular antiprism {3,4,5,6,7,8} split into 4 tets.
            // This corresponds to inserting the internal diagonal (3,7).
            std::array<int, 4>{3, 4, 6, 7},
            std::array<int, 4>{3, 5, 8, 7},
            std::array<int, 4>{3, 4, 5, 7},
            std::array<int, 4>{3, 6, 8, 7},
        };
    }
    else
    {
        // Delaunay-style split on the canonical prism with one midpoint on the
        // non-root edge (0,3).
        canonical_tets = {
            std::array<int, 4>{6, 2, 1, 0},
            std::array<int, 4>{6, 5, 4, 3},
            std::array<int, 4>{6, 4, 2, 1},
            std::array<int, 4>{6, 5, 4, 2},
        };
    }

    result.tets.clear();
    result.tets.reserve(canonical_tets.size());

    for (std::array<int, 4> tet : canonical_tets)
    {
        std::array<int, 4> local_ids = {
            map_canonical_local_id(tet[0]),
            map_canonical_local_id(tet[1]),
            map_canonical_local_id(tet[2]),
            map_canonical_local_id(tet[3]),
        };
        orient_tet_positive(
            std::span<const T>(all_coords.data(), all_coords.size()),
            coord_dim,
            local_ids);

        std::array<int, 4> token_tet = {
            all_tokens[static_cast<std::size_t>(local_ids[0])],
            all_tokens[static_cast<std::size_t>(local_ids[1])],
            all_tokens[static_cast<std::size_t>(local_ids[2])],
            all_tokens[static_cast<std::size_t>(local_ids[3])],
        };
        result.tets.push_back(token_tet);
    }

    return result;
}

/// Signed tetra determinant in local coordinates. Positive means the tet
/// orientation is consistent with the ambient local coordinate frame.
template <std::floating_point T>
T signed_tet_jacobian(std::span<const T> vertex_coords,
                      int coord_dim,
                      const std::array<int, 4>& tet_local_ids)
{
    if (coord_dim != 3)
        throw std::invalid_argument("signed_tet_jacobian: coord_dim must be 3");

    auto coord = [&](int local_id, int d) -> T
    {
        return vertex_coords[static_cast<std::size_t>(local_id * coord_dim + d)];
    };

    const T ax = coord(tet_local_ids[1], 0) - coord(tet_local_ids[0], 0);
    const T ay = coord(tet_local_ids[1], 1) - coord(tet_local_ids[0], 1);
    const T az = coord(tet_local_ids[1], 2) - coord(tet_local_ids[0], 2);

    const T bx = coord(tet_local_ids[2], 0) - coord(tet_local_ids[0], 0);
    const T by = coord(tet_local_ids[2], 1) - coord(tet_local_ids[0], 1);
    const T bz = coord(tet_local_ids[2], 2) - coord(tet_local_ids[0], 2);

    const T cx = coord(tet_local_ids[3], 0) - coord(tet_local_ids[0], 0);
    const T cy = coord(tet_local_ids[3], 1) - coord(tet_local_ids[0], 1);
    const T cz = coord(tet_local_ids[3], 2) - coord(tet_local_ids[0], 2);

    return ax * (by * cz - bz * cy)
         - ay * (bx * cz - bz * cx)
         + az * (bx * cy - by * cx);
}

/// Flip the last two vertices if needed so the tetra has non-negative
/// orientation in the local coordinate frame.
template <std::floating_point T>
void orient_tet_positive(std::span<const T> vertex_coords,
                         int coord_dim,
                         std::array<int, 4>& tet_local_ids)
{
    if (signed_tet_jacobian(vertex_coords, coord_dim, tet_local_ids) < T(0))
        std::swap(tet_local_ids[2], tet_local_ids[3]);
}

} // namespace cutcells::cell::prism_midpoint
