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

namespace cutcells::cell::quad_midpoint
{

struct LocalEdge
{
    int a = -1;
    int b = -1;
};

struct QuadRoleAnalysis
{
    std::array<int, 4> root_local_ids = {-1, -1, -1, -1};
    int num_roots = 0;
    std::array<int, 4> non_root_local_ids = {-1, -1, -1, -1};
    int num_non_roots = 0;
    std::array<int, 4> canonical_to_actual = {0, 1, 2, 3};
    LocalEdge midpoint_edge = {-1, -1};
};

bool is_root_token(int token);

QuadRoleAnalysis analyze_triangle_derived_quadrilateral(
    std::span<const int> quad_tokens);

inline int parent_triangle_edge_token(int token_a, int token_b)
{
    if (is_root_token(token_a) || is_root_token(token_b))
    {
        throw std::runtime_error(
            "parent_triangle_edge_token: midpoint endpoints must be original vertices");
    }

    const int a = token_a - 100;
    const int b = token_b - 100;
    if (a < 0 || a >= 3 || b < 0 || b >= 3)
    {
        throw std::runtime_error(
            "parent_triangle_edge_token: invalid original-vertex token");
    }

    constexpr std::array<LocalEdge, 3> parent_edges = {{
        {0, 1}, {1, 2}, {2, 0},
    }};
    for (std::size_t e = 0; e < parent_edges.size(); ++e)
    {
        const LocalEdge edge = parent_edges[e];
        if ((edge.a == a && edge.b == b) || (edge.a == b && edge.b == a))
            return static_cast<int>(e);
    }

    throw std::runtime_error("parent_triangle_edge_token: parent edge not found");
}

template <std::floating_point T>
struct MidpointInsertionResult
{
    QuadRoleAnalysis analysis;
    std::vector<T> added_vertex_coords;
    std::vector<int> added_vertex_tokens;
    std::vector<LocalEdge> added_vertex_edges;
    std::vector<std::array<int, 3>> triangles;
};

template <std::floating_point T>
MidpointInsertionResult<T> split_triangle_derived_quadrilateral(
    std::span<const T> quad_vertex_coords,
    int coord_dim,
    std::span<const int> quad_tokens,
    int& next_token_base)
{
    if (coord_dim <= 0)
    {
        throw std::invalid_argument(
            "split_triangle_derived_quadrilateral: coord_dim must be positive");
    }
    if (quad_vertex_coords.size() != static_cast<std::size_t>(4 * coord_dim))
    {
        throw std::invalid_argument(
            "split_triangle_derived_quadrilateral: expected 4 quadrilateral vertices in flat storage");
    }
    if (quad_tokens.size() != 4)
    {
        throw std::invalid_argument(
            "split_triangle_derived_quadrilateral: expected 4 quadrilateral tokens");
    }

    MidpointInsertionResult<T> result;
    result.analysis = analyze_triangle_derived_quadrilateral(quad_tokens);
    result.added_vertex_edges.push_back(result.analysis.midpoint_edge);
    const LocalEdge midpoint_edge = result.analysis.midpoint_edge;
    result.added_vertex_tokens.push_back(parent_triangle_edge_token(
        quad_tokens[static_cast<std::size_t>(midpoint_edge.a)],
        quad_tokens[static_cast<std::size_t>(midpoint_edge.b)]));
    (void) next_token_base;

    for (int d = 0; d < coord_dim; ++d)
    {
        const T xa = quad_vertex_coords[static_cast<std::size_t>(
            result.analysis.midpoint_edge.a * coord_dim + d)];
        const T xb = quad_vertex_coords[static_cast<std::size_t>(
            result.analysis.midpoint_edge.b * coord_dim + d)];
        result.added_vertex_coords.push_back(T(0.5) * (xa + xb));
    }

    auto map_canonical_local_id = [&](int canonical_id) -> int
    {
        if (canonical_id < 4)
            return result.analysis.canonical_to_actual[static_cast<std::size_t>(canonical_id)];
        if (canonical_id == 4)
            return 4;
        throw std::runtime_error(
            "split_triangle_derived_quadrilateral: unsupported canonical local id");
    };

    const std::vector<int> all_tokens = {
        quad_tokens[0],
        quad_tokens[1],
        quad_tokens[2],
        quad_tokens[3],
        result.added_vertex_tokens[0],
    };

    const std::array<std::array<int, 3>, 3> canonical_triangles = {{
        {2, 1, 4},
        {4, 1, 0},
        {3, 4, 0},
    }};

    result.triangles.reserve(canonical_triangles.size());
    for (const auto& tri : canonical_triangles)
    {
        const std::array<int, 3> local_ids = {
            map_canonical_local_id(tri[0]),
            map_canonical_local_id(tri[1]),
            map_canonical_local_id(tri[2]),
        };

        result.triangles.push_back({
            all_tokens[static_cast<std::size_t>(local_ids[0])],
            all_tokens[static_cast<std::size_t>(local_ids[1])],
            all_tokens[static_cast<std::size_t>(local_ids[2])],
        });
    }

    return result;
}

} // namespace cutcells::cell::quad_midpoint
