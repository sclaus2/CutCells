// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include "cell_types.h"

#include <array>
#include <concepts>
#include <span>
#include <stdexcept>
#include <vector>

namespace cutcells::cell
{

template <std::floating_point T>
inline std::vector<T> reference_vertices(type cell_type)
{
    switch (cell_type)
    {
        case type::interval:
            return {T(0), T(1)};
        case type::triangle:
            return {T(0), T(0), T(1), T(0), T(0), T(1)};
        case type::quadrilateral:
            return {T(0), T(0), T(1), T(0), T(0), T(1), T(1), T(1)};
        case type::tetrahedron:
            return {T(0), T(0), T(0),
                    T(1), T(0), T(0),
                    T(0), T(1), T(0),
                    T(0), T(0), T(1)};
        case type::hexahedron:
            return {T(0), T(0), T(0),
                    T(1), T(0), T(0),
                    T(0), T(1), T(0),
                    T(1), T(1), T(0),
                    T(0), T(0), T(1),
                    T(1), T(0), T(1),
                    T(0), T(1), T(1),
                    T(1), T(1), T(1)};
        case type::prism:
            return {T(0), T(0), T(0),
                    T(1), T(0), T(0),
                    T(0), T(1), T(0),
                    T(0), T(0), T(1),
                    T(1), T(0), T(1),
                    T(0), T(1), T(1)};
        case type::pyramid:
            return {T(0), T(0), T(0),
                    T(1), T(0), T(0),
                    T(0), T(1), T(0),
                    T(1), T(1), T(0),
                    T(0), T(0), T(1)};
        default:
            throw std::invalid_argument("reference_vertices: unsupported cell type");
    }
}

inline std::span<const int> vtk_to_basix_vertex_permutation(type cell_type)
{
    static constexpr std::array<int, 2> interval = {0, 1};
    static constexpr std::array<int, 3> triangle = {0, 1, 2};
    static constexpr std::array<int, 4> quadrilateral = {0, 1, 3, 2};
    static constexpr std::array<int, 4> tetrahedron = {0, 1, 2, 3};
    static constexpr std::array<int, 8> hexahedron = {0, 1, 3, 2, 4, 5, 7, 6};
    static constexpr std::array<int, 6> prism = {0, 1, 2, 3, 4, 5};
    static constexpr std::array<int, 5> pyramid = {0, 1, 3, 2, 4};

    switch (cell_type)
    {
        case type::interval: return std::span(interval);
        case type::triangle: return std::span(triangle);
        case type::quadrilateral: return std::span(quadrilateral);
        case type::tetrahedron: return std::span(tetrahedron);
        case type::hexahedron: return std::span(hexahedron);
        case type::prism: return std::span(prism);
        case type::pyramid: return std::span(pyramid);
        default:
            throw std::invalid_argument("vtk_to_basix_vertex_permutation: unsupported cell type");
    }
}

inline std::span<const int> basix_to_vtk_vertex_permutation(type cell_type)
{
    // For the supported first-order cells the inverse happens to equal the forward map.
    return vtk_to_basix_vertex_permutation(cell_type);
}

inline int vtk_to_basix_vertex(type cell_type, int vtk_vertex)
{
    const auto perm = vtk_to_basix_vertex_permutation(cell_type);
    if (vtk_vertex < 0 || vtk_vertex >= static_cast<int>(perm.size()))
        throw std::invalid_argument("vtk_to_basix_vertex: vertex id out of bounds");

    for (std::size_t i = 0; i < perm.size(); ++i)
    {
        if (perm[i] == vtk_vertex)
            return static_cast<int>(i);
    }

    throw std::invalid_argument("vtk_to_basix_vertex: invalid vtk vertex id");
}

inline int basix_to_vtk_vertex(type cell_type, int basix_vertex)
{
    const auto perm = basix_to_vtk_vertex_permutation(cell_type);
    if (basix_vertex < 0 || basix_vertex >= static_cast<int>(perm.size()))
        throw std::invalid_argument("basix_to_vtk_vertex: vertex id out of bounds");
    return perm[static_cast<std::size_t>(basix_vertex)];
}

inline std::span<const int> vtk_to_basix_edge_permutation(type cell_type)
{
    static constexpr std::array<int, 1> interval = {0};
    static constexpr std::array<int, 3> triangle = {2, 0, 1};
    static constexpr std::array<int, 4> quadrilateral = {0, 2, 3, 1};
    static constexpr std::array<int, 6> tetrahedron = {5, 2, 4, 3, 1, 0};
    static constexpr std::array<int, 12> hexahedron = {0, 3, 1, 2, 4, 7, 5, 6, 8, 9, 11, 10};
    static constexpr std::array<int, 9> prism = {0, 2, 1, 6, 8, 7, 3, 4, 5};
    static constexpr std::array<int, 8> pyramid = {0, 2, 3, 1, 4, 5, 6, 7};

    switch (cell_type)
    {
        case type::interval: return std::span(interval);
        case type::triangle: return std::span(triangle);
        case type::quadrilateral: return std::span(quadrilateral);
        case type::tetrahedron: return std::span(tetrahedron);
        case type::hexahedron: return std::span(hexahedron);
        case type::prism: return std::span(prism);
        case type::pyramid: return std::span(pyramid);
        default:
            throw std::invalid_argument("vtk_to_basix_edge_permutation: unsupported cell type");
    }
}

inline int vtk_to_basix_edge(type cell_type, int vtk_edge)
{
    const auto perm = vtk_to_basix_edge_permutation(cell_type);
    if (vtk_edge < 0 || vtk_edge >= static_cast<int>(perm.size()))
        throw std::invalid_argument("vtk_to_basix_edge: edge id out of bounds");
    return perm[static_cast<std::size_t>(vtk_edge)];
}

inline int basix_to_vtk_edge(type cell_type, int basix_edge)
{
    const auto perm = vtk_to_basix_edge_permutation(cell_type);
    if (basix_edge < 0 || basix_edge >= static_cast<int>(perm.size()))
        throw std::invalid_argument("basix_to_vtk_edge: edge id out of bounds");

    for (std::size_t i = 0; i < perm.size(); ++i)
    {
        if (perm[i] == basix_edge)
            return static_cast<int>(i);
    }

    throw std::invalid_argument("basix_to_vtk_edge: invalid basix edge id");
}

inline int vtk_parent_entity_token_to_basix(type cell_type, int token)
{
    if (token < 0)
        return token;
    if (token < 100)
        return vtk_to_basix_edge(cell_type, token);
    if (token < 200)
        return 100 + vtk_to_basix_vertex(cell_type, token - 100);
    return token;
}

template <std::floating_point T>
inline std::vector<T> permute_vertex_data(std::span<const T> data, int ncomp,
                                          std::span<const int> permutation)
{
    if (ncomp < 0)
        throw std::invalid_argument("permute_vertex_data: negative component count");
    if (data.size() != static_cast<std::size_t>(ncomp) * permutation.size())
        throw std::invalid_argument("permute_vertex_data: data size does not match permutation");

    std::vector<T> out(data.size());
    for (std::size_t i = 0; i < permutation.size(); ++i)
    {
        const int src = permutation[i];
        if (src < 0 || src >= static_cast<int>(permutation.size()))
            throw std::invalid_argument("permute_vertex_data: permutation entry out of bounds");

        for (int c = 0; c < ncomp; ++c)
        {
            out[i * static_cast<std::size_t>(ncomp) + static_cast<std::size_t>(c)]
                = data[static_cast<std::size_t>(src) * static_cast<std::size_t>(ncomp)
                     + static_cast<std::size_t>(c)];
        }
    }
    return out;
}

template <std::integral I>
inline std::vector<I> permute_vertex_ids(std::span<const I> ids,
                                         std::span<const int> permutation)
{
    if (ids.size() != permutation.size())
        throw std::invalid_argument("permute_vertex_ids: size does not match permutation");

    std::vector<I> out(ids.size());
    for (std::size_t i = 0; i < permutation.size(); ++i)
    {
        const int src = permutation[i];
        if (src < 0 || src >= static_cast<int>(permutation.size()))
            throw std::invalid_argument("permute_vertex_ids: permutation entry out of bounds");
        out[i] = ids[static_cast<std::size_t>(src)];
    }
    return out;
}

template <std::integral I, std::size_t N>
inline void reorder_subcell_vertices_from_vtk_to_basix(type cell_type,
                                                       std::array<I, N>& vertices,
                                                       int nvertices)
{
    const auto perm = vtk_to_basix_vertex_permutation(cell_type);
    if (nvertices != static_cast<int>(perm.size()))
        return;

    const auto reordered = permute_vertex_ids<I>(
        std::span<const I>(vertices.data(), static_cast<std::size_t>(nvertices)), perm);
    for (int i = 0; i < nvertices; ++i)
        vertices[static_cast<std::size_t>(i)] = reordered[static_cast<std::size_t>(i)];
}

template <std::size_t N>
inline std::array<int, N> remap_token_to_vertex_map_from_vtk_to_basix(
    type cell_type, const std::array<int, N>& vtk_token_to_vertex,
    int n_edges, int n_vertices, int n_special = 0)
{
    std::array<int, N> basix_token_to_vertex;
    basix_token_to_vertex.fill(-1);

    for (int token = 0; token < n_edges; ++token)
    {
        if (vtk_token_to_vertex[static_cast<std::size_t>(token)] >= 0)
        {
            basix_token_to_vertex[static_cast<std::size_t>(
                vtk_parent_entity_token_to_basix(cell_type, token))]
                = vtk_token_to_vertex[static_cast<std::size_t>(token)];
        }
    }

    for (int token = 100; token < 100 + n_vertices; ++token)
    {
        if (vtk_token_to_vertex[static_cast<std::size_t>(token)] >= 0)
        {
            basix_token_to_vertex[static_cast<std::size_t>(
                vtk_parent_entity_token_to_basix(cell_type, token))]
                = vtk_token_to_vertex[static_cast<std::size_t>(token)];
        }
    }

    for (int token = 200; token < 200 + n_special; ++token)
    {
        if (vtk_token_to_vertex[static_cast<std::size_t>(token)] >= 0)
            basix_token_to_vertex[static_cast<std::size_t>(token)]
                = vtk_token_to_vertex[static_cast<std::size_t>(token)];
    }

    return basix_token_to_vertex;
}

} // namespace cutcells::cell
