// Copyright (c) 2022 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include "cell_types.h"
#include <array>
#include <span>
#include <stdexcept>

namespace cutcells::cell
{
/// @brief Cell topology data for edge-based clipping algorithms.
/// Edge and face numbering follows VTK conventions.

//-----------------------------------------------------------------------------
// Edge definitions per cell type (VTK ordering)
//-----------------------------------------------------------------------------

/// Interval: 1 edge (0-1)
inline constexpr std::array<std::array<int, 2>, 1> interval_edges = {{
    {0, 1}
}};

/// Triangle: 3 edges
/// Basix vertex order: 0=(0,0), 1=(1,0), 2=(0,1)
/// Basix edge order: 0=(1,2), 1=(0,2), 2=(0,1)
inline constexpr std::array<std::array<int, 2>, 3> triangle_edges = {{
    {1, 2}, {0, 2}, {0, 1}
}};

/// Quadrilateral: 4 edges
/// Basix vertex order: 0=(0,0), 1=(1,0), 2=(0,1), 3=(1,1)
/// Basix edge order: 0=(0,1), 1=(0,2), 2=(1,3), 3=(2,3)
///
///   2 --- 3
///   |     |
///   0 --- 1
inline constexpr std::array<std::array<int, 2>, 4> quadrilateral_edges = {{
    {0, 1}, {0, 2}, {1, 3}, {2, 3}
}};

/// Tetrahedron: 4 triangular faces
/// Basix face order: opposite vertex 0, 1, 2, 3
inline constexpr std::array<std::array<int, 3>, 4> tetrahedron_faces = {{
    {1, 2, 3},  // face 0: opposite vertex 0
    {0, 2, 3},  // face 1: opposite vertex 1
    {0, 1, 3},  // face 2: opposite vertex 2
    {0, 1, 2}   // face 3: opposite vertex 3
}};
inline constexpr std::array<int, 4> tetrahedron_face_sizes = {3, 3, 3, 3};

/// Tetrahedron: 6 edges
/// Basix vertex order: 0=(0,0,0), 1=(1,0,0), 2=(0,1,0), 3=(0,0,1)
/// Basix edge order: 0=(2,3), 1=(1,3), 2=(1,2), 3=(0,3), 4=(0,2), 5=(0,1)
inline constexpr std::array<std::array<int, 2>, 6> tetrahedron_edges = {{
    {2, 3}, {1, 3}, {1, 2}, {0, 3}, {0, 2}, {0, 1}
}};

/// Hexahedron: 12 edges
/// Basix vertex numbering:
///   0=(0,0,0), 1=(1,0,0), 2=(0,1,0), 3=(1,1,0),
///   4=(0,0,1), 5=(1,0,1), 6=(0,1,1), 7=(1,1,1)
///
///        6 -------- 7
///       /|         /|
///      / |        / |
///     4 -------- 5  |
///     |  2 ------|-- 3
///     | /        | /
///     |/         |/
///     0 -------- 1
///
/// Basix edge order:
///   Bottom face: (0,1),(0,2),(1,3),(2,3)
///   Verticals:   (0,4),(1,5),(2,6),(3,7)
///   Top face:    (4,5),(4,6),(5,7),(6,7)
inline constexpr std::array<std::array<int, 2>, 12> hexahedron_edges = {{
    // Bottom face edges (0-3)
    {0, 1}, {0, 2}, {1, 3}, {2, 3},
    // Vertical edges (4-7)
    {0, 4}, {1, 5}, {2, 6}, {3, 7},
    // Top face edges (8-11)
    {4, 5}, {4, 6}, {5, 7}, {6, 7}
}};

/// Prism (Wedge): 9 edges
/// Basix vertex numbering (same as VTK for prism):
///   0=(0,0,0), 1=(1,0,0), 2=(0,1,0), 3=(0,0,1), 4=(1,0,1), 5=(0,1,1)
///
/// Basix edge order:
///   Bottom triangle: (0,1),(0,2),(1,2)
///   Verticals:       (0,3),(1,4),(2,5)
///   Top triangle:    (3,4),(3,5),(4,5)
inline constexpr std::array<std::array<int, 2>, 9> prism_edges = {{
    // Bottom triangle (0-2)
    {0, 1}, {0, 2}, {1, 2},
    // Vertical edges (3-5)
    {0, 3}, {1, 4}, {2, 5},
    // Top triangle (6-8)
    {3, 4}, {3, 5}, {4, 5}
}};

/// Pyramid: 8 edges
/// Basix vertex numbering:
///   0=(0,0,0), 1=(1,0,0), 2=(0,1,0), 3=(1,1,0), 4=(0,0,1) (apex)
///
///         4
///        /|\
///       / | \
///      /  |  \
///     /   |   \
///    2----|----3
///    |    |    |
///    0---------1
///
/// Basix edge order:
///   Base: (0,1),(0,2),(1,3),(2,3)
///   Apex connections: (0,4),(1,4),(3,4),(2,4)
inline constexpr std::array<std::array<int, 2>, 8> pyramid_edges = {{
    // Base edges (0-3)
    {0, 1}, {0, 2}, {1, 3}, {2, 3},
    // Apex edges (4-7)
    {0, 4}, {1, 4}, {3, 4}, {2, 4}
}};

//-----------------------------------------------------------------------------
// Edge count per cell type
//-----------------------------------------------------------------------------

inline constexpr int num_edges(type cell_type)
{
    switch (cell_type)
    {
    case type::point:         return 0;
    case type::interval:      return 1;
    case type::triangle:      return 3;
    case type::quadrilateral: return 4;
    case type::tetrahedron:   return 6;
    case type::hexahedron:    return 12;
    case type::prism:         return 9;
    case type::pyramid:       return 8;
    default:
        throw std::invalid_argument("Unknown cell type in num_edges");
    }
}

//-----------------------------------------------------------------------------
// Edge accessor (returns span to edge array for cell type)
//-----------------------------------------------------------------------------

/// Get edges for a given cell type.
/// @param cell_type The cell type
/// @return Span of edge vertex pairs [v0, v1] for each edge
inline std::span<const std::array<int, 2>> edges(type cell_type)
{
    switch (cell_type)
    {
    case type::interval:
        return std::span(interval_edges);
    case type::triangle:
        return std::span(triangle_edges);
    case type::quadrilateral:
        return std::span(quadrilateral_edges);
    case type::tetrahedron:
        return std::span(tetrahedron_edges);
    case type::hexahedron:
        return std::span(hexahedron_edges);
    case type::prism:
        return std::span(prism_edges);
    case type::pyramid:
        return std::span(pyramid_edges);
    default:
        throw std::invalid_argument("Unknown cell type in edges()");
    }
}

//-----------------------------------------------------------------------------
// Face definitions (for interface stitching, optional future use)
//-----------------------------------------------------------------------------

/// Hexahedron: 6 faces (quads)
/// Face order: -X, +X, -Y, +Y, -Z, +Z (left, right, front, back, bottom, top)
inline constexpr std::array<std::array<int, 4>, 6> hexahedron_faces = {{
    {0, 3, 7, 4},  // -X (left)
    {1, 2, 6, 5},  // +X (right)
    {0, 1, 5, 4},  // -Y (front)
    {3, 2, 6, 7},  // +Y (back)
    {0, 1, 2, 3},  // -Z (bottom)
    {4, 5, 6, 7}   // +Z (top)
}};

/// Prism: 5 faces (2 triangles + 3 quads)
/// Represented as max 4 vertices per face, with -1 padding for triangles
inline constexpr std::array<std::array<int, 4>, 5> prism_faces = {{
    {0, 1, 2, -1},  // bottom triangle
    {3, 4, 5, -1},  // top triangle
    {0, 1, 4, 3},   // quad face
    {1, 2, 5, 4},   // quad face
    {2, 0, 3, 5}    // quad face
}};

/// Pyramid: 5 faces (1 quad base + 4 triangles)
inline constexpr std::array<std::array<int, 4>, 5> pyramid_faces = {{
    {0, 1, 2, 3},   // base quad
    {0, 1, 4, -1},  // triangle
    {1, 2, 4, -1},  // triangle
    {2, 3, 4, -1},  // triangle
    {3, 0, 4, -1}   // triangle
}};

//-----------------------------------------------------------------------------
// Face vertex count per face (for mixed cells)
//-----------------------------------------------------------------------------

inline constexpr std::array<int, 6> hexahedron_face_sizes = {4, 4, 4, 4, 4, 4};
inline constexpr std::array<int, 5> prism_face_sizes = {3, 3, 4, 4, 4};
inline constexpr std::array<int, 5> pyramid_face_sizes = {4, 3, 3, 3, 3};

/// Get number of faces for a cell type
inline constexpr int num_faces(type cell_type)
{
    switch (cell_type)
    {
    case type::point:         return 0;
    case type::interval:      return 0;
    case type::triangle:      return 1;
    case type::quadrilateral: return 1;
    case type::tetrahedron:   return 4;
    case type::hexahedron:    return 6;
    case type::prism:         return 5;
    case type::pyramid:       return 5;
    default:
        throw std::invalid_argument("Unknown cell type in num_faces");
    }
}

/// Get the number of vertices per face for a cell type.
inline std::span<const int> face_sizes(type cell_type)
{
    switch (cell_type)
    {
    case type::tetrahedron:   return std::span(tetrahedron_face_sizes);
    case type::hexahedron:    return std::span(hexahedron_face_sizes);
    case type::prism:         return std::span(prism_face_sizes);
    case type::pyramid:       return std::span(pyramid_face_sizes);
    default:
        throw std::invalid_argument("Unknown cell type in face_sizes");
    }
}

/// Get the vertex indices for face `face_id` of a tetrahedron.
/// Returns a span of 3 vertex indices.
inline std::span<const int, 3> tet_face(int face_id)
{
    return std::span<const int, 3>(tetrahedron_faces[static_cast<std::size_t>(face_id)]);
}

/// Get the cell type of a specific face of a cell.
inline type face_type(type cell_type, int face_id)
{
    switch (cell_type)
    {
    case type::tetrahedron:
        return type::triangle;
    case type::hexahedron:
        return type::quadrilateral;
    case type::prism:
        return (face_id < 2) ? type::triangle : type::quadrilateral;
    case type::pyramid:
        return (face_id == 0) ? type::quadrilateral : type::triangle;
    default:
        throw std::invalid_argument("Unknown cell type in face_type");
    }
}

/// Get the vertex indices for face `face_id` of a 3D cell.
/// The returned span has `face_sizes(cell_type)[face_id]` entries.
inline std::span<const int> face_vertices(type cell_type, int face_id)
{
    const int fsize = face_sizes(cell_type)[static_cast<std::size_t>(face_id)];
    switch (cell_type)
    {
    case type::tetrahedron:
        return std::span<const int>(
            tetrahedron_faces[static_cast<std::size_t>(face_id)].data(),
            static_cast<std::size_t>(fsize));
    case type::hexahedron:
        return std::span<const int>(
            hexahedron_faces[static_cast<std::size_t>(face_id)].data(),
            static_cast<std::size_t>(fsize));
    case type::prism:
        return std::span<const int>(
            prism_faces[static_cast<std::size_t>(face_id)].data(),
            static_cast<std::size_t>(fsize));
    case type::pyramid:
        return std::span<const int>(
            pyramid_faces[static_cast<std::size_t>(face_id)].data(),
            static_cast<std::size_t>(fsize));
    default:
        throw std::invalid_argument("Unknown cell type in face_vertices");
    }
}

} // namespace cutcells::cell
