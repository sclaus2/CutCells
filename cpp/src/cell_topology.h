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
/// Edge and face numbering follows Basix conventions.

//-----------------------------------------------------------------------------
// Edge definitions per cell type (Basix ordering)
//-----------------------------------------------------------------------------

/// Interval: 1 edge (0-1)
inline constexpr std::array<std::array<int, 2>, 1> interval_edges = {{
    {0, 1}
}};

/// Triangle: 3 edges
/// Basix edge order: (1,2), (0,2), (0,1)
inline constexpr std::array<std::array<int, 2>, 3> triangle_edges = {{
    {1, 2}, {0, 2}, {0, 1}
}};

/// Quadrilateral: 4 edges
/// Basix edge order: (0,1), (0,2), (1,3), (2,3)
///
///   3 --- 2
///   |     |
///   0 --- 1
inline constexpr std::array<std::array<int, 2>, 4> quadrilateral_edges = {{
    {0, 1}, {0, 2}, {1, 3}, {2, 3}
}};

/// Tetrahedron: 6 edges
/// Basix edge order: (2,3), (1,3), (1,2), (0,3), (0,2), (0,1)
inline constexpr std::array<std::array<int, 2>, 6> tetrahedron_edges = {{
    {2, 3}, {1, 3}, {1, 2}, {0, 3}, {0, 2}, {0, 1}
}};

/// Hexahedron: 12 edges
/// Basix edge order:
///   (0,1), (0,2), (0,4), (1,3), (1,5), (2,3),
///   (2,6), (3,7), (4,5), (4,6), (5,7), (6,7)
inline constexpr std::array<std::array<int, 2>, 12> hexahedron_edges = {{
    {0, 1}, {0, 2}, {0, 4}, {1, 3},
    {1, 5}, {2, 3}, {2, 6}, {3, 7},
    {4, 5}, {4, 6}, {5, 7}, {6, 7}
}};

/// Prism (Wedge): 9 edges
/// Basix edge order:
///   (0,1), (0,2), (0,3), (1,2), (1,4),
///   (2,5), (3,4), (3,5), (4,5)
inline constexpr std::array<std::array<int, 2>, 9> prism_edges = {{
    {0, 1}, {0, 2}, {0, 3},
    {1, 2}, {1, 4}, {2, 5},
    {3, 4}, {3, 5}, {4, 5}
}};

/// Pyramid: 8 edges
/// Basix edge order:
///   (0,1), (0,2), (0,4), (1,3), (1,4), (2,3), (2,4), (3,4)
inline constexpr std::array<std::array<int, 2>, 8> pyramid_edges = {{
    {0, 1}, {0, 2}, {0, 4}, {1, 3},
    {1, 4}, {2, 3}, {2, 4}, {3, 4}
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

/// Hexahedron: 6 faces (quads), Basix face order.
inline constexpr std::array<std::array<int, 4>, 6> hexahedron_faces = {{
    {0, 1, 2, 3},
    {0, 1, 4, 5},
    {0, 2, 4, 6},
    {1, 3, 5, 7},
    {2, 3, 6, 7},
    {4, 5, 6, 7}
}};

/// Prism: 5 faces (2 triangles + 3 quads), Basix face order.
/// Stored with max 4 vertices per face, using -1 padding for triangles.
inline constexpr std::array<std::array<int, 4>, 5> prism_faces = {{
    {0, 1, 2, -1},
    {0, 1, 3, 4},
    {0, 2, 3, 5},
    {1, 2, 4, 5},
    {3, 4, 5, -1}
}};

/// Pyramid: 5 faces (1 quad + 4 triangles), Basix face order.
inline constexpr std::array<std::array<int, 4>, 5> pyramid_faces = {{
    {0, 1, 2, 3},
    {0, 1, 4, -1},
    {0, 2, 4, -1},
    {1, 3, 4, -1},
    {2, 3, 4, -1}
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

} // namespace cutcells::cell
