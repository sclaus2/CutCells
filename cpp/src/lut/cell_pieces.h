// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <span>
#include <vector>

#include "../cell_types.h"
#include "triangulation.h"

/// The lookup-table backend: a cut cell is subdivided by a Pk-iso-P1 template,
/// and each sub-cell is cut by the lookup tables (cut_<cell>) with the P1
/// interpolants of the level sets.
namespace cutcells::lut
{

/// Options of the lookup-table backend.
struct Options
{
    /// Order k of the Pk-iso-P1 template that subdivides a cut cell, 1 to 4;
    /// 0: the degree of the cell's Pk level sets, 2 for analytic ones.
    int template_order = 0;
    /// Split the cut pieces into simplices.
    bool triangulate = false;
    /// How, when triangulate is set: classical or midpoint (cell::cut).
    cell::TriangulationStrategy triangulation = cell::TriangulationStrategy::classical;
};

/// Straight pieces of one cell, each with the side of every cutting level set
/// it lies on.
template <std::floating_point T>
struct Pieces
{
    int tdim = 0;                        ///< of the cell; also the stride of vertices
    std::vector<T> vertices;             ///< reference coordinates in the cell
    std::vector<int> offsets = {0};      ///< piece p has vertices offsets[p] to offsets[p + 1] - 1
    std::vector<cell::type> types;
    std::vector<std::uint64_t> negative; ///< bit i: the piece lies where level set i < 0
    std::vector<std::uint64_t> positive; ///< bit i: where level set i > 0 (or touches 0 from above)
    std::vector<std::uint64_t> zero;     ///< bit i: the piece lies in level set i's zero set (0: a volume piece)

    int n_pieces() const { return static_cast<int>(types.size()); }
};

/// @brief The vertices of the Pk-iso-P1 template of order @p template_order
/// (1 to 4) in the reference cell, tdim coordinates each: where cut_cell takes
/// the level sets' values.
std::span<const double> template_vertices(cell::type cell_type, int template_order);

/// @brief The pieces of a cell that n level sets cut, on the Pk-iso-P1
/// template of order @p template_order.
///
/// Each template sub-cell gets the level sets' P1 interpolants of @p values
/// (multilinear on quadrilaterals), is classified by their signs at its
/// vertices, and is cut by the lookup tables one level set after the other,
/// both sides kept. Hexahedra on which the values are affine are cut whole,
/// others as their Kuhn tetrahedra: for multilinear values the hexahedron's
/// tables of both sides do not fit together. A value within 64 eps max|v| of 0 on a cell counts
/// as positive: the level set touches 0 from above there, and intersections
/// land on the vertex. With bit i of @p zero_sets, the pieces of level set i's
/// zero set are added, cut by the other level sets; with @p curves also the
/// pieces where two of those level sets vanish.
///
/// @param values         the level sets at the template's vertices,
///                       values[i * n_template_vertices + v]
/// @param triangulation  how the tables split cut pieces into simplices
///                       (none: they keep prisms, pyramids, ...)
/// @param out            pieces in the cell's reference coordinates (cleared first)
template <std::floating_point T>
void cut_cell(cell::type cell_type, int template_order, std::span<const T> values, int n_level_sets,
              std::uint64_t zero_sets, bool curves, cell::TriangulationStrategy triangulation, Pieces<T>& out);

} // namespace cutcells::lut
