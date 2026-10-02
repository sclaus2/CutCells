// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <span>
#include <vector>

#include "../cell_types.h"

namespace cutcells::lut
{

/// The map of a cell from its reference cell to physical coordinates, given by
/// the cell's vertices in Basix order: affine on simplices, multilinear on
/// quadrilaterals and hexahedra.
template <std::floating_point T>
struct CellMap
{
    cell::type type = cell::type::point;
    int gdim = 0;
    std::vector<T> vertices; ///< physical coordinates, gdim per vertex

    int tdim() const { return cell::get_tdim(type); }
};

/// @brief Whether a quadrilateral or hexahedron is a parallelogram or a
/// parallelepiped: vertex v lies at vertex 0 plus the edges from vertex 0 to
/// vertices 2^d for the bits d of v (Basix order), up to rounding.
/// @param x     the vertices, dim coordinates each
/// @param tdim  2 for a quadrilateral, 3 for a hexahedron
template <std::floating_point T>
bool is_parallelotope(std::span<const T> x, int tdim, int dim);

/// @brief The physical points of reference points.
/// @param xi  reference coordinates, tdim per point
/// @param x   physical coordinates, gdim per point (resized)
template <std::floating_point T>
void push_forward(const CellMap<T>& map, std::span<const T> xi, std::vector<T>& x);

/// @brief Append a rule on a straight piece of a cell: the reference rule of
/// the piece's type exact for degree @p degree, mapped into the cell. Prisms,
/// pyramids, and quadrilaterals and hexahedra other than parallelograms and
/// parallelepipeds are split into simplices, so that on pieces with planar
/// faces in affine cells the rule stays exact for polynomials of that degree.
///
/// Points are in the cell's reference coordinates. Weights are physical: they
/// sum to the measure of the piece's image under @p map, a volume for pieces of
/// the cell's dimension, an area or a length for lower ones, 1 for a point.
///
/// @param piece_vertices  the piece's vertices in the cell's reference
///                        coordinates, Basix order, tdim per vertex
template <std::floating_point T>
void append_piece_rule(const CellMap<T>& map, cell::type piece_type, std::span<const T> piece_vertices, int degree,
                       std::vector<T>& points, std::vector<T>& weights);

} // namespace cutcells::lut
