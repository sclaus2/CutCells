// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <span>
#include <vector>

#include "cell_source.h"

namespace cutcells::part
{

/// @brief Bernstein coefficients of a degree-n polynomial on the two halves of
/// a tetrahedron split at the midpoint of edge (p, q), by de Casteljau's
/// algorithm along the edge.
///
/// Coefficients are in CutCells' simplex order (bernstein.h): multi-index
/// (n - i - j - k, i, j, k) at position i + offsets of (j, k), with vertex 0
/// carrying the first entry.
/// @param keep_p  the half that keeps vertex p (vertex q moves to the midpoint)
/// @param keep_q  the half that keeps vertex q (vertex p moves to the midpoint)
template <std::floating_point T>
void bisect_tetrahedron(std::span<const T> coeffs, int n, int p, int q, std::vector<T>& keep_p,
                        std::vector<T>& keep_q);

/// @brief The sign of a level set on a tetrahedron or hexahedron, for
/// classifying cells: +1 or -1 if bounds prove it (up to touching zero), 0 if
/// it changes sign or the bounds cannot tell.
///
/// Hexahedra: bounds over the box and its halves (quadrays::cell_sign).
/// Tetrahedra: bounds over the tetrahedron itself and its halves, split at the
/// midpoint of their longest edge: the Bernstein coefficients of a Pk level
/// set, or Taylor models over the parallelepiped of a piece's edges whose
/// linear part is bounded at the piece's vertices. Either way down to
/// @p max_depth splits.
template <std::floating_point T, std::integral I>
int cell_sign(const CellSource<T, I>& cs, cell::type type, int max_depth);

} // namespace cutcells::part
