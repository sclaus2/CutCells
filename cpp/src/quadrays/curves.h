// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <span>

#include "clipped_box.h"
#include "engine.h"
#include "source.h"

namespace cutcells::quadrays
{

/// @brief Quadrature points of the curve where level sets @p a and @p b vanish
/// in a cell of dimension 3, appended to @p out: no height lines, no boxes.
///
/// The curve's points on the cell's faces come from Bernstein or Taylor bounds
/// and Newton's method. From each, the curve is followed along
/// grad a x grad b to where it leaves the cell, which pairs them into arcs;
/// closed curves inside the cell are found from their extreme points along an
/// axis. Arcs are split where the gradients become parallel (two branches
/// cross, and each goes straight on) and where a further level set the terms
/// name vanishes, and they take Gauss-Legendre points in physical arc length.
/// The measure is thus right to rounding at any @p q. Each arc takes @p q points
/// with Options::moment_tolerance off or Options::tolerance below 0; otherwise as
/// many more as a moment of the mass matrix asks for, or it is split.
///
/// Values and gradients are all it reads of the curve, so polynomial and
/// analytic level sets take the same path. A crossing on the cell's boundary
/// ends the arcs that reach it.
///
/// @param cell   the cell, tdim == 3
/// @param phis   the cell's level sets
/// @param terms  the selection, as for integrate
/// @param a,b    the two level sets that vanish
/// @param q      Gauss-Legendre points per arc at least
/// @param opt    Options::moment_tolerance and Options::tolerance are read
/// @param out    the points in box coordinates, with physical weights
/// @param stats  Stats::curve_arcs and Stats::curve_fallbacks are counted
/// @return false if the curve lies in a face of the cell or could not be
///         followed: @p out is then unchanged, and the caller integrates the
///         cell another way
template <std::floating_point T>
bool curve_points(const ClippedBox<T>& cell, std::span<const Source<T>> phis, std::span<const SelectionTerm> terms,
                  int a, int b, int q, const Options& opt, CellPoints<T>& out, Stats& stats);

} // namespace cutcells::quadrays
