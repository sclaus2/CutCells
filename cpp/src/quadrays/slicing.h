// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <concepts>
#include <cstdint>
#include <span>
#include <vector>

#include "../selection_expr.h"
#include "box_bernstein.h"
#include "clipped_box.h"
#include "engine.h"

/// Slicing: a cell of dimension 3 integrated slice by slice.
///
/// A sweep gives plane slices of the cell. Each slice is a cell of dimension 2
/// lying in 3D, a rectangle of its plane clipped by the cell's faces, and the
/// level sets restricted to it are Bernstein forms again (restrict_affine,
/// exact). The engine integrates the slice as it does any cell of dimension 2,
/// with its own height directions, so that the edge where two zero sets meet is
/// a corner of the slice and never the projection of a curve onto a base. The
/// rule of the cell is the sum over the slices,
///
///   int_part f dx = int ( int_{slice at t} f / |grad tau| dA ) dt,
///
/// tau the sweep's coordinate: 1 / |grad tau| is 1 for parallel planes and the
/// distance to the line for planes turning about a line, never negative.
///
/// This header builds the slices and sums their rules for positions and weights
/// along the sweep that the caller gives. Where the pattern of a slice changes
/// with t (the sweep's events) the positions have to be split; finding the
/// events is not part of this module yet.
namespace cutcells::quadrays
{

/// The family of slices.
enum class SweepKind : std::uint8_t
{
    parallel = 0, ///< planes orthogonal to the direction, at distance t from the origin along it
    turning = 1   ///< half-planes bounded by the line through the origin along the direction, at angle t
};

/// A sweep in physical coordinates.
template <std::floating_point T>
struct Sweep
{
    SweepKind kind = SweepKind::parallel;
    Vec3<T> origin = {0, 0, 0};
    Vec3<T> direction = {0, 0, 1}; ///< any length
};

/// One slice of a cell: a cell of dimension 2 whose box is the bounding
/// rectangle of the slice in its plane, clipped by the cell's faces.
template <std::floating_point T>
struct Slice
{
    ClippedBox<T> box;                ///< tdim 2; reference coordinates are the box's own
    Vec3<T> origin = {0, 0, 0};       ///< u = origin + matrix s: the cell's box coordinates of the
    std::array<T, 6> matrix = {};     ///< slice's box coordinates s (3 x 2, row-major)
    T radius_origin = 0;              ///< turning: distance to the line is radius_origin + radius_step s_0
    T radius_step = 0;
};

/// @brief The slice of @p cell at position @p t of @p sweep.
/// @return false if the plane misses the cell (or meets it in less than an area)
template <std::floating_point T>
bool make_slice(const ClippedBox<T>& cell, const Sweep<T>& sweep, T t, Slice<T>& slice);

/// @brief Quadrature points of the part of a cell of dimension 3 that @p terms
/// select, as the sum of the rules of its slices, appended to @p out.
///
/// The points are in the cell's box coordinates and the weights physical, as
/// integrate gives them. All terms are volume terms (no zero set asked for).
///
/// @param cell       the cell (tdim 3)
/// @param forms      its level sets as Bernstein forms on its box
/// @param terms      the selection; bit l of a term refers to forms[l]
/// @param sweep      the family of slices
/// @param positions  positions t of the slices along the sweep
/// @param weights    their weights (a rule in t), one per position
/// @param q          Gauss-Legendre points per segment in the slices
/// @throws std::invalid_argument if the cell is not of dimension 3, a term asks
///         for a zero set, or positions and weights differ in number
template <std::floating_point T>
void integrate_slices(const ClippedBox<T>& cell, std::span<const BoxBernstein<T>> forms,
                      std::span<const SelectionTerm> terms, const Sweep<T>& sweep, std::span<const T> positions,
                      std::span<const T> weights, int q, const Options& opt, CellPoints<T>& out, Stats& stats);

} // namespace cutcells::quadrays
