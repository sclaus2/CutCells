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
/// The slice integral is smooth in t except at the sweep's events, where the
/// pattern of a slice changes: the slice touches a zero set, an edge between
/// two zero sets (or a zero set and a face of the cell) is tangent to it, or it
/// passes a point where three of them meet. sweep_events finds them for the
/// planes u_k = const of a box, and integrate_by_slices cuts the sweep there
/// and checks each piece, for the cells it takes.
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

/// @brief The positions along box axis @p axis at which the pattern of the
/// slices u_axis = const of a part may change (events), sorted, in (0, 1).
///
/// The events are the values of u_axis at the points of the part's closure
/// where u_axis is stationary on a zero set (both other partial derivatives
/// vanish), where it is stationary on the curve in which two zero sets meet,
/// or a zero set meets a lateral face of the box, where three zero sets meet,
/// two meet on a lateral face, or one crosses a lateral edge. Each is a
/// polynomial system; boxes are bisected and dropped where the Bernstein
/// coefficients of one of its equations, or of a level set of the selection,
/// have one sign, so that no event lies outside the positions returned. A
/// position may hold no event (the bounds are not sharp), which costs slices
/// and nothing else, and events closer than 1e-5 come back as one position.
///
/// @param forms      the level sets on the box (dimension 3, no clips)
/// @param terms      the volume selection; bit l of a term refers to forms[l]
/// @param max_boxes  boxes the search may visit
/// @return false if the search does not finish within @p max_boxes, as where
///         an edge or a zero set lies in one slice, or two zero sets touch
///         along a curve; @p events is then meaningless
template <std::floating_point T>
bool sweep_events(std::span<const BoxBernstein<T>> forms, std::span<const SelectionTerm> terms, int axis,
                  std::vector<T>& events, int max_boxes = 60000);

/// @brief True if two zero sets may meet on the part's closure in the box:
/// false proves that they do not.
template <std::floating_point T>
bool may_hold_ridge(std::span<const BoxBernstein<T>> forms, std::span<const SelectionTerm> terms);

/// @brief Gauss-Legendre positions and weights on the piece [a, b] of a
/// sweep, appended: through t = a + L u^2 if @p flat_a (the end a is an event),
/// which turns the half-integer powers the slice integral has at an event into
/// integer ones, through t = b - L (1 - u)^2 if @p flat_b, plain otherwise. A
/// piece between two events is to be halved by the caller.
template <std::floating_point T>
void piece_rule(T a, bool flat_a, T b, bool flat_b, int m, std::vector<T>& positions, std::vector<T>& weights);

/// @brief integrate for a volume part of a cell that holds an edge between two
/// zero sets, by slices along the box axis with the fewest events.
///
/// Takes boxes without clips (hexahedra) with two or more level sets as
/// Bernstein forms. The sweep is cut at its events into pieces, and the slices
/// of each piece are increased, from q by half each time, until the piece's
/// integrals of the monomials of degree up to 2 in the box coordinates move by
/// less than 10^-(q + 4) of the part's volume.
/// @return false, with nothing appended, if it does not take the cell, no two
///         zero sets meet on the part, the events of no axis are found, or a
///         piece does not settle with 64 slices (as at a saddle of a zero set):
///         the caller then integrates as usual
template <std::floating_point T>
bool integrate_by_slices(const ClippedBox<T>& cell, std::span<const BoxBernstein<T>> forms,
                         std::span<const SelectionTerm> terms, int q, const Options& opt, CellPoints<T>& out,
                         Stats& stats);

} // namespace cutcells::quadrays
