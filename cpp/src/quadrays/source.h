// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <concepts>
#include <span>
#include <vector>

#include "analytic.h"
#include "box_bernstein.h"
#include "clipped_box.h"

namespace cutcells::quadrays
{

/// The level set of one cell as the engine reads it: a Bernstein form on the
/// cell's box, or an analytic level set in physical coordinates with the
/// cell's map x = origin + jacobian u. Exactly one of bernstein and analytic
/// is set; both are borrowed.
///
/// Values and gradients come from either kind. For bounds and roots the engine
/// restricts a Bernstein form itself (restrict_affine); an analytic level set
/// gives them through affine_bounds and line_roots below.
template <std::floating_point T>
struct Source
{
    const BoxBernstein<T>* bernstein = nullptr;
    const AnalyticLevelSet* analytic = nullptr;
    Vec3<T> origin{};
    Mat3<T> jacobian{};

    bool is_analytic() const { return analytic != nullptr; }
};

/// Bounds of psi(s) = phi(origin + matrix s) over s in [0, 1]^m.
template <std::floating_point T>
struct AffineBounds
{
    int sign = 0;                 ///< +1 or -1 if psi has that sign everywhere, 0 if it may vanish
    T magnitude = 0;              ///< upper bound of |psi|
    bool has_derivatives = false; ///< false: lower and upper are unknown
    std::array<T, 3> lower{};     ///< lower bounds of |d psi / d s_j|, 0 if d_j psi may vanish
    std::array<T, 3> upper{};     ///< upper bounds of |d psi / d s_j|
};

/// @brief A Bernstein form on the cell's box as a source.
template <std::floating_point T>
Source<T> bernstein_source(const BoxBernstein<T>& phi);

/// @brief An analytic level set on @p cell as a source.
/// @throws std::invalid_argument if @p phi lacks value, gradient or box_bounds
template <std::floating_point T>
Source<T> analytic_source(const AnalyticLevelSet& phi, const ClippedBox<T>& cell);

/// @brief phi at box coordinates @p u.
template <std::floating_point T>
T evaluate(const Source<T>& phi, std::span<const T> u);

/// @brief Gradient of phi with respect to the box coordinates @p u, into @p g (3 entries).
template <std::floating_point T>
void gradient(const Source<T>& phi, std::span<const T> u, std::span<T> g);

/// @brief Size of phi on the cell, for tolerances relative to it: the largest
/// Bernstein coefficient, or the largest |phi| at the corners of the box.
template <std::floating_point T>
T reference_magnitude(const Source<T>& phi);

/// @brief Bounds of an analytic level set restricted to the affine image
/// u = origin + matrix s, s in [0, 1]^m, of a box (matrix 3 x m, row-major).
/// @return false if no bound holds (the value's is lost too); @p out is then
///         meaningless
template <std::floating_point T>
bool affine_bounds(const Source<T>& phi, std::span<const T> origin, std::span<const T> matrix, int m,
                   AffineBounds<T>& out);

/// @brief Roots of t -> phi(origin + t direction) on (a, b) for an analytic
/// level set, appended to @p roots: none where a Taylor model certifies the
/// value's sign or bounds |phi| by @p zero (the line lies in the zero set or
/// touches it), one where it certifies the derivative's sign and the ends
/// differ in sign; bisection otherwise, down to 2^-40 of the segment and at
/// most 4096 intervals per line.
template <std::floating_point T>
void line_roots(const Source<T>& phi, std::span<const T> origin, std::span<const T> direction, T a, T b, T zero,
                std::vector<T>& roots);

/// @brief The root of t -> phi(origin + t direction) in (a, b) for an analytic
/// level set, given its values @p ga and @p gb of opposite signs at the ends.
template <std::floating_point T>
T line_root(const Source<T>& phi, std::span<const T> origin, std::span<const T> direction, T a, T b, T ga,
            T gb);

/// @brief +1 or -1 if Bernstein coefficients prove that sign (coefficients
/// that are 0 up to rounding allowed: the polynomial may touch zero), else 0.
template <std::floating_point T>
int coefficient_sign(std::span<const T> coeffs);

/// @brief The sign of a level set on a cell, for classifying cells as inside,
/// outside or cut.
///
/// Bounds over the cell's box prove a sign (Bernstein coefficients of the
/// restriction, or Taylor models), or its corners show both; otherwise the box
/// is bisected, along its longest edge, down to @p max_depth levels, and
/// sub-boxes outside the clips are skipped.
/// @return +1 or -1 if phi has that sign on the cell (up to touching zero), 0 if
///         it changes sign or the bounds cannot tell (where it touches zero
///         inside, or where a bound is lost)
template <std::floating_point T>
int cell_sign(const ClippedBox<T>& cell, const Source<T>& phi, int max_depth = 12);

/// @brief cell_sign for an analytic level set.
template <std::floating_point T>
int cell_sign(const ClippedBox<T>& cell, const AnalyticLevelSet& phi, int max_depth = 12);

} // namespace cutcells::quadrays
