// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <concepts>
#include <limits>
#include <span>
#include <vector>

#include "../cell_types.h"

namespace cutcells::quadrays
{

/// Tolerance that equals @p tol in double precision and scales with the machine
/// epsilon of T otherwise.
template <std::floating_point T>
constexpr T scaled_tolerance(double tol)
{
    return static_cast<T>(tol * (static_cast<double>(std::numeric_limits<T>::epsilon())
                                 / std::numeric_limits<double>::epsilon()));
}

/// Largest degree per variable of a BoxBernstein form.
inline constexpr int max_box_degree = 64;

/// Tensor-product Bernstein polynomial on the unit box [0, 1]^dim, dim <= 3:
///
///   p(s) = sum_a coeffs[a_0 + n_0 (a_1 + n_1 a_2)]
///              B^{degree_0}_{a_0}(s_0) B^{degree_1}_{a_1}(s_1) B^{degree_2}_{a_2}(s_2),
///
/// with n_j = degree_j + 1. Variables j >= dim have degree 0.
template <std::floating_point T>
struct BoxBernstein
{
    int dim = 0;
    std::array<int, 3> degree = {0, 0, 0};
    std::vector<T> coeffs;

    int size() const { return (degree[0] + 1) * (degree[1] + 1) * (degree[2] + 1); }
};

/// @brief The level set of a cell as an exact tensor Bernstein form on the cell's box.
///
/// The box of a simplex is the unit box of its reference coordinates; the
/// polynomial is extended to the whole box. Simplices are converted by a fixed
/// matrix per degree, computed once in exact rational arithmetic; tensor cells
/// are only reordered.
///
/// @param cell_type  triangle, tetrahedron, quadrilateral or hexahedron
/// @param degree     polynomial degree (at most 12 for simplices)
/// @param coeffs     Bernstein coefficients in CutCells' order (bernstein.h)
/// @param out        form of degree @p degree in each variable
template <std::floating_point T>
void cell_bernstein_on_box(cell::type cell_type, int degree, std::span<const T> coeffs,
                           BoxBernstein<T>& out);

/// @brief Value of @p p at the point @p s (any point, not only in the box).
template <std::floating_point T>
T evaluate(const BoxBernstein<T>& p, std::span<const T> s);

/// @brief Gradient of @p p at @p s; @p grad has p.dim entries.
template <std::floating_point T>
void gradient(const BoxBernstein<T>& p, std::span<const T> s, std::span<T> grad);

/// @brief Bernstein form of the partial derivative of @p p along variable @p direction.
template <std::floating_point T>
void derivative(const BoxBernstein<T>& p, int direction, BoxBernstein<T>& out);

/// @brief Bernstein form on [0, 1]^m of s -> p(origin + matrix s).
///
/// Exact (no interpolation): each variable of @p p is an affine function of s,
/// substituted by de Casteljau's algorithm with an affine argument. The degree
/// of the result in s_j is the sum of p.degree[i] over the rows i with
/// matrix(i, j) != 0.
///
/// @param p       form in p.dim variables
/// @param origin  p.dim entries
/// @param matrix  p.dim x m entries, row-major
/// @param m       number of new variables, 0 <= m <= 3
/// @param out     result
/// @param work    scratch buffer
template <std::floating_point T>
void restrict_affine(const BoxBernstein<T>& p, std::span<const T> origin,
                     std::span<const T> matrix, int m, BoxBernstein<T>& out,
                     std::vector<T>& work);

/// @brief Bernstein form of @p p on the sub-box [lo, hi], rescaled to the unit box.
template <std::floating_point T>
void subdivide(const BoxBernstein<T>& p, std::span<const T> lo, std::span<const T> hi,
               BoxBernstein<T>& out, std::vector<T>& work);

/// @brief True if the coefficients take both signs: the polynomial may change
/// sign in the box.
///
/// Coefficients that are all >= 0 (or all <= 0) give a polynomial of one sign;
/// it can only touch zero, which needs no breakpoint and no bisection.
template <std::floating_point T>
bool may_vanish(std::span<const T> coeffs);

/// @brief Largest absolute coefficient.
template <std::floating_point T>
T max_abs(std::span<const T> coeffs);

/// @brief Euclidean norm without underflow or overflow.
template <std::floating_point T>
T scaled_norm(std::span<const T> v);

/// @brief Certified direction margins of @p p on a box of side lengths @p lengths.
///
/// For each direction k: a lower bound of |d_k p| over the box, divided by an
/// upper bound of |grad p|, both from Bernstein coefficients of the
/// derivatives (0 if d_k p may change sign). Derivatives are taken with
/// respect to the box's own coordinates y = lo + lengths * s.
///
/// @param ratio  p.dim entries
template <std::floating_point T>
void margins(const BoxBernstein<T>& p, std::span<const T> lengths, std::span<T> ratio,
             BoxBernstein<T>& work);

/// @brief Value at s of the univariate Bernstein polynomial with coefficients @p c
/// on [0, 1] (de Casteljau).
template <std::floating_point T>
T evaluate_1d(std::span<const T> c, T s);

/// @brief Root of the univariate polynomial with Bernstein coefficients @p c on [a, b]
/// inside the bracket [lo, hi] subset of [a, b], where its values @p glo and
/// @p ghi have opposite signs: bisection with secant steps (Illinois).
template <std::floating_point T>
T bracketed_root(std::span<const T> c, T a, T b, T lo, T hi, T glo, T ghi);

/// @brief All roots in (a, b) of the univariate polynomial with Bernstein
/// coefficients @p c on [a, b], isolated by de Casteljau subdivision and
/// Descartes' rule of signs, then refined by bracketed_root. Appended to @p roots.
template <std::floating_point T>
void isolate_roots(std::span<const T> c, T a, T b, std::vector<T>& roots,
                   std::vector<T>& work);

} // namespace cutcells::quadrays
