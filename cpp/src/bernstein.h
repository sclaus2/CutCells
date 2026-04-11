// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <span>
#include <vector>

#include "cell_types.h"

namespace cutcells::bernstein
{

/// Number of Bernstein basis functions for a given cell type and degree.
///
/// Simplex cells (interval, triangle, tetrahedron): C(n+d, d)
/// Tensor-product cells (quadrilateral, hexahedron): (n+1)^d
int num_polynomials(cell::type ctype, int degree);

/// Evaluate a Bernstein expansion at a point xi in reference coordinates.
///
/// f(xi) = sum_i coeffs[i] * B_i(xi)
///
/// @param ctype   Cell type
/// @param degree  Polynomial degree
/// @param coeffs  Bernstein coefficients (size = num_polynomials(ctype, degree))
/// @param xi      Point in reference coordinates (size = tdim)
/// @return        Value of the expansion at xi
template <std::floating_point T>
T evaluate(cell::type ctype, int degree,
           std::span<const T> coeffs, std::span<const T> xi);

/// Gradient of a Bernstein expansion at xi in reference coordinates.
///
/// grad[k] = d/d(xi_k) sum_i coeffs[i] * B_i(xi)
///
/// @param ctype   Cell type
/// @param degree  Polynomial degree
/// @param coeffs  Bernstein coefficients (size = num_polynomials(ctype, degree))
/// @param xi      Point in reference coordinates (size = tdim)
/// @param[out] grad  Gradient vector (size = tdim)
template <std::floating_point T>
void gradient(cell::type ctype, int degree,
              std::span<const T> coeffs, std::span<const T> xi,
              std::span<T> grad);

/// Evaluate all N basis functions at a point.
///
/// out[i] = B_i(xi) for i = 0, ..., N-1
///
/// @param ctype   Cell type
/// @param degree  Polynomial degree
/// @param xi      Point in reference coordinates (size = tdim)
/// @param[out] out  Basis function values (size = num_polynomials(ctype, degree))
template <std::floating_point T>
void evaluate_basis(cell::type ctype, int degree,
                    std::span<const T> xi, std::span<T> out);

/// Convert Lagrange nodal values to Bernstein coefficients.
///
/// Builds the Bernstein evaluation matrix V (V_{ij} = B_j(ref_point_i)),
/// then solves V * coeffs = nodal_values.
///
/// @param ctype        Cell type
/// @param degree       Polynomial degree
/// @param ref_points   Reference coordinates of Lagrange nodes (flat: ndofs * tdim)
/// @param nodal_values Function values at Lagrange nodes (ndofs)
/// @param[out] coeffs  Bernstein coefficients (resized to num_polynomials)
template <std::floating_point T>
void lagrange_to_bernstein(cell::type ctype, int degree,
                           std::span<const T> ref_points,
                           std::span<const T> nodal_values,
                           std::vector<T>& coeffs);

/// Compute Bernstein coefficients of the partial derivative d/d(xi_direction)
/// of a degree-n polynomial on a simplex or tensor-product cell.
///
/// For simplex cells the result is a degree-(n-1) Bernstein expansion with
/// coefficients  c'_beta = n * (c_{beta + e_{dir+1}} - c_{beta + e_0}).
///
/// For tensor-product cells the result is a degree-(n-1) expansion in the
/// differentiated direction and degree-n in the others.
///
/// @param ctype      Cell type (interval, triangle, tetrahedron, quad, hex).
/// @param degree     Polynomial degree of the input.
/// @param coeffs     Bernstein coefficients (size = num_polynomials(ctype, degree)).
/// @param direction  Coordinate direction 0..tdim-1.
/// @param[out] deriv_coeffs  Derivative Bernstein coefficients.
template <std::floating_point T>
void derivative_coefficients(cell::type ctype, int degree,
                             std::span<const T> coeffs,
                             int direction,
                             std::vector<T>& deriv_coeffs);

// ---- Helper utilities exposed for testing ----

/// True for interval, triangle, tetrahedron.
inline bool is_simplex(cell::type ct)
{
    return ct == cell::type::interval
        || ct == cell::type::triangle
        || ct == cell::type::tetrahedron;
}

/// True for quadrilateral, hexahedron.
inline bool is_tensor_product(cell::type ct)
{
    return ct == cell::type::quadrilateral
        || ct == cell::type::hexahedron;
}

} // namespace cutcells::bernstein
