// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include "cell_types.h"
#include "edge_root.h"
#include "iso_refine.h"
#include "level_set.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace cutcells
{

enum class LocalLevelSetBackend : uint8_t
{
    nodal_signs = 0,
    bernstein = 1,
    analytical_callbacks = 2
};

template <std::floating_point T>
struct BernsteinCell
{
    cell::type cell_type = cell::type::point;
    int degree = 1;
    int tdim = 0;
    std::vector<T> coeffs;
};

template <std::floating_point T>
struct ParentPolynomialContext
{
    cell::type cell_type = cell::type::point;
    int degree = -1;
    BernsteinCell<T> bernstein;
};

int interval_bernstein_size(int p);
int triangle_bernstein_size(int p);
int tetrahedron_bernstein_size(int p);
int quadrilateral_bernstein_size(int p);
int hexahedron_bernstein_size(int p);
int bernstein_size(cell::type ct, int p);
int lagrange_node_count(cell::type ct, int degree);

int interval_i_to_flat(int p, int i);
int triangle_ijk_to_flat(int p, int i, int j, int k);
int tetrahedron_ijkl_to_flat(int p, int i, int j, int k, int l);
int quad_ij_to_flat(int p, int i, int j);
int hex_ijk_to_flat(int p, int i, int j, int k);

int infer_lagrange_order_from_num_nodes(cell::type ct, int n_vertices);
const std::vector<double>& lagrange_to_bernstein_inverse_matrix(cell::type ct, int order);

template <std::floating_point T>
inline T powi(T x, int n)
{
    if (n < 0)
        throw std::invalid_argument("powi: negative exponent");
    if (n == 0)
        return T(1);
    return std::pow(x, static_cast<T>(n));
}

template <std::floating_point T>
inline T binomial_coeff(int n, int k)
{
    if (k < 0 || k > n)
        return T(0);
    if (k == 0 || k == n)
        return T(1);
    if (k > n - k)
        k = n - k;
    T c = T(1);
    for (int i = 1; i <= k; ++i)
        c = c * static_cast<T>(n - (k - i)) / static_cast<T>(i);
    return c;
}

template <std::floating_point T>
inline T multinomial_coeff_triangle(int p, int i, int j, int k)
{
    return binomial_coeff<T>(p, i) * binomial_coeff<T>(p - i, j);
}

template <std::floating_point T>
inline T multinomial_coeff_tetrahedron(int p, int i, int j, int k, int l)
{
    return binomial_coeff<T>(p, i) * binomial_coeff<T>(p - i, j)
           * binomial_coeff<T>(p - i - j, k);
}

template <std::floating_point T>
inline T bernstein_basis_interval(int p, int i, T x)
{
    return binomial_coeff<T>(p, i) * powi(x, i) * powi(T(1) - x, p - i);
}

template <std::floating_point T>
inline std::vector<T> bernstein_coefficients_from_equispaced_nodal_1d(
    std::span<const T> nodal_values)
{
    const int p = static_cast<int>(nodal_values.size()) - 1;
    if (p < 0)
        return {};

    std::vector<T> a(static_cast<std::size_t>((p + 1) * (p + 1)), T(0));
    auto A = [&](int r, int c) -> T&
    {
        return a[static_cast<std::size_t>(r * (p + 1) + c)];
    };

    for (int j = 0; j <= p; ++j)
    {
        const T x = (p == 0) ? T(0) : static_cast<T>(j) / static_cast<T>(p);
        for (int i = 0; i <= p; ++i)
            A(j, i) = bernstein_basis_interval<T>(p, i, x);
    }

    std::vector<T> b(nodal_values.begin(), nodal_values.end());
    for (int k = 0; k <= p; ++k)
    {
        int piv = k;
        T maxv = std::abs(A(k, k));
        for (int r = k + 1; r <= p; ++r)
        {
            const T v = std::abs(A(r, k));
            if (v > maxv)
            {
                maxv = v;
                piv = r;
            }
        }
        if (maxv <= std::numeric_limits<T>::epsilon())
            return {};

        if (piv != k)
        {
            for (int c = k; c <= p; ++c)
                std::swap(A(k, c), A(piv, c));
            std::swap(b[static_cast<std::size_t>(k)], b[static_cast<std::size_t>(piv)]);
        }

        const T diag = A(k, k);
        for (int c = k; c <= p; ++c)
            A(k, c) /= diag;
        b[static_cast<std::size_t>(k)] /= diag;

        for (int r = k + 1; r <= p; ++r)
        {
            const T f = A(r, k);
            if (std::abs(f) <= std::numeric_limits<T>::epsilon())
                continue;
            for (int c = k; c <= p; ++c)
                A(r, c) -= f * A(k, c);
            b[static_cast<std::size_t>(r)] -= f * b[static_cast<std::size_t>(k)];
        }
    }

    std::vector<T> x(static_cast<std::size_t>(p + 1), T(0));
    for (int r = p; r >= 0; --r)
    {
        T s = b[static_cast<std::size_t>(r)];
        for (int c = r + 1; c <= p; ++c)
            s -= A(r, c) * x[static_cast<std::size_t>(c)];
        x[static_cast<std::size_t>(r)] = s;
    }
    return x;
}

template <std::floating_point T>
inline T bernstein_basis_triangle(int p, int i, int j, int k, T x, T y)
{
    const T l0 = T(1) - x - y;
    const T l1 = x;
    const T l2 = y;
    return multinomial_coeff_triangle<T>(p, i, j, k)
           * powi(l0, i) * powi(l1, j) * powi(l2, k);
}

template <std::floating_point T>
inline T bernstein_basis_tetrahedron(int p, int i, int j, int k, int l, T x, T y, T z)
{
    const T l0 = T(1) - x - y - z;
    const T l1 = x;
    const T l2 = y;
    const T l3 = z;
    return multinomial_coeff_tetrahedron<T>(p, i, j, k, l)
           * powi(l0, i) * powi(l1, j) * powi(l2, k) * powi(l3, l);
}

template <std::floating_point T>
inline T bernstein_basis_quadrilateral(int p, int i, int j, T x, T y)
{
    return bernstein_basis_interval<T>(p, i, x) * bernstein_basis_interval<T>(p, j, y);
}

template <std::floating_point T>
inline T bernstein_basis_hexahedron(int p, int i, int j, int k, T x, T y, T z)
{
    return bernstein_basis_interval<T>(p, i, x) * bernstein_basis_interval<T>(p, j, y)
           * bernstein_basis_interval<T>(p, k, z);
}

inline std::span<const double> bernstein_reference_nodes(cell::type ct, int order)
{
    if (order < 1)
        throw std::invalid_argument("bernstein_reference_nodes: order must be >= 1");

    if (ct == cell::type::interval)
    {
        switch (order)
        {
        case 1:
        {
            static const std::vector<double> x = {0.0, 1.0};
            return x;
        }
        case 2:
        case 3:
        case 4:
            return iso_p1_ref_coords(ct, order);
        default:
            break;
        }
    }

    if (order == 1)
        return p1_ref_coords(ct);
    return iso_p1_ref_coords(ct, order);
}

template <std::floating_point T>
inline void lagrange_to_bernstein_cell(
    cell::type ct,
    int order,
    std::span<const T> lagrange_values,
    std::span<T> bernstein_coeffs)
{
    const int n = bernstein_size(ct, order);
    if (static_cast<int>(lagrange_values.size()) != n)
        throw std::invalid_argument("lagrange_to_bernstein_cell: invalid lagrange_values size");
    if (static_cast<int>(bernstein_coeffs.size()) != n)
        throw std::invalid_argument("lagrange_to_bernstein_cell: invalid bernstein_coeffs size");

    const auto& inv = lagrange_to_bernstein_inverse_matrix(ct, order);
    for (int i = 0; i < n; ++i)
    {
        T v = T(0);
        for (int j = 0; j < n; ++j)
            v += static_cast<T>(inv[static_cast<std::size_t>(i * n + j)]) * lagrange_values[static_cast<std::size_t>(j)];
        bernstein_coeffs[static_cast<std::size_t>(i)] = v;
    }
}

template <std::floating_point T>
inline BernsteinCell<T> make_bernstein_cell(
    cell::type ct,
    int order,
    std::span<const T> lagrange_values)
{
    BernsteinCell<T> out;
    out.cell_type = ct;
    out.degree = order;
    out.tdim = cell::get_tdim(ct);
    out.coeffs.resize(static_cast<std::size_t>(bernstein_size(ct, order)));
    lagrange_to_bernstein_cell<T>(ct, order, lagrange_values,
                                  std::span<T>(out.coeffs.data(), out.coeffs.size()));
    return out;
}

template <std::floating_point T, std::integral I = int>
inline ParentPolynomialContext<T> make_parent_polynomial_context(
    cell::type                      parent_cell_type,
    const LevelSetFunction<T, I>&   level_set)
{
    if (!level_set.has_nodal_values())
        throw std::invalid_argument("make_parent_polynomial_context: nodal values are required");

    int degree = level_set.degree;
    if (degree < 1)
        degree = infer_lagrange_order_from_num_nodes(
            parent_cell_type, static_cast<int>(level_set.nodal_values.size()));

    const int expected_nodes = lagrange_node_count(parent_cell_type, degree);
    if (static_cast<int>(level_set.nodal_values.size()) != expected_nodes)
        throw std::invalid_argument(
            "make_parent_polynomial_context: nodal_values size does not match the parent degree");

    ParentPolynomialContext<T> out;
    out.cell_type = parent_cell_type;
    out.degree = degree;
    out.bernstein = make_bernstein_cell<T>(
        parent_cell_type, degree, level_set.nodal_values);
    return out;
}

template <std::floating_point T>
inline T evaluate_bernstein_cell(const BernsteinCell<T>& cell_poly, std::span<const T> x_ref)
{
    const int p = cell_poly.degree;
    if (cell_poly.tdim != static_cast<int>(x_ref.size()))
        throw std::invalid_argument("evaluate_bernstein_cell: point dimension mismatch");

    T value = T(0);
    switch (cell_poly.cell_type)
    {
    case cell::type::interval:
        for (int i = 0; i <= p; ++i)
        {
            const int idx = interval_i_to_flat(p, i);
            value += cell_poly.coeffs[static_cast<std::size_t>(idx)]
                     * bernstein_basis_interval<T>(p, i, x_ref[0]);
        }
        return value;
    case cell::type::triangle:
        for (int i = 0; i <= p; ++i)
        {
            for (int j = 0; j <= p - i; ++j)
            {
                const int k = p - i - j;
                const int idx = triangle_ijk_to_flat(p, i, j, k);
                value += cell_poly.coeffs[static_cast<std::size_t>(idx)]
                         * bernstein_basis_triangle<T>(p, i, j, k, x_ref[0], x_ref[1]);
            }
        }
        return value;
    case cell::type::quadrilateral:
        for (int j = 0; j <= p; ++j)
        {
            for (int i = 0; i <= p; ++i)
            {
                const int idx = quad_ij_to_flat(p, i, j);
                value += cell_poly.coeffs[static_cast<std::size_t>(idx)]
                         * bernstein_basis_quadrilateral<T>(p, i, j, x_ref[0], x_ref[1]);
            }
        }
        return value;
    case cell::type::tetrahedron:
        for (int i = 0; i <= p; ++i)
        {
            for (int j = 0; j <= p - i; ++j)
            {
                for (int k = 0; k <= p - i - j; ++k)
                {
                    const int l = p - i - j - k;
                    const int idx = tetrahedron_ijkl_to_flat(p, i, j, k, l);
                    value += cell_poly.coeffs[static_cast<std::size_t>(idx)]
                             * bernstein_basis_tetrahedron<T>(
                                 p, i, j, k, l, x_ref[0], x_ref[1], x_ref[2]);
                }
            }
        }
        return value;
    case cell::type::hexahedron:
        for (int k = 0; k <= p; ++k)
        {
            for (int j = 0; j <= p; ++j)
            {
                for (int i = 0; i <= p; ++i)
                {
                    const int idx = hex_ijk_to_flat(p, i, j, k);
                    value += cell_poly.coeffs[static_cast<std::size_t>(idx)]
                             * bernstein_basis_hexahedron<T>(p, i, j, k, x_ref[0], x_ref[1], x_ref[2]);
                }
            }
        }
        return value;
    default:
        throw std::invalid_argument("evaluate_bernstein_cell: unsupported cell type");
    }
}

template <std::floating_point T>
inline void extract_edge_bernstein_coeffs(
    const BernsteinCell<T>& cell_poly,
    int edge_id,
    std::vector<T>& edge_coeffs)
{
    const int p = cell_poly.degree;
    edge_coeffs.assign(static_cast<std::size_t>(p + 1), T(0));

    switch (cell_poly.cell_type)
    {
    case cell::type::interval:
        if (edge_id != 0)
            throw std::invalid_argument("extract_edge_bernstein_coeffs: invalid interval edge id");
        for (int i = 0; i <= p; ++i)
            edge_coeffs[static_cast<std::size_t>(i)] = cell_poly.coeffs[static_cast<std::size_t>(i)];
        return;
    case cell::type::triangle:
        switch (edge_id)
        {
        case 0:
            for (int k = 0; k <= p; ++k)
                edge_coeffs[static_cast<std::size_t>(k)]
                    = cell_poly.coeffs[static_cast<std::size_t>(triangle_ijk_to_flat(p, 0, p - k, k))];
            return;
        case 1:
            for (int k = 0; k <= p; ++k)
                edge_coeffs[static_cast<std::size_t>(k)]
                    = cell_poly.coeffs[static_cast<std::size_t>(triangle_ijk_to_flat(p, p - k, 0, k))];
            return;
        case 2:
            for (int k = 0; k <= p; ++k)
                edge_coeffs[static_cast<std::size_t>(k)]
                    = cell_poly.coeffs[static_cast<std::size_t>(triangle_ijk_to_flat(p, p - k, k, 0))];
            return;
        default:
            throw std::invalid_argument("extract_edge_bernstein_coeffs: invalid triangle edge id");
        }
    case cell::type::quadrilateral:
        switch (edge_id)
        {
        case 0:
            for (int i = 0; i <= p; ++i)
                edge_coeffs[static_cast<std::size_t>(i)]
                    = cell_poly.coeffs[static_cast<std::size_t>(quad_ij_to_flat(p, i, 0))];
            return;
        case 1:
            for (int j = 0; j <= p; ++j)
                edge_coeffs[static_cast<std::size_t>(j)]
                    = cell_poly.coeffs[static_cast<std::size_t>(quad_ij_to_flat(p, p, j))];
            return;
        case 2:
            for (int i = 0; i <= p; ++i)
                edge_coeffs[static_cast<std::size_t>(i)]
                    = cell_poly.coeffs[static_cast<std::size_t>(quad_ij_to_flat(p, i, p))];
            return;
        case 3:
            for (int j = 0; j <= p; ++j)
                edge_coeffs[static_cast<std::size_t>(j)]
                    = cell_poly.coeffs[static_cast<std::size_t>(quad_ij_to_flat(p, 0, j))];
            return;
        default:
            throw std::invalid_argument("extract_edge_bernstein_coeffs: invalid quadrilateral edge id");
        }
    case cell::type::tetrahedron:
    {
        const auto set = [&](int a, int b, int k)
        {
            std::array<int, 4> exp = {0, 0, 0, 0};
            exp[static_cast<std::size_t>(a)] = p - k;
            exp[static_cast<std::size_t>(b)] = k;
            return tetrahedron_ijkl_to_flat(p, exp[0], exp[1], exp[2], exp[3]);
        };
        switch (edge_id)
        {
        case 0:
            for (int k = 0; k <= p; ++k)
                edge_coeffs[static_cast<std::size_t>(k)] = cell_poly.coeffs[static_cast<std::size_t>(set(2, 3, k))];
            return;
        case 1:
            for (int k = 0; k <= p; ++k)
                edge_coeffs[static_cast<std::size_t>(k)] = cell_poly.coeffs[static_cast<std::size_t>(set(1, 3, k))];
            return;
        case 2:
            for (int k = 0; k <= p; ++k)
                edge_coeffs[static_cast<std::size_t>(k)] = cell_poly.coeffs[static_cast<std::size_t>(set(1, 2, k))];
            return;
        case 3:
            for (int k = 0; k <= p; ++k)
                edge_coeffs[static_cast<std::size_t>(k)] = cell_poly.coeffs[static_cast<std::size_t>(set(0, 3, k))];
            return;
        case 4:
            for (int k = 0; k <= p; ++k)
                edge_coeffs[static_cast<std::size_t>(k)] = cell_poly.coeffs[static_cast<std::size_t>(set(0, 2, k))];
            return;
        case 5:
            for (int k = 0; k <= p; ++k)
                edge_coeffs[static_cast<std::size_t>(k)] = cell_poly.coeffs[static_cast<std::size_t>(set(0, 1, k))];
            return;
        default:
            throw std::invalid_argument("extract_edge_bernstein_coeffs: invalid tetrahedron edge id");
        }
    }
    case cell::type::hexahedron:
        switch (edge_id)
        {
        case 0:
            for (int i = 0; i <= p; ++i)
                edge_coeffs[static_cast<std::size_t>(i)] = cell_poly.coeffs[static_cast<std::size_t>(hex_ijk_to_flat(p, i, 0, 0))];
            return;
        case 1:
            for (int j = 0; j <= p; ++j)
                edge_coeffs[static_cast<std::size_t>(j)] = cell_poly.coeffs[static_cast<std::size_t>(hex_ijk_to_flat(p, 0, j, 0))];
            return;
        case 2:
            for (int k = 0; k <= p; ++k)
                edge_coeffs[static_cast<std::size_t>(k)] = cell_poly.coeffs[static_cast<std::size_t>(hex_ijk_to_flat(p, 0, 0, k))];
            return;
        case 3:
            for (int j = 0; j <= p; ++j)
                edge_coeffs[static_cast<std::size_t>(j)] = cell_poly.coeffs[static_cast<std::size_t>(hex_ijk_to_flat(p, p, j, 0))];
            return;
        case 4:
            for (int k = 0; k <= p; ++k)
                edge_coeffs[static_cast<std::size_t>(k)] = cell_poly.coeffs[static_cast<std::size_t>(hex_ijk_to_flat(p, p, 0, k))];
            return;
        case 5:
            for (int i = 0; i <= p; ++i)
                edge_coeffs[static_cast<std::size_t>(i)] = cell_poly.coeffs[static_cast<std::size_t>(hex_ijk_to_flat(p, i, p, 0))];
            return;
        case 6:
            for (int k = 0; k <= p; ++k)
                edge_coeffs[static_cast<std::size_t>(k)] = cell_poly.coeffs[static_cast<std::size_t>(hex_ijk_to_flat(p, 0, p, k))];
            return;
        case 7:
            for (int k = 0; k <= p; ++k)
                edge_coeffs[static_cast<std::size_t>(k)] = cell_poly.coeffs[static_cast<std::size_t>(hex_ijk_to_flat(p, p, p, k))];
            return;
        case 8:
            for (int i = 0; i <= p; ++i)
                edge_coeffs[static_cast<std::size_t>(i)] = cell_poly.coeffs[static_cast<std::size_t>(hex_ijk_to_flat(p, i, 0, p))];
            return;
        case 9:
            for (int j = 0; j <= p; ++j)
                edge_coeffs[static_cast<std::size_t>(j)] = cell_poly.coeffs[static_cast<std::size_t>(hex_ijk_to_flat(p, 0, j, p))];
            return;
        case 10:
            for (int j = 0; j <= p; ++j)
                edge_coeffs[static_cast<std::size_t>(j)] = cell_poly.coeffs[static_cast<std::size_t>(hex_ijk_to_flat(p, p, j, p))];
            return;
        case 11:
            for (int i = 0; i <= p; ++i)
                edge_coeffs[static_cast<std::size_t>(i)] = cell_poly.coeffs[static_cast<std::size_t>(hex_ijk_to_flat(p, i, p, p))];
            return;
        default:
            throw std::invalid_argument("extract_edge_bernstein_coeffs: invalid hexahedron edge id");
        }
    default:
        throw std::invalid_argument("extract_edge_bernstein_coeffs: unsupported cell type");
    }
}

template <std::floating_point T>
inline void bernstein_derivative_1d(std::span<const T> coeffs, std::vector<T>& dcoeffs)
{
    if (coeffs.empty())
    {
        dcoeffs.clear();
        return;
    }

    const int p = static_cast<int>(coeffs.size()) - 1;
    dcoeffs.assign(static_cast<std::size_t>(std::max(0, p)), T(0));
    for (int i = 0; i < p; ++i)
        dcoeffs[static_cast<std::size_t>(i)] = static_cast<T>(p)
                                               * (coeffs[static_cast<std::size_t>(i + 1)]
                                                  - coeffs[static_cast<std::size_t>(i)]);
}

template <std::floating_point T>
inline T de_casteljau_eval_1d(std::span<const T> coeffs, T t)
{
    if (coeffs.empty())
        throw std::invalid_argument("de_casteljau_eval_1d: empty coefficient array");

    std::vector<T> work(coeffs.begin(), coeffs.end());
    const int p = static_cast<int>(work.size()) - 1;
    for (int r = 1; r <= p; ++r)
    {
        for (int i = 0; i <= p - r; ++i)
        {
            work[static_cast<std::size_t>(i)] = (T(1) - t) * work[static_cast<std::size_t>(i)]
                                                + t * work[static_cast<std::size_t>(i + 1)];
        }
    }
    return work[0];
}

template <std::floating_point T>
inline void de_casteljau_split_1d(
    std::span<const T> parent,
    std::vector<T>& left,
    std::vector<T>& right,
    T t_split = T(0.5))
{
    if (parent.empty())
        throw std::invalid_argument("de_casteljau_split_1d: empty coefficient array");

    const int p = static_cast<int>(parent.size()) - 1;
    std::vector<T> work(parent.begin(), parent.end());
    left.resize(static_cast<std::size_t>(p + 1));
    right.resize(static_cast<std::size_t>(p + 1));
    left[0] = work[0];
    right[static_cast<std::size_t>(p)] = work[static_cast<std::size_t>(p)];
    for (int r = 1; r <= p; ++r)
    {
        for (int i = 0; i <= p - r; ++i)
        {
            work[static_cast<std::size_t>(i)] = (T(1) - t_split) * work[static_cast<std::size_t>(i)]
                                                + t_split * work[static_cast<std::size_t>(i + 1)];
        }
        left[static_cast<std::size_t>(r)] = work[0];
        right[static_cast<std::size_t>(p - r)] = work[static_cast<std::size_t>(p - r)];
    }
}

template <std::floating_point T>
inline bool bernstein_all_strict_one_sign(std::span<const T> coeffs, T tol)
{
    bool all_pos = true;
    bool all_neg = true;
    for (const T c : coeffs)
    {
        all_pos = all_pos && (c > tol);
        all_neg = all_neg && (c < -tol);
    }
    return all_pos || all_neg;
}

template <std::floating_point T>
inline bool bernstein_has_near_zero(std::span<const T> coeffs, T tol)
{
    for (const T c : coeffs)
    {
        if (std::abs(c) <= tol)
            return true;
    }
    return false;
}

template <std::floating_point T>
inline bool bernstein_derivative_strict_one_sign(std::span<const T> coeffs, T tol)
{
    std::vector<T> dcoeffs;
    bernstein_derivative_1d(coeffs, dcoeffs);
    if (dcoeffs.empty())
        return true;
    return bernstein_all_strict_one_sign(std::span<const T>(dcoeffs.data(), dcoeffs.size()), tol);
}

template <std::floating_point T>
inline int bernstein_sign_variation_count(std::span<const T> coeffs, T tol)
{
    int count = 0;
    int prev_sign = 0;
    for (const T c : coeffs)
    {
        const int s = (c > tol) ? 1 : ((c < -tol) ? -1 : 0);
        if (s == 0)
            continue;
        if (prev_sign != 0 && s != prev_sign)
            ++count;
        prev_sign = s;
    }
    return count;
}

// ============================================================================
// Exact De Casteljau / blossom-based segment restriction helpers
// ============================================================================

/// Apply one multivariate blend step to a degree-r triangle Bernstein array.
/// Input:  work[], size (r+1)(r+2)/2, Bernstein coefficients at degree r.
/// Output: out[],  size r(r+1)/2,     Bernstein coefficients at degree r-1.
/// Formula: out[β] = Σ_k lambda[k] * work[β + e_k]   for |β| = r-1.
template <std::floating_point T>
inline void simplex_bernstein_blend_step_triangle(
    int                    r,
    std::span<const T>     work,
    const std::array<T,3>& lambda,
    std::vector<T>&        out)
{
    const int r1 = r - 1;
    out.resize(static_cast<std::size_t>(r1 + 1) * static_cast<std::size_t>(r1 + 2) / 2);
    for (int b0 = 0; b0 <= r1; ++b0)
        for (int b1 = 0; b1 <= r1 - b0; ++b1)
        {
            const int b2      = r1 - b0 - b1;
            const int idx_out = triangle_ijk_to_flat(r1, b0, b1, b2);
            out[static_cast<std::size_t>(idx_out)] =
                lambda[0] * work[static_cast<std::size_t>(triangle_ijk_to_flat(r, b0 + 1, b1,     b2    ))]
              + lambda[1] * work[static_cast<std::size_t>(triangle_ijk_to_flat(r, b0,     b1 + 1, b2    ))]
              + lambda[2] * work[static_cast<std::size_t>(triangle_ijk_to_flat(r, b0,     b1,     b2 + 1))];
        }
}

/// Apply one multivariate blend step to a degree-r tetrahedron Bernstein array.
/// Input:  work[], size (r+1)(r+2)(r+3)/6.
/// Output: out[],  size r(r+1)(r+2)/6.
template <std::floating_point T>
inline void simplex_bernstein_blend_step_tetrahedron(
    int                    r,
    std::span<const T>     work,
    const std::array<T,4>& lambda,
    std::vector<T>&        out)
{
    const int r1 = r - 1;
    out.resize(static_cast<std::size_t>(r1 + 1)
               * static_cast<std::size_t>(r1 + 2)
               * static_cast<std::size_t>(r1 + 3) / 6);
    for (int b0 = 0; b0 <= r1; ++b0)
        for (int b1 = 0; b1 <= r1 - b0; ++b1)
            for (int b2 = 0; b2 <= r1 - b0 - b1; ++b2)
            {
                const int b3      = r1 - b0 - b1 - b2;
                const int idx_out = tetrahedron_ijkl_to_flat(r1, b0, b1, b2, b3);
                out[static_cast<std::size_t>(idx_out)] =
                    lambda[0] * work[static_cast<std::size_t>(tetrahedron_ijkl_to_flat(r, b0 + 1, b1,     b2,     b3    ))]
                  + lambda[1] * work[static_cast<std::size_t>(tetrahedron_ijkl_to_flat(r, b0,     b1 + 1, b2,     b3    ))]
                  + lambda[2] * work[static_cast<std::size_t>(tetrahedron_ijkl_to_flat(r, b0,     b1,     b2 + 1, b3    ))]
                  + lambda[3] * work[static_cast<std::size_t>(tetrahedron_ijkl_to_flat(r, b0,     b1,     b2,     b3 + 1))];
            }
}

/// Exact restriction of a triangle Bernstein polynomial to an arbitrary segment.
///
/// Uses the multivariate blossom formula via De Casteljau blend steps.
/// Numerically exact (no Vandermonde solve): the j-th output Bernstein coefficient is
///   d[j] = blossom(lam_end^j, lam_start^{p-j})
/// computed by sharing the lam_end prefix across all j.
template <std::floating_point T>
inline void restrict_bernstein_triangle_to_segment_exact(
    const BernsteinCell<T>& cell_poly,
    std::span<const T>      x0_ref,     // length 2
    std::span<const T>      x1_ref,     // length 2
    std::vector<T>&         edge_coeffs)
{
    const int p = cell_poly.degree;
    edge_coeffs.resize(static_cast<std::size_t>(p + 1));

    // Barycentric coords: λ0 = 1-x-y, λ1 = x, λ2 = y
    const std::array<T, 3> lam_start = {
        T(1) - x0_ref[0] - x0_ref[1], x0_ref[0], x0_ref[1]};
    const std::array<T, 3> lam_end = {
        T(1) - x1_ref[0] - x1_ref[1], x1_ref[0], x1_ref[1]};

    // D_end[j] = j sequential lam_end blend steps from parent coefficients (degree p-j).
    // Built once and reused for each d[j].
    std::vector<std::vector<T>> D_end(static_cast<std::size_t>(p + 1));
    D_end[0].assign(cell_poly.coeffs.begin(), cell_poly.coeffs.end());
    for (int j = 1; j <= p; ++j)
        simplex_bernstein_blend_step_triangle<T>(
            p - (j - 1),
            std::span<const T>(D_end[j - 1].data(), D_end[j - 1].size()),
            lam_end, D_end[j]);

    // d[j] = apply (p-j) lam_start steps to D_end[j], reaching degree 0 (scalar)
    for (int j = 0; j <= p; ++j)
    {
        if (p - j == 0)
        {
            edge_coeffs[static_cast<std::size_t>(j)] = D_end[j][0];
            continue;
        }
        std::vector<T> work = D_end[j];
        std::vector<T> tmp;
        int cur_deg = p - j;
        for (int k = 0; k < p - j; ++k)
        {
            simplex_bernstein_blend_step_triangle<T>(
                cur_deg,
                std::span<const T>(work.data(), work.size()),
                lam_start, tmp);
            --cur_deg;
            std::swap(work, tmp);
        }
        edge_coeffs[static_cast<std::size_t>(j)] = work[0];
    }
}

/// Exact restriction of a tetrahedron Bernstein polynomial to an arbitrary segment.
/// Same blossom approach as the triangle case, using 4-component barycentric coords.
template <std::floating_point T>
inline void restrict_bernstein_tetrahedron_to_segment_exact(
    const BernsteinCell<T>& cell_poly,
    std::span<const T>      x0_ref,     // length 3
    std::span<const T>      x1_ref,     // length 3
    std::vector<T>&         edge_coeffs)
{
    const int p = cell_poly.degree;
    edge_coeffs.resize(static_cast<std::size_t>(p + 1));

    // Barycentric coords: λ0 = 1-x-y-z, λ1 = x, λ2 = y, λ3 = z
    const std::array<T, 4> lam_start = {
        T(1) - x0_ref[0] - x0_ref[1] - x0_ref[2],
        x0_ref[0], x0_ref[1], x0_ref[2]};
    const std::array<T, 4> lam_end = {
        T(1) - x1_ref[0] - x1_ref[1] - x1_ref[2],
        x1_ref[0], x1_ref[1], x1_ref[2]};

    std::vector<std::vector<T>> D_end(static_cast<std::size_t>(p + 1));
    D_end[0].assign(cell_poly.coeffs.begin(), cell_poly.coeffs.end());
    for (int j = 1; j <= p; ++j)
        simplex_bernstein_blend_step_tetrahedron<T>(
            p - (j - 1),
            std::span<const T>(D_end[j - 1].data(), D_end[j - 1].size()),
            lam_end, D_end[j]);

    for (int j = 0; j <= p; ++j)
    {
        if (p - j == 0)
        {
            edge_coeffs[static_cast<std::size_t>(j)] = D_end[j][0];
            continue;
        }
        std::vector<T> work = D_end[j];
        std::vector<T> tmp;
        int cur_deg = p - j;
        for (int k = 0; k < p - j; ++k)
        {
            simplex_bernstein_blend_step_tetrahedron<T>(
                cur_deg,
                std::span<const T>(work.data(), work.size()),
                lam_start, tmp);
            --cur_deg;
            std::swap(work, tmp);
        }
        edge_coeffs[static_cast<std::size_t>(j)] = work[0];
    }
}

/// Exact restriction of a quadrilateral Bernstein polynomial to an axis-aligned segment.
///
/// Red-refined quad sub-edges are always axis-aligned in reference coordinates.
/// Throws std::invalid_argument for non-axis-aligned segments.
///
/// Algorithm:
///   1. Collapse the constant direction by evaluating each 1-D slice with De Casteljau.
///   2. Restrict the varying direction to the sub-interval [a,b] with De Casteljau splits.
template <std::floating_point T>
inline void restrict_bernstein_quad_to_segment_exact(
    const BernsteinCell<T>& cell_poly,
    std::span<const T>      x0_ref,     // length 2
    std::span<const T>      x1_ref,     // length 2
    std::vector<T>&         edge_coeffs)
{
    const int p   = cell_poly.degree;
    const T   tol = T(1e-10);

    const T dx = x1_ref[0] - x0_ref[0];
    const T dy = x1_ref[1] - x0_ref[1];
    const bool axis_x = (std::abs(dy) <= tol) && (std::abs(dx) > tol);
    const bool axis_y = (std::abs(dx) <= tol) && (std::abs(dy) > tol);

    if (!axis_x && !axis_y)
        throw std::invalid_argument(
            "restrict_bernstein_quad_to_segment_exact: non-axis-aligned segment");

    // Step 1: collapse the constant direction
    std::vector<T> collapsed(static_cast<std::size_t>(p + 1));
    if (axis_x)
    {
        // y is constant; slice over i (x basis index), eval each y-column at y0
        const T y0 = x0_ref[1];
        for (int i = 0; i <= p; ++i)
        {
            std::vector<T> col(static_cast<std::size_t>(p + 1));
            for (int j = 0; j <= p; ++j)
                col[static_cast<std::size_t>(j)] =
                    cell_poly.coeffs[static_cast<std::size_t>(quad_ij_to_flat(p, i, j))];
            collapsed[static_cast<std::size_t>(i)] =
                de_casteljau_eval_1d<T>(std::span<const T>(col.data(), col.size()), y0);
        }
    }
    else
    {
        // x is constant; slice over j (y basis index), eval each x-row at x0
        const T x0c = x0_ref[0];
        for (int j = 0; j <= p; ++j)
        {
            std::vector<T> row(static_cast<std::size_t>(p + 1));
            for (int i = 0; i <= p; ++i)
                row[static_cast<std::size_t>(i)] =
                    cell_poly.coeffs[static_cast<std::size_t>(quad_ij_to_flat(p, i, j))];
            collapsed[static_cast<std::size_t>(j)] =
                de_casteljau_eval_1d<T>(std::span<const T>(row.data(), row.size()), x0c);
        }
    }

    // Step 2: restrict varying direction to sub-interval [a,b] (a <= b)
    const T a = axis_x ? std::min(x0_ref[0], x1_ref[0]) : std::min(x0_ref[1], x1_ref[1]);
    const T b = axis_x ? std::max(x0_ref[0], x1_ref[0]) : std::max(x0_ref[1], x1_ref[1]);

    std::vector<T> working(collapsed.begin(), collapsed.end());
    if (a > tol)
    {
        std::vector<T> left, right;
        de_casteljau_split_1d<T>(
            std::span<const T>(working.data(), working.size()), left, right, a);
        working = std::move(right);
    }
    const T remaining = T(1) - a;
    if (remaining > tol && b < T(1) - tol)
    {
        const T t_split = (b - a) / remaining;
        std::vector<T> left, right;
        de_casteljau_split_1d<T>(
            std::span<const T>(working.data(), working.size()), left, right, t_split);
        working = std::move(left);
    }

    edge_coeffs = std::move(working);
    // Orient coefficients to match the direction from x0 to x1
    if ((axis_x && x0_ref[0] > x1_ref[0]) || (axis_y && x0_ref[1] > x1_ref[1]))
        std::reverse(edge_coeffs.begin(), edge_coeffs.end());
}

/// Exact restriction of a hexahedron Bernstein polynomial to an axis-aligned segment.
///
/// Red-refined hex sub-edges are always axis-aligned in reference coordinates.
/// Throws std::invalid_argument for non-axis-aligned segments.
///
/// Algorithm: collapse the two constant directions with De Casteljau evaluation,
/// then restrict the varying direction to its sub-interval.
template <std::floating_point T>
inline void restrict_bernstein_hex_to_segment_exact(
    const BernsteinCell<T>& cell_poly,
    std::span<const T>      x0_ref,     // length 3
    std::span<const T>      x1_ref,     // length 3
    std::vector<T>&         edge_coeffs)
{
    const int p   = cell_poly.degree;
    const T   tol = T(1e-10);

    const T dx = x1_ref[0] - x0_ref[0];
    const T dy = x1_ref[1] - x0_ref[1];
    const T dz = x1_ref[2] - x0_ref[2];
    const bool axis_x = (std::abs(dy) <= tol) && (std::abs(dz) <= tol) && (std::abs(dx) > tol);
    const bool axis_y = (std::abs(dx) <= tol) && (std::abs(dz) <= tol) && (std::abs(dy) > tol);
    const bool axis_z = (std::abs(dx) <= tol) && (std::abs(dy) <= tol) && (std::abs(dz) > tol);

    if (!axis_x && !axis_y && !axis_z)
        throw std::invalid_argument(
            "restrict_bernstein_hex_to_segment_exact: non-axis-aligned segment");

    // Step 1: collapse the two constant directions
    std::vector<T> collapsed(static_cast<std::size_t>(p + 1));
    if (axis_x)
    {
        // Vary i; constant y0, z0
        const T y0 = x0_ref[1];
        const T z0 = x0_ref[2];
        for (int i = 0; i <= p; ++i)
        {
            // First collapse k-direction at z0 for each j
            std::vector<T> slice_j(static_cast<std::size_t>(p + 1));
            for (int j = 0; j <= p; ++j)
            {
                std::vector<T> col_k(static_cast<std::size_t>(p + 1));
                for (int k = 0; k <= p; ++k)
                    col_k[static_cast<std::size_t>(k)] =
                        cell_poly.coeffs[static_cast<std::size_t>(hex_ijk_to_flat(p, i, j, k))];
                slice_j[static_cast<std::size_t>(j)] =
                    de_casteljau_eval_1d<T>(std::span<const T>(col_k.data(), col_k.size()), z0);
            }
            // Then collapse j-direction at y0
            collapsed[static_cast<std::size_t>(i)] =
                de_casteljau_eval_1d<T>(std::span<const T>(slice_j.data(), slice_j.size()), y0);
        }
    }
    else if (axis_y)
    {
        // Vary j; constant x0, z0
        const T x0c = x0_ref[0];
        const T z0  = x0_ref[2];
        for (int j = 0; j <= p; ++j)
        {
            std::vector<T> slice_i(static_cast<std::size_t>(p + 1));
            for (int i = 0; i <= p; ++i)
            {
                std::vector<T> col_k(static_cast<std::size_t>(p + 1));
                for (int k = 0; k <= p; ++k)
                    col_k[static_cast<std::size_t>(k)] =
                        cell_poly.coeffs[static_cast<std::size_t>(hex_ijk_to_flat(p, i, j, k))];
                slice_i[static_cast<std::size_t>(i)] =
                    de_casteljau_eval_1d<T>(std::span<const T>(col_k.data(), col_k.size()), z0);
            }
            collapsed[static_cast<std::size_t>(j)] =
                de_casteljau_eval_1d<T>(std::span<const T>(slice_i.data(), slice_i.size()), x0c);
        }
    }
    else  // axis_z
    {
        // Vary k; constant x0, y0
        const T x0c = x0_ref[0];
        const T y0  = x0_ref[1];
        for (int k = 0; k <= p; ++k)
        {
            std::vector<T> slice_i(static_cast<std::size_t>(p + 1));
            for (int i = 0; i <= p; ++i)
            {
                std::vector<T> col_j(static_cast<std::size_t>(p + 1));
                for (int j = 0; j <= p; ++j)
                    col_j[static_cast<std::size_t>(j)] =
                        cell_poly.coeffs[static_cast<std::size_t>(hex_ijk_to_flat(p, i, j, k))];
                slice_i[static_cast<std::size_t>(i)] =
                    de_casteljau_eval_1d<T>(std::span<const T>(col_j.data(), col_j.size()), y0);
            }
            collapsed[static_cast<std::size_t>(k)] =
                de_casteljau_eval_1d<T>(std::span<const T>(slice_i.data(), slice_i.size()), x0c);
        }
    }

    // Step 2: restrict varying direction to sub-interval [a,b]
    T    a, b;
    bool reversed;
    if (axis_x)
    {
        a = std::min(x0_ref[0], x1_ref[0]);
        b = std::max(x0_ref[0], x1_ref[0]);
        reversed = (x0_ref[0] > x1_ref[0]);
    }
    else if (axis_y)
    {
        a = std::min(x0_ref[1], x1_ref[1]);
        b = std::max(x0_ref[1], x1_ref[1]);
        reversed = (x0_ref[1] > x1_ref[1]);
    }
    else
    {
        a = std::min(x0_ref[2], x1_ref[2]);
        b = std::max(x0_ref[2], x1_ref[2]);
        reversed = (x0_ref[2] > x1_ref[2]);
    }

    std::vector<T> working(collapsed.begin(), collapsed.end());
    if (a > tol)
    {
        std::vector<T> left, right;
        de_casteljau_split_1d<T>(
            std::span<const T>(working.data(), working.size()), left, right, a);
        working = std::move(right);
    }
    const T remaining = T(1) - a;
    if (remaining > tol && b < T(1) - tol)
    {
        const T t_split = (b - a) / remaining;
        std::vector<T> left, right;
        de_casteljau_split_1d<T>(
            std::span<const T>(working.data(), working.size()), left, right, t_split);
        working = std::move(left);
    }

    edge_coeffs = std::move(working);
    if (reversed)
        std::reverse(edge_coeffs.begin(), edge_coeffs.end());
}

// ============================================================================
// restrict_bernstein_to_segment_1d — dispatches to exact per-cell-type functions
// ============================================================================

template <std::floating_point T>
inline void restrict_bernstein_to_segment_1d(
    const BernsteinCell<T>& cell_poly,
    std::span<const T> x0_ref,
    std::span<const T> x1_ref,
    std::vector<T>& edge_coeffs)
{
    if (static_cast<int>(x0_ref.size()) != cell_poly.tdim
        || static_cast<int>(x1_ref.size()) != cell_poly.tdim)
        throw std::invalid_argument("restrict_bernstein_to_segment_1d: segment dimension mismatch");

    switch (cell_poly.cell_type)
    {
    case cell::type::triangle:
        restrict_bernstein_triangle_to_segment_exact<T>(cell_poly, x0_ref, x1_ref, edge_coeffs);
        return;
    case cell::type::tetrahedron:
        restrict_bernstein_tetrahedron_to_segment_exact<T>(cell_poly, x0_ref, x1_ref, edge_coeffs);
        return;
    case cell::type::quadrilateral:
        restrict_bernstein_quad_to_segment_exact<T>(cell_poly, x0_ref, x1_ref, edge_coeffs);
        return;
    case cell::type::hexahedron:
        restrict_bernstein_hex_to_segment_exact<T>(cell_poly, x0_ref, x1_ref, edge_coeffs);
        return;
    default:
        throw std::invalid_argument(
            "restrict_bernstein_to_segment_1d: unsupported cell type");
    }
}

template <std::floating_point T>
inline cell::edge_root::RootSolveInfo<T> bernstein_edge_root_info(
    std::span<const T> edge_coeffs,
    int max_depth,
    T xtol,
    T ftol)
{
    if (edge_coeffs.empty())
        throw std::invalid_argument("bernstein_edge_root_info: empty coefficient array");
    if (max_depth < 0)
        throw std::invalid_argument("bernstein_edge_root_info: max_depth must be >= 0");

    cell::edge_root::RootSolveInfo<T> out;
    out.t = T(0.5);
    out.residual = std::numeric_limits<T>::max();

    const T f0 = edge_coeffs.front();
    const T f1 = edge_coeffs.back();

    if (bernstein_all_strict_one_sign(edge_coeffs, ftol))
    {
        out.residual = std::min(std::abs(f0), std::abs(f1));
        return out;
    }

    auto bisect = [&](T a, T b) -> cell::edge_root::RootSolveInfo<T>
    {
        cell::edge_root::RootSolveInfo<T> info;
        T fa = de_casteljau_eval_1d<T>(edge_coeffs, a);
        T fb = de_casteljau_eval_1d<T>(edge_coeffs, b);
        info.evaluations += 2;

        if (std::abs(fa) <= ftol)
        {
            info.t = a;
            info.residual = std::abs(fa);
            info.converged = true;
            return info;
        }
        if (std::abs(fb) <= ftol)
        {
            info.t = b;
            info.residual = std::abs(fb);
            info.converged = true;
            return info;
        }

        for (int iter = 0; iter < std::max(1, max_depth); ++iter)
        {
            const T m = T(0.5) * (a + b);
            const T fm = de_casteljau_eval_1d<T>(edge_coeffs, m);
            ++info.evaluations;
            info.iterations = iter + 1;
            info.t = m;
            info.residual = std::abs(fm);
            if (std::abs(fm) <= ftol || std::abs(b - a) <= xtol)
            {
                info.converged = std::abs(fm) <= ftol || std::abs(b - a) <= xtol;
                return info;
            }

            if ((fa < T(0) && fm > T(0)) || (fa > T(0) && fm < T(0)))
            {
                b = m;
                fb = fm;
            }
            else
            {
                a = m;
                fa = fm;
            }
        }
        return info;
    };

    if (bernstein_derivative_strict_one_sign(edge_coeffs, ftol)
        && ((f0 < T(0) && f1 > T(0)) || (f0 > T(0) && f1 < T(0))))
    {
        return bisect(T(0), T(1));
    }

    struct Candidate
    {
        std::vector<T> coeffs;
        T a = T(0);
        T b = T(1);
        int depth = 0;
    };

    Candidate current{std::vector<T>(edge_coeffs.begin(), edge_coeffs.end()), T(0), T(1), 0};
    while (current.depth < max_depth)
    {
        if (bernstein_derivative_strict_one_sign(
                std::span<const T>(current.coeffs.data(), current.coeffs.size()), ftol)
            && ((current.coeffs.front() < T(0) && current.coeffs.back() > T(0))
                || (current.coeffs.front() > T(0) && current.coeffs.back() < T(0))))
        {
            const auto local = bisect(current.a, current.b);
            out = local;
            return out;
        }

        std::vector<T> left;
        std::vector<T> right;
        de_casteljau_split_1d<T>(
            std::span<const T>(current.coeffs.data(), current.coeffs.size()), left, right, T(0.5));
        const bool left_candidate = !bernstein_all_strict_one_sign(
            std::span<const T>(left.data(), left.size()), ftol);
        const bool right_candidate = !bernstein_all_strict_one_sign(
            std::span<const T>(right.data(), right.size()), ftol);

        const T mid = T(0.5) * (current.a + current.b);
        if (left_candidate && !right_candidate)
        {
            current = Candidate{std::move(left), current.a, mid, current.depth + 1};
        }
        else if (!left_candidate && right_candidate)
        {
            current = Candidate{std::move(right), mid, current.b, current.depth + 1};
        }
        else
        {
            break;
        }
    }

    out.t = T(0.5);
    out.residual = std::abs(de_casteljau_eval_1d<T>(edge_coeffs, out.t));
    out.evaluations = 1;
    out.iterations = current.depth;
    out.converged = false;
    return out;
}

} // namespace cutcells
