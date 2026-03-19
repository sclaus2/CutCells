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

template <std::floating_point T>
inline void restrict_bernstein_to_segment_1d(
    const BernsteinCell<T>& cell_poly,
    std::span<const T> x0_ref,
    std::span<const T> x1_ref,
    std::vector<T>& edge_coeffs)
{
    const int p = cell_poly.degree;
    if (static_cast<int>(x0_ref.size()) != cell_poly.tdim || static_cast<int>(x1_ref.size()) != cell_poly.tdim)
        throw std::invalid_argument("restrict_bernstein_to_segment_1d: segment dimension mismatch");

    std::vector<T> samples(static_cast<std::size_t>(p + 1), T(0));
    std::vector<T> x(static_cast<std::size_t>(cell_poly.tdim), T(0));
    for (int j = 0; j <= p; ++j)
    {
        const T t = (p == 0) ? T(0) : static_cast<T>(j) / static_cast<T>(p);
        for (int d = 0; d < cell_poly.tdim; ++d)
            x[static_cast<std::size_t>(d)] = (T(1) - t) * x0_ref[static_cast<std::size_t>(d)]
                                             + t * x1_ref[static_cast<std::size_t>(d)];
        samples[static_cast<std::size_t>(j)] = evaluate_bernstein_cell<T>(
            cell_poly, std::span<const T>(x.data(), x.size()));
    }

    edge_coeffs = bernstein_coefficients_from_equispaced_nodal_1d<T>(
        std::span<const T>(samples.data(), samples.size()));
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
