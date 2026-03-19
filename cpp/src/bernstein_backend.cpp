// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "bernstein_backend.h"

#include <map>

namespace cutcells
{
namespace
{

int check_order(int p)
{
    if (p < 0)
        throw std::invalid_argument("Bernstein degree must be non-negative");
    return p;
}

int invert_dense_square(std::vector<double>& a, int n)
{
    std::vector<double> inv(static_cast<std::size_t>(n * n), 0.0);
    for (int i = 0; i < n; ++i)
        inv[static_cast<std::size_t>(i * n + i)] = 1.0;

    for (int k = 0; k < n; ++k)
    {
        int piv = k;
        double maxv = std::abs(a[static_cast<std::size_t>(k * n + k)]);
        for (int r = k + 1; r < n; ++r)
        {
            const double v = std::abs(a[static_cast<std::size_t>(r * n + k)]);
            if (v > maxv)
            {
                maxv = v;
                piv = r;
            }
        }
        if (maxv <= std::numeric_limits<double>::epsilon())
            return -1;

        if (piv != k)
        {
            for (int c = 0; c < n; ++c)
            {
                std::swap(a[static_cast<std::size_t>(k * n + c)],
                          a[static_cast<std::size_t>(piv * n + c)]);
                std::swap(inv[static_cast<std::size_t>(k * n + c)],
                          inv[static_cast<std::size_t>(piv * n + c)]);
            }
        }

        const double diag = a[static_cast<std::size_t>(k * n + k)];
        for (int c = 0; c < n; ++c)
        {
            a[static_cast<std::size_t>(k * n + c)] /= diag;
            inv[static_cast<std::size_t>(k * n + c)] /= diag;
        }

        for (int r = 0; r < n; ++r)
        {
            if (r == k)
                continue;
            const double f = a[static_cast<std::size_t>(r * n + k)];
            if (std::abs(f) <= std::numeric_limits<double>::epsilon())
                continue;
            for (int c = 0; c < n; ++c)
            {
                a[static_cast<std::size_t>(r * n + c)] -= f * a[static_cast<std::size_t>(k * n + c)];
                inv[static_cast<std::size_t>(r * n + c)] -= f * inv[static_cast<std::size_t>(k * n + c)];
            }
        }
    }

    a = std::move(inv);
    return 0;
}

double evaluate_bernstein_basis_double(cell::type ct, int p, int basis_id, std::span<const double> x)
{
    switch (ct)
    {
    case cell::type::interval:
        return bernstein_basis_interval<double>(p, basis_id, x[0]);
    case cell::type::triangle:
        for (int i = 0; i <= p; ++i)
        {
            for (int j = 0; j <= p - i; ++j)
            {
                const int k = p - i - j;
                if (triangle_ijk_to_flat(p, i, j, k) == basis_id)
                    return bernstein_basis_triangle<double>(p, i, j, k, x[0], x[1]);
            }
        }
        break;
    case cell::type::quadrilateral:
        for (int j = 0; j <= p; ++j)
        {
            for (int i = 0; i <= p; ++i)
            {
                if (quad_ij_to_flat(p, i, j) == basis_id)
                    return bernstein_basis_quadrilateral<double>(p, i, j, x[0], x[1]);
            }
        }
        break;
    case cell::type::tetrahedron:
        for (int i = 0; i <= p; ++i)
        {
            for (int j = 0; j <= p - i; ++j)
            {
                for (int k = 0; k <= p - i - j; ++k)
                {
                    const int l = p - i - j - k;
                    if (tetrahedron_ijkl_to_flat(p, i, j, k, l) == basis_id)
                        return bernstein_basis_tetrahedron<double>(p, i, j, k, l, x[0], x[1], x[2]);
                }
            }
        }
        break;
    case cell::type::hexahedron:
        for (int k = 0; k <= p; ++k)
        {
            for (int j = 0; j <= p; ++j)
            {
                for (int i = 0; i <= p; ++i)
                {
                    if (hex_ijk_to_flat(p, i, j, k) == basis_id)
                        return bernstein_basis_hexahedron<double>(p, i, j, k, x[0], x[1], x[2]);
                }
            }
        }
        break;
    default:
        break;
    }

    throw std::invalid_argument("evaluate_bernstein_basis_double: unsupported basis id");
}

std::vector<double> build_inverse_matrix(cell::type ct, int order)
{
    const int n = bernstein_size(ct, order);
    const auto xref = bernstein_reference_nodes(ct, order);
    const int tdim = cell::get_tdim(ct);
    if (static_cast<int>(xref.size()) != n * tdim)
        throw std::invalid_argument("build_inverse_matrix: unexpected reference-node size");

    std::vector<double> a(static_cast<std::size_t>(n * n), 0.0);
    for (int row = 0; row < n; ++row)
    {
        const std::span<const double> x(
            xref.data() + static_cast<std::size_t>(row * tdim),
            static_cast<std::size_t>(tdim));
        for (int col = 0; col < n; ++col)
        {
            a[static_cast<std::size_t>(row * n + col)] = evaluate_bernstein_basis_double(ct, order, col, x);
        }
    }

    if (invert_dense_square(a, n) != 0)
        throw std::runtime_error("build_inverse_matrix: singular interpolation matrix");
    return a;
}

} // namespace

int interval_bernstein_size(int p)
{
    return check_order(p) + 1;
}

int triangle_bernstein_size(int p)
{
    check_order(p);
    return (p + 1) * (p + 2) / 2;
}

int tetrahedron_bernstein_size(int p)
{
    check_order(p);
    return (p + 1) * (p + 2) * (p + 3) / 6;
}

int quadrilateral_bernstein_size(int p)
{
    check_order(p);
    return (p + 1) * (p + 1);
}

int hexahedron_bernstein_size(int p)
{
    check_order(p);
    return (p + 1) * (p + 1) * (p + 1);
}

int bernstein_size(cell::type ct, int p)
{
    switch (ct)
    {
    case cell::type::interval:
        return interval_bernstein_size(p);
    case cell::type::triangle:
        return triangle_bernstein_size(p);
    case cell::type::quadrilateral:
        return quadrilateral_bernstein_size(p);
    case cell::type::tetrahedron:
        return tetrahedron_bernstein_size(p);
    case cell::type::hexahedron:
        return hexahedron_bernstein_size(p);
    case cell::type::prism:
    case cell::type::pyramid:
        throw std::invalid_argument("bernstein_size: prism and pyramid are not implemented yet");
    default:
        throw std::invalid_argument("bernstein_size: unsupported cell type");
    }
}

int lagrange_node_count(cell::type ct, int degree)
{
    if (degree < 1)
        throw std::invalid_argument("lagrange_node_count: degree must be >= 1");
    if (degree == 1)
        return cell::get_num_vertices(ct);
    return iso_p1_template(ct, degree).n_vertices;
}

int interval_i_to_flat(int p, int i)
{
    if (i < 0 || i > p)
        throw std::invalid_argument("interval_i_to_flat: index out of range");
    return i;
}

int triangle_ijk_to_flat(int p, int i, int j, int k)
{
    if (i < 0 || j < 0 || k < 0 || i + j + k != p)
        throw std::invalid_argument("triangle_ijk_to_flat: invalid barycentric index");

    int flat = 0;
    for (int ii = 0; ii < i; ++ii)
    {
        for (int jj = 0; jj <= p - ii; ++jj)
            ++flat;
    }
    for (int jj = 0; jj < j; ++jj)
        ++flat;
    return flat;
}

int tetrahedron_ijkl_to_flat(int p, int i, int j, int k, int l)
{
    if (i < 0 || j < 0 || k < 0 || l < 0 || i + j + k + l != p)
        throw std::invalid_argument("tetrahedron_ijkl_to_flat: invalid barycentric index");

    int flat = 0;
    for (int ii = 0; ii < i; ++ii)
    {
        for (int jj = 0; jj <= p - ii; ++jj)
        {
            for (int kk = 0; kk <= p - ii - jj; ++kk)
                ++flat;
        }
    }
    for (int jj = 0; jj < j; ++jj)
    {
        for (int kk = 0; kk <= p - i - jj; ++kk)
            ++flat;
    }
    for (int kk = 0; kk < k; ++kk)
        ++flat;
    return flat;
}

int quad_ij_to_flat(int p, int i, int j)
{
    if (i < 0 || i > p || j < 0 || j > p)
        throw std::invalid_argument("quad_ij_to_flat: index out of range");
    return j * (p + 1) + i;
}

int hex_ijk_to_flat(int p, int i, int j, int k)
{
    if (i < 0 || i > p || j < 0 || j > p || k < 0 || k > p)
        throw std::invalid_argument("hex_ijk_to_flat: index out of range");
    return (k * (p + 1) + j) * (p + 1) + i;
}

int infer_lagrange_order_from_num_nodes(cell::type ct, int n_vertices)
{
    if (n_vertices == cell::get_num_vertices(ct))
        return 1;

    for (int p = 2; p <= 4; ++p)
    {
        if (iso_p1_template(ct, p).n_vertices == n_vertices)
            return p;
    }

    throw std::invalid_argument("infer_lagrange_order_from_num_nodes: unsupported node count");
}

const std::vector<double>& lagrange_to_bernstein_inverse_matrix(cell::type ct, int order)
{
    if (order < 1 || order > 4)
        throw std::invalid_argument("lagrange_to_bernstein_inverse_matrix: supported orders are 1..4");

    using Key = std::pair<cell::type, int>;
    static std::map<Key, std::vector<double>> cache;

    const Key key{ct, order};
    const auto it = cache.find(key);
    if (it != cache.end())
        return it->second;

    auto [inserted_it, _] = cache.emplace(key, build_inverse_matrix(ct, order));
    return inserted_it->second;
}

} // namespace cutcells
