// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "piece_rules.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "../quadrature_tables.h"
#include "triangulation.h"

namespace cutcells::lut
{

namespace
{

bool is_simplex(cell::type type)
{
    return type == cell::type::interval || type == cell::type::triangle || type == cell::type::tetrahedron;
}

/// The point x(X) of a straight cell given by its vertices (Basix order, dim
/// coordinates each) and the derivative dx/dX (dim x tdim of the cell,
/// row-major); multilinear on quadrilaterals and hexahedra, the triangle's map
/// times the interval's on prisms, Basix's rational P1 map on pyramids (affine
/// where the base is a parallelogram).
template <std::floating_point T>
void map_point(cell::type type, const T* vertices, int dim, const T* X, T* x, T* jacobian)
{
    const int m = cell::get_tdim(type);
    if (type == cell::type::prism)
    {
        const T l[3] = {T(1) - X[0] - X[1], X[0], X[1]};
        const T dl[3][2] = {{-1, -1}, {1, 0}, {0, 1}};
        for (int i = 0; i < dim; ++i)
        {
            x[i] = T(0);
            for (int j = 0; j < 3; ++j)
                jacobian[i * 3 + j] = T(0);
            for (int a = 0; a < 3; ++a)
            {
                const T bottom = vertices[a * dim + i], top = vertices[(a + 3) * dim + i];
                x[i] += l[a] * ((T(1) - X[2]) * bottom + X[2] * top);
                for (int j = 0; j < 2; ++j)
                    jacobian[i * 3 + j] += dl[a][j] * ((T(1) - X[2]) * bottom + X[2] * top);
                jacobian[i * 3 + 2] += l[a] * (top - bottom);
            }
        }
        return;
    }
    if (type == cell::type::pyramid)
    {
        // x = v0 (1 - X - Y - Z) + v1 X + v2 Y + v4 Z + (v0 - v1 - v2 + v3) X Y / (1 - Z)
        const T r = std::max(T(1) - X[2], std::numeric_limits<T>::epsilon());
        for (int i = 0; i < dim; ++i)
        {
            const T v0 = vertices[i], v1 = vertices[dim + i], v2 = vertices[2 * dim + i], v3 = vertices[3 * dim + i],
                    v4 = vertices[4 * dim + i];
            const T w = v0 - v1 - v2 + v3;
            x[i] = v0 * (T(1) - X[0] - X[1] - X[2]) + v1 * X[0] + v2 * X[1] + v4 * X[2] + w * X[0] * X[1] / r;
            jacobian[i * 3] = v1 - v0 + w * X[1] / r;
            jacobian[i * 3 + 1] = v2 - v0 + w * X[0] / r;
            jacobian[i * 3 + 2] = v4 - v0 + w * X[0] * X[1] / (r * r);
        }
        return;
    }
    if (m == 0 || is_simplex(type))
    {
        for (int i = 0; i < dim; ++i)
        {
            x[i] = vertices[i];
            for (int j = 0; j < m; ++j)
            {
                const T e = vertices[(j + 1) * dim + i] - vertices[i];
                x[i] += X[j] * e;
                jacobian[i * m + j] = e;
            }
        }
        return;
    }
    if (type != cell::type::quadrilateral && type != cell::type::hexahedron)
        throw std::invalid_argument("lut: no map for this cell type");
    // vertex v of the unit square or cube has coordinate d = bit d of v
    std::array<T, 8> n{};
    std::array<T, 24> dn{};
    for (int v = 0; v < (1 << m); ++v)
    {
        n[v] = T(1);
        for (int d = 0; d < m; ++d)
        {
            n[v] *= ((v >> d) & 1) ? X[d] : T(1) - X[d];
            T g = ((v >> d) & 1) ? T(1) : T(-1);
            for (int e = 0; e < m; ++e)
                if (e != d)
                    g *= ((v >> e) & 1) ? X[e] : T(1) - X[e];
            dn[v * m + d] = g;
        }
    }
    for (int i = 0; i < dim; ++i)
    {
        x[i] = T(0);
        for (int j = 0; j < m; ++j)
            jacobian[i * m + j] = T(0);
        for (int v = 0; v < (1 << m); ++v)
        {
            const T c = vertices[v * dim + i];
            x[i] += n[v] * c;
            for (int j = 0; j < m; ++j)
                jacobian[i * m + j] += dn[v * m + j] * c;
        }
    }
}

template <std::floating_point T>
T determinant(const T* a, int n)
{
    if (n == 1)
        return a[0];
    if (n == 2)
        return a[0] * a[3] - a[1] * a[2];
    return a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6])
           + a[2] * (a[3] * a[7] - a[4] * a[6]);
}

/// The m-dimensional measure of the columns of b (rows x m, row-major): a
/// determinant, a length, or for two columns in 3D the norm of their cross
/// product (a Gram determinant would lose half the digits on slivers).
template <std::floating_point T>
T measure(const T* b, int rows, int m)
{
    if (m == 0)
        return T(1);
    if (rows == m)
        return std::abs(determinant(b, m));
    if (m == 1)
    {
        T s = T(0);
        for (int k = 0; k < rows; ++k)
            s += b[k] * b[k];
        return std::sqrt(s);
    }
    // two columns in 3D
    const T c0 = b[2] * b[5] - b[4] * b[3], c1 = b[4] * b[1] - b[0] * b[5], c2 = b[0] * b[3] - b[2] * b[1];
    return std::sqrt(c0 * c0 + c1 * c1 + c2 * c2);
}

/// Reference rules by cell type, for one degree.
template <std::floating_point T>
const quadrature::ReferenceQuadratureRule<T>& reference_rule(cell::type type, int degree)
{
    thread_local std::array<std::array<quadrature::ReferenceQuadratureRule<T>, 11>, 8> cache;
    if (degree < 1 || degree > 10)
        throw std::invalid_argument("lut: rules go from degree 1 to 10");
    quadrature::ReferenceQuadratureRule<T>& rule = cache[static_cast<int>(type)][degree];
    if (rule._weights.empty())
        rule = quadrature::get_reference_rule<T>(type, degree);
    return rule;
}

} // namespace

template <std::floating_point T>
bool is_parallelotope(std::span<const T> x, int tdim, int dim)
{
    T size = T(0), gap = T(0);
    for (int v = 1; v < (1 << tdim); ++v)
        for (int i = 0; i < dim; ++i)
        {
            T p = x[static_cast<std::size_t>(i)];
            for (int d = 0; d < tdim; ++d)
                if ((v >> d) & 1)
                    p += x[static_cast<std::size_t>((1 << d) * dim + i)] - x[static_cast<std::size_t>(i)];
            const T xv = x[static_cast<std::size_t>(v * dim + i)];
            gap = std::max(gap, std::abs(xv - p));
            size = std::max(size, std::abs(xv - x[static_cast<std::size_t>(i)]));
        }
    return gap <= T(64) * std::numeric_limits<T>::epsilon() * size;
}

template <std::floating_point T>
void push_forward(const CellMap<T>& map, std::span<const T> xi, std::vector<T>& x)
{
    const int tdim = map.tdim(), gdim = map.gdim;
    const std::size_t n = xi.size() / static_cast<std::size_t>(tdim);
    x.resize(n * static_cast<std::size_t>(gdim));
    std::array<T, 9> jacobian;
    for (std::size_t p = 0; p < n; ++p)
        map_point(map.type, map.vertices.data(), gdim, xi.data() + p * tdim, x.data() + p * gdim, jacobian.data());
}

template <std::floating_point T>
void append_piece_rule(const CellMap<T>& map, cell::type piece_type, std::span<const T> piece_vertices, int degree,
                       std::vector<T>& points, std::vector<T>& weights)
{
    const int tdim = map.tdim(), gdim = map.gdim;
    const int m = cell::get_tdim(piece_type);
    const bool tensor = piece_type == cell::type::quadrilateral || piece_type == cell::type::hexahedron;
    if (piece_type == cell::type::prism || piece_type == cell::type::pyramid
        || (tensor && !is_parallelotope(piece_vertices, m, tdim)))
    {
        // simplices keep the rule exact on pieces with planar faces
        int ids[8] = {0, 1, 2, 3, 4, 5, 6, 7};
        std::vector<std::vector<int>> simplices;
        cell::triangulation(piece_type, ids, simplices);
        const cell::type simplex = m == 3 ? cell::type::tetrahedron : cell::type::triangle;
        std::array<T, 12> s;
        for (const std::vector<int>& t : simplices)
        {
            for (int j = 0; j <= m; ++j)
                for (int d = 0; d < tdim; ++d)
                    s[j * tdim + d] = piece_vertices[static_cast<std::size_t>(t[j] * tdim + d)];
            append_piece_rule(map, simplex, std::span<const T>(s.data(), (m + 1) * tdim), degree, points, weights);
        }
        return;
    }
    if (m == 0)
    {
        points.insert(points.end(), piece_vertices.begin(), piece_vertices.begin() + tdim);
        weights.push_back(T(1));
        return;
    }
    const quadrature::ReferenceQuadratureRule<T>& rule = reference_rule<T>(piece_type, degree);
    std::array<T, 3> xi, x;
    std::array<T, 9> a, j, b;
    for (int q = 0; q < rule._num_points; ++q)
    {
        // the piece in the reference cell, then the cell in physical space
        map_point(piece_type, piece_vertices.data(), tdim, rule._points.data() + q * m, xi.data(), a.data());
        map_point(map.type, map.vertices.data(), gdim, xi.data(), x.data(), j.data());
        for (int r = 0; r < gdim; ++r)
            for (int c = 0; c < m; ++c)
            {
                T s = T(0);
                for (int k = 0; k < tdim; ++k)
                    s += j[r * tdim + k] * a[k * m + c];
                b[r * m + c] = s;
            }
        points.insert(points.end(), xi.begin(), xi.begin() + tdim);
        weights.push_back(rule._weights[q] * measure(b.data(), gdim, m));
    }
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template bool is_parallelotope<float>(std::span<const float>, int, int);
template bool is_parallelotope<double>(std::span<const double>, int, int);
template void push_forward<float>(const CellMap<float>&, std::span<const float>, std::vector<float>&);
template void push_forward<double>(const CellMap<double>&, std::span<const double>, std::vector<double>&);
template void append_piece_rule<float>(const CellMap<float>&, cell::type, std::span<const float>, int,
                                       std::vector<float>&, std::vector<float>&);
template void append_piece_rule<double>(const CellMap<double>&, cell::type, std::span<const double>, int,
                                        std::vector<double>&, std::vector<double>&);

} // namespace cutcells::lut
