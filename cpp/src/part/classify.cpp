// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "classify.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>

#include "../quadrays/analytic.h"

namespace cutcells::part
{

namespace
{

/// Positions of the multi-indices (n - i - j - k, i, j, k) in CutCells'
/// simplex order of degree n, at i + (n + 1) (j + (n + 1) k).
const std::vector<int>& simplex_positions(int n)
{
    thread_local std::vector<std::vector<int>> cache;
    if (static_cast<int>(cache.size()) <= n)
        cache.resize(n + 1);
    std::vector<int>& pos = cache[n];
    if (pos.empty())
    {
        pos.assign(static_cast<std::size_t>((n + 1) * (n + 1) * (n + 1)), -1);
        int index = 0;
        for (int k = 0; k <= n; ++k)
            for (int j = 0; j <= n - k; ++j)
                for (int i = 0; i <= n - k - j; ++i)
                    pos[static_cast<std::size_t>(i + (n + 1) * (j + (n + 1) * k))] = index++;
    }
    return pos;
}

template <std::floating_point T>
using Vertices = std::array<quadrays::Vec3<T>, 4>;

/// The edge of a piece with the longest physical length.
template <std::floating_point T>
std::pair<int, int> longest_edge(const quadrays::ClippedBox<T>& box, const Vertices<T>& v)
{
    std::pair<int, int> edge = {0, 1};
    T longest = -1;
    for (int a = 0; a < 4; ++a)
        for (int b = a + 1; b < 4; ++b)
        {
            T length = 0;
            for (int i = 0; i < 3; ++i)
            {
                T d = 0;
                for (int k = 0; k < 3; ++k)
                    d += box.jacobian[i][k] * (v[a][k] - v[b][k]);
                length += d * d;
            }
            if (length > longest)
            {
                longest = length;
                edge = {a, b};
            }
        }
    return edge;
}

template <std::floating_point T>
std::pair<Vertices<T>, Vertices<T>> halves(const Vertices<T>& v, int p, int q)
{
    quadrays::Vec3<T> mid;
    for (int i = 0; i < 3; ++i)
        mid[i] = T(0.5) * (v[p][i] + v[q][i]);
    Vertices<T> keep_p = v, keep_q = v;
    keep_p[q] = mid;
    keep_q[p] = mid;
    return {keep_p, keep_q};
}

/// Sign of a polynomial on a piece of a tetrahedron from its Bernstein
/// coefficients on the piece, bisecting the longest edge.
template <std::floating_point T>
int bernstein_tet_sign(const quadrays::ClippedBox<T>& box, std::span<const T> c, int n, const Vertices<T>& v,
                       int depth, int max_depth, std::vector<std::vector<T>>& scratch)
{
    const int s = quadrays::coefficient_sign(c);
    if (s != 0 || depth >= max_depth)
        return s;
    const auto [p, q] = longest_edge(box, v);
    std::vector<T>& a = scratch[static_cast<std::size_t>(2 * depth)];
    std::vector<T>& b = scratch[static_cast<std::size_t>(2 * depth + 1)];
    bisect_tetrahedron(c, n, p, q, a, b);
    const auto [va, vb] = halves(v, p, q);
    const int sa = bernstein_tet_sign(box, std::span<const T>(a), n, va, depth + 1, max_depth, scratch);
    if (sa == 0)
        return 0;
    const int sb = bernstein_tet_sign(box, std::span<const T>(b), n, vb, depth + 1, max_depth, scratch);
    return sa == sb ? sa : 0;
}

/// Sign of an analytic level set on a piece of a tetrahedron: Taylor models
/// over the parallelepiped of the piece's edges from its vertex 0, whose
/// linear part takes its extremes over the piece at the piece's vertices, or
/// without Taylor models, intervals over the piece's bounding box.
template <std::floating_point T>
int analytic_tet_sign(const quadrays::ClippedBox<T>& box, const quadrays::AnalyticLevelSet& phi, const Vertices<T>& v,
                      int depth, int max_depth)
{
    if (phi.taylor_bounds == nullptr)
    {
        double lo[3], hi[3], b[8];
        for (int i = 0; i < 3; ++i)
        {
            lo[i] = std::numeric_limits<double>::infinity();
            hi[i] = -lo[i];
            for (const quadrays::Vec3<T>& u : v)
            {
                double x = box.origin[i];
                for (int k = 0; k < 3; ++k)
                    x += static_cast<double>(box.jacobian[i][k]) * u[k];
                lo[i] = std::min(lo[i], x);
                hi[i] = std::max(hi[i], x);
            }
        }
        if (phi.box_bounds(lo, hi, b, phi.context) != 0)
        {
            const double tol = 10 * std::numeric_limits<double>::epsilon() * 0.5 * (b[1] - b[0]);
            if (b[0] > -tol)
                return 1;
            if (b[1] < tol)
                return -1;
        }
        if (depth >= max_depth)
            return 0;
        const auto [p, q] = longest_edge(box, v);
        const auto [va, vb] = halves(v, p, q);
        const int sa = analytic_tet_sign(box, phi, va, depth + 1, max_depth);
        if (sa == 0)
            return 0;
        const int sb = analytic_tet_sign(box, phi, vb, depth + 1, max_depth);
        return sa == sb ? sa : 0;
    }
    double centre[3], axes[9];
    for (int i = 0; i < 3; ++i)
    {
        double c = box.origin[i];
        for (int k = 0; k < 3; ++k)
        {
            double mid = v[0][k];
            for (int j = 1; j < 4; ++j)
                mid += 0.5 * (static_cast<double>(v[j][k]) - v[0][k]);
            c += static_cast<double>(box.jacobian[i][k]) * mid;
        }
        centre[i] = c;
        for (int j = 0; j < 3; ++j)
        {
            double a = 0;
            for (int k = 0; k < 3; ++k)
                a += static_cast<double>(box.jacobian[i][k]) * 0.5 * (static_cast<double>(v[j + 1][k]) - v[0][k]);
            axes[i * 3 + j] = a;
        }
    }
    double models[4 * 5];
    if (quadrays::parallelepiped_bounds(phi, centre, axes, 3, models) != 0)
    {
        // the piece's vertices at t = (-1, -1, -1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)
        const double alpha = models[0], eps = models[4];
        const double base = alpha - models[1] - models[2] - models[3];
        double lo = base, hi = base, deviation = eps;
        for (int j = 0; j < 3; ++j)
        {
            const double value = base + 2 * models[1 + j];
            lo = std::min(lo, value);
            hi = std::max(hi, value);
            deviation += std::abs(models[1 + j]);
        }
        // a few rounding errors of tolerance, as in quadrays::certain_sign
        const double tol = 10 * std::numeric_limits<double>::epsilon() * deviation;
        if (lo - eps > -tol)
            return 1;
        if (hi + eps < tol)
            return -1;
    }
    if (depth >= max_depth)
        return 0;
    const auto [p, q] = longest_edge(box, v);
    const auto [va, vb] = halves(v, p, q);
    const int sa = analytic_tet_sign(box, phi, va, depth + 1, max_depth);
    if (sa == 0)
        return 0;
    const int sb = analytic_tet_sign(box, phi, vb, depth + 1, max_depth);
    return sa == sb ? sa : 0;
}

} // namespace

template <std::floating_point T>
void bisect_tetrahedron(std::span<const T> coeffs, int n, int p, int q, std::vector<T>& keep_p,
                        std::vector<T>& keep_q)
{
    const std::vector<int>& pos = simplex_positions(n);
    auto position = [&](const std::array<int, 4>& a)
    { return pos[static_cast<std::size_t>(a[1] + (n + 1) * (a[2] + (n + 1) * a[3]))]; };
    keep_p.assign(coeffs.size(), T(0));
    keep_q.assign(coeffs.size(), T(0));
    std::vector<T> w(static_cast<std::size_t>(n + 1)), left(w.size()), right(w.size());
    for (int k = 0; k <= n; ++k)
        for (int j = 0; j <= n - k; ++j)
            for (int i = 0; i <= n - k - j; ++i)
            {
                // each multi-index without a power of vertex q starts one
                // sequence along the edge from vertex p to vertex q
                const std::array<int, 4> a = {n - i - j - k, i, j, k};
                if (a[q] != 0)
                    continue;
                const int m = a[p];
                std::array<int, 4> b = a;
                for (int s = 0; s <= m; ++s)
                {
                    b[p] = m - s;
                    b[q] = s;
                    w[s] = coeffs[static_cast<std::size_t>(position(b))];
                }
                left[0] = w[0];
                right[m] = w[m];
                for (int r = 1; r <= m; ++r)
                {
                    for (int t = 0; t + r <= m; ++t)
                        w[t] = T(0.5) * (w[t] + w[t + 1]);
                    left[r] = w[0];
                    right[m - r] = w[m - r];
                }
                for (int s = 0; s <= m; ++s)
                {
                    b[p] = m - s;
                    b[q] = s;
                    keep_p[static_cast<std::size_t>(position(b))] = left[s];
                    keep_q[static_cast<std::size_t>(position(b))] = right[s];
                }
            }
}

template <std::floating_point T, std::integral I>
int cell_sign(const CellSource<T, I>& cs, cell::type type, int max_depth)
{
    if (type == cell::type::hexahedron)
        return quadrays::cell_sign(cs.box, cs.source, max_depth);
    const Vertices<T> v = {quadrays::Vec3<T>{0, 0, 0}, quadrays::Vec3<T>{1, 0, 0}, quadrays::Vec3<T>{0, 1, 0},
                           quadrays::Vec3<T>{0, 0, 1}};
    if (cs.source.bernstein != nullptr)
    {
        const std::span<const T> c(cs.ls_cell.bernstein_coeffs);
        const int n = cs.ls_cell.bernstein_order;
        // the coefficients at the vertices are the values there
        const std::vector<int>& pos = simplex_positions(n);
        const T corners[4] = {c[static_cast<std::size_t>(pos[0])], c[static_cast<std::size_t>(pos[n])],
                              c[static_cast<std::size_t>(pos[(n + 1) * n])],
                              c[static_cast<std::size_t>(pos[(n + 1) * (n + 1) * n])]};
        const int whole = quadrays::coefficient_sign(c);
        if (whole != 0)
            return whole;
        const bool positive = std::any_of(corners, corners + 4, [](T x) { return x > T(0); });
        const bool negative = std::any_of(corners, corners + 4, [](T x) { return x < T(0); });
        if (positive && negative)
            return 0;
        thread_local std::vector<std::vector<T>> scratch;
        scratch.resize(static_cast<std::size_t>(2 * (max_depth + 1)));
        return bernstein_tet_sign(cs.box, c, n, v, 0, max_depth, scratch);
    }
    const int whole = analytic_tet_sign(cs.box, *cs.source.analytic, v, 0, 0);
    if (whole != 0)
        return whole;
    bool positive = false, negative = false;
    for (const quadrays::Vec3<T>& u : v)
    {
        const T value = quadrays::evaluate(cs.source, std::span<const T>(u));
        positive |= value > T(0);
        negative |= value < T(0);
    }
    if (positive && negative)
        return 0;
    return analytic_tet_sign(cs.box, *cs.source.analytic, v, 0, max_depth);
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template void bisect_tetrahedron<float>(std::span<const float>, int, int, int, std::vector<float>&,
                                        std::vector<float>&);
template void bisect_tetrahedron<double>(std::span<const double>, int, int, int, std::vector<double>&,
                                         std::vector<double>&);
template int cell_sign<float, int>(const CellSource<float, int>&, cell::type, int);
template int cell_sign<double, int>(const CellSource<double, int>&, cell::type, int);

} // namespace cutcells::part
