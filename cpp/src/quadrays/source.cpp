// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "source.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cutcells::quadrays
{

namespace
{

/// Physical point of box coordinates u.
template <std::floating_point T>
std::array<double, 3> physical(const Source<T>& phi, const T* u)
{
    std::array<double, 3> x;
    for (int i = 0; i < 3; ++i)
    {
        double s = phi.origin[i];
        for (int k = 0; k < 3; ++k)
            s += static_cast<double>(phi.jacobian[i][k]) * u[k];
        x[i] = s;
    }
    return x;
}

/// Physical image of a box direction.
template <std::floating_point T>
std::array<double, 3> physical_direction(const Source<T>& phi, const T* v)
{
    std::array<double, 3> x;
    for (int i = 0; i < 3; ++i)
    {
        double s = 0;
        for (int k = 0; k < 3; ++k)
            s += static_cast<double>(phi.jacobian[i][k]) * v[k];
        x[i] = s;
    }
    return x;
}

/// Sign of a Taylor model row (alpha, beta_0..beta_{m-1}, eps), with the
/// tolerance of certain_sign.
int row_sign(const double* row, int m, double& alpha, double& deviation)
{
    alpha = row[0];
    deviation = row[m + 1];
    for (int j = 0; j < m; ++j)
        deviation += std::abs(row[1 + j]);
    const double x = (1.0 - 10.0 * std::numeric_limits<double>::epsilon()) * deviation;
    return alpha > x ? 1 : (alpha < -x ? -1 : 0);
}

template <std::floating_point T>
T value_on_line(const Source<T>& phi, const std::array<double, 3>& p, const std::array<double, 3>& v, T t)
{
    const double x[3] = {p[0] + t * v[0], p[1] + t * v[1], p[2] + t * v[2]};
    return static_cast<T>(phi.analytic->value(x, phi.analytic->context));
}

/// Intervals of one line that line_roots may visit; beyond it, every interval
/// left is treated as one where the slope has a sign.
constexpr int max_line_intervals = 4096;

template <std::floating_point T>
void isolate(const Source<T>& phi, const std::array<double, 3>& p, const std::array<double, 3>& v, T a, T b,
             double zero, int depth, int& visited, std::vector<T>& roots)
{
    const double c = 0.5 * (static_cast<double>(a) + b), h = 0.5 * (static_cast<double>(b) - a);
    const double centre[3] = {p[0] + c * v[0], p[1] + c * v[1], p[2] + c * v[2]};
    const double axis[3] = {h * v[0], h * v[1], h * v[2]};
    double models[2 * 3];
    int value_sign = 0, slope_sign = 0;
    double magnitude = std::numeric_limits<double>::infinity();
    ++visited;
    const int status = parallelepiped_bounds(*phi.analytic, centre, axis, 1, models);
    if (status != 0)
    {
        double alpha, dev;
        value_sign = row_sign(models, 1, alpha, dev);
        magnitude = std::abs(alpha) + dev;
        if (status == 1)
            slope_sign = row_sign(models + 3, 1, alpha, dev);
    }
    if (value_sign != 0)
        return;
    if (slope_sign != 0 || depth >= 40 || visited >= max_line_intervals)
    {
        const T ga = value_on_line(phi, p, v, a), gb = value_on_line(phi, p, v, b);
        if (ga != T(0) && gb != T(0) && (ga > T(0)) != (gb > T(0)))
            roots.push_back(illinois_root([&](T t) { return value_on_line(phi, p, v, t); }, a, b, ga, gb));
        return;
    }
    // phi vanishes on the interval up to the tolerance: the line lies in the
    // zero set (a plane in a face, say) or touches it; no breakpoint
    if (magnitude <= zero)
        return;
    const T m = T(0.5) * (a + b);
    isolate(phi, p, v, a, m, zero, depth + 1, visited, roots);
    isolate(phi, p, v, m, b, zero, depth + 1, visited, roots);
}

/// Sign of phi on the sub-box [lo, hi] of a cell: +1, -1, 0 if it may vanish,
/// or 2 if the sub-box misses the clipped region.
template <std::floating_point T>
int box_sign(const ClippedBox<T>& cell, const Source<T>& phi, const Vec3<T>& lo, const Vec3<T>& hi, int depth,
             int max_depth)
{
    if (!may_meet_clips(cell, lo, hi))
        return 2;
    double centre[3], axes[9];
    for (int i = 0; i < 3; ++i)
    {
        double c = cell.origin[i];
        for (int k = 0; k < 3; ++k)
        {
            c += static_cast<double>(cell.jacobian[i][k]) * 0.5 * (static_cast<double>(lo[k]) + hi[k]);
            axes[i * 3 + k] = static_cast<double>(cell.jacobian[i][k]) * 0.5 * (static_cast<double>(hi[k]) - lo[k]);
        }
        centre[i] = c;
    }
    if (phi.bernstein != nullptr)
    {
        thread_local BoxBernstein<T> sub;
        thread_local std::vector<T> work;
        subdivide(*phi.bernstein, std::span<const T>(lo), std::span<const T>(hi), sub, work);
        const int s = coefficient_sign(std::span<const T>(sub.coeffs));
        if (s != 0)
            return s;
    }
    else
    {
        double models[4 * 5];
        if (parallelepiped_bounds(*phi.analytic, centre, axes, 3, models) != 0)
        {
            double alpha, dev;
            const int s = row_sign(models, 3, alpha, dev);
            if (s != 0)
                return s;
        }
    }
    if (depth >= max_depth)
        return 0;
    // halve the longest edge
    int axis = 0;
    double longest = -1;
    for (int k = 0; k < 3; ++k)
    {
        double length = 0;
        for (int i = 0; i < 3; ++i)
            length += axes[i * 3 + k] * axes[i * 3 + k];
        if (length > longest)
        {
            longest = length;
            axis = k;
        }
    }
    Vec3<T> first_hi = hi, second_lo = lo;
    first_hi[axis] = second_lo[axis] = T(0.5) * (lo[axis] + hi[axis]);
    const int a = box_sign(cell, phi, lo, first_hi, depth + 1, max_depth);
    if (a == 0)
        return 0;
    const int b = box_sign(cell, phi, second_lo, hi, depth + 1, max_depth);
    if (b == 0)
        return 0;
    if (a == 2)
        return b;
    if (b == 2)
        return a;
    return a == b ? a : 0;
}

} // namespace

template <std::floating_point T>
Source<T> bernstein_source(const BoxBernstein<T>& phi)
{
    Source<T> s;
    s.bernstein = &phi;
    return s;
}

template <std::floating_point T>
Source<T> analytic_source(const AnalyticLevelSet& phi, const ClippedBox<T>& cell)
{
    if (phi.value == nullptr || phi.gradient == nullptr || phi.box_bounds == nullptr)
    {
        throw std::invalid_argument(
            "quadrays: an analytic level set needs value, gradient and box_bounds");
    }
    Source<T> s;
    s.analytic = &phi;
    s.origin = cell.origin;
    s.jacobian = cell.jacobian;
    return s;
}

template <std::floating_point T>
T evaluate(const Source<T>& phi, std::span<const T> u)
{
    if (phi.bernstein != nullptr)
        return evaluate(*phi.bernstein, u);
    const std::array<double, 3> x = physical(phi, u.data());
    return static_cast<T>(phi.analytic->value(x.data(), phi.analytic->context));
}

template <std::floating_point T>
void gradient(const Source<T>& phi, std::span<const T> u, std::span<T> g)
{
    if (phi.bernstein != nullptr)
    {
        gradient(*phi.bernstein, u, g);
        return;
    }
    // grad_u phi = J^T grad_x phi
    const std::array<double, 3> x = physical(phi, u.data());
    double gx[3];
    phi.analytic->gradient(x.data(), gx, phi.analytic->context);
    for (int k = 0; k < 3; ++k)
    {
        double s = 0;
        for (int i = 0; i < 3; ++i)
            s += gx[i] * phi.jacobian[i][k];
        g[k] = static_cast<T>(s);
    }
}

template <std::floating_point T>
T reference_magnitude(const Source<T>& phi)
{
    if (phi.bernstein != nullptr)
        return max_abs(std::span<const T>(phi.bernstein->coeffs));
    T m = T(0);
    for (int c = 0; c < 8; ++c)
    {
        const std::array<T, 3> u = {T(c & 1), T((c >> 1) & 1), T((c >> 2) & 1)};
        m = std::max(m, std::abs(evaluate(phi, std::span<const T>(u))));
    }
    return m;
}

template <std::floating_point T>
bool affine_bounds(const Source<T>& phi, std::span<const T> origin, std::span<const T> matrix, int m,
                   AffineBounds<T>& out)
{
    if (phi.analytic == nullptr)
        throw std::invalid_argument("quadrays: affine_bounds needs an analytic level set");
    // the parallelepiped centre + axes t, t in [-1, 1]^m, with s = (t + 1) / 2
    std::array<T, 3> centre_u;
    std::array<T, 9> half{};
    for (int i = 0; i < 3; ++i)
    {
        T c = origin[i];
        for (int j = 0; j < m; ++j)
        {
            half[i * m + j] = T(0.5) * matrix[i * m + j];
            c += half[i * m + j];
        }
        centre_u[i] = c;
    }
    const std::array<double, 3> centre = physical(phi, centre_u.data());
    double axes[9];
    for (int j = 0; j < m; ++j)
    {
        const T column[3] = {half[j], half[m + j], half[2 * m + j]};
        const std::array<double, 3> a = physical_direction(phi, column);
        for (int i = 0; i < 3; ++i)
            axes[i * m + j] = a[i];
    }
    double models[4 * 5];
    const int status = parallelepiped_bounds(*phi.analytic, centre.data(), axes, m, models);
    if (status == 0)
        return false;
    double alpha, dev;
    out.sign = row_sign(models, m, alpha, dev);
    out.magnitude = static_cast<T>(std::abs(alpha) + dev);
    out.has_derivatives = status == 1;
    out.lower.fill(T(0));
    out.upper.fill(T(0));
    for (int j = 0; j < m && out.has_derivatives; ++j)
    {
        // d / d s_j = 2 d / d t_j
        const int s = row_sign(models + (m + 2) * (1 + j), m, alpha, dev);
        out.upper[j] = static_cast<T>(2 * (std::abs(alpha) + dev));
        out.lower[j] = s != 0 ? static_cast<T>(2 * (std::abs(alpha) - dev)) : T(0);
    }
    return true;
}

template <std::floating_point T>
void line_roots(const Source<T>& phi, std::span<const T> origin, std::span<const T> direction, T a, T b, T zero,
                std::vector<T>& roots)
{
    if (phi.analytic == nullptr)
        throw std::invalid_argument("quadrays: line_roots needs an analytic level set");
    int visited = 0;
    isolate(phi, physical(phi, origin.data()), physical_direction(phi, direction.data()), a, b,
            static_cast<double>(zero), 0, visited, roots);
}

template <std::floating_point T>
T line_root(const Source<T>& phi, std::span<const T> origin, std::span<const T> direction, T a, T b, T ga,
            T gb)
{
    if (phi.analytic == nullptr)
        throw std::invalid_argument("quadrays: line_root needs an analytic level set");
    const std::array<double, 3> p = physical(phi, origin.data()), v = physical_direction(phi, direction.data());
    return illinois_root([&](T t) { return value_on_line(phi, p, v, t); }, a, b, ga, gb);
}

template <std::floating_point T>
int coefficient_sign(std::span<const T> coeffs)
{
    const T m = max_abs(coeffs);
    if (!(m > T(0)))
        return 0;
    const T tol = T(64) * std::numeric_limits<T>::epsilon() * m;
    bool positive = true, negative = true;
    for (const T c : coeffs)
    {
        positive &= c >= -tol;
        negative &= c <= tol;
    }
    return positive ? 1 : (negative ? -1 : 0);
}

template <std::floating_point T>
int cell_sign(const ClippedBox<T>& cell, const AnalyticLevelSet& phi, int max_depth)
{
    return cell_sign(cell, analytic_source(phi, cell), max_depth);
}

template <std::floating_point T>
int cell_sign(const ClippedBox<T>& cell, const Source<T>& phi, int max_depth)
{
    const Vec3<T> lo = {0, 0, 0}, hi = {1, 1, 1};
    // the whole box first: most cells are far from the zero set
    const int s = box_sign(cell, phi, lo, hi, 0, 0);
    if (s == 1 || s == -1)
        return s;
    // corners of the clipped region with both signs
    int seen = 0;
    for (int c = 0; c < 8; ++c)
    {
        const Vec3<T> u = {T(c & 1), T((c >> 1) & 1), T((c >> 2) & 1)};
        if (!inside_clips(cell, u, scaled_tolerance<T>(1e-12)))
            continue;
        const T v = evaluate(phi, std::span<const T>(u));
        seen |= (v > T(0) ? 1 : 0) | (v < T(0) ? 2 : 0);
    }
    if (seen == 3)
        return 0;
    const int r = box_sign(cell, phi, lo, hi, 0, max_depth);
    return r == 2 ? 0 : r;
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template Source<float> bernstein_source<float>(const BoxBernstein<float>&);
template Source<double> bernstein_source<double>(const BoxBernstein<double>&);
template Source<float> analytic_source<float>(const AnalyticLevelSet&, const ClippedBox<float>&);
template Source<double> analytic_source<double>(const AnalyticLevelSet&, const ClippedBox<double>&);
template float evaluate<float>(const Source<float>&, std::span<const float>);
template double evaluate<double>(const Source<double>&, std::span<const double>);
template void gradient<float>(const Source<float>&, std::span<const float>, std::span<float>);
template void gradient<double>(const Source<double>&, std::span<const double>, std::span<double>);
template float reference_magnitude<float>(const Source<float>&);
template double reference_magnitude<double>(const Source<double>&);
template bool affine_bounds<float>(const Source<float>&, std::span<const float>, std::span<const float>, int,
                                   AffineBounds<float>&);
template bool affine_bounds<double>(const Source<double>&, std::span<const double>, std::span<const double>, int,
                                    AffineBounds<double>&);
template void line_roots<float>(const Source<float>&, std::span<const float>, std::span<const float>, float, float,
                                float, std::vector<float>&);
template void line_roots<double>(const Source<double>&, std::span<const double>, std::span<const double>, double,
                                 double, double, std::vector<double>&);
template float line_root<float>(const Source<float>&, std::span<const float>, std::span<const float>, float, float,
                                float, float);
template double line_root<double>(const Source<double>&, std::span<const double>, std::span<const double>, double,
                                  double, double, double);
template int coefficient_sign<float>(std::span<const float>);
template int coefficient_sign<double>(std::span<const double>);
template int cell_sign<float>(const ClippedBox<float>&, const Source<float>&, int);
template int cell_sign<double>(const ClippedBox<double>&, const Source<double>&, int);
template int cell_sign<float>(const ClippedBox<float>&, const AnalyticLevelSet&, int);
template int cell_sign<double>(const ClippedBox<double>&, const AnalyticLevelSet&, int);

} // namespace cutcells::quadrays
