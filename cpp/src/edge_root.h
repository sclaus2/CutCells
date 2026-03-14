// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cutcells::cell::edge_root
{
enum class method
{
    linear = 0,
    brent = 1,
    itp = 2,
    newton = 3
};

template <std::floating_point T>
struct RootSolveInfo
{
    T t = T(0.5);
    T residual = std::numeric_limits<T>::max();
    int iterations = 0;
    int evaluations = 0;
    bool converged = false;
};

template <std::floating_point T>
inline T clamp01(T t)
{
    return std::clamp(t, T(0), T(1));
}

template <std::floating_point T>
inline T linear_root_parameter(const T v0, const T v1, const T level = T(0))
{
    const T denom = v1 - v0;
    if (std::abs(denom) <= std::numeric_limits<T>::epsilon())
        return T(0.5);
    return clamp01<T>((level - v0) / denom);
}

template <std::floating_point T>
inline void interpolate_point(const std::span<const T> p0, const std::span<const T> p1,
                              const T t, const std::span<T> out)
{
    if (p0.size() != p1.size() || p0.size() != out.size())
        throw std::invalid_argument("interpolate_point: inconsistent point dimensions");

    for (std::size_t i = 0; i < p0.size(); ++i)
        out[i] = p0[i] + (p1[i] - p0[i]) * t;
}

template <std::floating_point T>
inline void linear_intersection_point(const std::span<const T> p0, const std::span<const T> p1,
                                      const T v0, const T v1, const std::span<T> out,
                                      const T level = T(0))
{
    const T t = linear_root_parameter(v0, v1, level);
    interpolate_point<T>(p0, p1, t, out);
}

template <std::floating_point T>
inline void linear_intersection_point(const std::span<const T> p0, const std::span<const T> p1,
                                      const T v0, const T v1, std::vector<T>& out,
                                      const int offset = 0, const T level = T(0))
{
    if (offset < 0)
        throw std::invalid_argument("linear_intersection_point: negative offset");
    if (out.size() < static_cast<std::size_t>(offset) + p0.size())
        throw std::invalid_argument("linear_intersection_point: output vector too small");

    linear_intersection_point<T>(p0, p1, v0, v1,
                                 std::span<T>(out.data() + offset, p0.size()), level);
}

template <std::floating_point T, typename Phi>
inline T evaluate_edge_function(Phi& phi, const std::span<const T> p0, const std::span<const T> p1,
                                const T t)
{
    std::array<T, 3> x = {T(0), T(0), T(0)};
    const std::size_t gdim = p0.size();
    if (gdim > x.size())
        throw std::invalid_argument("evaluate_edge_function: only gdim <= 3 supported");

    for (std::size_t j = 0; j < gdim; ++j)
        x[j] = p0[j] + (p1[j] - p0[j]) * t;

    return static_cast<T>(std::invoke(phi, std::span<const T>(x.data(), gdim)));
}

template <std::floating_point T, typename Eval>
inline T brent_parameter(Eval&& eval, T a, T b, T fa, T fb,
                         const int max_iter, const T xtol, const T ftol,
                         int* out_iterations = nullptr,
                         bool* out_converged = nullptr)
{
    if (out_iterations)
        *out_iterations = 0;
    if (out_converged)
        *out_converged = false;

    if (std::abs(fa) <= ftol)
    {
        if (out_converged)
            *out_converged = true;
        return a;
    }
    if (std::abs(fb) <= ftol)
    {
        if (out_converged)
            *out_converged = true;
        return b;
    }

    if (fa * fb > T(0))
        return clamp01<T>(T(0.5) * (a + b));

    T c = a;
    T fc = fa;
    T d = b - a;
    T e = d;

    for (int iter = 0; iter < max_iter; ++iter)
    {
        if ((fb > T(0) && fc > T(0)) || (fb < T(0) && fc < T(0)))
        {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }

        if (std::abs(fc) < std::abs(fb))
        {
            const T a_old = a;
            const T fa_old = fa;
            a = b;
            fa = fb;
            b = c;
            fb = fc;
            c = a_old;
            fc = fa_old;
        }

        const T eps = std::numeric_limits<T>::epsilon();
        const T tol1 = T(2) * eps * std::abs(b) + T(0.5) * xtol;
        const T xm = T(0.5) * (c - b);

        if (std::abs(xm) <= tol1 || std::abs(fb) <= ftol)
        {
            if (out_iterations)
                *out_iterations = iter + 1;
            if (out_converged)
                *out_converged = true;
            return clamp01<T>(b);
        }

        if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb))
        {
            T s = fb / fa;
            T p = T(0);
            T q = T(0);

            if (a == c)
            {
                p = T(2) * xm * s;
                q = T(1) - s;
            }
            else
            {
                const T q0 = fa / fc;
                const T r = fb / fc;
                p = s * (T(2) * xm * q0 * (q0 - r) - (b - a) * (r - T(1)));
                q = (q0 - T(1)) * (r - T(1)) * (s - T(1));
            }

            if (p > T(0))
                q = -q;
            p = std::abs(p);

            const T min1 = T(3) * xm * q - std::abs(tol1 * q);
            const T min2 = std::abs(e * q);
            if (T(2) * p < std::min(min1, min2))
            {
                e = d;
                d = p / q;
            }
            else
            {
                d = xm;
                e = d;
            }
        }
        else
        {
            d = xm;
            e = d;
        }

        a = b;
        fa = fb;

        if (std::abs(d) > tol1)
            b += d;
        else
            b += (xm > T(0) ? tol1 : -tol1);

        fb = eval(b);
    }

    if (out_iterations)
        *out_iterations = max_iter;
    return clamp01<T>(b);
}

template <std::floating_point T, typename Eval>
inline T itp_parameter(Eval&& eval, T a, T b, T fa, T fb, const T initial_guess,
                       const int max_iter, const T xtol, const T ftol,
                       int* out_iterations = nullptr,
                       bool* out_converged = nullptr)
{
    if (out_iterations)
        *out_iterations = 0;
    if (out_converged)
        *out_converged = false;

    if (std::abs(fa) <= ftol)
    {
        if (out_converged)
            *out_converged = true;
        return a;
    }
    if (std::abs(fb) <= ftol)
    {
        if (out_converged)
            *out_converged = true;
        return b;
    }

    if (fa * fb > T(0))
        return clamp01<T>(initial_guess);

    // Use the linear guess first to contract the bracket if possible.
    const T t0 = std::clamp(initial_guess, a, b);
    if (t0 > a && t0 < b)
    {
        const T f0 = eval(t0);
        if (std::abs(f0) <= ftol)
        {
            if (out_converged)
                *out_converged = true;
            return t0;
        }
        if (fa * f0 <= T(0))
        {
            b = t0;
            fb = f0;
        }
        else if (f0 * fb <= T(0))
        {
            a = t0;
            fa = f0;
        }
    }

    const T width0 = b - a;
    const int n_max = std::max(0, static_cast<int>(std::ceil(std::log2(width0 / xtol))));
    const T k1 = T(0.2) / std::max(width0, std::numeric_limits<T>::epsilon());
    const T k2 = T(2);
    const int iter_max = std::min(max_iter, n_max + 8);

    for (int iter = 0; iter < iter_max; ++iter)
    {
        const T width = b - a;
        if (width <= xtol)
        {
            if (out_iterations)
                *out_iterations = iter + 1;
            if (out_converged)
                *out_converged = true;
            return clamp01<T>(T(0.5) * (a + b));
        }

        const T m = T(0.5) * (a + b);
        const T xf = (a * fb - b * fa) / (fb - fa);
        const T sigma = (m >= xf) ? T(1) : T(-1);
        const T delta = k1 * std::pow(width, k2);
        const T xt = (std::abs(m - xf) <= delta) ? (xf + sigma * delta) : m;

        T r = xtol * std::ldexp(T(1), std::max(0, n_max - iter)) - T(0.5) * width;
        r = std::max(T(0), r);

        T x = (std::abs(xt - m) <= r) ? xt : (m - sigma * r);
        x = std::clamp(x, a, b);
        const T fx = eval(x);

        if (std::abs(fx) <= ftol)
        {
            if (out_iterations)
                *out_iterations = iter + 1;
            if (out_converged)
                *out_converged = true;
            return clamp01<T>(x);
        }

        if (fa * fx <= T(0))
        {
            b = x;
            fb = fx;
        }
        else
        {
            a = x;
            fa = fx;
        }
    }

    if (out_iterations)
        *out_iterations = iter_max;
    return clamp01<T>(T(0.5) * (a + b));
}

template <std::floating_point T, typename Eval>
inline T newton_parameter(Eval&& eval, T a, T b, T fa, T fb, const T initial_guess,
                          const int max_iter, const T xtol, const T ftol,
                          int* out_iterations = nullptr,
                          bool* out_converged = nullptr)
{
    if (out_iterations)
        *out_iterations = 0;
    if (out_converged)
        *out_converged = false;

    if (std::abs(fa) <= ftol)
    {
        if (out_converged)
            *out_converged = true;
        return a;
    }
    if (std::abs(fb) <= ftol)
    {
        if (out_converged)
            *out_converged = true;
        return b;
    }
    if (fa * fb > T(0))
        return clamp01<T>(initial_guess);

    T l = a;
    T r = b;
    T fl = fa;
    T fr = fb;
    T t = std::clamp(initial_guess, l, r);
    const T eps = std::numeric_limits<T>::epsilon();

    for (int iter = 0; iter < max_iter; ++iter)
    {
        if ((r - l) <= xtol)
        {
            if (out_iterations)
                *out_iterations = iter + 1;
            if (out_converged)
                *out_converged = true;
            return clamp01<T>(T(0.5) * (l + r));
        }

        const T f = eval(t);
        if (std::abs(f) <= ftol)
        {
            if (out_iterations)
                *out_iterations = iter + 1;
            if (out_converged)
                *out_converged = true;
            return clamp01<T>(t);
        }

        if (fl * f <= T(0))
        {
            r = t;
            fr = f;
        }
        else if (f * fr <= T(0))
        {
            l = t;
            fl = f;
        }

        const T width = r - l;
        T h = std::sqrt(eps) * std::max(T(1), std::abs(t));
        h = std::max(h, T(0.01) * width);
        h = std::min(h, T(0.25) * width);
        h = std::max(h, xtol);

        const T t_minus = std::max(l, t - h);
        const T t_plus = std::min(r, t + h);
        T t_candidate = T(0.5) * (l + r);

        if (t_plus > t_minus + xtol)
        {
            const T f_minus = eval(t_minus);
            const T f_plus = eval(t_plus);
            const T df = (f_plus - f_minus) / (t_plus - t_minus);

            if (std::isfinite(df) && std::abs(df) > T(10) * eps)
            {
                const T t_newton = t - f / df;
                if (std::isfinite(t_newton) && t_newton > l && t_newton < r)
                    t_candidate = t_newton;
            }
        }

        if (!(t_candidate > l && t_candidate < r))
            t_candidate = T(0.5) * (l + r);
        t = t_candidate;
    }

    if (out_iterations)
        *out_iterations = max_iter;
    return clamp01<T>(t);
}

template <std::floating_point T, typename Phi>
inline RootSolveInfo<T> find_root_parameter_info(const std::span<const T> p0,
                                                 const std::span<const T> p1,
                                                 Phi&& level_set,
                                                 const method root_method = method::linear,
                                                 const T level = T(0),
                                                 const int max_iter = 64,
                                                 const T xtol = T(1e-12),
                                                 const T ftol = T(1e-12))
{
    if (p0.size() != p1.size())
        throw std::invalid_argument("find_root_parameter: inconsistent point dimensions");

    RootSolveInfo<T> info;
    int eval_count = 0;
    auto eval = [&](const T t) -> T
    {
        ++eval_count;
        return evaluate_edge_function<T>(level_set, p0, p1, t) - level;
    };

    T a = T(0);
    T b = T(1);
    T fa = eval(a);
    T fb = eval(b);
    const T a0 = a;
    const T b0 = b;
    const T fa0 = fa;
    const T fb0 = fb;

    if (std::abs(fa) <= ftol)
    {
        info.t = T(0);
        info.residual = std::abs(fa);
        info.evaluations = eval_count;
        info.converged = true;
        return info;
    }
    if (std::abs(fb) <= ftol)
    {
        info.t = T(1);
        info.residual = std::abs(fb);
        info.evaluations = eval_count;
        info.converged = true;
        return info;
    }

    const T t_linear = linear_root_parameter<T>(fa, fb, T(0));
    if (root_method == method::linear)
    {
        const T f_linear = eval(t_linear);
        info.t = t_linear;
        info.residual = std::abs(f_linear);
        info.evaluations = eval_count;
        info.converged = (std::abs(f_linear) <= ftol);
        return info;
    }

    // Contract the bracket with the linear guess when possible.
    bool have_f_linear = false;
    T f_linear = T(0);
    if (t_linear > T(0) && t_linear < T(1))
    {
        f_linear = eval(t_linear);
        have_f_linear = true;
        if (std::abs(f_linear) <= ftol)
        {
            info.t = t_linear;
            info.residual = std::abs(f_linear);
            info.evaluations = eval_count;
            info.converged = true;
            return info;
        }
        if (fa * f_linear <= T(0))
        {
            b = t_linear;
            fb = f_linear;
        }
        else if (f_linear * fb <= T(0))
        {
            a = t_linear;
            fa = f_linear;
        }
    }

    if (fa * fb > T(0))
    {
        if (!have_f_linear)
            f_linear = eval(t_linear);
        info.t = t_linear;
        info.residual = std::abs(f_linear);
        info.evaluations = eval_count;
        return info;
    }

    int iterations = 0;
    bool converged = false;
    T t = t_linear;

    if (root_method == method::brent)
    {
        t = brent_parameter<T>(eval, a, b, fa, fb, max_iter, xtol, ftol, &iterations, &converged);
    }
    else if (root_method == method::itp)
    {
        t = itp_parameter<T>(eval, a, b, fa, fb, t_linear, max_iter, xtol, ftol, &iterations, &converged);
    }
    else
    {
        t = newton_parameter<T>(eval, a, b, fa, fb, t_linear, max_iter, xtol, ftol, &iterations, &converged);
    }

    const T f = eval(t);
    T f_abs = std::abs(f);

    // Safety fallback: if nonlinear solver returned a poor residual while we
    // have a valid sign-changing bracket, finish with robust bisection.
    auto bisect_on_bracket = [&](T left, T right, T f_left, T f_right) -> T
    {
        if (f_left * f_right > T(0))
            return t;

        T l = left;
        T r = right;
        T fl = f_left;
        T fr = f_right;
        T mid = T(0.5) * (l + r);
        for (int k = 0; k < max_iter; ++k)
        {
            mid = T(0.5) * (l + r);
            const T fm = eval(mid);
            ++iterations;
            if (std::abs(fm) <= ftol || std::abs(r - l) <= xtol)
                return mid;

            if (fl * fm <= T(0))
            {
                r = mid;
                fr = fm;
            }
            else
            {
                l = mid;
                fl = fm;
            }
        }
        (void)fr;
        return mid;
    };

    if (f_abs > ftol)
    {
        if (fa * fb <= T(0))
        {
            t = bisect_on_bracket(a, b, fa, fb);
        }
        else if (fa0 * fb0 <= T(0))
        {
            t = bisect_on_bracket(a0, b0, fa0, fb0);
        }
        f_abs = std::abs(eval(t));
    }

    info.t = t;
    info.residual = f_abs;
    info.iterations = iterations;
    info.evaluations = eval_count;
    info.converged = (f_abs <= ftol) || (converged && f_abs <= T(100) * ftol);
    return info;
}

template <std::floating_point T, typename Phi>
inline T find_root_parameter(const std::span<const T> p0, const std::span<const T> p1, Phi&& level_set,
                             const method root_method = method::linear,
                             const T level = T(0), const int max_iter = 64,
                             const T xtol = T(1e-12), const T ftol = T(1e-12))
{
    return find_root_parameter_info<T>(p0, p1, std::forward<Phi>(level_set),
                                       root_method, level, max_iter, xtol, ftol).t;
}

template <std::floating_point T, typename Phi>
inline void find_root_point(const std::span<const T> p0, const std::span<const T> p1,
                            const std::span<T> root_point, Phi&& level_set,
                            const method root_method = method::linear,
                            const T level = T(0), const int max_iter = 64,
                            const T xtol = T(1e-12), const T ftol = T(1e-12))
{
    const T t = find_root_parameter<T>(p0, p1, std::forward<Phi>(level_set),
                                       root_method, level, max_iter, xtol, ftol);
    interpolate_point<T>(p0, p1, t, root_point);
}

template <std::floating_point T, typename EdgeMap>
inline void compute_case_intersections_from_edge_ids(
    const std::span<const T> vertex_coordinates, const int gdim,
    const std::span<const T> ls_values,
    const int (*edges)[2], const int num_edges,
    const std::span<const int> intersected_edge_ids,
    std::vector<T>& intersection_points,
    EdgeMap& edge_ip_map)
{
    if (gdim <= 0 || gdim > 3)
        throw std::invalid_argument("compute_case_intersections_from_edge_ids: gdim must be 1..3");

    intersection_points.resize(intersected_edge_ids.size() * static_cast<std::size_t>(gdim));

    std::array<T, 3> p0 = {};
    std::array<T, 3> p1 = {};

    for (std::size_t ip = 0; ip < intersected_edge_ids.size(); ++ip)
    {
        const int edge_id = intersected_edge_ids[ip];
        if (edge_id < 0 || edge_id >= num_edges)
            throw std::invalid_argument("compute_case_intersections_from_edge_ids: invalid edge id");

        const int v0 = edges[edge_id][0];
        const int v1 = edges[edge_id][1];

        for (int j = 0; j < gdim; ++j)
        {
            p0[static_cast<std::size_t>(j)] = vertex_coordinates[v0 * gdim + j];
            p1[static_cast<std::size_t>(j)] = vertex_coordinates[v1 * gdim + j];
        }

        linear_intersection_point<T>(std::span<const T>(p0.data(), static_cast<std::size_t>(gdim)),
                                     std::span<const T>(p1.data(), static_cast<std::size_t>(gdim)),
                                     ls_values[v0], ls_values[v1],
                                     intersection_points, static_cast<int>(ip * gdim), T(0));
        edge_ip_map[edge_id] = static_cast<int>(ip);
    }
}

template <std::floating_point T, typename EdgeMap>
inline void compute_case_intersections_from_edge_mask(
    const std::span<const T> vertex_coordinates, const int gdim,
    const std::span<const T> ls_values,
    const int (*edges)[2], const int num_edges,
    const std::span<const int> intersected_edge_mask,
    std::vector<T>& intersection_points,
    EdgeMap& edge_ip_map)
{
    if (gdim <= 0 || gdim > 3)
        throw std::invalid_argument("compute_case_intersections_from_edge_mask: gdim must be 1..3");
    if (intersected_edge_mask.size() != static_cast<std::size_t>(num_edges))
        throw std::invalid_argument("compute_case_intersections_from_edge_mask: invalid mask size");

    intersection_points.clear();
    intersection_points.reserve(static_cast<std::size_t>(num_edges) * static_cast<std::size_t>(gdim));

    std::array<T, 3> p0 = {};
    std::array<T, 3> p1 = {};

    int ip = 0;
    for (int edge_id = 0; edge_id < num_edges; ++edge_id)
    {
        if (intersected_edge_mask[static_cast<std::size_t>(edge_id)] == 0)
            continue;

        const int v0 = edges[edge_id][0];
        const int v1 = edges[edge_id][1];

        for (int j = 0; j < gdim; ++j)
        {
            p0[static_cast<std::size_t>(j)] = vertex_coordinates[v0 * gdim + j];
            p1[static_cast<std::size_t>(j)] = vertex_coordinates[v1 * gdim + j];
        }

        const int offset = static_cast<int>(intersection_points.size());
        intersection_points.resize(intersection_points.size() + static_cast<std::size_t>(gdim));
        linear_intersection_point<T>(std::span<const T>(p0.data(), static_cast<std::size_t>(gdim)),
                                     std::span<const T>(p1.data(), static_cast<std::size_t>(gdim)),
                                     ls_values[v0], ls_values[v1], intersection_points, offset, T(0));
        edge_ip_map[edge_id] = ip++;
    }
}
} // namespace cutcells::cell::edge_root
