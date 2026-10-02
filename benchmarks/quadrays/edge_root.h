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
template <std::floating_point T>
inline T default_value_tolerance()
{
    return T(64) * std::numeric_limits<T>::epsilon();
}

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

/// Result of a ray root-finding solve: phi(x_s + t*d) = 0.
template <std::floating_point T>
struct RayRootResult
{
    bool valid = false;
    bool ambiguous_two_sided = false;
    bool used_gradient_dir = false;
    bool used_fallback_dir = false;

    T    t = T(0);
    T    phi_at_root = std::numeric_limits<T>::max();
    T    dist = std::numeric_limits<T>::max();
    T    residual = std::numeric_limits<T>::max();
    int  iterations = 0;
    int  evaluations = 0;
    bool converged = false;
    bool domain_exit = false;
    bool gradient_degenerate = false;

    std::array<T, 3> x_ref = {T(0), T(0), T(0)};

    bool pos_valid = false;
    bool neg_valid = false;
    T pos_t = T(0);
    T neg_t = T(0);
};

template <std::floating_point T>
struct RaySearchOptions
{
    int max_iter = 64;
    T xtol = T(1e-12);
    T ftol = default_value_tolerance<T>();
    T domain_tol = T(1e-10);
    std::array<T, 3> probe_scales = {T(1), T(2), T(4)};
    T max_probe_distance = std::numeric_limits<T>::infinity();
};

template <std::floating_point T>
inline T linear_root_parameter(const T v0, const T v1, const T level = T(0))
{
    const T denom = v1 - v0;
    if (std::abs(denom) <= std::numeric_limits<T>::epsilon())
        return T(0.5);
    return (level - v0) / denom;
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

/// General-purpose Brent solver on an arbitrary real interval [a, b].
///
/// Unlike brent_parameter, this does NOT clamp the result to [0, 1].
/// Use this when the parameter domain extends beyond the unit interval,
/// e.g. for ray-parameter root finding.
///
/// @param eval  callable T(T) — function whose root is sought
/// @param a     left bracket
/// @param b     right bracket
/// @param fa    eval(a)
/// @param fb    eval(b) — must satisfy fa * fb <= 0
/// @param max_iter  maximum iterations
/// @param xtol  tolerance on bracket width
/// @param ftol  tolerance on function value
template <std::floating_point T, typename Eval>
inline T brent_solve(Eval&& eval, T a, T b, T fa, T fb,
                     int max_iter, T xtol, T ftol,
                     int* out_iterations = nullptr,
                     bool* out_converged = nullptr)
{
    if (out_iterations)
        *out_iterations = 0;
    if (out_converged)
        *out_converged = false;

    if (std::abs(fa) <= ftol)
    {
        if (out_converged) *out_converged = true;
        return a;
    }
    if (std::abs(fb) <= ftol)
    {
        if (out_converged) *out_converged = true;
        return b;
    }
    if (fa * fb > T(0))
        return T(0.5) * (a + b);

    T c = a, fc = fa;
    T d = b - a, e = d;

    for (int iter = 0; iter < max_iter; ++iter)
    {
        if ((fb > T(0) && fc > T(0)) || (fb < T(0) && fc < T(0)))
        {
            c = a; fc = fa; d = b - a; e = d;
        }
        if (std::abs(fc) < std::abs(fb))
        {
            a = b; fa = fb; b = c; fb = fc; c = a; fc = fa;
        }

        const T eps  = std::numeric_limits<T>::epsilon();
        const T tol1 = T(2) * eps * std::abs(b) + T(0.5) * xtol;
        const T xm   = T(0.5) * (c - b);

        if (std::abs(xm) <= tol1 || std::abs(fb) <= ftol)
        {
            if (out_iterations) *out_iterations = iter + 1;
            if (out_converged)  *out_converged = true;
            return b;
        }

        if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb))
        {
            T s = fb / fa;
            T p, q;
            if (a == c)
            {
                p = T(2) * xm * s;
                q = T(1) - s;
            }
            else
            {
                T q0 = fa / fc, r = fb / fc;
                p = s * (T(2) * xm * q0 * (q0 - r) - (b - a) * (r - T(1)));
                q = (q0 - T(1)) * (r - T(1)) * (s - T(1));
            }
            if (p > T(0)) q = -q;
            p = std::abs(p);

            const T min1 = T(3) * xm * q - std::abs(tol1 * q);
            const T min2 = std::abs(e * q);
            if (T(2) * p < std::min(min1, min2))
            { e = d; d = p / q; }
            else
            { d = xm; e = d; }
        }
        else
        {
            d = xm; e = d;
        }

        a = b; fa = fb;
        if (std::abs(d) > tol1) b += d;
        else b += (xm > T(0) ? tol1 : -tol1);
        fb = eval(b);
    }

    if (out_iterations) *out_iterations = max_iter;
    return b;
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
        throw std::domain_error("brent_parameter: no sign change in bracket [a, b]");

    T c = a;
    T fc = fa;
    T d = b - a;
    T e = d;

    for (int iter = 0; iter < max_iter; ++iter)
    {
        if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0))
        {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }

        if (std::abs(fc) < std::abs(fb))
        {
            // Three-way rotation: new b = c_old (better), new a = b_old, new c = b_old.
            // Execute sequentially so that after  a = b,  'a' already holds b_old.
            a = b;
            fa = fb;
            b = c;
            fb = fc;
            c = a;   // a is already b_old — maintains bracket [b, c] with f(b)*f(c) ≤ 0
            fc = fa; // fa is already fb_old
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
            return b;
        }

        if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb))
        {
            T s = fb / fa;
            T p, q;
            if (a == c)
            {
                p = T(2) * xm * s;
                q = T(1) - s;
            }
            else
            {
                T r = fb / fc;
                T q0 = fa / fc;
                p = s * (T(2) * xm * q0 * (q0 - r) - (b - a) * (r - T(1)));
                q = (q0 - T(1)) * (r - T(1)) * (s - T(1));
            }

            if (p > T(0))
                q = -q;
            p = std::abs(p);
            T min1 = T(3) * xm * q - std::abs(tol1 * q);
            T min2 = std::abs(e * q);
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
            b += (xm > 0 ? tol1 : -tol1);
        fb = eval(b);
    }

    if (out_iterations)
        *out_iterations = max_iter;
    return b;
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
        throw std::domain_error("itp_parameter: no sign change in bracket [a, b]");

    // Use the linear guess first to contract the bracket if possible.
    const T t0 = std::clamp(initial_guess, a, b);
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
            return T(0.5) * (a + b);
        }

        const T m = T(0.5) * (a + b);
        const T xf = (a * fb - b * fa) / (fb - fa);
        const T sigma = (m >= xf) ? T(1) : T(-1);

        // Truncation: delta = min(k1 * width^k2, |m - xf|)
        // Then xt = xf + sigma * delta  (clamps between xf and m).
        const T delta = std::min(k1 * std::pow(width, k2), std::abs(m - xf));
        const T xt = xf + sigma * delta;

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
            return x;
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
    return T(0.5) * (a + b);
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
        throw std::domain_error("newton_parameter: no sign change in bracket [a, b]");

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
            return T(0.5) * (l + r);
        }

        const T f = eval(t);
        if (std::abs(f) <= ftol)
        {
            if (out_iterations)
                *out_iterations = iter + 1;
            if (out_converged)
                *out_converged = true;
            return t;
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
    return t;
}

template <std::floating_point T, typename Phi>
inline RootSolveInfo<T> find_root_parameter_info(const std::span<const T> p0,
                                                 const std::span<const T> p1,
                                                 Phi&& level_set,
                                                 const method root_method = method::linear,
                                                 const T level = T(0),
                                                 const int max_iter = 64,
                                                 const T xtol = T(1e-12),
                                                 const T ftol = default_value_tolerance<T>())
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
                             const T xtol = T(1e-12),
                             const T ftol = default_value_tolerance<T>())
{
    return find_root_parameter_info<T>(p0, p1, std::forward<Phi>(level_set),
                                       root_method, level, max_iter, xtol, ftol).t;
}

template <std::floating_point T, typename Phi>
inline void find_root_point(const std::span<const T> p0, const std::span<const T> p1,
                            const std::span<T> root_point, Phi&& level_set,
                            const method root_method = method::linear,
                            const T level = T(0), const int max_iter = 64,
                            const T xtol = T(1e-12),
                            const T ftol = default_value_tolerance<T>())
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

/// Check whether a reference point is inside the parent cell reference domain.
///
/// @param x_ref      reference coordinates, size = tdim
/// @param cell_type  cell type defining the reference domain
/// @param tol        tolerance for boundary proximity
template <std::floating_point T>
inline bool is_inside_reference_domain(
    std::span<const T> x_ref,
    cutcells::cell::type cell_type,
    T                  tol = T(1e-10))
{
    using cutcells::cell::type;
    switch (cell_type)
    {
    case type::interval:
    {
        return x_ref[0] >= -tol && x_ref[0] <= T(1) + tol;
    }
    case type::triangle:
    {
        const T x = x_ref[0];
        const T y = x_ref[1];
        return x >= -tol && y >= -tol && (x + y) <= T(1) + tol;
    }
    case type::tetrahedron:
    {
        const T x = x_ref[0];
        const T y = x_ref[1];
        const T z = x_ref[2];
        return x >= -tol && y >= -tol && z >= -tol && (x + y + z) <= T(1) + tol;
    }
    case type::quadrilateral:
    {
        return x_ref[0] >= -tol && x_ref[0] <= T(1) + tol
            && x_ref[1] >= -tol && x_ref[1] <= T(1) + tol;
    }
    case type::hexahedron:
    {
        return x_ref[0] >= -tol && x_ref[0] <= T(1) + tol
            && x_ref[1] >= -tol && x_ref[1] <= T(1) + tol
            && x_ref[2] >= -tol && x_ref[2] <= T(1) + tol;
    }
    case type::prism:
    {
        // Prism: triangle in (x,y) x interval in z
        const T x = x_ref[0];
        const T y = x_ref[1];
        const T z = x_ref[2];
        return x >= -tol && y >= -tol && (x + y) <= T(1) + tol
            && z >= -tol && z <= T(1) + tol;
    }
    case type::pyramid:
    {
        // Pyramid: x,y in [0, 1-z], z in [0, 1]
        const T x = x_ref[0];
        const T y = x_ref[1];
        const T z = x_ref[2];
        const T limit = T(1) - z;
        return z >= -tol && z <= T(1) + tol
            && x >= -tol && x <= limit + tol
            && y >= -tol && y <= limit + tol;
    }
    default:
        return false;
    }
}

/// Solve phi(x_s + t*d) = 0 for t using Newton → Brent → bisection fallback chain.
///
/// @param eval_phi   callable (const T* x_ref, int tdim) -> T
/// @param eval_grad  callable (const T* x_ref, int tdim, T* grad_out) -> void (reference-space gradient)
/// @param x_s        starting point in reference coordinates, size = tdim
/// @param d          ray direction in reference coordinates, size = tdim
/// @param cell_type  reference domain for the domain-exit check
/// @param tdim       topological dimension
/// @param h_local    local entity size for scale-aware bracket probing (segment length or face diameter)
/// @param tol        convergence tolerance
template <std::floating_point T, typename EvalPhi, typename EvalGrad>
inline RayRootResult<T> find_ray_root_one_direction(
    EvalPhi&&            eval_phi,
    EvalGrad&&           eval_grad,
    std::span<const T>   x_s,
    std::span<const T>   d,
    cutcells::cell::type cell_type,
    int                  tdim,
    T                    h_local,
    int                  sign,
    const RaySearchOptions<T>& opts)
{
    RayRootResult<T> result;
    if (sign != 1 && sign != -1)
        throw std::invalid_argument("find_ray_root_one_direction: sign must be +1 or -1");
    if (tdim <= 0 || tdim > 3)
        throw std::invalid_argument("find_ray_root_one_direction: tdim must be in [1, 3]");

    // Compute norm of direction
    T d_norm = T(0);
    for (int k = 0; k < tdim; ++k)
        d_norm += d[k] * d[k];
    d_norm = std::sqrt(d_norm);

    if (d_norm < opts.ftol)
    {
        result.gradient_degenerate = true;
        return result;
    }

    // Normalized direction
    std::vector<T> d_hat(static_cast<std::size_t>(tdim));
    for (int k = 0; k < tdim; ++k)
        d_hat[static_cast<std::size_t>(k)] = d[k] / d_norm;

    // Helper: evaluate phi at x_s + t * d_hat
    std::vector<T> x_tmp(static_cast<std::size_t>(tdim));
    const auto phi_at_t = [&](T t_val) -> T
    {
        for (int k = 0; k < tdim; ++k)
            x_tmp[static_cast<std::size_t>(k)] = x_s[k] + t_val * d_hat[static_cast<std::size_t>(k)];
        ++result.evaluations;
        return eval_phi(x_tmp.data(), tdim);
    };

    T t_a = T(0);
    T f_a = phi_at_t(T(0));

    // Check if already at root
    if (std::abs(f_a) <= opts.ftol)
    {
        result.t = T(0);
        result.phi_at_root = f_a;
        result.dist = T(0);
        result.residual = std::abs(f_a);
        result.converged = true;
        result.valid = true;
        for (int k = 0; k < tdim; ++k)
            result.x_ref[static_cast<std::size_t>(k)] = x_s[static_cast<std::size_t>(k)];
        return result;
    }

    // Find a bracket
    T t_b = T(0);
    T f_b = T(0);
    bool have_bracket = false;

    for (const T scale : opts.probe_scales)
    {
        const T unclamped = scale * h_local;
        const T max_probe = std::isfinite(opts.max_probe_distance)
            ? std::min(unclamped, opts.max_probe_distance)
            : unclamped;
        const T t_probe = static_cast<T>(sign) * max_probe;
        // Check that probe is inside reference domain
        for (int k = 0; k < tdim; ++k)
            x_tmp[static_cast<std::size_t>(k)] = x_s[k] + t_probe * d_hat[static_cast<std::size_t>(k)];
        if (!is_inside_reference_domain<T>(
                std::span<const T>(x_tmp.data(), static_cast<std::size_t>(tdim)),
                cell_type, opts.domain_tol))
            continue;

        const T f_probe = phi_at_t(t_probe);
        if (std::abs(f_probe) <= opts.ftol)
        {
            result.t = t_probe;
            result.phi_at_root = f_probe;
            result.dist = std::abs(t_probe);
            result.residual = std::abs(f_probe);
            result.converged = true;
            result.valid = true;
            for (int k = 0; k < tdim; ++k)
                result.x_ref[static_cast<std::size_t>(k)]
                    = x_s[static_cast<std::size_t>(k)] + t_probe * d_hat[static_cast<std::size_t>(k)];
            return result;
        }

        if (f_a * f_probe <= T(0))
        {
            t_b = t_probe;
            f_b = f_probe;
            have_bracket = true;
            break;
        }

        // Keep the smaller-magnitude endpoint as t_a
        if (std::abs(f_probe) < std::abs(f_a))
        {
            t_a = t_probe;
            f_a = f_probe;
        }
    }

    if (!have_bracket)
    {
        // No sign change found; return best estimate
        result.t = t_a;
        result.phi_at_root = f_a;
        result.dist = std::abs(t_a);
        result.residual = std::abs(f_a);
        result.converged = false;
        return result;
    }

    // Ensure t_a < t_b
    if (t_a > t_b)
    {
        std::swap(t_a, t_b);
        std::swap(f_a, f_b);
    }

    // Newton iterations with bracket safety
    const int max_iter = opts.max_iter;
    T t = T(0.5) * (t_a + t_b);
    // Try initial Newton step from t_a
    {
        std::vector<T> grad(static_cast<std::size_t>(tdim));
        for (int k = 0; k < tdim; ++k)
            x_tmp[static_cast<std::size_t>(k)] = x_s[k] + t_a * d_hat[static_cast<std::size_t>(k)];
        eval_grad(x_tmp.data(), tdim, grad.data());
        T dphidt = T(0);
        for (int k = 0; k < tdim; ++k)
            dphidt += grad[static_cast<std::size_t>(k)] * d_hat[static_cast<std::size_t>(k)];
        if (std::abs(dphidt) > T(1e-14))
        {
            const T t_newton = t_a - f_a / dphidt;
            if (t_newton > t_a && t_newton < t_b)
                t = t_newton;
        }
    }

    // Brent fallback (unclamped — ray parameter can be negative)
    int brent_iters = 0;
    bool brent_converged = false;
    const T t_brent = brent_solve<T>(
        phi_at_t, t_a, t_b, f_a, f_b,
        max_iter, opts.xtol, opts.ftol, &brent_iters, &brent_converged);
    result.iterations += brent_iters;
    t = t_brent;

    const T f_final = phi_at_t(t);
    result.t = t;
    result.phi_at_root = f_final;
    result.dist = std::abs(t);
    result.residual = std::abs(f_final);
    result.converged = (result.residual <= opts.ftol) || brent_converged;

    // Check domain exit
    for (int k = 0; k < tdim; ++k)
        x_tmp[static_cast<std::size_t>(k)] = x_s[k] + t * d_hat[static_cast<std::size_t>(k)];
    result.domain_exit = !is_inside_reference_domain<T>(
        std::span<const T>(x_tmp.data(), static_cast<std::size_t>(tdim)),
        cell_type, opts.domain_tol);
    result.valid = result.converged && !result.domain_exit;
    for (int k = 0; k < tdim; ++k)
        result.x_ref[static_cast<std::size_t>(k)] = x_tmp[static_cast<std::size_t>(k)];

    return result;
}

template <std::floating_point T, typename EvalPhi, typename EvalGrad>
inline RayRootResult<T> find_ray_root(
    EvalPhi&&            eval_phi,
    EvalGrad&&           eval_grad,
    std::span<const T>   x_s,
    std::span<const T>   d,
    cutcells::cell::type cell_type,
    int                  tdim,
    T                    h_local,
    const RaySearchOptions<T>& opts)
{
    auto pos = find_ray_root_one_direction<T>(
        eval_phi, eval_grad, x_s, d, cell_type, tdim, h_local, +1, opts);
    auto neg = find_ray_root_one_direction<T>(
        eval_phi, eval_grad, x_s, d, cell_type, tdim, h_local, -1, opts);

    RayRootResult<T> out;
    out.evaluations = pos.evaluations + neg.evaluations;
    out.iterations = pos.iterations + neg.iterations;
    out.gradient_degenerate = pos.gradient_degenerate && neg.gradient_degenerate;

    out.pos_valid = pos.valid;
    out.neg_valid = neg.valid;
    out.pos_t = pos.t;
    out.neg_t = neg.t;
    out.ambiguous_two_sided = pos.valid && neg.valid;

    if (pos.valid && neg.valid)
        out = (std::abs(pos.t) <= std::abs(neg.t)) ? pos : neg;
    else if (pos.valid)
        out = pos;
    else if (neg.valid)
        out = neg;
    else
        out = (pos.residual <= neg.residual) ? pos : neg;

    out.pos_valid = pos.valid;
    out.neg_valid = neg.valid;
    out.pos_t = pos.t;
    out.neg_t = neg.t;
    out.ambiguous_two_sided = pos.valid && neg.valid;
    out.evaluations = pos.evaluations + neg.evaluations;
    out.iterations = pos.iterations + neg.iterations;
    out.gradient_degenerate = pos.gradient_degenerate && neg.gradient_degenerate;
    return out;
}

} // namespace cutcells::cell::edge_root
