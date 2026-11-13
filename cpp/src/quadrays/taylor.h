// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <limits>
#include <type_traits>

/// Scalar types for evaluating a templated level-set functor
/// (template <typename T> T operator()(const std::array<T, 3>& x) const):
/// Taylor<T, N> for bounds over a region, Dual<V, N> for derivatives,
/// Dual<Taylor<T, N>, N> for bounds of the derivatives and
/// Dual<Dual<Taylor<T, N>, N>, N> for bounds of the second ones. Functors call sqrt,
/// exp, log, sin, cos, abs, min and max unqualified (after using std::sqrt and
/// so on), so that argument-dependent lookup finds the versions below.
namespace cutcells::quadrays
{

// ============================================================================
// First-order Taylor models
// ============================================================================

/// First-order Taylor model of a function of y in [-1, 1]^N:
///
///   f(y) in alpha + beta . y + [-eps, eps].
///
/// The arithmetic is Saye's (2015) and algoim's Interval<N>, with the variables
/// scaled to [-1, 1], so the half-widths need no storage. sqrt, log, exp and
/// 1 / x are concave or convex about alpha, so their remainders come from
/// their values at the ends of the range, which is tighter than the
/// second-order bound C/2 b^2 and keeps a linear model of sqrt where the range
/// reaches 0; a square's remainder is one-sided. A distance function
/// sqrt(sum (x_i - c_i)^2) then certifies boxes several times larger than with
/// C/2 b^2. An operation that cannot bound its result (sqrt of a negative
/// range, log of a range reaching 0, division by a range containing 0) returns
/// a lost model, of NaN entries, and so does every operation on one: bounded()
/// tells whether a result holds. algoim throws instead, but bounds are lost
/// often (around the centre of a distance function, say), and an exception
/// costs far more than the bound.
template <std::floating_point T, int N>
struct Taylor
{
    T alpha = 0;
    std::array<T, N> beta{};
    T eps = 0;

    Taylor() = default;
    /// The constant c.
    explicit Taylor(T c) : alpha(c) {}
};

/// Largest |f(y) - alpha| over y in [-1, 1]^N.
template <std::floating_point T, int N>
T deviation(const Taylor<T, N>& a)
{
    T b = a.eps;
    for (int k = 0; k < N; ++k)
        b += std::abs(a.beta[k]);
    return b;
}

/// The result of an operation that cannot bound it.
template <std::floating_point T, int N>
Taylor<T, N> lost_model()
{
    constexpr T nan = std::numeric_limits<T>::quiet_NaN();
    Taylor<T, N> r(nan);
    r.beta.fill(nan);
    r.eps = nan;
    return r;
}

/// True if the model bounds its function: no operation on the way lost the
/// bound (lost_model) and nothing overflowed.
template <std::floating_point T, int N>
bool bounded(const Taylor<T, N>& a)
{
    bool finite = std::isfinite(a.alpha) && std::isfinite(a.eps);
    for (int k = 0; k < N; ++k)
        finite = finite && std::isfinite(a.beta[k]);
    return finite;
}

/// +1 or -1 if the model has that sign everywhere, 0 otherwise. A tolerance of
/// a few rounding errors lets x -> x count as positive on [0, 1].
template <std::floating_point T, int N>
int certain_sign(const Taylor<T, N>& a)
{
    const T x = (T(1) - T(10) * std::numeric_limits<T>::epsilon()) * deviation(a);
    if (a.alpha > x)
        return 1;
    if (a.alpha < -x)
        return -1;
    return 0;
}

/// The constant model with values in [lo, hi].
template <std::floating_point T, int N>
Taylor<T, N> constant_range(T lo, T hi)
{
    Taylor<T, N> r(T(0.5) * (lo + hi));
    r.eps = T(0.5) * (hi - lo);
    return r;
}

template <std::floating_point T, int N>
Taylor<T, N> operator-(const Taylor<T, N>& a)
{
    Taylor<T, N> r = a;
    r.alpha = -a.alpha;
    for (int k = 0; k < N; ++k)
        r.beta[k] = -a.beta[k];
    return r;
}

template <std::floating_point T, int N>
Taylor<T, N> operator+(const Taylor<T, N>& a, const Taylor<T, N>& b)
{
    Taylor<T, N> r;
    r.alpha = a.alpha + b.alpha;
    for (int k = 0; k < N; ++k)
        r.beta[k] = a.beta[k] + b.beta[k];
    r.eps = a.eps + b.eps;
    return r;
}

template <std::floating_point T, int N>
Taylor<T, N> operator-(const Taylor<T, N>& a, const Taylor<T, N>& b)
{
    return a + (-b);
}

template <std::floating_point T, int N>
Taylor<T, N> operator*(const Taylor<T, N>& a, const Taylor<T, N>& b)
{
    // A square, the same model twice (x * x in a functor): (alpha + l)^2 =
    // alpha^2 + 2 alpha l + l^2 with l^2 in [0, d^2], one-sided, half the
    // remainder of a product. Two models with equal entries may bound
    // different functions, whose product is no square.
    if (&a == &b)
    {
        const T d = deviation(a);
        Taylor<T, N> r;
        r.alpha = a.alpha * a.alpha + T(0.5) * d * d;
        for (int k = 0; k < N; ++k)
            r.beta[k] = T(2) * a.alpha * a.beta[k];
        r.eps = T(2) * std::abs(a.alpha) * a.eps + T(0.5) * d * d;
        return r;
    }
    // (alpha_a + l_a)(alpha_b + l_b) with |l| <= deviation
    Taylor<T, N> r;
    r.eps = deviation(a) * deviation(b) + std::abs(a.alpha) * b.eps + std::abs(b.alpha) * a.eps;
    for (int k = 0; k < N; ++k)
        r.beta[k] = a.alpha * b.beta[k] + b.alpha * a.beta[k];
    r.alpha = a.alpha * b.alpha;
    return r;
}

/// 1 / b about alpha: 1 / x = 1 / alpha - (x - alpha) / alpha^2 + g(x) with g
/// convex and g(alpha) = 0 its minimum, so that on the range [alpha - d,
/// alpha + d], which must not contain 0, g lies between 0 and its larger end
/// value, d^2 / (alpha^2 (|alpha| - d)).
template <std::floating_point T, int N>
Taylor<T, N> reciprocal(const Taylor<T, N>& b)
{
    if (!bounded(b))
        return lost_model<T, N>();
    const T d = deviation(b), m = std::abs(b.alpha);
    if (!(d < m))
        return lost_model<T, N>(); // the range may contain 0
    const T inv = T(1) / b.alpha;
    // g on the end nearer 0, with the sign of 1 / alpha, and a few rounding errors
    const T gmax = d * d / (b.alpha * b.alpha * (m - d)) * (T(1) + T(8) * std::numeric_limits<T>::epsilon());
    const T g = b.alpha > T(0) ? gmax : -gmax;
    Taylor<T, N> r;
    r.alpha = inv + T(0.5) * g;
    for (int k = 0; k < N; ++k)
        r.beta[k] = -b.beta[k] * inv * inv;
    r.eps = b.eps * inv * inv + T(0.5) * gmax;
    return r;
}

template <std::floating_point T, int N>
Taylor<T, N> operator/(const Taylor<T, N>& a, const Taylor<T, N>& b)
{
    return a * reciprocal(b);
}

// with constants
template <std::floating_point T, int N>
Taylor<T, N> operator+(const Taylor<T, N>& a, std::type_identity_t<T> c)
{
    Taylor<T, N> r = a;
    r.alpha += c;
    return r;
}
template <std::floating_point T, int N>
Taylor<T, N> operator+(std::type_identity_t<T> c, const Taylor<T, N>& a)
{
    return a + c;
}
template <std::floating_point T, int N>
Taylor<T, N> operator-(const Taylor<T, N>& a, std::type_identity_t<T> c)
{
    return a + (-c);
}
template <std::floating_point T, int N>
Taylor<T, N> operator-(std::type_identity_t<T> c, const Taylor<T, N>& a)
{
    return (-a) + c;
}
template <std::floating_point T, int N>
Taylor<T, N> operator*(const Taylor<T, N>& a, std::type_identity_t<T> c)
{
    Taylor<T, N> r;
    r.alpha = a.alpha * c;
    for (int k = 0; k < N; ++k)
        r.beta[k] = a.beta[k] * c;
    r.eps = a.eps * std::abs(c);
    return r;
}
template <std::floating_point T, int N>
Taylor<T, N> operator*(std::type_identity_t<T> c, const Taylor<T, N>& a)
{
    return a * c;
}
template <std::floating_point T, int N>
Taylor<T, N> operator/(const Taylor<T, N>& a, std::type_identity_t<T> c)
{
    return a * (T(1) / c);
}
template <std::floating_point T, int N>
Taylor<T, N> operator/(std::type_identity_t<T> c, const Taylor<T, N>& a)
{
    return Taylor<T, N>(c) / a;
}

/// sqrt about alpha: sqrt(x) = s + (x - alpha) / (2 s) + g(x) with g concave
/// and g(alpha) = 0 its maximum, so that on the range [max(alpha - b, 0),
/// alpha + b] g lies between its smaller end value and 0. This bound never
/// exceeds the second-order remainder C/2 b^2 (C bounding |sqrt''| on the
/// range) and stays finite where the range reaches 0 (around the centre of a
/// distance function, where that remainder has no bound).
template <std::floating_point T, int N>
Taylor<T, N> sqrt(const Taylor<T, N>& a)
{
    if (!bounded(a))
        return lost_model<T, N>();
    const T b = deviation(a);
    const T hi = a.alpha + b;
    if (!(hi >= T(0)))
        return lost_model<T, N>(); // a negative range
    if (!(a.alpha > T(0)))
        return constant_range<T, N>(T(0), std::sqrt(hi));
    const T s = std::sqrt(a.alpha), lo = std::max(a.alpha - b, T(0));
    const auto g = [&](T x) { return std::sqrt(x) - s - (x - a.alpha) / (2 * s); };
    // a few rounding errors of the terms of g
    const T pad = T(8) * std::numeric_limits<T>::epsilon() * (s + b / (2 * s));
    const T gmin = std::min(g(lo), g(hi)) - pad;
    Taylor<T, N> r;
    r.alpha = s + T(0.5) * gmin;
    for (int k = 0; k < N; ++k)
        r.beta[k] = T(0.5) / s * a.beta[k];
    r.eps = T(0.5) / s * a.eps - T(0.5) * gmin;
    return r;
}

/// log about alpha: log(x) = log(alpha) + (x - alpha) / alpha + g(x) with g
/// concave and g(alpha) = 0 its maximum, so that on [alpha - b, alpha + b],
/// which must not reach 0, g lies between its smaller end value and 0.
template <std::floating_point T, int N>
Taylor<T, N> log(const Taylor<T, N>& a)
{
    if (!bounded(a))
        return lost_model<T, N>();
    const T b = deviation(a);
    if (!(b < a.alpha))
        return lost_model<T, N>(); // the range reaches 0
    const T inv = T(1) / a.alpha, l = std::log(a.alpha);
    const auto g = [&](T x) { return std::log(x) - l - (x - a.alpha) * inv; };
    const T pad = T(8) * std::numeric_limits<T>::epsilon() * (std::abs(l) + b * inv);
    const T gmin = std::min(g(a.alpha - b), g(a.alpha + b)) - pad;
    Taylor<T, N> r;
    r.alpha = l + T(0.5) * gmin;
    for (int k = 0; k < N; ++k)
        r.beta[k] = inv * a.beta[k];
    r.eps = inv * a.eps - T(0.5) * gmin;
    return r;
}

/// exp about alpha: exp(x) = e + e (x - alpha) + g(x) with g convex and
/// g(alpha) = 0 its minimum, so that on [alpha - b, alpha + b] g lies between
/// 0 and its larger end value.
template <std::floating_point T, int N>
Taylor<T, N> exp(const Taylor<T, N>& a)
{
    if (!bounded(a))
        return lost_model<T, N>();
    const T b = deviation(a);
    const T e = std::exp(a.alpha);
    const auto g = [&](T x) { return std::exp(x) - e - e * (x - a.alpha); };
    const T pad = T(8) * std::numeric_limits<T>::epsilon() * (std::exp(a.alpha + b) + e * b);
    const T gmax = std::max(g(a.alpha - b), g(a.alpha + b)) + pad;
    Taylor<T, N> r;
    r.alpha = e + T(0.5) * gmax;
    for (int k = 0; k < N; ++k)
        r.beta[k] = e * a.beta[k];
    r.eps = e * a.eps + T(0.5) * gmax;
    return r;
}

template <std::floating_point T, int N>
Taylor<T, N> sin(const Taylor<T, N>& a)
{
    const T b = deviation(a);
    const T c = std::cos(a.alpha);
    Taylor<T, N> r;
    r.alpha = std::sin(a.alpha);
    for (int k = 0; k < N; ++k)
        r.beta[k] = c * a.beta[k];
    r.eps = std::abs(c) * a.eps + T(0.5) * b * b;
    return r;
}

template <std::floating_point T, int N>
Taylor<T, N> cos(const Taylor<T, N>& a)
{
    const T b = deviation(a);
    const T s = -std::sin(a.alpha);
    Taylor<T, N> r;
    r.alpha = std::cos(a.alpha);
    for (int k = 0; k < N; ++k)
        r.beta[k] = s * a.beta[k];
    r.eps = std::abs(s) * a.eps + T(0.5) * b * b;
    return r;
}

/// Branch helpers: the branch that holds everywhere, or a constant model
/// enclosing both. A lost operand gives a lost result, which std::min and
/// std::max on its NaN entries would not.
template <std::floating_point T, int N>
bool certainly_less(const Taylor<T, N>& a, const Taylor<T, N>& b)
{
    return certain_sign(a - b) < 0;
}

template <std::floating_point T, int N>
Taylor<T, N> hull(const Taylor<T, N>& a, const Taylor<T, N>& b)
{
    if (!bounded(a) || !bounded(b))
        return lost_model<T, N>();
    const T da = deviation(a), db = deviation(b);
    return constant_range<T, N>(std::min(a.alpha - da, b.alpha - db), std::max(a.alpha + da, b.alpha + db));
}

template <std::floating_point T, int N>
Taylor<T, N> abs(const Taylor<T, N>& a)
{
    if (!bounded(a))
        return lost_model<T, N>();
    const int s = certain_sign(a);
    if (s > 0)
        return a;
    if (s < 0)
        return -a;
    const T d = deviation(a);
    return constant_range<T, N>(T(0), std::max(d - a.alpha, a.alpha + d));
}

template <std::floating_point T, int N>
Taylor<T, N> min(const Taylor<T, N>& a, const Taylor<T, N>& b)
{
    if (!bounded(a) || !bounded(b))
        return lost_model<T, N>();
    if (certainly_less(a, b))
        return a;
    if (certainly_less(b, a))
        return b;
    const T da = deviation(a), db = deviation(b);
    return constant_range<T, N>(std::min(a.alpha - da, b.alpha - db), std::min(a.alpha + da, b.alpha + db));
}

template <std::floating_point T, int N>
Taylor<T, N> max(const Taylor<T, N>& a, const Taylor<T, N>& b)
{
    if (!bounded(a) || !bounded(b))
        return lost_model<T, N>();
    if (certainly_less(a, b))
        return b;
    if (certainly_less(b, a))
        return a;
    const T da = deviation(a), db = deviation(b);
    return constant_range<T, N>(std::max(a.alpha - da, b.alpha - db), std::max(a.alpha + da, b.alpha + db));
}

// ============================================================================
// Forward-mode derivatives
// ============================================================================

/// Value and first derivatives with respect to N variables; V is a floating-point
/// type or a Taylor model.
template <typename V, int N>
struct Dual
{
    V v{};
    std::array<V, N> d{};

    Dual() = default;
    /// The constant c.
    explicit Dual(double c) : v(c) {}
};

/// Sign certainty, enclosures and bounds of plain numbers, for Dual<double, N>.
template <std::floating_point T>
int certain_sign(T v)
{
    return (v > T(0)) - (v < T(0));
}
template <std::floating_point T>
bool bounded(T v)
{
    return std::isfinite(v);
}
template <std::floating_point T>
bool certainly_less(T a, T b)
{
    return a < b;
}
/// Reached for plain numbers only on ties, where both branches agree.
template <std::floating_point T>
T hull(T a, T)
{
    return a;
}

/// Sign certainty of a Dual's value.
template <typename V, int N>
int certain_sign(const Dual<V, N>& a)
{
    return certain_sign(a.v);
}

/// Comparison of the values, for Duals of Duals (second derivatives).
template <typename V, int N>
bool certainly_less(const Dual<V, N>& a, const Dual<V, N>& b)
{
    return certainly_less(a.v, b.v);
}

/// True if the value and every derivative are bounded.
template <typename V, int N>
bool bounded(const Dual<V, N>& a)
{
    bool all = bounded(a.v);
    for (int k = 0; k < N; ++k)
        all = all && bounded(a.d[k]);
    return all;
}

template <typename V, int N>
Dual<V, N> operator-(const Dual<V, N>& a)
{
    Dual<V, N> r;
    r.v = -a.v;
    for (int k = 0; k < N; ++k)
        r.d[k] = -a.d[k];
    return r;
}

template <typename V, int N>
Dual<V, N> operator+(const Dual<V, N>& a, const Dual<V, N>& b)
{
    Dual<V, N> r;
    r.v = a.v + b.v;
    for (int k = 0; k < N; ++k)
        r.d[k] = a.d[k] + b.d[k];
    return r;
}

template <typename V, int N>
Dual<V, N> operator-(const Dual<V, N>& a, const Dual<V, N>& b)
{
    Dual<V, N> r;
    r.v = a.v - b.v;
    for (int k = 0; k < N; ++k)
        r.d[k] = a.d[k] - b.d[k];
    return r;
}

template <typename V, int N>
Dual<V, N> operator*(const Dual<V, N>& a, const Dual<V, N>& b)
{
    Dual<V, N> r;
    r.v = a.v * b.v;
    for (int k = 0; k < N; ++k)
        r.d[k] = a.v * b.d[k] + b.v * a.d[k];
    return r;
}

template <typename V, int N>
Dual<V, N> operator/(const Dual<V, N>& a, const Dual<V, N>& b)
{
    Dual<V, N> r;
    r.v = a.v / b.v;
    for (int k = 0; k < N; ++k)
        r.d[k] = (a.d[k] - r.v * b.d[k]) / b.v;
    return r;
}

// with constants
template <typename V, int N>
Dual<V, N> operator+(const Dual<V, N>& a, double c)
{
    Dual<V, N> r = a;
    r.v = a.v + V(c);
    return r;
}
template <typename V, int N>
Dual<V, N> operator+(double c, const Dual<V, N>& a)
{
    return a + c;
}
template <typename V, int N>
Dual<V, N> operator-(const Dual<V, N>& a, double c)
{
    return a + (-c);
}
template <typename V, int N>
Dual<V, N> operator-(double c, const Dual<V, N>& a)
{
    return (-a) + c;
}
template <typename V, int N>
Dual<V, N> operator*(const Dual<V, N>& a, double c)
{
    Dual<V, N> r;
    r.v = a.v * V(c);
    for (int k = 0; k < N; ++k)
        r.d[k] = a.d[k] * V(c);
    return r;
}
template <typename V, int N>
Dual<V, N> operator*(double c, const Dual<V, N>& a)
{
    return a * c;
}
template <typename V, int N>
Dual<V, N> operator/(const Dual<V, N>& a, double c)
{
    return a * (1.0 / c);
}
template <typename V, int N>
Dual<V, N> operator/(double c, const Dual<V, N>& a)
{
    return Dual<V, N>(c) / a;
}

/// True if a derivative is exactly 0: plain numbers, and models or Duals all of
/// whose entries are 0, so that the function does not vary over the region.
template <std::floating_point T>
bool identically_zero(T v)
{
    return v == T(0);
}
template <std::floating_point T, int N>
bool identically_zero(const Taylor<T, N>& a)
{
    return a.alpha == T(0) && a.eps == T(0) && std::all_of(a.beta.begin(), a.beta.end(), [](T b) { return b == T(0); });
}
template <typename V, int N>
bool identically_zero(const Dual<V, N>& a)
{
    return identically_zero(a.v)
           && std::all_of(a.d.begin(), a.d.end(), [](const V& d) { return identically_zero(d); });
}

/// sqrt(a), d sqrt(a) = d a / (2 sqrt(a)). A derivative of a that is exactly 0
/// stays 0, also where sqrt(a) vanishes and the quotient is 0 / 0: the length
/// of (max(d_0, 0), max(d_1, 0)) in a capped shape's distance, inside the
/// shape, does not vary, and neither point gradients nor Taylor models of it
/// are undefined. Where a's derivative does not vanish, as for a distance at
/// its centre, the derivative's bound is lost.
template <typename V, int N>
Dual<V, N> sqrt(const Dual<V, N>& a)
{
    using std::sqrt;
    Dual<V, N> r;
    r.v = sqrt(a.v);
    const V twice = r.v * V(2.0);
    for (int k = 0; k < N; ++k)
        r.d[k] = identically_zero(a.d[k]) ? V(0.0) : a.d[k] / twice;
    return r;
}

template <typename V, int N>
Dual<V, N> log(const Dual<V, N>& a)
{
    using std::log;
    Dual<V, N> r;
    r.v = log(a.v);
    for (int k = 0; k < N; ++k)
        r.d[k] = a.d[k] / a.v;
    return r;
}

template <typename V, int N>
Dual<V, N> exp(const Dual<V, N>& a)
{
    using std::exp;
    Dual<V, N> r;
    r.v = exp(a.v);
    for (int k = 0; k < N; ++k)
        r.d[k] = r.v * a.d[k];
    return r;
}

template <typename V, int N>
Dual<V, N> sin(const Dual<V, N>& a)
{
    using std::cos;
    using std::sin;
    Dual<V, N> r;
    r.v = sin(a.v);
    const V c = cos(a.v);
    for (int k = 0; k < N; ++k)
        r.d[k] = c * a.d[k];
    return r;
}

template <typename V, int N>
Dual<V, N> cos(const Dual<V, N>& a)
{
    using std::cos;
    using std::sin;
    Dual<V, N> r;
    r.v = cos(a.v);
    const V s = -sin(a.v);
    for (int k = 0; k < N; ++k)
        r.d[k] = s * a.d[k];
    return r;
}

template <typename V, int N>
Dual<V, N> hull(const Dual<V, N>& a, const Dual<V, N>& b)
{
    Dual<V, N> r;
    r.v = hull(a.v, b.v);
    for (int k = 0; k < N; ++k)
        r.d[k] = hull(a.d[k], b.d[k]);
    return r;
}

/// |a|; where the sign is not certain, the derivative encloses both branches.
template <typename V, int N>
Dual<V, N> abs(const Dual<V, N>& a)
{
    using std::abs;
    const int s = certain_sign(a.v);
    if (s > 0)
        return a;
    if (s < 0)
        return -a;
    Dual<V, N> r = hull(a, -a);
    r.v = abs(a.v);
    return r;
}

template <typename V, int N>
Dual<V, N> min(const Dual<V, N>& a, const Dual<V, N>& b)
{
    using std::min;
    if (certainly_less(a.v, b.v))
        return a;
    if (certainly_less(b.v, a.v))
        return b;
    Dual<V, N> r = hull(a, b);
    r.v = min(a.v, b.v);
    return r;
}

template <typename V, int N>
Dual<V, N> max(const Dual<V, N>& a, const Dual<V, N>& b)
{
    using std::max;
    if (certainly_less(a.v, b.v))
        return b;
    if (certainly_less(b.v, a.v))
        return a;
    Dual<V, N> r = hull(a, b);
    r.v = max(a.v, b.v);
    return r;
}

/// a^p for an integer p, by repeated squaring.
template <typename V>
V integer_power(const V& a, long p)
{
    V result(1.0), base = a;
    for (long e = p < 0 ? -p : p; e > 0; e >>= 1)
    {
        if (e & 1)
            result = result * base;
        base = base * base;
    }
    return p < 0 ? V(1.0) / result : result;
}

} // namespace cutcells::quadrays
