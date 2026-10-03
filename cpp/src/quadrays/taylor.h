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
#include <stdexcept>
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
/// scaled to [-1, 1], so the half-widths need no storage, and with the full
/// remainder for sqrt and log. An operation that cannot bound its result (sqrt
/// or log of a range reaching 0, division by a range containing 0) throws
/// std::domain_error.
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
    // (alpha_a + l_a)(alpha_b + l_b) with |l| <= deviation
    Taylor<T, N> r;
    r.eps = deviation(a) * deviation(b) + std::abs(a.alpha) * b.eps + std::abs(b.alpha) * a.eps;
    for (int k = 0; k < N; ++k)
        r.beta[k] = a.alpha * b.beta[k] + b.alpha * a.beta[k];
    r.alpha = a.alpha * b.alpha;
    return r;
}

template <std::floating_point T, int N>
Taylor<T, N> operator/(const Taylor<T, N>& a, const Taylor<T, N>& b)
{
    const T inv = T(1) / b.alpha;
    const T tau = deviation(b) * std::abs(inv);
    if (!(tau < T(1)))
        throw std::domain_error("Taylor: division by a range that may contain 0");
    const T rho = tau * tau / ((T(1) - tau) * (T(1) - tau) * (T(1) - tau));
    const T b_num = deviation(a);
    Taylor<T, N> r;
    r.eps = (std::abs(a.alpha * inv) * b.eps + (std::abs(a.alpha) + b_num) * rho + b_num * tau + a.eps) * std::abs(inv);
    for (int k = 0; k < N; ++k)
        r.beta[k] = a.beta[k] * inv - a.alpha * b.beta[k] * inv * inv;
    r.alpha = a.alpha * inv;
    return r;
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

/// sqrt with the full first-order remainder |f'(alpha)| eps + C/2 b^2, C bounding
/// |sqrt''| = x^(-3/2) / 4 on [alpha - b, alpha + b]. algoim's sqrt leaves out
/// the first term, which declared cut cells uncut for |x - c| - r. A range
/// reaching 0 has no first-order model, only the values [0, sqrt(alpha + b)]
/// (around the centre of a distance function, say); its derivative, which a
/// Dual divides by sqrt, then loses its bound.
template <std::floating_point T, int N>
Taylor<T, N> sqrt(const Taylor<T, N>& a)
{
    const T b = deviation(a);
    if (!(b < a.alpha))
    {
        if (!(a.alpha + b >= T(0)))
            throw std::domain_error("Taylor: sqrt of a negative range");
        return constant_range<T, N>(T(0), std::sqrt(a.alpha + b));
    }
    const T s = std::sqrt(a.alpha);
    const T lo = a.alpha - b;
    const T C = T(0.25) / (lo * std::sqrt(lo));
    Taylor<T, N> r;
    r.alpha = s;
    for (int k = 0; k < N; ++k)
        r.beta[k] = T(0.5) / s * a.beta[k];
    r.eps = T(0.5) / s * a.eps + T(0.5) * C * b * b;
    return r;
}

/// log, with |log''| <= 1 / (alpha - b)^2.
template <std::floating_point T, int N>
Taylor<T, N> log(const Taylor<T, N>& a)
{
    const T b = deviation(a);
    if (!(b < a.alpha))
        throw std::domain_error("Taylor: log of a range reaching 0");
    const T inv = T(1) / a.alpha;
    const T C = T(1) / ((a.alpha - b) * (a.alpha - b));
    Taylor<T, N> r;
    r.alpha = std::log(a.alpha);
    for (int k = 0; k < N; ++k)
        r.beta[k] = inv * a.beta[k];
    r.eps = inv * a.eps + T(0.5) * C * b * b;
    return r;
}

/// exp, with exp'' bounded by its value at alpha + b.
template <std::floating_point T, int N>
Taylor<T, N> exp(const Taylor<T, N>& a)
{
    const T b = deviation(a);
    const T e = std::exp(a.alpha);
    Taylor<T, N> r;
    r.alpha = e;
    for (int k = 0; k < N; ++k)
        r.beta[k] = e * a.beta[k];
    r.eps = e * a.eps + T(0.5) * std::exp(a.alpha + b) * b * b;
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
/// enclosing both.
template <std::floating_point T, int N>
bool certainly_less(const Taylor<T, N>& a, const Taylor<T, N>& b)
{
    return certain_sign(a - b) < 0;
}

template <std::floating_point T, int N>
Taylor<T, N> hull(const Taylor<T, N>& a, const Taylor<T, N>& b)
{
    const T da = deviation(a), db = deviation(b);
    return constant_range<T, N>(std::min(a.alpha - da, b.alpha - db), std::max(a.alpha + da, b.alpha + db));
}

template <std::floating_point T, int N>
Taylor<T, N> abs(const Taylor<T, N>& a)
{
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

/// Sign certainty and enclosures of plain numbers, for Dual<double, N>.
template <std::floating_point T>
int certain_sign(T v)
{
    return (v > T(0)) - (v < T(0));
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

template <typename V, int N>
Dual<V, N> sqrt(const Dual<V, N>& a)
{
    using std::sqrt;
    Dual<V, N> r;
    r.v = sqrt(a.v);
    const V twice = r.v * V(2.0);
    for (int k = 0; k < N; ++k)
        r.d[k] = a.d[k] / twice;
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
