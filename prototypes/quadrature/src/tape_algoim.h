// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

// The branch helpers and log() that tape.h's evaluate<T> needs for algoim's
// Taylor-model intervals. They live in namespace algoim so that argument-dependent
// lookup finds them. Include only in the translation unit that includes algoim
// (algoim's headers define non-inline functions).

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "interval.hpp"

namespace algoim
{

/// Values an interval can take: alpha -+ maxDeviation.
template <int N>
std::pair<real, real> range(const Interval<N>& a)
{
    const real b = a.maxDeviation();
    return {a.alpha - b, a.alpha + b};
}

/// The constant function with values in [lo, hi].
template <int N>
Interval<N> constant_range(real lo, real hi)
{
    return Interval<N>(0.5 * (lo + hi), uvector<real, N>(0.0), 0.5 * (hi - lo));
}

template <int N>
int certain_sign(const Interval<N>& a)
{
    return a.sign();
}

template <int N>
bool certainly_less(const Interval<N>& a, const Interval<N>& b)
{
    return (a - b).sign() < 0;
}

template <int N>
Interval<N> hull(const Interval<N>& a, const Interval<N>& b)
{
    const auto [al, ah] = range(a);
    const auto [bl, bh] = range(b);
    return constant_range<N>(std::min(al, bl), std::max(ah, bh));
}

template <int N>
Interval<N> tabs(const Interval<N>& a)
{
    const int s = a.sign();
    if (s > 0)
        return a;
    if (s < 0)
        return -a;
    const auto [lo, hi] = range(a);
    return constant_range<N>(0.0, std::max(-lo, hi));
}

template <int N>
Interval<N> tmin(const Interval<N>& a, const Interval<N>& b)
{
    if (certainly_less(a, b))
        return a;
    if (certainly_less(b, a))
        return b;
    const auto [al, ah] = range(a);
    const auto [bl, bh] = range(b);
    return constant_range<N>(std::min(al, bl), std::min(ah, bh));
}

template <int N>
Interval<N> tmax(const Interval<N>& a, const Interval<N>& b)
{
    if (certainly_less(a, b))
        return b;
    if (certainly_less(b, a))
        return a;
    const auto [al, ah] = range(a);
    const auto [bl, bh] = range(b);
    return constant_range<N>(std::max(al, bl), std::max(ah, bh));
}

/// sqrt with the full first-order remainder |f'(alpha)| eps + C/2 b^2, where C bounds
/// |sqrt''| = x^{-3/2} / 4 on [alpha - b, alpha + b]. algoim::sqrt leaves out the
/// first term, so its bound is too tight for a non-linear argument (eps > 0): for
/// |x - c| - r it declared cut cells uncut.
template <int N>
Interval<N> tsqrt(const Interval<N>& i)
{
    const real b = i.maxDeviation();
    if (b >= i.alpha)
        throw std::domain_error("Unable to compute sqrt() with supplied argument");
    const real s = std::sqrt(i.alpha);
    const real lo = i.alpha - b;
    const real C = 0.25 / (lo * std::sqrt(lo));
    return Interval<N>(s, (0.5 / s) * i.beta, (0.5 / s) * i.eps + 0.5 * C * b * b);
}

/// log with the first-order Taylor bound of interval.hpp: |log''| <= 1 / (alpha - b)^2.
template <int N>
Interval<N> log(const Interval<N>& i)
{
    const real b = i.maxDeviation();
    if (b >= i.alpha)
        throw std::domain_error("Unable to compute log() with supplied argument");
    const real inv = 1.0 / i.alpha;
    const real C = 1.0 / ((i.alpha - b) * (i.alpha - b));
    return Interval<N>(std::log(i.alpha), inv * i.beta, inv * i.eps + 0.5 * C * b * b);
}

} // namespace algoim
