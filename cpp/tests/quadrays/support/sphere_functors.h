// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// The test sphere as algoim-style functors of physical coordinates, for the
// analytic interface (quadrays/analytic.h) and for algoim's 2015 engine.
// Shared by the tests in cpp/tests/quadrays and the drivers in benchmarks/quadrays.

#pragma once

#include <array>
#include <cmath>

namespace cutcells::quadrays::support
{

/// |x - c| - r: the signed distance, not a polynomial.
struct SphereDistance
{
    std::array<double, 3> c = {0, 0, 0};
    double r = 0;

    template <typename V>
    V operator()(const std::array<V, 3>& x) const
    {
        using std::sqrt;
        const V dx = x[0] - c[0], dy = x[1] - c[1], dz = x[2] - c[2];
        return sqrt(dx * dx + dy * dy + dz * dz) - r;
    }
};

/// |x - c|^2 - r^2: the polynomial the Bernstein path interpolates exactly.
struct SphereQuadratic
{
    std::array<double, 3> c = {0, 0, 0};
    double r = 0;

    template <typename V>
    V operator()(const std::array<V, 3>& x) const
    {
        const V dx = x[0] - c[0], dy = x[1] - c[1], dz = x[2] - c[2];
        return dx * dx + dy * dy + dz * dz - r * r;
    }
};

} // namespace cutcells::quadrays::support
