// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <cmath>
#include <concepts>
#include <vector>

namespace cutcells::quadrays
{

/// Gauss-Legendre nodes (ascending) and weights on [0, 1], by Newton's method
/// on the Legendre polynomial in long double.
template <std::floating_point T>
void gauss_legendre(int q, std::vector<T>& x, std::vector<T>& w)
{
    x.resize(q);
    w.resize(q);
    const long double pi = 3.141592653589793238462643383279502884L;
    for (int i = 0; i < q; ++i)
    {
        long double z = std::cos(pi * (i + 0.75L) / (q + 0.5L));
        long double p1 = 0, dp = 0;
        for (int it = 0; it < 100; ++it)
        {
            long double p0 = 1;
            p1 = z;
            for (int j = 2; j <= q; ++j)
            {
                const long double p2 = ((2 * j - 1) * z * p1 - (j - 1) * p0) / j;
                p0 = p1;
                p1 = p2;
            }
            dp = q == 1 ? 1.0L : q * (z * p1 - p0) / (z * z - 1);
            const long double dz = p1 / dp;
            z -= dz;
            if (std::abs(dz) <= 1e-19L)
                break;
        }
        // derivative at the converged node
        long double p0 = 1;
        p1 = z;
        for (int j = 2; j <= q; ++j)
        {
            const long double p2 = ((2 * j - 1) * z * p1 - (j - 1) * p0) / j;
            p0 = p1;
            p1 = p2;
        }
        dp = q == 1 ? 1.0L : q * (z * p1 - p0) / (z * z - 1);
        x[i] = static_cast<T>((1 - z) / 2);
        w[i] = static_cast<T>(1 / ((1 - z * z) * dp * dp));
    }
}

} // namespace cutcells::quadrays
