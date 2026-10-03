// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "analytic.h"

#include <algorithm>
#include <cmath>

namespace cutcells::quadrays
{

int parallelepiped_bounds(const AnalyticLevelSet& phi, const double* centre, const double* axes, int m,
                          double* models)
{
    if (phi.taylor_bounds != nullptr)
        return phi.taylor_bounds(centre, axes, m, models, phi.context);
    if (phi.box_bounds == nullptr)
        return 0;

    // intervals over the bounding box of the parallelepiped
    double lo[3], hi[3], b[8];
    for (int i = 0; i < 3; ++i)
    {
        double r = 0;
        for (int j = 0; j < m; ++j)
            r += std::abs(axes[i * m + j]);
        lo[i] = centre[i] - r;
        hi[i] = centre[i] + r;
    }
    const int status = phi.box_bounds(lo, hi, b, phi.context);
    if (status == 0)
        return 0;

    const int width = m + 2;
    auto constant = [&](int row, double l, double h)
    {
        double* out = models + width * row;
        out[0] = 0.5 * (l + h);
        for (int j = 0; j < m; ++j)
            out[1 + j] = 0;
        out[m + 1] = 0.5 * (h - l);
    };
    constant(0, b[0], b[1]);
    if (status == 2)
        return 2;
    // d phi / d t_j = sum_i axes_ij d phi / d x_i, in interval arithmetic
    for (int j = 0; j < m; ++j)
    {
        double l = 0, h = 0;
        for (int i = 0; i < 3; ++i)
        {
            const double a = axes[i * m + j], gl = b[2 + 2 * i], gh = b[3 + 2 * i];
            l += std::min(a * gl, a * gh);
            h += std::max(a * gl, a * gh);
        }
        constant(1 + j, l, h);
    }
    return 1;
}

int parallelepiped_hessian(const AnalyticLevelSet& phi, const double* centre, const double* axes, int m,
                           double* bounds)
{
    if (phi.hessian_bounds == nullptr)
        return 0;
    return phi.hessian_bounds(centre, axes, m, bounds, phi.context);
}

} // namespace cutcells::quadrays
