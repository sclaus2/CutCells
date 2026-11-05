// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "analytic.h"

#include <limits>

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

namespace
{
double derivative_value(const double* x, void* context)
{
    const DirectionalDerivative& d = *static_cast<const DirectionalDerivative*>(context);
    double g[3];
    d.phi->gradient(x, g, d.phi->context);
    return g[0] * d.direction[0] + g[1] * d.direction[1] + g[2] * d.direction[2];
}

double derivative_gradient(const double* x, double* grad, void* context)
{
    const DirectionalDerivative& d = *static_cast<const DirectionalDerivative*>(context);
    // the Hessian at x: its bounds over a box too small to widen them
    const double h = 1e-4;
    const double axes[9] = {h, 0, 0, 0, h, 0, 0, 0, h};
    double bounds[18];
    grad[0] = grad[1] = grad[2] = std::numeric_limits<double>::quiet_NaN();
    if (d.phi->hessian_bounds(x, axes, 3, bounds, d.phi->context) != 0)
        for (int i = 0; i < 3; ++i)
        {
            double sum = 0;
            for (int j = 0; j < 3; ++j)
                sum += 0.5 * (bounds[2 * (i * 3 + j)] + bounds[2 * (i * 3 + j) + 1]) / (h * h) * d.direction[j];
            grad[i] = sum;
        }
    return derivative_value(x, context);
}

int derivative_box_bounds(const double* lo, const double* hi, double* b, void* context)
{
    const DirectionalDerivative& d = *static_cast<const DirectionalDerivative*>(context);
    double parent[8];
    if (d.phi->box_bounds(lo, hi, parent, d.phi->context) != 1)
        return 0;
    b[0] = b[1] = 0;
    for (int i = 0; i < 3; ++i)
    {
        const double a0 = parent[2 + 2 * i] * d.direction[i], a1 = parent[3 + 2 * i] * d.direction[i];
        b[0] += std::min(a0, a1);
        b[1] += std::max(a0, a1);
    }
    double centre[3], half[3], axes[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}, bounds[18];
    for (int i = 0; i < 3; ++i)
    {
        centre[i] = 0.5 * (lo[i] + hi[i]);
        half[i] = std::max(0.5 * (hi[i] - lo[i]), 1e-8 * (1 + std::abs(centre[i])));
        axes[i * 3 + i] = half[i];
    }
    if (d.phi->hessian_bounds(centre, axes, 3, bounds, d.phi->context) == 0)
        return 2;
    for (int i = 0; i < 3; ++i)
    {
        double low = 0, high = 0;
        for (int j = 0; j < 3; ++j)
        {
            const double scale = d.direction[j] / (half[i] * half[j]);
            const double a0 = bounds[2 * (i * 3 + j)] * scale, a1 = bounds[2 * (i * 3 + j) + 1] * scale;
            low += std::min(a0, a1);
            high += std::max(a0, a1);
        }
        b[2 + 2 * i] = low;
        b[3 + 2 * i] = high;
    }
    return 1;
}
} // namespace

namespace
{
double shifted_value(const double* x, void* context)
{
    const ShiftedLevelSet& s = *static_cast<const ShiftedLevelSet*>(context);
    return s.phi->value(x, s.phi->context) + s.shift;
}

double shifted_gradient(const double* x, double* grad, void* context)
{
    const ShiftedLevelSet& s = *static_cast<const ShiftedLevelSet*>(context);
    return s.phi->gradient(x, grad, s.phi->context) + s.shift;
}

int shifted_box_bounds(const double* lo, const double* hi, double* b, void* context)
{
    const ShiftedLevelSet& s = *static_cast<const ShiftedLevelSet*>(context);
    const int ok = s.phi->box_bounds(lo, hi, b, s.phi->context);
    if (ok != 0)
    {
        b[0] += s.shift;
        b[1] += s.shift;
    }
    return ok;
}

int shifted_taylor_bounds(const double* centre, const double* axes, int m, double* models, void* context)
{
    const ShiftedLevelSet& s = *static_cast<const ShiftedLevelSet*>(context);
    const int ok = s.phi->taylor_bounds(centre, axes, m, models, s.phi->context);
    if (ok != 0)
        models[0] += s.shift;
    return ok;
}

int shifted_hessian_bounds(const double* centre, const double* axes, int m, double* bounds, void* context)
{
    const ShiftedLevelSet& s = *static_cast<const ShiftedLevelSet*>(context);
    return s.phi->hessian_bounds(centre, axes, m, bounds, s.phi->context);
}
} // namespace

AnalyticLevelSet shifted_level_set(const ShiftedLevelSet& s)
{
    AnalyticLevelSet out;
    out.context = const_cast<ShiftedLevelSet*>(&s);
    out.value = shifted_value;
    out.gradient = shifted_gradient;
    out.box_bounds = shifted_box_bounds;
    out.taylor_bounds = s.phi->taylor_bounds != nullptr ? shifted_taylor_bounds : nullptr;
    out.hessian_bounds = s.phi->hessian_bounds != nullptr ? shifted_hessian_bounds : nullptr;
    return out;
}

AnalyticLevelSet directional_derivative(const DirectionalDerivative& d)
{
    AnalyticLevelSet out;
    out.context = const_cast<DirectionalDerivative*>(&d);
    out.value = derivative_value;
    out.gradient = derivative_gradient;
    out.box_bounds = derivative_box_bounds;
    return out;
}

} // namespace cutcells::quadrays
