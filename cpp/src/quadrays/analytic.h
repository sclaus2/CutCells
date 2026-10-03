// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <stdexcept>

#include "taylor.h"

namespace cutcells::quadrays
{

/// An analytic level set as quadrays asks for it, after algoim: values and
/// gradients at points, and bounds of both over regions, in physical
/// coordinates. A plain struct of function pointers and a context pointer,
/// so any library can fill it, from C++ or through a Python capsule, without
/// CutCells' templates.
///
/// box_bounds is required. taylor_bounds is optional: Taylor models over the
/// slanted boxes of simplices are tighter than intervals over their bounding
/// boxes, and save bisections. Both return 1, or 2 if only the value's bounds
/// hold (the derivatives', e.g. of a distance around its centre, are lost), or
/// 0 if no bound holds.
struct AnalyticLevelSet
{
    /// Passed to every function; owned by whoever fills the struct.
    void* context = nullptr;

    /// phi(x), for x with 3 coordinates.
    double (*value)(const double* x, void* context) = nullptr;

    /// phi(x); the gradient goes to grad (3 entries).
    double (*gradient)(const double* x, double* grad, void* context) = nullptr;

    /// Bounds over the box [lo, hi] (3 entries each): phi in [b[0], b[1]] and
    /// d phi / d x_i in [b[2 + 2 i], b[3 + 2 i]].
    int (*box_bounds)(const double* lo, const double* hi, double* b, void* context) = nullptr;

    /// Optional. First-order Taylor models over the parallelepiped
    /// {centre + axes t : t in [-1, 1]^m}, axes 3 x m row-major, 1 <= m <= 3:
    /// models[(m + 2) r .. (m + 2) (r + 1)) = (alpha, beta_0, ..., beta_{m-1}, eps)
    /// for phi (r = 0) and for d phi / d t_j (r = 1 + j), meaning values in
    /// alpha + beta . t + [-eps, eps].
    int (*taylor_bounds)(const double* centre, const double* axes, int m, double* models,
                         void* context) = nullptr;

    /// Optional. Bounds of the second derivatives over the same parallelepiped:
    /// d^2 phi / dt_i dt_j in [bounds[2 (i m + j)], bounds[2 (i m + j) + 1]].
    /// Returns 1, or 0 if no bound holds. Without it, a level set with two
    /// sheets in a cell is not certified with two roots per height line.
    int (*hessian_bounds)(const double* centre, const double* axes, int m, double* bounds,
                          void* context) = nullptr;
};

/// @brief Taylor models over a parallelepiped, in the layout of taylor_bounds:
/// from taylor_bounds if @p phi has it, otherwise constant models from
/// box_bounds over the parallelepiped's bounding box.
/// @return 1, 2 if only the value's model (the first row) holds, 0 if none does
int parallelepiped_bounds(const AnalyticLevelSet& phi, const double* centre, const double* axes, int m,
                          double* models);

/// @brief Bounds of the second derivatives over a parallelepiped, in the layout
/// of hessian_bounds.
/// @return 1, or 0 if @p phi has no hessian_bounds or no bound holds
int parallelepiped_hessian(const AnalyticLevelSet& phi, const double* centre, const double* axes, int m,
                           double* bounds);

// ============================================================================
// Adapter for templated functors
// ============================================================================

/// Taylor models of a functor and of its derivatives along the axes of a
/// parallelepiped, in the layout of AnalyticLevelSet::taylor_bounds.
template <typename F, int M>
int functor_taylor_models(const F& f, const double* centre, const double* axes, double* models)
{
    using TM = Taylor<double, M>;
    std::array<Dual<TM, M>, 3> x;
    for (int i = 0; i < 3; ++i)
    {
        x[i].v = TM(centre[i]);
        for (int j = 0; j < M; ++j)
        {
            x[i].v.beta[j] = axes[i * M + j];
            x[i].d[j] = TM(axes[i * M + j]);
        }
    }
    auto store = [&](const TM& t, int row)
    {
        double* out = models + (M + 2) * row;
        out[0] = t.alpha;
        for (int j = 0; j < M; ++j)
            out[1 + j] = t.beta[j];
        out[M + 1] = t.eps;
    };
    try
    {
        const Dual<TM, M> r = f(x);
        store(r.v, 0);
        for (int j = 0; j < M; ++j)
            store(r.d[j], 1 + j);
        return 1;
    }
    catch (const std::domain_error&)
    {
    }
    // the derivatives have no bound: the value alone
    try
    {
        const std::array<TM, 3> xv = {x[0].v, x[1].v, x[2].v};
        store(f(xv), 0);
        return 2;
    }
    catch (const std::domain_error&)
    {
        return 0;
    }
}

/// Bounds of the second derivatives of a functor along the axes of a
/// parallelepiped, in the layout of AnalyticLevelSet::hessian_bounds.
template <typename F, int M>
int functor_hessian_bounds(const F& f, const double* centre, const double* axes, double* bounds)
{
    using TM = Taylor<double, M>;
    using D1 = Dual<TM, M>;
    std::array<Dual<D1, M>, 3> x;
    for (int i = 0; i < 3; ++i)
    {
        x[i].v.v = TM(centre[i]);
        for (int j = 0; j < M; ++j)
        {
            x[i].v.v.beta[j] = axes[i * M + j];
            x[i].v.d[j] = TM(axes[i * M + j]);
            x[i].d[j].v = TM(axes[i * M + j]);
        }
    }
    try
    {
        const Dual<D1, M> r = f(x);
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < M; ++j)
            {
                const TM& h = r.d[i].d[j];
                const double dev = deviation(h);
                bounds[2 * (i * M + j)] = h.alpha - dev;
                bounds[2 * (i * M + j) + 1] = h.alpha + dev;
            }
        return 1;
    }
    catch (const std::domain_error&)
    {
        return 0;
    }
}

/// Bounds of a functor and of its gradient over a box, in the layout of
/// AnalyticLevelSet::box_bounds.
template <typename F>
int functor_box_bounds(const F& f, const double* lo, const double* hi, double* b)
{
    using TM = Taylor<double, 3>;
    std::array<Dual<TM, 3>, 3> x;
    for (int i = 0; i < 3; ++i)
    {
        x[i].v = TM(0.5 * (lo[i] + hi[i]));
        x[i].v.beta[i] = 0.5 * (hi[i] - lo[i]);
        x[i].d[i] = TM(1.0);
    }
    auto range = [](const TM& t, double* out)
    {
        const double dev = deviation(t);
        out[0] = t.alpha - dev;
        out[1] = t.alpha + dev;
    };
    try
    {
        const Dual<TM, 3> r = f(x);
        range(r.v, b);
        for (int i = 0; i < 3; ++i)
            range(r.d[i], b + 2 + 2 * i);
        return 1;
    }
    catch (const std::domain_error&)
    {
    }
    // the gradient has no bound: the value alone
    try
    {
        const std::array<TM, 3> xv = {x[0].v, x[1].v, x[2].v};
        range(f(xv), b);
        return 2;
    }
    catch (const std::domain_error&)
    {
        return 0;
    }
}

/// @brief The interface of an algoim-style functor
///   template <typename T> T operator()(const std::array<T, 3>& x) const,
/// in physical coordinates. It is evaluated with double for values,
/// Dual<double, 3> for gradients, Dual<Taylor<double, m>, m> for bounds
/// (Taylor<double, m> where the derivatives have none) and
/// Dual<Dual<Taylor<double, m>, m>, m> for second derivatives, so it must call sqrt,
/// exp, log, sin, cos, abs, min and max unqualified (with using std::sqrt and
/// so on). The functor must outlive the result.
template <typename F>
AnalyticLevelSet analytic_level_set(const F& functor)
{
    AnalyticLevelSet phi;
    phi.context = const_cast<void*>(static_cast<const void*>(&functor)); // only read
    phi.value = [](const double* x, void* context) -> double
    { return (*static_cast<const F*>(context))(std::array<double, 3>{x[0], x[1], x[2]}); };
    phi.gradient = [](const double* x, double* grad, void* context) -> double
    {
        std::array<Dual<double, 3>, 3> xd;
        for (int i = 0; i < 3; ++i)
        {
            xd[i].v = x[i];
            xd[i].d[i] = 1.0;
        }
        const Dual<double, 3> r = (*static_cast<const F*>(context))(xd);
        for (int i = 0; i < 3; ++i)
            grad[i] = r.d[i];
        return r.v;
    };
    phi.box_bounds = [](const double* lo, const double* hi, double* b, void* context) -> int
    { return functor_box_bounds(*static_cast<const F*>(context), lo, hi, b); };
    phi.taylor_bounds = [](const double* centre, const double* axes, int m, double* models, void* context) -> int
    {
        const F& f = *static_cast<const F*>(context);
        switch (m)
        {
        case 1:
            return functor_taylor_models<F, 1>(f, centre, axes, models);
        case 2:
            return functor_taylor_models<F, 2>(f, centre, axes, models);
        case 3:
            return functor_taylor_models<F, 3>(f, centre, axes, models);
        default:
            return 0;
        }
    };
    phi.hessian_bounds = [](const double* centre, const double* axes, int m, double* bounds, void* context) -> int
    {
        const F& f = *static_cast<const F*>(context);
        switch (m)
        {
        case 1:
            return functor_hessian_bounds<F, 1>(f, centre, axes, bounds);
        case 2:
            return functor_hessian_bounds<F, 2>(f, centre, axes, bounds);
        case 3:
            return functor_hessian_bounds<F, 3>(f, centre, axes, bounds);
        default:
            return 0;
        }
    };
    return phi;
}

} // namespace cutcells::quadrays
