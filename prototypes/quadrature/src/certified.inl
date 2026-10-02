// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Included by generators.cpp: algoim's headers define non-inline functions, so
// all code using them must live in one translation unit.

#include "certified.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <map>
#include <vector>

#include "quadrature_multipoly.hpp"

namespace cutcells::proto
{

namespace certify_detail
{
using algoim::real;
using algoim::SparkStack; // needed by the algoim_spark_alloc macro

template <int D>
using VecD = std::array<double, D>;

constexpr double tiny = 1e-14;
constexpr double infinity = std::numeric_limits<double>::infinity();

/// Euclidean norm without underflow or overflow: phi may be scaled by any factor.
template <std::size_t N>
double scaled_norm(const std::array<double, N>& v)
{
    double m = 0;
    for (double x : v)
        m = std::max(m, std::abs(x));
    if (m == 0 || !std::isfinite(m))
        return m;
    double s = 0;
    for (double x : v)
        s += (x / m) * (x / m);
    return m * std::sqrt(s);
}

struct Context
{
    const algoim::xarray<real, 3>* phi = nullptr; ///< Bernstein form of phi on the cell's unit box
    int degree = 2;                               ///< degree of phi in each variable
    int q = 3;
    CertifyOptions opt;
    CertifyStats* stats = nullptr;
    int vis_order = 0;                     ///< > 0: leaf-cell nodes (vis_order + 1 per segment) instead of Gauss points
    mutable std::array<int, 3> box_ids{}; ///< counters of certified boxes, per level
    double phi_scale = 0;                  ///< largest |Bernstein coefficient| of phi on the cell
    mutable int bisections = 0;            ///< bisections so far in this cell (CertifyOptions::max_bisections)
    const Tape* tape = nullptr;            ///< analytic level set: evaluate it instead of phi's Bernstein form
    Vec3 x0{};                             ///< with a tape: the cell's map x = x0 + jac u
    Mat3 jac{};
};

/// psi(y) = phi(A y + b) (curved), or a . y + c (linear), on level coordinates y.
template <int D>
struct Func
{
    std::uint32_t origin = 1; ///< diagnostics: chain of 4-bit codes, see origin_name
    bool linear = false;
    std::array<VecD<D>, 3> A{};
    Vec3 b{};
    VecD<D> a{};
    double c = 0;
};

/// c . y <= d
template <int D>
struct Half
{
    VecD<D> c{};
    double d = 0;
};

template <int D>
struct Problem
{
    VecD<D> lo{}, hi{};
    std::vector<Func<D>> funcs;
    std::vector<Half<D>> clips;
    bool rotated = false; ///< already in the diagonal frame (level 2)
};

/// y_k = alpha + beta . y' on the base coordinates y'
template <int D>
struct Bound
{
    double alpha = 0;
    VecD<D - 1> beta{};
    int kind = 7; ///< diagnostics: 2/3 box face below/above, 4/5 clip below/above, 6 linear root
};

/// Origin codes: 1 phi; restricted to 2 lower box face, 3 upper box face, 4 lower
/// clip, 5 upper clip, 6 linear root, 7 unchanged; new linear functions: 8 switch
/// of lower bounds, 9 switch of upper bounds, 10 difference of linear roots.
inline std::string origin_name(std::uint32_t origin)
{
    static const char* names[] = {"?", "phi", "box-lo", "box-hi", "clip-lo", "clip-up", "lin-root", "same",
                                  "switch-lo", "switch-up", "root-diff"};
    std::vector<std::string> parts;
    for (; origin != 0; origin >>= 4)
        parts.push_back(names[std::min<std::uint32_t>(origin & 15u, 10u)]);
    std::string s;
    for (auto it = parts.rbegin(); it != parts.rend(); ++it)
        s += (s.empty() ? "" : "|") + *it;
    return s;
}

/// Where an emitted point sits in the decomposition. Per level (index = number of
/// free coordinates - 1): the certified box, the segment along its height line and
/// the node within the segment. Leaf cells are assembled from it.
struct Tag
{
    std::array<int, 3> box{{-1, -1, -1}};
    std::array<int, 3> segment{{0, 0, 0}};
    std::array<int, 3> node{{0, 0, 0}};
};

template <int D>
using Emit = std::function<void(const VecD<D>&, double, const Tag&)>;

enum class LineRule
{
    segments, ///< volume: Gauss-Legendre on the segments between bounds and roots
    interface ///< top level only: the root of phi, weighted by |grad phi| / |d_k phi|
};

// ============================================================================
// Evaluation
// ============================================================================

algoim::uvector<real, 3> uv(const Vec3& u) { return algoim::uvector<real, 3>(u[0], u[1], u[2]); }

template <int D>
Vec3 to_u(const Func<D>& f, const VecD<D>& y)
{
    Vec3 u = f.b;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < D; ++j)
            u[i] += f.A[i][j] * y[j];
    return u;
}

/// Physical point of box coordinates u (tape level sets).
Vec3 to_x(const Context& ctx, const Vec3& u)
{
    Vec3 x = ctx.x0;
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            x[i] += ctx.jac[i][k] * u[k];
    return x;
}

/// Gradient of phi with respect to the box coordinates u (tape level sets).
Vec3 tape_gradient_u(const Context& ctx, const Vec3& u)
{
    Vec3 gx;
    value_gradient(*ctx.tape, to_x(ctx, u), gx);
    Vec3 gu = {0, 0, 0};
    for (int i = 0; i < 3; ++i)
        for (int m = 0; m < 3; ++m)
            gu[i] += gx[m] * ctx.jac[m][i];
    return gu;
}

template <int D>
double value(const Context& ctx, const Func<D>& f, const VecD<D>& y)
{
    if (f.linear)
    {
        double s = f.c;
        for (int j = 0; j < D; ++j)
            s += f.a[j] * y[j];
        return s;
    }
    if (ctx.tape)
        return tape_value(*ctx.tape, to_x(ctx, to_u<D>(f, y)));
    return algoim::bernstein::evalBernsteinPoly(*ctx.phi, uv(to_u<D>(f, y)));
}

template <int D>
double derivative(const Context& ctx, const Func<D>& f, const VecD<D>& y, int k)
{
    if (f.linear)
        return f.a[k];
    if (ctx.tape)
    {
        const Vec3 gu = tape_gradient_u(ctx, to_u<D>(f, y));
        double s = 0;
        for (int i = 0; i < 3; ++i)
            s += gu[i] * f.A[i][k];
        return s;
    }
    const auto g = algoim::bernstein::evalBernsteinPolyGradient(*ctx.phi, uv(to_u<D>(f, y)));
    double s = 0;
    for (int i = 0; i < 3; ++i)
        s += g(i) * f.A[i][k];
    return s;
}

template <int D>
VecD<D> insert(const VecD<D - 1>& yb, int k, double t)
{
    VecD<D> y{};
    for (int j = 0, jb = 0; j < D; ++j)
        y[j] = j == k ? t : yb[jb++];
    return y;
}

template <int D>
double bound_value(const Bound<D>& b, const VecD<D - 1>& yb)
{
    double s = b.alpha;
    for (int j = 0; j < D - 1; ++j)
        s += b.beta[j] * yb[j];
    return s;
}

/// f restricted to y_k = alpha + beta . y'
template <int D>
Func<D - 1> restrict_to(const Func<D>& f, int k, const Bound<D>& b)
{
    Func<D - 1> r;
    r.origin = (f.origin << 4) | static_cast<std::uint32_t>(b.kind);
    r.linear = f.linear;
    for (int jb = 0; jb < D - 1; ++jb)
    {
        const int j = jb < k ? jb : jb + 1;
        if (f.linear)
            r.a[jb] = f.a[j] + f.a[k] * b.beta[jb];
        else
            for (int i = 0; i < 3; ++i)
                r.A[i][jb] = f.A[i][j] + f.A[i][k] * b.beta[jb];
    }
    if (f.linear)
        r.c = f.c + f.a[k] * b.alpha;
    else
        for (int i = 0; i < 3; ++i)
            r.b[i] = f.b[i] + f.A[i][k] * b.alpha;
    return r;
}

/// Linear base function b1(y') - b2(y')
template <int D>
Func<D - 1> difference(const Bound<D>& b1, const Bound<D>& b2)
{
    Func<D - 1> r;
    r.origin = b1.kind == 6 ? 10u : (b1.kind == 2 || b1.kind == 4 ? 8u : 9u);
    r.linear = true;
    r.c = b1.alpha - b2.alpha;
    for (int j = 0; j < D - 1; ++j)
        r.a[j] = b1.beta[j] - b2.beta[j];
    return r;
}

// ============================================================================
// Bernstein forms, bounds and roots
// ============================================================================

/// Bernstein coefficients of a curved function on the box [lo, hi]. Its degree in
/// y_j is the degree of phi times the number of coordinates u_i depending on y_j.
template <int D>
void bernstein_form(const Context& ctx, const Func<D>& f, const VecD<D>& lo, const VecD<D>& hi,
                    std::vector<real>& buffer, algoim::uvector<int, D>& ext)
{
    int size = 1;
    for (int j = 0; j < D; ++j)
    {
        int deg = 0;
        for (int i = 0; i < 3; ++i)
            if (std::abs(f.A[i][j]) > tiny)
                deg += ctx.degree;
        ext(j) = std::max(deg, 1) + 1;
        size *= ext(j);
    }
    buffer.assign(size, 0.0);
    algoim::xarray<real, D> p(buffer.data(), ext);
    algoim::bernstein::bernsteinInterpolate<D>(
        [&](const algoim::uvector<real, D>& s)
        {
            VecD<D> y;
            for (int j = 0; j < D; ++j)
                y[j] = lo[j] + s(j) * (hi[j] - lo[j]);
            return value<D>(ctx, f, y);
        },
        p);
}

bool may_vanish(const std::vector<real>& coeffs)
{
    bool pos = false, neg = false;
    for (real v : coeffs)
    {
        pos |= v >= 0;
        neg |= v <= 0;
    }
    return pos && neg;
}

/// For each direction k: certified lower bound of |d_k psi| over the box divided by
/// an upper bound of |grad psi| (0 if d_k psi may change sign).
template <int D>
VecD<D> margins(std::vector<real>& coeffs, const algoim::uvector<int, D>& ext, const VecD<D>& lo, const VecD<D>& hi)
{
    algoim::xarray<real, D> p(coeffs.data(), ext);
    VecD<D> lower{}, upper{};
    std::vector<real> buffer;
    for (int k = 0; k < D; ++k)
    {
        if (ext(k) < 2)
            continue;
        algoim::uvector<int, D> e = ext;
        e(k) -= 1;
        int size = 1;
        for (int j = 0; j < D; ++j)
            size *= e(j);
        buffer.assign(size, 0.0);
        algoim::xarray<real, D> dk(buffer.data(), e);
        algoim::bernstein::bernsteinDerivative(p, k, dk);
        bool pos = true, neg = true;
        double amin = infinity, amax = 0;
        for (int i = 0; i < size; ++i)
        {
            pos &= buffer[i] > 0;
            neg &= buffer[i] < 0;
            amin = std::min(amin, std::abs(buffer[i]));
            amax = std::max(amax, std::abs(buffer[i]));
        }
        const double scale = 1.0 / (hi[k] - lo[k]);
        lower[k] = (pos || neg) ? amin * scale : 0.0;
        upper[k] = amax * scale;
    }
    const double norm = scaled_norm(upper);
    VecD<D> ratio{};
    for (int k = 0; k < D; ++k)
        ratio[k] = norm > 0 ? lower[k] / norm : 0.0;
    return ratio;
}

/// Margins from M^D sub-cells. Only sub-cells that may meet the clipped region and on
/// which psi may vanish count. If d_k psi has one strict sign on all of them, psi has
/// at most one root on every height line in the region (consecutive simple roots have
/// derivatives of opposite sign), and ratio[k] is the smallest local
/// min|d_k psi| / max|grad psi|. Returns false if no sub-cell counts: the function
/// does not vanish in the cell and can be dropped.
template <int D>
bool local_margins(std::vector<real>& coeffs, const algoim::uvector<int, D>& ext, const VecD<D>& lo, const VecD<D>& hi,
                   const std::vector<Half<D>>& clips, int M, VecD<D>& ratio)
{
    algoim::xarray<real, D> p(coeffs.data(), ext);
    int size = 1;
    for (int j = 0; j < D; ++j)
        size *= ext(j);
    std::vector<real> sub_buffer(size), d_buffer;
    algoim::xarray<real, D> sub(sub_buffer.data(), ext);
    VecD<D> local, scale;
    std::array<int, D> sign{};
    local.fill(infinity);
    for (int j = 0; j < D; ++j)
        scale[j] = M / (hi[j] - lo[j]);
    bool relevant = false;
    int cells = 1;
    for (int j = 0; j < D; ++j)
        cells *= M;
    for (int n = 0; n < cells; ++n)
    {
        algoim::uvector<real, D> a, b;
        VecD<D> ylo, yhi;
        for (int j = 0, m = n; j < D; ++j, m /= M)
        {
            a(j) = real(m % M) / M;
            b(j) = real(m % M + 1) / M;
            ylo[j] = lo[j] + (hi[j] - lo[j]) * a(j);
            yhi[j] = lo[j] + (hi[j] - lo[j]) * b(j);
        }
        bool meets = true;
        for (const Half<D>& h : clips)
        {
            double vmin = 0;
            for (int j = 0; j < D; ++j)
                vmin += h.c[j] * (h.c[j] >= 0 ? ylo[j] : yhi[j]);
            meets &= vmin <= h.d + 1e-12;
        }
        if (!meets)
            continue;
        algoim::bernstein::deCasteljau(p, a, b, sub);
        if (!may_vanish(sub_buffer))
            continue;
        relevant = true;
        VecD<D> dmin{}, dmax{};
        std::array<int, D> s{};
        for (int k = 0; k < D; ++k)
        {
            if (ext(k) < 2)
                continue;
            algoim::uvector<int, D> e = ext;
            e(k) -= 1;
            int dsize = 1;
            for (int j = 0; j < D; ++j)
                dsize *= e(j);
            d_buffer.assign(dsize, 0.0);
            algoim::xarray<real, D> dk(d_buffer.data(), e);
            algoim::bernstein::bernsteinDerivative(sub, k, dk);
            bool pos = true, neg = true;
            double amin = infinity, amax = 0;
            for (int i = 0; i < dsize; ++i)
            {
                pos &= d_buffer[i] > 0;
                neg &= d_buffer[i] < 0;
                amin = std::min(amin, std::abs(d_buffer[i]));
                amax = std::max(amax, std::abs(d_buffer[i]));
            }
            s[k] = pos ? 1 : (neg ? -1 : 0);
            dmin[k] = (pos || neg) ? amin * scale[k] : 0.0;
            dmax[k] = amax * scale[k];
        }
        const double norm = scaled_norm(dmax);
        for (int k = 0; k < D; ++k)
        {
            // one strict sign shared by all counted sub-cells
            if (s[k] == 0 || (sign[k] != 0 && sign[k] != s[k]))
                sign[k] = 2;
            else if (sign[k] == 0)
                sign[k] = s[k];
            local[k] = std::min(local[k], norm > 0 ? dmin[k] / norm : 0.0);
        }
    }
    for (int k = 0; k < D; ++k)
        ratio[k] = (relevant && (sign[k] == 1 || sign[k] == -1)) ? local[k] : 0.0;
    return relevant;
}

/// Root of g in (a, b) given a sign change; bisection with secant steps (Illinois).
double bracketed_root(const std::function<double(double)>& g, double a, double b, double ga, double gb)
{
    int side = 0;
    for (int it = 0; it < 200; ++it)
    {
        double x = (a * gb - b * ga) / (gb - ga);
        if (!(x > a && x < b))
            x = 0.5 * (a + b);
        const double gx = g(x);
        if (gx == 0.0 || b - a <= 4e-16 * std::max(1.0, std::abs(x)))
            return x;
        if ((gx > 0) == (gb > 0))
        {
            b = x;
            gb = gx;
            if (side == -1)
                ga *= 0.5;
            side = -1;
        }
        else
        {
            a = x;
            ga = gx;
            if (side == 1)
                gb *= 0.5;
            side = 1;
        }
    }
    return 0.5 * (a + b);
}

/// All roots in (a, b) of the polynomial with Bernstein coefficients c on [a, b],
/// isolated by de Casteljau subdivision and Descartes' rule of signs.
void isolate_roots(const std::vector<double>& c, double a, double b, const std::function<double(double)>& g,
                   std::vector<double>& out, int depth = 0)
{
    int changes = 0, last = 0;
    for (double v : c)
    {
        const int s = (v > 0) - (v < 0);
        if (s == 0)
            continue;
        if (last != 0 && s != last)
            ++changes;
        last = s;
    }
    if (changes == 0)
        return;
    if (changes == 1 || depth >= 48)
    {
        const double ga = g(a), gb = g(b);
        if (ga != 0.0 && gb != 0.0 && (ga > 0) != (gb > 0))
            out.push_back(bracketed_root(g, a, b, ga, gb));
        return;
    }
    const int n = static_cast<int>(c.size());
    std::vector<double> left(n), right(n), work = c;
    for (int r = 0; r < n; ++r)
    {
        left[r] = work[0];
        right[n - 1 - r] = work[n - 1 - r];
        for (int i = 0; i < n - 1 - r; ++i)
            work[i] = 0.5 * (work[i] + work[i + 1]);
    }
    const double mid = 0.5 * (a + b);
    isolate_roots(left, a, mid, g, out, depth + 1);
    isolate_roots(right, mid, b, g, out, depth + 1);
}

/// All roots in (a, b) of a curved function along a line, of degree `deg` in t.
void line_roots(const Context& ctx, int deg, double a, double b, const std::function<double(double)>& g,
                std::vector<double>& out)
{
    deg = std::max(deg, 1);
    std::vector<real> c(deg + 1);
    algoim::xarray<real, 1> p(c.data(), algoim::uvector<int, 1>(deg + 1));
    algoim::bernstein::bernsteinInterpolate<1>([&](const algoim::uvector<real, 1>& s) { return g(a + s(0) * (b - a)); },
                                               p);
    (void)ctx;
    isolate_roots(c, a, b, g, out);
}

// ============================================================================
// Tape level sets: bounds from first-order Taylor models
// ============================================================================

/// Bounds of psi(y) = phi(x0 + jac (A y + b)) and of its derivatives on [lo, hi],
/// from algoim::Interval<D> (value at the centre, gradient, remainder): the affine
/// map from the box to physical space is represented exactly. Sets may_vanish, the
/// direction margins (as margins() does from Bernstein coefficients) and the
/// largest |psi|. A lost bound (sqrt or division near zero) makes nothing certain.
template <int D>
void tape_bounds(const Context& ctx, const Func<D>& f, const VecD<D>& lo, const VecD<D>& hi, bool& may_vanish,
                 VecD<D>& ratio, double& magnitude)
{
    using TM = algoim::Interval<D>;
    for (int j = 0; j < D; ++j)
        TM::delta(j) = 0.5 * (hi[j] - lo[j]);
    std::array<Dual<TM, D>, 3> x;
    for (int i = 0; i < 3; ++i)
    {
        double centre = ctx.x0[i];
        for (int a = 0; a < 3; ++a)
            centre += ctx.jac[i][a] * f.b[a];
        algoim::uvector<real, D> beta;
        for (int j = 0; j < D; ++j)
        {
            double m = 0;
            for (int a = 0; a < 3; ++a)
                m += ctx.jac[i][a] * f.A[a][j];
            beta(j) = m;
            centre += m * 0.5 * (lo[j] + hi[j]);
        }
        x[i].v = TM(centre, beta);
        for (int j = 0; j < D; ++j)
            x[i].d[j] = TM(beta(j));
    }
    thread_local std::vector<Dual<TM, D>> regs;
    try
    {
        const Dual<TM, D> r = evaluate(*ctx.tape, x, regs);
        may_vanish = r.v.sign() == 0;
        magnitude = std::abs(r.v.alpha) + r.v.maxDeviation();
        VecD<D> lower{}, upper{};
        for (int k = 0; k < D; ++k)
        {
            const double dev = r.d[k].maxDeviation();
            upper[k] = std::abs(r.d[k].alpha) + dev;
            lower[k] = r.d[k].sign() != 0 ? std::abs(r.d[k].alpha) - dev : 0.0;
        }
        const double norm = scaled_norm(upper);
        for (int k = 0; k < D; ++k)
            ratio[k] = norm > 0 ? lower[k] / norm : 0.0;
    }
    catch (const std::domain_error&)
    {
        may_vanish = true;
        ratio.fill(0.0);
        magnitude = infinity;
    }
}

/// Roots of t -> phi(p + t v) on (a, b), isolated with Taylor models: none where
/// the value has a certain sign, at most one where the derivative has; bisection
/// otherwise, down to 2^-40 of the segment.
void tape_line_roots(const Context& ctx, const Vec3& p, const Vec3& v, double a, double b, std::vector<double>& out,
                     int depth = 0)
{
    using TM = algoim::Interval<1>;
    TM::delta(0) = 0.5 * (b - a);
    const double c = 0.5 * (a + b);
    std::array<Dual<TM, 1>, 3> x;
    for (int i = 0; i < 3; ++i)
    {
        x[i].v = TM(p[i] + v[i] * c, algoim::uvector<real, 1>(v[i]));
        x[i].d[0] = TM(v[i]);
    }
    thread_local std::vector<Dual<TM, 1>> regs;
    int value_sign = 0, slope_sign = 0;
    try
    {
        const Dual<TM, 1> r = evaluate(*ctx.tape, x, regs);
        value_sign = r.v.sign();
        slope_sign = r.d[0].sign();
    }
    catch (const std::domain_error&)
    {
    }
    if (value_sign != 0)
        return;
    const auto g = [&](double t) { return tape_value(*ctx.tape, {p[0] + t * v[0], p[1] + t * v[1], p[2] + t * v[2]}); };
    if (slope_sign != 0 || depth >= 40)
    {
        const double ga = g(a), gb = g(b);
        if (ga != 0.0 && gb != 0.0 && (ga > 0) != (gb > 0))
            out.push_back(bracketed_root(g, a, b, ga, gb));
        return;
    }
    const double m = 0.5 * (a + b);
    tape_line_roots(ctx, p, v, a, m, out, depth + 1);
    tape_line_roots(ctx, p, v, m, b, out, depth + 1);
}

/// The line y = y0 + t e_k of a level as p + t v in physical space.
template <int D>
void physical_line(const Context& ctx, const Func<D>& f, const VecD<D>& y0, int k, Vec3& p, Vec3& v)
{
    p = to_x(ctx, to_u<D>(f, y0));
    for (int i = 0; i < 3; ++i)
    {
        v[i] = 0;
        for (int a = 0; a < 3; ++a)
            v[i] += ctx.jac[i][a] * f.A[a][k];
    }
}


// ============================================================================
// Clipped boxes
// ============================================================================

template <int D>
double determinant(const std::array<VecD<D>, D>& m)
{
    if constexpr (D == 1)
        return m[0][0];
    else if constexpr (D == 2)
        return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    else
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
               + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

/// Vertices of {y in [lo, hi] : clips} by enumeration of plane triples (D <= 3).
template <int D>
std::vector<VecD<D>> polytope_vertices(const VecD<D>& lo, const VecD<D>& hi, const std::vector<Half<D>>& clips)
{
    std::vector<VecD<D>> vertices;
    std::vector<Half<D>> planes;
    for (int j = 0; j < D; ++j)
    {
        Half<D> l, u;
        l.c[j] = -1.0;
        l.d = -lo[j];
        u.c[j] = 1.0;
        u.d = hi[j];
        planes.push_back(l);
        planes.push_back(u);
    }
    planes.insert(planes.end(), clips.begin(), clips.end());
    const int np = static_cast<int>(planes.size());
    std::array<int, D> idx{};
    std::function<void(int, int)> choose = [&](int start, int level)
    {
        if (level == D)
        {
            std::array<VecD<D>, D> m;
            for (int r = 0; r < D; ++r)
                m[r] = planes[idx[r]].c;
            const double det = determinant<D>(m);
            if (std::abs(det) < 1e-13)
                return;
            VecD<D> x;
            for (int col = 0; col < D; ++col)
            {
                auto mc = m;
                for (int r = 0; r < D; ++r)
                    mc[r][col] = planes[idx[r]].d;
                x[col] = determinant<D>(mc) / det;
            }
            for (const Half<D>& h : planes)
            {
                double s = 0;
                for (int j = 0; j < D; ++j)
                    s += h.c[j] * x[j];
                if (s > h.d + 1e-11)
                    return;
            }
            vertices.push_back(x);
            return;
        }
        for (int i = start; i < np; ++i)
        {
            idx[level] = i;
            choose(i + 1, level + 1);
        }
    };
    choose(0, 0);
    return vertices;
}

/// Shrink [lo, hi] to the bounding box of {y in [lo, hi] : clips}; false if empty.
template <int D>
bool tighten(VecD<D>& lo, VecD<D>& hi, const std::vector<Half<D>>& clips)
{
    if (clips.empty())
        return true;
    const std::vector<VecD<D>> vertices = polytope_vertices<D>(lo, hi, clips);
    if (vertices.empty())
        return false;
    VecD<D> nlo, nhi;
    nlo.fill(infinity);
    nhi.fill(-infinity);
    for (const VecD<D>& x : vertices)
        for (int j = 0; j < D; ++j)
        {
            nlo[j] = std::min(nlo[j], x[j]);
            nhi[j] = std::max(nhi[j], x[j]);
        }
    for (int j = 0; j < D; ++j)
    {
        lo[j] = std::max(lo[j], nlo[j]);
        hi[j] = std::min(hi[j], nhi[j]);
        if (!(hi[j] - lo[j] > 1e-13))
            return false;
    }
    return true;
}

template <int D>
bool linear_may_vanish(const Func<D>& f, const VecD<D>& lo, const VecD<D>& hi)
{
    double vmin = f.c, vmax = f.c;
    for (int j = 0; j < D; ++j)
    {
        vmin += std::min(f.a[j] * lo[j], f.a[j] * hi[j]);
        vmax += std::max(f.a[j] * lo[j], f.a[j] * hi[j]);
    }
    return vmin <= 0 && vmax >= 0;
}

// ============================================================================
// Integration by dimension reduction
// ============================================================================

/// Why did certification fail? Sample the clipped region on a 9^D grid and look at
/// the function with the smallest margin for the chosen axis near its zero set.
template <int D>
void diagnose_failure(const Context& ctx, const Problem<D>& p, const std::vector<int>& curved,
                      const std::vector<VecD<D>>& ratios, int k, double best, int depth)
{
    std::size_t blocker = 0;
    for (std::size_t i = 1; i < ratios.size(); ++i)
        if (ratios[i][k] < ratios[blocker][k])
            blocker = i;
    const Func<D>& f = p.funcs[curved[blocker]];
    const int G = 9;
    std::vector<VecD<D>> points;
    std::array<int, D> idx{};
    const int total = static_cast<int>(std::pow(G, D));
    for (int n = 0; n < total; ++n)
    {
        int m = n;
        VecD<D> y;
        for (int j = 0; j < D; ++j)
        {
            idx[j] = m % G;
            m /= G;
            y[j] = p.lo[j] + (p.hi[j] - p.lo[j]) * (idx[j] + 0.5) / G;
        }
        bool inside = true;
        for (const Half<D>& h : p.clips)
        {
            double sum = 0;
            for (int j = 0; j < D; ++j)
                sum += h.c[j] * y[j];
            inside &= sum <= h.d;
        }
        if (inside)
            points.push_back(y);
    }
    std::vector<double> values(points.size());
    double vmax = 0;
    bool pos = false, neg = false;
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        values[i] = value<D>(ctx, f, points[i]);
        vmax = std::max(vmax, std::abs(values[i]));
        pos |= values[i] > 0;
        neg |= values[i] < 0;
    }
    std::string category;
    if (!(pos && neg))
        category = "zero set outside the cell";
    else
    {
        VecD<D> rmin;
        rmin.fill(infinity);
        for (std::size_t i = 0; i < points.size(); ++i)
        {
            if (std::abs(values[i]) > 0.2 * vmax)
                continue;
            VecD<D> g;
            for (int j = 0; j < D; ++j)
                g[j] = derivative<D>(ctx, f, points[i], j);
            const double norm = scaled_norm(g);
            for (int j = 0; j < D && norm > 0; ++j)
                rmin[j] = std::min(rmin[j], std::abs(g[j]) / norm);
        }
        double other = 0;
        for (int j = 0; j < D; ++j)
            other = std::max(other, rmin[j]);
        if (rmin[k] >= ctx.opt.margin)
            category = "margin holds near the zero set in the cell";
        else if (other >= ctx.opt.margin)
            category = "another axis has the margin there";
        else
            category = "margin fails near the zero set in the cell";
    }
    const std::string key = "level " + std::to_string(D) + " | " + origin_name(f.origin) + " | "
                            + (best > 0 ? "bound below margin" : "no certain sign") + " | " + category;
    ++ctx.stats->causes[key];
    ++ctx.stats->causes["level " + std::to_string(D) + " | depth " + std::to_string(depth)];

    // From the certified bounds alone: does one function fail on every axis, or does
    // each function have an axis but no axis suits all (a conflict)?
    std::string axes;
    for (std::size_t i = 0; i < ratios.size() && axes.empty(); ++i)
    {
        double top = 0;
        for (int j = 0; j < D; ++j)
            top = std::max(top, ratios[i][j]);
        if (top < ctx.opt.margin)
            axes = "one function fails on every axis: " + origin_name(p.funcs[curved[i]].origin);
    }
    if (axes.empty())
    {
        std::vector<std::string> names;
        for (int j = 0; j < D; ++j)
        {
            std::size_t b = 0;
            for (std::size_t i = 1; i < ratios.size(); ++i)
                if (ratios[i][j] < ratios[b][j])
                    b = i;
            names.push_back(origin_name(p.funcs[curved[b]].origin));
        }
        std::sort(names.begin(), names.end());
        names.erase(std::unique(names.begin(), names.end()), names.end());
        axes = "conflict:";
        for (const std::string& n : names)
            axes += " " + n;
    }
    // is the box cut by the boundary of the clipped region?
    bool cut = false;
    for (const Half<D>& h : p.clips)
    {
        double vmin = 0, vmax = 0;
        for (int j = 0; j < D; ++j)
        {
            vmin += h.c[j] * (h.c[j] >= 0 ? p.lo[j] : p.hi[j]);
            vmax += h.c[j] * (h.c[j] >= 0 ? p.hi[j] : p.lo[j]);
        }
        cut |= vmin < h.d - 1e-12 && vmax > h.d + 1e-12;
    }
    ++ctx.stats->causes["axes | level " + std::to_string(D) + (depth >= ctx.opt.max_depth ? " | at depth limit" : " | before limit")
                        + " | " + axes + (cut ? " | box cut by a face" : " | box inside")];
}

/// Nodes on a segment [a, b]: Gauss-Legendre points and weights for quadrature, or
/// vis_order + 1 equispaced points, pulled 1e-5 inside the ends, for leaf cells (closer
/// to a breakpoint, a root may fall on either side of a bound by rounding).
template <typename F>
void segment_nodes(const Context& ctx, double a, double b, F&& f)
{
    if (ctx.vis_order > 0)
    {
        const double eps = 1e-5;
        for (int j = 0; j <= ctx.vis_order; ++j)
            f(j, a + (b - a) * (eps + (1.0 - 2.0 * eps) * j / ctx.vis_order), 0.0);
        return;
    }
    for (int j = 0; j < ctx.q; ++j)
        f(j, a + (b - a) * algoim::GaussQuad::x(ctx.q, j), (b - a) * algoim::GaussQuad::w(ctx.q, j));
}

template <int D>
void integrate(const Context& ctx, Problem<D> p, const Emit<D>& emit, int depth, LineRule rule);

/// One free coordinate: Gauss-Legendre on the segments between bounds and roots.
template <>
void integrate<1>(const Context& ctx, Problem<1> p, const Emit<1>& emit, int, LineRule)
{
    double L = p.lo[0], U = p.hi[0];
    for (const Half<1>& h : p.clips)
    {
        if (h.c[0] > tiny)
            U = std::min(U, h.d / h.c[0]);
        else if (h.c[0] < -tiny)
            L = std::max(L, h.d / h.c[0]);
        else if (h.d < 0)
            return;
    }
    if (!(U > L))
        return;
    std::vector<double> nodes = {L, U};
    for (const Func<1>& f : p.funcs)
    {
        if (f.linear)
        {
            if (std::abs(f.a[0]) > tiny)
            {
                const double t = -f.c / f.a[0];
                if (t > L && t < U)
                    nodes.push_back(t);
            }
            continue;
        }
        int deg = 0;
        for (int i = 0; i < 3; ++i)
            if (std::abs(f.A[i][0]) > tiny)
                deg += ctx.degree;
        if (ctx.tape)
        {
            Vec3 p0, v;
            physical_line<1>(ctx, f, {0.0}, 0, p0, v);
            tape_line_roots(ctx, p0, v, L, U, nodes);
            continue;
        }
        line_roots(ctx, deg, L, U, [&](double t) { return value<1>(ctx, f, {t}); }, nodes);
    }
    std::sort(nodes.begin(), nodes.end());
    Tag tag;
    tag.box[0] = ctx.box_ids[0]++;
    int segment = 0;
    for (std::size_t s = 0; s + 1 < nodes.size(); ++s)
    {
        const double a = nodes[s], b = nodes[s + 1];
        if (b - a <= 1e-15)
            continue;
        tag.segment[0] = segment++;
        segment_nodes(ctx, a, b,
                      [&](int j, double t, double w)
                      {
                          tag.node[0] = j;
                          emit({t}, w, tag);
                      });
    }
}

/// The functions of a (tightened) problem that may vanish in it, their direction
/// margins and the best height direction.
template <int D>
struct Analysis
{
    std::vector<Func<D>> funcs;
    std::vector<VecD<D>> ratios;
    std::vector<int> curved; ///< index in funcs of each entry of ratios
    int k = 0;
    double best = -1;
    bool certified = false;
};

template <int D>
Analysis<D> analyse(const Context& ctx, const Problem<D>& p)
{
    Analysis<D> an;
    std::vector<Func<D>>& funcs = an.funcs;
    std::vector<VecD<D>>& ratios = an.ratios;
    std::vector<int>& curved = an.curved;
    std::vector<real> coeffs;
    algoim::uvector<int, D> ext;
    for (const Func<D>& f : p.funcs)
    {
        if (f.linear)
        {
            if (linear_may_vanish<D>(f, p.lo, p.hi))
                funcs.push_back(f);
            continue;
        }
        if (ctx.tape)
        {
            bool may = false;
            double magnitude = 0;
            VecD<D> ratio{};
            tape_bounds<D>(ctx, f, p.lo, p.hi, may, ratio, magnitude);
            if (!may || magnitude <= 1e-12 * ctx.phi_scale)
                continue; // no zero in the box, or phi = 0 on a face
            curved.push_back(static_cast<int>(funcs.size()));
            funcs.push_back(f);
            ratios.push_back(ratio);
            continue;
        }
        bernstein_form<D>(ctx, f, p.lo, p.hi, coeffs, ext);
        // phi vanishes identically here (phi = 0 on a face, say): only rounding noise is
        // left, which would block certification; it has no root inside the region
        double cmax = 0;
        for (real v : coeffs)
            cmax = std::max(cmax, std::abs(v));
        if (cmax <= 1e-12 * ctx.phi_scale)
            continue;
        VecD<D> ratio{};
        if (ctx.opt.mask_subdivisions > 1)
        {
            if (!local_margins<D>(coeffs, ext, p.lo, p.hi, p.clips, ctx.opt.mask_subdivisions, ratio))
                continue; // no zero of this function inside the cell
        }
        else
        {
            if (!may_vanish(coeffs))
                continue;
            ratio = margins<D>(coeffs, ext, p.lo, p.hi);
        }
        curved.push_back(static_cast<int>(funcs.size()));
        funcs.push_back(f);
        ratios.push_back(ratio);
    }
    // height direction: the best certified margin over all curved functions
    for (int kk = 0; kk < D; ++kk)
    {
        double score = 1.0;
        for (const VecD<D>& r : ratios)
            score = std::min(score, r[kk]);
        if (score > an.best || (score == an.best && p.hi[kk] - p.lo[kk] > p.hi[an.k] - p.lo[an.k]))
        {
            an.best = score;
            an.k = kk;
        }
    }
    an.certified = ratios.empty() || (an.best > 0 && an.best >= ctx.opt.margin);
    return an;
}

/// A level-2 problem in the diagonal frame y = T z, T = [[1/2, -1/2], [1/2, 1/2]]:
/// z_0 runs along (1, 1) and z_1 along (-1, 1). The box becomes four clips.
inline Problem<2> rotate_diagonal(const Problem<2>& p)
{
    auto transpose_times = [](const VecD<2>& c) { return VecD<2>{0.5 * (c[0] + c[1]), 0.5 * (c[1] - c[0])}; };
    Problem<2> r;
    r.rotated = true;
    for (const Half<2>& h : p.clips)
        r.clips.push_back({transpose_times(h.c), h.d});
    for (int j = 0; j < 2; ++j)
    {
        VecD<2> e{};
        e[j] = 1.0;
        r.clips.push_back({transpose_times(e), p.hi[j]});
        e[j] = -1.0;
        r.clips.push_back({transpose_times(e), -p.lo[j]});
    }
    // z_0 = y_0 + y_1, z_1 = y_1 - y_0
    r.lo = {p.lo[0] + p.lo[1], p.lo[1] - p.hi[0]};
    r.hi = {p.hi[0] + p.hi[1], p.hi[1] - p.lo[0]};
    for (const Func<2>& f : p.funcs)
    {
        Func<2> g = f;
        if (f.linear)
            g.a = transpose_times(f.a);
        else
            for (int i = 0; i < 3; ++i)
                g.A[i] = transpose_times(f.A[i]);
        r.funcs.push_back(g);
    }
    return r;
}

template <int D>
void integrate(const Context& ctx, Problem<D> p, const Emit<D>& emit, int depth, LineRule rule)
{
    if (!tighten<D>(p.lo, p.hi, p.clips))
        return;

    // keep the functions that may vanish on the box, with their direction margins
    Analysis<D> an = analyse<D>(ctx, p);
    if (rule == LineRule::interface && an.funcs.empty())
        return; // the level set does not cross this box
    p.funcs = an.funcs;
    const std::vector<VecD<D>>& ratios = an.ratios;
    const std::vector<int>& curved = an.curved;
    const int k = an.k;
    const double best = an.best;
    const bool certified = an.certified;

    // Level 2: two functions whose zero curves meet (on a face of a tet, say) may each
    // need a different axis, and no bisection separates them. A diagonal frame often
    // suits both.
    if constexpr (D == 2)
    {
        if (!certified && ctx.opt.diagonal_frames && !p.rotated)
        {
            Problem<2> r = rotate_diagonal(p);
            if (tighten<2>(r.lo, r.hi, r.clips) && analyse<2>(ctx, r).certified)
            {
                ++ctx.stats->rotations;
                const Emit<2> back = [&emit](const VecD<2>& z, double w, const Tag& tag)
                { emit({0.5 * (z[0] - z[1]), 0.5 * (z[0] + z[1])}, 0.5 * w, tag); };
                integrate<2>(ctx, std::move(r), back, depth, rule);
                return;
            }
        }
    }
    static const bool trace = std::getenv("CERTIFY_TRACE") != nullptr;
    if (trace)
    {
        std::fprintf(stderr, "%*slevel %d depth %d box", 2 * (3 - D), "", D, depth);
        for (int j = 0; j < D; ++j)
            std::fprintf(stderr, " [%.4f, %.4f]", p.lo[j], p.hi[j]);
        std::fprintf(stderr, " k %d %s, %zu clips\n", k, certified ? "certified" : "FAILED", p.clips.size());
        for (std::size_t i = 0; i < ratios.size(); ++i)
        {
            std::fprintf(stderr, "%*s  %-22s ratios", 2 * (3 - D), "", origin_name(p.funcs[curved[i]].origin).c_str());
            for (int j = 0; j < D; ++j)
                std::fprintf(stderr, " %.3f", ratios[i][j]);
            std::fprintf(stderr, "\n");
        }
        for (const Func<D>& f : p.funcs)
            if (f.linear)
                std::fprintf(stderr, "%*s  %-22s linear\n", 2 * (3 - D), "", origin_name(f.origin).c_str());
    }
    if (!certified && ctx.opt.diagnose)
        diagnose_failure<D>(ctx, p, curved, ratios, k, best, depth);
    if (!certified)
    {
        if (depth < ctx.opt.max_depth && ctx.bisections < ctx.opt.max_bisections)
        {
            ++ctx.bisections;
            ++ctx.stats->bisections;
            int axis = 0;
            for (int j = 1; j < D; ++j)
                if (p.hi[j] - p.lo[j] > p.hi[axis] - p.lo[axis])
                    axis = j;
            const double mid = 0.5 * (p.lo[axis] + p.hi[axis]);
            Problem<D> first = p, second = p;
            first.hi[axis] = mid;
            second.lo[axis] = mid;
            integrate<D>(ctx, std::move(first), emit, depth + 1, rule);
            integrate<D>(ctx, std::move(second), emit, depth + 1, rule);
            return;
        }
        ++ctx.stats->uncertified;
    }

    // bounds of the height lines: box faces and clip planes
    std::vector<Bound<D>> lowers, uppers;
    Problem<D - 1> base;
    for (int jb = 0; jb < D - 1; ++jb)
    {
        const int j = jb < k ? jb : jb + 1;
        base.lo[jb] = p.lo[j];
        base.hi[jb] = p.hi[j];
    }
    {
        Bound<D> l, u;
        l.alpha = p.lo[k];
        u.alpha = p.hi[k];
        l.kind = 2;
        u.kind = 3;
        lowers.push_back(l);
        uppers.push_back(u);
    }
    for (const Half<D>& h : p.clips)
    {
        const double ck = h.c[k];
        if (std::abs(ck) <= tiny)
        {
            Half<D - 1> hb;
            for (int jb = 0; jb < D - 1; ++jb)
                hb.c[jb] = h.c[jb < k ? jb : jb + 1];
            hb.d = h.d;
            base.clips.push_back(hb);
            continue;
        }
        Bound<D> b;
        b.alpha = h.d / ck;
        b.kind = ck > 0 ? 5 : 4;
        for (int jb = 0; jb < D - 1; ++jb)
            b.beta[jb] = -h.c[jb < k ? jb : jb + 1] / ck;
        (ck > 0 ? uppers : lowers).push_back(b);
    }
    // Fourier-Motzkin: the height line is non-empty where every lower <= every upper
    for (std::size_t i = 0; i < lowers.size(); ++i)
        for (std::size_t j = 0; j < uppers.size(); ++j)
        {
            if (i == 0 && j == 0)
                continue;
            Half<D - 1> hb;
            for (int jb = 0; jb < D - 1; ++jb)
                hb.c[jb] = lowers[i].beta[jb] - uppers[j].beta[jb];
            hb.d = uppers[j].alpha - lowers[i].alpha;
            base.clips.push_back(hb);
        }

    // Drop bounds that are never active on the base region: a lower bound below
    // another lower bound at every vertex of the base polytope (bounds are affine),
    // or an upper bound above another upper bound. Their restrictions and switch
    // functions would only cause needless bisections.
    if (ctx.opt.prune_bounds && lowers.size() + uppers.size() > 2)
    {
        const std::vector<VecD<D - 1>> corners = polytope_vertices<D - 1>(base.lo, base.hi, base.clips);
        auto prune = [&](std::vector<Bound<D>>& list, bool lower)
        {
            std::vector<Bound<D>> kept;
            for (std::size_t i = 0; i < list.size(); ++i)
            {
                bool dominated = false;
                for (std::size_t j = 0; j < list.size() && !dominated; ++j)
                {
                    if (j == i)
                        continue;
                    bool all = !corners.empty();
                    for (const VecD<D - 1>& y : corners)
                    {
                        const double di = bound_value<D>(list[i], y), dj = bound_value<D>(list[j], y);
                        // j at least as tight everywhere; ties keep the earlier bound
                        const bool tighter = lower ? (dj > di + 1e-12 || (std::abs(dj - di) <= 1e-12 && j < i))
                                                   : (dj < di - 1e-12 || (std::abs(dj - di) <= 1e-12 && j < i));
                        all &= tighter || std::abs(dj - di) <= 1e-12;
                        if (!all)
                            break;
                    }
                    // all corners tied: keep the earlier one
                    if (all)
                    {
                        bool strictly = false;
                        for (const VecD<D - 1>& y : corners)
                        {
                            const double di = bound_value<D>(list[i], y), dj = bound_value<D>(list[j], y);
                            strictly |= lower ? dj > di + 1e-12 : dj < di - 1e-12;
                        }
                        dominated = strictly || j < i;
                    }
                }
                if (!dominated)
                    kept.push_back(list[i]);
            }
            list = kept;
        };
        prune(lowers, true);
        prune(uppers, false);
    }

    // base functions: where the integrand along the line changes form
    std::vector<Bound<D>> bounds = lowers;
    bounds.insert(bounds.end(), uppers.begin(), uppers.end());
    std::vector<Bound<D>> linear_roots;
    for (const Func<D>& f : p.funcs)
        if (f.linear && std::abs(f.a[k]) > tiny)
        {
            Bound<D> r;
            r.kind = 6;
            r.alpha = -f.c / f.a[k];
            for (int jb = 0; jb < D - 1; ++jb)
                r.beta[jb] = -f.a[jb < k ? jb : jb + 1] / f.a[k];
            linear_roots.push_back(r);
        }
    for (const Func<D>& f : p.funcs)
    {
        if (f.linear && std::abs(f.a[k]) <= tiny)
        {
            base.funcs.push_back(restrict_to<D>(f, k, Bound<D>{})); // independent of y_k
            continue;
        }
        for (const Bound<D>& b : bounds) // a root reaches a bound
            base.funcs.push_back(restrict_to<D>(f, k, b));
        if (!f.linear)
            for (const Bound<D>& r : linear_roots) // a curved root meets a linear one
                base.funcs.push_back(restrict_to<D>(f, k, r));
    }
    for (std::size_t i = 0; i < linear_roots.size(); ++i)
        for (std::size_t j = i + 1; j < linear_roots.size(); ++j)
            base.funcs.push_back(difference<D>(linear_roots[i], linear_roots[j]));
    for (std::size_t i = 0; i < lowers.size(); ++i) // the active lower bound changes
        for (std::size_t j = i + 1; j < lowers.size(); ++j)
            base.funcs.push_back(difference<D>(lowers[i], lowers[j]));
    for (std::size_t i = 0; i < uppers.size(); ++i)
        for (std::size_t j = i + 1; j < uppers.size(); ++j)
            base.funcs.push_back(difference<D>(uppers[i], uppers[j]));

    // integrand of the base: the rule along the height line
    const int box_id = ctx.box_ids[D - 1]++;
    const Emit<D - 1> line = [&, k, certified, box_id](const VecD<D - 1>& yb, double w, const Tag& below)
    {
        double L = -infinity, U = infinity;
        for (const Bound<D>& b : lowers)
            L = std::max(L, bound_value<D>(b, yb));
        for (const Bound<D>& b : uppers)
            U = std::min(U, bound_value<D>(b, yb));
        if (!(U > L))
            return;
        std::vector<double> nodes;
        for (const Func<D>& f : p.funcs)
        {
            if (f.linear)
            {
                if (std::abs(f.a[k]) > tiny)
                {
                    double t = f.c;
                    for (int jb = 0; jb < D - 1; ++jb)
                        t += f.a[jb < k ? jb : jb + 1] * yb[jb];
                    t = -t / f.a[k];
                    if (t > L && t < U)
                        nodes.push_back(t);
                }
                continue;
            }
            const auto g = [&](double t) { return value<D>(ctx, f, insert<D>(yb, k, t)); };
            if (certified)
            {
                const double gl = g(L), gu = g(U);
                if (gl != 0.0 && gu != 0.0 && (gl > 0) != (gu > 0))
                    nodes.push_back(bracketed_root(g, L, U, gl, gu));
            }
            else
            {
                if (ctx.tape)
                {
                    Vec3 p0, v;
                    physical_line<D>(ctx, f, insert<D>(yb, k, 0.0), k, p0, v);
                    tape_line_roots(ctx, p0, v, L, U, nodes);
                    continue;
                }
                int deg = 0;
                for (int i = 0; i < 3; ++i)
                    if (std::abs(f.A[i][k]) > tiny)
                        deg += ctx.degree;
                line_roots(ctx, deg, L, U, g, nodes);
            }
        }
        if (rule == LineRule::interface)
        {
            // one curved function at the top level: phi; one root per certified line
            Tag tag = below;
            tag.box[D - 1] = box_id;
            for (std::size_t r = 0; r < nodes.size(); ++r)
            {
                const double t = nodes[r];
                tag.segment[D - 1] = static_cast<int>(r);
                const VecD<D> y = insert<D>(yb, k, t);
                VecD<D> grad;
                for (int j = 0; j < D; ++j)
                    grad[j] = derivative<D>(ctx, p.funcs.front(), y, j);
                const double dk = std::abs(grad[k]);
                if (dk > 0)
                    emit(y, w * scaled_norm(grad) / dk, tag);
            }
            return;
        }
        nodes.push_back(L);
        nodes.push_back(U);
        std::sort(nodes.begin(), nodes.end());
        Tag tag = below;
        tag.box[D - 1] = box_id;
        for (std::size_t s = 0; s + 1 < nodes.size(); ++s)
        {
            const double a = nodes[s], b = nodes[s + 1];
            // leaf cells keep degenerate segments so segment numbers stay aligned
            if (ctx.vis_order == 0 && b - a <= 1e-15)
                continue;
            tag.segment[D - 1] = static_cast<int>(s);
            segment_nodes(ctx, a, b,
                          [&](int j, double t, double wt)
                          {
                              tag.node[D - 1] = j;
                              emit(insert<D>(yb, k, t), w * wt, tag);
                          });
        }
    };
    integrate<D - 1>(ctx, std::move(base), line, 0, LineRule::segments);
}

Vec3 phi_gradient(const algoim::xarray<real, 3>& phi, const Vec3& u)
{
    const auto g = algoim::bernstein::evalBernsteinPolyGradient(phi, uv(u));
    return {g(0), g(1), g(2)};
}

double surface_scale(const Mat3& inv, const Vec3& g)
{
    // the ratio does not depend on the scale of g; normalise it first
    const double m = scaled_norm(g);
    if (!(m > 0))
        return 0.0;
    Vec3 y = {0, 0, 0};
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            y[i] += inv[k][i] * g[k] / m;
    return scaled_norm(y);
}
/// Bernstein form of the level set on the cell's unit box.
struct CellLevelSet
{
    std::vector<real> buffer;
    algoim::xarray<real, 3> phi;

    CellLevelSet(const ClippedBox& cell, const LevelSet& ls)
        : buffer(ls.tape ? 1 : static_cast<std::size_t>((ls.degree + 1) * (ls.degree + 1) * (ls.degree + 1))),
          phi(buffer.data(), algoim::uvector<int, 3>(ls.tape ? 1 : ls.degree + 1))
    {
        if (ls.tape)
        {
            // a tape is evaluated directly; its scale is the largest |phi| at the box corners
            buffer[0] = 0;
            for (int c = 0; c < 8; ++c)
                buffer[0] = std::max(buffer[0], std::abs(tape_value(*ls.tape, physical_point(cell, {double(c & 1),
                                                                         double((c >> 1) & 1), double((c >> 2) & 1)}))));
            return;
        }
        algoim::bernstein::bernsteinInterpolate<3>(
            [&](const algoim::uvector<real, 3>& u) { return ls.value(physical_point(cell, {u(0), u(1), u(2)})); }, phi);
    }

    double scale() const
    {
        double m = 0;
        for (real v : buffer)
            m = std::max(m, std::abs(v));
        return m;
    }
};

/// Context for a cell: phi's Bernstein form, or the tape and the cell's map.
Context cell_context(const ClippedBox& cell, const LevelSet& ls, const CellLevelSet& cls, const CertifyOptions& opt,
                     CertifyStats& stats)
{
    Context ctx;
    ctx.phi = &cls.phi;
    ctx.phi_scale = cls.scale();
    ctx.degree = ls.degree;
    ctx.opt = opt;
    ctx.stats = &stats;
    ctx.tape = ls.tape;
    ctx.x0 = cell.origin;
    ctx.jac = cell.jacobian;
    return ctx;
}

double phi_at(const Context& ctx, const CellLevelSet& cls, const Vec3& u)
{
    return ctx.tape ? tape_value(*ctx.tape, to_x(ctx, u)) : algoim::bernstein::evalBernsteinPoly(cls.phi, uv(u));
}

Vec3 gradient_at(const Context& ctx, const CellLevelSet& cls, const Vec3& u)
{
    return ctx.tape ? tape_gradient_u(ctx, u) : phi_gradient(cls.phi, u);
}

/// The top level: the unit box with the cell's clips and phi itself.
Problem<3> top_problem(const ClippedBox& cell)
{
    Problem<3> top;
    top.lo = {0, 0, 0};
    top.hi = {1, 1, 1};
    Func<3> f;
    for (int i = 0; i < 3; ++i)
        f.A[i][i] = 1.0;
    top.funcs.push_back(f);
    for (const HalfSpace& h : cell.clips)
    {
        Half<3> c;
        c.c = h.c;
        c.d = h.d;
        top.clips.push_back(c);
    }
    return top;
}
} // namespace certify_detail

void certified_bisection(const ClippedBox& cell, const LevelSet& ls, const SelectionTerm& term, int q,
                         const CertifyOptions& opt, Rule& rule, CertifyStats& stats)
{
    using namespace certify_detail;
    const CellLevelSet cls(cell, ls);
    Context ctx = cell_context(cell, ls, cls, opt, stats);
    ctx.q = q;

    const PartKind kind = part_kind(term);
    const double detj = std::abs(jacobian_determinant(cell));
    auto append = [&](const VecD<3>& u, double w)
    {
        const Vec3 xi = reference_point(cell, u);
        rule.points.insert(rule.points.end(), xi.begin(), xi.end());
        rule.weights.push_back(w);
    };
    if (kind == PartKind::interface)
    {
        const Mat3 inv = inverse_jacobian(cell);
        integrate<3>(
            ctx, top_problem(cell),
            [&](const VecD<3>& u, double w, const Tag&)
            { append(u, w * detj * surface_scale(inv, gradient_at(ctx, cls, u))); },
            0, LineRule::interface);
        return;
    }
    integrate<3>(
        ctx, top_problem(cell),
        [&](const VecD<3>& u, double w, const Tag&)
        {
            if (kind != PartKind::whole)
            {
                const double v = phi_at(ctx, cls, u);
                if ((kind == PartKind::negative && !(v < 0)) || (kind == PartKind::positive && !(v > 0)))
                    return;
            }
            append(u, w * detj);
        },
        0, LineRule::segments);
}

void certified_leaves(const ClippedBox& cell, const LevelSet& ls, const SelectionTerm& term, int degree,
                      const CertifyOptions& opt, std::int32_t parent, LeafMesh& mesh, CertifyStats& stats)
{
    using namespace certify_detail;
    if (degree < 1)
        throw std::runtime_error("certified_leaves: degree must be at least 1");
    const CellLevelSet cls(cell, ls);
    Context ctx = cell_context(cell, ls, cls, opt, stats);
    ctx.q = 1;
    ctx.vis_order = degree;

    const PartKind kind = part_kind(term);
    const bool surface = kind == PartKind::interface;
    const int p1 = degree + 1;
    const int n_nodes = surface ? p1 * p1 : p1 * p1 * p1;

    struct Leaf
    {
        std::vector<Vec3> u;
        std::vector<char> set;
        double phi_sum = 0;
    };
    std::map<std::array<int, 6>, Leaf> leaves;
    integrate<3>(
        ctx, top_problem(cell),
        [&](const VecD<3>& u, double, const Tag& tag)
        {
            const std::array<int, 6> key
                = {tag.box[0], tag.box[1], tag.box[2], tag.segment[0], tag.segment[1], tag.segment[2]};
            Leaf& leaf = leaves[key];
            if (leaf.u.empty())
            {
                leaf.u.resize(n_nodes);
                leaf.set.assign(n_nodes, 0);
            }
            const int index = surface ? tag.node[0] + p1 * tag.node[1]
                                      : tag.node[0] + p1 * (tag.node[1] + p1 * tag.node[2]);
            if (index < 0 || index >= n_nodes || leaf.set[index])
                return;
            leaf.u[index] = u;
            leaf.set[index] = 1;
            leaf.phi_sum += phi_at(ctx, cls, u);
        },
        0, surface ? LineRule::interface : LineRule::segments);

    std::vector<std::int32_t> conn(n_nodes);
    for (const auto& [key, leaf] : leaves)
    {
        if (std::count(leaf.set.begin(), leaf.set.end(), char(1)) != n_nodes)
        {
            ++stats.incomplete_leaves; // node counts differed across the leaf
            static const bool debug = std::getenv("LEAF_DEBUG") != nullptr;
            if (debug)
            {
                std::fprintf(stderr, "incomplete leaf: boxes %d %d %d segments %d %d %d, missing nodes:", key[0], key[1],
                             key[2], key[3], key[4], key[5]);
                for (int index = 0; index < n_nodes; ++index)
                    if (!leaf.set[index])
                        std::fprintf(stderr, " (%d,%d,%d)", index % p1, (index / p1) % p1, index / (p1 * p1));
                std::fprintf(stderr, "\n");
            }
            continue;
        }
        if (!surface && kind != PartKind::whole)
        {
            // a volume leaf lies on one side of the level set; nodes on the interface are ~0
            if ((kind == PartKind::negative && !(leaf.phi_sum < 0)) || (kind == PartKind::positive && !(leaf.phi_sum > 0)))
                continue;
        }
        const std::int32_t first = mesh.n_points();
        std::vector<Vec3> x(n_nodes);
        for (int index = 0; index < n_nodes; ++index)
        {
            x[index] = physical_point(cell, leaf.u[index]);
            mesh.points.insert(mesh.points.end(), x[index].begin(), x[index].end());
        }
        // Orientation: hexahedra with a positive Jacobian, interface quadrilaterals with
        // their normal along grad phi. Summing over all corners keeps the sign reliable
        // for collapsed leaves.
        auto node = [&](int i, int j, int kk) -> const Vec3& { return x[surface ? i + p1 * j : i + p1 * (j + p1 * kk)]; };
        auto sub = [](const Vec3& a, const Vec3& b) { return Vec3{a[0] - b[0], a[1] - b[1], a[2] - b[2]}; };
        auto cross = [](const Vec3& a, const Vec3& b)
        { return Vec3{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; };
        double orientation = 0;
        if (surface)
        {
            Vec3 normal = {0, 0, 0};
            for (int c = 0; c < 4; ++c)
            {
                const int i0 = (c & 1) * degree, j0 = ((c >> 1) & 1) * degree;
                const double s = (i0 ? -1.0 : 1.0) * (j0 ? -1.0 : 1.0);
                const Vec3 n = cross(sub(node(degree - i0, j0, 0), node(i0, j0, 0)), sub(node(i0, degree - j0, 0), node(i0, j0, 0)));
                for (int d = 0; d < 3; ++d)
                    normal[d] += s * n[d];
            }
            // physical gradient of phi at the leaf centre: J^{-T} grad_u phi
            const Mat3 inv = inverse_jacobian(cell);
            const Vec3 gu = gradient_at(ctx, cls, leaf.u[n_nodes / 2]);
            for (int d = 0; d < 3; ++d)
            {
                double gx = 0;
                for (int e = 0; e < 3; ++e)
                    gx += inv[e][d] * gu[e];
                orientation += normal[d] * gx;
            }
        }
        else
            for (int c = 0; c < 8; ++c)
            {
                const int i0 = (c & 1) * degree, j0 = ((c >> 1) & 1) * degree, k0 = ((c >> 2) & 1) * degree;
                const double s = (i0 ? -1.0 : 1.0) * (j0 ? -1.0 : 1.0) * (k0 ? -1.0 : 1.0);
                const Vec3 o = node(i0, j0, k0);
                const Vec3 n = cross(sub(node(degree - i0, j0, k0), o), sub(node(i0, degree - j0, k0), o));
                const Vec3 e = sub(node(i0, j0, degree - k0), o);
                orientation += s * (n[0] * e[0] + n[1] * e[1] + n[2] * e[2]);
            }
        const bool flip = orientation < 0;
        if (surface)
            for (int j = 0; j < p1; ++j)
                for (int i = 0; i < p1; ++i)
                    conn[vtk_lagrange_quad_index(flip ? degree - i : i, j, degree)] = first + i + p1 * j;
        else
            for (int kk = 0; kk < p1; ++kk)
                for (int j = 0; j < p1; ++j)
                    for (int i = 0; i < p1; ++i)
                        conn[vtk_lagrange_hex_index(flip ? degree - i : i, j, kk, degree)] = first + i + p1 * (j + p1 * kk);
        mesh.connectivity.insert(mesh.connectivity.end(), conn.begin(), conn.end());
        mesh.offsets.push_back(static_cast<std::int32_t>(mesh.connectivity.size()));
        mesh.types.push_back(surface ? vtk_lagrange_quadrilateral : vtk_lagrange_hexahedron);
        mesh.parent.push_back(parent);
        mesh.degree.push_back(degree);
    }
}

} // namespace cutcells::proto
