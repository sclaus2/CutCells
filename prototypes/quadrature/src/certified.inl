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
#include <functional>
#include <limits>
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

struct Context
{
    const algoim::xarray<real, 3>* phi = nullptr; ///< Bernstein form of phi on the cell's unit box
    int degree = 2;                               ///< degree of phi in each variable
    int q = 3;
    CertifyOptions opt;
    CertifyStats* stats = nullptr;
};

/// psi(y) = phi(A y + b) (curved), or a . y + c (linear), on level coordinates y.
template <int D>
struct Func
{
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
};

/// y_k = alpha + beta . y' on the base coordinates y'
template <int D>
struct Bound
{
    double alpha = 0;
    VecD<D - 1> beta{};
};

template <int D>
using Emit = std::function<void(const VecD<D>&, double)>;

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
    return algoim::bernstein::evalBernsteinPoly(*ctx.phi, uv(to_u<D>(f, y)));
}

template <int D>
double derivative(const Context& ctx, const Func<D>& f, const VecD<D>& y, int k)
{
    if (f.linear)
        return f.a[k];
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
    double norm = 0;
    for (int k = 0; k < D; ++k)
        norm += upper[k] * upper[k];
    norm = std::sqrt(norm);
    VecD<D> ratio{};
    for (int k = 0; k < D; ++k)
        ratio[k] = norm > 0 ? lower[k] / norm : 0.0;
    return ratio;
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
        line_roots(ctx, deg, L, U, [&](double t) { return value<1>(ctx, f, {t}); }, nodes);
    }
    std::sort(nodes.begin(), nodes.end());
    for (std::size_t s = 0; s + 1 < nodes.size(); ++s)
    {
        const double a = nodes[s], b = nodes[s + 1];
        if (b - a <= 1e-15)
            continue;
        for (int j = 0; j < ctx.q; ++j)
            emit({a + (b - a) * algoim::GaussQuad::x(ctx.q, j)}, (b - a) * algoim::GaussQuad::w(ctx.q, j));
    }
}

template <int D>
void integrate(const Context& ctx, Problem<D> p, const Emit<D>& emit, int depth, LineRule rule)
{
    if (!tighten<D>(p.lo, p.hi, p.clips))
        return;

    // keep the functions that may vanish on the box, with their direction margins
    std::vector<Func<D>> funcs;
    std::vector<VecD<D>> ratios;
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
        bernstein_form<D>(ctx, f, p.lo, p.hi, coeffs, ext);
        if (!may_vanish(coeffs))
            continue;
        funcs.push_back(f);
        ratios.push_back(margins<D>(coeffs, ext, p.lo, p.hi));
    }
    if (rule == LineRule::interface && funcs.empty())
        return; // the level set does not cross this box
    p.funcs = std::move(funcs);

    // height direction: the best certified margin over all curved functions
    int k = 0;
    double best = -1;
    for (int kk = 0; kk < D; ++kk)
    {
        double score = 1.0;
        for (const VecD<D>& r : ratios)
            score = std::min(score, r[kk]);
        if (score > best || (score == best && p.hi[kk] - p.lo[kk] > p.hi[k] - p.lo[k]))
        {
            best = score;
            k = kk;
        }
    }
    const bool certified = ratios.empty() || (best > 0 && best >= ctx.opt.margin);
    if (!certified)
    {
        if (depth < ctx.opt.max_depth)
        {
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
    const Emit<D - 1> line = [&, k, certified](const VecD<D - 1>& yb, double w)
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
            for (double t : nodes)
            {
                const VecD<D> y = insert<D>(yb, k, t);
                double gn = 0;
                for (int j = 0; j < D; ++j)
                {
                    const double dj = derivative<D>(ctx, p.funcs.front(), y, j);
                    gn += dj * dj;
                }
                const double dk = std::abs(derivative<D>(ctx, p.funcs.front(), y, k));
                if (dk > 0)
                    emit(y, w * std::sqrt(gn) / dk);
            }
            return;
        }
        nodes.push_back(L);
        nodes.push_back(U);
        std::sort(nodes.begin(), nodes.end());
        for (std::size_t s = 0; s + 1 < nodes.size(); ++s)
        {
            const double a = nodes[s], b = nodes[s + 1];
            if (b - a <= 1e-15)
                continue;
            for (int j = 0; j < ctx.q; ++j)
                emit(insert<D>(yb, k, a + (b - a) * algoim::GaussQuad::x(ctx.q, j)),
                     w * (b - a) * algoim::GaussQuad::w(ctx.q, j));
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
    Vec3 y = {0, 0, 0};
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            y[i] += inv[k][i] * g[k];
    return std::sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]) / std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
}
} // namespace certify_detail

void certified_bisection(const ClippedBox& cell, const LevelSet& ls, const SelectionTerm& term, int q,
                         const CertifyOptions& opt, Rule& rule, CertifyStats& stats)
{
    using namespace certify_detail;
    const int extent = ls.degree + 1;
    std::vector<real> buffer(extent * extent * extent);
    algoim::xarray<real, 3> phi(buffer.data(), algoim::uvector<int, 3>(extent));
    algoim::bernstein::bernsteinInterpolate<3>(
        [&](const algoim::uvector<real, 3>& u) { return ls.value(physical_point(cell, {u(0), u(1), u(2)})); }, phi);

    Context ctx;
    ctx.phi = &phi;
    ctx.degree = ls.degree;
    ctx.q = q;
    ctx.opt = opt;
    ctx.stats = &stats;

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
            ctx, top, [&](const VecD<3>& u, double w) { append(u, w * detj * surface_scale(inv, phi_gradient(phi, u))); },
            0, LineRule::interface);
        return;
    }
    integrate<3>(
        ctx, top,
        [&](const VecD<3>& u, double w)
        {
            if (kind != PartKind::whole)
            {
                const double v = algoim::bernstein::evalBernsteinPoly(phi, uv(u));
                if ((kind == PartKind::negative && !(v < 0)) || (kind == PartKind::positive && !(v > 0)))
                    return;
            }
            append(u, w * detj);
        },
        0, LineRule::segments);
}

} // namespace cutcells::proto
