// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "engine.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cutcells::quadrays
{

Part part_of(const SelectionTerm& term, int level_set)
{
    if (level_set < 0 || level_set >= 64)
        throw std::invalid_argument("quadrays: level-set index out of range");
    const std::uint64_t bit = std::uint64_t(1) << level_set;
    if ((term.zero_required | term.negative_required | term.positive_required) & ~bit)
    {
        throw std::invalid_argument(
            "quadrays: the selection term constrains another level set; "
            "one level set per term is supported");
    }
    const bool zero = term.zero_required & bit, negative = term.negative_required & bit,
               positive = term.positive_required & bit;
    if (int(zero) + int(negative) + int(positive) > 1)
        throw std::invalid_argument("quadrays: the selection term asks for contradictory signs");
    if (zero)
        return Part::interface;
    if (negative)
        return Part::negative;
    if (positive)
        return Part::positive;
    return Part::whole;
}

namespace
{

template <typename T, int D>
using VecD = std::array<T, D>;

// Tolerances of the prototype in double precision, scaled for float.
template <std::floating_point T>
constexpr T tiny = scaled_tolerance<T>(1e-14);
template <std::floating_point T>
constexpr T vertex_det_tol = scaled_tolerance<T>(1e-13);
template <std::floating_point T>
constexpr T feasible_tol = scaled_tolerance<T>(1e-11);
template <std::floating_point T>
constexpr T extent_tol = scaled_tolerance<T>(1e-13);
template <std::floating_point T>
constexpr T compare_tol = scaled_tolerance<T>(1e-12);
template <std::floating_point T>
constexpr T segment_tol = scaled_tolerance<T>(1e-15);
template <std::floating_point T>
constexpr T zero_function_tol = scaled_tolerance<T>(1e-12);

template <std::floating_point T>
constexpr T infinity = std::numeric_limits<T>::infinity();

/// Entries of a level map below tiny count as zero, as for the degrees.
template <std::floating_point T>
T snap(T a)
{
    return std::abs(a) > tiny<T> ? a : T(0);
}

// ============================================================================
// Problems at one level
// ============================================================================

/// psi(y) = phi(A y + b) (curved), or a . y + c (linear), on level coordinates y.
template <typename T, int D>
struct Func
{
    std::uint32_t origin = 1; ///< diagnostics: chain of 4-bit codes, see origin_name
    bool linear = false;
    std::array<VecD<T, D>, 3> A{};
    Vec3<T> b{};
    VecD<T, D> a{};
    T c = 0;
};

/// c . y <= d
template <typename T, int D>
struct Half
{
    VecD<T, D> c{};
    T d = 0;
};

template <typename T, int D>
struct Problem
{
    VecD<T, D> lo{}, hi{};
    std::vector<Func<T, D>> funcs;
    std::vector<Half<T, D>> clips;
};

/// y_k = alpha + beta . y' on the base coordinates y'
template <typename T, int D>
struct Bound
{
    T alpha = 0;
    VecD<T, D - 1> beta{};
    int kind = 7; ///< diagnostics: 2/3 box face below/above, 4/5 clip below/above, 6 linear root
};

/// Origin codes: 1 phi; restricted to 2 lower box face, 3 upper box face, 4 lower
/// clip, 5 upper clip, 6 linear root, 7 unchanged; new linear functions: 8 switch
/// of lower bounds, 9 switch of upper bounds, 10 difference of linear roots.
std::string origin_name(std::uint32_t origin)
{
    static const char* names[] = {"?",       "phi",      "box-lo",    "box-hi",   "clip-lo", "clip-up",
                                  "lin-root", "same",    "switch-lo", "switch-up", "root-diff"};
    std::vector<std::string> parts;
    for (; origin != 0; origin >>= 4)
        parts.push_back(names[std::min<std::uint32_t>(origin & 15u, 10u)]);
    std::string s;
    for (auto it = parts.rbegin(); it != parts.rend(); ++it)
        s += (s.empty() ? "" : "|") + *it;
    return s;
}

// ============================================================================
// Context: the cell, the options and scratch buffers per level
// ============================================================================

template <std::floating_point T>
struct Context
{
    Source<T> phi;  ///< the level set on the cell's unit box
    T phi_scale = 0; ///< size of phi (reference_magnitude)
    Options opt;
    T margin = 0;
    Stats* stats = nullptr;
    int vis_order = 0; ///< > 0: leaf nodes (vis_order + 1 per segment) instead of Gauss points
    Part part = Part::whole;
    T detj = 0;       ///< |det| of the box-to-physical map
    Mat3<T> inv = {}; ///< its inverse (interface weights)
    std::array<int, 3> box_ids{};
    int bisections = 0; ///< bisections so far in this cell

    // Gauss-Legendre rule on [0, 1]
    int q = 0;
    std::vector<T> gauss_x, gauss_w;

    // scratch, indexed by level (number of free coordinates)
    std::array<BoxBernstein<T>, 4> form, deriv, sub, line;
    std::array<std::vector<T>, 4> form_work, line_work, root_work, nodes;
};

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

// ============================================================================
// Evaluation
// ============================================================================

template <typename T, int D>
Vec3<T> to_u(const Func<T, D>& f, const VecD<T, D>& y)
{
    Vec3<T> u = f.b;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < D; ++j)
            u[i] += f.A[i][j] * y[j];
    return u;
}

template <std::floating_point T, int D>
T value(const Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& y)
{
    if (f.linear)
    {
        T s = f.c;
        for (int j = 0; j < D; ++j)
            s += f.a[j] * y[j];
        return s;
    }
    const Vec3<T> u = to_u<T, D>(f, y);
    return evaluate(ctx.phi, std::span<const T>(u));
}

template <std::floating_point T, int D>
T derivative_along(const Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& y, int k)
{
    if (f.linear)
        return f.a[k];
    const Vec3<T> u = to_u<T, D>(f, y);
    Vec3<T> g;
    gradient(ctx.phi, std::span<const T>(u), std::span<T>(g));
    T s = T(0);
    for (int i = 0; i < 3; ++i)
        s += g[i] * f.A[i][k];
    return s;
}

template <typename T, int D>
VecD<T, D> insert(const VecD<T, D - 1>& yb, int k, T t)
{
    VecD<T, D> y{};
    for (int j = 0, jb = 0; j < D; ++j)
        y[j] = j == k ? t : yb[jb++];
    return y;
}

template <typename T, int D>
T bound_value(const Bound<T, D>& b, const VecD<T, D - 1>& yb)
{
    T s = b.alpha;
    for (int j = 0; j < D - 1; ++j)
        s += b.beta[j] * yb[j];
    return s;
}

/// f restricted to y_k = alpha + beta . y'
template <typename T, int D>
Func<T, D - 1> restrict_to(const Func<T, D>& f, int k, const Bound<T, D>& b)
{
    Func<T, D - 1> r;
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
template <typename T, int D>
Func<T, D - 1> difference(const Bound<T, D>& b1, const Bound<T, D>& b2)
{
    Func<T, D - 1> r;
    r.origin = b1.kind == 6 ? 10u : (b1.kind == 2 || b1.kind == 4 ? 8u : 9u);
    r.linear = true;
    r.c = b1.alpha - b2.alpha;
    for (int j = 0; j < D - 1; ++j)
        r.a[j] = b1.beta[j] - b2.beta[j];
    return r;
}

// ============================================================================
// Bernstein forms of curved functions
// ============================================================================

/// Bernstein form of a curved function on the box [lo, hi] of its level: exact
/// restriction of phi to the affine image of the box.
template <std::floating_point T, int D>
void bernstein_form(Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi,
                    BoxBernstein<T>& out)
{
    std::array<T, 3> origin;
    std::array<T, 3 * D> matrix;
    for (int i = 0; i < 3; ++i)
    {
        T o = f.b[i];
        for (int j = 0; j < D; ++j)
        {
            const T a = snap(f.A[i][j]);
            o += a * lo[j];
            matrix[i * D + j] = a * (hi[j] - lo[j]);
        }
        origin[i] = o;
    }
    restrict_affine(*ctx.phi.bernstein, std::span<const T>(origin), std::span<const T>(matrix), D, out,
                    ctx.form_work[D]);
}

/// Bernstein form on [L, U] of t -> psi(y0 + (t - L) e_k), y0 on the line at t = L.
template <std::floating_point T, int D>
void line_form(Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& y0, int k, T L, T U,
               BoxBernstein<T>& out)
{
    std::array<T, 3> origin, column;
    for (int i = 0; i < 3; ++i)
    {
        T o = f.b[i];
        for (int j = 0; j < D; ++j)
            o += snap(f.A[i][j]) * y0[j];
        origin[i] = o;
        column[i] = snap(f.A[i][k]) * (U - L);
    }
    restrict_affine(*ctx.phi.bernstein, std::span<const T>(origin), std::span<const T>(column), 1, out,
                    ctx.line_work[D]);
}

/// The line y = y0 + t e_k of a level as u = u0 + t dir in box coordinates.
template <typename T, int D>
void box_line(const Func<T, D>& f, const VecD<T, D>& y0, int k, Vec3<T>& u0, Vec3<T>& dir)
{
    u0 = to_u<T, D>(f, y0);
    for (int i = 0; i < 3; ++i)
        dir[i] = f.A[i][k];
}

/// Bounds of a curved function on [lo, hi] for an analytic level set, from
/// Taylor models of phi on the affine image of the box (the map to physical
/// space is affine, so the models represent it exactly). Sets the sign of the
/// function if it is certain, the direction margins as margins() does from
/// Bernstein coefficients, and the largest |psi|. Without derivative bounds
/// no direction is certified; without any bound nothing is certain.
template <std::floating_point T, int D>
void analytic_bounds(const Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi,
                     int& sign, VecD<T, D>& ratio, T& magnitude)
{
    std::array<T, 3> origin;
    std::array<T, 3 * D> matrix;
    VecD<T, D> lengths;
    for (int j = 0; j < D; ++j)
        lengths[j] = hi[j] - lo[j];
    for (int i = 0; i < 3; ++i)
    {
        T o = f.b[i];
        for (int j = 0; j < D; ++j)
        {
            o += f.A[i][j] * lo[j];
            matrix[i * D + j] = f.A[i][j] * lengths[j];
        }
        origin[i] = o;
    }
    AffineBounds<T> b;
    if (!affine_bounds(ctx.phi, std::span<const T>(origin), std::span<const T>(matrix), D, b))
    {
        sign = 0;
        ratio.fill(T(0));
        magnitude = infinity<T>;
        return;
    }
    sign = b.sign;
    magnitude = b.magnitude;
    if (!b.has_derivatives)
    {
        ratio.fill(T(0));
        return;
    }
    // derivatives with respect to y_j = lo_j + lengths_j s_j
    VecD<T, D> lower, upper;
    for (int j = 0; j < D; ++j)
    {
        lower[j] = b.lower[j] / lengths[j];
        upper[j] = b.upper[j] / lengths[j];
    }
    const T norm = scaled_norm(std::span<const T>(upper));
    for (int j = 0; j < D; ++j)
        ratio[j] = norm > T(0) ? lower[j] / norm : T(0);
}

/// Margins from M^D sub-cells. Only sub-cells that may meet the clipped region and
/// on which psi may vanish count. If d_k psi has one strict sign on all of them,
/// psi has at most one root on every height line in the region, and ratio[k] is the
/// smallest local min|d_k psi| / max|grad psi|. Returns false if no sub-cell
/// counts: the function does not vanish in the cell and can be dropped.
template <std::floating_point T, int D>
bool local_margins(Context<T>& ctx, const BoxBernstein<T>& p, const VecD<T, D>& lo, const VecD<T, D>& hi,
                   const std::vector<Half<T, D>>& clips, int M, VecD<T, D>& ratio)
{
    VecD<T, D> local, scale;
    std::array<int, D> sign{};
    local.fill(infinity<T>);
    for (int j = 0; j < D; ++j)
        scale[j] = M / (hi[j] - lo[j]);
    bool relevant = false;
    int cells = 1;
    for (int j = 0; j < D; ++j)
        cells *= M;
    BoxBernstein<T>& sub = ctx.sub[D];
    BoxBernstein<T>& dk = ctx.deriv[D];
    for (int n = 0; n < cells; ++n)
    {
        VecD<T, D> a, b, ylo, yhi;
        for (int j = 0, m = n; j < D; ++j, m /= M)
        {
            a[j] = T(m % M) / M;
            b[j] = T(m % M + 1) / M;
            ylo[j] = lo[j] + (hi[j] - lo[j]) * a[j];
            yhi[j] = lo[j] + (hi[j] - lo[j]) * b[j];
        }
        bool meets = true;
        for (const Half<T, D>& h : clips)
        {
            T vmin = T(0);
            for (int j = 0; j < D; ++j)
                vmin += h.c[j] * (h.c[j] >= T(0) ? ylo[j] : yhi[j]);
            meets &= vmin <= h.d + compare_tol<T>;
        }
        if (!meets)
            continue;
        subdivide(p, std::span<const T>(a), std::span<const T>(b), sub, ctx.form_work[D]);
        if (!may_vanish(std::span<const T>(sub.coeffs)))
            continue;
        relevant = true;
        VecD<T, D> dmin{}, dmax{};
        std::array<int, D> s{};
        for (int k = 0; k < D; ++k)
        {
            if (sub.degree[k] == 0)
                continue;
            derivative(sub, k, dk);
            bool pos = true, neg = true;
            T amin = infinity<T>, amax = T(0);
            for (const T v : dk.coeffs)
            {
                pos &= v > T(0);
                neg &= v < T(0);
                amin = std::min(amin, std::abs(v));
                amax = std::max(amax, std::abs(v));
            }
            s[k] = pos ? 1 : (neg ? -1 : 0);
            dmin[k] = (pos || neg) ? amin * scale[k] : T(0);
            dmax[k] = amax * scale[k];
        }
        const T norm = scaled_norm(std::span<const T>(dmax));
        for (int k = 0; k < D; ++k)
        {
            // one strict sign shared by all counted sub-cells
            if (s[k] == 0 || (sign[k] != 0 && sign[k] != s[k]))
                sign[k] = 2;
            else if (sign[k] == 0)
                sign[k] = s[k];
            local[k] = std::min(local[k], norm > T(0) ? dmin[k] / norm : T(0));
        }
    }
    for (int k = 0; k < D; ++k)
        ratio[k] = (relevant && (sign[k] == 1 || sign[k] == -1)) ? local[k] : T(0);
    return relevant;
}

// ============================================================================
// Clipped boxes at one level
// ============================================================================

template <typename T, int D>
T determinant(const std::array<VecD<T, D>, D>& m)
{
    if constexpr (D == 1)
        return m[0][0];
    else if constexpr (D == 2)
        return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    else
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
               - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
               + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

/// Vertices of {y in [lo, hi] : clips} by enumeration of plane D-tuples.
template <std::floating_point T, int D>
void polytope_vertices(const VecD<T, D>& lo, const VecD<T, D>& hi, const std::vector<Half<T, D>>& clips,
                       std::vector<VecD<T, D>>& vertices)
{
    vertices.clear();
    std::vector<Half<T, D>> planes;
    planes.reserve(2 * D + clips.size());
    for (int j = 0; j < D; ++j)
    {
        Half<T, D> l, u;
        l.c[j] = T(-1);
        l.d = -lo[j];
        u.c[j] = T(1);
        u.d = hi[j];
        planes.push_back(l);
        planes.push_back(u);
    }
    planes.insert(planes.end(), clips.begin(), clips.end());
    const int np = static_cast<int>(planes.size());

    auto vertex = [&](const std::array<int, D>& idx)
    {
        std::array<VecD<T, D>, D> m;
        for (int r = 0; r < D; ++r)
            m[r] = planes[idx[r]].c;
        const T det = determinant<T, D>(m);
        if (std::abs(det) < vertex_det_tol<T>)
            return;
        VecD<T, D> x;
        for (int col = 0; col < D; ++col)
        {
            auto mc = m;
            for (int r = 0; r < D; ++r)
                mc[r][col] = planes[idx[r]].d;
            x[col] = determinant<T, D>(mc) / det;
        }
        for (const Half<T, D>& h : planes)
        {
            T s = T(0);
            for (int j = 0; j < D; ++j)
                s += h.c[j] * x[j];
            if (s > h.d + feasible_tol<T>)
                return;
        }
        vertices.push_back(x);
    };
    if constexpr (D == 1)
    {
        for (int a = 0; a < np; ++a)
            vertex({a});
    }
    else if constexpr (D == 2)
    {
        for (int a = 0; a < np; ++a)
            for (int b = a + 1; b < np; ++b)
                vertex({a, b});
    }
    else
    {
        for (int a = 0; a < np; ++a)
            for (int b = a + 1; b < np; ++b)
                for (int c = b + 1; c < np; ++c)
                    vertex({a, b, c});
    }
}

/// Shrink [lo, hi] to the bounding box of {y in [lo, hi] : clips}; false if empty.
template <std::floating_point T, int D>
bool tighten(VecD<T, D>& lo, VecD<T, D>& hi, const std::vector<Half<T, D>>& clips)
{
    if (clips.empty())
        return true;
    std::vector<VecD<T, D>> vertices;
    polytope_vertices<T, D>(lo, hi, clips, vertices);
    if (vertices.empty())
        return false;
    VecD<T, D> nlo, nhi;
    nlo.fill(infinity<T>);
    nhi.fill(-infinity<T>);
    for (const VecD<T, D>& x : vertices)
        for (int j = 0; j < D; ++j)
        {
            nlo[j] = std::min(nlo[j], x[j]);
            nhi[j] = std::max(nhi[j], x[j]);
        }
    for (int j = 0; j < D; ++j)
    {
        lo[j] = std::max(lo[j], nlo[j]);
        hi[j] = std::min(hi[j], nhi[j]);
        if (!(hi[j] - lo[j] > extent_tol<T>))
            return false;
    }
    return true;
}

template <std::floating_point T, int D>
bool linear_may_vanish(const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi)
{
    T vmin = f.c, vmax = f.c;
    for (int j = 0; j < D; ++j)
    {
        vmin += std::min(f.a[j] * lo[j], f.a[j] * hi[j]);
        vmax += std::max(f.a[j] * lo[j], f.a[j] * hi[j]);
    }
    return vmin <= T(0) && vmax >= T(0);
}

// ============================================================================
// Analysis of one box: which functions may vanish, and in which direction
// ============================================================================

template <typename T, int D>
struct Analysis
{
    std::vector<Func<T, D>> funcs;
    std::vector<VecD<T, D>> ratios;
    std::vector<int> curved; ///< index in funcs of each entry of ratios
    int k = 0;
    T best = -1;
    bool certified = false;
    int dropped_sign = 0; ///< sign of a curved function dropped for having one sign on the box
};

template <std::floating_point T, int D>
Analysis<T, D> analyse(Context<T>& ctx, const Problem<T, D>& p)
{
    Analysis<T, D> an;
    BoxBernstein<T>& form = ctx.form[D];
    VecD<T, D> lengths;
    for (int j = 0; j < D; ++j)
        lengths[j] = p.hi[j] - p.lo[j];
    for (const Func<T, D>& f : p.funcs)
    {
        if (f.linear)
        {
            if (linear_may_vanish<T, D>(f, p.lo, p.hi))
                an.funcs.push_back(f);
            continue;
        }
        if (ctx.phi.is_analytic())
        {
            int sign = 0;
            VecD<T, D> ratio{};
            T magnitude = T(0);
            analytic_bounds<T, D>(ctx, f, p.lo, p.hi, sign, ratio, magnitude);
            // phi = 0 on a face, say: no root inside the region
            if (magnitude <= zero_function_tol<T> * ctx.phi_scale)
                continue;
            if (sign != 0)
            {
                an.dropped_sign = sign;
                continue;
            }
            an.curved.push_back(static_cast<int>(an.funcs.size()));
            an.funcs.push_back(f);
            an.ratios.push_back(ratio);
            continue;
        }
        bernstein_form<T, D>(ctx, f, p.lo, p.hi, form);
        const std::span<const T> coeffs(form.coeffs);
        // phi vanishes identically here (phi = 0 on a face, say): only rounding
        // noise is left, which would block certification; it has no root inside
        if (max_abs(coeffs) <= zero_function_tol<T> * ctx.phi_scale)
            continue;
        VecD<T, D> ratio{};
        if (ctx.opt.mask_subdivisions > 1)
        {
            if (!local_margins<T, D>(ctx, form, p.lo, p.hi, p.clips, ctx.opt.mask_subdivisions, ratio))
                continue; // no zero of this function inside the cell
        }
        else
        {
            if (!may_vanish(coeffs))
            {
                // one sign on the box, zeros at most where it touches zero
                const T largest = *std::max_element(coeffs.begin(), coeffs.end(), [](T a, T b)
                                                    { return std::abs(a) < std::abs(b); });
                an.dropped_sign = largest > T(0) ? 1 : -1;
                continue;
            }
            margins(form, std::span<const T>(lengths), std::span<T>(ratio), ctx.deriv[D]);
        }
        an.curved.push_back(static_cast<int>(an.funcs.size()));
        an.funcs.push_back(f);
        an.ratios.push_back(ratio);
    }
    // height direction: the best certified margin over all curved functions
    for (int kk = 0; kk < D; ++kk)
    {
        T score = T(1);
        for (const VecD<T, D>& r : an.ratios)
            score = std::min(score, r[kk]);
        if (score > an.best || (score == an.best && p.hi[kk] - p.lo[kk] > p.hi[an.k] - p.lo[an.k]))
        {
            an.best = score;
            an.k = kk;
        }
    }
    an.certified = an.ratios.empty() || (an.best > T(0) && an.best >= ctx.margin);
    return an;
}

/// Why did certification fail? Sample the clipped region on a 9^D grid and look
/// at the function with the smallest margin for the chosen axis near its zero set.
template <std::floating_point T, int D>
void diagnose_failure(Context<T>& ctx, const Problem<T, D>& p, const Analysis<T, D>& an, int depth)
{
    const std::vector<VecD<T, D>>& ratios = an.ratios;
    const int k = an.k;
    std::size_t blocker = 0;
    for (std::size_t i = 1; i < ratios.size(); ++i)
        if (ratios[i][k] < ratios[blocker][k])
            blocker = i;
    const Func<T, D>& f = p.funcs[an.curved[blocker]];
    const int G = 9;
    std::vector<VecD<T, D>> points;
    int total = 1;
    for (int j = 0; j < D; ++j)
        total *= G;
    for (int n = 0; n < total; ++n)
    {
        int m = n;
        VecD<T, D> y;
        for (int j = 0; j < D; ++j)
        {
            y[j] = p.lo[j] + (p.hi[j] - p.lo[j]) * (T(m % G) + T(0.5)) / G;
            m /= G;
        }
        bool inside = true;
        for (const Half<T, D>& h : p.clips)
        {
            T sum = T(0);
            for (int j = 0; j < D; ++j)
                sum += h.c[j] * y[j];
            inside &= sum <= h.d;
        }
        if (inside)
            points.push_back(y);
    }
    std::vector<T> values(points.size());
    T vmax = T(0);
    bool pos = false, neg = false;
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        values[i] = value<T, D>(ctx, f, points[i]);
        vmax = std::max(vmax, std::abs(values[i]));
        pos |= values[i] > T(0);
        neg |= values[i] < T(0);
    }
    std::string category;
    if (!(pos && neg))
        category = "zero set outside the cell";
    else
    {
        VecD<T, D> rmin;
        rmin.fill(infinity<T>);
        for (std::size_t i = 0; i < points.size(); ++i)
        {
            if (std::abs(values[i]) > T(0.2) * vmax)
                continue;
            VecD<T, D> g;
            for (int j = 0; j < D; ++j)
                g[j] = derivative_along<T, D>(ctx, f, points[i], j);
            const T norm = scaled_norm(std::span<const T>(g));
            for (int j = 0; j < D && norm > T(0); ++j)
                rmin[j] = std::min(rmin[j], std::abs(g[j]) / norm);
        }
        T other = T(0);
        for (int j = 0; j < D; ++j)
            other = std::max(other, rmin[j]);
        if (rmin[k] >= ctx.margin)
            category = "margin holds near the zero set in the cell";
        else if (other >= ctx.margin)
            category = "another axis has the margin there";
        else
            category = "margin fails near the zero set in the cell";
    }
    std::map<std::string, std::int64_t>& causes = ctx.stats->causes;
    const std::string level = "level " + std::to_string(D);
    ++causes[level + " | " + origin_name(f.origin) + " | "
             + (an.best > T(0) ? "bound below margin" : "no certain sign") + " | " + category];
    ++causes[level + " | depth " + std::to_string(depth)];

    // From the certified bounds alone: does one function fail on every axis, or
    // does each function have an axis but no axis suits all (a conflict)?
    std::string axes;
    for (std::size_t i = 0; i < ratios.size() && axes.empty(); ++i)
    {
        T top = T(0);
        for (int j = 0; j < D; ++j)
            top = std::max(top, ratios[i][j]);
        if (top < ctx.margin)
            axes = "one function fails on every axis: " + origin_name(p.funcs[an.curved[i]].origin);
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
            names.push_back(origin_name(p.funcs[an.curved[b]].origin));
        }
        std::sort(names.begin(), names.end());
        names.erase(std::unique(names.begin(), names.end()), names.end());
        axes = "conflict:";
        for (const std::string& n : names)
            axes += " " + n;
    }
    // is the box cut by the boundary of the clipped region?
    bool cut = false;
    for (const Half<T, D>& h : p.clips)
    {
        T vmin = T(0), vmax_c = T(0);
        for (int j = 0; j < D; ++j)
        {
            vmin += h.c[j] * (h.c[j] >= T(0) ? p.lo[j] : p.hi[j]);
            vmax_c += h.c[j] * (h.c[j] >= T(0) ? p.hi[j] : p.lo[j]);
        }
        cut |= vmin < h.d - compare_tol<T> && vmax_c > h.d + compare_tol<T>;
    }
    ++causes["axes | " + level + (depth >= ctx.opt.max_depth ? " | at depth limit" : " | before limit") + " | "
             + axes + (cut ? " | box cut by a face" : " | box inside")];
}

/// A level-2 problem in the diagonal frame y = T z, T = [[1/2, -1/2], [1/2, 1/2]]:
/// z_0 runs along (1, 1) and z_1 along (-1, 1). The box becomes four clips.
template <typename T>
Problem<T, 2> rotate_diagonal(const Problem<T, 2>& p)
{
    auto transpose_times = [](const VecD<T, 2>& c)
    { return VecD<T, 2>{T(0.5) * (c[0] + c[1]), T(0.5) * (c[1] - c[0])}; };
    Problem<T, 2> r;
    for (const Half<T, 2>& h : p.clips)
        r.clips.push_back({transpose_times(h.c), h.d});
    for (int j = 0; j < 2; ++j)
    {
        VecD<T, 2> e{};
        e[j] = T(1);
        r.clips.push_back({transpose_times(e), p.hi[j]});
        e[j] = T(-1);
        r.clips.push_back({transpose_times(e), -p.lo[j]});
    }
    // z_0 = y_0 + y_1, z_1 = y_1 - y_0
    r.lo = {p.lo[0] + p.lo[1], p.lo[1] - p.hi[0]};
    r.hi = {p.hi[0] + p.hi[1], p.hi[1] - p.lo[0]};
    for (const Func<T, 2>& f : p.funcs)
    {
        Func<T, 2> g = f;
        if (f.linear)
            g.a = transpose_times(f.a);
        else
            for (int i = 0; i < 3; ++i)
                g.A[i] = transpose_times(f.A[i]);
        r.funcs.push_back(g);
    }
    return r;
}

// ============================================================================
// Integration by dimension reduction
// ============================================================================

/// Nodes on a segment [a, b]: Gauss-Legendre points and weights for quadrature,
/// or vis_order + 1 equispaced points pulled 1e-5 inside the ends for leaf cells
/// (closer to a breakpoint, a root may fall on either side of a bound by rounding).
template <std::floating_point T, typename F>
void segment_nodes(const Context<T>& ctx, T a, T b, F&& f)
{
    if (ctx.vis_order > 0)
    {
        const T eps = T(1e-5);
        for (int j = 0; j <= ctx.vis_order; ++j)
            f(j, a + (b - a) * (eps + (T(1) - 2 * eps) * T(j) / T(ctx.vis_order)), T(0));
        return;
    }
    for (int j = 0; j < ctx.q; ++j)
        f(j, a + (b - a) * ctx.gauss_x[j], (b - a) * ctx.gauss_w[j]);
}

/// One free coordinate: Gauss-Legendre on the segments between bounds and roots.
template <std::floating_point T, typename Emit>
void integrate_line(Context<T>& ctx, const Problem<T, 1>& p, const Emit& emit)
{
    T L = p.lo[0], U = p.hi[0];
    for (const Half<T, 1>& h : p.clips)
    {
        if (h.c[0] > tiny<T>)
            U = std::min(U, h.d / h.c[0]);
        else if (h.c[0] < -tiny<T>)
            L = std::max(L, h.d / h.c[0]);
        else if (h.d < T(0))
            return;
    }
    if (!(U > L))
        return;
    std::vector<T>& nodes = ctx.nodes[1];
    nodes.clear();
    nodes.push_back(L);
    nodes.push_back(U);
    for (const Func<T, 1>& f : p.funcs)
    {
        if (f.linear)
        {
            if (std::abs(f.a[0]) > tiny<T>)
            {
                const T t = -f.c / f.a[0];
                if (t > L && t < U)
                    nodes.push_back(t);
            }
            continue;
        }
        if (ctx.phi.is_analytic())
        {
            Vec3<T> u0, dir;
            box_line<T, 1>(f, VecD<T, 1>{T(0)}, 0, u0, dir);
            line_roots(ctx.phi, std::span<const T>(u0), std::span<const T>(dir), L, U,
                       zero_function_tol<T> * ctx.phi_scale, nodes);
            continue;
        }
        line_form<T, 1>(ctx, f, VecD<T, 1>{L}, 0, L, U, ctx.line[1]);
        isolate_roots(std::span<const T>(ctx.line[1].coeffs), L, U, nodes, ctx.root_work[1]);
    }
    std::sort(nodes.begin(), nodes.end());
    NodeTag tag;
    tag.box[0] = ctx.box_ids[0]++;
    int segment = 0;
    for (std::size_t s = 0; s + 1 < nodes.size(); ++s)
    {
        const T a = nodes[s], b = nodes[s + 1];
        if (b - a <= segment_tol<T>)
            continue;
        tag.segment[0] = segment++;
        segment_nodes(ctx, a, b,
                      [&](int j, T t, T w)
                      {
                          tag.node[0] = j;
                          emit(VecD<T, 1>{t}, w, tag);
                      });
    }
}

template <std::floating_point T, int D, bool Rotated, typename Emit>
void integrate_box(Context<T>& ctx, Problem<T, D> p, const Emit& emit, int depth, bool interface);

/// The dimension reduction of a tightened and analysed box: bisect if it is not
/// certified, otherwise integrate along height lines over the base level.
template <std::floating_point T, int D, bool Rotated, typename Emit>
void reduce(Context<T>& ctx, Problem<T, D> p, const Analysis<T, D>& an, const Emit& emit, int depth,
            bool interface)
{
    const int k = an.k;
    const bool certified = an.certified;
    if (!certified && ctx.opt.diagnose)
        diagnose_failure<T, D>(ctx, p, an, depth);
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
            const T mid = T(0.5) * (p.lo[axis] + p.hi[axis]);
            Problem<T, D> first = p, second = std::move(p);
            first.hi[axis] = mid;
            second.lo[axis] = mid;
            integrate_box<T, D, Rotated>(ctx, std::move(first), emit, depth + 1, interface);
            integrate_box<T, D, Rotated>(ctx, std::move(second), emit, depth + 1, interface);
            return;
        }
        ++ctx.stats->uncertified;
    }

    // bounds of the height lines: box faces and clip planes
    std::vector<Bound<T, D>> lowers, uppers;
    Problem<T, D - 1> base;
    for (int jb = 0; jb < D - 1; ++jb)
    {
        const int j = jb < k ? jb : jb + 1;
        base.lo[jb] = p.lo[j];
        base.hi[jb] = p.hi[j];
    }
    {
        Bound<T, D> l, u;
        l.alpha = p.lo[k];
        u.alpha = p.hi[k];
        l.kind = 2;
        u.kind = 3;
        lowers.push_back(l);
        uppers.push_back(u);
    }
    for (const Half<T, D>& h : p.clips)
    {
        const T ck = h.c[k];
        if (std::abs(ck) <= tiny<T>)
        {
            Half<T, D - 1> hb;
            for (int jb = 0; jb < D - 1; ++jb)
                hb.c[jb] = h.c[jb < k ? jb : jb + 1];
            hb.d = h.d;
            base.clips.push_back(hb);
            continue;
        }
        Bound<T, D> b;
        b.alpha = h.d / ck;
        b.kind = ck > T(0) ? 5 : 4;
        for (int jb = 0; jb < D - 1; ++jb)
            b.beta[jb] = -h.c[jb < k ? jb : jb + 1] / ck;
        (ck > T(0) ? uppers : lowers).push_back(b);
    }
    // Fourier-Motzkin: the height line is non-empty where every lower <= every upper
    for (std::size_t i = 0; i < lowers.size(); ++i)
        for (std::size_t j = 0; j < uppers.size(); ++j)
        {
            if (i == 0 && j == 0)
                continue;
            Half<T, D - 1> hb;
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
        std::vector<VecD<T, D - 1>> corners;
        polytope_vertices<T, D - 1>(base.lo, base.hi, base.clips, corners);
        const T tol = compare_tol<T>;
        auto prune = [&](std::vector<Bound<T, D>>& list, bool lower)
        {
            std::vector<Bound<T, D>> kept;
            for (std::size_t i = 0; i < list.size(); ++i)
            {
                bool dominated = false;
                for (std::size_t j = 0; j < list.size() && !dominated; ++j)
                {
                    if (j == i)
                        continue;
                    bool all = !corners.empty();
                    for (const VecD<T, D - 1>& y : corners)
                    {
                        const T di = bound_value<T, D>(list[i], y), dj = bound_value<T, D>(list[j], y);
                        // j at least as tight everywhere; ties keep the earlier bound
                        const bool tighter = lower ? (dj > di + tol || (std::abs(dj - di) <= tol && j < i))
                                                   : (dj < di - tol || (std::abs(dj - di) <= tol && j < i));
                        all &= tighter || std::abs(dj - di) <= tol;
                        if (!all)
                            break;
                    }
                    // all corners tied: keep the earlier one
                    if (all)
                    {
                        bool strictly = false;
                        for (const VecD<T, D - 1>& y : corners)
                        {
                            const T di = bound_value<T, D>(list[i], y), dj = bound_value<T, D>(list[j], y);
                            strictly |= lower ? dj > di + tol : dj < di - tol;
                        }
                        dominated = strictly || j < i;
                    }
                }
                if (!dominated)
                    kept.push_back(list[i]);
            }
            list = std::move(kept);
        };
        prune(lowers, true);
        prune(uppers, false);
    }

    // base functions: where the integrand along the line changes form
    std::vector<Bound<T, D>> bounds = lowers;
    bounds.insert(bounds.end(), uppers.begin(), uppers.end());
    std::vector<Bound<T, D>> linear_roots;
    for (const Func<T, D>& f : p.funcs)
        if (f.linear && std::abs(f.a[k]) > tiny<T>)
        {
            Bound<T, D> r;
            r.kind = 6;
            r.alpha = -f.c / f.a[k];
            for (int jb = 0; jb < D - 1; ++jb)
                r.beta[jb] = -f.a[jb < k ? jb : jb + 1] / f.a[k];
            linear_roots.push_back(r);
        }
    for (const Func<T, D>& f : p.funcs)
    {
        if (f.linear && std::abs(f.a[k]) <= tiny<T>)
        {
            base.funcs.push_back(restrict_to<T, D>(f, k, Bound<T, D>{})); // independent of y_k
            continue;
        }
        for (const Bound<T, D>& b : bounds) // a root reaches a bound
            base.funcs.push_back(restrict_to<T, D>(f, k, b));
        if (!f.linear)
            for (const Bound<T, D>& r : linear_roots) // a curved root meets a linear one
                base.funcs.push_back(restrict_to<T, D>(f, k, r));
    }
    for (std::size_t i = 0; i < linear_roots.size(); ++i)
        for (std::size_t j = i + 1; j < linear_roots.size(); ++j)
            base.funcs.push_back(difference<T, D>(linear_roots[i], linear_roots[j]));
    for (std::size_t i = 0; i < lowers.size(); ++i) // the active lower bound changes
        for (std::size_t j = i + 1; j < lowers.size(); ++j)
            base.funcs.push_back(difference<T, D>(lowers[i], lowers[j]));
    for (std::size_t i = 0; i < uppers.size(); ++i)
        for (std::size_t j = i + 1; j < uppers.size(); ++j)
            base.funcs.push_back(difference<T, D>(uppers[i], uppers[j]));

    // the rule along the height lines: the integrand of the base level
    const int box_id = ctx.box_ids[D - 1]++;
    const int dropped_sign = an.dropped_sign;
    const auto line = [&, k, certified, box_id](const VecD<T, D - 1>& yb, T w, const NodeTag& below)
    {
        T L = -infinity<T>, U = infinity<T>;
        for (const Bound<T, D>& b : lowers)
            L = std::max(L, bound_value<T, D>(b, yb));
        for (const Bound<T, D>& b : uppers)
            U = std::min(U, bound_value<T, D>(b, yb));
        if (!(U > L))
            return;
        std::vector<T>& nodes = ctx.nodes[D];
        nodes.clear();
        const VecD<T, D> yL = insert<T, D>(yb, k, L);
        bool curved_line = false; // at the top level, for Bernstein forms: phi on this line is in ctx.line[3]
        for (const Func<T, D>& f : p.funcs)
        {
            if (f.linear)
            {
                if (std::abs(f.a[k]) > tiny<T>)
                {
                    T t = f.c;
                    for (int jb = 0; jb < D - 1; ++jb)
                        t += f.a[jb < k ? jb : jb + 1] * yb[jb];
                    t = -t / f.a[k];
                    if (t > L && t < U)
                        nodes.push_back(t);
                }
                continue;
            }
            curved_line = true;
            if (ctx.phi.is_analytic())
            {
                Vec3<T> u0, dir;
                box_line<T, D>(f, insert<T, D>(yb, k, T(0)), k, u0, dir);
                if (certified)
                {
                    const T gl = value<T, D>(ctx, f, yL), gu = value<T, D>(ctx, f, insert<T, D>(yb, k, U));
                    if (gl != T(0) && gu != T(0) && (gl > T(0)) != (gu > T(0)))
                        nodes.push_back(
                            line_root(ctx.phi, std::span<const T>(u0), std::span<const T>(dir), L, U, gl, gu));
                }
                else
                    line_roots(ctx.phi, std::span<const T>(u0), std::span<const T>(dir), L, U,
                               zero_function_tol<T> * ctx.phi_scale, nodes);
                continue;
            }
            line_form<T, D>(ctx, f, yL, k, L, U, ctx.line[D]);
            const std::span<const T> c(ctx.line[D].coeffs);
            if (certified)
            {
                const T gl = c.front(), gu = c.back();
                if (gl != T(0) && gu != T(0) && (gl > T(0)) != (gu > T(0)))
                    nodes.push_back(bracketed_root(c, L, U, L, U, gl, gu));
            }
            else
                isolate_roots(c, L, U, nodes, ctx.root_work[D]);
        }

        NodeTag tag = below;
        tag.box[D - 1] = box_id;
        if constexpr (D == 3)
        {
            if (interface)
            {
                // the root of phi on each certified line, weighted by the physical
                // surface measure: |det J| |J^-T grad phi| / |d_k phi|
                for (std::size_t r = 0; r < nodes.size(); ++r)
                {
                    tag.segment[D - 1] = static_cast<int>(r);
                    const VecD<T, D> y = insert<T, D>(yb, k, nodes[r]);
                    Vec3<T> g;
                    gradient(ctx.phi, std::span<const T>(y), std::span<T>(g));
                    const T dk = std::abs(g[k]);
                    if (!(dk > T(0)))
                        continue;
                    const T gnorm = scaled_norm(std::span<const T>(g));
                    Vec3<T> gx = {0, 0, 0};
                    for (int i = 0; i < 3; ++i)
                        for (int m = 0; m < 3; ++m)
                            gx[i] += ctx.inv[m][i] * g[m] / gnorm;
                    emit(y, w * gnorm / dk * ctx.detj * scaled_norm(std::span<const T>(gx)), tag);
                }
                return;
            }
        }
        nodes.push_back(L);
        nodes.push_back(U);
        std::sort(nodes.begin(), nodes.end());
        for (std::size_t s = 0; s + 1 < nodes.size(); ++s)
        {
            const T a = nodes[s], b = nodes[s + 1];
            // leaf cells keep degenerate segments so segment numbers stay aligned
            if (ctx.vis_order == 0 && b - a <= segment_tol<T>)
                continue;
            tag.segment[D - 1] = static_cast<int>(s);
            segment_nodes(ctx, a, b,
                          [&](int j, T t, T wt)
                          {
                              tag.node[D - 1] = j;
                              const VecD<T, D> y = insert<T, D>(yb, k, t);
                              if constexpr (D == 3)
                              {
                                  if (ctx.vis_order == 0)
                                  {
                                      // top level: keep the points of the selected side
                                      if (ctx.part == Part::negative || ctx.part == Part::positive)
                                      {
                                          T v;
                                          if (curved_line && !ctx.phi.is_analytic())
                                              v = evaluate_1d(std::span<const T>(ctx.line[3].coeffs), (t - L) / (U - L));
                                          else if (!curved_line && dropped_sign != 0)
                                              v = T(dropped_sign);
                                          else
                                              v = evaluate(ctx.phi, std::span<const T>(y));
                                          if ((ctx.part == Part::negative && !(v < T(0)))
                                              || (ctx.part == Part::positive && !(v > T(0))))
                                              return;
                                      }
                                      emit(y, w * wt * ctx.detj, tag);
                                      return;
                                  }
                              }
                              emit(y, w * wt, tag);
                          });
        }
    };
    if constexpr (D == 2)
        integrate_line<T>(ctx, base, line);
    else
        integrate_box<T, D - 1, false>(ctx, std::move(base), line, 0, false);
}

template <std::floating_point T, int D, bool Rotated, typename Emit>
void integrate_box(Context<T>& ctx, Problem<T, D> p, const Emit& emit, int depth, bool interface)
{
    if (!tighten<T, D>(p.lo, p.hi, p.clips))
        return;

    // keep the functions that may vanish on the box, with their direction margins
    Analysis<T, D> an = analyse<T, D>(ctx, p);
    if (interface && an.funcs.empty())
        return; // the level set does not cross this box
    p.funcs = std::move(an.funcs);

    // Level 2: two functions whose zero curves meet (on a face of a tet, say) may
    // each need a different axis, and no bisection separates them. A diagonal
    // frame often suits both.
    if constexpr (D == 2 && !Rotated)
    {
        if (!an.certified && ctx.opt.diagonal_frames)
        {
            Problem<T, 2> r = rotate_diagonal<T>(p);
            if (tighten<T, 2>(r.lo, r.hi, r.clips))
            {
                Analysis<T, 2> ar = analyse<T, 2>(ctx, r);
                if (ar.certified)
                {
                    ++ctx.stats->rotations;
                    r.funcs = std::move(ar.funcs);
                    const auto back = [&emit](const VecD<T, 2>& z, T w, const NodeTag& tag)
                    { emit(VecD<T, 2>{T(0.5) * (z[0] - z[1]), T(0.5) * (z[0] + z[1])}, T(0.5) * w, tag); };
                    reduce<T, 2, true>(ctx, std::move(r), ar, back, depth, interface);
                    return;
                }
            }
        }
    }
    reduce<T, D, Rotated>(ctx, std::move(p), an, emit, depth, interface);
}

/// The engine on one cell: the top level is the unit box with the cell's clips
/// and phi itself.
template <std::floating_point T>
void run(const ClippedBox<T>& cell, const Source<T>& phi, Part part, int q, int vis_order,
         const Options& opt, CellPoints<T>& out, Stats& stats)
{
    if ((phi.bernstein == nullptr) == (phi.analytic == nullptr))
        throw std::invalid_argument("quadrays: a source is a Bernstein form or an analytic level set");
    if (phi.bernstein != nullptr && phi.bernstein->dim != 3)
        throw std::invalid_argument("quadrays: the level set must be a form in 3 variables");
    if (vis_order == 0 && q < 1)
        throw std::invalid_argument("quadrays: at least one Gauss point per segment is needed");

    thread_local Context<T> ctx;
    ctx.phi = phi;
    ctx.phi_scale = reference_magnitude(phi);
    ctx.opt = opt;
    ctx.margin = static_cast<T>(opt.margin);
    ctx.stats = &stats;
    ctx.vis_order = vis_order;
    ctx.part = part;
    ctx.detj = std::abs(jacobian_determinant(cell));
    ctx.inv = part == Part::interface ? inverse_jacobian(cell) : Mat3<T>{};
    ctx.box_ids = {0, 0, 0};
    ctx.bisections = 0;
    if (vis_order == 0 && ctx.q != q)
    {
        gauss_legendre<T>(q, ctx.gauss_x, ctx.gauss_w);
        ctx.q = q;
    }

    Problem<T, 3> top;
    top.lo = {0, 0, 0};
    top.hi = {1, 1, 1};
    Func<T, 3> f;
    for (int i = 0; i < 3; ++i)
        f.A[i][i] = T(1);
    top.funcs.push_back(f);
    for (const HalfSpace<T>& h : cell.clips)
        top.clips.push_back({h.c, h.d});

    const auto emit = [&](const VecD<T, 3>& u, T w, const NodeTag& tag)
    {
        out.points.insert(out.points.end(), u.begin(), u.end());
        out.weights.push_back(w);
        if (vis_order > 0)
            out.tags.push_back(tag);
    };
    integrate_box<T, 3, false>(ctx, std::move(top), emit, 0, part == Part::interface);
}

} // namespace

template <std::floating_point T>
void integrate(const ClippedBox<T>& cell, const Source<T>& phi, Part part, int q, const Options& opt,
               CellPoints<T>& out, Stats& stats)
{
    run(cell, phi, part, q, 0, opt, out, stats);
}

template <std::floating_point T>
void integrate(const ClippedBox<T>& cell, const BoxBernstein<T>& phi, Part part, int q,
               const Options& opt, CellPoints<T>& out, Stats& stats)
{
    run(cell, bernstein_source(phi), part, q, 0, opt, out, stats);
}

template <std::floating_point T>
void leaf_nodes(const ClippedBox<T>& cell, const Source<T>& phi, Part part, int degree, const Options& opt,
                CellPoints<T>& out, Stats& stats)
{
    if (degree < 1)
        throw std::invalid_argument("quadrays: leaf cells need degree 1 or more");
    run(cell, phi, part, 0, degree, opt, out, stats);
}

template <std::floating_point T>
void leaf_nodes(const ClippedBox<T>& cell, const BoxBernstein<T>& phi, Part part, int degree,
                const Options& opt, CellPoints<T>& out, Stats& stats)
{
    leaf_nodes(cell, bernstein_source(phi), part, degree, opt, out, stats);
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template void integrate<float>(const ClippedBox<float>&, const Source<float>&, Part, int, const Options&,
                               CellPoints<float>&, Stats&);
template void integrate<double>(const ClippedBox<double>&, const Source<double>&, Part, int, const Options&,
                                CellPoints<double>&, Stats&);
template void leaf_nodes<float>(const ClippedBox<float>&, const Source<float>&, Part, int, const Options&,
                                CellPoints<float>&, Stats&);
template void leaf_nodes<double>(const ClippedBox<double>&, const Source<double>&, Part, int, const Options&,
                                 CellPoints<double>&, Stats&);
template void integrate<float>(const ClippedBox<float>&, const BoxBernstein<float>&, Part, int,
                               const Options&, CellPoints<float>&, Stats&);
template void integrate<double>(const ClippedBox<double>&, const BoxBernstein<double>&, Part, int,
                                const Options&, CellPoints<double>&, Stats&);
template void leaf_nodes<float>(const ClippedBox<float>&, const BoxBernstein<float>&, Part, int,
                                const Options&, CellPoints<float>&, Stats&);
template void leaf_nodes<double>(const ClippedBox<double>&, const BoxBernstein<double>&, Part, int,
                                 const Options&, CellPoints<double>&, Stats&);

} // namespace cutcells::quadrays
