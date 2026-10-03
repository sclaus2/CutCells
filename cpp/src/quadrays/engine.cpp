// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "engine.h"

#include <algorithm>
#include <bit>
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

SelectionTerm term_of(Part part, int level_set)
{
    if (level_set < 0 || level_set >= 64)
        throw std::invalid_argument("quadrays: level-set index out of range");
    const std::uint64_t bit = std::uint64_t(1) << level_set;
    SelectionTerm term;
    switch (part)
    {
    case Part::negative:
        term.negative_required = bit;
        break;
    case Part::positive:
        term.positive_required = bit;
        break;
    case Part::interface:
        term.zero_required = bit;
        break;
    default:
        break;
    }
    return term;
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

/// a and b non-zero with opposite signs: a root lies between them.
template <std::floating_point T>
bool opposite(T a, T b)
{
    return a != T(0) && b != T(0) && (a > T(0)) != (b > T(0));
}

// ============================================================================
// Problems at one level
// ============================================================================

/// psi(y) = phi_ls(A y + b) (curved), or a . y + c (linear), on level coordinates y.
template <typename T, int D>
struct Func
{
    std::uint32_t origin = 1; ///< diagnostics: chain of 4-bit codes, see origin_name
    bool linear = false;
    int ls = 0; ///< the level set of a curved function
    std::array<VecD<T, D>, 3> A{};
    Vec3<T> b{};
    VecD<T, D> a{};
    T c = 0;
};

/// A level set restricted to another's zero set along the height lines of the
/// level above: s(y) = on(y, r(y)), with r(y) the root of t -> under(y, t). It
/// vanishes where the roots of under and on cross. The inner functions are
/// curved functions of the level above, whose coordinate k is t; under is
/// monotone in t on the box's height range [t_lo, t_hi]. The root is followed
/// half that range beyond it, so that s stays smooth where the root leaves the
/// box (it would have a kink there if clamped at the range, which no box
/// along that curve could certify); beyond, it is clamped. Crossings outside
/// [t_lo, t_hi] are not in the box: they only add breakpoints.
template <typename T, int D>
struct Surface
{
    std::uint32_t origin = 11;
    Func<T, D + 1> under, on;
    int k = 0;
    T t_lo = 0, t_hi = 0;
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
    std::vector<Surface<T, D>> surfaces;
    std::vector<Half<T, D>> clips;
    /// Top level: the map from the level's coordinates to box coordinates (the
    /// identity, or a rotation of it).
    Func<T, D> frame;
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
/// of lower bounds, 9 switch of upper bounds, 10 difference of linear roots; 11
/// a surface function (one level set on another's zero set).
std::string origin_name(std::uint32_t origin)
{
    static const char* names[] = {"?",        "phi",       "box-lo",    "box-hi",   "clip-lo", "clip-up",
                                  "lin-root", "same",      "switch-lo", "switch-up", "root-diff", "surface",
                                  "corner"};
    std::vector<std::string> parts;
    for (; origin != 0; origin >>= 4)
        parts.push_back(names[std::min<std::uint32_t>(origin & 15u, 12u)]);
    std::string s;
    for (auto it = parts.rbegin(); it != parts.rend(); ++it)
        s += (s.empty() ? "" : "|") + *it;
    return s;
}

/// The sign conditions of a selection term on the cell's level sets.
struct Requirement
{
    std::uint64_t negative = 0, positive = 0, zero = 0;
};

// ============================================================================
// Context: the cell, the options and scratch buffers per level
// ============================================================================

template <std::floating_point T>
struct Context
{
    std::vector<Source<T>> phis; ///< the level sets on the cell's unit box
    std::vector<BoxBernstein<T>> linear_forms; ///< of analytic level sets linear on the cell
    std::vector<T> scales;       ///< size of each (reference_magnitude)
    std::vector<Requirement> terms;
    std::uint64_t used = 0; ///< level sets some term constrains
    int surface = -1;       ///< interface parts: the level set whose zero set is integrated
    Options opt;
    T margin = 0;
    Stats* stats = nullptr;
    int vis_order = 0; ///< > 0: leaf nodes (vis_order + 1 per segment) instead of Gauss points
    T detj = 0;        ///< |det| of the box-to-physical map
    Mat3<T> inv = {};  ///< its inverse (interface weights)
    std::array<int, 3> box_ids{};
    int bisections = 0; ///< bisections so far in this cell

    // Gauss-Legendre rule on [0, 1]
    int q = 0;
    std::vector<T> gauss_x, gauss_w;

    // scratch, indexed by level (number of free coordinates)
    std::array<BoxBernstein<T>, 4> form, deriv, second, sub, line;
    std::array<std::vector<T>, 4> form_work, line_work, root_work, nodes, slope_line;
    // scratch of surface functions, indexed by the level of their inner functions
    std::array<BoxBernstein<T>, 4> inner_form, inner_deriv, inner_line;
    std::array<std::vector<T>, 4> inner_work;
    // the level sets' Bernstein forms on the current top-level height line
    std::vector<BoxBernstein<T>> top_lines;
    std::vector<char> top_line_set;
    std::vector<T> top_values;
};

template <std::floating_point T>
bool is_analytic(const Context<T>& ctx, int ls)
{
    return ctx.phis[static_cast<std::size_t>(ls)].is_analytic();
}

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

/// Box coordinates of a top-level point (the identity map, u2 = 0 on 2D cells).
template <typename T, int D>
Vec3<T> pad(const VecD<T, D>& y)
{
    Vec3<T> u = {0, 0, 0};
    for (int j = 0; j < D; ++j)
        u[j] = y[j];
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
    return evaluate(ctx.phis[static_cast<std::size_t>(f.ls)], std::span<const T>(u));
}

template <std::floating_point T, int D>
T derivative_along(const Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& y, int k)
{
    if (f.linear)
        return f.a[k];
    const Vec3<T> u = to_u<T, D>(f, y);
    Vec3<T> g = {0, 0, 0};
    gradient(ctx.phis[static_cast<std::size_t>(f.ls)], std::span<const T>(u), std::span<T>(g));
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
    r.ls = f.ls;
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

/// A surface function restricted to y_k = alpha + beta . y': its inner
/// functions restricted along the same coordinate.
template <typename T, int D>
Surface<T, D - 1> restrict_surface(const Surface<T, D>& s, int k, const Bound<T, D>& b)
{
    Surface<T, D - 1> r;
    r.origin = (s.origin << 4) | static_cast<std::uint32_t>(b.kind);
    r.t_lo = s.t_lo;
    r.t_hi = s.t_hi;
    // y_k among the inner coordinates insert(y, s.k, t), and t after removing it
    const int K = k < s.k ? k : k + 1;
    r.k = k < s.k ? s.k - 1 : s.k;
    Bound<T, D + 1> inner;
    inner.alpha = b.alpha;
    inner.kind = b.kind;
    for (int jb = 0; jb < D - 1; ++jb)
        inner.beta[jb < r.k ? jb : jb + 1] = b.beta[jb];
    inner.beta[r.k] = T(0);
    r.under = restrict_to<T, D + 1>(s.under, K, inner);
    r.on = restrict_to<T, D + 1>(s.on, K, inner);
    return r;
}

/// A surface function along the line y = insert(yb, k, tau), as a surface
/// function of tau: inner coordinates (tau, t).
template <typename T, int D>
Surface<T, 1> surface_on_line(const Surface<T, D>& s, const VecD<T, D - 1>& yb, int k)
{
    Surface<T, 1> r;
    r.origin = s.origin;
    r.t_lo = s.t_lo;
    r.t_hi = s.t_hi;
    r.k = 1;
    const int K = k < s.k ? k : k + 1; // tau among the inner coordinates
    // the inner coordinates at tau = 0 and t = 0
    const VecD<T, D + 1> y0 = insert<T, D + 1>(insert<T, D>(yb, k, T(0)), s.k, T(0));
    auto project = [&](const Func<T, D + 1>& f)
    {
        Func<T, 2> g;
        g.origin = f.origin;
        g.linear = f.linear;
        g.ls = f.ls;
        if (f.linear)
        {
            g.a = {f.a[K], f.a[s.k]};
            g.c = f.c;
            for (int j = 0; j < D + 1; ++j)
                g.c += f.a[j] * y0[j];
        }
        else
        {
            g.b = to_u<T, D + 1>(f, y0);
            for (int i = 0; i < 3; ++i)
                g.A[i] = {f.A[i][K], f.A[i][s.k]};
        }
        return g;
    };
    r.under = project(s.under);
    r.on = project(s.on);
    return r;
}

// ============================================================================
// Bernstein forms of curved functions
// ============================================================================

/// Bernstein form of a curved function on the box [lo, hi] of its level: exact
/// restriction of its level set to the affine image of the box.
template <std::floating_point T, int D>
void bernstein_form(const Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi,
                    BoxBernstein<T>& out, std::vector<T>& work)
{
    const BoxBernstein<T>& phi = *ctx.phis[static_cast<std::size_t>(f.ls)].bernstein;
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
    // a form of a 2D cell reads the first two rows
    restrict_affine(phi, std::span<const T>(origin.data(), static_cast<std::size_t>(phi.dim)),
                    std::span<const T>(matrix.data(), static_cast<std::size_t>(phi.dim * D)), D, out, work);
}

/// Bernstein form on [L, U] of t -> psi(y0 + (t - L) e_k), y0 on the line at t = L.
template <std::floating_point T, int D>
void line_form(const Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& y0, int k, T L, T U,
               BoxBernstein<T>& out, std::vector<T>& work)
{
    const BoxBernstein<T>& phi = *ctx.phis[static_cast<std::size_t>(f.ls)].bernstein;
    std::array<T, 3> origin, column;
    for (int i = 0; i < 3; ++i)
    {
        T o = f.b[i];
        for (int j = 0; j < D; ++j)
            o += snap(f.A[i][j]) * y0[j];
        origin[i] = o;
        column[i] = snap(f.A[i][k]) * (U - L);
    }
    restrict_affine(phi, std::span<const T>(origin.data(), static_cast<std::size_t>(phi.dim)),
                    std::span<const T>(column.data(), static_cast<std::size_t>(phi.dim)), 1, out, work);
}

/// The line y = y0 + t e_k of a level as u = u0 + t dir in box coordinates.
template <typename T, int D>
void box_line(const Func<T, D>& f, const VecD<T, D>& y0, int k, Vec3<T>& u0, Vec3<T>& dir)
{
    u0 = to_u<T, D>(f, y0);
    for (int i = 0; i < 3; ++i)
        dir[i] = f.A[i][k];
}

/// The affine image of [lo, hi] for an analytic level set: origin and matrix
/// (3 x D, row-major) of the map s -> A (lo + diag(hi - lo) s) + b.
template <typename T, int D>
void affine_image(const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi, std::array<T, 3>& origin,
                  std::array<T, 3 * D>& matrix)
{
    for (int i = 0; i < 3; ++i)
    {
        T o = f.b[i];
        for (int j = 0; j < D; ++j)
        {
            o += f.A[i][j] * lo[j];
            matrix[i * D + j] = f.A[i][j] * (hi[j] - lo[j]);
        }
        origin[i] = o;
    }
}

/// Bounds of a curved function on [lo, hi] for an analytic level set, from
/// Taylor models of phi on the affine image of the box (the map to physical
/// space is affine, so the models represent it exactly). Sets the sign of the
/// function if it is certain, the direction margins as margins() does from
/// Bernstein coefficients, and the largest |psi|. Without derivative bounds
/// no direction is certified; without any bound nothing is certain.
template <std::floating_point T, int D>
void analytic_bounds(const Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi,
                     int& sign, VecD<T, D>& ratio, T& magnitude, VecD<T, D>& upper)
{
    upper.fill(infinity<T>);
    std::array<T, 3> origin;
    std::array<T, 3 * D> matrix;
    affine_image<T, D>(f, lo, hi, origin, matrix);
    VecD<T, D> lengths;
    for (int j = 0; j < D; ++j)
        lengths[j] = hi[j] - lo[j];
    AffineBounds<T> b;
    if (!affine_bounds(ctx.phis[static_cast<std::size_t>(f.ls)], std::span<const T>(origin),
                       std::span<const T>(matrix), D, b))
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
    VecD<T, D> lower;
    for (int j = 0; j < D; ++j)
    {
        lower[j] = b.lower[j] / lengths[j];
        upper[j] = b.upper[j] / lengths[j];
    }
    const T norm = scaled_norm(std::span<const T>(upper));
    for (int j = 0; j < D; ++j)
        ratio[j] = norm > T(0) ? lower[j] / norm : T(0);
}

/// analytic_bounds with margins from Taylor models over the M^D sub-boxes of
/// [lo, hi] that may meet the clipped region. The sign and the size come from
/// the whole box: a zero set on the boundary between sub-boxes touches each of
/// them with one sign, and would be lost. The margins count the sub-boxes on
/// which psi may vanish, as local_margins does for Bernstein forms: if d_k psi
/// has one strict sign on all of them, psi has at most one root on every
/// height line in the region. Without such a sub-box (the zero set on their
/// boundaries) or a bound, the whole box's margins stand.
template <std::floating_point T, int D>
void analytic_local_bounds(const Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi,
                           const std::vector<Half<T, D>>& clips, int M, int& sign, VecD<T, D>& ratio,
                           T& magnitude, VecD<T, D>& upper_out)
{
    analytic_bounds<T, D>(ctx, f, lo, hi, sign, ratio, magnitude, upper_out);
    if (sign != 0 || !std::isfinite(magnitude))
        return;
    int cells = 1;
    for (int j = 0; j < D; ++j)
        cells *= M;
    bool relevant = false;
    std::array<int, D> dsign{};
    VecD<T, D> lower, upper;
    lower.fill(infinity<T>);
    upper.fill(T(0));
    const Source<T>& phi = ctx.phis[static_cast<std::size_t>(f.ls)];
    for (int n = 0; n < cells; ++n)
    {
        VecD<T, D> ylo, yhi;
        for (int j = 0, m = n; j < D; ++j, m /= M)
        {
            ylo[j] = lo[j] + (hi[j] - lo[j]) * T(m % M) / T(M);
            yhi[j] = lo[j] + (hi[j] - lo[j]) * T(m % M + 1) / T(M);
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
        std::array<T, 3> origin;
        std::array<T, 3 * D> matrix;
        affine_image<T, D>(f, ylo, yhi, origin, matrix);
        AffineBounds<T> b;
        if (!affine_bounds(phi, std::span<const T>(origin), std::span<const T>(matrix), D, b) || !b.has_derivatives)
            return; // the whole box's margins stand
        if (b.sign != 0)
            continue; // no root inside (at most on its boundary)
        relevant = true;
        for (int j = 0; j < D; ++j)
        {
            const T length = yhi[j] - ylo[j];
            const int s = b.dmin[j] > T(0) ? 1 : (b.dmax[j] < T(0) ? -1 : 0);
            // one strict sign shared by all counted sub-boxes
            if (s == 0 || (dsign[j] != 0 && dsign[j] != s))
                dsign[j] = 2;
            else if (dsign[j] == 0)
                dsign[j] = s;
            lower[j] = std::min(lower[j], s != 0 ? b.lower[j] / length : T(0));
            upper[j] = std::max(upper[j], b.upper[j] / length);
        }
    }
    if (!relevant)
        return; // the zero set on sub-box boundaries: the whole box's margins
    const T norm = scaled_norm(std::span<const T>(upper));
    VecD<T, D> local;
    for (int j = 0; j < D; ++j)
        local[j] = (norm > T(0) && (dsign[j] == 1 || dsign[j] == -1)) ? lower[j] / norm : T(0);
    ratio = local;
    upper_out = upper;
}

/// Margins from M^D sub-cells. Only sub-cells that may meet the clipped region and
/// on which psi may vanish count. If d_k psi has one strict sign on all of them,
/// psi has at most one root on every height line in the region, and ratio[k] is the
/// smallest local min|d_k psi| / max|grad psi|. Returns false if no sub-cell
/// counts: the function does not vanish in the cell and can be dropped.
template <std::floating_point T, int D>
bool local_margins(Context<T>& ctx, const BoxBernstein<T>& p, const VecD<T, D>& lo, const VecD<T, D>& hi,
                   const std::vector<Half<T, D>>& clips, int M, VecD<T, D>& ratio, VecD<T, D>& upper)
{
    upper.fill(T(0));
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
            upper[k] = std::max(upper[k], dmax[k]);
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
// Bounds of curved functions and surface functions as intervals
// ============================================================================

/// Bounds of a curved function over a box: the sign of its value if certain,
/// and intervals of its derivatives with respect to the level's coordinates.
template <typename T, int D>
struct Ranges
{
    bool valid = false; ///< false: no bound holds
    int sign = 0;
    bool has_derivatives = false;
    VecD<T, D> dmin{}, dmax{};
};

template <std::floating_point T, int D>
Ranges<T, D> curved_ranges(Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi,
                           BoxBernstein<T>& form, BoxBernstein<T>& deriv, std::vector<T>& work)
{
    Ranges<T, D> r;
    if (is_analytic(ctx, f.ls))
    {
        std::array<T, 3> origin;
        std::array<T, 3 * D> matrix;
        affine_image<T, D>(f, lo, hi, origin, matrix);
        AffineBounds<T> b;
        if (!affine_bounds(ctx.phis[static_cast<std::size_t>(f.ls)], std::span<const T>(origin),
                           std::span<const T>(matrix), D, b))
            return r;
        r.valid = true;
        r.sign = b.sign;
        r.has_derivatives = b.has_derivatives;
        for (int j = 0; j < D && r.has_derivatives; ++j)
        {
            r.dmin[j] = b.dmin[j] / (hi[j] - lo[j]);
            r.dmax[j] = b.dmax[j] / (hi[j] - lo[j]);
        }
        return r;
    }
    bernstein_form<T, D>(ctx, f, lo, hi, form, work);
    r.valid = true;
    r.sign = coefficient_sign(std::span<const T>(form.coeffs));
    r.has_derivatives = true;
    for (int j = 0; j < D; ++j)
    {
        if (form.degree[j] == 0)
            continue; // constant along y_j
        derivative(form, j, deriv);
        T vmin = infinity<T>, vmax = -infinity<T>;
        for (const T v : deriv.coeffs)
        {
            vmin = std::min(vmin, v);
            vmax = std::max(vmax, v);
        }
        r.dmin[j] = vmin / (hi[j] - lo[j]);
        r.dmax[j] = vmax / (hi[j] - lo[j]);
    }
    return r;
}

template <std::floating_point T>
std::array<T, 2> interval_mul(T a, T b, T c, T d)
{
    const T p[4] = {a * c, a * d, b * c, b * d};
    return {*std::min_element(p, p + 4), *std::max_element(p, p + 4)};
}

/// The range in which a surface function follows the root of under.
template <typename T, int D>
std::array<T, 2> surface_range(const Surface<T, D>& s)
{
    const T e = T(0.5) * (s.t_hi - s.t_lo);
    return {s.t_lo - e, s.t_hi + e};
}

/// under at the inner point Y.
template <std::floating_point T, int D>
T under_value(const Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D + 1>& Y)
{
    return value<T, D + 1>(ctx, s.under, Y);
}

/// The root of t -> under(y, t) in [a, b], where its values @p fa and @p fb
/// have opposite signs.
template <std::floating_point T, int D>
T under_root(Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D>& y, T a, T b, T fa, T fb)
{
    if (is_analytic(ctx, s.under.ls))
    {
        Vec3<T> u0, dir;
        box_line<T, D + 1>(s.under, insert<T, D + 1>(y, s.k, T(0)), s.k, u0, dir);
        return line_root(ctx.phis[static_cast<std::size_t>(s.under.ls)], std::span<const T>(u0),
                         std::span<const T>(dir), a, b, fa, fb);
    }
    BoxBernstein<T>& line = ctx.inner_line[D + 1];
    line_form<T, D + 1>(ctx, s.under, insert<T, D + 1>(y, s.k, a), s.k, a, b, line, ctx.inner_work[D + 1]);
    return bracketed_root(std::span<const T>(line.coeffs), a, b, a, b, fa, fb);
}

/// The root of t -> under(y, t): in [t_lo, t_hi], where under is monotone,
/// else beyond the end where |under| is smaller, within surface_range, else
/// that end of the range.
template <std::floating_point T, int D>
T surface_root(Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D>& y)
{
    const T flo = under_value<T, D>(ctx, s, insert<T, D + 1>(y, s.k, s.t_lo)),
            fhi = under_value<T, D>(ctx, s, insert<T, D + 1>(y, s.k, s.t_hi));
    if (flo == T(0))
        return s.t_lo;
    if (fhi == T(0))
        return s.t_hi;
    if (opposite(flo, fhi))
        return under_root<T, D>(ctx, s, y, s.t_lo, s.t_hi, flo, fhi);
    const std::array<T, 2> range = surface_range(s);
    if (std::abs(flo) <= std::abs(fhi))
    {
        const T fa = under_value<T, D>(ctx, s, insert<T, D + 1>(y, s.k, range[0]));
        return opposite(fa, flo) ? under_root<T, D>(ctx, s, y, range[0], s.t_lo, fa, flo) : range[0];
    }
    const T fb = under_value<T, D>(ctx, s, insert<T, D + 1>(y, s.k, range[1]));
    return opposite(fhi, fb) ? under_root<T, D>(ctx, s, y, s.t_hi, range[1], fhi, fb) : range[1];
}

/// The value of on at the root of under: the surface function at y.
template <std::floating_point T, int D>
T surface_value(Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D>& y)
{
    return value<T, D + 1>(ctx, s.on, insert<T, D + 1>(y, s.k, surface_root<T, D>(ctx, s, y)));
}

/// The gradient of a surface function at y: d_j on - d_t on d_j under / d_t
/// under at the root, or d_j on where the root is clamped.
template <std::floating_point T, int D>
void surface_gradient(Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D>& y, VecD<T, D>& g)
{
    const T r = surface_root<T, D>(ctx, s, y);
    const VecD<T, D + 1> yr = insert<T, D + 1>(y, s.k, r);
    const std::array<T, 2> range = surface_range(s);
    const T ut = derivative_along<T, D + 1>(ctx, s.under, yr, s.k);
    const T ot = derivative_along<T, D + 1>(ctx, s.on, yr, s.k);
    const bool interior = r > range[0] && r < range[1] && ut != T(0);
    for (int j = 0; j < D; ++j)
    {
        const int J = j < s.k ? j : j + 1;
        g[j] = derivative_along<T, D + 1>(ctx, s.on, yr, J);
        if (interior)
            g[j] -= ot * derivative_along<T, D + 1>(ctx, s.under, yr, J) / ut;
    }
}

/// What bounds say about a surface function on a box.
enum class SurfaceCase
{
    none,     ///< it has no root in the box (its on function has one sign), or
              ///< under has none, so its zeros are those of on on a box face
    bounded,  ///< it may vanish; dmin, dmax bound its derivatives
    unbounded ///< it may vanish; no derivative bound
};

/// Bounds of a surface function on [lo, hi] from those of its inner
/// functions. The root r of under lies, over the box, within r(centre) plus
/// its slope bound |d_j r| <= max |d_j under| / min |d_t under| (bounds on
/// [lo, hi] x [t_lo, t_hi], or on a slab about r(centre) that holds the root
/// where d_t under changes sign on that range) times the half widths; the
/// inner functions are then bounded on that slab, which bisection of the base
/// thins. Where r is
/// interior, d_j s = d_j on - d_t on d_j under / d_t under; where it is
/// clamped, d_j s = d_j on. Intervals of both enclose every one-sided
/// derivative. @p centre_value receives s at the centre.
template <std::floating_point T, int D>
SurfaceCase surface_bounds(Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D>& lo, const VecD<T, D>& hi,
                           VecD<T, D>& dmin, VecD<T, D>& dmax, T& centre_value)
{
    VecD<T, D> centre;
    for (int j = 0; j < D; ++j)
        centre[j] = T(0.5) * (lo[j] + hi[j]);
    const T rc = surface_root<T, D>(ctx, s, centre);
    centre_value = value<T, D + 1>(ctx, s.on, insert<T, D + 1>(centre, s.k, rc));

    // the slab of the root over the box
    const std::array<T, 2> range = surface_range(s);
    VecD<T, D + 1> ilo = insert<T, D + 1>(lo, s.k, range[0]), ihi = insert<T, D + 1>(hi, s.k, range[1]);
    const auto under_ranges = [&](const VecD<T, D + 1>& a, const VecD<T, D + 1>& b)
    {
        return curved_ranges<T, D + 1>(ctx, s.under, a, b, ctx.inner_form[D + 1], ctx.inner_deriv[D + 1],
                                       ctx.inner_work[D + 1]);
    };
    // how far the root strays from rc over the box, from bounds of under on a
    // slab where d_t under has one sign: |d_j r| <= max |d_j under| / min |d_t under|
    const auto slab_reach = [&](const Ranges<T, D + 1>& r, T& reach)
    {
        const T wl = r.dmin[s.k], wh = r.dmax[s.k];
        if (!(wl > T(0) || wh < T(0)))
            return false;
        const T slope_min = std::min(std::abs(wl), std::abs(wh));
        reach = T(0);
        for (int j = 0; j < D; ++j)
        {
            const int J = j < s.k ? j : j + 1;
            reach += std::max(std::abs(r.dmin[J]), std::abs(r.dmax[J])) / slope_min * T(0.5) * (hi[j] - lo[j]);
        }
        return true;
    };
    const Ranges<T, D + 1> wide = under_ranges(ilo, ihi);
    if (wide.valid && wide.sign != 0)
        return SurfaceCase::none; // no root: zeros of on beyond the box
    if (!wide.valid)
        return SurfaceCase::unbounded;
    T reach = T(0);
    if (!wide.has_derivatives || !slab_reach(wide, reach))
    {
        // d_t under changes sign on the range (the box above was certified on
        // its parts where under may vanish only, Taylor sub-boxes), or its
        // derivatives have no bounds there (a distance function, whose Taylor
        // models fail near its centre): a slab about rc on which d_t under has
        // one sign holds the root over the box if the reach its bounds give is
        // at most its half width. Grow it from a thin one.
        if (!(rc > range[0] && rc < range[1]))
        {
            // no root at the centre: none if under keeps one sign on thinner
            // slabs of the range, whose bounds are tighter
            const int pieces = 8;
            int sign = 0;
            for (int i = 0; i < pieces; ++i)
            {
                ilo[s.k] = range[0] + (range[1] - range[0]) * T(i) / T(pieces);
                ihi[s.k] = range[0] + (range[1] - range[0]) * T(i + 1) / T(pieces);
                const Ranges<T, D + 1> piece = under_ranges(ilo, ihi);
                if (!piece.valid || piece.sign == 0 || (sign != 0 && piece.sign != sign))
                    return SurfaceCase::unbounded;
                sign = piece.sign;
            }
            return SurfaceCase::none;
        }
        T half = std::max(T(1e-12), tiny<T>) * (s.t_hi - s.t_lo);
        bool held = false;
        for (int it = 0; it < 6 && !held; ++it)
        {
            ilo[s.k] = std::max(range[0], rc - half);
            ihi[s.k] = std::min(range[1], rc + half);
            const Ranges<T, D + 1> local = under_ranges(ilo, ihi);
            T r = T(0);
            if (!local.valid || !local.has_derivatives || !slab_reach(local, r))
                return SurfaceCase::unbounded;
            held = r <= half;
            if (held)
                reach = r;
            else
                half = T(2) * r;
        }
        if (!held)
            return SurfaceCase::unbounded;
    }
    // a root outside the box's range everywhere: its crossings are not in the box
    if (!(rc - reach <= s.t_hi && rc + reach >= s.t_lo))
        return SurfaceCase::none;
    const T rlo = std::max(range[0], rc - reach), rhi = std::min(range[1], rc + reach);
    // the root stays inside the range it is followed in: no clamping
    const bool interior = rc - reach > range[0] && rc + reach < range[1];
    ilo[s.k] = rlo;
    ihi[s.k] = std::max(rhi, rlo + std::max(T(1e-12), tiny<T>) * (s.t_hi - s.t_lo));

    const Ranges<T, D + 1> on = curved_ranges<T, D + 1>(ctx, s.on, ilo, ihi, ctx.inner_form[D + 1],
                                                        ctx.inner_deriv[D + 1], ctx.inner_work[D + 1]);
    if (on.valid && on.sign != 0)
        return SurfaceCase::none;
    const Ranges<T, D + 1> under = under_ranges(ilo, ihi);
    if (!on.valid || !under.valid || !on.has_derivatives || !under.has_derivatives)
        return SurfaceCase::unbounded;
    const T ul = under.dmin[s.k], uh = under.dmax[s.k];
    if (!(ul > T(0) || uh < T(0)))
        return SurfaceCase::unbounded;
    for (int j = 0; j < D; ++j)
    {
        const int J = j < s.k ? j : j + 1;
        // d_j under / d_t under, times d_t on
        const std::array<T, 2> q = interval_mul(under.dmin[J], under.dmax[J], T(1) / uh, T(1) / ul);
        const std::array<T, 2> p = interval_mul(on.dmin[s.k], on.dmax[s.k], q[0], q[1]);
        dmin[j] = on.dmin[J] - p[1];
        dmax[j] = on.dmax[J] - p[0];
        if (!interior)
        {
            dmin[j] = std::min(dmin[j], on.dmin[J]);
            dmax[j] = std::max(dmax[j], on.dmax[J]);
        }
    }
    return SurfaceCase::bounded;
}

/// Margins of a surface function on [lo, hi] as margins() gives them from
/// derivative bounds; false if it has no root there. With derivative bounds,
/// its value at the centre and the mean-value bound may show one sign.
template <std::floating_point T, int D>
bool surface_margins(Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D>& lo, const VecD<T, D>& hi,
                     VecD<T, D>& ratio)
{
    ratio.fill(T(0));
    VecD<T, D> dmin{}, dmax{};
    T centre_value = T(0);
    const SurfaceCase sc = surface_bounds<T, D>(ctx, s, lo, hi, dmin, dmax, centre_value);
    if (sc == SurfaceCase::none)
        return false;
    if (sc == SurfaceCase::unbounded)
        return true;
    VecD<T, D> upper;
    T reach = T(0);
    for (int j = 0; j < D; ++j)
    {
        upper[j] = std::max(std::abs(dmin[j]), std::abs(dmax[j]));
        reach += upper[j] * T(0.5) * (hi[j] - lo[j]);
    }
    if (std::abs(centre_value) > reach)
        return false;
    const T norm = scaled_norm(std::span<const T>(upper));
    for (int j = 0; j < D; ++j)
    {
        const bool strict = dmin[j] > T(0) || dmax[j] < T(0);
        ratio[j] = (strict && norm > T(0)) ? std::min(std::abs(dmin[j]), std::abs(dmax[j])) / norm : T(0);
    }
    return true;
}

/// Roots of a surface function of one coordinate on (a, b), appended to
/// @p roots: none where bounds show one sign, one where they show a monotone
/// function whose ends differ in sign, bisection otherwise (down to 2^-24 of
/// the interval, at most @p budget intervals).
template <std::floating_point T>
void surface_line_roots(Context<T>& ctx, const Surface<T, 1>& s, T a, T b, T sa, T sb, int depth, int& budget,
                        std::vector<T>& roots)
{
    --budget;
    VecD<T, 1> dmin{}, dmax{};
    T centre_value = T(0);
    const SurfaceCase sc = surface_bounds<T, 1>(ctx, s, VecD<T, 1>{a}, VecD<T, 1>{b}, dmin, dmax, centre_value);
    if (sc == SurfaceCase::none)
        return;
    if (sc == SurfaceCase::bounded
        && std::abs(centre_value) > std::max(std::abs(dmin[0]), std::abs(dmax[0])) * T(0.5) * (b - a))
        return; // one sign on (a, b)
    const bool monotone = sc == SurfaceCase::bounded && (dmin[0] > T(0) || dmax[0] < T(0));
    if (monotone || depth >= 24 || budget <= 0)
    {
        if (opposite(sa, sb))
            roots.push_back(illinois_root([&](T t) { return surface_value<T, 1>(ctx, s, VecD<T, 1>{t}); }, a, b, sa,
                                          sb));
        return;
    }
    const T m = T(0.5) * (a + b);
    const T sm = surface_value<T, 1>(ctx, s, VecD<T, 1>{m});
    if (sm == T(0))
        roots.push_back(m);
    surface_line_roots<T>(ctx, s, a, m, sa, sm, depth + 1, budget, roots);
    surface_line_roots<T>(ctx, s, m, b, sm, sb, depth + 1, budget, roots);
}

template <std::floating_point T>
void surface_line_roots(Context<T>& ctx, const Surface<T, 1>& s, T a, T b, std::vector<T>& roots)
{
    int budget = 256;
    const T sa = surface_value<T, 1>(ctx, s, VecD<T, 1>{a}), sb = surface_value<T, 1>(ctx, s, VecD<T, 1>{b});
    surface_line_roots<T>(ctx, s, a, b, sa, sb, 0, budget, roots);
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

/// True if a clip plane cuts the box [lo, hi].
template <std::floating_point T, int D>
bool box_clipped(const VecD<T, D>& lo, const VecD<T, D>& hi, const std::vector<Half<T, D>>& clips)
{
    for (const Half<T, D>& h : clips)
    {
        T vmax = T(0);
        for (int j = 0; j < D; ++j)
            vmax += h.c[j] * (h.c[j] >= T(0) ? hi[j] : lo[j]);
        if (vmax > h.d + compare_tol<T>)
            return true;
    }
    return false;
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
    std::vector<Surface<T, D>> surfaces;
    std::vector<VecD<T, D>> ratios;         ///< margins of each curved function in curved
    std::vector<VecD<T, D>> uppers;         ///< upper bounds of its |d_j psi|
    std::vector<VecD<T, D>> surface_ratios; ///< margins of each surface function
    std::vector<int> curved;                ///< index in funcs of each entry of ratios
    std::vector<char> two_roots;            ///< per entry of funcs: two roots per line along k
    std::vector<char> independent;          ///< per entry of funcs: constant along k
    /// per entry of funcs: top level of an interface part, another level set
    /// than the interface's; only evaluated at the interface's roots, so it
    /// needs no margin along k
    std::vector<char> passive;
    int k = 0;
    T best = -1;
    bool certified = false;
    /// Top level: each level set's sign on the box (0: may vanish or unknown),
    /// the terms left once those signs are known, and whether they decide the
    /// box (+1 all of it, -1 none of it, 0 the level sets decide).
    std::vector<int> signs;
    std::vector<Requirement> terms;
    int decided = 0;
};

/// The terms that the level sets' signs on a box leave: terms contradicted by
/// a known sign go, conditions met by one are dropped.
inline int residual_terms(const std::vector<Requirement>& terms, const std::vector<int>& signs, int surface,
                          bool surface_vanishes, std::vector<Requirement>& out)
{
    out.clear();
    if (surface >= 0 && !surface_vanishes)
        return -1; // the zero set does not cross the box
    for (const Requirement& t : terms)
    {
        Requirement r;
        bool holds = true;
        for (std::uint64_t bits = t.negative | t.positive; bits != 0 && holds; bits &= bits - 1)
        {
            const int l = std::countr_zero(bits);
            const std::uint64_t bit = std::uint64_t(1) << l;
            const int s = signs[static_cast<std::size_t>(l)];
            const int wanted = (t.negative & bit) ? -1 : 1;
            if (s == 0)
                (wanted < 0 ? r.negative : r.positive) |= bit;
            else
                holds = s == wanted;
        }
        if (!holds)
            continue;
        r.zero = t.zero;
        if ((r.negative | r.positive) == 0)
        {
            out.assign(1, r);
            return 1;
        }
        out.push_back(r);
    }
    return out.empty() ? -1 : 0;
}

/// Second-derivative margin of a curved function along k on [lo, hi]:
/// min |d_k d_k psi| / max |grad d_k psi|, from Bernstein coefficients or
/// Hessian bounds; 0 if d_k d_k psi may vanish or nothing bounds it.
template <std::floating_point T, int D>
T second_margin(Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi, int k,
                int& curvature)
{
    curvature = 0;
    VecD<T, D> lower{}, upper{};
    if (is_analytic(ctx, f.ls))
    {
        std::array<T, 3> origin;
        std::array<T, 3 * D> matrix;
        affine_image<T, D>(f, lo, hi, origin, matrix);
        std::array<T, 9> hlo, hhi;
        if (!affine_hessian(ctx.phis[static_cast<std::size_t>(f.ls)], std::span<const T>(origin),
                            std::span<const T>(matrix), D, hlo, hhi))
            return T(0);
        for (int j = 0; j < D; ++j)
        {
            const T scale = T(1) / ((hi[j] - lo[j]) * (hi[k] - lo[k]));
            const T a = hlo[k * D + j] * scale, b = hhi[k * D + j] * scale;
            upper[j] = std::max(std::abs(a), std::abs(b));
            if (j == k)
            {
                curvature = a > T(0) ? 1 : (b < T(0) ? -1 : 0);
                lower[j] = curvature != 0 ? std::min(std::abs(a), std::abs(b)) : T(0);
            }
        }
    }
    else
    {
        bernstein_form<T, D>(ctx, f, lo, hi, ctx.form[D], ctx.form_work[D]);
        if (ctx.form[D].degree[k] < 2)
            return T(0);
        derivative(ctx.form[D], k, ctx.deriv[D]);
        for (int j = 0; j < D; ++j)
        {
            if (ctx.deriv[D].degree[j] == 0)
                continue;
            derivative(ctx.deriv[D], j, ctx.second[D]);
            bool pos = true, neg = true;
            T amin = infinity<T>, amax = T(0);
            for (const T v : ctx.second[D].coeffs)
            {
                pos &= v > T(0);
                neg &= v < T(0);
                amin = std::min(amin, std::abs(v));
                amax = std::max(amax, std::abs(v));
            }
            const T scale = T(1) / ((hi[j] - lo[j]) * (hi[k] - lo[k]));
            upper[j] = amax * scale;
            if (j == k)
            {
                curvature = pos ? 1 : (neg ? -1 : 0);
                lower[j] = curvature != 0 ? amin * scale : T(0);
            }
        }
    }
    const T norm = scaled_norm(std::span<const T>(upper));
    return norm > T(0) && curvature != 0 ? lower[k] / norm : T(0);
}

/// The sign of psi on the plane y_k = beta . y' + alpha over the base box
/// [blo, bhi] (Bernstein coefficients, strictly, or Taylor models); 0 if the
/// bounds do not show one.
template <std::floating_point T, int D>
int plane_sign(Context<T>& ctx, const Func<T, D>& f, int k, const Bound<T, D>& plane, const VecD<T, D - 1>& blo,
               const VecD<T, D - 1>& bhi)
{
    const Func<T, D - 1> slice = restrict_to<T, D>(f, k, plane);
    if (is_analytic(ctx, f.ls))
    {
        std::array<T, 3> origin;
        std::array<T, 3 * (D - 1)> matrix;
        affine_image<T, D - 1>(slice, blo, bhi, origin, matrix);
        AffineBounds<T> b;
        if (!affine_bounds(ctx.phis[static_cast<std::size_t>(f.ls)], std::span<const T>(origin),
                           std::span<const T>(matrix), D - 1, b))
            return 0;
        return b.sign;
    }
    bernstein_form<T, D - 1>(ctx, slice, blo, bhi, ctx.form[D - 1], ctx.form_work[D - 1]);
    bool pos = true, neg = true;
    for (const T v : ctx.form[D - 1].coeffs)
    {
        pos &= v > T(0);
        neg &= v < T(0);
    }
    return pos ? 1 : (neg ? -1 : 0);
}

/// Two roots per line along k for a function psi whose derivative along k is
/// monotone (second-derivative margin; curvature its sign), judged from the
/// roots themselves rather than by bisection. On a grid of height lines over
/// the base, the extreme point of psi (the root of d_k psi, or the box end
/// nearer to it) lies between the two sheets; the centre line must cross
/// both sheets inside the box (else psi has one root there, a case for the
/// margin: -1). Over each cell of the grid, a
/// plane through the extreme points of its corners does too. If psi has the
/// sign -curvature on every such plane (Bernstein coefficients or Taylor
/// models), the extreme value of psi on every height line has that sign, so
/// its two roots never merge in the box. Returns the smallest
/// |d_k psi| / |grad psi| at the roots on the lines through the base's
/// corners and centre (1 if none; 0 where another axis is much better there,
/// near a silhouette along k), or -1 if a plane does not show the sign.
template <std::floating_point T, int D>
T two_root_score(Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi, int k,
                 int curvature)
{
    // |d_k psi| / |grad psi| at a root, or 0 where another axis is much better
    auto slope_ratio = [&](const VecD<T, D>& y)
    {
        VecD<T, D> g;
        for (int j = 0; j < D; ++j)
            g[j] = derivative_along<T, D>(ctx, f, y, j);
        const T norm = scaled_norm(std::span<const T>(g));
        if (!(norm > T(0)))
            return T(0);
        T best = T(0);
        for (int j = 0; j < D; ++j)
            best = std::max(best, std::abs(g[j]));
        return std::abs(g[k]) >= T(0.75) * best ? std::abs(g[k]) / norm : T(0);
    };
    auto line_point = [&](const VecD<T, D - 1>& yb, T t)
    {
        VecD<T, D> y = insert<T, D>(yb, k, t);
        return y;
    };
    // the extreme point on the line through yb, and the roots on either side
    // (their slopes lower score)
    T score = T(1);
    int found = 0; // roots found on the last line examined
    auto extreme_on = [&](const VecD<T, D - 1>& yb, bool roots) -> T
    {
        found = 0;
        auto slope = [&](T t) { return derivative_along<T, D>(ctx, f, line_point(yb, t), k); };
        auto val = [&](T t) { return value<T, D>(ctx, f, line_point(yb, t)); };
        const T dl = slope(lo[k]), dh = slope(hi[k]);
        // without an extreme point in the box, the end nearer to it (a NaN
        // marks the line for the caller)
        if (!opposite(dl, dh))
            return std::numeric_limits<T>::quiet_NaN();
        const T c = illinois_root(slope, lo[k], hi[k], dl, dh);
        if (roots)
        {
            const std::array<T, 3> ends = {lo[k], c, hi[k]};
            for (int piece = 0; piece < 2; ++piece)
            {
                const T a = ends[static_cast<std::size_t>(piece)], b = ends[static_cast<std::size_t>(piece + 1)];
                if (!(b > a))
                    continue;
                const T fa = val(a), fb = val(b);
                if (opposite(fa, fb))
                {
                    score = std::min(score, slope_ratio(line_point(yb, illinois_root(val, a, b, fa, fb))));
                    ++found;
                }
            }
        }
        return c;
    };
    // the grid of lines over the base: M cells per base direction
    constexpr int M = 4;
    constexpr int nodes = D == 3 ? (M + 1) * (M + 1) : M + 1;
    VecD<T, D - 1> blo, bhi;
    for (int jb = 0; jb < D - 1; ++jb)
    {
        blo[static_cast<std::size_t>(jb)] = lo[jb < k ? jb : jb + 1];
        bhi[static_cast<std::size_t>(jb)] = hi[jb < k ? jb : jb + 1];
    }
    std::array<VecD<T, D - 1>, nodes> base{};
    std::array<T, nodes> extreme{};
    for (int n = 0; n < nodes; ++n)
    {
        int m = n;
        for (int jb = 0; jb < D - 1; ++jb, m /= (M + 1))
            base[static_cast<std::size_t>(n)][static_cast<std::size_t>(jb)]
                = blo[static_cast<std::size_t>(jb)]
                  + (bhi[static_cast<std::size_t>(jb)] - blo[static_cast<std::size_t>(jb)]) * T(m % (M + 1)) / T(M);
        // roots on the lines through the base's corners and centre
        bool corner = true;
        m = n;
        for (int jb = 0; jb < D - 1; ++jb, m /= (M + 1))
            corner &= m % (M + 1) == 0 || m % (M + 1) == M;
        bool centre = true;
        m = n;
        for (int jb = 0; jb < D - 1; ++jb, m /= (M + 1))
            centre &= m % (M + 1) == M / 2;
        T c = extreme_on(base[static_cast<std::size_t>(n)], corner || centre);
        // two sheets: both crossed by the centre line inside the box; else one
        // root there, a case for the margin
        if (centre && (!std::isfinite(c) || found < 2))
            return T(-1);
        if (!std::isfinite(c))
        {
            const VecD<T, D - 1>& yb = base[static_cast<std::size_t>(n)];
            const T dl = derivative_along<T, D>(ctx, f, line_point(yb, lo[k]), k),
                    dh = derivative_along<T, D>(ctx, f, line_point(yb, hi[k]), k);
            c = std::abs(dl) <= std::abs(dh) ? lo[k] : hi[k];
        }
        extreme[static_cast<std::size_t>(n)] = c;
    }
    // a plane per grid cell through the extreme points of three of its corners
    constexpr int cells = D == 3 ? M * M : M;
    for (int cell = 0; cell < cells; ++cell)
    {
        const int i = cell % M, j = D == 3 ? cell / M : 0;
        const int n00 = D == 3 ? i + (M + 1) * j : i;
        VecD<T, D - 1> clo = base[static_cast<std::size_t>(n00)], chi;
        Bound<T, D> plane;
        plane.kind = 7;
        if constexpr (D == 3)
        {
            const int n10 = n00 + 1, n01 = n00 + (M + 1);
            chi = {base[static_cast<std::size_t>(n10)][0], base[static_cast<std::size_t>(n01)][1]};
            plane.beta[0] = (extreme[static_cast<std::size_t>(n10)] - extreme[static_cast<std::size_t>(n00)])
                            / (chi[0] - clo[0]);
            plane.beta[1] = (extreme[static_cast<std::size_t>(n01)] - extreme[static_cast<std::size_t>(n00)])
                            / (chi[1] - clo[1]);
        }
        else
        {
            const int n10 = n00 + 1;
            chi = {base[static_cast<std::size_t>(n10)][0]};
            plane.beta[0] = (extreme[static_cast<std::size_t>(n10)] - extreme[static_cast<std::size_t>(n00)])
                            / (chi[0] - clo[0]);
        }
        plane.alpha = extreme[static_cast<std::size_t>(n00)];
        for (int jb = 0; jb < D - 1; ++jb)
            plane.alpha -= plane.beta[static_cast<std::size_t>(jb)] * clo[static_cast<std::size_t>(jb)];
        if (plane_sign<T, D>(ctx, f, k, plane, clo, chi) != -curvature)
            return T(-1);
    }
    return score;
}

template <std::floating_point T, int D>
Analysis<T, D> analyse(Context<T>& ctx, const Problem<T, D>& p, bool top, int depth)
{
    Analysis<T, D> an;
    BoxBernstein<T>& form = ctx.form[D];
    VecD<T, D> lengths;
    for (int j = 0; j < D; ++j)
        lengths[j] = p.hi[j] - p.lo[j];
    if (top)
        an.signs.assign(ctx.phis.size(), 0);
    const bool clipped = ctx.opt.taylor_subdivisions > 1 && box_clipped<T, D>(p.lo, p.hi, p.clips);
    bool surface_vanishes = false;
    for (const Func<T, D>& f : p.funcs)
    {
        if (f.linear)
        {
            if (linear_may_vanish<T, D>(f, p.lo, p.hi))
                an.funcs.push_back(f);
            continue;
        }
        const T scale = ctx.scales[static_cast<std::size_t>(f.ls)];
        if (is_analytic(ctx, f.ls))
        {
            int sign = 0;
            VecD<T, D> ratio{}, upper{};
            T magnitude = T(0);
            if (clipped)
                analytic_local_bounds<T, D>(ctx, f, p.lo, p.hi, p.clips, ctx.opt.taylor_subdivisions, sign, ratio,
                                            magnitude, upper);
            else
                analytic_bounds<T, D>(ctx, f, p.lo, p.hi, sign, ratio, magnitude, upper);
            // phi = 0 on a face, say: no root inside the region
            if (magnitude <= zero_function_tol<T> * scale)
                continue;
            if (sign != 0)
            {
                if (top)
                    an.signs[static_cast<std::size_t>(f.ls)] = sign;
                continue;
            }
            an.curved.push_back(static_cast<int>(an.funcs.size()));
            an.funcs.push_back(f);
            an.ratios.push_back(ratio);
            an.uppers.push_back(upper);
            surface_vanishes |= f.ls == ctx.surface;
            continue;
        }
        bernstein_form<T, D>(ctx, f, p.lo, p.hi, form, ctx.form_work[D]);
        const std::span<const T> coeffs(form.coeffs);
        // phi vanishes identically here (phi = 0 on a face, say): only rounding
        // noise is left, which would block certification; it has no root inside
        if (max_abs(coeffs) <= zero_function_tol<T> * scale)
            continue;
        VecD<T, D> ratio{}, upper{};
        if (ctx.opt.mask_subdivisions > 1)
        {
            if (!local_margins<T, D>(ctx, form, p.lo, p.hi, p.clips, ctx.opt.mask_subdivisions, ratio, upper))
                continue; // no zero of this function inside the cell
        }
        else
        {
            if (!may_vanish(coeffs))
            {
                // one sign on the box, zeros at most where it touches zero
                const T largest = *std::max_element(coeffs.begin(), coeffs.end(), [](T a, T b)
                                                    { return std::abs(a) < std::abs(b); });
                if (top)
                    an.signs[static_cast<std::size_t>(f.ls)] = largest > T(0) ? 1 : -1;
                continue;
            }
            margins(form, std::span<const T>(lengths), std::span<T>(ratio), ctx.deriv[D], std::span<T>(upper));
        }
        an.curved.push_back(static_cast<int>(an.funcs.size()));
        an.funcs.push_back(f);
        an.ratios.push_back(ratio);
        an.uppers.push_back(upper);
        surface_vanishes |= f.ls == ctx.surface;
    }

    // Top level: the level sets' signs may decide the terms on the box, and
    // level sets the remaining terms do not name need no breakpoints.
    if (top)
    {
        an.decided = residual_terms(ctx.terms, an.signs, ctx.surface, surface_vanishes, an.terms);
        std::uint64_t needed = 0;
        if (an.decided == 0)
            for (const Requirement& r : an.terms)
                needed |= r.negative | r.positive;
        if (ctx.surface >= 0)
            needed |= std::uint64_t(1) << ctx.surface;
        if (an.decided == -1)
            needed = 0;
        std::vector<Func<T, D>> funcs;
        std::vector<VecD<T, D>> ratios, uppers;
        std::vector<int> curved;
        for (std::size_t i = 0, c = 0; i < an.funcs.size(); ++i)
        {
            const bool is_curved = c < an.curved.size() && an.curved[c] == static_cast<int>(i);
            if (!is_curved || ((needed >> an.funcs[i].ls) & 1))
            {
                if (is_curved)
                {
                    curved.push_back(static_cast<int>(funcs.size()));
                    ratios.push_back(an.ratios[c]);
                    uppers.push_back(an.uppers[c]);
                }
                funcs.push_back(an.funcs[i]);
            }
            c += is_curved;
        }
        an.funcs = std::move(funcs);
        an.ratios = std::move(ratios);
        an.uppers = std::move(uppers);
        an.curved = std::move(curved);
    }

    // A function constant along y_k has no root on its height lines and needs no
    // margin there; its restrictions to the bounds carry it to the base. Not the
    // level set of an interface, whose roots the top level integrates.
    const auto is_passive = [&](std::size_t c)
    { return top && ctx.surface >= 0 && an.funcs[static_cast<std::size_t>(an.curved[c])].ls != ctx.surface; };
    const auto margin_of = [&](std::size_t c, int kk)
    {
        if (is_passive(c))
            return T(1);
        const VecD<T, D>& upper = an.uppers[c];
        const T norm = scaled_norm(std::span<const T>(upper));
        const bool constant = norm > T(0) && norm < infinity<T> && upper[kk] <= tiny<T> * norm
                              && !(top && an.funcs[static_cast<std::size_t>(an.curved[c])].ls == ctx.surface);
        return constant ? T(1) : an.ratios[c][kk];
    };

    // surface functions: dropped where they have no root
    if constexpr (D <= 2)
    {
        for (const Surface<T, D>& s : p.surfaces)
        {
            VecD<T, D> ratio{};
            if (!surface_margins<T, D>(ctx, s, p.lo, p.hi, ratio))
                continue;
            an.surfaces.push_back(s);
            an.surface_ratios.push_back(ratio);
        }
    }

    // height direction: the best certified margin over all curved functions
    for (int kk = 0; kk < D; ++kk)
    {
        T best = T(1);
        for (std::size_t c = 0; c < an.ratios.size(); ++c)
            best = std::min(best, margin_of(c, kk));
        for (const VecD<T, D>& r : an.surface_ratios)
            best = std::min(best, r[kk]);
        if (best > an.best || (best == an.best && p.hi[kk] - p.lo[kk] > p.hi[an.k] - p.lo[an.k]))
        {
            an.best = best;
            an.k = kk;
        }
    }
    an.two_roots.assign(an.funcs.size(), 0);
    an.independent.assign(an.funcs.size(), 0);
    an.passive.assign(an.funcs.size(), 0);
    const auto mark_independent = [&]()
    {
        for (std::size_t c = 0; c < an.ratios.size(); ++c)
        {
            const std::size_t i = static_cast<std::size_t>(an.curved[c]);
            an.passive[i] = is_passive(c);
            an.independent[i] = !an.passive[i] && margin_of(c, an.k) == T(1) && an.ratios[c][an.k] < T(1);
        }
    };
    an.certified = (an.ratios.empty() && an.surface_ratios.empty()) || (an.best > T(0) && an.best >= ctx.margin);
    if (an.certified || depth < ctx.opt.two_roots_depth || an.ratios.empty())
    {
        mark_independent();
        return an;
    }

    // Two roots per line, where no direction gives one: bisection rarely
    // separates two sheets of one level set, so the roots themselves decide. A
    // direction qualifies where every curved function has a margin, or a
    // monotone derivative along it (second-derivative margin) with roots that
    // never merge in the box and are not steep (two_root_score). Each line
    // then splits at the extreme point.
    T best2 = -1;
    int k2 = -1;
    std::vector<int> chosen;
    for (int kk = 0; kk < D; ++kk)
    {
        T score = T(1);
        for (const VecD<T, D>& r : an.surface_ratios)
            score = std::min(score, r[kk]);
        if (!(score >= ctx.margin))
            continue;
        std::vector<int> two;
        for (std::size_t c = 0; c < an.ratios.size() && score >= ctx.margin; ++c)
        {
            if (margin_of(c, kk) >= ctx.margin)
            {
                score = std::min(score, margin_of(c, kk));
                continue;
            }
            const Func<T, D>& f = an.funcs[static_cast<std::size_t>(an.curved[c])];
            int cv = 0;
            if (!(second_margin<T, D>(ctx, f, p.lo, p.hi, kk, cv) >= ctx.margin))
            {
                score = T(-1);
                break;
            }
            score = std::min(score, two_root_score<T, D>(ctx, f, p.lo, p.hi, kk, cv));
            two.push_back(static_cast<int>(c));
        }
        if (score >= ctx.margin && score > best2)
        {
            best2 = score;
            k2 = kk;
            chosen = two;
        }
    }
    if (k2 < 0)
        return an;
    for (const int c : chosen)
        an.two_roots[static_cast<std::size_t>(an.curved[static_cast<std::size_t>(c)])] = 1;
    an.k = k2;
    an.best = best2;
    an.certified = true;
    mark_independent();
    ++ctx.stats->two_roots;
    return an;
}

/// Why did certification fail? Sample the clipped region on a 9^D grid and look
/// at the function with the smallest margin for the chosen axis near its zero set.
template <std::floating_point T, int D>
void diagnose_failure(Context<T>& ctx, const Problem<T, D>& p, const Analysis<T, D>& an, int depth)
{
    std::map<std::string, std::int64_t>& causes = ctx.stats->causes;
    const std::string level = "level " + std::to_string(D);
    const std::vector<VecD<T, D>>& ratios = an.ratios;
    if (ratios.empty())
    {
        ++causes[level + " | surface function | depth " + std::to_string(depth)];
        return;
    }
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
    // surface functions: the two columns of y among the inner coordinates
    for (const Surface<T, 2>& s : p.surfaces)
    {
        Surface<T, 2> g = s;
        const int j0 = s.k == 0 ? 1 : 0, j1 = s.k == 2 ? 1 : 2;
        for (Func<T, 3>* f : {&g.under, &g.on})
            for (int i = 0; i < 3; ++i)
            {
                const VecD<T, 2> c = transpose_times(VecD<T, 2>{f->A[i][j0], f->A[i][j1]});
                f->A[i][j0] = c[0];
                f->A[i][j1] = c[1];
            }
        r.surfaces.push_back(g);
    }
    return r;
}

/// A problem in a rotated frame y = R z, R orthogonal with the new axes as its
/// columns. The box becomes 2 D clips; the functions, the clips and the frame
/// map are composed with R.
template <typename T, int D>
Problem<T, D> rotate_frame(const Problem<T, D>& p, const std::array<VecD<T, D>, D>& R)
{
    // a row vector c of y, as a row vector of z: c R
    auto times = [&R](const VecD<T, D>& c)
    {
        VecD<T, D> r{};
        for (int m = 0; m < D; ++m)
            for (int j = 0; j < D; ++j)
                r[m] += c[j] * R[j][m];
        return r;
    };
    Problem<T, D> r;
    r.frame = p.frame;
    for (int i = 0; i < 3; ++i)
        r.frame.A[i] = times(p.frame.A[i]);
    for (const Half<T, D>& h : p.clips)
        r.clips.push_back({times(h.c), h.d});
    for (int j = 0; j < D; ++j)
    {
        VecD<T, D> e{};
        e[j] = T(1);
        r.clips.push_back({times(e), p.hi[j]});
        e[j] = T(-1);
        r.clips.push_back({times(e), -p.lo[j]});
    }
    // the bounding box of the box's corners, z = R^T y
    r.lo.fill(infinity<T>);
    r.hi.fill(-infinity<T>);
    for (int corner = 0; corner < (1 << D); ++corner)
        for (int m = 0; m < D; ++m)
        {
            T z = T(0);
            for (int j = 0; j < D; ++j)
                z += R[j][m] * (((corner >> j) & 1) ? p.hi[j] : p.lo[j]);
            r.lo[m] = std::min(r.lo[m], z);
            r.hi[m] = std::max(r.hi[m], z);
        }
    for (const Func<T, D>& f : p.funcs)
    {
        Func<T, D> g = f;
        if (f.linear)
            g.a = times(f.a);
        else
            for (int i = 0; i < 3; ++i)
                g.A[i] = times(f.A[i]);
        r.funcs.push_back(g);
    }
    if constexpr (D <= 2)
    {
        // surface functions: the columns of y among the inner coordinates
        for (const Surface<T, D>& s : p.surfaces)
        {
            Surface<T, D> g = s;
            for (Func<T, D + 1>* f : {&g.under, &g.on})
                for (int i = 0; i < 3; ++i)
                {
                    VecD<T, D> c;
                    for (int j = 0; j < D; ++j)
                        c[j] = f->A[i][j < s.k ? j : j + 1];
                    c = times(c);
                    for (int j = 0; j < D; ++j)
                        f->A[i][j < s.k ? j : j + 1] = c[j];
                }
            r.surfaces.push_back(g);
        }
    }
    return r;
}

/// A rotated frame for a box where no axis suits all curved functions, as
/// where the zero sets of two level sets meet: the height direction among the
/// unit normals of the functions (curved ones and surface functions) at the
/// box centre and the bisectors of pairs of them that makes the smallest
/// |n_i . d| largest. False if none reaches the margin there. The curved
/// functions are those of @p p that @p an lists.
template <std::floating_point T, int D>
bool choose_frame(Context<T>& ctx, const Problem<T, D>& p, const Analysis<T, D>& an,
                  std::array<VecD<T, D>, D>& R)
{
    VecD<T, D> centre;
    for (int j = 0; j < D; ++j)
        centre[j] = T(0.5) * (p.lo[j] + p.hi[j]);
    std::vector<VecD<T, D>> normals;
    auto add = [&normals](VecD<T, D> g)
    {
        const T norm = scaled_norm(std::span<const T>(g));
        if (!(norm > T(0)) || !std::isfinite(norm))
            return;
        for (int j = 0; j < D; ++j)
            g[j] /= norm;
        normals.push_back(g);
    };
    for (const int c : an.curved)
    {
        const Func<T, D>& f = p.funcs[static_cast<std::size_t>(c)];
        VecD<T, D> g;
        for (int j = 0; j < D; ++j)
            g[j] = derivative_along<T, D>(ctx, f, centre, j);
        add(g);
    }
    if constexpr (D <= 2)
    {
        for (const Surface<T, D>& s : p.surfaces)
        {
            VecD<T, D> g;
            surface_gradient<T, D>(ctx, s, centre, g);
            add(g);
        }
    }
    if (normals.size() < 2)
        return false;
    std::vector<VecD<T, D>> candidates = normals;
    for (std::size_t a = 0; a < normals.size(); ++a)
        for (std::size_t b = a + 1; b < normals.size(); ++b)
            for (const T sign : {T(1), T(-1)})
            {
                VecD<T, D> d;
                for (int j = 0; j < D; ++j)
                    d[j] = normals[a][j] + sign * normals[b][j];
                const T norm = scaled_norm(std::span<const T>(d));
                if (!(norm > T(0.1)))
                    continue;
                for (int j = 0; j < D; ++j)
                    d[j] /= norm;
                candidates.push_back(d);
            }
    T best = -1;
    VecD<T, D> dir{};
    for (const VecD<T, D>& d : candidates)
    {
        T score = T(1);
        for (const VecD<T, D>& n : normals)
        {
            T dot = T(0);
            for (int j = 0; j < D; ++j)
                dot += n[j] * d[j];
            score = std::min(score, std::abs(dot));
        }
        if (score > best)
        {
            best = score;
            dir = d;
        }
    }
    if (!(best >= ctx.margin))
        return false;
    // an orthonormal frame with dir as its first axis (columns of R)
    for (int j = 0; j < D; ++j)
        R[j][0] = dir[j];
    if constexpr (D == 2)
    {
        R[0][1] = -dir[1];
        R[1][1] = dir[0];
    }
    else if constexpr (D == 3)
    {
        int e = 0;
        for (int j = 1; j < 3; ++j)
            if (std::abs(dir[j]) < std::abs(dir[e]))
                e = j;
        VecD<T, 3> a{};
        a[e] = T(1);
        T dot = dir[e];
        T norm = T(0);
        for (int j = 0; j < 3; ++j)
        {
            a[j] -= dot * dir[j];
            norm += a[j] * a[j];
        }
        norm = std::sqrt(norm);
        for (int j = 0; j < 3; ++j)
            R[j][1] = a[j] / norm;
        const VecD<T, 3> b = {dir[1] * R[2][1] - dir[2] * R[1][1], dir[2] * R[0][1] - dir[0] * R[2][1],
                              dir[0] * R[1][1] - dir[1] * R[0][1]};
        for (int j = 0; j < 3; ++j)
            R[j][2] = b[j];
    }
    else
        return false;
    return true;
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

/// Roots in (L, U) of a curved function on the height line y = insert(yb, k, t),
/// appended to @p roots: on a certified box the one root the ends bracket, or
/// with two roots per line the roots on either side of the extreme point;
/// otherwise all roots, isolated. Bernstein forms leave the line's form in
/// @p form.
template <std::floating_point T, int D>
void curved_roots(Context<T>& ctx, const Func<T, D>& f, const VecD<T, D - 1>& yb, int k, T L, T U, bool certified,
                  bool two, BoxBernstein<T>& form, std::vector<T>& roots)
{
    const Source<T>& phi = ctx.phis[static_cast<std::size_t>(f.ls)];
    const VecD<T, D> yL = insert<T, D>(yb, k, L);
    if (phi.is_analytic())
    {
        Vec3<T> u0, dir;
        box_line<T, D>(f, insert<T, D>(yb, k, T(0)), k, u0, dir);
        if (!certified)
        {
            line_roots(phi, std::span<const T>(u0), std::span<const T>(dir), L, U,
                       zero_function_tol<T> * ctx.scales[static_cast<std::size_t>(f.ls)], roots);
            return;
        }
        const T gl = value<T, D>(ctx, f, yL), gu = value<T, D>(ctx, f, insert<T, D>(yb, k, U));
        T c = U, gc = gu;
        if (two)
        {
            // the extreme point: the root of d_k psi, which is monotone
            auto slope = [&](T t) { return derivative_along<T, D>(ctx, f, insert<T, D>(yb, k, t), k); };
            const T dl = slope(L), du = slope(U);
            if (opposite(dl, du))
            {
                c = illinois_root(slope, L, U, dl, du);
                gc = value<T, D>(ctx, f, insert<T, D>(yb, k, c));
                if (opposite(gl, gc))
                    roots.push_back(line_root(phi, std::span<const T>(u0), std::span<const T>(dir), L, c, gl, gc));
                if (opposite(gc, gu))
                    roots.push_back(line_root(phi, std::span<const T>(u0), std::span<const T>(dir), c, U, gc, gu));
                return;
            }
        }
        if (opposite(gl, gu))
            roots.push_back(line_root(phi, std::span<const T>(u0), std::span<const T>(dir), L, U, gl, gu));
        return;
    }
    line_form<T, D>(ctx, f, yL, k, L, U, form, ctx.line_work[D]);
    const std::span<const T> c(form.coeffs);
    if (!certified)
    {
        isolate_roots(c, L, U, roots, ctx.root_work[D]);
        return;
    }
    const T gl = c.front(), gu = c.back();
    if (two && c.size() > 2)
    {
        // the extreme point: the root of the line's derivative, which is monotone
        std::vector<T>& dc = ctx.slope_line[D];
        const int n = static_cast<int>(c.size()) - 1;
        dc.resize(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i)
            dc[static_cast<std::size_t>(i)] = T(n) * (c[static_cast<std::size_t>(i + 1)] - c[static_cast<std::size_t>(i)]);
        const T dl = dc.front(), du = dc.back();
        if (opposite(dl, du))
        {
            const T x = bracketed_root(std::span<const T>(dc), L, U, L, U, dl, du);
            const T gx = evaluate_1d(c, (x - L) / (U - L));
            if (opposite(gl, gx))
                roots.push_back(bracketed_root(c, L, U, L, x, gl, gx));
            if (opposite(gx, gu))
                roots.push_back(bracketed_root(c, L, U, x, U, gx, gu));
            return;
        }
    }
    if (opposite(gl, gu))
        roots.push_back(bracketed_root(c, L, U, L, U, gl, gu));
}

/// Does a top-level point satisfy one of the terms left on its box? Values of
/// the level sets come from their forms on the current line where they exist.
template <std::floating_point T>
bool holds(Context<T>& ctx, const std::vector<Requirement>& terms, const Vec3<T>& u, T s)
{
    std::uint64_t evaluated = 0;
    for (const Requirement& r : terms)
    {
        bool ok = true;
        for (std::uint64_t bits = r.negative | r.positive; bits != 0 && ok; bits &= bits - 1)
        {
            const int l = std::countr_zero(bits);
            const std::uint64_t bit = std::uint64_t(1) << l;
            T& v = ctx.top_values[static_cast<std::size_t>(l)];
            if (!(evaluated & bit))
            {
                evaluated |= bit;
                if (ctx.top_line_set[static_cast<std::size_t>(l)])
                    v = evaluate_1d(std::span<const T>(ctx.top_lines[static_cast<std::size_t>(l)].coeffs), s);
                else
                    v = evaluate(ctx.phis[static_cast<std::size_t>(l)], std::span<const T>(u));
            }
            ok = (r.negative & bit) ? v < T(0) : v > T(0);
        }
        if (ok)
            return true;
    }
    return false;
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
        if (is_analytic(ctx, f.ls))
        {
            Vec3<T> u0, dir;
            box_line<T, 1>(f, VecD<T, 1>{T(0)}, 0, u0, dir);
            line_roots(ctx.phis[static_cast<std::size_t>(f.ls)], std::span<const T>(u0), std::span<const T>(dir), L,
                       U, zero_function_tol<T> * ctx.scales[static_cast<std::size_t>(f.ls)], nodes);
            continue;
        }
        line_form<T, 1>(ctx, f, VecD<T, 1>{L}, 0, L, U, ctx.line[1], ctx.line_work[1]);
        isolate_roots(std::span<const T>(ctx.line[1].coeffs), L, U, nodes, ctx.root_work[1]);
    }
    for (const Surface<T, 1>& s : p.surfaces)
        surface_line_roots<T>(ctx, s, L, U, nodes);
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

/// Where the roots of two surface functions on the height lines of a level-2
/// box cross: three level sets meet there (a corner), and the integrand of the
/// line below has a kink. Both are monotone along k on the certified box; the
/// difference of their roots is sampled along the base and its sign changes
/// are refined. The base coordinates of the crossings are appended.
template <std::floating_point T>
void surface_crossings(Context<T>& ctx, const Surface<T, 2>& a, const Surface<T, 2>& b, int k, T lo, T hi, T blo,
                       T bhi, std::vector<T>& out)
{
    auto root = [&](const Surface<T, 2>& s, T y) -> T
    {
        const T fl = surface_value<T, 2>(ctx, s, insert<T, 2>(VecD<T, 1>{y}, k, lo)),
                fh = surface_value<T, 2>(ctx, s, insert<T, 2>(VecD<T, 1>{y}, k, hi));
        if (!opposite(fl, fh))
            return std::numeric_limits<T>::quiet_NaN();
        return illinois_root([&](T t) { return surface_value<T, 2>(ctx, s, insert<T, 2>(VecD<T, 1>{y}, k, t)); }, lo, hi,
                             fl, fh);
    };
    auto gap = [&](T y) { return root(a, y) - root(b, y); };
    constexpr int samples = 9;
    T y0 = blo, g0 = gap(blo);
    for (int i = 1; i < samples; ++i)
    {
        const T y1 = blo + (bhi - blo) * T(i) / T(samples - 1), g1 = gap(y1);
        if (std::isfinite(g0) && std::isfinite(g1) && opposite(g0, g1))
        {
            // the roots may vanish inside: bisect while both exist
            T l = y0, h = y1, gl = g0;
            for (int it = 0; it < 60 && h - l > segment_tol<T> * (bhi - blo); ++it)
            {
                const T m = T(0.5) * (l + h), gm = gap(m);
                if (!std::isfinite(gm))
                    break;
                if (opposite(gl, gm))
                    h = m;
                else
                {
                    l = m;
                    gl = gm;
                }
            }
            out.push_back(T(0.5) * (l + h));
        }
        y0 = y1;
        g0 = g1;
    }
}

template <std::floating_point T, int D, bool Rotated, int Top, typename Emit>
void integrate_box(Context<T>& ctx, Problem<T, D> p, const Emit& emit, int depth, bool interface);

/// True if the interface's level set vanishes on the whole plane y_axis = mid
/// of the box: halves split there would share its zero set as a face, on which
/// neither has a root.
template <std::floating_point T, int D>
bool interface_on_plane(Context<T>& ctx, const Problem<T, D>& p, int axis, T mid)
{
    if constexpr (D < 2)
        return false;
    else
    {
        thread_local BoxBernstein<T> form;
        thread_local std::vector<T> work;
        for (const Func<T, D>& f : p.funcs)
        {
            if (f.linear || f.ls != ctx.surface)
                continue;
            Bound<T, D> plane;
            plane.alpha = mid;
            const Func<T, D - 1> g = restrict_to<T, D>(f, axis, plane);
            VecD<T, D - 1> lo, hi;
            for (int jb = 0; jb < D - 1; ++jb)
            {
                const int j = jb < axis ? jb : jb + 1;
                lo[jb] = p.lo[j];
                hi[jb] = p.hi[j];
            }
            const T tol = zero_function_tol<T> * ctx.scales[static_cast<std::size_t>(f.ls)];
            if (is_analytic(ctx, f.ls))
            {
                int sign = 0;
                VecD<T, D - 1> ratio{}, upper{};
                T magnitude = T(0);
                analytic_bounds<T, D - 1>(ctx, g, lo, hi, sign, ratio, magnitude, upper);
                return magnitude <= tol;
            }
            bernstein_form<T, D - 1>(ctx, g, lo, hi, form, work);
            return max_abs(std::span<const T>(form.coeffs)) <= tol;
        }
        return false;
    }
}

/// The dimension reduction of a tightened and analysed box: bisect if it is not
/// certified, otherwise integrate along height lines over the base level.
template <std::floating_point T, int D, bool Rotated, int Top, typename Emit>
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
            T mid = T(0.5) * (p.lo[axis] + p.hi[axis]);
            if constexpr (D == Top)
                if (interface && interface_on_plane<T, D>(ctx, p, axis, mid))
                    mid = p.lo[axis] + T(0.625) * (p.hi[axis] - p.lo[axis]);
            Problem<T, D> first = p, second = std::move(p);
            first.hi[axis] = mid;
            second.lo[axis] = mid;
            integrate_box<T, D, Rotated, Top>(ctx, std::move(first), emit, depth + 1, interface);
            integrate_box<T, D, Rotated, Top>(ctx, std::move(second), emit, depth + 1, interface);
            return;
        }
        ++ctx.stats->uncertified;
    }

    // bounds of the height lines: box faces and clip planes
    std::vector<Bound<T, D>> lowers, uppers;
    Problem<T, D - 1> base_box;
    for (int jb = 0; jb < D - 1; ++jb)
    {
        const int j = jb < k ? jb : jb + 1;
        base_box.lo[jb] = p.lo[j];
        base_box.hi[jb] = p.hi[j];
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
            base_box.clips.push_back(hb);
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
            base_box.clips.push_back(hb);
        }

    // Drop bounds that are never active on the base region: a lower bound below
    // another lower bound at every vertex of the base polytope (bounds are affine),
    // or an upper bound above another upper bound. Their restrictions and switch
    // functions would only cause needless bisections.
    if (ctx.opt.prune_bounds && lowers.size() + uppers.size() > 2)
    {
        std::vector<VecD<T, D - 1>> corners;
        polytope_vertices<T, D - 1>(base_box.lo, base_box.hi, base_box.clips, corners);
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

    // Duplicates (a clip in a box face, say) would be counted twice below.
    auto dedupe = [](std::vector<Bound<T, D>>& list)
    {
        std::vector<Bound<T, D>> kept;
        for (const Bound<T, D>& b : list)
        {
            bool same = false;
            for (const Bound<T, D>& c : kept)
            {
                bool equal = std::abs(b.alpha - c.alpha) <= compare_tol<T>;
                for (int jb = 0; jb < D - 1; ++jb)
                    equal &= std::abs(b.beta[jb] - c.beta[jb]) <= compare_tol<T>;
                same |= equal;
            }
            if (!same)
                kept.push_back(b);
        }
        list = std::move(kept);
    };
    dedupe(lowers);
    dedupe(uppers);

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

    // corners: the roots of two surface functions crossing (three level sets)
    std::vector<T> corners;
    if constexpr (D == 2)
    {
        if (certified)
            for (std::size_t i = 0; i < p.surfaces.size(); ++i)
                for (std::size_t j = i + 1; j < p.surfaces.size(); ++j)
                    surface_crossings<T>(ctx, p.surfaces[i], p.surfaces[j], k, p.lo[k], p.hi[k], base_box.lo[0],
                                         base_box.hi[0], corners);
    }

    // The base of the height lines between the given lower and upper bounds, on
    // the base polytope cut by the given clips: its functions are where the
    // integrand along the lines changes form.
    const auto integrate_base = [&](const std::vector<Bound<T, D>>& lows, const std::vector<Bound<T, D>>& ups,
                                    const std::vector<Half<T, D - 1>>& clips)
    {
        Problem<T, D - 1> base;
        base.lo = base_box.lo;
        base.hi = base_box.hi;
        base.clips = base_box.clips;
        base.clips.insert(base.clips.end(), clips.begin(), clips.end());
        std::vector<Bound<T, D>> bounds = lows;
        bounds.insert(bounds.end(), ups.begin(), ups.end());
        // A passive function has no roots that matter here, but its restrictions
        // to the bounds pair with the interface's: where the interface's zero set
        // leaves through a face and crosses the other zero set there, the
        // integrand below has a kink, which only their surface function marks.
        for (std::size_t i = 0; i < p.funcs.size(); ++i)
        {
            const Func<T, D>& f = p.funcs[i];
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
        for (std::size_t i = 0; i < lows.size(); ++i) // the active lower bound changes
            for (std::size_t j = i + 1; j < lows.size(); ++j)
                base.funcs.push_back(difference<T, D>(lows[i], lows[j]));
        for (std::size_t i = 0; i < ups.size(); ++i)
            for (std::size_t j = i + 1; j < ups.size(); ++j)
                base.funcs.push_back(difference<T, D>(ups[i], ups[j]));
        if constexpr (D <= 2)
        {
            for (const Surface<T, D>& s : p.surfaces)
            {
                for (const Bound<T, D>& b : bounds)
                    base.surfaces.push_back(restrict_surface<T, D>(s, k, b));
                for (const Bound<T, D>& r : linear_roots)
                    base.surfaces.push_back(restrict_surface<T, D>(s, k, r));
            }
        }
        if constexpr (D == 2)
        {
            for (const T y : corners)
            {
                Func<T, 1> f;
                f.origin = 12;
                f.linear = true;
                f.a[0] = T(1);
                f.c = -y;
                base.funcs.push_back(f);
            }
        }
        // Where two level sets cut a certified box, their roots on the height lines
        // cross on a ridge: the base splits there, where one level set vanishes on
        // the other's zero set. The surface runs along one with a single root per line.
        if (certified)
        {
            for (std::size_t i = 0; i < p.funcs.size(); ++i)
                for (std::size_t j = 0; j < p.funcs.size(); ++j)
                {
                    const Func<T, D>& f = p.funcs[i];
                    const Func<T, D>& g = p.funcs[j];
                    // a function constant along y_k: no roots there to cross; a
                    // passive one is only crossed
                    if (f.linear || g.linear || f.ls == g.ls || an.two_roots[i] || an.independent[i]
                        || an.independent[j] || an.passive[i])
                        continue;
                    // one surface per pair: along the first with one root per line
                    if (j < i && !an.two_roots[j] && !an.passive[j])
                        continue;
                    Surface<T, D - 1> s;
                    s.under = f;
                    s.on = g;
                    s.k = k;
                    s.t_lo = p.lo[k];
                    s.t_hi = p.hi[k];
                    base.surfaces.push_back(s);
                    ++ctx.stats->surfaces;
                }
        }

        // the rule along the height lines: the integrand of the base level
        const int box_id = ctx.box_ids[D - 1]++;
        const auto line = [&, k, certified, box_id](const VecD<T, D - 1>& yb, T w, const NodeTag& below)
        {
            T L = -infinity<T>, U = infinity<T>;
            for (const Bound<T, D>& b : lows)
                L = std::max(L, bound_value<T, D>(b, yb));
            for (const Bound<T, D>& b : ups)
                U = std::min(U, bound_value<T, D>(b, yb));
            if (!(U > L))
                return;
            std::vector<T>& nodes = ctx.nodes[D];
            nodes.clear();
            if constexpr (D == Top)
                std::fill(ctx.top_line_set.begin(), ctx.top_line_set.end(), char(0));
            for (std::size_t i = 0; i < p.funcs.size(); ++i)
            {
                const Func<T, D>& f = p.funcs[i];
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
                BoxBernstein<T>* form = &ctx.line[D];
                if constexpr (D == Top)
                {
                    // the interface: the roots of its level set; the others are
                    // evaluated at them
                    if (interface && f.ls != ctx.surface)
                        continue;
                    if (!is_analytic(ctx, f.ls))
                    {
                        form = &ctx.top_lines[static_cast<std::size_t>(f.ls)];
                        ctx.top_line_set[static_cast<std::size_t>(f.ls)] = 1;
                    }
                }
                curved_roots<T, D>(ctx, f, yb, k, L, U, certified, an.two_roots[i] != 0, *form, nodes);
            }
            if constexpr (D <= 2)
            {
                for (const Surface<T, D>& s : p.surfaces)
                {
                    if (certified)
                    {
                        const T sl = surface_value<T, D>(ctx, s, insert<T, D>(yb, k, L)),
                                su = surface_value<T, D>(ctx, s, insert<T, D>(yb, k, U));
                        if (opposite(sl, su))
                            nodes.push_back(illinois_root(
                                [&](T t) { return surface_value<T, D>(ctx, s, insert<T, D>(yb, k, t)); }, L, U, sl, su));
                    }
                    else
                        surface_line_roots<T>(ctx, surface_on_line<T, D>(s, yb, k), L, U, nodes);
                }
            }

            NodeTag tag = below;
            tag.box[D - 1] = box_id;
            if constexpr (D == Top)
            {
                if (interface)
                {
                    // the roots of the level set on each line, weighted by the physical
                    // surface measure: |det J| |J^-T grad phi| / |d_k phi|
                    std::sort(nodes.begin(), nodes.end());
                    const Source<T>& phi = ctx.phis[static_cast<std::size_t>(ctx.surface)];
                    for (std::size_t r = 0; r < nodes.size(); ++r)
                    {
                        tag.segment[D - 1] = static_cast<int>(r);
                        const VecD<T, D> y = insert<T, D>(yb, k, nodes[r]);
                        const Vec3<T> u = to_u<T, D>(p.frame, y);
                        if (ctx.vis_order == 0 && an.decided == 0 && !holds<T>(ctx, an.terms, u, T(0)))
                            continue;
                        Vec3<T> g = {0, 0, 0};
                        gradient(phi, std::span<const T>(u), std::span<T>(g));
                        // d_k phi in the level's (possibly rotated) coordinates
                        T dk = T(0);
                        for (int i = 0; i < 3; ++i)
                            dk += p.frame.A[i][k] * g[i];
                        dk = std::abs(dk);
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
                                  if constexpr (D == Top)
                                  {
                                      if (ctx.vis_order == 0)
                                      {
                                          // top level: keep the points of the selected part
                                          if (an.decided == 0
                                              && !holds<T>(ctx, an.terms, to_u<T, D>(p.frame, y), (t - L) / (U - L)))
                                              return;
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
            integrate_box<T, D - 1, false, Top>(ctx, std::move(base), line, 0, false);
    };

    // With several active bounds, the base splits into the regions where one
    // lower and one upper bound are active: each needs only the restrictions to
    // those two, which are fewer functions to certify together. Where several
    // level sets meet, or in a rotated frame (whose box faces are clips), that
    // keeps the levels below tractable; with one level set the joint base
    // bisects more and is a little more accurate.
    bool several = !p.surfaces.empty();
    for (std::size_t i = 1; i < p.funcs.size() && !several; ++i)
        several = !p.funcs[i].linear && !p.funcs[0].linear && p.funcs[i].ls != p.funcs[0].ls;
    for (std::size_t i = 0; i < p.funcs.size() && !several; ++i)
        for (std::size_t j = i + 1; j < p.funcs.size() && !several; ++j)
            several = !p.funcs[i].linear && !p.funcs[j].linear && p.funcs[i].ls != p.funcs[j].ls;
    if (!ctx.opt.split_bounds || !(several || Rotated) || lowers.size() * uppers.size() == 1)
    {
        integrate_base(lowers, uppers, {});
        return;
    }
    std::vector<Half<T, D - 1>> clips;
    for (std::size_t i = 0; i < lowers.size(); ++i)
        for (std::size_t j = 0; j < uppers.size(); ++j)
        {
            clips.clear();
            // lower_m <= lower_i and upper_j <= upper_m
            for (std::size_t m = 0; m < lowers.size(); ++m)
                if (m != i)
                {
                    Half<T, D - 1> h;
                    for (int jb = 0; jb < D - 1; ++jb)
                        h.c[jb] = lowers[m].beta[jb] - lowers[i].beta[jb];
                    h.d = lowers[i].alpha - lowers[m].alpha;
                    clips.push_back(h);
                }
            for (std::size_t m = 0; m < uppers.size(); ++m)
                if (m != j)
                {
                    Half<T, D - 1> h;
                    for (int jb = 0; jb < D - 1; ++jb)
                        h.c[jb] = uppers[j].beta[jb] - uppers[m].beta[jb];
                    h.d = uppers[m].alpha - uppers[j].alpha;
                    clips.push_back(h);
                }
            integrate_base({lowers[i]}, {uppers[j]}, clips);
        }
}

template <std::floating_point T, int D, bool Rotated, int Top, typename Emit>
void integrate_box(Context<T>& ctx, Problem<T, D> p, const Emit& emit, int depth, bool interface)
{
    if (!tighten<T, D>(p.lo, p.hi, p.clips))
        return;

    // keep the functions that may vanish on the box, with their direction margins
    Analysis<T, D> an = analyse<T, D>(ctx, p, D == Top, depth);
    if constexpr (D == Top)
    {
        if (an.decided < 0)
            return; // the selection holds nowhere in the box
        if (an.decided > 0 && !interface)
        {
            // it holds on the whole box: no breakpoints needed
            an.funcs.clear();
            an.ratios.clear();
            an.uppers.clear();
            an.curved.clear();
            an.two_roots.clear();
            an.independent.clear();
            an.passive.clear();
            an.certified = true;
        }
    }
    else
    {
        if (interface && an.funcs.empty())
            return;
    }
    p.funcs = std::move(an.funcs);
    p.surfaces = std::move(an.surfaces);

    // Level 2: two functions whose zero curves meet (on a face of a tet, say) may
    // each need a different axis, and no bisection separates them. A diagonal
    // frame often suits both. Not at the top level of 2D cells, whose points
    // are the box's own coordinates.
    if constexpr (D == 2 && !Rotated && Top != 2)
    {
        if (!an.certified && ctx.opt.diagonal_frames)
        {
            Problem<T, 2> r = rotate_diagonal<T>(p);
            if (tighten<T, 2>(r.lo, r.hi, r.clips))
            {
                Analysis<T, 2> ar = analyse<T, 2>(ctx, r, false, depth);
                if (ar.certified)
                {
                    ++ctx.stats->rotations;
                    r.funcs = std::move(ar.funcs);
                    r.surfaces = std::move(ar.surfaces);
                    const auto back = [&emit](const VecD<T, 2>& z, T w, const NodeTag& tag)
                    { emit(VecD<T, 2>{T(0.5) * (z[0] - z[1]), T(0.5) * (z[0] + z[1])}, T(0.5) * w, tag); };
                    reduce<T, 2, true, Top>(ctx, std::move(r), ar, back, depth, interface);
                    return;
                }
            }
        }
    }
    // Where the zero sets of two level sets meet, their normals may differ so
    // much that no axis suits both: a height direction between them does.
    if constexpr (!Rotated)
    {
        // two level sets: curved functions of both, or surface functions
        bool two = !p.surfaces.empty();
        for (std::size_t c = 1; c < an.curved.size() && !two; ++c)
            two = p.funcs[static_cast<std::size_t>(an.curved[c])].ls
                  != p.funcs[static_cast<std::size_t>(an.curved[0])].ls;
        std::array<VecD<T, D>, D> R;
        // in 3D only once bisection has not helped: a rotated box has its faces
        // as clip planes, which make the levels below costly
        const bool allowed = D < 3 || depth >= ctx.opt.rotation_depth;
        if (!an.certified && two && allowed && choose_frame<T, D>(ctx, p, an, R))
        {
            Problem<T, D> r = rotate_frame<T, D>(p, R);
            if (tighten<T, D>(r.lo, r.hi, r.clips))
            {
                Analysis<T, D> ar = analyse<T, D>(ctx, r, D == Top, depth);
                if constexpr (D == Top)
                {
                    // the rotated box fits the region differently: its bounds may decide it
                    if (ar.decided < 0)
                        return;
                    if (ar.decided > 0 && !interface)
                    {
                        ar.funcs.clear();
                        ar.ratios.clear();
                        ar.uppers.clear();
                        ar.curved.clear();
                        ar.two_roots.clear();
                        ar.independent.clear();
                        ar.passive.clear();
                        ar.certified = true;
                    }
                }
                if (ar.certified)
                {
                    ++ctx.stats->rotations;
                    r.funcs = std::move(ar.funcs);
                    r.surfaces = std::move(ar.surfaces);
                    const auto back = [&emit, R](const VecD<T, D>& z, T w, const NodeTag& tag)
                    {
                        VecD<T, D> y{};
                        for (int j = 0; j < D; ++j)
                            for (int m = 0; m < D; ++m)
                                y[j] += R[j][m] * z[m];
                        emit(y, w, tag);
                    };
                    reduce<T, D, true, Top>(ctx, std::move(r), ar, back, depth, interface);
                    return;
                }
            }
        }
    }
    reduce<T, D, Rotated, Top>(ctx, std::move(p), an, emit, depth, interface);
}

/// The top level: the unit box with the cell's clips and one function per
/// level set the terms name.
template <std::floating_point T, int Top>
void run_top(Context<T>& ctx, const ClippedBox<T>& cell, CellPoints<T>& out, bool interface)
{
    Problem<T, Top> top;
    top.lo.fill(T(0));
    top.hi.fill(T(1));
    for (int i = 0; i < Top; ++i)
        top.frame.A[i][i] = T(1);
    for (int l = 0; l < static_cast<int>(ctx.phis.size()); ++l)
    {
        if (!((ctx.used >> l) & 1))
            continue;
        Func<T, Top> f;
        f.ls = l;
        for (int i = 0; i < Top; ++i)
            f.A[i][i] = T(1);
        top.funcs.push_back(f);
    }
    for (const HalfSpace<T>& h : cell.clips)
    {
        Half<T, Top> c;
        for (int j = 0; j < Top; ++j)
            c.c[j] = h.c[j];
        c.d = h.d;
        top.clips.push_back(c);
    }
    const int vis_order = ctx.vis_order;
    const auto emit = [&out, vis_order](const VecD<T, Top>& y, T w, const NodeTag& tag)
    {
        const Vec3<T> u = pad<T, Top>(y);
        out.points.insert(out.points.end(), u.begin(), u.end());
        out.weights.push_back(w);
        if (vis_order > 0)
            out.tags.push_back(tag);
    };
    integrate_box<T, Top, false, Top>(ctx, std::move(top), emit, 0, interface);
}

/// The engine on one cell.
template <std::floating_point T>
void run(const ClippedBox<T>& cell, std::span<const Source<T>> phis, std::span<const SelectionTerm> terms, int q,
         int vis_order, const Options& opt, CellPoints<T>& out, Stats& stats)
{
    if (phis.empty() || phis.size() > 64)
        throw std::invalid_argument("quadrays: give 1 to 64 level sets");
    if (cell.tdim != 2 && cell.tdim != 3)
        throw std::invalid_argument("quadrays: cells of dimension 2 or 3");
    for (const Source<T>& phi : phis)
    {
        if ((phi.bernstein == nullptr) == (phi.analytic == nullptr))
            throw std::invalid_argument("quadrays: a source is a Bernstein form or an analytic level set");
        if (phi.bernstein != nullptr && phi.bernstein->dim != cell.tdim)
            throw std::invalid_argument("quadrays: the level set must be a form in as many variables as the cell has");
        if (phi.analytic != nullptr && phi.tdim != cell.tdim)
            throw std::invalid_argument("quadrays: the analytic level set belongs to a cell of another dimension");
    }
    if (vis_order == 0 && q < 1)
        throw std::invalid_argument("quadrays: at least one Gauss point per segment is needed");

    thread_local Context<T> ctx;
    const std::size_t n = phis.size();
    const std::uint64_t all = n == 64 ? ~std::uint64_t(0) : (std::uint64_t(1) << n) - 1;
    ctx.terms.clear();
    ctx.used = 0;
    std::uint64_t zero = 0;
    bool first = true;
    for (const SelectionTerm& t : terms)
    {
        const Requirement r = {t.negative_required, t.positive_required, t.zero_required};
        if ((r.negative | r.positive | r.zero) & ~all)
            throw std::invalid_argument("quadrays: a term names a level set the cell does not have");
        if (std::popcount(r.zero) > 1)
            throw std::invalid_argument("quadrays: curves where two level sets vanish are not integrated");
        if (!first && r.zero != zero)
            throw std::invalid_argument("quadrays: the terms ask for different zero sets");
        zero = r.zero;
        first = false;
        // a term that asks for two signs of one level set selects nothing
        if ((r.negative & r.positive) | (r.zero & (r.negative | r.positive)))
            continue;
        ctx.terms.push_back(r);
        ctx.used |= r.negative | r.positive | r.zero;
    }
    if (ctx.terms.empty())
        return;
    ctx.surface = zero != 0 ? std::countr_zero(zero) : -1;
    ctx.phis.assign(phis.begin(), phis.end());
    // An analytic level set that is linear on the cell (a plane) is read as its
    // form of degree 1, as a level set of degree 1 is: its zero set then needs
    // no surface function where it meets another's, whose bounds would come
    // from Taylor models.
    ctx.linear_forms.resize(n);
    for (std::size_t l = 0; l < n; ++l)
        if (((ctx.used >> l) & 1) && phis[l].is_analytic() && linear_form(phis[l], ctx.linear_forms[l]))
            ctx.phis[l] = bernstein_source(ctx.linear_forms[l]);
    ctx.scales.resize(n);
    for (std::size_t l = 0; l < n; ++l)
        ctx.scales[l] = ((ctx.used >> l) & 1) ? reference_magnitude(ctx.phis[l]) : T(0);
    ctx.top_lines.resize(n);
    ctx.top_line_set.assign(n, 0);
    ctx.top_values.assign(n, T(0));
    ctx.opt = opt;
    ctx.margin = static_cast<T>(opt.margin);
    ctx.stats = &stats;
    ctx.vis_order = vis_order;
    const bool interface = ctx.surface >= 0;
    ctx.detj = std::abs(jacobian_determinant(cell));
    ctx.inv = interface ? inverse_jacobian(cell) : Mat3<T>{};
    ctx.box_ids = {0, 0, 0};
    ctx.bisections = 0;
    if (vis_order == 0 && ctx.q != q)
    {
        gauss_legendre<T>(q, ctx.gauss_x, ctx.gauss_w);
        ctx.q = q;
    }
    if (cell.tdim == 3)
        run_top<T, 3>(ctx, cell, out, interface);
    else
        run_top<T, 2>(ctx, cell, out, interface);
}

} // namespace

template <std::floating_point T>
void integrate(const ClippedBox<T>& cell, std::span<const Source<T>> phis, std::span<const SelectionTerm> terms,
               int q, const Options& opt, CellPoints<T>& out, Stats& stats)
{
    run(cell, phis, terms, q, 0, opt, out, stats);
}

template <std::floating_point T>
void integrate(const ClippedBox<T>& cell, const Source<T>& phi, Part part, int q, const Options& opt,
               CellPoints<T>& out, Stats& stats)
{
    const SelectionTerm term = term_of(part);
    run(cell, std::span<const Source<T>>(&phi, 1), std::span<const SelectionTerm>(&term, 1), q, 0, opt, out, stats);
}

template <std::floating_point T>
void integrate(const ClippedBox<T>& cell, const BoxBernstein<T>& phi, Part part, int q,
               const Options& opt, CellPoints<T>& out, Stats& stats)
{
    integrate(cell, bernstein_source(phi), part, q, opt, out, stats);
}

template <std::floating_point T>
void leaf_nodes(const ClippedBox<T>& cell, std::span<const Source<T>> phis, std::span<const SelectionTerm> terms,
                int degree, const Options& opt, CellPoints<T>& out, Stats& stats)
{
    if (degree < 1)
        throw std::invalid_argument("quadrays: leaf cells need degree 1 or more");
    run(cell, phis, terms, 0, degree, opt, out, stats);
}

template <std::floating_point T>
void leaf_nodes(const ClippedBox<T>& cell, const Source<T>& phi, Part part, int degree, const Options& opt,
                CellPoints<T>& out, Stats& stats)
{
    const SelectionTerm term = term_of(part);
    leaf_nodes(cell, std::span<const Source<T>>(&phi, 1), std::span<const SelectionTerm>(&term, 1), degree, opt, out,
               stats);
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

template void integrate<float>(const ClippedBox<float>&, std::span<const Source<float>>,
                               std::span<const SelectionTerm>, int, const Options&, CellPoints<float>&, Stats&);
template void integrate<double>(const ClippedBox<double>&, std::span<const Source<double>>,
                                std::span<const SelectionTerm>, int, const Options&, CellPoints<double>&, Stats&);
template void integrate<float>(const ClippedBox<float>&, const Source<float>&, Part, int, const Options&,
                               CellPoints<float>&, Stats&);
template void integrate<double>(const ClippedBox<double>&, const Source<double>&, Part, int, const Options&,
                                CellPoints<double>&, Stats&);
template void integrate<float>(const ClippedBox<float>&, const BoxBernstein<float>&, Part, int,
                               const Options&, CellPoints<float>&, Stats&);
template void integrate<double>(const ClippedBox<double>&, const BoxBernstein<double>&, Part, int,
                                const Options&, CellPoints<double>&, Stats&);
template void leaf_nodes<float>(const ClippedBox<float>&, std::span<const Source<float>>,
                                std::span<const SelectionTerm>, int, const Options&, CellPoints<float>&, Stats&);
template void leaf_nodes<double>(const ClippedBox<double>&, std::span<const Source<double>>,
                                 std::span<const SelectionTerm>, int, const Options&, CellPoints<double>&, Stats&);
template void leaf_nodes<float>(const ClippedBox<float>&, const Source<float>&, Part, int, const Options&,
                                CellPoints<float>&, Stats&);
template void leaf_nodes<double>(const ClippedBox<double>&, const Source<double>&, Part, int, const Options&,
                                 CellPoints<double>&, Stats&);
template void leaf_nodes<float>(const ClippedBox<float>&, const BoxBernstein<float>&, Part, int,
                                const Options&, CellPoints<float>&, Stats&);
template void leaf_nodes<double>(const ClippedBox<double>&, const BoxBernstein<double>&, Part, int,
                                 const Options&, CellPoints<double>&, Stats&);

} // namespace cutcells::quadrays
