// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "curves.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "box_bernstein.h"
#include "gauss.h"

namespace cutcells::quadrays
{

namespace
{

// ============================================================================
// Small vectors
// ============================================================================

template <std::floating_point T>
T dot(const Vec3<T>& a, const Vec3<T>& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <std::floating_point T>
Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

template <std::floating_point T>
T norm(const Vec3<T>& a)
{
    return std::sqrt(dot(a, a));
}

template <std::floating_point T>
Vec3<T> add(const Vec3<T>& a, T s, const Vec3<T>& b)
{
    return {a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2]};
}

template <std::floating_point T>
Vec3<T> sub(const Vec3<T>& a, const Vec3<T>& b)
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

template <std::floating_point T>
Vec3<T> times(T s, const Vec3<T>& a)
{
    return {s * a[0], s * a[1], s * a[2]};
}

template <std::floating_point T>
Vec3<T> apply(const Mat3<T>& m, const Vec3<T>& a)
{
    return {m[0][0] * a[0] + m[0][1] * a[1] + m[0][2] * a[2], m[1][0] * a[0] + m[1][1] * a[1] + m[1][2] * a[2],
            m[2][0] * a[0] + m[2][1] * a[1] + m[2][2] * a[2]};
}

// ============================================================================
// The curve of one cell
// ============================================================================

template <std::floating_point T>
struct Plane
{
    Vec3<T> n{}; ///< outward
    T d = 0;
};

template <std::floating_point T>
struct Face
{
    Plane<T> plane;
    Vec3<T> origin{};
    std::array<T, 6> matrix{}; ///< 3 x 2: the rectangle s -> origin + matrix s holds the face
};

template <std::floating_point T>
using Path = std::vector<Vec3<T>>;

template <std::floating_point T>
struct Curve
{
    const Source<T>* a = nullptr;
    const Source<T>* b = nullptr;
    std::vector<Plane<T>> planes;
    Mat3<T> jacobian{};
    Vec3<T> ref_origin{};     ///< box coordinates to the cell's reference coordinates
    Mat3<T> ref_jacobian{};
    T diameter = 0;
    T x_tol = 0;      ///< Newton's last step
    T accept_tol = 0; ///< its last step where it stalls at rounding
    T inside_tol = 0; ///< a point counts as in the cell
    T seed_tol = 0;   ///< two points on a face are one
};

template <std::floating_point T>
void value_gradient(const Source<T>& phi, const Vec3<T>& u, T& v, Vec3<T>& g)
{
    v = evaluate(phi, std::span<const T>(u));
    gradient(phi, std::span<const T>(u), std::span<T>(g));
}

template <std::floating_point T>
Vec3<T> tangent(const Curve<T>& c, const Vec3<T>& x)
{
    T va, vb;
    Vec3<T> ga, gb;
    value_gradient(*c.a, x, va, ga);
    value_gradient(*c.b, x, vb, gb);
    return cross(ga, gb);
}

template <std::floating_point T>
bool inside(const Curve<T>& c, const Vec3<T>& x, T tol)
{
    for (const Plane<T>& p : c.planes)
        if (dot(p.n, x) - p.d > tol)
            return false;
    return true;
}

/// Newton's method onto a = b = 0: in the plane normal . x = offset if there is
/// one, else by the smallest step. False where the gradients are parallel.
template <std::floating_point T>
bool correct(const Curve<T>& c, Vec3<T>& x, const Vec3<T>* normal, T offset, int iterations = 12)
{
    const T singular = scaled_tolerance<T>(1e-14);
    T step = std::numeric_limits<T>::infinity();
    for (int it = 0; it < iterations; ++it)
    {
        T va, vb;
        Vec3<T> ga, gb;
        value_gradient(*c.a, x, va, ga);
        value_gradient(*c.b, x, vb, gb);
        const T na = norm(ga), nb = norm(gb);
        if (!(na > T(0)) || !(nb > T(0)) || !std::isfinite(va) || !std::isfinite(vb))
            return false;
        Vec3<T> dx;
        if (normal != nullptr)
        {
            const Vec3<T> c1 = cross(gb, *normal), c2 = cross(*normal, ga), c3 = cross(ga, gb);
            const T det = dot(ga, c1);
            if (!(std::abs(det) > singular * na * nb * norm(*normal)))
                return false;
            const T r = dot(*normal, x) - offset;
            for (int i = 0; i < 3; ++i)
                dx[i] = -(va * c1[i] + vb * c2[i] + r * c3[i]) / det;
        }
        else
        {
            const T aa = dot(ga, ga), ab = dot(ga, gb), bb = dot(gb, gb);
            const T det = aa * bb - ab * ab;
            if (!(det > singular * aa * bb))
                return false;
            const T la = (bb * va - ab * vb) / det, lb = (aa * vb - ab * va) / det;
            for (int i = 0; i < 3; ++i)
                dx[i] = -(la * ga[i] + lb * gb[i]);
        }
        x = add(x, T(1), dx);
        step = norm(dx);
        if (step <= c.x_tol)
            return true;
    }
    return step <= c.accept_tol;
}

/// Newton's method onto a = b = f = 0 from the point @p x.
template <std::floating_point T>
bool correct_three(const Curve<T>& c, const Source<T>& f, Vec3<T>& x)
{
    const T singular = scaled_tolerance<T>(1e-14);
    T step = std::numeric_limits<T>::infinity();
    for (int it = 0; it < 12; ++it)
    {
        T va, vb, vf;
        Vec3<T> ga, gb, gf;
        value_gradient(*c.a, x, va, ga);
        value_gradient(*c.b, x, vb, gb);
        value_gradient(f, x, vf, gf);
        const Vec3<T> c1 = cross(gb, gf), c2 = cross(gf, ga), c3 = cross(ga, gb);
        const T det = dot(ga, c1);
        if (!(std::abs(det) > singular * norm(ga) * norm(gb) * norm(gf)))
            return false;
        Vec3<T> dx;
        for (int i = 0; i < 3; ++i)
            dx[i] = -(va * c1[i] + vb * c2[i] + vf * c3[i]) / det;
        x = add(x, T(1), dx);
        step = norm(dx);
        if (step <= c.x_tol)
            return true;
    }
    return step <= c.accept_tol;
}

// ============================================================================
// Bounds of a level set on a face, or on boxes of the cell
// ============================================================================

template <std::floating_point T>
struct Restricted
{
    const Source<T>* phi = nullptr;
    int m = 0;
    Vec3<T> origin{};
    std::array<T, 9> matrix{};
    BoxBernstein<T> form;                 ///< polynomial: on [0, 1]^m
    std::array<BoxBernstein<T>, 3> deriv; ///< polynomial: its derivatives, degree 0 where it is constant
    T tol = 0;                            ///< a value within this may be a zero
};

template <std::floating_point T>
struct Ranges
{
    T lo = 0, hi = 0;
    bool has_derivatives = false;
    std::array<T, 3> dlo{}, dhi{};

    bool may_vanish() const { return lo <= T(0) && hi >= T(0); }
};

template <std::floating_point T>
void restrict(const Source<T>& phi, const Vec3<T>& origin, std::span<const T> matrix, int m, Restricted<T>& r,
              std::vector<T>& work)
{
    r.phi = &phi;
    r.m = m;
    r.origin = origin;
    std::copy(matrix.begin(), matrix.end(), r.matrix.begin());
    r.tol = scaled_tolerance<T>(1e-13) * reference_magnitude(phi);
    if (phi.is_analytic())
        return;
    restrict_affine(*phi.bernstein, std::span<const T>(origin), matrix, m, r.form, work);
    for (int j = 0; j < m; ++j)
        if (r.form.degree[static_cast<std::size_t>(j)] > 0)
            derivative(r.form, j, r.deriv[static_cast<std::size_t>(j)]);
}

template <std::floating_point T>
void coefficient_range(std::span<const T> coeffs, T& lo, T& hi)
{
    lo = std::numeric_limits<T>::infinity();
    hi = -lo;
    for (const T v : coeffs)
    {
        lo = std::min(lo, v);
        hi = std::max(hi, v);
    }
}

/// Bounds on the part [lo, hi] of [0, 1]^m. Derivatives are with respect to s in
/// [0, 1]^m, whatever the part.
template <std::floating_point T>
void ranges(const Restricted<T>& r, const std::array<T, 3>& lo, const std::array<T, 3>& hi, Ranges<T>& out,
            BoxBernstein<T>& scratch, std::vector<T>& work)
{
    out.has_derivatives = false;
    if (r.phi->is_analytic())
    {
        Vec3<T> origin = r.origin;
        std::array<T, 9> matrix{};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < r.m; ++j)
            {
                const T e = r.matrix[static_cast<std::size_t>(i * r.m + j)];
                origin[i] += e * lo[j];
                matrix[static_cast<std::size_t>(i * r.m + j)] = e * (hi[j] - lo[j]);
            }
        AffineBounds<T> bounds;
        const std::size_t size = static_cast<std::size_t>(3 * r.m);
        if (!affine_bounds<T>(*r.phi, std::span<const T>(origin), std::span<const T>(matrix.data(), size), r.m, bounds))
        {
            out.lo = -std::numeric_limits<T>::infinity();
            out.hi = std::numeric_limits<T>::infinity();
            return;
        }
        out.lo = bounds.sign > 0 ? T(1) : -bounds.magnitude - r.tol;
        out.hi = bounds.sign < 0 ? T(-1) : bounds.magnitude + r.tol;
        if (bounds.has_derivatives)
        {
            out.has_derivatives = true;
            for (int j = 0; j < r.m; ++j)
            {
                const T len = hi[j] - lo[j];
                out.dlo[j] = bounds.dmin[j] / len;
                out.dhi[j] = bounds.dmax[j] / len;
            }
        }
        return;
    }
    const int m = r.m;
    const std::span<const T> l(lo.data(), static_cast<std::size_t>(m)), h(hi.data(), static_cast<std::size_t>(m));
    subdivide(r.form, l, h, scratch, work);
    coefficient_range<T>(scratch.coeffs, out.lo, out.hi);
    out.lo -= r.tol;
    out.hi += r.tol;
    out.has_derivatives = true;
    for (int j = 0; j < m; ++j)
    {
        if (r.form.degree[static_cast<std::size_t>(j)] == 0)
        {
            out.dlo[j] = out.dhi[j] = T(0);
            continue;
        }
        subdivide(r.deriv[static_cast<std::size_t>(j)], l, h, scratch, work);
        coefficient_range<T>(scratch.coeffs, out.dlo[j], out.dhi[j]);
    }
}

template <std::floating_point T>
void interval_product(T a, T b, T c, T d, T& lo, T& hi)
{
    const T p[4] = {a * c, a * d, b * c, b * d};
    lo = *std::min_element(p, p + 4);
    hi = *std::max_element(p, p + 4);
}

/// Whether u v - w z, with each in its interval, keeps away from zero.
template <std::floating_point T>
bool difference_definite(T ul, T uh, T vl, T vh, T wl, T wh, T zl, T zh)
{
    T lo1, hi1, lo2, hi2;
    interval_product(ul, uh, vl, vh, lo1, hi1);
    interval_product(wl, wh, zl, zh, lo2, hi2);
    return lo1 - hi2 > T(0) || hi1 - lo2 < T(0);
}

// ============================================================================
// The curve's points on the faces
// ============================================================================

enum class Search : std::uint8_t
{
    done,
    degenerate, ///< a level set vanishes on the face
    failed
};

template <std::floating_point T>
void add_seed(const Curve<T>& c, const Vec3<T>& x, std::vector<Vec3<T>>& seeds)
{
    for (const Vec3<T>& s : seeds)
        if (norm(sub(s, x)) <= c.seed_tol)
            return;
    seeds.push_back(x);
}

/// All a = b = 0 on a face: the rectangle that holds the face is subdivided until
/// the bounds rule a part out, or show at most one point in it.
template <std::floating_point T>
Search face_seeds(const Curve<T>& c, const Face<T>& f, std::vector<Vec3<T>>& seeds, int& budget)
{
    thread_local Restricted<T> ra, rb;
    thread_local BoxBernstein<T> scratch;
    thread_local std::vector<T> work;
    restrict(*c.a, f.origin, std::span<const T>(f.matrix), 2, ra, work);
    restrict(*c.b, f.origin, std::span<const T>(f.matrix), 2, rb, work);
    // a level set that vanishes on the whole face: the curve, if any, lies in it
    if (!c.a->is_analytic() && max_abs<T>(ra.form.coeffs) <= ra.tol)
        return Search::degenerate;
    if (!c.b->is_analytic() && max_abs<T>(rb.form.coeffs) <= rb.tol)
        return Search::degenerate;

    struct Part
    {
        std::array<T, 3> lo, hi;
        int depth;
    };
    constexpr int max_depth = 13;
    std::vector<Part> stack = {{{T(0), T(0), T(0)}, {T(1), T(1), T(0)}, 0}};
    while (!stack.empty())
    {
        if (--budget < 0)
            return Search::failed;
        const Part p = stack.back();
        stack.pop_back();
        const T mx = T(0.5) * (p.lo[0] + p.hi[0]), my = T(0.5) * (p.lo[1] + p.hi[1]);
        const auto point = [&f](T s, T t)
        {
            return Vec3<T>{f.origin[0] + f.matrix[0] * s + f.matrix[1] * t, f.origin[1] + f.matrix[2] * s + f.matrix[3] * t,
                           f.origin[2] + f.matrix[4] * s + f.matrix[5] * t};
        };
        // outside the face: its corners all beyond one plane
        const std::array<Vec3<T>, 4> corners
            = {point(p.lo[0], p.lo[1]), point(p.hi[0], p.lo[1]), point(p.lo[0], p.hi[1]), point(p.hi[0], p.hi[1])};
        bool outside = false;
        for (const Plane<T>& q : c.planes)
        {
            bool all = true;
            for (const Vec3<T>& x : corners)
                all = all && dot(q.n, x) - q.d > c.inside_tol;
            outside = outside || all;
        }
        if (outside)
            continue;
        Ranges<T> ga, gb;
        ranges(ra, p.lo, p.hi, ga, scratch, work);
        if (!ga.may_vanish())
            continue;
        ranges(rb, p.lo, p.hi, gb, scratch, work);
        if (!gb.may_vanish())
            continue;
        // the 2 x 2 Jacobian keeps its sign: at most one point here
        const bool unique = ga.has_derivatives && gb.has_derivatives
                            && difference_definite(ga.dlo[0], ga.dhi[0], gb.dlo[1], gb.dhi[1], ga.dlo[1], ga.dhi[1],
                                                   gb.dlo[0], gb.dhi[0]);
        bool found = false;
        if (unique || p.depth >= max_depth)
        {
            Vec3<T> x = point(mx, my);
            if (correct(c, x, &f.plane.n, f.plane.d) && inside(c, x, c.inside_tol))
            {
                add_seed(c, x, seeds);
                // in this part: no other
                const Vec3<T> o = sub(x, f.origin);
                const T e1 = f.matrix[0] * f.matrix[0] + f.matrix[2] * f.matrix[2] + f.matrix[4] * f.matrix[4];
                const T e2 = f.matrix[1] * f.matrix[1] + f.matrix[3] * f.matrix[3] + f.matrix[5] * f.matrix[5];
                const T s = (o[0] * f.matrix[0] + o[1] * f.matrix[2] + o[2] * f.matrix[4]) / e1;
                const T t = (o[0] * f.matrix[1] + o[1] * f.matrix[3] + o[2] * f.matrix[5]) / e2;
                found = s >= p.lo[0] && s <= p.hi[0] && t >= p.lo[1] && t <= p.hi[1];
            }
        }
        if ((unique && found) || p.depth >= max_depth)
            continue;
        const int axis = p.hi[0] - p.lo[0] >= p.hi[1] - p.lo[1] ? 0 : 1;
        // off-centre: a point on a dyadic line would lie on the part's edge
        const T mid = p.lo[axis] + T(0.53) * (p.hi[axis] - p.lo[axis]);
        Part first = p, second = p;
        first.hi[axis] = mid;
        second.lo[axis] = mid;
        first.depth = second.depth = p.depth + 1;
        stack.push_back(first);
        stack.push_back(second);
    }
    return Search::done;
}

// ============================================================================
// Following the curve
// ============================================================================

enum class Walk : std::uint8_t
{
    failed,
    left,   ///< reached a face
    closed, ///< came back to its start
    tiny    ///< ends where it starts, within rounding
};

/// From @p start along @p t, inside the cell: a predictor step, Newton in the
/// plane across it, and steps small enough that the tangent turns little. A
/// branch goes straight through a crossing, where the raw tangent flips.
template <std::floating_point T>
Walk follow(const Curve<T>& c, const Vec3<T>& start, Vec3<T> t, bool may_close, Path<T>& path, Vec3<T>& last_tangent)
{
    constexpr int max_steps = 5000;
    const T ds_max = T(0.15) * c.diameter, ds_min = scaled_tolerance<T>(1e-9) * c.diameter;
    const T cos_turn = T(0.9659);
    const Vec3<T> t0 = t;
    path.assign(1, start);
    Vec3<T> x = start;
    T ds = T(0.05) * c.diameter;
    for (int step = 0; step < max_steps; ++step)
    {
        if (ds < ds_min)
            return path.size() == 1 ? Walk::tiny : Walk::failed;
        const Vec3<T> pred = add(x, ds, t);
        Vec3<T> y = pred;
        if (!correct(c, y, &t, dot(t, pred)) || norm(sub(y, x)) > T(2) * ds)
        {
            ds *= T(0.5);
            continue;
        }
        const Vec3<T> raw = tangent(c, y);
        const T n = norm(raw);
        Vec3<T> ty = n > T(0) ? times(T(1) / n, raw) : t;
        if (dot(ty, t) < T(0))
            ty = times(T(-1), ty);
        if (dot(ty, t) < cos_turn)
        {
            ds *= T(0.5);
            continue;
        }
        int leaves = -1;
        T frac = std::numeric_limits<T>::infinity();
        bool back = false;
        for (std::size_t m = 0; m < c.planes.size(); ++m)
        {
            const Plane<T>& p = c.planes[m];
            const T vy = dot(p.n, y) - p.d;
            if (vy <= T(0))
                continue;
            const T vx = dot(p.n, x) - p.d;
            if (vx >= -c.inside_tol && path.size() == 1)
            {
                back = true; // out again through the plane it started on: a smaller step
                break;
            }
            const T f = vx >= T(0) ? T(0) : -vx / (vy - vx);
            if (f < frac)
            {
                frac = f;
                leaves = static_cast<int>(m);
            }
        }
        if (back)
        {
            ds *= T(0.5);
            continue;
        }
        if (leaves >= 0)
        {
            const Plane<T>& p = c.planes[static_cast<std::size_t>(leaves)];
            Vec3<T> e = add(x, frac, sub(y, x));
            if (!correct(c, e, &p.n, p.d) || !inside(c, e, c.inside_tol))
            {
                ds *= T(0.5);
                continue;
            }
            path.push_back(e);
            last_tangent = ty;
            return Walk::left;
        }
        if (may_close && path.size() > 2)
        {
            const Vec3<T> sx = sub(start, x), yx = sub(y, x);
            const T f = dot(sx, yx) / dot(yx, yx);
            if (f > T(0) && f <= T(1) && norm(add(sx, -f, yx)) <= T(0.1) * ds && dot(t, t0) > T(0.5))
            {
                path.push_back(start);
                last_tangent = t;
                return Walk::closed;
            }
        }
        path.push_back(y);
        x = y;
        t = ty;
        ds = std::min(T(1.5) * ds, ds_max);
    }
    return Walk::failed;
}

/// Whether the curve goes from @p x along @p dir into the cell.
template <std::floating_point T>
bool enters(const Curve<T>& c, const Vec3<T>& x, const Vec3<T>& dir)
{
    const T active = scaled_tolerance<T>(1e-9);
    bool touches = false;
    for (const Plane<T>& p : c.planes)
    {
        if (std::abs(dot(p.n, x) - p.d) > active)
            continue;
        const T s = dot(p.n, dir);
        if (s > active)
            return false;
        touches = touches || s >= -active;
    }
    if (!touches)
        return true;
    Vec3<T> y = add(x, scaled_tolerance<T>(1e-5) * c.diameter, dir);
    return correct(c, y, &dir, dot(dir, y)) && inside(c, y, c.seed_tol);
}

// ============================================================================
// Closed curves
// ============================================================================

/// Newton's method on a = b = T_k = 0, an extreme point of the curve along axis k,
/// the Jacobian's last row by differences of the gradients.
template <std::floating_point T>
bool extreme_point(const Curve<T>& c, int k, Vec3<T>& x)
{
    const T h = std::sqrt(std::numeric_limits<T>::epsilon());
    for (int it = 0; it < 20; ++it)
    {
        T va, vb;
        Vec3<T> ga, gb;
        value_gradient(*c.a, x, va, ga);
        value_gradient(*c.b, x, vb, gb);
        const T r3 = cross(ga, gb)[static_cast<std::size_t>(k)];
        Vec3<T> gt;
        for (int j = 0; j < 3; ++j)
        {
            Vec3<T> p = x, m = x;
            p[j] += h;
            m[j] -= h;
            gt[j] = (tangent(c, p)[static_cast<std::size_t>(k)] - tangent(c, m)[static_cast<std::size_t>(k)]) / (T(2) * h);
        }
        const Vec3<T> c1 = cross(gb, gt), c2 = cross(gt, ga), c3 = cross(ga, gb);
        const T det = dot(ga, c1);
        if (!(std::abs(det) > T(0)))
            return false;
        Vec3<T> dx;
        for (int i = 0; i < 3; ++i)
            dx[i] = -(va * c1[i] + vb * c2[i] + r3 * c3[i]) / det;
        x = add(x, T(1), dx);
        if (norm(dx) <= c.accept_tol)
            return true;
    }
    return false;
}

/// Points of the curve from which closed curves inside the cell can be followed.
/// A closed curve has an extreme point along every axis, where T_k = 0: where the
/// bounds keep T_k from vanishing wherever a and b may, for some k, there is none.
template <std::floating_point T>
void closed_curve_seeds(const Curve<T>& c, std::vector<Vec3<T>>& seeds)
{
    thread_local Restricted<T> ra, rb;
    thread_local BoxBernstein<T> scratch;
    thread_local std::vector<T> work;
    const std::array<T, 9> identity = {T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1)};
    restrict(*c.a, Vec3<T>{T(0), T(0), T(0)}, std::span<const T>(identity), 3, ra, work);
    restrict(*c.b, Vec3<T>{T(0), T(0), T(0)}, std::span<const T>(identity), 3, rb, work);

    struct Part
    {
        std::array<T, 3> lo, hi;
        unsigned alive;
    };
    constexpr int max_depth = 5;
    std::vector<Part> level = {{{T(0), T(0), T(0)}, {T(1), T(1), T(1)}, 7u}};
    std::vector<Part> next;
    std::array<int, 3> count{};
    int k_best = 0;
    for (int depth = 0; depth <= max_depth; ++depth)
    {
        next.clear();
        count = {0, 0, 0};
        for (const Part& p : level)
        {
            std::array<Vec3<T>, 8> corners;
            for (int i = 0; i < 8; ++i)
                corners[static_cast<std::size_t>(i)] = {(i & 1) ? p.hi[0] : p.lo[0], (i & 2) ? p.hi[1] : p.lo[1], (i & 4) ? p.hi[2] : p.lo[2]};
            bool outside = false;
            for (const Plane<T>& q : c.planes)
            {
                bool all = true;
                for (const Vec3<T>& x : corners)
                    all = all && dot(q.n, x) - q.d > c.inside_tol;
                outside = outside || all;
            }
            if (outside)
                continue;
            Ranges<T> ga, gb;
            ranges(ra, p.lo, p.hi, ga, scratch, work);
            if (!ga.may_vanish())
                continue;
            ranges(rb, p.lo, p.hi, gb, scratch, work);
            if (!gb.may_vanish())
                continue;
            unsigned alive = p.alive;
            if (ga.has_derivatives && gb.has_derivatives)
                for (int k = 0; k < 3; ++k)
                {
                    const int i = (k + 1) % 3, j = (k + 2) % 3;
                    if ((alive & (1u << k))
                        && difference_definite(ga.dlo[i], ga.dhi[i], gb.dlo[j], gb.dhi[j], ga.dlo[j], ga.dhi[j], gb.dlo[i], gb.dhi[i]))
                        alive &= ~(1u << k);
                }
            if (alive == 0)
                continue;
            for (int k = 0; k < 3; ++k)
                count[static_cast<std::size_t>(k)] += (alive >> k) & 1u;
            Part q = p;
            q.alive = alive;
            next.push_back(q);
        }
        k_best = static_cast<int>(std::min_element(count.begin(), count.end()) - count.begin());
        if (count[static_cast<std::size_t>(k_best)] == 0)
            return;
        if (depth == max_depth)
            break;
        level.clear();
        for (const Part& p : next)
            for (int i = 0; i < 8; ++i)
            {
                Part q = p;
                for (int j = 0; j < 3; ++j)
                {
                    const T mid = T(0.5) * (p.lo[j] + p.hi[j]);
                    (((i >> j) & 1) ? q.lo[j] : q.hi[j]) = mid;
                }
                level.push_back(q);
            }
    }
    for (const Part& p : next)
    {
        if (!((p.alive >> k_best) & 1u))
            continue;
        Vec3<T> x = {T(0.5) * (p.lo[0] + p.hi[0]), T(0.5) * (p.lo[1] + p.hi[1]), T(0.5) * (p.lo[2] + p.hi[2])};
        if (extreme_point(c, k_best, x) && inside(c, x, c.inside_tol))
            add_seed(c, x, seeds);
    }
}

template <std::floating_point T>
bool on_paths(const std::vector<Path<T>>& arcs, const Vec3<T>& x)
{
    for (const Path<T>& path : arcs)
        for (std::size_t i = 0; i + 1 < path.size(); ++i)
        {
            const Vec3<T> d = sub(path[i + 1], path[i]), r = sub(x, path[i]);
            const T dd = dot(d, d);
            if (!(dd > T(0)))
                continue;
            const T f = std::clamp(dot(r, d) / dd, T(0), T(1));
            if (norm(add(r, -f, d)) <= T(0.15) * std::sqrt(dd))
                return true;
        }
    return false;
}

// ============================================================================
// Splits: crossings, further level sets, and arcs that turn too far
// ============================================================================

/// The point between @p x0 and @p x1 where the two gradients are parallel:
/// grad a = mu grad b, a = 0, Newton in (x, mu), the Hessians by differences.
template <std::floating_point T>
Vec3<T> crossing_point(const Curve<T>& c, const Vec3<T>& x0, const Vec3<T>& x1)
{
    const Vec3<T> middle = add(x0, T(0.5), sub(x1, x0));
    const T h = std::sqrt(std::numeric_limits<T>::epsilon());
    Vec3<T> x = middle;
    T va, vb;
    Vec3<T> ga, gb;
    value_gradient(*c.b, x, vb, gb);
    value_gradient(*c.a, x, va, ga);
    T mu = dot(ga, gb) / std::max(dot(gb, gb), std::numeric_limits<T>::min());
    for (int it = 0; it < 30; ++it)
    {
        value_gradient(*c.a, x, va, ga);
        value_gradient(*c.b, x, vb, gb);
        std::array<std::array<T, 5>, 4> m{};
        for (int j = 0; j < 3; ++j)
        {
            Vec3<T> p = x, q = x, gap, gbp, gaq, gbq;
            T v;
            p[j] += h;
            q[j] -= h;
            value_gradient(*c.a, p, v, gap);
            value_gradient(*c.b, p, v, gbp);
            value_gradient(*c.a, q, v, gaq);
            value_gradient(*c.b, q, v, gbq);
            for (int i = 0; i < 3; ++i)
                m[i][j] = ((gap[i] - mu * gbp[i]) - (gaq[i] - mu * gbq[i])) / (T(2) * h);
            m[3][j] = ga[j];
        }
        for (int i = 0; i < 3; ++i)
        {
            m[i][3] = -gb[i];
            m[i][4] = -(ga[i] - mu * gb[i]);
        }
        m[3][4] = -va;
        // Gaussian elimination with pivoting
        bool singular = false;
        for (int col = 0; col < 4 && !singular; ++col)
        {
            int pivot = col;
            for (int i = col + 1; i < 4; ++i)
                if (std::abs(m[i][col]) > std::abs(m[pivot][col]))
                    pivot = i;
            std::swap(m[pivot], m[col]);
            if (!(std::abs(m[col][col]) > T(0)))
            {
                singular = true;
                break;
            }
            for (int i = col + 1; i < 4; ++i)
            {
                const T f = m[i][col] / m[col][col];
                for (int j = col; j < 5; ++j)
                    m[i][j] -= f * m[col][j];
            }
        }
        if (singular)
            return middle;
        std::array<T, 4> z{};
        for (int i = 3; i >= 0; --i)
        {
            T s = m[i][4];
            for (int j = i + 1; j < 4; ++j)
                s -= m[i][j] * z[j];
            z[i] = s / m[i][i];
        }
        for (int i = 0; i < 3; ++i)
            x[i] += z[i];
        mu += z[3];
        const T step = std::sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);
        if (step <= c.accept_tol)
            return norm(sub(x, middle)) <= T(2) * norm(sub(x1, x0)) ? x : middle;
    }
    return middle;
}

/// A path as the pieces between the crossings on it.
template <std::floating_point T>
void split_at_crossings(const Curve<T>& c, const Path<T>& path, std::vector<Path<T>>& pieces)
{
    Path<T> current = {path[0]};
    Vec3<T> previous = tangent(c, path[0]);
    for (std::size_t i = 1; i < path.size(); ++i)
    {
        const Vec3<T> raw = tangent(c, path[i]);
        if (dot(raw, previous) < T(0))
        {
            const Vec3<T> x = crossing_point(c, path[i - 1], path[i]);
            current.push_back(x);
            pieces.push_back(current);
            current.assign(1, x);
        }
        current.push_back(path[i]);
        previous = raw;
    }
    pieces.push_back(current);
}

/// A path as the pieces between the points where @p f vanishes on it.
template <std::floating_point T>
void split_at_zeros(const Curve<T>& c, const Source<T>& f, const Path<T>& path, std::vector<Path<T>>& pieces)
{
    Path<T> current = {path[0]};
    T before = evaluate(f, std::span<const T>(path[0]));
    for (std::size_t i = 1; i < path.size(); ++i)
    {
        // the middle of the step too: a level set can cross and come back within one
        const Vec3<T> mid = add(path[i - 1], T(0.5), sub(path[i], path[i - 1]));
        for (const Vec3<T>* at : {&mid, &path[i]})
        {
            const T v = evaluate(f, std::span<const T>(*at));
            if ((before < T(0)) != (v < T(0)))
            {
                Vec3<T> x = add(current.back(), before / (before - v), sub(*at, current.back()));
                if (correct_three(c, f, x) && norm(sub(x, *at)) <= T(2) * norm(sub(path[i], path[i - 1])))
                {
                    current.push_back(x);
                    pieces.push_back(current);
                    current.assign(1, x);
                }
            }
            before = v;
        }
        current.push_back(path[i]);
    }
    pieces.push_back(current);
}

template <std::floating_point T>
bool term_holds(std::span<const Source<T>> phis, std::span<const SelectionTerm> terms, const Vec3<T>& u)
{
    for (const SelectionTerm& t : terms)
    {
        if ((t.negative_required & t.positive_required) | (t.zero_required & (t.negative_required | t.positive_required)))
            continue;
        bool ok = true;
        for (std::uint64_t bits = t.negative_required | t.positive_required; bits != 0 && ok; bits &= bits - 1)
        {
            const int l = std::countr_zero(bits);
            const T v = evaluate(phis[static_cast<std::size_t>(l)], std::span<const T>(u));
            ok = (t.negative_required >> l) & 1u ? v < T(0) : v > T(0);
        }
        if (ok)
            return true;
    }
    return false;
}

/// Pieces whose chord every step follows within 60 degrees: the chord's coordinate
/// is then a good parameter.
template <std::floating_point T>
void monotone_pieces(const Path<T>& path, std::vector<Path<T>>& pieces)
{
    std::vector<Path<T>> stack = {path};
    while (!stack.empty())
    {
        Path<T> p = std::move(stack.back());
        stack.pop_back();
        if (p.size() < 2)
            continue;
        const Vec3<T> chord = sub(p.back(), p.front());
        const T L = norm(chord);
        bool good = L > T(0);
        for (std::size_t i = 0; i + 1 < p.size() && good; ++i)
        {
            const Vec3<T> step = sub(p[i + 1], p[i]);
            good = dot(step, chord) > T(0.5) * norm(step) * L;
        }
        if (good || p.size() == 2)
        {
            pieces.push_back(std::move(p));
            continue;
        }
        const std::size_t m = p.size() / 2;
        stack.emplace_back(p.begin() + static_cast<std::ptrdiff_t>(m), p.end());
        stack.emplace_back(p.begin(), p.begin() + static_cast<std::ptrdiff_t>(m) + 1);
    }
}

// ============================================================================
// The rule on an arc
// ============================================================================

template <std::floating_point T>
struct ArcFit
{
    Vec3<T> origin{}, chord{};
    T length = 0; ///< of the chord, in box coordinates
    T total = 0;  ///< of the arc, physical
    std::vector<T> tau;
    std::vector<T> rate;     ///< Chebyshev series of ds / dtau on the chord
    std::vector<T> integral; ///< of its integral, from the chord's start
};

template <std::floating_point T>
T chebyshev(const std::vector<T>& a, T x)
{
    T b1 = 0, b2 = 0;
    for (std::size_t k = a.size(); k-- > 1;)
    {
        const T b = T(2) * x * b1 - b2 + a[k];
        b2 = b1;
        b1 = b;
    }
    return x * b1 - b2 + a[0];
}

template <std::floating_point T>
bool arc_point(const Curve<T>& c, const Path<T>& path, const ArcFit<T>& f, T tau, Vec3<T>& x)
{
    const std::size_t n = f.tau.size();
    std::size_t i = static_cast<std::size_t>(std::upper_bound(f.tau.begin(), f.tau.end(), tau) - f.tau.begin());
    i = i == 0 ? 0 : std::min(i - 1, n - 2);
    const T span = f.tau[i + 1] - f.tau[i];
    const T s = span > T(0) ? (tau - f.tau[i]) / span : T(0);
    const Vec3<T> guess = add(path[i], s, sub(path[i + 1], path[i]));
    x = guess;
    if (!correct(c, x, &f.chord, dot(f.chord, f.origin) + tau))
        return false;
    return norm(sub(x, guess)) <= T(2) * norm(sub(path[i + 1], path[i])) + c.seed_tol;
}

/// ds / dtau from a Chebyshev fit in the coordinate tau along the chord, and its
/// integral, the arc's length. False if the series has not settled: split the arc.
template <std::floating_point T>
bool fit_arc(const Curve<T>& c, const Path<T>& path, ArcFit<T>& f)
{
    constexpr int n = 16;
    f.origin = path.front();
    const Vec3<T> chord = sub(path.back(), path.front());
    f.length = norm(chord);
    if (!(f.length > T(0)))
        return false;
    f.chord = times(T(1) / f.length, chord);
    f.tau.resize(path.size());
    for (std::size_t i = 0; i < path.size(); ++i)
        f.tau[i] = dot(f.chord, sub(path[i], f.origin));
    for (std::size_t i = 1; i < path.size(); ++i)
        if (!(f.tau[i] >= f.tau[i - 1]))
            return false;
    const T pi = std::acos(T(-1));
    std::array<T, n> values{};
    for (int i = 0; i < n; ++i)
    {
        const T theta = pi * (T(i) + T(0.5)) / T(n);
        const T tau = T(0.5) * f.length * (std::cos(theta) + T(1));
        Vec3<T> x;
        if (!arc_point(c, path, f, tau, x))
            return false;
        const Vec3<T> t = tangent(c, x);
        const T along = std::abs(dot(t, f.chord)), nt = norm(t);
        if (!(along > T(0.2) * nt))
            return false;
        values[static_cast<std::size_t>(i)] = norm(apply(c.jacobian, t)) / along;
    }
    f.rate.assign(n, T(0));
    for (int k = 0; k < n; ++k)
    {
        T sum = 0;
        for (int i = 0; i < n; ++i)
            sum += values[static_cast<std::size_t>(i)] * std::cos(T(k) * pi * (T(i) + T(0.5)) / T(n));
        f.rate[static_cast<std::size_t>(k)] = T(2) * sum / T(n);
    }
    f.rate[0] *= T(0.5);
    // integral of the series against x in [-1, 1], times dtau / dx
    f.integral.assign(n + 1, T(0));
    const auto a = [&f](int k) { return k < n ? f.rate[static_cast<std::size_t>(k)] : T(0); };
    f.integral[1] = a(0) - T(0.5) * a(2);
    for (int k = 2; k <= n; ++k)
        f.integral[static_cast<std::size_t>(k)] = (a(k - 1) - a(k + 1)) / T(2 * k);
    for (T& v : f.integral)
        v *= T(0.5) * f.length;
    f.integral[0] = -chebyshev(f.integral, T(-1));
    f.total = chebyshev(f.integral, T(1));
    const T tail = std::abs(f.rate[n - 1]) + std::abs(f.rate[n - 2]);
    return f.total > T(0) && tail <= scaled_tolerance<T>(1e-10) * std::abs(f.rate[0]);
}

/// Gauss-Legendre points in arc length, from the fit.
template <std::floating_point T>
bool arc_nodes(const Curve<T>& c, const Path<T>& path, const ArcFit<T>& f, int q, std::vector<T>& x, std::vector<T>& w)
{
    thread_local std::vector<std::vector<T>> xs, ws;
    if (static_cast<int>(xs.size()) <= q)
    {
        xs.resize(static_cast<std::size_t>(q) + 1);
        ws.resize(static_cast<std::size_t>(q) + 1);
    }
    std::vector<T>& gx = xs[static_cast<std::size_t>(q)];
    std::vector<T>& gw = ws[static_cast<std::size_t>(q)];
    if (gx.empty())
        gauss_legendre<T>(q, gx, gw);
    x.resize(3 * static_cast<std::size_t>(q));
    w.resize(static_cast<std::size_t>(q));
    for (int j = 0; j < q; ++j)
    {
        const std::size_t js = static_cast<std::size_t>(j);
        const T target = f.total * gx[js];
        T lo = T(-1), hi = T(1);
        T s = T(2) * gx[js] - T(1);
        for (int it = 0; it < 40; ++it)
        {
            const T r = chebyshev(f.integral, s) - target;
            if (std::abs(r) <= scaled_tolerance<T>(1e-15) * f.total)
                break;
            (r > T(0) ? hi : lo) = s;
            const T slope = T(0.5) * f.length * chebyshev(f.rate, s);
            T next = s - r / slope;
            if (!(next > lo && next < hi))
                next = T(0.5) * (lo + hi);
            s = next;
        }
        Vec3<T> u;
        if (!arc_point(c, path, f, T(0.5) * f.length * (s + T(1)), u))
            return false;
        for (std::size_t i = 0; i < 3; ++i)
            x[3 * js + i] = u[i];
        w[js] = f.total * gw[js];
    }
    return true;
}

/// A moment of the mass matrix: the sum over the cell's reference coordinates, which
/// assembly sees, of a Legendre polynomial squared (Options::moment_tolerance). Box
/// coordinates would not do: a cell's box may be sheared against them.
template <std::floating_point T>
T moment(const Curve<T>& c, const std::vector<T>& x, const std::vector<T>& w, int degree)
{
    T sum = 0;
    for (std::size_t i = 0; i < w.size(); ++i)
    {
        T shape = 0;
        for (std::size_t j = 0; j < 3; ++j)
        {
            T xi = c.ref_origin[j];
            for (std::size_t k = 0; k < 3; ++k)
                xi += c.ref_jacobian[j][k] * x[3 * i + k];
            const T t = T(2) * xi - T(1);
            T low = T(1), high = t;
            for (int m = 1; m < degree; ++m)
            {
                const T next = (T(2 * m + 1) * t * high - T(m) * low) / T(m + 1);
                low = high;
                high = next;
            }
            shape += high * high;
        }
        sum += w[i] * shape;
    }
    return sum;
}

template <std::floating_point T>
void emit(const std::vector<T>& x, const std::vector<T>& w, CellPoints<T>& out)
{
    out.points.insert(out.points.end(), x.begin(), x.end());
    out.weights.insert(out.weights.end(), w.begin(), w.end());
}

/// The rule on one arc: @p q points or more until the moment settles; an arc whose
/// fit or moment does not settle at all is split, as far as @p depth allows.
template <std::floating_point T>
bool integrate_arc(const Curve<T>& c, const Path<T>& path, int q, T moment_relative, int degree, int depth, CellPoints<T>& out)
{
    constexpr int max_depth = 5;
    if (path.size() < 2)
        return true;
    const auto halves = [&]()
    {
        Path<T> p = path;
        if (p.size() == 2)
        {
            Vec3<T> m = add(p[0], T(0.5), sub(p[1], p[0]));
            const Vec3<T> d = sub(p[1], p[0]);
            if (!correct(c, m, &d, dot(d, m)))
                return false;
            p.insert(p.begin() + 1, m);
        }
        const std::size_t m = p.size() / 2;
        return integrate_arc(c, Path<T>(p.begin(), p.begin() + static_cast<std::ptrdiff_t>(m) + 1), q, moment_relative, degree, depth + 1, out)
               && integrate_arc(c, Path<T>(p.begin() + static_cast<std::ptrdiff_t>(m), p.end()), q, moment_relative, degree, depth + 1, out);
    };
    if (norm(sub(path.back(), path.front())) <= c.x_tol)
        return true;
    ArcFit<T> f;
    if (!fit_arc(c, path, f) && depth < max_depth)
        return halves();
    if (!(f.total > T(0)))
        return false;
    thread_local std::vector<T> x0, w0, x1, w1;
    if (!arc_nodes(c, path, f, q, x0, w0))
        return depth < max_depth && halves();
    if (moment_relative > T(0))
    {
        const int q_max = 2 * q + 4;
        T m0 = moment(c, x0, w0, degree);
        bool settled = false;
        for (int r = q; r < q_max; ++r)
        {
            if (!arc_nodes(c, path, f, r + 1, x1, w1))
                return depth < max_depth && halves();
            const T m1 = moment(c, x1, w1, degree);
            if (std::abs(m1 - m0) <= moment_relative * f.total)
            {
                settled = true;
                break;
            }
            x0.swap(x1);
            w0.swap(w1);
            m0 = m1;
        }
        if (!settled && depth < max_depth)
            return halves();
    }
    emit(x0, w0, out);
    return true;
}

template <std::floating_point T>
void make_faces(const Polytope<T>& poly, std::vector<Plane<T>>& planes, std::vector<Face<T>>& faces)
{
    planes.clear();
    faces.clear();
    for (int k = 0; k < poly.n_faces(); ++k)
    {
        const Vec3<T>& n = poly.face_normals[static_cast<std::size_t>(k)];
        const std::size_t first = static_cast<std::size_t>(poly.face_offsets[static_cast<std::size_t>(k)]);
        const std::size_t last = static_cast<std::size_t>(poly.face_offsets[static_cast<std::size_t>(k) + 1]);
        if (first == last)
            continue;
        const std::array<int, 2>& e0 = poly.edges[static_cast<std::size_t>(poly.face_edges[first])];
        const Vec3<T>& v0 = poly.vertices[static_cast<std::size_t>(e0[0])];
        const Vec3<T> along = sub(poly.vertices[static_cast<std::size_t>(e0[1])], v0);
        const T len = norm(along);
        if (!(len > T(0)))
            continue;
        const Vec3<T> e1 = times(T(1) / len, along), e2 = cross(n, e1);
        T lo1 = 0, hi1 = 0, lo2 = 0, hi2 = 0;
        for (std::size_t i = first; i < last; ++i)
            for (const int v : poly.edges[static_cast<std::size_t>(poly.face_edges[i])])
            {
                const Vec3<T> r = sub(poly.vertices[static_cast<std::size_t>(v)], v0);
                lo1 = std::min(lo1, dot(r, e1));
                hi1 = std::max(hi1, dot(r, e1));
                lo2 = std::min(lo2, dot(r, e2));
                hi2 = std::max(hi2, dot(r, e2));
            }
        Face<T> f;
        f.plane = {n, dot(n, v0)};
        f.origin = add(add(v0, lo1, e1), lo2, e2);
        for (std::size_t i = 0; i < 3; ++i)
        {
            f.matrix[2 * i] = (hi1 - lo1) * e1[i];
            f.matrix[2 * i + 1] = (hi2 - lo2) * e2[i];
        }
        planes.push_back(f.plane);
        faces.push_back(f);
    }
}

} // namespace

template <std::floating_point T>
bool curve_points(const ClippedBox<T>& cell, std::span<const Source<T>> phis, std::span<const SelectionTerm> terms,
                  int a, int b, int q, const Options& opt, CellPoints<T>& out, Stats& stats)
{
    const std::size_t start_points = out.points.size(), start_weights = out.weights.size();
    const auto failed = [&](bool count)
    {
        out.points.resize(start_points);
        out.weights.resize(start_weights);
        stats.curve_fallbacks += count;
        return false;
    };

    thread_local Polytope<T> poly;
    thread_local std::vector<Face<T>> faces;
    Curve<T> c;
    c.a = &phis[static_cast<std::size_t>(a)];
    c.b = &phis[static_cast<std::size_t>(b)];
    clipped_polytope(cell, poly);
    make_faces(poly, c.planes, faces);
    c.jacobian = cell.jacobian;
    c.ref_origin = cell.ref_origin;
    c.ref_jacobian = cell.ref_jacobian;
    for (int i = 0; i < poly.n_vertices(); ++i)
        for (int j = 0; j < i; ++j)
            c.diameter = std::max(c.diameter, norm(sub(poly.vertices[static_cast<std::size_t>(i)], poly.vertices[static_cast<std::size_t>(j)])));
    if (!(c.diameter > T(0)))
        return failed(true);
    c.x_tol = scaled_tolerance<T>(1e-14);
    c.accept_tol = scaled_tolerance<T>(1e-11);
    c.inside_tol = scaled_tolerance<T>(1e-10);
    c.seed_tol = scaled_tolerance<T>(1e-9);

    // the curve's points on the faces
    thread_local std::vector<Vec3<T>> seeds;
    seeds.clear();
    int budget = 6000;
    for (const Face<T>& f : faces)
    {
        const Search s = face_seeds(c, f, seeds, budget);
        if (s == Search::degenerate)
            return failed(false);
        if (s == Search::failed)
            return failed(true);
    }

    // follow the curve from each, into the cell
    thread_local std::vector<Path<T>> arcs;
    arcs.clear();
    std::vector<Path<T>> consumed(seeds.size());
    Path<T> path;
    Vec3<T> last;
    const T match = scaled_tolerance<T>(1e-7);
    for (std::size_t s = 0; s < seeds.size(); ++s)
    {
        const Vec3<T> raw = tangent(c, seeds[s]);
        const T n = norm(raw);
        if (!(n > T(0)))
            continue;
        const Vec3<T> t = times(T(1) / n, raw);
        for (const T sign : {T(1), T(-1)})
        {
            const Vec3<T> dir = times(sign, t);
            bool taken = false;
            for (const Vec3<T>& d : consumed[s])
                taken = taken || dot(d, dir) > T(0.9);
            if (taken || !enters(c, seeds[s], dir))
                continue;
            const Walk w = follow(c, seeds[s], dir, false, path, last);
            if (w == Walk::failed)
                return failed(true);
            if (w == Walk::tiny)
                continue;
            consumed[s].push_back(dir);
            for (std::size_t r = 0; r < seeds.size(); ++r)
                if (norm(sub(seeds[r], path.back())) <= match)
                    consumed[r].push_back(times(T(-1), last));
            arcs.push_back(path);
        }
    }

    // closed curves inside the cell
    thread_local std::vector<Vec3<T>> extremes;
    extremes.clear();
    closed_curve_seeds(c, extremes);
    for (const Vec3<T>& x : extremes)
    {
        if (on_paths(arcs, x))
            continue;
        const Vec3<T> raw = tangent(c, x);
        const T n = norm(raw);
        if (!(n > T(0)))
            continue;
        const Vec3<T> t = times(T(1) / n, raw);
        const Walk w = follow(c, x, t, true, path, last);
        if (w == Walk::failed)
            return failed(true);
        if (w == Walk::tiny)
            continue;
        if (w == Walk::left)
        {
            // reached a face the search did not find: the other way too
            Path<T> other;
            const Walk w2 = follow(c, x, times(T(-1), t), false, other, last);
            if (w2 == Walk::failed)
                return failed(true);
            if (w2 == Walk::left)
            {
                std::reverse(other.begin(), other.end());
                other.insert(other.end(), path.begin() + 1, path.end());
                path.swap(other);
            }
        }
        arcs.push_back(path);
    }

    // level sets other than the two that the terms name
    std::uint64_t others = 0;
    for (const SelectionTerm& term : terms)
        others |= term.negative_required | term.positive_required;
    others &= ~((std::uint64_t(1) << a) | (std::uint64_t(1) << b));

    const bool check = opt.moment_tolerance > 0 && opt.tolerance >= 0;
    const T moment_relative = check ? static_cast<T>(std::pow(opt.moment_tolerance, std::min(q, 6) / 6.0)) : T(0);
    const int degree = std::max(1, (2 * q - 1) / 2);

    std::vector<Path<T>> split, current, next, monotone;
    for (const Path<T>& arc : arcs)
    {
        current.clear();
        split_at_crossings(c, arc, current);
        for (std::uint64_t bits = others; bits != 0; bits &= bits - 1)
        {
            const Source<T>& f = phis[static_cast<std::size_t>(std::countr_zero(bits))];
            next.clear();
            for (const Path<T>& p : current)
                split_at_zeros(c, f, p, next);
            current.swap(next);
        }
        for (const Path<T>& p : current)
        {
            const Vec3<T> mid = p.size() > 2 ? p[p.size() / 2] : add(p[0], T(0.5), sub(p[1], p[0]));
            if (!term_holds<T>(phis, terms, mid))
                continue;
            monotone.clear();
            monotone_pieces(p, monotone);
            for (const Path<T>& piece : monotone)
            {
                if (!integrate_arc(c, piece, q, moment_relative, degree, 0, out))
                    return failed(true);
                ++stats.curve_arcs;
            }
        }
    }
    return true;
}

template bool curve_points<float>(const ClippedBox<float>&, std::span<const Source<float>>, std::span<const SelectionTerm>,
                                  int, int, int, const Options&, CellPoints<float>&, Stats&);
template bool curve_points<double>(const ClippedBox<double>&, std::span<const Source<double>>, std::span<const SelectionTerm>,
                                   int, int, int, const Options&, CellPoints<double>&, Stats&);

} // namespace cutcells::quadrays
