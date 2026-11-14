// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "slicing.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "gauss.h"
#include "source.h"

namespace cutcells::quadrays
{

namespace
{

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
Vec3<T> unit(const Vec3<T>& a)
{
    const T n = std::sqrt(dot(a, a));
    if (!(n > T(0)))
        throw std::invalid_argument("quadrays: a sweep needs a direction");
    return {a[0] / n, a[1] / n, a[2] / n};
}

template <std::floating_point T>
Vec3<T> times(const Mat3<T>& m, const Vec3<T>& v)
{
    return {dot(m[0], v), dot(m[1], v), dot(m[2], v)};
}

/// Two unit vectors that complete the unit vector @p n to a right-handed
/// orthonormal frame (a, b, n), the same for the same n.
template <std::floating_point T>
void frame_of(const Vec3<T>& n, Vec3<T>& a, Vec3<T>& b)
{
    int least = 0;
    for (int i = 1; i < 3; ++i)
        if (std::abs(n[i]) < std::abs(n[least]))
            least = i;
    Vec3<T> axis = {0, 0, 0};
    axis[least] = T(1);
    a = unit(cross(n, axis));
    b = cross(n, a);
}

/// Half-plane a p + b q <= d of the slice's plane coordinates.
template <std::floating_point T>
struct HalfPlane
{
    T a = 0, b = 0, d = 0;
};

/// The part of the convex polygon @p poly (p, q pairs) in @p h, into @p out.
template <std::floating_point T>
void clip_polygon(const std::vector<std::array<T, 2>>& poly, const HalfPlane<T>& h,
                  std::vector<std::array<T, 2>>& out)
{
    out.clear();
    const std::size_t n = poly.size();
    for (std::size_t k = 0; k < n; ++k)
    {
        const std::array<T, 2>& p = poly[k];
        const std::array<T, 2>& r = poly[(k + 1) % n];
        const T fp = h.a * p[0] + h.b * p[1] - h.d;
        const T fr = h.a * r[0] + h.b * r[1] - h.d;
        if (fp <= T(0))
            out.push_back(p);
        if ((fp < T(0) && fr > T(0)) || (fr < T(0) && fp > T(0)))
        {
            const T s = fp / (fp - fr);
            out.push_back({p[0] + s * (r[0] - p[0]), p[1] + s * (r[1] - p[1])});
        }
    }
}

/// The zero sets the terms ask for, as bits: none for a volume, one for a
/// surface, two for a curve.
/// @throws std::invalid_argument if they ask for more than two, or differ
inline std::uint64_t zero_sets_of(std::span<const SelectionTerm> terms)
{
    if (terms.empty())
        throw std::invalid_argument("quadrays: slices need a selection");
    const std::uint64_t zero = terms[0].zero_required;
    for (const SelectionTerm& term : terms)
        if (term.zero_required != zero)
            throw std::invalid_argument("quadrays: the terms of a sliced part ask for the same zero sets");
    if (std::popcount(zero) > 2)
        throw std::invalid_argument("quadrays: slices integrate volumes, surfaces and curves");
    return zero;
}

/// The lowest level set of @p mask, and the next: -1 where there is none.
inline std::array<int, 2> level_sets_of(std::uint64_t mask)
{
    std::array<int, 2> l = {-1, -1};
    for (int k = 0, n = 0; k < 64 && n < 2; ++k)
        if ((mask >> k) & 1)
            l[static_cast<std::size_t>(n++)] = k;
    return l;
}

/// The rectangle [lo, hi] of a slice's box coordinates that holds the part's
/// closure in the slice, with some room around it: boxes are bisected and
/// dropped where a level set the terms ask a sign of has the other one, a zero
/// set the part lies on cannot vanish, or a clip excludes them.
/// @return 0 if the part misses the slice, 1 if its rectangle has more than a
///         tenth of the box's area (@p lo and @p hi are then the whole box),
///         2 otherwise
template <std::floating_point T>
int fit_to_part(const ClippedBox<T>& box, std::span<const Source<T>> phis, std::span<const SelectionTerm> terms,
                std::uint64_t zero, std::array<T, 2>& lo, std::array<T, 2>& hi)
{
    struct Rect
    {
        std::array<T, 2> lo, hi;
    };
    thread_local std::vector<Rect> kept, next;
    thread_local std::vector<T> size, work;
    thread_local std::vector<std::array<T, 2>> range;
    thread_local BoxBernstein<T> scratch;
    const std::size_t n = phis.size();
    size.resize(n);
    range.resize(n);
    for (std::size_t l = 0; l < n; ++l)
        size[l] = reference_magnitude(phis[l]);
    const T eps = T(256) * std::numeric_limits<T>::epsilon();
    const auto may_meet = [&](const Rect& r)
    {
        for (const HalfSpace<T>& h : box.clips)
            if (std::min(h.c[0] * r.lo[0], h.c[0] * r.hi[0]) + std::min(h.c[1] * r.lo[1], h.c[1] * r.hi[1]) > h.d + eps)
                return false;
        const T w0 = r.hi[0] - r.lo[0], w1 = r.hi[1] - r.lo[1];
        for (std::size_t l = 0; l < n; ++l)
        {
            const T all = std::numeric_limits<T>::infinity();
            range[l] = {-all, all};
            if (phis[l].is_analytic())
            {
                const std::array<T, 3> origin = {r.lo[0], r.lo[1], T(0)};
                const std::array<T, 6> matrix = {w0, 0, 0, w1, 0, 0};
                AffineBounds<T> bounds;
                if (affine_bounds(phis[l], std::span<const T>(origin), std::span<const T>(matrix), 2, bounds)
                    && bounds.sign != 0)
                    range[l] = bounds.sign > 0 ? std::array<T, 2>{all, all} : std::array<T, 2>{-all, -all};
            }
            else
            {
                const std::array<T, 4> matrix = {w0, 0, 0, w1};
                restrict_affine(*phis[l].bernstein, std::span<const T>(r.lo), std::span<const T>(matrix), 2, scratch, work);
                const auto [mn, mx] = std::minmax_element(scratch.coeffs.begin(), scratch.coeffs.end());
                range[l] = {*mn, *mx};
            }
            const T tol = eps * size[l];
            if (((zero >> l) & 1) && (range[l][0] > tol || range[l][1] < -tol))
                return false;
        }
        for (const SelectionTerm& term : terms)
        {
            bool possible = true;
            for (std::size_t l = 0; l < n && possible; ++l)
            {
                const std::uint64_t bit = std::uint64_t(1) << l;
                const T tol = eps * size[l];
                possible = !((term.negative_required & bit) && range[l][0] > tol)
                           && !((term.positive_required & bit) && range[l][1] < -tol);
            }
            if (possible)
                return true;
        }
        return false;
    };
    kept.assign(1, Rect{{T(0), T(0)}, {T(1), T(1)}});
    for (int depth = 0; depth < 16 && kept.size() <= 160; ++depth)
    {
        next.clear();
        for (const Rect& r : kept)
        {
            const T m0 = T(0.5) * (r.lo[0] + r.hi[0]), m1 = T(0.5) * (r.lo[1] + r.hi[1]);
            const Rect children[4] = {{{r.lo[0], r.lo[1]}, {m0, m1}}, {{m0, r.lo[1]}, {r.hi[0], m1}},
                                      {{r.lo[0], m1}, {m0, r.hi[1]}}, {{m0, m1}, {r.hi[0], r.hi[1]}}};
            for (const Rect& child : children)
                if (may_meet(child))
                    next.push_back(child);
        }
        kept.swap(next);
        if (kept.empty())
            return 0;
    }
    lo = {T(1), T(1)};
    hi = {T(0), T(0)};
    for (const Rect& r : kept)
        for (std::size_t j = 0; j < 2; ++j)
        {
            lo[j] = std::min(lo[j], r.lo[j]);
            hi[j] = std::max(hi[j], r.hi[j]);
        }
    // room around the part: the engine's boxes should not hug its boundary
    const T leaf = kept[0].hi[0] - kept[0].lo[0];
    for (std::size_t j = 0; j < 2; ++j)
    {
        const T room = std::max(T(0.2) * (hi[j] - lo[j]), leaf);
        lo[j] = std::max(T(0), lo[j] - room);
        hi[j] = std::min(T(1), hi[j] + room);
    }
    // only a part that is small in its slice gets its own cell: the engine's rule depends on its
    // cell, and a cell that follows the part changes from slice to slice, which the rule along the
    // sweep sees as noise
    if ((hi[0] - lo[0]) * (hi[1] - lo[1]) > T(0.1))
    {
        lo = {T(0), T(0)};
        hi = {T(1), T(1)};
        return 1;
    }
    return 2;
}

} // namespace

template <std::floating_point T>
bool make_slice(const ClippedBox<T>& cell, const Sweep<T>& sweep, T t, Slice<T>& slice)
{
    if (cell.tdim != 3)
        throw std::invalid_argument("quadrays: slices are of cells of dimension 3");
    // the plane: x = base + p e1 + q e2, and for a turning sweep p >= 0 only
    const Vec3<T> n = unit(sweep.direction);
    Vec3<T> a, b;
    frame_of(n, a, b);
    Vec3<T> base = sweep.origin, e1 = a, e2 = b;
    const bool turning = sweep.kind == SweepKind::turning;
    if (turning)
    {
        const T c = std::cos(t), s = std::sin(t);
        e1 = {c * a[0] + s * b[0], c * a[1] + s * b[1], c * a[2] + s * b[2]};
        e2 = n;
    }
    else
        base = {base[0] + t * n[0], base[1] + t * n[1], base[2] + t * n[2]};

    // the plane in the cell's box coordinates: u = u0 + p b1 + q b2
    const Mat3<T> inverse = inverse_jacobian(cell);
    const Vec3<T> u0 = times(inverse, Vec3<T>{base[0] - cell.origin[0], base[1] - cell.origin[1], base[2] - cell.origin[2]});
    const Vec3<T> b1 = times(inverse, e1), b2 = times(inverse, e2);

    // a rectangle of the plane that holds the cell's box, from the box's corners
    T lo[2] = {std::numeric_limits<T>::max(), std::numeric_limits<T>::max()};
    T hi[2] = {std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest()};
    for (int corner = 0; corner < 8; ++corner)
    {
        const Vec3<T> x = physical_point(cell, Vec3<T>{T(corner & 1), T((corner >> 1) & 1), T((corner >> 2) & 1)});
        const Vec3<T> r = {x[0] - base[0], x[1] - base[1], x[2] - base[2]};
        const T pq[2] = {dot(r, e1), dot(r, e2)};
        for (int j = 0; j < 2; ++j)
        {
            lo[j] = std::min(lo[j], pq[j]);
            hi[j] = std::max(hi[j], pq[j]);
        }
    }
    if (turning)
        lo[0] = std::max(lo[0], T(0));
    const T size = std::max(hi[0] - lo[0], hi[1] - lo[1]);
    if (!(hi[0] > lo[0]) || !(hi[1] > lo[1]))
        return false;

    // the cell in the plane: the faces of its box and its clips
    thread_local std::vector<HalfPlane<T>> faces;
    thread_local std::vector<std::array<T, 2>> poly, clipped;
    faces.clear();
    const auto add = [&](const Vec3<T>& c, T d)
    {
        const T norm = std::sqrt(dot(c, c));
        faces.push_back({dot(c, b1) / norm, dot(c, b2) / norm, (d - dot(c, u0)) / norm});
    };
    for (int i = 0; i < 3; ++i)
    {
        Vec3<T> c = {0, 0, 0};
        c[i] = T(1);
        add(c, T(1));
        c[i] = T(-1);
        add(c, T(0));
    }
    for (const HalfSpace<T>& h : cell.clips)
        add(h.c, h.d);
    poly = {{lo[0], lo[1]}, {hi[0], lo[1]}, {hi[0], hi[1]}, {lo[0], hi[1]}};
    for (const HalfPlane<T>& h : faces)
    {
        clip_polygon(poly, h, clipped);
        poly.swap(clipped);
        if (poly.size() < 3)
            return false;
    }
    T area = 0;
    for (std::size_t k = 0; k < poly.size(); ++k)
    {
        const std::array<T, 2>& p = poly[k];
        const std::array<T, 2>& r = poly[(k + 1) % poly.size()];
        area += p[0] * r[1] - p[1] * r[0];
    }
    if (!(std::abs(area) > scaled_tolerance<T>(1e-13) * size * size))
        return false;

    // the slice's box: the bounding rectangle of that polygon
    T p0 = poly[0][0], p1 = poly[0][0], q0 = poly[0][1], q1 = poly[0][1];
    for (const std::array<T, 2>& p : poly)
    {
        p0 = std::min(p0, p[0]);
        p1 = std::max(p1, p[0]);
        q0 = std::min(q0, p[1]);
        q1 = std::max(q1, p[1]);
    }
    const T dp = p1 - p0, dq = q1 - q0;
    if (!(dp > scaled_tolerance<T>(1e-12) * size) || !(dq > scaled_tolerance<T>(1e-12) * size))
        return false;
    for (int i = 0; i < 3; ++i)
    {
        slice.origin[i] = u0[i] + p0 * b1[i] + q0 * b2[i];
        slice.matrix[2 * i] = dp * b1[i];
        slice.matrix[2 * i + 1] = dq * b2[i];
    }
    slice.radius_origin = turning ? p0 : T(0);
    slice.radius_step = turning ? dp : T(0);
    const Vec3<T> normal = cross(e1, e2);
    ClippedBox<T>& box = slice.box;
    box.tdim = 2;
    box.ref_origin = {0, 0, 0};
    for (int i = 0; i < 3; ++i)
    {
        box.origin[i] = base[i] + p0 * e1[i] + q0 * e2[i];
        box.jacobian[i] = {dp * e1[i], dq * e2[i], normal[i]};
        box.ref_jacobian[i] = {T(i == 0), T(i == 1), T(i == 2)};
    }
    // the faces that cut the rectangle, in its box coordinates
    box.clips.clear();
    const T tol = scaled_tolerance<T>(1e-13) * size;
    for (const HalfPlane<T>& h : faces)
    {
        const T c0 = h.a * dp, c1 = h.b * dq, d = h.d - h.a * p0 - h.b * q0;
        if (std::max(c0, T(0)) + std::max(c1, T(0)) <= d + tol)
            continue;
        const T norm = std::sqrt(c0 * c0 + c1 * c1);
        box.clips.push_back({{c0 / norm, c1 / norm, T(0)}, d / norm});
    }
    return true;
}

template <std::floating_point T>
void integrate_slices(const ClippedBox<T>& cell, std::span<const Source<T>> phis,
                      std::span<const SelectionTerm> terms, const Sweep<T>& sweep, std::span<const T> positions,
                      std::span<const T> weights, int q, const Options& opt, CellPoints<T>& out, Stats& stats)
{
    if (positions.size() != weights.size())
        throw std::invalid_argument("quadrays: slices need one weight per position");
    const std::uint64_t zero_sets = zero_sets_of(terms);
    const std::array<int, 2> zero = level_sets_of(zero_sets);
    const Mat3<T> inverse = inverse_jacobian(cell);
    // the physical gradient of level set l at the box coordinates u
    const auto physical_gradient = [&](int l, const Vec3<T>& u)
    {
        Vec3<T> du;
        gradient(phis[static_cast<std::size_t>(l)], std::span<const T>(u), std::span<T>(du));
        Vec3<T> dx = {0, 0, 0};
        for (std::size_t r = 0; r < 3; ++r)
            for (std::size_t c = 0; c < 3; ++c)
                dx[c] += inverse[r][c] * du[r];
        return dx;
    };
    thread_local Slice<T> slice;
    thread_local std::vector<BoxBernstein<T>> restricted, fitted;
    thread_local std::vector<Source<T>> sources, fitted_sources;
    thread_local ClippedBox<T> fitted_box;
    thread_local CellPoints<T> points;
    thread_local std::vector<T> work;
    restricted.resize(phis.size());
    for (std::size_t k = 0; k < positions.size(); ++k)
    {
        if (!make_slice(cell, sweep, positions[k], slice))
            continue;
        sources.clear();
        for (std::size_t l = 0; l < phis.size(); ++l)
        {
            // an analytic level set is read on the slice as on any cell of dimension 2 in 3D
            if (phis[l].is_analytic())
            {
                sources.push_back(analytic_source(*phis[l].analytic, slice.box));
                continue;
            }
            restrict_affine(*phis[l].bernstein, std::span<const T>(slice.origin), std::span<const T>(slice.matrix), 2,
                            restricted[l], work);
            sources.push_back(bernstein_source(restricted[l]));
        }
        // a part that is small in the slice gets its own rectangle as the engine's cell: tolerances
        // are relative to the cell, and a small part in a large one would be coarse
        std::array<T, 2> lo, hi;
        const int fit = fit_to_part(slice.box, std::span<const Source<T>>(sources), terms, zero_sets, lo, hi);
        if (fit == 0)
            continue;
        points.points.clear();
        points.weights.clear();
        points.tags.clear();
        if (fit == 1)
            integrate(slice.box, std::span<const Source<T>>(sources), terms, q, opt, points, stats);
        else
        {
            fitted_box = sub_box(slice.box, Vec3<T>{lo[0], lo[1], T(0)}, Vec3<T>{hi[0], hi[1], T(1)});
            fitted.resize(phis.size());
            fitted_sources.clear();
            for (std::size_t l = 0; l < phis.size(); ++l)
            {
                if (phis[l].is_analytic())
                {
                    fitted_sources.push_back(analytic_source(*phis[l].analytic, fitted_box));
                    continue;
                }
                subdivide(restricted[l], std::span<const T>(lo), std::span<const T>(hi), fitted[l], work);
                fitted_sources.push_back(bernstein_source(fitted[l]));
            }
            integrate(fitted_box, std::span<const Source<T>>(fitted_sources), terms, q, opt, points, stats);
        }
        for (int i = 0; i < points.n_points(); ++i)
        {
            // the slice's box coordinates of the point
            const T s0 = lo[0] + (hi[0] - lo[0]) * points.points[3 * i], s1 = lo[1] + (hi[1] - lo[1]) * points.points[3 * i + 1];
            for (int j = 0; j < 3; ++j)
                out.points.push_back(slice.origin[j] + slice.matrix[2 * j] * s0 + slice.matrix[2 * j + 1] * s1);
            T factor = sweep.kind == SweepKind::turning ? slice.radius_origin + slice.radius_step * s0 : T(1);
            if (zero[0] >= 0)
            {
                const std::size_t first = out.points.size() - 3;
                const Vec3<T> u = {out.points[first], out.points[first + 1], out.points[first + 2]};
                const Vec3<T> normal = {slice.box.jacobian[0][2], slice.box.jacobian[1][2], slice.box.jacobian[2][2]};
                const Vec3<T> da = physical_gradient(zero[0], u);
                T cosine; // of the angle the slices are crossed at
                if (zero[1] < 0)
                {
                    // a surface: its slice is a curve, whose length becomes the surface's measure by
                    // the angle between the zero set's normal and the slice's
                    const Vec3<T> across = cross(da, normal);
                    cosine = std::sqrt(dot(across, across) / dot(da, da));
                }
                else
                {
                    // a curve: its slice is points, which become its length by the angle between
                    // its tangent and the slice's normal
                    const Vec3<T> tangent = cross(da, physical_gradient(zero[1], u));
                    cosine = std::abs(dot(tangent, normal)) / std::sqrt(dot(tangent, tangent));
                }
                factor /= std::max(cosine, std::numeric_limits<T>::min());
            }
            out.weights.push_back(points.weights[i] * weights[k] * factor);
        }
    }
}

// ============================================================================
// Events of a sweep along a box axis
// ============================================================================

namespace
{

template <std::floating_point T>
struct Range
{
    T lo = 0, hi = 0;
};

/// A level set on the box: for a Bernstein form its partial derivatives as
/// forms too, for an analytic level set the source, whose bounds over a box
/// come from affine_bounds. size: the functions' sizes, for rounding.
template <std::floating_point T>
struct Jet
{
    const Source<T>* source = nullptr;
    std::array<BoxBernstein<T>, 4> form; ///< value, d/du_0, d/du_1, d/du_2 (Bernstein forms only)
    std::array<T, 4> size = {};
};

/// A face of the cell in box coordinates: c . u <= d with |c| = 1.
template <std::floating_point T>
struct Face
{
    Vec3<T> c = {0, 0, 0};
    T d = 0;
};

/// The cell as the search reads it.
template <std::floating_point T>
struct Geometry
{
    std::vector<Jet<T>> jets;
    std::vector<Face<T>> faces;
    std::vector<std::array<int, 2>> edges; ///< the two faces an edge lies in
    std::vector<Vec3<T>> vertices;
};

enum class SystemKind : std::uint8_t
{
    critical,      ///< f_a = 0, d_p f_a = 0, d_q f_a = 0
    tangency,      ///< f_a = 0, f_b = 0, d_p f_a d_q f_b - d_q f_a d_p f_b = 0
    triple,        ///< f_a = f_b = f_c = 0
    face_tangency, ///< on face b: f_a = 0, d_p f_a c_q - d_q f_a c_p = 0 (c the face's normal)
    face_pair,     ///< on face c: f_a = f_b = 0
    edge_single    ///< on edge b: f_a = 0
};

struct System
{
    SystemKind kind;
    int a = 0, b = 0, c = 0;
};

template <std::floating_point T>
struct BoxRanges
{
    const Geometry<T>* cell = nullptr;
    Vec3<T> lo{}, hi{};
    std::vector<std::array<Range<T>, 4>> range;
    std::vector<std::array<bool, 4>> known;
    BoxBernstein<T> scratch;
    std::vector<T> work;

    void reset(const Geometry<T>& g, const Vec3<T>& l, const Vec3<T>& h)
    {
        cell = &g;
        lo = l;
        hi = h;
        range.assign(g.jets.size(), {});
        known.assign(g.jets.size(), {false, false, false, false});
    }
};

/// The range of function d (0: value, 1 + i: d/du_i) of level set l on the
/// box. Of an analytic level set the value's range is only as good as its
/// sign (all that is asked of it), and a bound that is lost is the whole line.
template <std::floating_point T>
const Range<T>& range_of(BoxRanges<T>& r, int l, int d)
{
    const std::size_t ls = static_cast<std::size_t>(l), ds = static_cast<std::size_t>(d);
    if (r.known[ls][ds])
        return r.range[ls][ds];
    const Jet<T>& jet = r.cell->jets[ls];
    const std::array<T, 9> matrix = {r.hi[0] - r.lo[0], 0, 0, 0, r.hi[1] - r.lo[1], 0, 0, 0, r.hi[2] - r.lo[2]};
    if (!jet.source->is_analytic())
    {
        restrict_affine(jet.form[ds], std::span<const T>(r.lo), std::span<const T>(matrix), 3, r.scratch, r.work);
        const auto [mn, mx] = std::minmax_element(r.scratch.coeffs.begin(), r.scratch.coeffs.end());
        r.range[ls][ds] = {*mn, *mx};
        r.known[ls][ds] = true;
        return r.range[ls][ds];
    }
    const T all = std::numeric_limits<T>::infinity(), far = std::numeric_limits<T>::max();
    AffineBounds<T> bounds;
    std::array<Range<T>, 4>& out = r.range[ls];
    out = {Range<T>{-all, all}, Range<T>{-all, all}, Range<T>{-all, all}, Range<T>{-all, all}};
    if (affine_bounds(*jet.source, std::span<const T>(r.lo), std::span<const T>(matrix), 3, bounds))
    {
        out[0] = bounds.sign > 0 ? Range<T>{far, far}
                                 : (bounds.sign < 0 ? Range<T>{-far, -far} : Range<T>{-bounds.magnitude, bounds.magnitude});
        // the bounds are of the derivatives along the box's edges
        for (std::size_t j = 0; j < 3 && bounds.has_derivatives; ++j)
            out[1 + j] = {bounds.dmin[j] / matrix[4 * j], bounds.dmax[j] / matrix[4 * j]};
    }
    r.known[ls] = {true, true, true, true};
    return out[ds];
}

/// True if both ends of the range are numbers: a bound that holds.
template <std::floating_point T>
bool bounded(const Range<T>& x)
{
    return std::isfinite(x.lo) && std::isfinite(x.hi);
}

template <std::floating_point T>
constexpr T rounding = T(256) * std::numeric_limits<T>::epsilon();

/// True if function d of level set l may vanish on the box.
template <std::floating_point T>
bool may_vanish_on(BoxRanges<T>& r, int l, int d)
{
    const Range<T>& x = range_of(r, l, d);
    const T eps = rounding<T> * r.cell->jets[static_cast<std::size_t>(l)].size[static_cast<std::size_t>(d)];
    return x.lo <= eps && x.hi >= -eps;
}

/// The range of c . u - d of face j on the box.
template <std::floating_point T>
Range<T> face_range(const BoxRanges<T>& r, int j)
{
    const Face<T>& f = r.cell->faces[static_cast<std::size_t>(j)];
    Range<T> x = {-f.d, -f.d};
    for (std::size_t i = 0; i < 3; ++i)
    {
        const T a = f.c[i] * r.lo[i], b = f.c[i] * r.hi[i];
        x.lo += std::min(a, b);
        x.hi += std::max(a, b);
    }
    return x;
}

template <std::floating_point T>
bool may_touch_face(const BoxRanges<T>& r, int j)
{
    const Range<T> x = face_range(r, j);
    return x.lo <= rounding<T> && x.hi >= -rounding<T>;
}

/// True if the signs the term asks of the level sets may hold where they
/// take the values @p value (a range per level set).
template <std::floating_point T, typename Value>
bool term_may_hold(const SelectionTerm& term, const std::vector<Jet<T>>& jets, const Value& value)
{
    for (int l = 0; l < static_cast<int>(jets.size()); ++l)
    {
        const std::uint64_t bit = std::uint64_t(1) << l;
        if (!((term.negative_required | term.positive_required) & bit))
            continue;
        const Range<T> x = value(l);
        const T eps = rounding<T> * jets[static_cast<std::size_t>(l)].size[0];
        if (((term.negative_required & bit) && x.lo > eps) || ((term.positive_required & bit) && x.hi < -eps))
            return false;
    }
    return true;
}

/// True if the box may meet the cell and, in it, the closure of the part the
/// terms select.
template <std::floating_point T>
bool may_meet_part(BoxRanges<T>& r, std::span<const SelectionTerm> terms)
{
    for (int j = 0; j < static_cast<int>(r.cell->faces.size()); ++j)
        if (face_range(r, j).lo > rounding<T>)
            return false;
    for (const SelectionTerm& term : terms)
        if (term_may_hold(term, r.cell->jets, [&](int l) { return range_of(r, l, 0); }))
            return true;
    return false;
}

/// True if the system may have a solution on the box; p and q are the box
/// axes other than the sweep's.
template <std::floating_point T>
bool may_hold(BoxRanges<T>& r, const System& s, int p, int q)
{
    const std::vector<Jet<T>>& jets = r.cell->jets;
    const std::size_t sp = static_cast<std::size_t>(1 + p), sq = static_cast<std::size_t>(1 + q);
    switch (s.kind)
    {
    case SystemKind::critical:
        return may_vanish_on(r, s.a, 0) && may_vanish_on(r, s.a, 1 + p) && may_vanish_on(r, s.a, 1 + q);
    case SystemKind::tangency:
    {
        if (!may_vanish_on(r, s.a, 0) || !may_vanish_on(r, s.b, 0))
            return false;
        // interval of d_p f_a d_q f_b - d_q f_a d_p f_b
        const auto product = [](const Range<T>& x, const Range<T>& y)
        {
            const T v[4] = {x.lo * y.lo, x.lo * y.hi, x.hi * y.lo, x.hi * y.hi};
            return Range<T>{*std::min_element(v, v + 4), *std::max_element(v, v + 4)};
        };
        if (!bounded(range_of(r, s.a, 1 + p)) || !bounded(range_of(r, s.a, 1 + q))
            || !bounded(range_of(r, s.b, 1 + p)) || !bounded(range_of(r, s.b, 1 + q)))
            return true;
        const Range<T> first = product(range_of(r, s.a, 1 + p), range_of(r, s.b, 1 + q));
        const Range<T> second = product(range_of(r, s.a, 1 + q), range_of(r, s.b, 1 + p));
        const std::size_t a = static_cast<std::size_t>(s.a), b = static_cast<std::size_t>(s.b);
        const T eps = rounding<T> * (jets[a].size[sp] * jets[b].size[sq] + jets[a].size[sq] * jets[b].size[sp]);
        return first.lo - second.hi <= eps && first.hi - second.lo >= -eps;
    }
    case SystemKind::triple:
        return may_vanish_on(r, s.a, 0) && may_vanish_on(r, s.b, 0) && may_vanish_on(r, s.c, 0);
    case SystemKind::face_tangency:
    {
        if (!may_vanish_on(r, s.a, 0) || !may_touch_face(r, s.b))
            return false;
        // interval of d_p f_a c_q - d_q f_a c_p
        const Face<T>& face = r.cell->faces[static_cast<std::size_t>(s.b)];
        if (!bounded(range_of(r, s.a, 1 + p)) || !bounded(range_of(r, s.a, 1 + q)))
            return true;
        const auto scaled = [](const Range<T>& x, T c)
        { return c >= T(0) ? Range<T>{c * x.lo, c * x.hi} : Range<T>{c * x.hi, c * x.lo}; };
        const Range<T> first = scaled(range_of(r, s.a, 1 + p), face.c[static_cast<std::size_t>(q)]);
        const Range<T> second = scaled(range_of(r, s.a, 1 + q), face.c[static_cast<std::size_t>(p)]);
        const std::size_t a = static_cast<std::size_t>(s.a);
        const T eps = rounding<T> * (jets[a].size[sp] + jets[a].size[sq]);
        return first.lo - second.hi <= eps && first.hi - second.lo >= -eps;
    }
    case SystemKind::face_pair:
        return may_vanish_on(r, s.a, 0) && may_vanish_on(r, s.b, 0) && may_touch_face(r, s.c);
    case SystemKind::edge_single:
    {
        const std::array<int, 2>& edge = r.cell->edges[static_cast<std::size_t>(s.b)];
        return may_vanish_on(r, s.a, 0) && may_touch_face(r, edge[0]) && may_touch_face(r, edge[1]);
    }
    }
    return false;
}

/// A box the search kept: its range along the sweep's axis, its centre, and
/// the first system that may hold on it (-1: a vertex of the cell).
template <std::floating_point T>
struct Hit
{
    T lo = 0, hi = 0;
    Vec3<T> centre = {0, 0, 0};
    int system = -1;
};

/// Bisect the unit box and keep the boxes on which some system may hold and
/// which may meet the part, down to @p width along @p axis; they go to
/// @p found.
/// @return false if more than @p budget boxes were visited
template <std::floating_point T>
bool search(const Geometry<T>& cell, std::span<const SelectionTerm> terms, const std::vector<System>& systems,
            int axis, int p, int q, T width, std::vector<Hit<T>>& found, int& budget)
{
    struct Box
    {
        Vec3<T> lo, hi;
    };
    thread_local std::vector<Box> stack;
    thread_local BoxRanges<T> ranges;
    stack.clear();
    stack.push_back({Vec3<T>{0, 0, 0}, Vec3<T>{1, 1, 1}});
    const std::size_t k = static_cast<std::size_t>(axis);
    while (!stack.empty())
    {
        const Box box = stack.back();
        stack.pop_back();
        if (--budget < 0)
            return false;
        ranges.reset(cell, box.lo, box.hi);
        if (!may_meet_part(ranges, terms))
            continue;
        int held = -1;
        for (std::size_t i = 0; i < systems.size() && held < 0; ++i)
            if (may_hold(ranges, systems[i], p, q))
                held = static_cast<int>(i);
        if (held < 0)
            continue;
        if (box.hi[k] - box.lo[k] <= width)
        {
            found.push_back({box.lo[k], box.hi[k],
                             Vec3<T>{T(0.5) * (box.lo[0] + box.hi[0]), T(0.5) * (box.lo[1] + box.hi[1]),
                                     T(0.5) * (box.lo[2] + box.hi[2])},
                             held});
            continue;
        }
        std::size_t widest = 0;
        for (std::size_t j = 1; j < 3; ++j)
            if (box.hi[j] - box.lo[j] > box.hi[widest] - box.lo[widest])
                widest = j;
        Box left = box, right = box;
        left.hi[widest] = right.lo[widest] = T(0.5) * (box.lo[widest] + box.hi[widest]);
        stack.push_back(left);
        stack.push_back(right);
    }
    return true;
}

/// The level sets with their derivatives, and the faces, edges and vertices of
/// the cell (the clipped box's polytope; faces in one plane count once).
template <std::floating_point T>
void make_geometry(const ClippedBox<T>& cell, std::span<const Source<T>> phis, Geometry<T>& g)
{
    g.jets.resize(phis.size());
    for (std::size_t l = 0; l < phis.size(); ++l)
    {
        Jet<T>& jet = g.jets[l];
        jet.source = &phis[l];
        if (phis[l].is_analytic())
        {
            jet.size.fill(reference_magnitude(phis[l]));
            continue;
        }
        if (phis[l].bernstein->dim != 3)
            throw std::invalid_argument("quadrays: the events of a sweep need level sets on a box of dimension 3");
        jet.form[0] = *phis[l].bernstein;
        for (int d = 0; d < 3; ++d)
            derivative(jet.form[0], d, jet.form[static_cast<std::size_t>(1 + d)]);
        for (std::size_t d = 0; d < 4; ++d)
            jet.size[d] = max_abs(std::span<const T>(jet.form[d].coeffs));
    }
    thread_local Polytope<T> poly;
    clipped_polytope(cell, poly);
    g.vertices = poly.vertices;
    g.faces.clear();
    g.edges.clear();
    const T same = scaled_tolerance<T>(1e-9);
    for (int f = 0; f < poly.n_faces(); ++f)
    {
        const Vec3<T>& c = poly.face_normals[static_cast<std::size_t>(f)];
        const int e = poly.face_edges[static_cast<std::size_t>(poly.face_offsets[static_cast<std::size_t>(f)])];
        const T d = dot(c, poly.vertices[static_cast<std::size_t>(poly.edges[static_cast<std::size_t>(e)][0])]);
        bool repeated = false;
        for (const Face<T>& other : g.faces)
            repeated = repeated || (std::abs(dot(c, other.c) - T(1)) < same && std::abs(d - other.d) < same);
        if (!repeated)
            g.faces.push_back({c, d});
    }
    const auto lies_on = [&](const Vec3<T>& v, const Face<T>& f) { return std::abs(dot(f.c, v) - f.d) < same; };
    for (int i = 0; i < static_cast<int>(g.faces.size()); ++i)
        for (int j = i + 1; j < static_cast<int>(g.faces.size()); ++j)
        {
            int shared = 0;
            for (const Vec3<T>& v : g.vertices)
                shared += lies_on(v, g.faces[static_cast<std::size_t>(i)]) && lies_on(v, g.faces[static_cast<std::size_t>(j)]);
            if (shared >= 2)
                g.edges.push_back({i, j});
        }
}

/// The solution of a system next to @p u, by Newton's method with a Jacobian
/// from differences: the search finds an event's position to its boxes' width,
/// and this to rounding.
/// @return false if the iteration does not settle within @p reach of its start
template <std::floating_point T>
bool polish(const Geometry<T>& cell, const System& s, int p, int q, T reach, Vec3<T>& u)
{
    const auto value = [&](int l, const Vec3<T>& x) { return evaluate(*cell.jets[static_cast<std::size_t>(l)].source, std::span<const T>(x)); };
    const auto slope = [&](int l, const Vec3<T>& x)
    {
        Vec3<T> g;
        gradient(*cell.jets[static_cast<std::size_t>(l)].source, std::span<const T>(x), std::span<T>(g));
        return g;
    };
    const auto face = [&](int j, const Vec3<T>& x)
    {
        const Face<T>& f = cell.faces[static_cast<std::size_t>(j)];
        return dot(f.c, x) - f.d;
    };
    const std::size_t sp = static_cast<std::size_t>(p), sq = static_cast<std::size_t>(q);
    const auto residual = [&](const Vec3<T>& x) -> Vec3<T>
    {
        switch (s.kind)
        {
        case SystemKind::critical:
        {
            const Vec3<T> g = slope(s.a, x);
            return {value(s.a, x), g[sp], g[sq]};
        }
        case SystemKind::tangency:
        {
            const Vec3<T> ga = slope(s.a, x), gb = slope(s.b, x);
            return {value(s.a, x), value(s.b, x), ga[sp] * gb[sq] - ga[sq] * gb[sp]};
        }
        case SystemKind::triple:
            return {value(s.a, x), value(s.b, x), value(s.c, x)};
        case SystemKind::face_tangency:
        {
            const Vec3<T> g = slope(s.a, x);
            const Face<T>& f = cell.faces[static_cast<std::size_t>(s.b)];
            return {value(s.a, x), face(s.b, x), g[sp] * f.c[sq] - g[sq] * f.c[sp]};
        }
        case SystemKind::face_pair:
            return {value(s.a, x), value(s.b, x), face(s.c, x)};
        case SystemKind::edge_single:
        {
            const std::array<int, 2>& edge = cell.edges[static_cast<std::size_t>(s.b)];
            return {value(s.a, x), face(edge[0], x), face(edge[1], x)};
        }
        }
        return {0, 0, 0};
    };
    const Vec3<T> start = u;
    const T h = std::cbrt(std::numeric_limits<T>::epsilon());
    for (int it = 0; it < 30; ++it)
    {
        const Vec3<T> r = residual(u);
        Mat3<T> jac; // jac[i][j] = d r_i / d u_j
        for (std::size_t j = 0; j < 3; ++j)
        {
            Vec3<T> up = u, down = u;
            up[j] += h;
            down[j] -= h;
            const Vec3<T> a = residual(up), b = residual(down);
            for (std::size_t i = 0; i < 3; ++i)
                jac[i][j] = (a[i] - b[i]) / (T(2) * h);
        }
        const T det = dot(jac[0], cross(jac[1], jac[2]));
        if (!(std::abs(det) > T(0)) || !std::isfinite(det))
            return false;
        // Cramer's rule for jac step = r
        Vec3<T> step;
        for (std::size_t j = 0; j < 3; ++j)
        {
            Mat3<T> m = jac;
            for (std::size_t i = 0; i < 3; ++i)
                m[i][j] = r[i];
            step[j] = dot(m[0], cross(m[1], m[2])) / det;
        }
        T moved = 0, away = 0;
        for (std::size_t j = 0; j < 3; ++j)
        {
            u[j] -= step[j];
            moved = std::max(moved, std::abs(step[j]));
            away = std::max(away, std::abs(u[j] - start[j]));
        }
        if (!(away <= reach))
            return false;
        if (moved <= T(8) * std::numeric_limits<T>::epsilon())
            return true;
    }
    return false;
}

} // namespace

template <std::floating_point T>
bool sweep_events(const ClippedBox<T>& cell, std::span<const Source<T>> phis,
                  std::span<const SelectionTerm> terms, int axis, std::vector<T>& events, int max_boxes)
{
    if (axis < 0 || axis > 2)
        throw std::invalid_argument("quadrays: a sweep runs along box axis 0, 1 or 2");
    if (cell.tdim != 3)
        throw std::invalid_argument("quadrays: slices are of cells of dimension 3");
    const std::uint64_t zero = zero_sets_of(terms);
    thread_local Geometry<T> g;
    thread_local std::vector<System> systems;
    thread_local std::vector<Hit<T>> found;
    make_geometry(cell, phis, g);
    const int n = static_cast<int>(phis.size());
    const int p = (axis + 1) % 3, q = (axis + 2) % 3;
    const T width = scaled_tolerance<T>(1e-7);

    // the systems that hold all the zero sets the part asks for: every one for a volume
    systems.clear();
    const auto bit = [](int a) { return std::uint64_t(1) << a; };
    const auto holds = [&](std::uint64_t sets) { return (zero & ~sets) == 0; };
    const int n_faces = static_cast<int>(g.faces.size());
    for (int a = 0; a < n; ++a)
    {
        if (holds(bit(a)))
            systems.push_back({SystemKind::critical, a, 0, 0});
        for (int b = a + 1; b < n; ++b)
        {
            if (holds(bit(a) | bit(b)))
            {
                systems.push_back({SystemKind::tangency, a, b, 0});
                for (int f = 0; f < n_faces; ++f)
                    systems.push_back({SystemKind::face_pair, a, b, f});
            }
            for (int c = b + 1; c < n; ++c)
                if (holds(bit(a) | bit(b) | bit(c)))
                    systems.push_back({SystemKind::triple, a, b, c});
        }
        if (!holds(bit(a)))
            continue;
        for (int f = 0; f < n_faces; ++f)
        {
            // on a face that lies in a slice the whole curve is at the face's position: its vertices give it
            const Face<T>& face = g.faces[static_cast<std::size_t>(f)];
            if (std::abs(face.c[static_cast<std::size_t>(p)]) + std::abs(face.c[static_cast<std::size_t>(q)])
                > scaled_tolerance<T>(1e-9))
                systems.push_back({SystemKind::face_tangency, a, f, 0});
        }
        for (int e = 0; e < static_cast<int>(g.edges.size()); ++e)
            systems.push_back({SystemKind::edge_single, a, e, 0});
    }
    found.clear();
    int budget = max_boxes;
    if (!search(g, terms, systems, axis, p, q, width, found, budget))
        return false;
    // an event's position enters a surface's or a curve's measure to first order or worse (where a
    // slice first touches a surface, its slice integral jumps): Newton's method takes it from the
    // boxes' width to rounding
    const std::size_t k = static_cast<std::size_t>(axis);
    for (Hit<T>& hit : found)
    {
        Vec3<T> u = hit.centre;
        if (polish(g, systems[static_cast<std::size_t>(hit.system)], p, q, T(100) * width, u))
            hit.lo = hit.hi = u[k];
    }

    // the cell's vertices on the part: the pattern of the slice polygon changes there
    thread_local std::vector<T> values;
    for (const Vec3<T>& v : g.vertices)
    {
        values.resize(phis.size());
        for (std::size_t l = 0; l < phis.size(); ++l)
            values[l] = evaluate(phis[l], std::span<const T>(v));
        const auto at = [&](int l) { return Range<T>{values[static_cast<std::size_t>(l)], values[static_cast<std::size_t>(l)]}; };
        bool on_part = false;
        for (const SelectionTerm& term : terms)
            on_part = on_part || term_may_hold(term, g.jets, at);
        for (const int l : level_sets_of(zero)) // a vertex of the cell on the part's zero sets
            if (l >= 0)
                on_part = on_part && std::abs(values[static_cast<std::size_t>(l)])
                                         <= scaled_tolerance<T>(1e-9) * g.jets[static_cast<std::size_t>(l)].size[0];
        if (on_part)
            found.push_back({v[k], v[k], v, -1});
    }

    // boxes that touch are one event; a longer run of them is kept as its two ends
    std::sort(found.begin(), found.end(), [](const Hit<T>& x, const Hit<T>& y) { return x.lo < y.lo; });
    events.clear();
    const T together = scaled_tolerance<T>(1e-5), margin = scaled_tolerance<T>(1e-6);
    const auto add = [&](T t)
    {
        if (t > margin && t < T(1) - margin && (events.empty() || t - events.back() > margin))
            events.push_back(t);
    };
    for (std::size_t i = 0; i < found.size();)
    {
        T lo = found[i].lo, hi = found[i].hi;
        std::size_t j = i + 1;
        for (; j < found.size() && found[j].lo <= hi + width; ++j)
            hi = std::max(hi, found[j].hi);
        if (hi - lo <= together)
            add(T(0.5) * (lo + hi));
        else
        {
            add(lo);
            add(hi);
        }
        i = j;
    }
    return true;
}

template <std::floating_point T>
void piece_rule(T a, bool flat_a, T b, bool flat_b, int m, std::vector<T>& positions, std::vector<T>& weights)
{
    thread_local std::vector<T> x, w;
    gauss_legendre(m, x, w);
    const T length = b - a;
    for (int i = 0; i < m; ++i)
    {
        const T u = x[static_cast<std::size_t>(i)], wu = w[static_cast<std::size_t>(i)];
        if (flat_a)
        {
            positions.push_back(a + length * u * u);
            weights.push_back(wu * T(2) * length * u);
        }
        else if (flat_b)
        {
            positions.push_back(b - length * (T(1) - u) * (T(1) - u));
            weights.push_back(wu * T(2) * length * (T(1) - u));
        }
        else
        {
            positions.push_back(a + length * u);
            weights.push_back(wu * length);
        }
    }
}

template <std::floating_point T>
bool integrate_by_slices(const ClippedBox<T>& cell, std::span<const Source<T>> phis,
                         std::span<const SelectionTerm> terms, int q, const Options& opt, CellPoints<T>& out,
                         Stats& stats, unsigned axes)
{
    if (cell.tdim != 3 || phis.empty() || terms.empty())
        return false;
    for (const Source<T>& phi : phis)
        if (!phi.is_analytic() && phi.bernstein->dim != 3)
            return false;
    const std::uint64_t zero = terms[0].zero_required;
    if (std::popcount(zero) > 2)
        return false;
    for (const SelectionTerm& term : terms)
        if (term.zero_required != zero)
            return false;
    thread_local std::vector<T> events, best, positions, weights;
    int axis = -1;
    for (int k = 0; k < 3; ++k)
        if (((axes >> k) & 1U) && sweep_events(cell, phis, terms, k, events)
            && (axis < 0 || events.size() < best.size()))
        {
            axis = k;
            best = events;
        }
    if (axis < 0)
        return false;

    // the planes u_axis = const in physical space: their normal is row axis of the inverse map
    const Mat3<T> inverse = inverse_jacobian(cell);
    const Vec3<T> normal = inverse[static_cast<std::size_t>(axis)];
    const T step = T(1) / std::sqrt(dot(normal, normal));
    Sweep<T> sweep;
    sweep.kind = SweepKind::parallel;
    sweep.origin = cell.origin;
    sweep.direction = normal;

    // the pieces of the sweep between events; one between two events is halved
    struct Piece
    {
        T a, b;
        bool flat_a, flat_b;
        CellPoints<T> points;
        std::array<T, 10> moments;
    };
    thread_local std::vector<Piece> pieces;
    pieces.clear();
    T a = 0;
    for (std::size_t i = 0; i <= best.size(); ++i)
    {
        const T b = i < best.size() ? best[i] : T(1);
        const bool flat_a = i > 0, flat_b = i < best.size();
        if (flat_a && flat_b)
        {
            const T mid = T(0.5) * (a + b);
            pieces.push_back({a, mid, true, false, {}, {}});
            pieces.push_back({mid, b, false, true, {}, {}});
        }
        else
            pieces.push_back({a, b, flat_a, flat_b, {}, {}});
        a = b;
    }
    // the rule of a piece with m slices, and its integrals of 1, u_i and u_i u_j
    const auto run = [&](Piece& piece, int m)
    {
        positions.clear();
        weights.clear();
        piece_rule(piece.a, piece.flat_a, piece.b, piece.flat_b, m, positions, weights);
        for (std::size_t i = 0; i < positions.size(); ++i)
        {
            positions[i] *= step;
            weights[i] *= step;
        }
        piece.points.points.clear();
        piece.points.weights.clear();
        integrate_slices(cell, phis, terms, sweep, std::span<const T>(positions), std::span<const T>(weights), q, opt,
                         piece.points, stats);
        piece.moments.fill(T(0));
        for (int i = 0; i < piece.points.n_points(); ++i)
        {
            const T* u = &piece.points.points[static_cast<std::size_t>(3 * i)];
            const T w = piece.points.weights[static_cast<std::size_t>(i)];
            piece.moments[0] += w;
            int k = 4;
            for (int c = 0; c < 3; ++c)
            {
                piece.moments[static_cast<std::size_t>(1 + c)] += w * u[c];
                for (int d = c; d < 3; ++d)
                    piece.moments[static_cast<std::size_t>(k++)] += w * u[c] * u[d];
            }
        }
    };
    const int first_flat = q, first_plain = std::max(3, q - 2), most = 64;
    T volume = 0;
    for (Piece& piece : pieces)
    {
        run(piece, piece.flat_a || piece.flat_b ? first_flat : first_plain);
        volume += piece.moments[0];
    }
    if (!(volume > T(0))) // nothing found; the measure is an area for a surface, a length for a curve
        return false;
    const T tolerance = std::max(scaled_tolerance<T>(1e-13), std::pow(T(10), -T(q + 4))) * volume;
    for (Piece& piece : pieces)
    {
        int m = piece.flat_a || piece.flat_b ? first_flat : first_plain;
        for (;;)
        {
            if (m >= most)
                return false;
            const std::array<T, 10> before = piece.moments;
            m = std::min(most, m + m / 2 + 1);
            run(piece, m);
            T moved = 0;
            for (std::size_t i = 0; i < 10; ++i)
                moved = std::max(moved, std::abs(piece.moments[i] - before[i]));
            if (moved <= tolerance)
                break;
        }
    }
    for (const Piece& piece : pieces)
    {
        out.points.insert(out.points.end(), piece.points.points.begin(), piece.points.points.end());
        out.weights.insert(out.weights.end(), piece.points.weights.begin(), piece.points.weights.end());
    }
    return true;
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template bool sweep_events<float>(const ClippedBox<float>&, std::span<const Source<float>>,
                                  std::span<const SelectionTerm>, int, std::vector<float>&, int);
template bool sweep_events<double>(const ClippedBox<double>&, std::span<const Source<double>>,
                                   std::span<const SelectionTerm>, int, std::vector<double>&, int);
template void piece_rule<float>(float, bool, float, bool, int, std::vector<float>&, std::vector<float>&);
template void piece_rule<double>(double, bool, double, bool, int, std::vector<double>&, std::vector<double>&);
template bool integrate_by_slices<float>(const ClippedBox<float>&, std::span<const Source<float>>,
                                         std::span<const SelectionTerm>, int, const Options&, CellPoints<float>&,
                                         Stats&, unsigned);
template bool integrate_by_slices<double>(const ClippedBox<double>&, std::span<const Source<double>>,
                                          std::span<const SelectionTerm>, int, const Options&,
                                          CellPoints<double>&, Stats&, unsigned);

template bool make_slice<float>(const ClippedBox<float>&, const Sweep<float>&, float, Slice<float>&);
template bool make_slice<double>(const ClippedBox<double>&, const Sweep<double>&, double, Slice<double>&);
template void integrate_slices<float>(const ClippedBox<float>&, std::span<const Source<float>>,
                                      std::span<const SelectionTerm>, const Sweep<float>&, std::span<const float>,
                                      std::span<const float>, int, const Options&, CellPoints<float>&, Stats&);
template void integrate_slices<double>(const ClippedBox<double>&, std::span<const Source<double>>,
                                       std::span<const SelectionTerm>, const Sweep<double>&,
                                       std::span<const double>, std::span<const double>, int, const Options&,
                                       CellPoints<double>&, Stats&);

} // namespace cutcells::quadrays
