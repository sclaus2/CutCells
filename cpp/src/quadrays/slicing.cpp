// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "slicing.h"

#include <algorithm>
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
void integrate_slices(const ClippedBox<T>& cell, std::span<const BoxBernstein<T>> forms,
                      std::span<const SelectionTerm> terms, const Sweep<T>& sweep, std::span<const T> positions,
                      std::span<const T> weights, int q, const Options& opt, CellPoints<T>& out, Stats& stats)
{
    if (positions.size() != weights.size())
        throw std::invalid_argument("quadrays: slices need one weight per position");
    for (const SelectionTerm& term : terms)
        if (term.zero_required != 0)
            throw std::invalid_argument("quadrays: slices integrate volume parts");
    thread_local Slice<T> slice;
    thread_local std::vector<BoxBernstein<T>> restricted;
    thread_local std::vector<Source<T>> sources;
    thread_local CellPoints<T> points;
    thread_local std::vector<T> work;
    restricted.resize(forms.size());
    for (std::size_t k = 0; k < positions.size(); ++k)
    {
        if (!make_slice(cell, sweep, positions[k], slice))
            continue;
        sources.clear();
        for (std::size_t l = 0; l < forms.size(); ++l)
        {
            restrict_affine(forms[l], std::span<const T>(slice.origin), std::span<const T>(slice.matrix), 2,
                            restricted[l], work);
            sources.push_back(bernstein_source(restricted[l]));
        }
        points.points.clear();
        points.weights.clear();
        points.tags.clear();
        integrate(slice.box, std::span<const Source<T>>(sources), terms, q, opt, points, stats);
        for (int i = 0; i < points.n_points(); ++i)
        {
            const T s0 = points.points[3 * i], s1 = points.points[3 * i + 1];
            for (int j = 0; j < 3; ++j)
                out.points.push_back(slice.origin[j] + slice.matrix[2 * j] * s0 + slice.matrix[2 * j + 1] * s1);
            const T factor = sweep.kind == SweepKind::turning ? slice.radius_origin + slice.radius_step * s0 : T(1);
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

/// A level set and its partial derivatives on the box, with their sizes.
template <std::floating_point T>
struct Jet
{
    std::array<BoxBernstein<T>, 4> form; ///< value, d/du_0, d/du_1, d/du_2
    std::array<T, 4> size = {};
};

enum class SystemKind : std::uint8_t
{
    critical,      ///< f_a = 0, d_p f_a = 0, d_q f_a = 0
    tangency,      ///< f_a = 0, f_b = 0, d_p f_a d_q f_b - d_q f_a d_p f_b = 0
    triple,        ///< f_a = f_b = f_c = 0
    face_tangency, ///< on a lateral face: f_a = 0, d_p f_a = 0
    pair,          ///< f_a = f_b = 0
    single         ///< f_a = 0
};

struct System
{
    SystemKind kind;
    int a = 0, b = 0, c = 0;
};

template <std::floating_point T>
struct BoxRanges
{
    const std::vector<Jet<T>>* jets = nullptr;
    Vec3<T> lo{}, hi{};
    std::vector<std::array<Range<T>, 4>> range;
    std::vector<std::array<bool, 4>> known;
    BoxBernstein<T> scratch;
    std::vector<T> work;

    void reset(const std::vector<Jet<T>>& j, const Vec3<T>& l, const Vec3<T>& h)
    {
        jets = &j;
        lo = l;
        hi = h;
        range.assign(j.size(), {});
        known.assign(j.size(), {false, false, false, false});
    }
};

/// The range of function d (0: value, 1 + i: d/du_i) of level set l on the box.
template <std::floating_point T>
const Range<T>& range_of(BoxRanges<T>& r, int l, int d)
{
    const std::size_t ls = static_cast<std::size_t>(l), ds = static_cast<std::size_t>(d);
    if (!r.known[ls][ds])
    {
        // the form on the box; a dimension of width 0 (a face, an edge) drops out
        const std::array<T, 9> matrix = {r.hi[0] - r.lo[0], 0, 0, 0, r.hi[1] - r.lo[1], 0, 0, 0, r.hi[2] - r.lo[2]};
        restrict_affine((*r.jets)[ls].form[ds], std::span<const T>(r.lo), std::span<const T>(matrix), 3, r.scratch,
                        r.work);
        const auto [mn, mx] = std::minmax_element(r.scratch.coeffs.begin(), r.scratch.coeffs.end());
        r.range[ls][ds] = {*mn, *mx};
        r.known[ls][ds] = true;
    }
    return r.range[ls][ds];
}

template <std::floating_point T>
constexpr T rounding = T(256) * std::numeric_limits<T>::epsilon();

/// True if function d of level set l may vanish on the box.
template <std::floating_point T>
bool may_vanish_on(BoxRanges<T>& r, int l, int d)
{
    const Range<T>& x = range_of(r, l, d);
    const T eps = rounding<T> * (*r.jets)[static_cast<std::size_t>(l)].size[static_cast<std::size_t>(d)];
    return x.lo <= eps && x.hi >= -eps;
}

/// True if the box may meet the closure of the part the terms select.
template <std::floating_point T>
bool may_meet_part(BoxRanges<T>& r, std::span<const SelectionTerm> terms)
{
    const int n = static_cast<int>(r.jets->size());
    for (const SelectionTerm& term : terms)
    {
        bool possible = true;
        for (int l = 0; l < n && possible; ++l)
        {
            const std::uint64_t bit = std::uint64_t(1) << l;
            if (!((term.negative_required | term.positive_required) & bit))
                continue;
            const Range<T>& x = range_of(r, l, 0);
            const T eps = rounding<T> * (*r.jets)[static_cast<std::size_t>(l)].size[0];
            if ((term.negative_required & bit) && x.lo > eps)
                possible = false;
            if ((term.positive_required & bit) && x.hi < -eps)
                possible = false;
        }
        if (possible)
            return true;
    }
    return false;
}

/// True if the system may have a solution on the box; p and q are the box
/// axes other than the sweep's (on a lateral face, p is the free one).
template <std::floating_point T>
bool may_hold(BoxRanges<T>& r, const System& s, int p, int q)
{
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
        const Range<T> first = product(range_of(r, s.a, 1 + p), range_of(r, s.b, 1 + q));
        const Range<T> second = product(range_of(r, s.a, 1 + q), range_of(r, s.b, 1 + p));
        const std::vector<Jet<T>>& jets = *r.jets;
        const std::size_t a = static_cast<std::size_t>(s.a), b = static_cast<std::size_t>(s.b);
        const std::size_t sp = static_cast<std::size_t>(1 + p), sq = static_cast<std::size_t>(1 + q);
        const T eps = rounding<T> * (jets[a].size[sp] * jets[b].size[sq] + jets[a].size[sq] * jets[b].size[sp]);
        return first.lo - second.hi <= eps && first.hi - second.lo >= -eps;
    }
    case SystemKind::triple:
        return may_vanish_on(r, s.a, 0) && may_vanish_on(r, s.b, 0) && may_vanish_on(r, s.c, 0);
    case SystemKind::face_tangency:
        return may_vanish_on(r, s.a, 0) && may_vanish_on(r, s.a, 1 + p);
    case SystemKind::pair:
        return may_vanish_on(r, s.a, 0) && may_vanish_on(r, s.b, 0);
    case SystemKind::single:
        return may_vanish_on(r, s.a, 0);
    }
    return false;
}

/// Bisect the box [lo, hi] (dimensions of width 0 stay) and keep the boxes on
/// which some system may hold and which may meet the part, down to @p width
/// along @p axis, or to @p size in every direction if @p axis is -1; the
/// ranges of the boxes kept along @p axis go to @p found.
/// @return false if more than @p budget boxes were visited
template <std::floating_point T>
bool search(const std::vector<Jet<T>>& jets, std::span<const SelectionTerm> terms,
            const std::vector<System>& systems, const Vec3<T>& lo, const Vec3<T>& hi, int axis, int p, int q,
            T width, std::vector<Range<T>>& found, int& budget)
{
    struct Box
    {
        Vec3<T> lo, hi;
    };
    thread_local std::vector<Box> stack;
    thread_local BoxRanges<T> ranges;
    stack.clear();
    stack.push_back({lo, hi});
    while (!stack.empty())
    {
        const Box box = stack.back();
        stack.pop_back();
        if (--budget < 0)
            return false;
        ranges.reset(jets, box.lo, box.hi);
        if (!may_meet_part(ranges, terms))
            continue;
        bool any = false;
        for (const System& s : systems)
            if (may_hold(ranges, s, p, q))
            {
                any = true;
                break;
            }
        if (!any)
            continue;
        int widest = 0;
        for (int j = 1; j < 3; ++j)
            if (box.hi[j] - box.lo[j] > box.hi[widest] - box.lo[widest])
                widest = j;
        const T extent = axis >= 0 ? box.hi[axis] - box.lo[axis] : box.hi[widest] - box.lo[widest];
        if (extent <= width)
        {
            found.push_back(axis >= 0 ? Range<T>{box.lo[axis], box.hi[axis]} : Range<T>{T(0), T(1)});
            if (axis < 0)
                return true; // one box is enough
            continue;
        }
        Box left = box, right = box;
        left.hi[widest] = right.lo[widest] = T(0.5) * (box.lo[widest] + box.hi[widest]);
        stack.push_back(left);
        stack.push_back(right);
    }
    return true;
}

template <std::floating_point T>
void make_jets(std::span<const BoxBernstein<T>> forms, std::vector<Jet<T>>& jets)
{
    jets.resize(forms.size());
    for (std::size_t l = 0; l < forms.size(); ++l)
    {
        if (forms[l].dim != 3)
            throw std::invalid_argument("quadrays: the events of a sweep need level sets on a box of dimension 3");
        jets[l].form[0] = forms[l];
        for (int d = 0; d < 3; ++d)
            derivative(forms[l], d, jets[l].form[static_cast<std::size_t>(1 + d)]);
        for (std::size_t d = 0; d < 4; ++d)
            jets[l].size[d] = max_abs(std::span<const T>(jets[l].form[d].coeffs));
    }
}

} // namespace

template <std::floating_point T>
bool sweep_events(std::span<const BoxBernstein<T>> forms, std::span<const SelectionTerm> terms, int axis,
                  std::vector<T>& events, int max_boxes)
{
    if (axis < 0 || axis > 2)
        throw std::invalid_argument("quadrays: a sweep runs along box axis 0, 1 or 2");
    thread_local std::vector<Jet<T>> jets;
    thread_local std::vector<System> systems;
    thread_local std::vector<Range<T>> found;
    make_jets(forms, jets);
    const int n = static_cast<int>(forms.size());
    const int p = (axis + 1) % 3, q = (axis + 2) % 3;
    const T width = scaled_tolerance<T>(1e-7);
    int budget = max_boxes;
    found.clear();

    // inside the box
    systems.clear();
    for (int a = 0; a < n; ++a)
    {
        systems.push_back({SystemKind::critical, a, 0, 0});
        for (int b = a + 1; b < n; ++b)
        {
            systems.push_back({SystemKind::tangency, a, b, 0});
            for (int c = b + 1; c < n; ++c)
                systems.push_back({SystemKind::triple, a, b, c});
        }
    }
    if (!search(jets, terms, systems, Vec3<T>{0, 0, 0}, Vec3<T>{1, 1, 1}, axis, p, q, width, found, budget))
        return false;

    // on the four lateral faces: the face fixes one of the two other axes
    for (int side = 0; side < 4; ++side)
    {
        const int fixed = side < 2 ? p : q, open = side < 2 ? q : p;
        Vec3<T> lo = {0, 0, 0}, hi = {1, 1, 1};
        lo[fixed] = hi[fixed] = T(side % 2);
        systems.clear();
        for (int a = 0; a < n; ++a)
        {
            systems.push_back({SystemKind::face_tangency, a, 0, 0});
            for (int b = a + 1; b < n; ++b)
                systems.push_back({SystemKind::pair, a, b, 0});
        }
        if (!search(jets, terms, systems, lo, hi, axis, open, open, width, found, budget))
            return false;
    }
    // on the four lateral edges
    systems.clear();
    for (int a = 0; a < n; ++a)
        systems.push_back({SystemKind::single, a, 0, 0});
    for (int corner = 0; corner < 4; ++corner)
    {
        Vec3<T> lo = {0, 0, 0}, hi = {1, 1, 1};
        lo[p] = hi[p] = T(corner & 1);
        lo[q] = hi[q] = T(corner >> 1);
        if (!search(jets, terms, systems, lo, hi, axis, p, q, width, found, budget))
            return false;
    }

    // boxes that touch are one event; a longer run of them is kept as its two ends
    std::sort(found.begin(), found.end(), [](const Range<T>& x, const Range<T>& y) { return x.lo < y.lo; });
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
bool may_hold_ridge(std::span<const BoxBernstein<T>> forms, std::span<const SelectionTerm> terms)
{
    thread_local std::vector<Jet<T>> jets;
    thread_local std::vector<System> systems;
    thread_local std::vector<Range<T>> found;
    make_jets(forms, jets);
    systems.clear();
    const int n = static_cast<int>(forms.size());
    for (int a = 0; a < n; ++a)
        for (int b = a + 1; b < n; ++b)
            systems.push_back({SystemKind::pair, a, b, 0});
    if (systems.empty())
        return false;
    found.clear();
    int budget = 20000;
    const bool finished = search(jets, terms, systems, Vec3<T>{0, 0, 0}, Vec3<T>{1, 1, 1}, -1, 0, 0, T(1) / T(64),
                                 found, budget);
    return !finished || !found.empty();
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
bool integrate_by_slices(const ClippedBox<T>& cell, std::span<const BoxBernstein<T>> forms,
                         std::span<const SelectionTerm> terms, int q, const Options& opt, CellPoints<T>& out,
                         Stats& stats)
{
    if (cell.tdim != 3 || !cell.clips.empty() || forms.size() < 2)
        return false;
    for (const BoxBernstein<T>& form : forms)
        if (form.dim != 3)
            return false;
    for (const SelectionTerm& term : terms)
        if (term.zero_required != 0)
            return false;
    if (!may_hold_ridge(forms, terms))
        return false;
    thread_local std::vector<T> events, best, positions, weights;
    int axis = -1;
    for (int k = 0; k < 3; ++k)
        if (sweep_events(forms, terms, k, events) && (axis < 0 || events.size() < best.size()))
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
        integrate_slices(cell, forms, terms, sweep, std::span<const T>(positions), std::span<const T>(weights), q, opt,
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
    if (!(volume > T(0)))
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

template bool sweep_events<float>(std::span<const BoxBernstein<float>>, std::span<const SelectionTerm>, int,
                                  std::vector<float>&, int);
template bool sweep_events<double>(std::span<const BoxBernstein<double>>, std::span<const SelectionTerm>, int,
                                   std::vector<double>&, int);
template bool may_hold_ridge<float>(std::span<const BoxBernstein<float>>, std::span<const SelectionTerm>);
template bool may_hold_ridge<double>(std::span<const BoxBernstein<double>>, std::span<const SelectionTerm>);
template void piece_rule<float>(float, bool, float, bool, int, std::vector<float>&, std::vector<float>&);
template void piece_rule<double>(double, bool, double, bool, int, std::vector<double>&, std::vector<double>&);
template bool integrate_by_slices<float>(const ClippedBox<float>&, std::span<const BoxBernstein<float>>,
                                         std::span<const SelectionTerm>, int, const Options&, CellPoints<float>&,
                                         Stats&);
template bool integrate_by_slices<double>(const ClippedBox<double>&, std::span<const BoxBernstein<double>>,
                                          std::span<const SelectionTerm>, int, const Options&,
                                          CellPoints<double>&, Stats&);

template bool make_slice<float>(const ClippedBox<float>&, const Sweep<float>&, float, Slice<float>&);
template bool make_slice<double>(const ClippedBox<double>&, const Sweep<double>&, double, Slice<double>&);
template void integrate_slices<float>(const ClippedBox<float>&, std::span<const BoxBernstein<float>>,
                                      std::span<const SelectionTerm>, const Sweep<float>&, std::span<const float>,
                                      std::span<const float>, int, const Options&, CellPoints<float>&, Stats&);
template void integrate_slices<double>(const ClippedBox<double>&, std::span<const BoxBernstein<double>>,
                                       std::span<const SelectionTerm>, const Sweep<double>&,
                                       std::span<const double>, std::span<const double>, int, const Options&,
                                       CellPoints<double>&, Stats&);

} // namespace cutcells::quadrays
