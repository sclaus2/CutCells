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

/// The level set whose zero set the terms ask for: -1 for a volume part.
/// @throws std::invalid_argument if they ask for more than one, or differ in it
inline int zero_set_of(std::span<const SelectionTerm> terms)
{
    if (terms.empty())
        throw std::invalid_argument("quadrays: slices need a selection");
    const std::uint64_t zero = terms[0].zero_required;
    for (const SelectionTerm& term : terms)
        if (term.zero_required != zero)
            throw std::invalid_argument("quadrays: the terms of a sliced part ask for the same zero set");
    if (zero == 0)
        return -1;
    if (zero & (zero - 1))
        throw std::invalid_argument("quadrays: slices integrate volumes and the zero set of one level set");
    int l = 0;
    while (!((zero >> l) & 1))
        ++l;
    return l;
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
    const int zero = zero_set_of(terms);
    const Mat3<T> inverse = inverse_jacobian(cell);
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
            T factor = sweep.kind == SweepKind::turning ? slice.radius_origin + slice.radius_step * s0 : T(1);
            if (zero >= 0)
            {
                // a curve of the zero set in the slice: its length to the surface's measure, by the
                // angle between the zero set's normal and the slice's
                const std::size_t first = out.points.size() - 3;
                const Vec3<T> u = {out.points[first], out.points[first + 1], out.points[first + 2]};
                Vec3<T> du;
                gradient(forms[static_cast<std::size_t>(zero)], std::span<const T>(u), std::span<T>(du));
                Vec3<T> dx = {0, 0, 0};
                for (std::size_t r = 0; r < 3; ++r)
                    for (std::size_t c = 0; c < 3; ++c)
                        dx[c] += inverse[r][c] * du[r];
                const Vec3<T> normal = {slice.box.jacobian[0][2], slice.box.jacobian[1][2], slice.box.jacobian[2][2]};
                const T cosine = dot(dx, normal) / std::sqrt(dot(dx, dx));
                factor /= std::sqrt(std::max(T(1) - cosine * cosine, std::numeric_limits<T>::min()));
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

/// A level set and its partial derivatives on the box, with their sizes.
template <std::floating_point T>
struct Jet
{
    std::array<BoxBernstein<T>, 4> form; ///< value, d/du_0, d/du_1, d/du_2
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

/// The range of function d (0: value, 1 + i: d/du_i) of level set l on the box.
template <std::floating_point T>
const Range<T>& range_of(BoxRanges<T>& r, int l, int d)
{
    const std::size_t ls = static_cast<std::size_t>(l), ds = static_cast<std::size_t>(d);
    if (!r.known[ls][ds])
    {
        const std::array<T, 9> matrix = {r.hi[0] - r.lo[0], 0, 0, 0, r.hi[1] - r.lo[1], 0, 0, 0, r.hi[2] - r.lo[2]};
        restrict_affine(r.cell->jets[ls].form[ds], std::span<const T>(r.lo), std::span<const T>(matrix), 3, r.scratch,
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

/// Bisect the unit box and keep the boxes on which some system may hold and
/// which may meet the part, down to @p width along @p axis; their ranges
/// along @p axis go to @p found.
/// @return false if more than @p budget boxes were visited
template <std::floating_point T>
bool search(const Geometry<T>& cell, std::span<const SelectionTerm> terms, const std::vector<System>& systems,
            int axis, int p, int q, T width, std::vector<Range<T>>& found, int& budget)
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
        bool any = false;
        for (const System& s : systems)
            if (may_hold(ranges, s, p, q))
            {
                any = true;
                break;
            }
        if (!any)
            continue;
        if (box.hi[k] - box.lo[k] <= width)
        {
            found.push_back({box.lo[k], box.hi[k]});
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
void make_geometry(const ClippedBox<T>& cell, std::span<const BoxBernstein<T>> forms, Geometry<T>& g)
{
    g.jets.resize(forms.size());
    for (std::size_t l = 0; l < forms.size(); ++l)
    {
        if (forms[l].dim != 3)
            throw std::invalid_argument("quadrays: the events of a sweep need level sets on a box of dimension 3");
        g.jets[l].form[0] = forms[l];
        for (int d = 0; d < 3; ++d)
            derivative(forms[l], d, g.jets[l].form[static_cast<std::size_t>(1 + d)]);
        for (std::size_t d = 0; d < 4; ++d)
            g.jets[l].size[d] = max_abs(std::span<const T>(g.jets[l].form[d].coeffs));
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

} // namespace

template <std::floating_point T>
bool sweep_events(const ClippedBox<T>& cell, std::span<const BoxBernstein<T>> forms,
                  std::span<const SelectionTerm> terms, int axis, std::vector<T>& events, int max_boxes)
{
    if (axis < 0 || axis > 2)
        throw std::invalid_argument("quadrays: a sweep runs along box axis 0, 1 or 2");
    if (cell.tdim != 3)
        throw std::invalid_argument("quadrays: slices are of cells of dimension 3");
    const int zero = zero_set_of(terms);
    thread_local Geometry<T> g;
    thread_local std::vector<System> systems;
    thread_local std::vector<Range<T>> found;
    make_geometry(cell, forms, g);
    const int n = static_cast<int>(forms.size());
    const int p = (axis + 1) % 3, q = (axis + 2) % 3;
    const T width = scaled_tolerance<T>(1e-7);

    // the systems: of all zero sets for a volume, of those with the part's own for a surface
    systems.clear();
    const auto takes = [&](int a) { return zero < 0 || a == zero; };
    const int n_faces = static_cast<int>(g.faces.size());
    for (int a = 0; a < n; ++a)
    {
        if (takes(a))
            systems.push_back({SystemKind::critical, a, 0, 0});
        for (int b = a + 1; b < n; ++b)
        {
            if (!takes(a) && !takes(b))
                continue;
            systems.push_back({SystemKind::tangency, a, b, 0});
            for (int c = b + 1; c < n; ++c)
                systems.push_back({SystemKind::triple, a, b, c});
            for (int f = 0; f < n_faces; ++f)
                systems.push_back({SystemKind::face_pair, a, b, f});
        }
        if (!takes(a))
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

    // the cell's vertices on the part: the pattern of the slice polygon changes there
    thread_local std::vector<T> values;
    const std::size_t k = static_cast<std::size_t>(axis);
    for (const Vec3<T>& v : g.vertices)
    {
        values.resize(forms.size());
        for (std::size_t l = 0; l < forms.size(); ++l)
            values[l] = evaluate(forms[l], std::span<const T>(v));
        const auto at = [&](int l) { return Range<T>{values[static_cast<std::size_t>(l)], values[static_cast<std::size_t>(l)]}; };
        bool on_part = false;
        for (const SelectionTerm& term : terms)
            on_part = on_part || term_may_hold(term, g.jets, at);
        if (zero >= 0)
            on_part = on_part && std::abs(values[static_cast<std::size_t>(zero)])
                                     <= scaled_tolerance<T>(1e-9) * g.jets[static_cast<std::size_t>(zero)].size[0];
        if (on_part)
            found.push_back({v[k], v[k]});
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
    if (cell.tdim != 3 || forms.empty() || terms.empty())
        return false;
    for (const BoxBernstein<T>& form : forms)
        if (form.dim != 3)
            return false;
    const std::uint64_t zero = terms[0].zero_required;
    if (zero & (zero - 1))
        return false;
    for (const SelectionTerm& term : terms)
        if (term.zero_required != zero)
            return false;
    thread_local std::vector<T> events, best, positions, weights;
    int axis = -1;
    for (int k = 0; k < 3; ++k)
        if (sweep_events(cell, forms, terms, k, events) && (axis < 0 || events.size() < best.size()))
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
    if (!(volume > T(0))) // nothing found; the measure is an area for a surface
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

template bool sweep_events<float>(const ClippedBox<float>&, std::span<const BoxBernstein<float>>,
                                  std::span<const SelectionTerm>, int, std::vector<float>&, int);
template bool sweep_events<double>(const ClippedBox<double>&, std::span<const BoxBernstein<double>>,
                                   std::span<const SelectionTerm>, int, std::vector<double>&, int);
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
