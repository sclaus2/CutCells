// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "slicing.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

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
// Explicit instantiations
// ============================================================================

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
