// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "clipped_box.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "box_bernstein.h"

namespace cutcells::quadrays
{

namespace
{
template <std::floating_point T>
T det3(const Mat3<T>& a)
{
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
           - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
           + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

template <std::floating_point T>
Vec3<T> affine(const Vec3<T>& origin, const Mat3<T>& m, const Vec3<T>& u)
{
    Vec3<T> y = origin;
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            y[i] += m[i][k] * u[k];
    return y;
}

template <std::floating_point T>
Mat3<T> identity()
{
    Mat3<T> m = {};
    for (int i = 0; i < 3; ++i)
        m[i][i] = T(1);
    return m;
}

template <std::floating_point T>
Mat3<T> inverse(const Mat3<T>& a)
{
    const T det = det3(a);
    Mat3<T> inv;
    inv[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
    inv[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det;
    inv[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
    inv[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
    inv[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det;
    inv[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;
    inv[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
    inv[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / det;
    inv[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;
    return inv;
}

template <std::floating_point T>
Mat3<T> product(const Mat3<T>& a, const Mat3<T>& b)
{
    Mat3<T> c = {};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                c[i][j] += a[i][k] * b[k][j];
    return c;
}

template <std::floating_point T>
T norm(const Vec3<T>& a)
{
    return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

/// Entries within a few ulps of 0 or +-1 snapped, so that frames along the
/// physical axes are exact.
template <std::floating_point T>
void snap(Mat3<T>& f)
{
    const T eps = T(16) * std::numeric_limits<T>::epsilon();
    for (auto& row : f)
        for (T& x : row)
        {
            if (std::abs(x) < eps)
                x = T(0);
            else if (std::abs(std::abs(x) - T(1)) < eps)
                x = x > T(0) ? T(1) : T(-1);
        }
}

/// The bounding box of @p points (relative to their first) in the frame
/// @p f (rows the axes) along its first @p naxes axes: lower ends and lengths.
template <std::floating_point T>
T box_in_frame(const Mat3<T>& f, std::span<const Vec3<T>> points, int naxes, Vec3<T>& lo, Vec3<T>& len)
{
    T measure = T(1);
    for (int k = 0; k < naxes; ++k)
    {
        T a = std::numeric_limits<T>::infinity(), b = -std::numeric_limits<T>::infinity();
        for (const Vec3<T>& x : points)
        {
            T p = T(0);
            for (int i = 0; i < 3; ++i)
                p += f[k][i] * (x[i] - points[0][i]);
            a = std::min(a, p);
            b = std::max(b, p);
        }
        lo[k] = a;
        len[k] = b - a;
        measure *= b - a;
    }
    return measure;
}

/// Of the candidate frames, the one whose bounding box of @p points is
/// smallest along the first @p naxes axes; ties keep the earlier frame.
template <std::floating_point T>
std::size_t smallest_box(std::span<Mat3<T>> frames, std::span<const Vec3<T>> points, int naxes)
{
    std::size_t best = 0;
    T best_measure = std::numeric_limits<T>::infinity();
    for (std::size_t f = 0; f < frames.size(); ++f)
    {
        snap(frames[f]);
        Vec3<T> lo, len;
        const T measure = box_in_frame(frames[f], points, naxes, lo, len);
        if (measure < best_measure * (T(1) - T(1e-9)))
        {
            best = f;
            best_measure = measure;
        }
    }
    return best;
}

/// The frame of a unit vector @p e0 and a second direction @p b
/// (Gram-Schmidt), right-handed; false if b is nearly parallel to e0.
template <std::floating_point T>
bool frame_from(const Vec3<T>& e0, const Vec3<T>& b, Mat3<T>& f)
{
    const T dot = b[0] * e0[0] + b[1] * e0[1] + b[2] * e0[2];
    Vec3<T> w = {b[0] - dot * e0[0], b[1] - dot * e0[1], b[2] - dot * e0[2]};
    const T lw = norm(w);
    if (!(lw > T(1e-6) * norm(b)))
        return false;
    for (T& c : w)
        c /= lw;
    f = {e0, w, Vec3<T>{e0[1] * w[2] - e0[2] * w[1], e0[2] * w[0] - e0[0] * w[2], e0[0] * w[1] - e0[1] * w[0]}};
    return true;
}

/// The reference map xi -> v0 + a xi of a cell, the box maps x(u) set in
/// @p box, and the cell's facets c . xi <= d in reference coordinates: the
/// box's reference map, and the facets as clips unless they hold on the
/// whole box.
template <std::floating_point T>
void finish_box(const Vec3<T>& v0, const Mat3<T>& a, std::span<const HalfSpace<T>> facets, ClippedBox<T>& box)
{
    const Mat3<T> ainv = inverse(a);
    box.ref_jacobian = product(ainv, box.jacobian);
    box.ref_origin = {0, 0, 0};
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            box.ref_origin[i] += ainv[i][k] * (box.origin[k] - v0[k]);
    if (box.tdim == 2)
    {
        box.ref_origin[2] = T(0);
        box.ref_jacobian[2] = {T(0), T(0), T(1)};
        for (int i = 0; i < 2; ++i)
            box.ref_jacobian[i][2] = T(0);
    }
    box.clips.clear();
    const T tol = scaled_tolerance<T>(1e-12);
    for (const HalfSpace<T>& g : facets)
    {
        // c . (ref_origin + ref_jacobian u) <= d
        HalfSpace<T> h;
        h.d = g.d;
        for (int i = 0; i < 3; ++i)
        {
            h.d -= g.c[i] * box.ref_origin[i];
            for (int k = 0; k < box.tdim; ++k)
                h.c[k] += g.c[i] * box.ref_jacobian[i][k];
        }
        T m = T(0);
        for (int k = 0; k < box.tdim; ++k)
            m = std::max(m, std::abs(h.c[k]));
        if (!(m > T(0)))
            continue;
        for (int k = 0; k < box.tdim; ++k)
            h.c[k] /= m;
        h.d /= m;
        T vmax = T(0); // the largest value of c . u over the box
        for (int k = 0; k < box.tdim; ++k)
            vmax += std::max(h.c[k], T(0));
        if (vmax > h.d + tol)
            box.clips.push_back(h);
    }
}

/// The reference map of a triangle, tetrahedron or prism (vertex 0 and the
/// edges to vertices 1, 2 and 3 as columns), and false for a degenerate cell.
template <std::floating_point T>
bool reference_map(int tdim, std::span<const Vec3<T>> v, Mat3<T>& a)
{
    a = {};
    T scale = T(1);
    for (int k = 0; k < tdim; ++k)
    {
        Vec3<T> e;
        for (int i = 0; i < 3; ++i)
            a[i][k] = e[i] = v[k + 1][i] - v[0][i];
        scale *= norm(e);
    }
    if (tdim == 2)
        a[2][2] = T(1);
    return std::abs(det3(a)) > T(1e-10) * scale;
}

/// The box of a triangle or tetrahedron in BoxFrame::orthogonal; false (box
/// untouched) for a degenerate simplex.
template <std::floating_point T>
bool orthogonal_simplex_box(int tdim, std::span<const T> x, ClippedBox<T>& box)
{
    const int nv = tdim + 1;
    std::array<Vec3<T>, 4> v{};
    for (int j = 0; j < nv; ++j)
        for (int i = 0; i < tdim; ++i)
            v[j][i] = x[j * tdim + i];
    const std::span<const Vec3<T>> vs(v.data(), static_cast<std::size_t>(nv));
    Mat3<T> a;
    if (!reference_map<T>(tdim, vs, a))
        return false;

    // candidate frames, rows the axes: the physical axes first, then those two
    // edges span (in 2D one edge and its normal)
    std::array<Vec3<T>, 6> edges;
    int ne = 0;
    for (int i = 0; i < nv; ++i)
        for (int j = i + 1; j < nv; ++j)
            edges[ne++] = {v[j][0] - v[i][0], v[j][1] - v[i][1], v[j][2] - v[i][2]};
    std::array<Mat3<T>, 31> frames;
    int nf = 0;
    frames[nf++] = identity<T>();
    for (int p = 0; p < ne; ++p)
    {
        const T lp = norm(edges[p]);
        const Vec3<T> e0 = {edges[p][0] / lp, edges[p][1] / lp, edges[p][2] / lp};
        if (tdim == 2)
            frames[nf++] = {e0, Vec3<T>{-e0[1], e0[0], T(0)}, Vec3<T>{T(0), T(0), T(1)}};
        else
            for (int q = 0; q < ne; ++q)
                if (q != p && frame_from(e0, edges[q], frames[nf]))
                    ++nf;
    }
    const Mat3<T>& r = frames[smallest_box<T>(std::span<Mat3<T>>(frames.data(), static_cast<std::size_t>(nf)), vs, tdim)];
    Vec3<T> lo, len;
    box_in_frame<T>(r, vs, tdim, lo, len);

    // x(u) = v0 + r^T (lo + diag(len) u)
    box.tdim = tdim;
    box.origin = v[0];
    box.jacobian = {};
    for (int k = 0; k < tdim; ++k)
        for (int i = 0; i < 3; ++i)
        {
            box.origin[i] += r[k][i] * lo[k];
            box.jacobian[i][k] = r[k][i] * len[k];
        }
    if (tdim == 2)
        box.jacobian[2][2] = T(1);
    // facets: xi_k >= 0 and xi_0 + ... <= 1
    std::array<HalfSpace<T>, 4> facets{};
    for (int k = 0; k < tdim; ++k)
    {
        facets[k].c[k] = T(-1);
        facets[tdim].c[k] = T(1);
    }
    facets[tdim].d = T(1);
    finish_box<T>(v[0], a, std::span<const HalfSpace<T>>(facets.data(), static_cast<std::size_t>(tdim + 1)), box);
    return true;
}

/// The box of a prism in BoxFrame::orthogonal: its bottom triangle's
/// smallest bounding rectangle in an orthonormal frame of its plane (the
/// frames of its edges, which contain the physical axes when an edge lies
/// along one), times the lateral edge v3 - v0, clipped by the triangle's
/// sides. Its level set keeps degree n in each box variable. False (box
/// untouched) for a degenerate prism.
template <std::floating_point T>
bool orthogonal_prism_box(std::span<const T> x, ClippedBox<T>& box)
{
    std::array<Vec3<T>, 6> v;
    for (int j = 0; j < 6; ++j)
        for (int i = 0; i < 3; ++i)
            v[j][i] = x[j * 3 + i];
    Mat3<T> a;
    if (!reference_map<T>(3, std::span<const Vec3<T>>(v), a))
        return false;
    const std::span<const Vec3<T>> bottom(v.data(), 3);
    const Vec3<T> e1 = {v[1][0] - v[0][0], v[1][1] - v[0][1], v[1][2] - v[0][2]},
                  e2 = {v[2][0] - v[0][0], v[2][1] - v[0][1], v[2][2] - v[0][2]};
    Vec3<T> nrm = {e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]};
    const T ln = norm(nrm);
    for (T& c : nrm)
        c /= ln;
    // in-plane frames: each edge and the normal's cross product with it
    std::array<Mat3<T>, 3> frames;
    for (int p = 0; p < 3; ++p)
    {
        const Vec3<T>& from = v[p];
        const Vec3<T>& to = v[(p + 1) % 3];
        Vec3<T> e0 = {to[0] - from[0], to[1] - from[1], to[2] - from[2]};
        const T l0 = norm(e0);
        for (T& c : e0)
            c /= l0;
        frames[p] = {e0, Vec3<T>{nrm[1] * e0[2] - nrm[2] * e0[1], nrm[2] * e0[0] - nrm[0] * e0[2],
                                 nrm[0] * e0[1] - nrm[1] * e0[0]},
                     nrm};
    }
    const Mat3<T>& r = frames[smallest_box<T>(std::span<Mat3<T>>(frames), bottom, 2)];
    Vec3<T> lo, len;
    box_in_frame<T>(r, bottom, 2, lo, len);

    // x(u) = v0 + r0 (lo0 + len0 u0) + r1 (lo1 + len1 u1) + (v3 - v0) u2
    box.tdim = 3;
    box.origin = v[0];
    box.jacobian = {};
    for (int i = 0; i < 3; ++i)
    {
        for (int k = 0; k < 2; ++k)
        {
            box.origin[i] += r[k][i] * lo[k];
            box.jacobian[i][k] = r[k][i] * len[k];
        }
        box.jacobian[i][2] = v[3][i] - v[0][i];
    }
    // facets: the sides xi0 >= 0, xi1 >= 0, xi0 + xi1 <= 1; the top and bottom are box faces
    std::array<HalfSpace<T>, 3> facets{};
    facets[0].c[0] = T(-1);
    facets[1].c[1] = T(-1);
    facets[2].c = {T(1), T(1), T(0)};
    facets[2].d = T(1);
    finish_box<T>(v[0], a, std::span<const HalfSpace<T>>(facets), box);
    return true;
}

/// Planes of the clipped region as c . u <= d: the six box faces, then the clips.
template <std::floating_point T>
std::vector<HalfSpace<T>> region_planes(const ClippedBox<T>& box)
{
    std::vector<HalfSpace<T>> planes;
    for (int k = 0; k < 3; ++k)
    {
        HalfSpace<T> lower, upper;
        lower.c[k] = T(-1);
        lower.d = T(0);
        upper.c[k] = T(1);
        upper.d = T(1);
        planes.push_back(lower);
        planes.push_back(upper);
    }
    planes.insert(planes.end(), box.clips.begin(), box.clips.end());
    return planes;
}
} // namespace

template <std::floating_point T>
void make_clipped_box(cell::type cell_type, std::span<const T> vertex_coords, int gdim,
                      ClippedBox<T>& box, BoxFrame frame)
{
    if (!supported_cell(cell_type))
        throw std::invalid_argument("quadrays: unsupported cell type " + cell::cell_type_to_str(cell_type));
    const int tdim = cell::get_tdim(cell_type);
    if (gdim != tdim)
    {
        throw std::invalid_argument("quadrays: a " + cell::cell_type_to_str(cell_type) + " must lie in "
                                    + std::to_string(tdim) + "D (gdim = " + std::to_string(tdim) + ")");
    }
    if (static_cast<int>(vertex_coords.size()) != cell::get_num_vertices(cell_type) * gdim)
        throw std::invalid_argument("quadrays: wrong number of vertex coordinates");
    if (frame == BoxFrame::orthogonal)
    {
        if ((cell_type == cell::type::triangle || cell_type == cell::type::tetrahedron)
            && orthogonal_simplex_box(tdim, vertex_coords, box))
            return;
        if (cell_type == cell::type::prism && orthogonal_prism_box(vertex_coords, box))
            return;
    }
    // Basix vertices spanning the box from vertex 0
    std::array<int, 3> axes = {1, 2, 0};
    if (cell_type == cell::type::tetrahedron || cell_type == cell::type::prism)
        axes[2] = 3;
    else if (cell_type == cell::type::hexahedron || cell_type == cell::type::pyramid)
        axes[2] = 4;

    box.tdim = tdim;
    box.origin = {0, 0, 0};
    box.jacobian = {};
    for (int i = 0; i < gdim; ++i)
    {
        box.origin[i] = vertex_coords[i];
        for (int k = 0; k < tdim; ++k)
            box.jacobian[i][k] = vertex_coords[axes[k] * gdim + i] - vertex_coords[i];
    }
    if (tdim == 2)
        box.jacobian[2][2] = T(1);
    box.ref_origin = {0, 0, 0};
    box.ref_jacobian = identity<T>();
    box.clips.clear();
    switch (cell_type)
    {
    case cell::type::triangle:
    case cell::type::prism:
        box.clips.push_back({{T(1), T(1), T(0)}, T(1)});
        break;
    case cell::type::tetrahedron:
        box.clips.push_back({{T(1), T(1), T(1)}, T(1)});
        break;
    case cell::type::pyramid:
        box.clips.push_back({{T(1), T(0), T(1)}, T(1)});
        box.clips.push_back({{T(0), T(1), T(1)}, T(1)});
        break;
    default:
        break;
    }
}

template <std::floating_point T>
T jacobian_determinant(const ClippedBox<T>& box)
{
    return det3(box.jacobian);
}

template <std::floating_point T>
T reference_jacobian_determinant(const ClippedBox<T>& box)
{
    // dx/dxi = (dx/du) (dxi/du)^-1
    return det3(box.jacobian) / det3(box.ref_jacobian);
}

template <std::floating_point T>
Mat3<T> inverse_jacobian(const ClippedBox<T>& box)
{
    if (det3(box.jacobian) == T(0))
        throw std::runtime_error("quadrays: singular box map");
    return inverse(box.jacobian);
}

template <std::floating_point T>
Vec3<T> physical_point(const ClippedBox<T>& box, const Vec3<T>& u)
{
    return affine(box.origin, box.jacobian, u);
}

template <std::floating_point T>
Vec3<T> reference_point(const ClippedBox<T>& box, const Vec3<T>& u)
{
    return affine(box.ref_origin, box.ref_jacobian, u);
}

template <std::floating_point T>
Vec3<T> box_point(const ClippedBox<T>& box, const Vec3<T>& xi)
{
    if (reference_frame(box))
        return xi;
    const Mat3<T> inv = inverse(box.ref_jacobian);
    Vec3<T> u = {0, 0, 0};
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            u[i] += inv[i][k] * (xi[k] - box.ref_origin[k]);
    if (box.tdim == 2)
        u[2] = T(0);
    return u;
}

template <std::floating_point T>
bool reference_frame(const ClippedBox<T>& box)
{
    for (int i = 0; i < 3; ++i)
    {
        if (box.ref_origin[i] != T(0))
            return false;
        for (int k = 0; k < 3; ++k)
            if (box.ref_jacobian[i][k] != (i == k ? T(1) : T(0)))
                return false;
    }
    return true;
}

template <std::floating_point T>
ClippedBox<T> sub_box(const ClippedBox<T>& box, const Vec3<T>& lo, const Vec3<T>& hi)
{
    ClippedBox<T> child;
    child.tdim = box.tdim;
    child.origin = physical_point(box, lo);
    child.ref_origin = reference_point(box, lo);
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
        {
            child.jacobian[i][k] = box.jacobian[i][k] * (hi[k] - lo[k]);
            child.ref_jacobian[i][k] = box.ref_jacobian[i][k] * (hi[k] - lo[k]);
        }
    for (const HalfSpace<T>& h : box.clips)
    {
        // c . (lo + D u') <= d  ->  (c o D) . u' <= d - c . lo
        HalfSpace<T> g;
        g.d = h.d;
        for (int k = 0; k < 3; ++k)
        {
            g.c[k] = h.c[k] * (hi[k] - lo[k]);
            g.d -= h.c[k] * lo[k];
        }
        child.clips.push_back(g);
    }
    return child;
}

template <std::floating_point T>
int longest_axis(const ClippedBox<T>& box)
{
    int best = 0;
    T best_len = T(-1);
    for (int k = 0; k < box.tdim; ++k)
    {
        T len = T(0);
        for (int i = 0; i < 3; ++i)
            len += box.jacobian[i][k] * box.jacobian[i][k];
        if (len > best_len)
        {
            best_len = len;
            best = k;
        }
    }
    return best;
}

template <std::floating_point T>
bool may_meet_clips(const ClippedBox<T>& box, const Vec3<T>& lo, const Vec3<T>& hi, T tol)
{
    for (const HalfSpace<T>& h : box.clips)
    {
        // the smallest value of c . u over the sub-box is attained at a corner
        T vmin = T(0);
        for (int k = 0; k < 3; ++k)
            vmin += h.c[k] * (h.c[k] >= T(0) ? lo[k] : hi[k]);
        if (vmin > h.d + tol)
            return false;
    }
    return true;
}

template <std::floating_point T>
bool inside_clips(const ClippedBox<T>& box, const Vec3<T>& u, T tol)
{
    for (const HalfSpace<T>& h : box.clips)
        if (h.c[0] * u[0] + h.c[1] * u[1] + h.c[2] * u[2] > h.d + tol)
            return false;
    return true;
}

template <std::floating_point T>
void clipped_polytope(const ClippedBox<T>& box, Polytope<T>& poly)
{
    const std::vector<HalfSpace<T>> planes = region_planes(box);
    const int np = static_cast<int>(planes.size());
    const T tol = scaled_tolerance<T>(1e-11);
    const T det_tol = scaled_tolerance<T>(1e-14);
    const T same_tol = scaled_tolerance<T>(1e-10);

    poly = Polytope<T>{};
    std::vector<std::vector<int>> active; // planes active at each vertex
    for (int a = 0; a < np; ++a)
        for (int b = a + 1; b < np; ++b)
            for (int c = b + 1; c < np; ++c)
            {
                const Mat3<T> m = {planes[a].c, planes[b].c, planes[c].c};
                const T det = det3(m);
                if (std::abs(det) < det_tol)
                    continue;
                Vec3<T> x;
                for (int col = 0; col < 3; ++col)
                {
                    Mat3<T> mc = m;
                    mc[0][col] = planes[a].d;
                    mc[1][col] = planes[b].d;
                    mc[2][col] = planes[c].d;
                    x[col] = det3(mc) / det;
                }
                bool feasible = true;
                for (const HalfSpace<T>& h : planes)
                    if (h.c[0] * x[0] + h.c[1] * x[1] + h.c[2] * x[2] > h.d + tol)
                    {
                        feasible = false;
                        break;
                    }
                if (!feasible)
                    continue;
                bool duplicate = false;
                for (const Vec3<T>& v : poly.vertices)
                    if (std::abs(v[0] - x[0]) + std::abs(v[1] - x[1]) + std::abs(v[2] - x[2]) < same_tol)
                    {
                        duplicate = true;
                        break;
                    }
                if (duplicate)
                    continue;
                std::vector<int> act;
                for (int p = 0; p < np; ++p)
                {
                    const HalfSpace<T>& h = planes[p];
                    if (std::abs(h.c[0] * x[0] + h.c[1] * x[1] + h.c[2] * x[2] - h.d) < tol)
                        act.push_back(p);
                }
                poly.vertices.push_back(x);
                active.push_back(act);
            }

    auto shared = [&](int v, int w)
    {
        int count = 0;
        for (int p : active[v])
            for (int q : active[w])
                count += p == q;
        return count;
    };
    const int nv = poly.n_vertices();
    for (int v = 0; v < nv; ++v)
        for (int w = v + 1; w < nv; ++w)
            if (shared(v, w) >= 2)
                poly.edges.push_back({v, w});

    for (int p = 0; p < np; ++p)
    {
        std::vector<int> on;
        for (int v = 0; v < nv; ++v)
            for (int q : active[v])
                if (q == p)
                    on.push_back(v);
        if (on.size() < 3)
            continue;
        for (int e = 0; e < static_cast<int>(poly.edges.size()); ++e)
        {
            bool a_on = false, b_on = false;
            for (int v : on)
            {
                a_on |= v == poly.edges[e][0];
                b_on |= v == poly.edges[e][1];
            }
            if (a_on && b_on)
                poly.face_edges.push_back(e);
        }
        poly.face_offsets.push_back(static_cast<int>(poly.face_edges.size()));
        const Vec3<T>& c = planes[p].c;
        const T n = std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
        poly.face_normals.push_back({c[0] / n, c[1] / n, c[2] / n});
    }
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template void make_clipped_box<float>(cell::type, std::span<const float>, int, ClippedBox<float>&, BoxFrame);
template Vec3<float> box_point<float>(const ClippedBox<float>&, const Vec3<float>&);
template bool reference_frame<float>(const ClippedBox<float>&);
template float jacobian_determinant<float>(const ClippedBox<float>&);
template float reference_jacobian_determinant<float>(const ClippedBox<float>&);
template Mat3<float> inverse_jacobian<float>(const ClippedBox<float>&);
template Vec3<float> physical_point<float>(const ClippedBox<float>&, const Vec3<float>&);
template Vec3<float> reference_point<float>(const ClippedBox<float>&, const Vec3<float>&);
template ClippedBox<float> sub_box<float>(const ClippedBox<float>&, const Vec3<float>&, const Vec3<float>&);
template int longest_axis<float>(const ClippedBox<float>&);
template bool may_meet_clips<float>(const ClippedBox<float>&, const Vec3<float>&, const Vec3<float>&, float);
template bool inside_clips<float>(const ClippedBox<float>&, const Vec3<float>&, float);
template void clipped_polytope<float>(const ClippedBox<float>&, Polytope<float>&);

template void make_clipped_box<double>(cell::type, std::span<const double>, int, ClippedBox<double>&, BoxFrame);
template Vec3<double> box_point<double>(const ClippedBox<double>&, const Vec3<double>&);
template bool reference_frame<double>(const ClippedBox<double>&);
template double jacobian_determinant<double>(const ClippedBox<double>&);
template double reference_jacobian_determinant<double>(const ClippedBox<double>&);
template Mat3<double> inverse_jacobian<double>(const ClippedBox<double>&);
template Vec3<double> physical_point<double>(const ClippedBox<double>&, const Vec3<double>&);
template Vec3<double> reference_point<double>(const ClippedBox<double>&, const Vec3<double>&);
template ClippedBox<double> sub_box<double>(const ClippedBox<double>&, const Vec3<double>&, const Vec3<double>&);
template int longest_axis<double>(const ClippedBox<double>&);
template bool may_meet_clips<double>(const ClippedBox<double>&, const Vec3<double>&, const Vec3<double>&, double);
template bool inside_clips<double>(const ClippedBox<double>&, const Vec3<double>&, double);
template void clipped_polytope<double>(const ClippedBox<double>&, Polytope<double>&);

} // namespace cutcells::quadrays
