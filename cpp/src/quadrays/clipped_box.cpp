// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "clipped_box.h"

#include <cmath>
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
                      ClippedBox<T>& box)
{
    if (gdim != 3)
        throw std::invalid_argument("quadrays: cells must lie in 3D (gdim = 3)");
    // Basix vertices spanning the box from vertex 0
    std::array<int, 3> axes;
    switch (cell_type)
    {
    case cell::type::tetrahedron:
        axes = {1, 2, 3};
        break;
    case cell::type::hexahedron:
        axes = {1, 2, 4};
        break;
    default:
        throw std::invalid_argument("quadrays: unsupported cell type "
                                    + cell::cell_type_to_str(cell_type));
    }
    if (static_cast<int>(vertex_coords.size()) != cell::get_num_vertices(cell_type) * gdim)
        throw std::invalid_argument("quadrays: wrong number of vertex coordinates");

    for (int i = 0; i < 3; ++i)
    {
        box.origin[i] = vertex_coords[i];
        for (int k = 0; k < 3; ++k)
            box.jacobian[i][k] = vertex_coords[axes[k] * gdim + i] - vertex_coords[i];
    }
    box.ref_origin = {0, 0, 0};
    box.ref_jacobian = identity<T>();
    box.clips.clear();
    if (cell_type == cell::type::tetrahedron)
        box.clips.push_back({{T(1), T(1), T(1)}, T(1)});
}

template <std::floating_point T>
T jacobian_determinant(const ClippedBox<T>& box)
{
    return det3(box.jacobian);
}

template <std::floating_point T>
Mat3<T> inverse_jacobian(const ClippedBox<T>& box)
{
    const Mat3<T>& a = box.jacobian;
    const T det = det3(a);
    if (det == T(0))
        throw std::runtime_error("quadrays: singular box map");
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
ClippedBox<T> sub_box(const ClippedBox<T>& box, const Vec3<T>& lo, const Vec3<T>& hi)
{
    ClippedBox<T> child;
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
    for (int k = 0; k < 3; ++k)
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

template void make_clipped_box<float>(cell::type, std::span<const float>, int, ClippedBox<float>&);
template float jacobian_determinant<float>(const ClippedBox<float>&);
template Mat3<float> inverse_jacobian<float>(const ClippedBox<float>&);
template Vec3<float> physical_point<float>(const ClippedBox<float>&, const Vec3<float>&);
template Vec3<float> reference_point<float>(const ClippedBox<float>&, const Vec3<float>&);
template ClippedBox<float> sub_box<float>(const ClippedBox<float>&, const Vec3<float>&, const Vec3<float>&);
template int longest_axis<float>(const ClippedBox<float>&);
template bool may_meet_clips<float>(const ClippedBox<float>&, const Vec3<float>&, const Vec3<float>&, float);
template bool inside_clips<float>(const ClippedBox<float>&, const Vec3<float>&, float);
template void clipped_polytope<float>(const ClippedBox<float>&, Polytope<float>&);

template void make_clipped_box<double>(cell::type, std::span<const double>, int, ClippedBox<double>&);
template double jacobian_determinant<double>(const ClippedBox<double>&);
template Mat3<double> inverse_jacobian<double>(const ClippedBox<double>&);
template Vec3<double> physical_point<double>(const ClippedBox<double>&, const Vec3<double>&);
template Vec3<double> reference_point<double>(const ClippedBox<double>&, const Vec3<double>&);
template ClippedBox<double> sub_box<double>(const ClippedBox<double>&, const Vec3<double>&, const Vec3<double>&);
template int longest_axis<double>(const ClippedBox<double>&);
template bool may_meet_clips<double>(const ClippedBox<double>&, const Vec3<double>&, const Vec3<double>&, double);
template bool inside_clips<double>(const ClippedBox<double>&, const Vec3<double>&, double);
template void clipped_polytope<double>(const ClippedBox<double>&, Polytope<double>&);

} // namespace cutcells::quadrays
