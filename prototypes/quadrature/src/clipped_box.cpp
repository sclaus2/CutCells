// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "clipped_box.h"

#include <cmath>
#include <stdexcept>

namespace cutcells::proto
{

namespace
{
double det3(const Mat3& a)
{
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
           - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
           + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

Vec3 affine(const Vec3& origin, const Mat3& m, const Vec3& u)
{
    Vec3 y = origin;
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            y[i] += m[i][k] * u[k];
    return y;
}

Mat3 identity()
{
    Mat3 m = {};
    for (int i = 0; i < 3; ++i)
        m[i][i] = 1.0;
    return m;
}

/// Planes of the clipped region as a . u <= b: the six box faces, then the clips.
std::vector<HalfSpace> region_planes(const ClippedBox& box)
{
    std::vector<HalfSpace> planes;
    for (int k = 0; k < 3; ++k)
    {
        HalfSpace lower, upper;
        lower.c[k] = -1.0;
        lower.d = 0.0;
        upper.c[k] = 1.0;
        upper.d = 1.0;
        planes.push_back(lower);
        planes.push_back(upper);
    }
    planes.insert(planes.end(), box.clips.begin(), box.clips.end());
    return planes;
}
} // namespace

ClippedBox hex_cell(const Vec3& lo, double h)
{
    ClippedBox box;
    box.origin = lo;
    for (int i = 0; i < 3; ++i)
        box.jacobian[i][i] = h;
    box.ref_jacobian = identity();
    return box;
}

ClippedBox tet_cell(const std::array<Vec3, 4>& X)
{
    ClippedBox box;
    box.origin = X[0];
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            box.jacobian[i][k] = X[k + 1][i] - X[0][i];
    box.ref_jacobian = identity();
    box.clips.push_back({{1.0, 1.0, 1.0}, 1.0});
    return box;
}

double jacobian_determinant(const ClippedBox& box) { return det3(box.jacobian); }

Mat3 inverse_jacobian(const ClippedBox& box)
{
    const Mat3& a = box.jacobian;
    const double det = det3(a);
    if (det == 0.0)
        throw std::runtime_error("inverse_jacobian: singular box map");
    Mat3 inv;
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

Vec3 physical_point(const ClippedBox& box, const Vec3& u) { return affine(box.origin, box.jacobian, u); }

Vec3 reference_point(const ClippedBox& box, const Vec3& u) { return affine(box.ref_origin, box.ref_jacobian, u); }

ClippedBox sub_box(const ClippedBox& box, const Vec3& lo, const Vec3& hi)
{
    ClippedBox child;
    child.origin = physical_point(box, lo);
    child.ref_origin = reference_point(box, lo);
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
        {
            child.jacobian[i][k] = box.jacobian[i][k] * (hi[k] - lo[k]);
            child.ref_jacobian[i][k] = box.ref_jacobian[i][k] * (hi[k] - lo[k]);
        }
    for (const HalfSpace& h : box.clips)
    {
        // c . (lo + D u') <= d  ->  (c o D) . u' <= d - c . lo
        HalfSpace g;
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

int longest_axis(const ClippedBox& box)
{
    int best = 0;
    double best_len = -1.0;
    for (int k = 0; k < 3; ++k)
    {
        double len = 0.0;
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

bool may_meet_clips(const ClippedBox& box, const Vec3& lo, const Vec3& hi, double tol)
{
    for (const HalfSpace& h : box.clips)
    {
        // smallest value of c . u over the sub-box is attained at a corner
        double vmin = 0.0;
        for (int k = 0; k < 3; ++k)
            vmin += h.c[k] * (h.c[k] >= 0.0 ? lo[k] : hi[k]);
        if (vmin > h.d + tol)
            return false;
    }
    return true;
}

bool inside_clips(const ClippedBox& box, const Vec3& u)
{
    for (const HalfSpace& h : box.clips)
        if (h.c[0] * u[0] + h.c[1] * u[1] + h.c[2] * u[2] > h.d)
            return false;
    return true;
}

Polytope clipped_polytope(const ClippedBox& box)
{
    const std::vector<HalfSpace> planes = region_planes(box);
    const int np = static_cast<int>(planes.size());
    const double tol = 1e-11;

    Polytope poly;
    std::vector<std::vector<int>> active; // planes active at each vertex
    for (int a = 0; a < np; ++a)
        for (int b = a + 1; b < np; ++b)
            for (int c = b + 1; c < np; ++c)
            {
                const Mat3 m = {planes[a].c, planes[b].c, planes[c].c};
                const double det = det3(m);
                if (std::abs(det) < 1e-14)
                    continue;
                Vec3 x;
                for (int col = 0; col < 3; ++col)
                {
                    Mat3 mc = m;
                    mc[0][col] = planes[a].d;
                    mc[1][col] = planes[b].d;
                    mc[2][col] = planes[c].d;
                    x[col] = det3(mc) / det;
                }
                bool feasible = true;
                for (const HalfSpace& h : planes)
                    if (h.c[0] * x[0] + h.c[1] * x[1] + h.c[2] * x[2] > h.d + tol)
                    {
                        feasible = false;
                        break;
                    }
                if (!feasible)
                    continue;
                bool duplicate = false;
                for (const Vec3& v : poly.vertices)
                    if (std::abs(v[0] - x[0]) + std::abs(v[1] - x[1]) + std::abs(v[2] - x[2]) < 1e-10)
                    {
                        duplicate = true;
                        break;
                    }
                if (duplicate)
                    continue;
                std::vector<int> act;
                for (int p = 0; p < np; ++p)
                {
                    const HalfSpace& h = planes[p];
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
        std::vector<int> face_edges;
        for (int e = 0; e < static_cast<int>(poly.edges.size()); ++e)
        {
            bool a_on = false, b_on = false;
            for (int v : on)
            {
                a_on |= v == poly.edges[e][0];
                b_on |= v == poly.edges[e][1];
            }
            if (a_on && b_on)
                face_edges.push_back(e);
        }
        const Vec3& c = planes[p].c;
        const double n = std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
        poly.face_normals.push_back({c[0] / n, c[1] / n, c[2] / n});
        poly.face_edges.push_back(face_edges);
    }
    return poly;
}

} // namespace cutcells::proto
