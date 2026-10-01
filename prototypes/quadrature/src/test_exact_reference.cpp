// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Checks the exact references: a spherical cap against its closed form, and the
// per-cell values summed over hexahedral and Kuhn-tetrahedral meshes against the
// sphere's area and the ball's volume, also for spheres through grid vertices and
// tangent to grid planes; plane cuts of a cube and of meshes against closed forms.
// Exits non-zero on failure.

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "exact_reference.h"

using namespace cutcells::proto::exact;

int main()
{
    const double r = 0.7;
    int failures = 0;
    auto check = [&](const char* what, double value, double expected, double tol)
    {
        const double rel = std::abs(value - expected) / std::abs(expected);
        std::printf("%-34s rel. error %.1e %s\n", what, rel, rel < tol ? "" : "FAILED");
        failures += rel >= tol;
    };

    // cap above z = 0.5 inside a box that contains it
    const auto cap_box = make_faces({{V3{-0.8, -0.8, 0.5}, V3{0.8, -0.8, 0.5}, V3{0.8, 0.8, 0.5}, V3{-0.8, 0.8, 0.5}},
                                     {V3{-0.8, -0.8, 0.8}, V3{0.8, -0.8, 0.8}, V3{0.8, 0.8, 0.8}, V3{-0.8, 0.8, 0.8}},
                                     {V3{-0.8, -0.8, 0.5}, V3{-0.8, 0.8, 0.5}, V3{-0.8, 0.8, 0.8}, V3{-0.8, -0.8, 0.8}},
                                     {V3{0.8, -0.8, 0.5}, V3{0.8, 0.8, 0.5}, V3{0.8, 0.8, 0.8}, V3{0.8, -0.8, 0.8}},
                                     {V3{-0.8, -0.8, 0.5}, V3{0.8, -0.8, 0.5}, V3{0.8, -0.8, 0.8}, V3{-0.8, -0.8, 0.8}},
                                     {V3{-0.8, 0.8, 0.5}, V3{0.8, 0.8, 0.5}, V3{0.8, 0.8, 0.8}, V3{-0.8, 0.8, 0.8}}});
    const double cap_area = sphere_area(cap_box, r);
    check("cap area", cap_area, 2 * M_PI * r * 0.2, 1e-13);
    check("cap volume", ball_volume(cap_box, r, cap_area), M_PI * 0.04 * (3 * r - 0.2) / 3, 1e-13);

    // Per-cell values summed over meshes of [-1, 1]^3 with n cells per side; f is
    // called with the faces of every hex and Kuhn tet (relative to origin).
    auto for_cells = [](int n, const V3& origin, auto&& f)
    {
        const double h = 2.0 / n;
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b)
                for (int d = 0; d < n; ++d)
                {
                    const V3 lo = {-1 + h * a - origin[0], -1 + h * b - origin[1], -1 + h * d - origin[2]};
                    f(box_faces(lo, h), true);
                    std::array<int, 3> p = {0, 1, 2};
                    do
                    {
                        std::array<V3, 4> X;
                        X[0] = lo;
                        for (int k = 0; k < 3; ++k)
                        {
                            X[k + 1] = X[k];
                            X[k + 1][p[k]] += h;
                        }
                        f(tet_faces(X), false);
                    } while (std::next_permutation(p.begin(), p.end()));
                }
    };
    auto sphere_sums = [&](const char* name, int n, const V3& centre, double radius)
    {
        double hex_area = 0, hex_volume = 0, tet_area = 0, tet_volume = 0;
        for_cells(n, centre,
                  [&](const std::vector<Face>& faces, bool hex)
                  {
                      const double a = sphere_area(faces, radius), v = ball_volume(faces, radius, a);
                      (hex ? hex_area : tet_area) += a;
                      (hex ? hex_volume : tet_volume) += v;
                  });
        const double area = 4 * M_PI * radius * radius, volume = 4.0 / 3.0 * M_PI * radius * radius * radius;
        char what[128];
        std::snprintf(what, sizeof what, "%s, hex areas", name);
        check(what, hex_area, area, 1e-12);
        std::snprintf(what, sizeof what, "%s, hex volumes", name);
        check(what, hex_volume, volume, 1e-12);
        std::snprintf(what, sizeof what, "%s, tet areas", name);
        check(what, tet_area, area, 1e-12);
        std::snprintf(what, sizeof what, "%s, tet volumes", name);
        check(what, tet_volume, volume, 1e-12);
    };
    sphere_sums("off the grid", 13, {0.013, -0.021, 0.007}, r);
    // n = 16, h = 0.125: through 6 vertices and tangent to grid planes there; through
    // 12 vertices; tangent to the plane x = 0.75 inside a face; a cap of height 1e-6
    sphere_sums("r = 4h, centre on a vertex", 16, {0, 0, 0}, 0.5);
    sphere_sums("r = sqrt(32) h, centre on a vertex", 16, {0, 0, 0}, std::sqrt(0.5));
    sphere_sums("tangent to x = 0.75", 16, {0.0123, -0.0371, 0.0217}, 0.75 - 0.0123);
    sphere_sums("cap of height 1e-6", 16, {0.0123, -0.0371, 0.0217}, 0.75 - 0.0123 + 1e-6);

    // plane cuts of the unit cube
    const auto cube = box_faces({0, 0, 0}, 1.0);
    const PlaneCut c1 = plane_cut(cube, {1, 0, 0}, 0.3);
    check("cube cut by x = 0.3, volume", c1.volume_below, 0.3, 1e-14);
    check("cube cut by x = 0.3, area", c1.cut_area, 1.0, 1e-14);
    const PlaneCut c2 = plane_cut(cube, {1, 1, 1}, 1.5);
    check("cube cut by x + y + z = 1.5, volume", c2.volume_below, 0.5, 1e-14);
    check("cube cut by x + y + z = 1.5, area", c2.cut_area, 0.75 * std::sqrt(3.0), 1e-14);
    const PlaneCut c3 = plane_cut(cube, {1, 0, 0}, 1.0);
    check("cube with face in x = 1, volume", c3.volume_below, 1.0, 1e-14);
    check("cube with face in x = 1, face", c3.face_in_plane, 1.0, 1e-14);

    // plane cuts summed over meshes (n = 16): through vertices, and on faces
    auto plane_sums = [&](const char* name, const V3& a, double b, double volume, double area)
    {
        double hv = 0, ha = 0, tv = 0, ta = 0;
        for_cells(16, {0, 0, 0},
                  [&](const std::vector<Face>& faces, bool hex)
                  {
                      const PlaneCut c = plane_cut(faces, a, b);
                      (hex ? hv : tv) += c.volume_below;
                      (hex ? ha : ta) += c.cut_area + 0.5 * c.face_in_plane; // a face in the plane is shared
                  });
        char what[128];
        std::snprintf(what, sizeof what, "%s, hex volumes", name);
        check(what, hv, volume, 1e-12);
        std::snprintf(what, sizeof what, "%s, hex areas", name);
        check(what, ha, area, 1e-12);
        std::snprintf(what, sizeof what, "%s, tet volumes", name);
        check(what, tv, volume, 1e-12);
        std::snprintf(what, sizeof what, "%s, tet areas", name);
        check(what, ta, area, 1e-12);
    };
    // x + y + z < t in [-1, 1]^3, 0 <= t <= 1: volume 8 (1/2 + (3t - t^3/3) / 8), area (3 - t^2) sqrt(3)
    const double t = 0.25;
    plane_sums("x + y + z = 0.25", {1, 1, 1}, t, 8 * (0.5 + (3 * t - t * t * t / 3) / 8), (3 - t * t) * std::sqrt(3.0));
    plane_sums("x = 0.25 (on faces)", {1, 0, 0}, 0.25, 5.0, 4.0);
    // x - y > 0.125: a triangular prism with legs 1.875 and length 2
    plane_sums("x - y = 0.125 (on tet faces)", {1, -1, 0}, 0.125, 8.0 - 1.875 * 1.875, 2.0 * std::sqrt(2.0) * 1.875);
    return failures == 0 ? 0 : 1;
}
