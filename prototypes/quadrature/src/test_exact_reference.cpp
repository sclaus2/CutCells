// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Checks the exact reference: a spherical cap against its closed form, and the
// per-cell values summed over hexahedral and Kuhn-tetrahedral meshes against the
// sphere's area and the ball's volume. Exits non-zero on failure.

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

    // per-cell values summed over meshes of [-1, 1]^3, sphere centre off the grid
    const V3 centre = {0.013, -0.021, 0.007};
    const int n = 13;
    const double h = 2.0 / n;
    double hex_area = 0, hex_volume = 0, tet_area = 0, tet_volume = 0;
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b)
            for (int d = 0; d < n; ++d)
            {
                const V3 lo = {-1 + h * a - centre[0], -1 + h * b - centre[1], -1 + h * d - centre[2]};
                const auto hex = box_faces(lo, h);
                const double ha = sphere_area(hex, r);
                hex_area += ha;
                hex_volume += ball_volume(hex, r, ha);
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
                    const auto tet = tet_faces(X);
                    const double ta = sphere_area(tet, r);
                    tet_area += ta;
                    tet_volume += ball_volume(tet, r, ta);
                } while (std::next_permutation(p.begin(), p.end()));
            }
    const double area = 4 * M_PI * r * r, volume = 4.0 / 3.0 * M_PI * r * r * r;
    check("hex mesh, summed areas", hex_area, area, 1e-12);
    check("hex mesh, summed volumes", hex_volume, volume, 1e-12);
    check("Kuhn tet mesh, summed areas", tet_area, area, 1e-12);
    check("Kuhn tet mesh, summed volumes", tet_volume, volume, 1e-12);
    return failures == 0 ? 0 : 1;
}
