// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// The front end of part/ on [-1, 1]^3 meshed with hexahedra or Kuhn tetrahedra:
//  - bisect_tetrahedron reproduces the polynomial on both halves;
//  - classification by bounds finds every cell the sphere cuts, for Pk and
//    analytic level sets, and only few others;
//  - faces lying in a zero set are owned once, by the negative side, so plane
//    parts through mesh faces integrate exactly;
//  - parts of the sphere, and of a sphere and a plane, against exact totals;
//  - two planes crossing in cells: their corner, union and faces, exact.
// Exits non-zero on failure.

#include <cutcells/bernstein.h>
#include <cutcells/level_set.h>
#include <cutcells/part/classify.h>
#include <cutcells/part/cut_result.h>
#include <cutcells/part/mesh_part.h>
#include <cutcells/part/output.h>
#include <cutcells/quadrays/analytic.h>

#include <cmath>
#include <cstdio>
#include <functional>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "../quadrays/support/exact_reference.h"
#include "../quadrays/support/sphere_functors.h"
#include "support/box_mesh.h"

using namespace cutcells;
using namespace cutcells::part;
using namespace cutcells::part::support;

namespace
{

int failures = 0;

void check(bool ok, const std::string& what)
{
    if (!ok)
    {
        std::printf("FAILED: %s\n", what.c_str());
        ++failures;
    }
}

const V3 centre = {0.0123, -0.0371, 0.0217};
const double radius = 0.7;

/// The plane x_axis = offset.
struct AxisPlane
{
    int axis = 0;
    double offset = 0;

    template <typename V>
    V operator()(const std::array<V, 3>& x) const
    {
        return x[axis] - offset;
    }
};

/// A Pk level set interpolating @p f on @p mesh.
LevelSetFunction<double, int> pk_level_set(const MeshView<double, int>& mesh, int degree,
                                           const std::function<double(const V3&)>& f, const std::string& name)
{
    LevelSetMeshData<double, int> data = create_level_set_mesh_data<double, int>(mesh, degree);
    std::vector<double> values(static_cast<std::size_t>(data.num_dofs()));
    for (int d = 0; d < data.num_dofs(); ++d)
    {
        const double* x = data.dof_coordinate(d);
        values[static_cast<std::size_t>(d)] = f({x[0], x[1], x[2]});
    }
    return create_level_set_function<double, int>(std::move(data), std::span<const double>(values), name);
}

template <typename F>
LevelSetFunction<double, int> analytic_ls(const F& functor, const std::string& name)
{
    auto phi = std::make_shared<const quadrays::AnalyticLevelSet>(quadrays::analytic_level_set(functor));
    return create_level_set_function<double, int>(phi, 3, name);
}

double total(const MeshPart<double, int>& part, int order, bool full)
{
    const auto rules = quadrature_rules(part, order, full);
    double s = 0;
    for (const double w : rules._weights)
        s += w;
    return s;
}

void test_bisection()
{
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const std::array<std::array<double, 3>, 4> ref = {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    double worst = 0;
    for (int n = 1; n <= 4; ++n)
    {
        std::vector<double> c(static_cast<std::size_t>(bernstein::num_polynomials(cell::type::tetrahedron, n)));
        for (double& x : c)
            x = 2 * uniform(rng) - 1;
        std::vector<double> a, b;
        for (int p = 0; p < 4; ++p)
            for (int q = p + 1; q < 4; ++q)
            {
                bisect_tetrahedron<double>(c, n, p, q, a, b);
                for (int half = 0; half < 2; ++half)
                {
                    auto v = ref;
                    for (int i = 0; i < 3; ++i)
                        (half == 0 ? v[q] : v[p])[i] = 0.5 * (ref[p][i] + ref[q][i]);
                    const std::vector<double>& h = half == 0 ? a : b;
                    for (int s = 0; s < 20; ++s)
                    {
                        // a point of the half by its own reference coordinates
                        double l[4] = {uniform(rng), uniform(rng), uniform(rng), uniform(rng)};
                        const double sum = l[0] + l[1] + l[2] + l[3];
                        std::array<double, 3> xi_local = {l[1] / sum, l[2] / sum, l[3] / sum}, x = {0, 0, 0};
                        for (int r = 0; r < 4; ++r)
                            for (int i = 0; i < 3; ++i)
                                x[i] += l[r] / sum * v[r][i];
                        const double on_half = bernstein::evaluate<double>(cell::type::tetrahedron, n, h, xi_local);
                        const double original = bernstein::evaluate<double>(cell::type::tetrahedron, n, c, x);
                        worst = std::max(worst, std::abs(on_half - original));
                    }
                }
            }
    }
    check(worst < 1e-13, "bisect_tetrahedron misses the polynomial by " + std::to_string(worst));
}

/// Classification of the sphere against the cells it cuts exactly.
void test_classification()
{
    const quadrays::support::SphereDistance distance = {centre, radius};
    auto quadratic = [](const V3& x)
    { return std::pow(x[0] - centre[0], 2) + std::pow(x[1] - centre[1], 2) + std::pow(x[2] - centre[2], 2) - radius * radius; };
    for (const char* kind : {"hex", "tet"})
    {
        BoxMesh mesh;
        make_box_mesh(kind, 8, centre, mesh);
        std::vector<char> exact(mesh.cells.size());
        int n_exact = 0;
        for (std::size_t c = 0; c < mesh.cells.size(); ++c)
        {
            exact[c] = quadrays::support::exact::sphere_area(mesh.cells[c].faces, radius) > 0;
            n_exact += exact[c];
        }
        for (int source = 0; source < 2; ++source)
        {
            const std::vector<LevelSetFunction<double, int>> ls
                = {source == 0 ? pk_level_set(mesh.view, 2, quadratic, "phi") : analytic_ls(distance, "phi")};
            const CutResult<double, int> r = cut<double, int>(mesh.view, ls);
            int missed = 0, extra = 0, wrong_side = 0;
            for (std::size_t c = 0; c < mesh.cells.size(); ++c)
            {
                const cell::domain d = r.domain(0, static_cast<int>(c));
                const bool cut_cell = d == cell::domain::intersected;
                missed += exact[c] && !cut_cell;
                extra += !exact[c] && cut_cell;
                if (!cut_cell)
                {
                    const bool inside = quadratic(mesh.cells[c].centroid) < 0;
                    wrong_side += inside != (d == cell::domain::inside);
                }
            }
            std::printf("classification %s %s: %d cut cells, %d missed, %d extra, %d on the wrong side\n", kind,
                        source == 0 ? "P2" : "distance", n_exact, missed, extra, wrong_side);
            check(missed == 0 && wrong_side == 0 && extra <= n_exact / 50,
                  std::string("classification of the sphere, ") + kind);
        }
    }
}

/// The plane x = 0.25 lies in mesh faces (h = 0.25).
void test_zero_faces()
{
    const AxisPlane plane = {0, 0.25};
    for (const char* kind : {"hex", "tet"})
    {
        BoxMesh mesh;
        make_box_mesh(kind, 8, {0, 0, 0}, mesh);
        for (int source = 0; source < 2; ++source)
        {
            const std::vector<LevelSetFunction<double, int>> ls
                = {source == 0 ? pk_level_set(mesh.view, 1, [](const V3& x) { return x[0] - 0.25; }, "phi")
                               : analytic_ls(plane, "phi")};
            const CutResult<double, int> r = cut<double, int>(mesh.view, ls);
            const int expected = std::string(kind) == "hex" ? 64 : 128;
            bool owners_negative = true;
            for (int z = 0; z < r.n_zero_faces(); ++z)
                owners_negative &= r.domain(0, r.zero_face_cells[static_cast<std::size_t>(z)]) == cell::domain::inside;
            const double area = total(select(r, "phi = 0"), 3, false);
            const double below = total(select(r, "phi < 0"), 3, true), above = total(select(r, "phi > 0"), 3, true);
            std::printf("zero faces %s %s: %d faces (%d expected), interface %.15f, volumes %.15f %.15f\n", kind,
                        source == 0 ? "P1" : "analytic", r.n_zero_faces(), expected, area, below, above);
            check(r.n_zero_faces() == expected && owners_negative && r.cut_cells.empty(),
                  std::string("zero faces of the plane, ") + kind);
            // sums over thousands of cells: exact up to their rounding
            check(std::abs(area - 4) < 1e-12 * 4 && std::abs(below - 5) < 1e-12 * 5 && std::abs(above - 3) < 1e-12 * 3,
                  std::string("plane parts through faces, ") + kind);
        }
    }
}

/// Parts of the sphere, and of a smaller sphere with a plane beside it.
void test_parts()
{
    const quadrays::support::SphereDistance distance = {centre, radius};
    const double ball = 4.0 / 3.0 * M_PI * radius * radius * radius, sphere = 4.0 * M_PI * radius * radius;
    const V3 c2 = {-0.45, 0.02, -0.03};
    const double r2 = 0.5;
    const quadrays::support::SphereDistance small = {c2, r2};
    const AxisPlane plane = {0, 0.3};
    for (const char* kind : {"hex", "tet"})
    {
        BoxMesh mesh;
        make_box_mesh(kind, 8, centre, mesh);
        const std::vector<LevelSetFunction<double, int>> one = {analytic_ls(distance, "phi")};
        const CutResult<double, int> r = cut<double, int>(mesh.view, one);
        const double volume = total(select(r, "phi < 0"), 5, true), outside = total(select(r, "phi > 0"), 5, true);
        const double area = total(select(r, "phi = 0"), 5, false), both = total(select(r, "phi < 0 or phi > 0"), 5, true);
        std::printf("sphere %s: volume %.1e, outside %.1e, area %.1e, both sides %.1e\n", kind, volume / ball - 1,
                    outside / (8 - ball) - 1, area / sphere - 1, both / 8 - 1);
        check(std::abs(volume / ball - 1) < 1e-7 && std::abs(outside / (8 - ball) - 1) < 1e-7
                  && std::abs(area / sphere - 1) < 1e-6 && std::abs(both / 8 - 1) < 1e-12,
              std::string("parts of the sphere, ") + kind);

        const std::vector<LevelSetFunction<double, int>> two = {analytic_ls(small, "phi1"), analytic_ls(plane, "phi2")};
        const CutResult<double, int> r12 = cut<double, int>(mesh.view, two);
        const double b2 = 4.0 / 3.0 * M_PI * r2 * r2 * r2, s2 = 4.0 * M_PI * r2 * r2, slab = 4 * 1.3;
        const double in = total(select(r12, "phi1 < 0 and phi2 < 0"), 5, true);
        const double out = total(select(r12, "phi1 > 0 and phi2 < 0"), 5, true);
        const double face = total(select(r12, "phi1 = 0 and phi2 < 0"), 5, false);
        const double cut_plane = total(select(r12, "phi2 = 0 and phi1 > 0"), 5, false);
        std::printf("sphere and plane %s: ball %.1e, slab without ball %.1e, sphere %.1e, plane %.1e\n", kind,
                    in / b2 - 1, out / (slab - b2) - 1, face / s2 - 1, cut_plane / 4 - 1);
        check(std::abs(in / b2 - 1) < 1e-7 && std::abs(out / (slab - b2) - 1) < 1e-7 && std::abs(face / s2 - 1) < 1e-6
                  && std::abs(cut_plane / 4 - 1) < 1e-13,
              std::string("parts of a sphere and a plane, ") + kind);

        // two planes crossing in cells: x < 0.3 and y < 0.1, their union and
        // the faces of the corner, exact up to rounding
        const AxisPlane plane_y = {1, 0.1};
        const std::vector<LevelSetFunction<double, int>> crossing = {analytic_ls(plane, "phi1"),
                                                                     analytic_ls(plane_y, "phi2")};
        const CutResult<double, int> rc = cut<double, int>(mesh.view, crossing);
        const double corner = total(select(rc, "phi1 < 0 and phi2 < 0"), 3, true);
        const double either = total(select(rc, "phi1 < 0 or phi2 < 0"), 3, true);
        const double face1 = total(select(rc, "phi1 = 0 and phi2 < 0"), 3, false);
        const double face2 = total(select(rc, "phi2 = 0 and phi1 < 0"), 3, false);
        const double faces = total(select(rc, "phi1 = 0 and phi2 < 0 or phi2 = 0 and phi1 < 0"), 3, false);
        std::printf("two planes %s: corner %.1e, union %.1e, faces %.1e %.1e, both %.1e\n", kind, corner / 2.86 - 1,
                    either / 6.74 - 1, face1 / 2.2 - 1, face2 / 2.6 - 1, faces / 4.8 - 1);
        // sums over thousands of cells: exact up to their rounding
        check(std::abs(corner / 2.86 - 1) < 1e-12 && std::abs(either / 6.74 - 1) < 1e-12
                  && std::abs(face1 / 2.2 - 1) < 1e-12 && std::abs(face2 / 2.6 - 1) < 1e-12
                  && std::abs(faces / 4.8 - 1) < 1e-12,
              std::string("two planes crossing in cells, ") + kind);
    }
}

} // namespace

int main()
{
    test_bisection();
    test_classification();
    test_zero_faces();
    test_parts();
    return failures == 0 ? 0 : 1;
}
