// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Analytic level sets through quadrays' interface (analytic.h):
//  - Taylor models and their adapter enclose sampled values and derivatives,
//    for square roots (the term algoim's sqrt leaves out), exp, log, sin, cos,
//    division and the branch helpers;
//  - planes are integrated exactly, per cell to 1e-13;
//  - the signed-distance sphere |x - c| - r reproduces the polynomial sphere's
//    errors (test_sphere.cpp) on hexahedra and tetrahedra, and the quadratic
//    sphere as a functor gives the polynomial path's errors on hexahedra;
//  - a hand-written ShapeForest tape of the distance gives the functor's rules;
//  - around the centre of a distance, where only the value has bounds, the ball
//    is integrated as accurately as with the polynomial;
//  - cells are classified by the level set's own bounds, caps that contain no
//    vertex included.
// Exits non-zero on failure.

#include <cutcells/quadrays/adapters/shapeforest_tape.h>
#include <cutcells/quadrays/analytic.h>
#include <cutcells/quadrays/rules.h>
#include <cutcells/selection_expr.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "support/exact_reference.h"
#include "support/sphere_functors.h"
#include "support/test_mesh.h"

using namespace cutcells;
using namespace cutcells::quadrays;
using namespace cutcells::quadrays::support;

namespace
{

const V3 centre = {0.0123, -0.0371, 0.0217};
const double radius = 0.7;

struct Plane
{
    V3 a;
    double b = 0;

    template <typename V>
    V operator()(const std::array<V, 3>& x) const
    {
        return x[0] * a[0] + x[1] * a[1] + x[2] * a[2] - b;
    }
};

struct Wavy
{
    template <typename V>
    V operator()(const std::array<V, 3>& x) const
    {
        using std::cos;
        using std::exp;
        using std::log;
        using std::sin;
        return sin(3.0 * x[0]) * exp(0.5 * x[1]) + log(2.0 + x[2] * x[2]) - cos(x[0] * x[1]) / (1.5 + x[2] * x[2])
               - 0.3;
    }
};

struct Branches
{
    template <typename V>
    V operator()(const std::array<V, 3>& x) const
    {
        using std::abs;
        using std::max;
        using std::min;
        return abs(x[0]) + max(x[1], 0.5 * x[2]) - min(x[2], x[0] * x[1]) - 0.4;
    }
};

int failures = 0;

void check(bool ok, const std::string& what)
{
    if (!ok)
    {
        std::printf("FAILED: %s\n", what.c_str());
        ++failures;
    }
}

// ============================================================================
// Enclosures
// ============================================================================

/// Taylor models over parallelepipeds and box bounds enclose values and
/// derivatives sampled on a 9^3 grid.
void test_enclosures(const AnalyticLevelSet& phi, const std::string& name, const V3& centre_x, double size,
                     int expect_status)
{
    const int G = 9;
    const double axes_all[9] = {0.8 * size, 0.3 * size, -0.2 * size, //
                                0.1 * size, 0.7 * size, 0.4 * size,  //
                                -0.3 * size, 0.2 * size, 0.9 * size};
    for (int m = 1; m <= 3; ++m)
    {
        double axes[9];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < m; ++j)
                axes[i * m + j] = axes_all[i * 3 + j];
        double models[4 * 5];
        const int status = parallelepiped_bounds(phi, centre_x.data(), axes, m, models);
        check(status == expect_status, name + ": Taylor models with status " + std::to_string(status));
        if (status == 0)
            continue;
        double worst = 0;
        const int total = m == 1 ? G : (m == 2 ? G * G : G * G * G);
        for (int n = 0; n < total; ++n)
        {
            double t[3] = {0, 0, 0};
            for (int j = 0, k = n; j < m; ++j, k /= G)
                t[j] = -1.0 + 2.0 * (k % G) / (G - 1);
            double x[3], g[3];
            for (int i = 0; i < 3; ++i)
            {
                x[i] = centre_x[i];
                for (int j = 0; j < m; ++j)
                    x[i] += axes[i * m + j] * t[j];
            }
            const double v = phi.gradient(x, g, phi.context);
            for (int row = 0; row <= (status == 1 ? m : 0); ++row)
            {
                // value, then d/dt_j = grad . axes_j
                double f = v;
                if (row > 0)
                {
                    f = 0;
                    for (int i = 0; i < 3; ++i)
                        f += g[i] * axes[i * m + row - 1];
                }
                const double* model = models + (m + 2) * row;
                double centre_value = model[0];
                for (int j = 0; j < m; ++j)
                    centre_value += model[1 + j] * t[j];
                const double excess = std::abs(f - centre_value) - model[m + 1];
                worst = std::max(worst, excess / (1.0 + std::abs(f)));
            }
        }
        check(worst <= 1e-12, name + ": Taylor model misses samples by " + std::to_string(worst) + " (m = "
                                  + std::to_string(m) + ")");
    }

    // box bounds over the bounding box of the parallelepiped with m = 3
    double lo[3], hi[3], b[8];
    for (int i = 0; i < 3; ++i)
    {
        double r = 0;
        for (int j = 0; j < 3; ++j)
            r += std::abs(axes_all[i * 3 + j]);
        lo[i] = centre_x[i] - r;
        hi[i] = centre_x[i] + r;
    }
    const int status = phi.box_bounds(lo, hi, b, phi.context);
    check(status == expect_status, name + ": box bounds with status " + std::to_string(status));
    if (status == 0)
        return;
    double worst = 0;
    for (int n = 0; n < G * G * G; ++n)
    {
        double x[3], g[3];
        for (int i = 0, k = n; i < 3; ++i, k /= G)
            x[i] = lo[i] + (hi[i] - lo[i]) * (k % G) / (G - 1);
        const double v = phi.gradient(x, g, phi.context);
        auto outside = [](double f, double l, double h) { return std::max({0.0, l - f, f - h}) / (1.0 + std::abs(f)); };
        worst = std::max(worst, outside(v, b[0], b[1]));
        for (int i = 0; i < 3 && status == 1; ++i)
            worst = std::max(worst, outside(g[i], b[2 + 2 * i], b[3 + 2 * i]));
    }
    check(worst <= 1e-12, name + ": box bounds miss samples by " + std::to_string(worst));
}

/// Values and gradients of the adapter against the functor and central differences.
template <typename F>
void test_point_values(const F& f, const AnalyticLevelSet& phi, const std::string& name)
{
    const double x[3] = {0.31, -0.42, 0.27};
    double g[3];
    const double v = phi.gradient(x, g, phi.context);
    const double expect = f(std::array<double, 3>{x[0], x[1], x[2]});
    const double tol = 1e-14 * (1 + std::abs(expect));
    check(std::abs(v - expect) <= tol && std::abs(phi.value(x, phi.context) - expect) <= tol, name + ": value");
    const double h = 1e-6;
    for (int i = 0; i < 3; ++i)
    {
        std::array<double, 3> xp = {x[0], x[1], x[2]}, xm = xp;
        xp[i] += h;
        xm[i] -= h;
        const double fd = (f(xp) - f(xm)) / (2 * h);
        check(std::abs(fd - g[i]) <= 1e-7 * (1 + std::abs(g[i])), name + ": gradient component " + std::to_string(i));
    }
}

// ============================================================================
// Integration
// ============================================================================

struct Errors
{
    double l1 = 0, worst = 0;
    int bisections = 0, uncertified = 0;
};

/// Per-cell errors of the sphere of radius 0.7 about centre, given by @p phi,
/// against exact values on [-1, 1]^3 with n^3 cubes or their Kuhn tetrahedra.
/// With coefficients_of non-null the cells get Bernstein forms of it instead.
Errors sphere_errors(const AnalyticLevelSet* phi, const std::string& mesh, int n, int q, const std::string& part_text,
                     const Options& opt = {})
{
    SelectionExpr expr = parse_selection_expr(part_text);
    compile_selection_expr(expr, {"phi"});
    const SelectionTerm& term = expr.terms.front();
    const bool surface = part_of(term) == Part::interface;
    const double h = 2.0 / n, r = radius, ball = 4.0 / 3.0 * M_PI * r * r * r, sphere = 4.0 * M_PI * r * r;
    auto quadratic = [](const V3& x)
    {
        return (x[0] - centre[0]) * (x[0] - centre[0]) + (x[1] - centre[1]) * (x[1] - centre[1])
               + (x[2] - centre[2]) * (x[2] - centre[2]) - radius * radius;
    };
    Errors e;
    Stats stats;
    double sum_abs = 0;
    std::vector<double> coeffs;
    for (int i0 = 0; i0 < n; ++i0)
        for (int i1 = 0; i1 < n; ++i1)
            for (int i2 = 0; i2 < n; ++i2)
            {
                const V3 lo = {-1 + h * i0, -1 + h * i1, -1 + h * i2};
                for (const TestCell& cell : grid_cells(mesh, lo, h, centre))
                {
                    double dc = 0;
                    for (int d = 0; d < 3; ++d)
                        dc += (cell.centroid[d] - centre[d]) * (cell.centroid[d] - centre[d]);
                    dc = std::sqrt(dc);
                    if (dc - cell.radius >= r || dc + cell.radius <= r)
                        continue;
                    const double area = exact::sphere_area(cell.faces, r);
                    if (area <= 0.0)
                        continue;
                    const double exact_value = surface ? area : exact::ball_volume(cell.faces, r, area);
                    quadrature::QuadratureRules<double> rules;
                    if (phi != nullptr)
                        append_cell_rules<double>(cell.type, cell.vertices, *phi, term, 0, q, opt, 0, rules, stats);
                    else
                    {
                        cell_coefficients(cell, 2, quadratic, coeffs);
                        append_cell_rules<double>(cell.type, cell.vertices, 2, coeffs, term, 0, q, opt, 0, rules,
                                                  stats);
                    }
                    double value = 0;
                    for (double w : rules._weights)
                        value += w;
                    sum_abs += std::abs(value - exact_value);
                    const double floor = surface ? 1e-3 * h * h : 1e-3 * h * h * h;
                    if (exact_value > floor)
                        e.worst = std::max(e.worst, std::abs(value - exact_value) / exact_value);
                }
            }
    e.l1 = sum_abs / (surface ? sphere : ball);
    e.bisections = stats.bisections;
    e.uncertified = stats.uncertified;
    return e;
}

void test_plane_exactness()
{
    const std::vector<Plane> planes = {{{1.0, 0.3, -0.2}, 0.0}, {{-0.45, 1.0, 0.7}, 0.1}};
    double worst_overall = 0;
    for (const Plane& plane : planes)
    {
        const AnalyticLevelSet phi = analytic_level_set(plane);
        for (const char* mesh : {"hex", "tet"})
            for (const char* part_text : {"phi < 0", "phi > 0", "phi = 0"})
            {
                SelectionExpr expr = parse_selection_expr(part_text);
                compile_selection_expr(expr, {"phi"});
                const SelectionTerm& term = expr.terms.front();
                const Part part = part_of(term);
                const int n = 5;
                const double h = 2.0 / n;
                const double cell_scale = part == Part::interface ? h * h : h * h * h;
                double total = 0, exact_total = 0, worst = 0;
                for (int i0 = 0; i0 < n; ++i0)
                    for (int i1 = 0; i1 < n; ++i1)
                        for (int i2 = 0; i2 < n; ++i2)
                        {
                            const V3 lo = {-1 + h * i0, -1 + h * i1, -1 + h * i2};
                            for (const TestCell& cell : grid_cells(mesh, lo, h, {0, 0, 0}))
                            {
                                const exact::PlaneCut cut = exact::plane_cut(cell.faces, plane.a, plane.b);
                                const double exact_value = part == Part::negative   ? cut.volume_below
                                                           : part == Part::positive ? cell.volume - cut.volume_below
                                                                                    : cut.cut_area;
                                quadrature::QuadratureRules<double> rules;
                                Stats stats;
                                append_cell_rules<double>(cell.type, cell.vertices, phi, term, 0, 3, Options{}, 0,
                                                          rules, stats);
                                double value = 0;
                                for (double w : rules._weights)
                                    value += w;
                                total += value;
                                exact_total += exact_value;
                                worst = std::max(worst, std::abs(value - exact_value) / cell_scale);
                            }
                        }
                const double rel_total = std::abs(total - exact_total) / exact_total;
                worst_overall = std::max({worst_overall, worst, rel_total});
                check(worst <= 1e-13 && rel_total <= 1e-13,
                      std::string("plane exactness, ") + mesh + ", " + part_text + ": worst cell "
                          + std::to_string(worst) + ", total " + std::to_string(rel_total));
            }
    }
    std::printf("planes: largest error %.1e (per cell relative to h^3 or h^2, totals relative)\n", worst_overall);

    // float: the interface is double, the engine float
    const Plane plane = planes.front();
    SelectionExpr expr = parse_selection_expr("phi < 0");
    compile_selection_expr(expr, {"phi"});
    const std::vector<float> vertices = {0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f};
    quadrature::QuadratureRules<float> rules;
    Stats stats;
    const Plane shifted = {plane.a, 0.25};
    const AnalyticLevelSet phi_shifted = analytic_level_set(shifted);
    append_cell_rules<float>(cell::type::tetrahedron, vertices, phi_shifted, expr.terms.front(), 0, 3, Options{}, 0,
                             rules, stats);
    double value = 0;
    for (float w : rules._weights)
        value += w;
    std::vector<exact::V3> tet = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    const double exact_value
        = exact::plane_cut(exact::tet_faces({tet[0], tet[1], tet[2], tet[3]}), shifted.a, shifted.b).volume_below;
    check(std::abs(value - exact_value) <= 1e-6, "plane exactness in float: " + std::to_string(value) + " against "
                                                    + std::to_string(exact_value));
}

/// cell_sign on cells around the sphere, and on the neighbours of the cell
/// [0, 0.5]^3 when the sphere about its centre reaches 0.01 beyond its faces:
/// caps in the middle of the neighbours' faces, away from all their vertices.
void test_cell_sign()
{
    auto sign_of = [](const AnalyticLevelSet& phi, cell::type type, const std::vector<double>& vertices)
    {
        ClippedBox<double> box;
        make_clipped_box<double>(type, vertices, 3, box);
        return cell_sign(box, phi);
    };
    auto hex = [](const V3& lo, double h)
    {
        std::vector<double> v;
        for (int k = 0; k < 8; ++k) // Basix order
            for (int d = 0; d < 3; ++d)
                v.push_back(lo[d] + h * ((k >> d) & 1));
        return v;
    };
    const SphereDistance distance = {centre, radius};
    const AnalyticLevelSet phi = analytic_level_set(distance);
    check(sign_of(phi, cell::type::hexahedron, hex({0.5, 0.5, 0.5}, 0.25)) == 1, "cell_sign: outside");
    check(sign_of(phi, cell::type::hexahedron, hex({0.0, 0.0, 0.0}, 0.25)) == -1, "cell_sign: inside");
    check(sign_of(phi, cell::type::hexahedron, hex({0.5, 0.0, 0.0}, 0.25)) == 0, "cell_sign: cut");

    const SphereDistance small = {{0.25, 0.25, 0.25}, 0.26};
    const AnalyticLevelSet caps = analytic_level_set(small);
    int cut_hexes = 0, cut_tets = 0, wrong = 0;
    for (const V3 lo : {V3{0.5, 0, 0}, V3{-0.5, 0, 0}, V3{0, 0.5, 0}, V3{0, -0.5, 0}, V3{0, 0, 0.5}, V3{0, 0, -0.5}})
    {
        const int s = sign_of(caps, cell::type::hexahedron, hex(lo, 0.5));
        cut_hexes += s == 0;
        wrong += s == -1;
        for (const TestCell& cell : grid_cells("tet", lo, 0.5, {0, 0, 0}))
        {
            const int t = sign_of(caps, cell::type::tetrahedron, cell.vertices);
            cut_tets += t == 0;
            wrong += t == -1;
        }
    }
    check(cut_hexes == 6 && wrong == 0, "cell_sign: caps in hexahedra (" + std::to_string(cut_hexes) + " of 6 cut)");
    check(cut_tets >= 6 && wrong == 0, "cell_sign: caps in tetrahedra (" + std::to_string(cut_tets) + " cut)");
}

/// A ShapeForest tape of |x - c| - r as ShapeForest lowers it.
shapeforest::Tape distance_tape(const V3& c, double r)
{
    using shapeforest::Op;
    shapeforest::Tape t;
    auto add = [&t](Op op, int a, int b, int cc, int out, double imm)
    {
        t.op.push_back(static_cast<std::uint8_t>(op));
        t.a.push_back(a);
        t.b.push_back(b);
        t.c.push_back(cc);
        t.out.push_back(out);
        t.imm.push_back(imm);
    };
    // registers 0..2: x, y, z; 3..5: centre; 6..8: differences; 9: length; 10: radius; 11: result
    t.n_registers = 12;
    t.inputs = {0, 1, 2};
    add(Op::constant, -1, -1, -1, 3, c[0]);
    add(Op::constant, -1, -1, -1, 4, c[1]);
    add(Op::constant, -1, -1, -1, 5, c[2]);
    add(Op::sub, 0, 3, -1, 6, 0);
    add(Op::sub, 1, 4, -1, 7, 0);
    add(Op::sub, 2, 5, -1, 8, 0);
    add(Op::length3, 6, 7, 8, 9, 0);
    add(Op::constant, -1, -1, -1, 10, r);
    add(Op::sub, 9, 10, -1, 11, 0);
    t.output = 11;
    shapeforest::prepare_tape(t);
    return t;
}

} // namespace

int main()
{
    const SphereDistance distance = {centre, radius};
    const SphereQuadratic quadratic = {centre, radius};
    const Wavy wavy;
    const Branches branches;
    const AnalyticLevelSet phi_distance = analytic_level_set(distance), phi_quadratic = analytic_level_set(quadratic),
                           phi_wavy = analytic_level_set(wavy), phi_branches = analytic_level_set(branches);

    // enclosures
    test_point_values(distance, phi_distance, "distance");
    test_point_values(wavy, phi_wavy, "wavy");
    test_point_values(branches, phi_branches, "branches");
    test_enclosures(phi_distance, "distance", {0.5, 0.4, -0.3}, 0.1, 1);
    test_enclosures(phi_distance, "distance at its centre", centre, 0.1, 2); // no gradient bound
    test_enclosures(phi_quadratic, "quadratic", {0.5, 0.4, -0.3}, 0.3, 1);
    test_enclosures(phi_wavy, "wavy", {0.2, -0.1, 0.3}, 0.2, 1);
    test_enclosures(phi_branches, "branches", {0.02, 0.1, 0.05}, 0.1, 1);
    test_enclosures(phi_branches, "branches, kinks inside", {0.0, 0.0, 0.0}, 0.2, 1);

    test_plane_exactness();
    test_cell_sign();

    // The distance against the polynomial sphere's numbers at n = 8
    // (test_sphere.cpp, the prototype's Bernstein path), 25% slack, with Taylor
    // models over whole boxes as in phase 2.
    struct Expected
    {
        const char* mesh;
        int q;
        const char* part;
        double l1, worst;
    };
    const Expected polynomial[] = {
        {"tet", 3, "phi < 0", 3.4e-7, 5.2e-5}, {"tet", 3, "phi = 0", 7.5e-6, 3.0e-4},
        {"tet", 5, "phi < 0", 1.1e-9, 3.9e-7}, {"tet", 5, "phi = 0", 5.5e-8, 3.9e-6},
        {"hex", 3, "phi < 0", 6.5e-6, 2.1e-4}, {"hex", 3, "phi = 0", 4.8e-5, 7.2e-4},
        {"hex", 5, "phi < 0", 2.2e-8, 9.4e-7}, {"hex", 5, "phi = 0", 3.7e-7, 8.6e-6},
    };
    Options whole_boxes;
    whole_boxes.taylor_subdivisions = 1;
    int whole_tet_bisections = 0;
    for (const Expected& p : polynomial)
    {
        const Errors e = sphere_errors(&phi_distance, p.mesh, 8, p.q, p.part, whole_boxes);
        const bool ok = e.l1 <= 1.25 * p.l1 && e.worst <= 1.25 * p.worst;
        std::printf("distance %s q = %d %-8s L1 %.1e (polynomial %.1e), worst %.1e (%.1e), bisections %d %s\n", p.mesh,
                    p.q, p.part, e.l1, p.l1, e.worst, p.worst, e.bisections, ok ? "" : "FAILED");
        failures += !ok;
        if (std::string(p.mesh) == "tet")
            whole_tet_bisections = e.bisections;
    }

    // Taylor models over the sub-boxes that meet a tetrahedron (the default):
    // at most half the whole boxes' bisections, and at q = 5 the errors of the
    // polynomial path in the same (orthogonal) frame, 25% slack. Whole boxes
    // bisect more, and so integrate more accurately at the same q.
    for (const Expected& p : polynomial)
    {
        if (std::string(p.mesh) != "tet" || p.q != 5)
            continue;
        const Errors bernstein = sphere_errors(nullptr, "tet", 8, 5, p.part);
        const Errors e = sphere_errors(&phi_distance, "tet", 8, 5, p.part);
        const bool ok = e.l1 <= 1.25 * bernstein.l1 && e.worst <= 1.25 * bernstein.worst
                        && 2 * e.bisections <= whole_tet_bisections;
        std::printf("distance tet sub-boxes q = 5 %-8s L1 %.1e (Bernstein %.1e), worst %.1e (%.1e), bisections %d "
                    "(whole boxes %d) %s\n",
                    p.part, e.l1, bernstein.l1, e.worst, bernstein.worst, e.bisections, whole_tet_bisections,
                    ok ? "" : "FAILED");
        failures += !ok;
    }

    // The quadratic sphere as a functor: the polynomial path's errors on hexahedra.
    for (const char* part : {"phi < 0", "phi = 0"})
    {
        const Errors a = sphere_errors(&phi_quadratic, "hex", 8, 5, part);
        const Errors b = sphere_errors(nullptr, "hex", 8, 5, part);
        const bool ok = std::abs(a.l1 - b.l1) <= 1e-2 * b.l1 && std::abs(a.worst - b.worst) <= 1e-2 * b.worst;
        std::printf("quadratic functor hex q = 5 %-8s L1 %.3e (Bernstein %.3e), worst %.3e (%.3e) %s\n", part, a.l1, b.l1,
                    a.worst, b.worst, ok ? "" : "FAILED");
        failures += !ok;
    }

    // The distance as a ShapeForest tape: the functor's rules.
    {
        const shapeforest::Tape tape = distance_tape(centre, radius);
        const shapeforest::TapeLevelSet tape_functor = {&tape};
        const AnalyticLevelSet phi_tape = analytic_level_set(tape_functor);
        for (const char* mesh : {"hex", "tet"})
            for (const char* part : {"phi < 0", "phi = 0"})
            {
                const Errors a = sphere_errors(&phi_tape, mesh, 8, 3, part);
                const Errors b = sphere_errors(&phi_distance, mesh, 8, 3, part);
                const bool ok = std::abs(a.l1 - b.l1) <= 1e-6 * b.l1 && a.bisections == b.bisections;
                std::printf("tape %s q = 3 %-8s L1 %.3e (functor %.3e), bisections %d (%d) %s\n", mesh, part, a.l1,
                            b.l1, a.bisections, b.bisections, ok ? "" : "FAILED");
                failures += !ok;
            }
    }

    // The distance about a mesh vertex, where its square root reaches 0: there
    // the value keeps a bound and the gradient loses its own. The ball is
    // integrated as accurately as the polynomial |x|^2 - r^2 on the same mesh.
    {
        const double r = 0.45;
        const SphereDistance on_vertex = {{0.0, 0.0, 0.0}, r};
        const AnalyticLevelSet phi = analytic_level_set(on_vertex);
        auto polynomial = [r](const V3& x) { return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - r * r; };
        SelectionExpr expr = parse_selection_expr("phi < 0");
        compile_selection_expr(expr, {"phi"});
        for (const char* mesh : {"hex", "tet"})
        {
            const int n = 4;
            const double h = 2.0 / n;
            double total[2] = {0, 0};
            bool positive = true;
            Stats stats[2];
            std::vector<double> coeffs;
            for (int i0 = 0; i0 < n; ++i0)
                for (int i1 = 0; i1 < n; ++i1)
                    for (int i2 = 0; i2 < n; ++i2)
                    {
                        const V3 lo = {-1 + h * i0, -1 + h * i1, -1 + h * i2};
                        for (const TestCell& cell : grid_cells(mesh, lo, h, {0, 0, 0}))
                        {
                            quadrature::QuadratureRules<double> rules[2];
                            append_cell_rules<double>(cell.type, cell.vertices, phi, expr.terms.front(), 0, 5,
                                                      Options{}, 0, rules[0], stats[0]);
                            cell_coefficients(cell, 2, polynomial, coeffs);
                            append_cell_rules<double>(cell.type, cell.vertices, 2, coeffs, expr.terms.front(), 0, 5,
                                                      Options{}, 0, rules[1], stats[1]);
                            for (int k = 0; k < 2; ++k)
                                for (double w : rules[k]._weights)
                                {
                                    total[k] += w;
                                    positive &= w > 0 && std::isfinite(w);
                                }
                        }
                    }
            const double exact_value = 4.0 / 3.0 * M_PI * r * r * r;
            const double error = std::abs(total[0] - exact_value) / exact_value;
            const double reference = std::abs(total[1] - exact_value) / exact_value;
            const bool ok = positive && error <= 1e-5 && error <= 3 * reference + 1e-9;
            std::printf("distance about a vertex, %s: volume error %.1e (polynomial %.1e), bisections %d (%d), "
                        "uncertified %d (%d) %s\n",
                        mesh, error, reference, stats[0].bisections, stats[1].bisections, stats[0].uncertified,
                        stats[1].uncertified, ok ? "" : "FAILED");
            failures += !ok;
        }
    }

    return failures == 0 ? 0 : 1;
}
