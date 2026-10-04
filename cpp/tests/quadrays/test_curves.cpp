// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// The curves where two level sets vanish ("a = 0 and b = 0"; points on 2D
// cells), with one source per level set, as Bernstein coefficients and as
// analytic level sets, on hexahedra and tetrahedra of [-1, 1]^3 (n = 8) and
// quadrilaterals and triangles of [-1, 1]^2 (n = 16):
//  - two planes: their line, per cell against its length in the cell, to
//    rounding, and its first moments;
//  - the rim of a lens of two balls, one circle, per cell against the arc in
//    the cell, for q = 3, 5, 8;
//  - the rims of a napkin ring (a ball outside a cylinder), two circles, per
//    cell, 4 pi r in all;
//  - the edges of a Steinmetz bicylinder, two ellipses crossing where the
//    cylinders touch, 8 sqrt(2) r E(1 / sqrt(2)) in all (E by the AGM);
//  - corners: those edges inside and outside a third cylinder (the
//    tricylinder's edges), 8 r int_0^{pi/4} sqrt(1 + cos^2 t) dt inside, the
//    two parts adding up to the edges;
//  - degenerate placements: a small ball about a cell's centre with a plane
//    through it, and two small balls symmetric about a cell's mid-plane (the
//    curve lies where the cell's box is bisected);
//  - on 2D cells, two disks crossing at two points, each of weight 1 at its
//    exact place; three disks, the crossings kept by the sign of the third;
//    two lines crossing in a cell, and in a vertex, where every cell sharing
//    it finds the point at most once; the degenerate placements of 3D;
//  - leaves: Lagrange curves (VTK 68) in 3D, vertices (VTK 1) on 2D cells,
//    every node on both zero sets up to 1e-5 h.
// Every rule must have finite points and weights, no negative weight, no point
// outside its cell or off the curve (beyond 1e-9 h of either zero set), and
// Stats::curve_lost must stay 0.
// The argument "bernstein" or "analytic" runs one kind of level set only.
// Exits non-zero on failure.

#include <cutcells/quadrays/rules.h>
#include <cutcells/selection_expr.h>

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "support/exact_reference.h"
#include "support/fields.h"
#include "support/test_mesh.h"

using namespace cutcells;
using namespace cutcells::quadrays;
using namespace cutcells::quadrays::support;

namespace
{

int failures = 0;

void report(const std::string& what, const Run& r, double error, double bound)
{
    std::vector<std::string> problems;
    if (r.problems() > 0)
        problems.push_back("fail/neg/out/side " + std::to_string(r.fail) + "/" + std::to_string(r.negative) + "/"
                           + std::to_string(r.outside) + "/" + std::to_string(r.side));
    if (r.max_bisections > Options{}.max_bisections)
        problems.push_back(std::to_string(r.max_bisections) + " bisections in one cell");
    if (r.curve_lost > 0)
        problems.push_back(std::to_string(r.curve_lost) + " boxes lost the curve");
    if (!(error <= bound))
        problems.push_back("error above the bound");
    std::string message;
    for (const std::string& p : problems)
        message += " " + p + ";";
    std::printf("%s %.1e (bound %.1e), max bisections per cell %d%s%s\n", what.c_str(), error, bound,
                r.max_bisections, problems.empty() ? "" : " FAILED:", message.c_str());
    failures += !problems.empty();
}

void check(bool ok, const std::string& what)
{
    std::printf("%s%s\n", what.c_str(), ok ? "" : " FAILED");
    failures += !ok;
}

/// Length of the line p + s t inside a convex polytope; its parameters there
/// go to [s_lo, s_hi].
double line_length_in(const std::vector<exact::Face>& faces, const V3& p, const V3& t, double& s_lo, double& s_hi)
{
    s_lo = -1e300;
    s_hi = 1e300;
    for (const exact::Face& f : faces)
    {
        const double ft = exact::dot(f.n, t), rest = f.d - exact::dot(f.n, p);
        if (std::abs(ft) < 1e-300)
        {
            if (rest < 0)
                return 0.0;
            continue;
        }
        if (ft > 0)
            s_hi = std::min(s_hi, rest / ft);
        else
            s_lo = std::max(s_lo, rest / ft);
    }
    return s_hi > s_lo ? (s_hi - s_lo) * exact::norm(t) : 0.0;
}

/// Length of the circle of radius rho about @p centre, in the plane through it
/// with unit normal @p e, inside a convex polytope.
double circle_length_in(const std::vector<exact::Face>& faces, const V3& centre, const V3& e, double rho)
{
    const std::vector<exact::V3> section = exact::plane_section(faces, e, exact::dot(e, centre));
    if (section.size() < 3)
        return 0.0;
    const V3 t0 = std::abs(e[0]) < 0.9 ? V3{1, 0, 0} : V3{0, 1, 0};
    const V3 e1 = exact::unit(exact::cross(e, t0)), e2 = exact::cross(e, e1);
    std::vector<exact::V3> loop;
    for (const exact::V3& v : section)
    {
        const V3 w = exact::sub(v, centre);
        loop.push_back({exact::dot(w, e1), exact::dot(w, e2), 0});
    }
    return exact::circle_polygon_length(loop, {0, 0, 0}, rho);
}

/// The complete elliptic integral of the second kind E(m), m = k^2, by the
/// arithmetic-geometric mean.
double elliptic_e(double m)
{
    double a = 1, b = std::sqrt(1 - m), c = std::sqrt(m), sum = 0.5 * c * c, power = 0.5;
    for (int it = 0; it < 30 && c > 1e-17; ++it)
    {
        const double an = 0.5 * (a + b);
        c = 0.5 * (a - b);
        b = std::sqrt(a * b);
        a = an;
        power *= 2;
        sum += power * c * c;
    }
    return M_PI / (2 * a) * (1 - sum);
}

/// The points where the circles |x - ca| = ra and |x - cb| = rb cross (z = 0).
std::vector<V3> circle_crossings(const V3& ca, double ra, const V3& cb, double rb)
{
    const double dx = cb[0] - ca[0], dy = cb[1] - ca[1], d = std::hypot(dx, dy);
    const double x1 = (d * d + ra * ra - rb * rb) / (2 * d), y1 = std::sqrt(ra * ra - x1 * x1);
    const V3 e = {dx / d, dy / d, 0};
    return {V3{ca[0] + x1 * e[0] - y1 * e[1], ca[1] + x1 * e[1] + y1 * e[0], 0},
            V3{ca[0] + x1 * e[0] + y1 * e[1], ca[1] + x1 * e[1] - y1 * e[0], 0}};
}

/// Do the points of @p r, each of weight 1, match @p exact one to one within
/// @p tol?
bool match_points(const Run& r, const std::vector<V3>& exact, double tol)
{
    if (r.points.size() != exact.size())
        return false;
    std::vector<char> used(exact.size(), 0);
    for (std::size_t p = 0; p < r.points.size(); ++p)
    {
        if (r.weights[p] != 1.0)
            return false;
        bool found = false;
        for (std::size_t e = 0; e < exact.size() && !found; ++e)
            if (!used[e] && exact::norm(exact::sub(r.points[p], exact[e])) <= tol)
                found = used[e] = 1;
        if (!found)
            return false;
    }
    return true;
}

std::string source_name(bool analytic) { return analytic ? "analytic" : "P2/P1"; }

} // namespace

int main(int argc, char** argv)
{
    const std::string only = argc > 1 ? argv[1] : "";
    std::vector<bool> sources; // analytic or not
    if (only != "analytic")
        sources.push_back(false);
    if (only != "bernstein")
        sources.push_back(true);
    const V3 c = {0.0123, -0.0371, 0.0217};
    const char* curve = "a = 0 and b = 0";
    char what[200];

    // ---- per cell: two planes, the lens's rim, the napkin ring's rims ----
    {
        // two planes in general position, their line p0 + s t
        const V3 na = {1.0, 0.3, -0.2}, nb = {-0.25, 1.0, 0.4};
        const double da = 0.1, db = -0.05;
        const V3 t = exact::cross(na, nb);
        const double tt = exact::dot(t, t);
        const V3 p0 = exact::mul(1.0 / tt, exact::add(exact::mul(da, exact::cross(nb, t)),
                                                       exact::mul(db, exact::cross(t, na))));
        // two balls: the rim is the circle where the spheres meet
        const V3 ca = {c[0] - 0.25, c[1], c[2]}, cb = {c[0] + 0.3, c[1] + 0.05, c[2]};
        const double R1 = 0.6, R2 = 0.5;
        const double d = exact::norm(exact::sub(cb, ca)), x1 = (d * d + R1 * R1 - R2 * R2) / (2 * d);
        const V3 e = exact::mul(1.0 / d, exact::sub(cb, ca)), rim = exact::add(ca, exact::mul(x1, e));
        const double rho = std::sqrt(R1 * R1 - x1 * x1);
        // a ball outside a cylinder: rims of radius r at heights +- h / 2
        const double R = 0.8, r = 0.45, half = std::sqrt(R * R - r * r);
        struct Shape
        {
            const char* name;
            std::vector<Field> fields;
            std::function<double(const TestCell&)> exact;
            std::vector<std::pair<int, double>> bounds; ///< q, per-cell L1 relative to the total
        };
        const std::vector<Shape> shapes = {
            {"two planes",
             {plane(na, da), plane(nb, db)},
             [&](const TestCell& cell)
             {
                 double lo = 0, hi = 0;
                 return line_length_in(cell.faces, p0, t, lo, hi);
             },
             {{1, 5e-15}, {3, 5e-15}}},
            {"lens rim",
             {ball(ca, R1), ball(cb, R2)},
             [&](const TestCell& cell) { return circle_length_in(cell.faces, rim, e, rho); },
             {{3, 4e-7}, {5, 6e-11}, {8, 1e-14}}},
            {"napkin rims",
             {ball(c, R), cylinder(c, r, 2)},
             [&](const TestCell& cell)
             {
                 return circle_length_in(cell.faces, {c[0], c[1], c[2] + half}, {0, 0, 1}, r)
                        + circle_length_in(cell.faces, {c[0], c[1], c[2] - half}, {0, 0, 1}, r);
             },
             {{3, 6e-5}, {5, 2e-7}, {8, 6e-11}}}};
        for (const std::string mesh : {"hex", "tet"})
        {
            const std::vector<TestCell> cells = mesh_3d(mesh, 8, {0, 0, 0});
            for (const Shape& shape : shapes)
            {
                std::vector<double> exact;
                for (const TestCell& cell : cells)
                    exact.push_back(shape.exact(cell));
                for (const bool analytic : sources)
                    for (const auto& [q, bound] : shape.bounds)
                    {
                        const Run run_ = run(shape.fields, curve, cells, 0.25, q, analytic, exact);
                        std::snprintf(what, sizeof what, "%-11s %-9s %s q = %d per-cell L1", shape.name,
                                      source_name(analytic).c_str(), mesh.c_str(), q);
                        report(what, run_, run_.l1 / run_.exact_total, bound);
                    }
            }
            // the line's length in the box and its first moments
            double length = 0;
            V3 moment = {0, 0, 0};
            for (const TestCell& cell : cells)
            {
                double lo = 0, hi = 0;
                const double l = line_length_in(cell.faces, p0, t, lo, hi);
                length += l;
                moment = exact::add(moment, exact::mul(l, exact::add(p0, exact::mul(0.5 * (lo + hi), t))));
            }
            for (const bool analytic : sources)
            {
                const Run run_ = run(shapes[0].fields, curve, cells, 0.25, 2, analytic, {});
                V3 m = {0, 0, 0};
                for (std::size_t p = 0; p < run_.points.size(); ++p)
                    m = exact::add(m, exact::mul(run_.weights[p], run_.points[p]));
                const double error = std::max({std::abs(run_.total - length), std::abs(m[0] - moment[0]),
                                               std::abs(m[1] - moment[1]), std::abs(m[2] - moment[2])})
                                     / length;
                std::snprintf(what, sizeof what, "two planes  %-9s %s q = 2 length and first moments",
                              source_name(analytic).c_str(), mesh.c_str());
                report(what, run_, error, 2e-15);
            }
        }
    }

    // ---- totals: the bicylinder's edges, and the tricylinder's (corners) ----
    {
        const double r = 0.6;
        const double edges = 8 * std::sqrt(2.0) * r * elliptic_e(0.5);
        const auto speed = [](double s) { return std::sqrt(1 + std::cos(s) * std::cos(s)); };
        const double inside = 8 * r * exact::tanh_sinh(speed, 0.0, M_PI / 4);
        const std::vector<Field> two = {cylinder(c, r, 2), cylinder(c, r, 1)};
        const std::vector<Field> three = {cylinder(c, r, 2), cylinder(c, r, 1), cylinder(c, r, 0)};
        for (const std::string mesh : {"hex", "tet"})
        {
            const std::vector<TestCell> cells = mesh_3d(mesh, 8, {0, 0, 0});
            for (const bool analytic : sources)
            {
                const Run all = run(two, curve, cells, 0.25, 5, analytic, {});
                std::snprintf(what, sizeof what, "bicylinder  %-9s %s q = 5 edges, total error",
                              source_name(analytic).c_str(), mesh.c_str());
                report(what, all, std::abs(all.total / edges - 1), 3e-7);
                const Run in = run(three, "a = 0 and b = 0 and c < 0", cells, 0.25, 5, analytic, {});
                std::snprintf(what, sizeof what, "tricylinder %-9s %s q = 5 edges in c, total error",
                              source_name(analytic).c_str(), mesh.c_str());
                report(what, in, std::abs(in.total / inside - 1), 2e-12);
                const Run out = run(three, "a = 0 and b = 0 and c > 0", cells, 0.25, 5, analytic, {});
                std::snprintf(what, sizeof what, "tricylinder %-9s %s q = 5 edges out of c, total error",
                              source_name(analytic).c_str(), mesh.c_str());
                report(what, out, std::abs(out.total / (edges - inside) - 1), 1e-9);
                std::snprintf(what, sizeof what, "tricylinder %-9s %s q = 5 in and out add up to the edges",
                              source_name(analytic).c_str(), mesh.c_str());
                Run both = in;
                accumulate(both, out);
                report(what, both, std::abs((in.total + out.total) / all.total - 1), 3e-7);
            }
        }
    }

    // ---- degenerate placements: the curve in a cell's mid-plane ----
    {
        const double h = 0.25;
        const V3 m = {0.125, 0.125, 0.125}; // the centre of a cell of n = 8
        const double rho = 0.1, offset = 0.05, rho2 = std::sqrt(rho * rho - offset * offset);
        struct Case
        {
            const char* name;
            std::vector<Field> fields;
            double exact;
            double bound; ///< of the total's relative error (a circle 0.4 h across resolved at q = 5)
        };
        const std::vector<Case> cases = {
            {"ball, plane through its centre", {ball(m, rho), plane({1, 0, 0}, m[0])}, 2 * M_PI * rho, 5e-5},
            {"balls symmetric about a mid-plane",
             {ball({m[0] - offset, m[1], m[2]}, rho), ball({m[0] + offset, m[1], m[2]}, rho)},
             2 * M_PI * rho2,
             1e-11}};
        for (const Case& cs : cases)
            for (const std::string mesh : {"hex", "tet"})
            {
                const std::vector<TestCell> cells = mesh_3d(mesh, 8, {0, 0, 0});
                for (const bool analytic : sources)
                {
                    const Run r = run(cs.fields, curve, cells, h, 5, analytic, {});
                    std::snprintf(what, sizeof what, "%-34s %-9s %s q = 5 total error", cs.name,
                                  source_name(analytic).c_str(), mesh.c_str());
                    report(what, r, std::abs(r.total / cs.exact - 1), cs.bound);
                }
            }
    }

    // ---- 2D: crossing points ----
    {
        const double h = 0.125;
        const V3 c2 = {c[0], c[1], 0};
        // two disks
        const V3 ca = {c2[0] - 0.2, c2[1], 0}, cb = {c2[0] + 0.3, c2[1] + 0.1, 0};
        const double r1 = 0.62, r2 = 0.5;
        const std::vector<Field> two = {ball(ca, r1), ball(cb, r2)};
        const std::vector<V3> crossings = circle_crossings(ca, r1, cb, r2);
        // three disks: the crossings of a and b on either side of c
        const V3 cc = {c2[0] + 0.02, c2[1] + 0.3, 0};
        const double r3 = 0.48;
        const std::vector<Field> three = {ball(ca, r1), ball(cb, r2), ball(cc, r3)};
        std::vector<V3> in_c, out_c;
        for (const V3& x : crossings)
            (value(three[2], x) < 0 ? in_c : out_c).push_back(x);
        // two lines crossing in a cell, and in a vertex; through the vertex both
        // slopes lie in (0, 1), so that both lines enter the cells north-east and
        // south-west of it (squares and triangles, which split along x = y).
        // Where they enter different cells, each only touches the others at the
        // vertex and no cell counts as cut by both: the point is not found.
        const V3 la = {1.0, 0.37, 0}, lb = {-0.6, 1.0, 0}, va = {0.37, -1.0, 0}, vb = {0.8, -1.0, 0};
        const V3 p_cell = {0.0123 + 0.03, -0.0371 + 0.02, 0}, p_vertex = {0.125, -0.25, 0};
        auto lines = [](const V3& na, const V3& nb, const V3& p)
        { return std::vector<Field>{plane(na, exact::dot(na, p)), plane(nb, exact::dot(nb, p))}; };
        // degenerate: the crossings on a cell's mid-line
        const V3 m = {0.0625, 0.0625, 0};
        const double rho = 0.04, offset = 0.02, y2 = std::sqrt(rho * rho - offset * offset);
        for (const std::string mesh : {"quad", "tri"})
        {
            const std::vector<TestCell2D> cells = mesh_2d(mesh, 16);
            for (const bool analytic : sources)
            {
                const std::string tag = source_name(analytic) + " " + mesh;
                const Run r = run(two, curve, cells, h, 3, analytic, {});
                report("two disks " + tag + ": 2 points of weight 1 at the crossings", r,
                       match_points(r, crossings, 1e-12) ? 0.0 : 1.0, 0.0);
                const Run rin = run(three, "a = 0 and b = 0 and c < 0", cells, h, 3, analytic, {});
                report("three disks " + tag + ": crossings of a and b in c", rin,
                       match_points(rin, in_c, 1e-12) ? 0.0 : 1.0, 0.0);
                const Run rout = run(three, "a = 0 and b = 0 and c > 0", cells, h, 3, analytic, {});
                report("three disks " + tag + ": crossings of a and b out of c", rout,
                       match_points(rout, out_c, 1e-12) ? 0.0 : 1.0, 0.0);
                const Run rl = run(lines(la, lb, p_cell), curve, cells, h, 3, analytic, {});
                report("two lines crossing in a cell " + tag, rl, match_points(rl, {p_cell}, 1e-12) ? 0.0 : 1.0, 0.0);
                // in a vertex: each cell sharing it finds the point at most once
                int found = 0, most = 0;
                bool placed = true;
                Run all;
                for (const TestCell2D& cell : cells)
                {
                    const Run one
                        = run(lines(va, vb, p_vertex), curve, std::vector<TestCell2D>{cell}, h, 3, analytic, {});
                    found += static_cast<int>(one.points.size());
                    most = std::max(most, static_cast<int>(one.points.size()));
                    for (std::size_t p = 0; p < one.points.size(); ++p)
                        placed &= exact::norm(exact::sub(one.points[p], p_vertex)) <= 1e-12 && one.weights[p] == 1.0;
                    accumulate(all, one);
                }
                std::snprintf(what, sizeof what,
                              "two lines crossing in a vertex %s: found by %d cells, at most once each", tag.c_str(),
                              found);
                report(what, all, found >= 1 && most == 1 && placed ? 0.0 : 1.0, 0.0);
                const std::vector<Field> centred = {ball(m, rho), plane({1, 0, 0}, m[0])};
                const Run rc = run(centred, curve, cells, h, 3, analytic, {});
                report("disk, line through its centre " + tag, rc,
                       match_points(rc, {V3{m[0], m[1] - rho, 0}, V3{m[0], m[1] + rho, 0}}, 1e-12) ? 0.0 : 1.0, 0.0);
                const std::vector<Field> symmetric = {ball({m[0] - offset, m[1], 0}, rho),
                                                      ball({m[0] + offset, m[1], 0}, rho)};
                const Run rs = run(symmetric, curve, cells, h, 3, analytic, {});
                report("disks symmetric about a mid-line " + tag, rs,
                       match_points(rs, {V3{m[0], m[1] - y2, 0}, V3{m[0], m[1] + y2, 0}}, 1e-12) ? 0.0 : 1.0, 0.0);
            }
        }
    }

    // ---- leaves ----
    {
        const V3 ca = {c[0] - 0.25, c[1], c[2]}, cb = {c[0] + 0.3, c[1] + 0.05, c[2]};
        for (const std::string mesh : {"hex", "tet"})
        {
            const std::vector<TestCell> cells = mesh_3d(mesh, 8, {0, 0, 0});
            for (const bool analytic : sources)
                for (const auto& [name, fields] :
                     {std::pair{"lens rim", std::vector<Field>{ball(ca, 0.6), ball(cb, 0.5)}},
                      std::pair{"bicylinder edges", std::vector<Field>{cylinder(c, 0.6, 2), cylinder(c, 0.6, 1)}}})
                {
                    const LeafRun r = leaves_of(fields, curve, cells, 0.25, analytic);
                    const bool ok = r.leaves > 0 && r.incomplete <= 0.03 * r.leaves && r.outside == 0
                                    && r.types.size() == 1 && r.types.count(vtk_lagrange_curve) == 1;
                    std::snprintf(what, sizeof what,
                                  "leaves %-16s %-9s %s: %ld Lagrange curves, %ld dropped, %ld nodes off the curve",
                                  name, source_name(analytic).c_str(), mesh.c_str(), r.leaves, r.incomplete, r.outside);
                    check(ok, what);
                }
        }
        const V3 c2 = {c[0], c[1], 0};
        for (const std::string mesh : {"quad", "tri"})
        {
            const std::vector<TestCell2D> cells = mesh_2d(mesh, 16);
            for (const bool analytic : sources)
            {
                const std::vector<Field> two = {ball({c2[0] - 0.2, c2[1], 0}, 0.62),
                                                ball({c2[0] + 0.3, c2[1] + 0.1, 0}, 0.5)};
                const LeafRun r = leaves_of(two, curve, cells, 0.125, analytic);
                const bool ok
                    = r.leaves == 2 && r.outside == 0 && r.types.size() == 1 && r.types.count(vtk_vertex) == 1;
                std::snprintf(what, sizeof what, "leaves two disks %-9s %s: %ld vertices, %ld nodes off the points",
                              source_name(analytic).c_str(), mesh.c_str(), r.leaves, r.outside);
                check(ok, what);
            }
        }
    }
    return failures == 0 ? 0 : 1;
}
