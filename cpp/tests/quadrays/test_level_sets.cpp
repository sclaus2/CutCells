// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Several level sets in one cell (append_rules with one source per level set),
// as Bernstein coefficients and as analytic level sets.
//  - A ball and a half-space against exact per-cell values on hexahedra and
//    tetrahedra of [-1, 1]^3 (n = 8): the ball above the plane, the sphere
//    above it, the disk the plane cuts from the ball, and the union of ball and
//    half-space; a disk and a half-plane likewise on quadrilaterals and
//    triangles of [-1, 1]^2 (n = 16). The per-cell L1, relative to the total,
//    stays below the bounds given for q = 3, 5, 8.
//  - Exact totals of a lens (two balls) and a Steinmetz bicylinder (two
//    cylinders) on n = 8 with q = 5: volumes of intersections and unions, and
//    interfaces bounded by the other level set, to 1e-6 (the bicylinder 1e-5:
//    its cylinders touch where their intersection curves cross).
//  - Robustness, q = 3 on n = 8: a plane through vertices with a plane 1e-3
//    off the tetrahedra's slanted faces; a sphere with its tangent plane; two
//    balls touching at an off-grid point; a sphere through 12 vertices with a
//    plane in grid faces. Every rule must have finite points and weights, no
//    negative weight, no point outside its cell or outside its part (beyond
//    1e-9 h of a zero set), and at most Options::max_bisections bisections; the
//    per-cell L1 must stay below 1e-3 of the exact total (spheres of radius 2h
//    and 1.2h on hexahedra reach 2e-4 at q = 3, as with one level set), and
//    parts of measure zero below 1e-12 in total.
// Exits non-zero on failure.

#include <cutcells/quadrays/analytic.h>
#include <cutcells/quadrays/rules.h>
#include <cutcells/quadrays/taylor.h>
#include <cutcells/selection_expr.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <functional>
#include <string>
#include <type_traits>
#include <vector>

#include "support/exact_reference.h"
#include "support/test_mesh.h"

using namespace cutcells;
using namespace cutcells::quadrays;
using namespace cutcells::quadrays::support;

namespace
{
/// A level set: ball (or disk, with c[2] = 0), plane (or line) or cylinder.
struct Field
{
    enum class Kind
    {
        ball,
        plane,
        cylinder
    };
    Kind kind = Kind::ball;
    V3 c = {0, 0, 0};
    double r = 0;
    V3 n = {0, 0, 1}; ///< plane: n . x - d
    double d = 0;
    int axis = 2; ///< cylinder

    template <typename V>
    V operator()(const std::array<V, 3>& x) const
    {
        switch (kind)
        {
        case Kind::plane:
            return x[0] * n[0] + x[1] * n[1] + x[2] * n[2] - d;
        case Kind::cylinder:
        {
            const int i = (axis + 1) % 3, j = (axis + 2) % 3;
            const V a = x[i] - c[i], b = x[j] - c[j];
            return a * a + b * b - r * r;
        }
        default:
        {
            const V a = x[0] - c[0], b = x[1] - c[1], e = x[2] - c[2];
            return a * a + b * b + e * e - r * r;
        }
        }
    }

    int degree() const { return kind == Kind::plane ? 1 : 2; }
};

Field ball(const V3& c, double r)
{
    Field f;
    f.c = c;
    f.r = r;
    return f;
}

Field plane(const V3& n, double d)
{
    Field f;
    f.kind = Field::Kind::plane;
    f.n = n;
    f.d = d;
    return f;
}

Field cylinder(const V3& c, double r, int axis)
{
    Field f;
    f.kind = Field::Kind::cylinder;
    f.c = c;
    f.r = r;
    f.axis = axis;
    return f;
}

double value(const Field& f, const V3& x) { return f(std::array<double, 3>{x[0], x[1], x[2]}); }

double gradient_norm(const Field& f, const V3& x)
{
    std::array<Dual<double, 3>, 3> X;
    for (int k = 0; k < 3; ++k)
    {
        X[k] = Dual<double, 3>(x[k]);
        X[k].d[k] = 1;
    }
    const Dual<double, 3> g = f(X);
    return std::hypot(g.d[0], g.d[1], g.d[2]);
}

/// Does x lie in the part of one of the terms, up to 1e-9 h of each zero set?
bool in_part(const std::vector<Field>& fields, const std::vector<SelectionTerm>& terms, const V3& x, double h)
{
    for (const SelectionTerm& t : terms)
    {
        bool ok = true;
        for (std::size_t l = 0; l < fields.size() && ok; ++l)
        {
            const std::uint64_t bit = std::uint64_t(1) << l;
            const double f = value(fields[l], x), reach = 1e-9 * h * gradient_norm(fields[l], x);
            if (t.zero_required & bit)
                ok = std::abs(f) <= reach;
            else if (t.negative_required & bit)
                ok = f <= reach;
            else if (t.positive_required & bit)
                ok = f >= -reach;
        }
        if (ok)
            return true;
    }
    return false;
}

struct Run
{
    long fail = 0, negative = 0, outside = 0, side = 0;
    int max_bisections = 0;
    double total = 0, exact_total = 0, l1 = 0;

    long problems() const { return fail + negative + outside + side; }
};

/// The rules of the part @p text of the level sets @p fields (named a, b, c) on
/// every cell, their checks, and the per-cell errors if @p exact (per cell) is
/// given.
template <typename Cell>
Run run(const std::vector<Field>& fields, const char* text, const std::vector<Cell>& cells, double h, int q,
        bool analytic, const std::vector<double>& exact)
{
    constexpr bool planar = std::is_same_v<Cell, TestCell2D>;
    constexpr int tdim = planar ? 2 : 3;
    const std::vector<std::string> names = {"a", "b", "c"};
    SelectionExpr expr = parse_selection_expr(text);
    compile_selection_expr(expr, std::vector<std::string>(names.begin(), names.begin() + fields.size()));
    std::vector<AnalyticLevelSet> als;
    for (const Field& f : fields)
        als.push_back(analytic_level_set(f));
    std::vector<std::vector<double>> coeffs(fields.size());
    std::vector<BoxBernstein<double>> forms(fields.size());
    const Options opt;
    Run r;
    for (std::size_t i = 0; i < cells.size(); ++i)
    {
        const Cell& cell = cells[i];
        ClippedBox<double> box;
        make_clipped_box<double>(cell.type, cell.vertices, tdim, box);
        std::vector<Source<double>> sources;
        for (std::size_t l = 0; l < fields.size(); ++l)
        {
            if (analytic)
            {
                sources.push_back(analytic_source(als[l], box));
                continue;
            }
            const Field& f = fields[l];
            cell_coefficients(cell, f.degree(), [&f](const V3& x) { return value(f, x); }, coeffs[l]);
            cell_bernstein_on_box<double>(cell.type, f.degree(), coeffs[l], forms[l]);
            sources.push_back(bernstein_source(forms[l]));
        }
        quadrature::QuadratureRules<double> rule;
        Stats stats;
        try
        {
            append_rules<double>(box, std::span<const Source<double>>(sources),
                                 std::span<const SelectionTerm>(expr.terms), q, opt, 0, rule, stats);
        }
        catch (const std::exception&)
        {
            ++r.fail;
            continue;
        }
        r.max_bisections = std::max(r.max_bisections, stats.bisections);
        double sum = 0;
        for (std::size_t p = 0; p < rule._weights.size(); ++p)
        {
            const double w = rule._weights[p];
            const std::span<const double> xi(rule._points.data() + tdim * p, tdim);
            bool finite = std::isfinite(w);
            for (double x : xi)
                finite &= std::isfinite(x);
            if (!finite)
            {
                ++r.fail;
                continue;
            }
            sum += w;
            r.negative += w < 0;
            r.outside += !in_reference_cell(cell.type, xi, 1e-12);
            V3 x;
            if constexpr (planar)
                x = physical(cell, std::array<double, 2>{xi[0], xi[1]});
            else
                x = physical(cell, V3{xi[0], xi[1], xi[2]});
            r.side += !in_part(fields, expr.terms, x, h);
        }
        r.total += sum;
        if (!exact.empty())
        {
            r.exact_total += exact[i];
            r.l1 += std::abs(sum - exact[i]);
        }
    }
    return r;
}

std::vector<TestCell> mesh_3d(const std::string& mesh, int n, const V3& origin)
{
    const double h = 2.0 / n;
    std::vector<TestCell> cells;
    for (int i0 = 0; i0 < n; ++i0)
        for (int i1 = 0; i1 < n; ++i1)
            for (int i2 = 0; i2 < n; ++i2)
                for (const TestCell& cell : grid_cells(mesh, {-1 + h * i0, -1 + h * i1, -1 + h * i2}, h, origin))
                    cells.push_back(cell);
    return cells;
}

std::vector<TestCell2D> mesh_2d(const std::string& mesh, int n)
{
    const double h = 2.0 / n;
    std::vector<TestCell2D> cells;
    for (int i0 = 0; i0 < n; ++i0)
        for (int i1 = 0; i1 < n; ++i1)
            for (const TestCell2D& cell : grid_cells_2d(mesh, {-1 + h * i0, -1 + h * i1, 0}, h))
                cells.push_back(cell);
    return cells;
}

/// Volume of a convex polytope (empty: 0).
double volume(const std::vector<exact::Face>& faces)
{
    double v = 0;
    for (const exact::Face& f : faces)
        v += f.d * exact::loop_area(f.loop);
    return v / 3.0;
}

/// The volume and area of the ball of radius r about the origin in a polytope.
std::array<double, 2> ball_parts(const std::vector<exact::Face>& faces, double r)
{
    if (faces.empty())
        return {0.0, 0.0};
    const double area = exact::sphere_area(faces, r);
    return {exact::ball_volume(faces, r, area), area};
}

/// Faces moved so that @p c is the origin.
std::vector<exact::Face> shifted(std::vector<exact::Face> faces, const V3& c)
{
    for (exact::Face& f : faces)
    {
        for (exact::V3& v : f.loop)
            v = exact::sub(v, c);
        f.d -= exact::dot(f.n, c);
    }
    return faces;
}

/// The volume below the plane n . x = d in a polytope (empty: 0), and the area of the cut.
std::array<double, 2> below_plane(const std::vector<exact::Face>& faces, const V3& n, double d)
{
    if (faces.empty())
        return {0.0, 0.0};
    const exact::PlaneCut cut = exact::plane_cut(faces, n, d);
    return {cut.volume_below, cut.cut_area};
}

struct Target
{
    const char* text;
    std::vector<double> exact; ///< per cell
};

int failures = 0;

void report(const std::string& what, const Run& r, double l1, double bound)
{
    std::vector<std::string> problems;
    if (r.problems() > 0)
        problems.push_back("fail/neg/out/side " + std::to_string(r.fail) + "/" + std::to_string(r.negative) + "/"
                           + std::to_string(r.outside) + "/" + std::to_string(r.side));
    if (r.max_bisections > Options{}.max_bisections)
        problems.push_back(std::to_string(r.max_bisections) + " bisections in one cell");
    if (!(l1 <= bound))
        problems.push_back("error above the bound");
    std::string message;
    for (const std::string& p : problems)
        message += " " + p + ";";
    std::printf("%s %.1e (bound %.1e), max bisections per cell %d%s%s\n", what.c_str(), l1, bound, r.max_bisections,
                problems.empty() ? "" : " FAILED:", message.c_str());
    failures += !problems.empty();
}
} // namespace

int main()
{
    // ---- a ball and a half-space, per cell ----
    {
        const V3 c = {0.0123, -0.0371, 0.0217};
        const double R = 0.8, a = 0.3, rho = std::sqrt(R * R - a * a);
        const std::vector<Field> fields = {ball(c, R), plane({0, 0, 1}, c[2] + a)};
        struct Bound
        {
            int q;
            double l1[4]; ///< by part
        };
        // about three times the largest L1 of the meshes and sources
        const Bound bounds[] = {{3, {8e-6, 4e-5, 2e-7, 2.5e-6}},
                                {5, {7e-9, 1.3e-7, 2e-11, 5e-9}},
                                {8, {1.5e-12, 1.5e-10, 1e-14, 6e-12}}};
        for (const std::string mesh : {"hex", "tet"})
        {
            const std::vector<TestCell> cells = mesh_3d(mesh, 8, c); // faces relative to c
            std::vector<Target> parts = {{"a < 0 and b > 0", {}}, {"a = 0 and b > 0", {}}, {"b = 0 and a < 0", {}},
                                       {"a < 0 or b > 0", {}}};
            for (const TestCell& cell : cells)
            {
                const auto above = ball_parts(exact::half_space(cell.faces, {0, 0, -1}, -a), R);
                const auto whole = ball_parts(cell.faces, R);
                const double up = cell.volume - exact::plane_cut(cell.faces, {0, 0, 1}, a).volume_below;
                parts[0].exact.push_back(above[0]);
                parts[1].exact.push_back(above[1]);
                parts[2].exact.push_back(exact::polygon_disk_area(exact::plane_section(cell.faces, {0, 0, 1}, a),
                                                                  {0, 0, 1}, {0, 0, a}, rho));
                parts[3].exact.push_back(whole[0] + up - above[0]);
            }
            for (const bool analytic : {false, true})
                for (const Bound& b : bounds)
                    for (std::size_t k = 0; k < parts.size(); ++k)
                    {
                        const Run r = run(fields, parts[k].text, cells, 0.25, b.q, analytic, parts[k].exact);
                        char what[160];
                        std::snprintf(what, sizeof what, "ball, plane %-9s %s q = %d %-16s L1", analytic ? "analytic" : "P2/P1",
                                      mesh.c_str(), b.q, parts[k].text);
                        report(what, r, r.l1 / r.exact_total, b.l1[k]);
                    }
        }
    }

    // ---- a disk and a half-plane, per cell ----
    {
        const V3 c = {0.0123, -0.0371, 0};
        const double r = 0.7;
        const V3 n = {0.3, 1, 0};
        const double d = 0.3 * c[0] + c[1] + 0.25;
        const std::vector<Field> fields = {ball(c, r), plane(n, d)};
        struct Bound
        {
            int q;
            double l1[4];
        };
        const Bound bounds[] = {{3, {3e-8, 3e-7, 1e-14, 2e-7}},
                                {5, {2e-12, 5e-11, 1e-14, 2.5e-10}},
                                {8, {2e-14, 3e-14, 1e-14, 3e-14}}};
        for (const std::string mesh : {"quad", "tri"})
        {
            const std::vector<TestCell2D> cells = mesh_2d(mesh, 16);
            std::vector<Target> parts = {{"a < 0 and b > 0", {}}, {"a = 0 and b > 0", {}}, {"b = 0 and a < 0", {}},
                                       {"a < 0 or b > 0", {}}};
            const V3 m = {-n[0], -n[1], 0};
            for (const TestCell2D& cell : cells)
            {
                const std::vector<exact::V3> above = exact::clip_polygon(cell.loop, m, -d);
                const double disk_above = exact::disk_polygon_area(above, c, r);
                parts[0].exact.push_back(disk_above);
                parts[1].exact.push_back(exact::circle_polygon_length(above, c, r));
                parts[2].exact.push_back(exact::line_length(cell.loop, n, d, c, r));
                parts[3].exact.push_back(exact::disk_polygon_area(cell.loop, c, r) + exact::polygon_area(above)
                                         - disk_above);
            }
            for (const bool analytic : {false, true})
                for (const Bound& b : bounds)
                    for (std::size_t k = 0; k < parts.size(); ++k)
                    {
                        const Run run_ = run(fields, parts[k].text, cells, 0.125, b.q, analytic, parts[k].exact);
                        char what[160];
                        std::snprintf(what, sizeof what, "disk, line  %-9s %s q = %d %-16s L1", analytic ? "analytic" : "P2/P1",
                                      mesh.c_str(), b.q, parts[k].text);
                        report(what, run_, run_.l1 / run_.exact_total, b.l1[k]);
                    }
        }
    }

    // ---- exact totals: a lens and a bicylinder ----
    {
        const V3 c = {0.0123, -0.0371, 0.0217};
        struct Shape
        {
            const char* name;
            std::vector<Field> fields;
            std::vector<std::pair<const char*, double>> parts; ///< text, exact total
            double bound;                                      ///< relative error of the totals
        };
        std::vector<Shape> shapes;
        {
            // two balls: the lens, its faces on either sphere, the union
            const double R1 = 0.6, R2 = 0.5, d = std::sqrt(0.55 * 0.55 + 0.05 * 0.05);
            const double V = M_PI * (R1 + R2 - d) * (R1 + R2 - d) * (d * d + 2 * d * (R1 + R2) - 3 * (R1 - R2) * (R1 - R2))
                             / (12 * d);
            const double h1 = (R2 - R1 + d) * (R2 + R1 - d) / (2 * d), h2 = (R1 - R2 + d) * (R1 + R2 - d) / (2 * d);
            const double V1 = 4.0 / 3 * M_PI * R1 * R1 * R1, V2 = 4.0 / 3 * M_PI * R2 * R2 * R2;
            shapes.push_back({"lens",
                              {ball({c[0] - 0.25, c[1], c[2]}, R1), ball({c[0] + 0.3, c[1] + 0.05, c[2]}, R2)},
                              {{"a < 0 and b < 0", V},
                               {"a = 0 and b < 0", 2 * M_PI * R1 * h1},
                               {"b = 0 and a < 0", 2 * M_PI * R2 * h2},
                               {"a < 0 or b < 0", V1 + V2 - V}},
                              1e-6});
        }
        {
            // two cylinders: the Steinmetz solid and its faces
            const double r = 0.6;
            shapes.push_back({"bicylinder",
                              {cylinder(c, r, 2), cylinder(c, r, 1)},
                              {{"a < 0 and b < 0", 16 * r * r * r / 3},
                               {"a = 0 and b < 0", 8 * r * r},
                               {"b = 0 and a < 0", 8 * r * r}},
                              1e-5}); // the cylinders touch where their intersection curves cross
        }
        for (const Shape& shape : shapes)
            for (const std::string mesh : {"hex", "tet"})
            {
                const std::vector<TestCell> cells = mesh_3d(mesh, 8, {0, 0, 0});
                for (const bool analytic : {false, true})
                    for (const auto& [text, exact_total] : shape.parts)
                    {
                        const Run r = run(shape.fields, text, cells, 0.25, 5, analytic, {});
                        char what[160];
                        std::snprintf(what, sizeof what, "%-10s %-9s %s q = 5 %-16s total error", shape.name,
                                      analytic ? "analytic" : "P2", mesh.c_str(), text);
                        report(what, r, std::abs(r.total / exact_total - 1), shape.bound);
                    }
            }
    }

    // ---- robustness ----
    {
        const double h16 = 2.0 / 16;
        const V3 off = {0.0123, -0.0371, 0.0217};
        struct Case
        {
            const char* name;
            std::vector<Field> fields;
            V3 origin; ///< of the cells' faces
            std::vector<const char*> parts;
            std::function<double(const TestCell&, int)> exact; ///< of part i in a cell
        };
        std::vector<Case> cases;
        {
            // through vertices, and 1e-3 off the tetrahedra's slanted faces x - y = 2 h16
            const V3 na = {1, 1, 1}, nb = {1, -1, 0};
            const double da = 2 * h16, db = 2 * h16 + 1e-3;
            cases.push_back({"planes",
                             {plane(na, da), plane(nb, db)},
                             {0, 0, 0},
                             {"a < 0 and b < 0", "a = 0 and b < 0", "b = 0 and a < 0"},
                             [=](const TestCell& cell, int i)
                             {
                                 if (i == 2)
                                     return below_plane(exact::half_space(cell.faces, na, da), nb, db)[1];
                                 return below_plane(exact::half_space(cell.faces, nb, db), na, da)[i];
                             }});
        }
        {
            // a sphere and its tangent plane at the top
            const double R = 0.5;
            cases.push_back({"sphere, tangent plane",
                             {ball(off, R), plane({0, 0, 1}, off[2] + R)},
                             off,
                             {"a < 0 and b < 0", "a = 0 and b < 0", "b = 0 and a < 0"},
                             [=](const TestCell& cell, int i)
                             { return i == 2 ? 0.0 : ball_parts(cell.faces, R)[i]; }});
        }
        {
            // two balls touching at an off-grid point
            const V3 ca = {off[0] - 0.3, off[1], off[2]}, cb = {off[0] + 0.3, off[1], off[2]};
            cases.push_back({"balls touching",
                             {ball(ca, 0.3), ball(cb, 0.3)},
                             {0, 0, 0},
                             {"a < 0 or b < 0", "a < 0 and b < 0", "a = 0 and b > 0"},
                             [=](const TestCell& cell, int i)
                             {
                                 const auto a = ball_parts(shifted(cell.faces, ca), 0.3);
                                 if (i == 1)
                                     return 0.0;
                                 if (i == 2)
                                     return a[1];
                                 return a[0] + ball_parts(shifted(cell.faces, cb), 0.3)[0];
                             }});
        }
        {
            // a sphere through 12 vertices, a plane in grid faces
            const double R = std::sqrt(32.0) * h16;
            cases.push_back({"sphere, plane in faces",
                             {ball({0, 0, 0}, R), plane({1, 0, 0}, 2 * h16)},
                             {0, 0, 0},
                             {"a < 0 and b < 0", "a = 0 and b < 0"},
                             [=](const TestCell& cell, int i)
                             { return ball_parts(exact::half_space(cell.faces, {1, 0, 0}, 2 * h16), R)[i]; }});
        }
        for (const Case& c : cases)
            for (const std::string mesh : {"hex", "tet"})
            {
                const std::vector<TestCell> cells = mesh_3d(mesh, 8, c.origin);
                for (std::size_t i = 0; i < c.parts.size(); ++i)
                {
                    std::vector<double> exact;
                    for (const TestCell& cell : cells)
                        exact.push_back(c.exact(cell, static_cast<int>(i)));
                    for (const bool analytic : {false, true})
                    {
                        const Run r = run(c.fields, c.parts[i], cells, 0.25, 3, analytic, exact);
                        const bool zero = r.exact_total == 0;
                        char what[160];
                        std::snprintf(what, sizeof what, "%-22s %-9s %s %-16s %s", c.name,
                                      analytic ? "analytic" : "P2/P1", mesh.c_str(), c.parts[i],
                                      zero ? "total" : "L1");
                        report(what, r, zero ? std::abs(r.total) : r.l1 / r.exact_total, zero ? 1e-12 : 1e-3);
                    }
                }
            }
    }
    return failures == 0 ? 0 : 1;
}
