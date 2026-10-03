// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// quadrays on triangles and quadrilaterals of [-1, 1]^2, as Bernstein
// coefficients and as analytic level sets, against exact per-cell values.
//  - Lines (degree 1 and 2, and analytic) are exact on every cell, cut or not:
//    per cell relative to h^2 or h, and in total, to 1e-13.
//  - The circle of radius 0.7 off the grid's symmetry converges with q: on
//    n = 8, the per-cell L1 of the disk's area and of the circle's length,
//    relative to their totals, stays below the bounds given for q = 3, 5, 8.
//  - The robustness cases in 2D, on n = 16 with q = 3: circles through
//    vertices, tangent to edges, cutting caps of 1e-6 and 1e-12, scaled by
//    1e+-150 and 1e+-200; lines in edges, 1e-12 off them, through vertices; two
//    circles touching or 1e-6 and 1e-10 apart; an annulus 1e-3 wide (two roots
//    on a line); two lines crossing in a vertex and in a cell. Every rule must
//    have finite points and weights, no negative weight, no point outside its
//    cell or on the wrong side, and at most Options::max_bisections bisections;
//    scaled level sets must give the unscaled totals, and every case a per-cell
//    L1 below 1e-4 of the exact total.
// Exits non-zero on failure.

#include <cutcells/quadrays/analytic.h>
#include <cutcells/quadrays/rules.h>
#include <cutcells/quadrays/taylor.h>
#include <cutcells/selection_expr.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <map>
#include <string>
#include <vector>

#include "support/exact_reference.h"
#include "support/test_mesh.h"

using namespace cutcells;
using namespace cutcells::quadrays;
using namespace cutcells::quadrays::support;

namespace
{
enum class Shape
{
    circle,      ///< |x - centre|^2 - radius^2
    line,        ///< normal . x - offset
    two_circles, ///< product of two circles with disjoint disks
    annulus,     ///< (|x - centre|^2 - radius^2)(|x - centre|^2 - radius2^2), radius < radius2
    cross        ///< (x - cx)^2 - (y - cy)^2 / 4: two lines crossing at centre
};

/// A level set in the plane z = 0, multiplied by scale.
struct Case
{
    std::string name;
    Shape shape = Shape::circle;
    V3 centre = {0, 0, 0};
    double radius = 0;
    V3 centre2 = {0, 0, 0};
    double radius2 = 0;
    V3 normal = {1, 0, 0};
    double offset = 0;
    double scale = 1;
};

int degree(const Case& c)
{
    switch (c.shape)
    {
    case Shape::line:
        return 1;
    case Shape::two_circles:
    case Shape::annulus:
        return 4;
    default:
        return 2;
    }
}

template <typename V>
V case_value(const Case& c, const std::array<V, 3>& x)
{
    auto sq = [&x](const V3& p)
    {
        const V a = x[0] - p[0], b = x[1] - p[1];
        return a * a + b * b;
    };
    V v;
    switch (c.shape)
    {
    case Shape::circle:
        v = sq(c.centre) - c.radius * c.radius;
        break;
    case Shape::line:
        v = x[0] * c.normal[0] + x[1] * c.normal[1] - c.offset;
        break;
    case Shape::two_circles:
        v = (sq(c.centre) - c.radius * c.radius) * (sq(c.centre2) - c.radius2 * c.radius2);
        break;
    case Shape::annulus:
        v = (sq(c.centre) - c.radius * c.radius) * (sq(c.centre) - c.radius2 * c.radius2);
        break;
    case Shape::cross:
    {
        const V a = x[0] - c.centre[0], b = x[1] - c.centre[1];
        v = a * a - 0.25 * (b * b);
        break;
    }
    }
    return v * c.scale;
}

double phi(const Case& c, const V3& x) { return case_value(c, std::array<double, 3>{x[0], x[1], x[2]}); }

double gradient_norm(const Case& c, const V3& x)
{
    std::array<Dual<double, 3>, 3> X;
    for (int k = 0; k < 3; ++k)
    {
        X[k] = Dual<double, 3>(x[k]);
        X[k].d[k] = 1;
    }
    const Dual<double, 3> f = case_value(c, X);
    return std::hypot(f.d[0], f.d[1]); // no underflow for phi scaled by 1e-200
}

struct CaseLevelSet
{
    const Case* c = nullptr;

    template <typename V>
    V operator()(const std::array<V, 3>& x) const
    {
        return case_value(*c, x);
    }
};

std::vector<Case> all_cases()
{
    const V3 off = {0.0123, -0.0371, 0};
    const double h = 2.0 / 16;
    std::vector<Case> cases;
    auto circle = [&](const char* name, const V3& centre, double radius, double scale = 1.0)
    {
        Case c;
        c.name = name;
        c.centre = centre;
        c.radius = radius;
        c.scale = scale;
        cases.push_back(c);
    };
    auto line = [&](const char* name, const V3& normal, double offset)
    {
        Case c;
        c.name = name;
        c.shape = Shape::line;
        c.normal = normal;
        c.offset = offset;
        cases.push_back(c);
    };
    circle("circle", off, 0.7);
    circle("circle-vertex-tangent", {0, 0, 0}, 4 * h); // through 4 vertices, tangent to grid lines there
    circle("circle-vertices", {0, 0, 0}, 5 * h);       // through 12 vertices: (3h, 4h), ...
    circle("circle-edge-tangent", off, 6 * h - off[0]);
    circle("circle-cap-1e-6", off, 6 * h - off[0] + 1e-6);
    circle("circle-cap-1e-12", off, 6 * h - off[0] + 1e-12);
    circle("circle-scale-1e-150", off, 0.7, 1e-150);
    circle("circle-scale-1e+150", off, 0.7, 1e150);
    circle("circle-scale-1e-200", off, 0.7, 1e-200);
    circle("circle-scale-1e+200", off, 0.7, 1e200);
    line("line", {1, 0.3, 0}, 0.01);
    line("line-edges", {1, 0, 0}, 2 * h);
    line("line-near-edges", {1, 0, 0}, 2 * h + 1e-12);
    line("line-vertices", {1, 1, 0}, 2 * h);
    line("line-diagonals", {1, -1, 0}, h); // in the triangles' slanted edges
    for (const auto& [name, gap] : {std::pair{"circles-touching", 0.0}, std::pair{"circles-gap-1e-6", 1e-6},
                                    std::pair{"circles-gap-1e-10", 1e-10}})
    {
        Case c;
        c.name = name;
        c.shape = Shape::two_circles;
        c.centre = {off[0] - 0.3, off[1], 0};
        c.centre2 = {off[0] + 0.3 + gap, off[1], 0};
        c.radius = c.radius2 = 0.3;
        cases.push_back(c);
    }
    Case annulus;
    annulus.name = "annulus-1e-3";
    annulus.shape = Shape::annulus;
    annulus.centre = off;
    annulus.radius = 0.5;
    annulus.radius2 = 0.5 + 1e-3;
    cases.push_back(annulus);
    for (const auto& [name, centre] : {std::pair{"cross-vertex", V3{0, 0, 0}}, std::pair{"cross-cell", off}})
    {
        Case c;
        c.name = name;
        c.shape = Shape::cross;
        c.centre = centre;
        cases.push_back(c);
    }
    return cases;
}

/// Exact part measures of one cell.
struct Reference
{
    double area = 0;        ///< phi < 0
    double length = 0;      ///< phi = 0 through the cell
    double edge_length = 0; ///< phi = 0 on an edge of the cell (no owner)
};

Reference disk(const TestCell2D& cell, const V3& c, double r)
{
    return {exact::disk_polygon_area(cell.loop, c, r), exact::circle_polygon_length(cell.loop, c, r), 0.0};
}

/// n . x < d in the cell, and the line n . x = d.
Reference half_plane(const TestCell2D& cell, const V3& n, double d)
{
    Reference out;
    const std::size_t m = cell.loop.size();
    for (std::size_t i = 0; i < m; ++i)
    {
        const V3& a = cell.loop[i];
        const V3& b = cell.loop[(i + 1) % m];
        if (std::abs(exact::dot(n, a) - d) <= 1e-13 && std::abs(exact::dot(n, b) - d) <= 1e-13)
        {
            out.edge_length = exact::norm(exact::sub(b, a));
            out.area = exact::dot(n, cell.loop[(i + 2) % m]) < d ? cell.area : 0.0;
            return out;
        }
    }
    out.area = exact::polygon_area(exact::clip_polygon(cell.loop, n, d));
    out.length = exact::line_length(cell.loop, n, d, {0, 0, 0}, 0.0);
    return out;
}

Reference reference(const Case& c, const TestCell2D& cell)
{
    switch (c.shape)
    {
    case Shape::circle:
        return disk(cell, c.centre, c.radius);
    case Shape::line:
        return half_plane(cell, c.normal, c.offset);
    case Shape::two_circles:
    {
        const Reference a = disk(cell, c.centre, c.radius), b = disk(cell, c.centre2, c.radius2);
        return {a.area + b.area, a.length + b.length, 0.0};
    }
    case Shape::annulus:
    {
        const Reference a = disk(cell, c.centre, c.radius), b = disk(cell, c.centre, c.radius2);
        return {b.area - a.area, a.length + b.length, 0.0};
    }
    case Shape::cross:
    {
        // phi < 0 where |x - cx| < |y - cy| / 2: two wedges, above and below the centre
        const double cx = c.centre[0], cy = c.centre[1];
        const std::vector<exact::V3> up
            = exact::clip_polygon(exact::clip_polygon(cell.loop, {1, -0.5, 0}, cx - 0.5 * cy), {-1, -0.5, 0}, -cx - 0.5 * cy);
        const std::vector<exact::V3> down
            = exact::clip_polygon(exact::clip_polygon(cell.loop, {1, 0.5, 0}, cx + 0.5 * cy), {-1, 0.5, 0}, -cx + 0.5 * cy);
        Reference out;
        out.area = exact::polygon_area(up) + exact::polygon_area(down);
        out.length = exact::line_length(cell.loop, {1, -0.5, 0}, cx - 0.5 * cy, {0, 0, 0}, 0.0)
                     + exact::line_length(cell.loop, {1, 0.5, 0}, cx + 0.5 * cy, {0, 0, 0}, 0.0);
        return out;
    }
    }
    return {};
}

struct Run
{
    long fail = 0, negative = 0, outside = 0, side = 0;
    int max_bisections = 0;
    double total = 0, exact_total = 0, l1 = 0;
};

/// One case on the n x n grid of [-1, 1]^2: rules of every cell, their checks
/// and the per-cell errors.
Run run_case(const Case& c, const std::string& mesh, int n, const SelectionTerm& term, int q, const Options& opt,
             bool analytic)
{
    const double h = 2.0 / n;
    const bool surface = part_of(term) == Part::interface;
    const CaseLevelSet functor = {&c};
    const AnalyticLevelSet phi_analytic = analytic_level_set(functor);
    std::vector<double> coeffs;
    Run r;
    for (int i0 = 0; i0 < n; ++i0)
        for (int i1 = 0; i1 < n; ++i1)
            for (const TestCell2D& cell : grid_cells_2d(mesh, {-1 + h * i0, -1 + h * i1, 0}, h))
            {
                const Reference ref = reference(c, cell);
                r.exact_total += surface ? ref.length + 0.5 * ref.edge_length : ref.area;
                quadrature::QuadratureRules<double> rule;
                Stats stats;
                try
                {
                    if (analytic)
                        append_cell_rules<double>(cell.type, cell.vertices, phi_analytic, term, 0, q, opt, 0, rule,
                                                  stats);
                    else
                    {
                        cell_coefficients(cell, degree(c), [&c](const V3& x) { return phi(c, x); }, coeffs);
                        append_cell_rules<double>(cell.type, cell.vertices, degree(c), coeffs, term, 0, q, opt, 0,
                                                  rule, stats);
                    }
                }
                catch (const std::exception&)
                {
                    ++r.fail;
                    continue;
                }
                r.max_bisections = std::max(r.max_bisections, stats.bisections);
                double value = 0;
                for (std::size_t p = 0; p < rule._weights.size(); ++p)
                {
                    const std::array<double, 2> xi = {rule._points[2 * p], rule._points[2 * p + 1]};
                    const double w = rule._weights[p];
                    if (!std::isfinite(w) || !std::isfinite(xi[0]) || !std::isfinite(xi[1]))
                    {
                        ++r.fail;
                        continue;
                    }
                    value += w;
                    r.negative += w < 0;
                    r.outside += !in_reference_cell(cell.type, xi, 1e-12);
                    const V3 x = physical(cell, xi);
                    const double f = phi(c, x), reach = 1e-9 * h * gradient_norm(c, x);
                    r.side += surface ? !(std::abs(f) <= reach) : !(f <= reach);
                }
                r.total += value;
                if (!(surface && ref.edge_length > 0))
                    r.l1 += std::abs(value - (surface ? ref.length : ref.area));
            }
    return r;
}

SelectionTerm term_of(const char* text)
{
    SelectionExpr expr = parse_selection_expr(text);
    compile_selection_expr(expr, {"phi"});
    return expr.terms.front();
}
} // namespace

int main()
{
    int failures = 0;
    const Options opt;

    // lines: exact on every cell
    {
        double worst_overall = 0;
        const std::vector<std::pair<V3, double>> lines = {{{1.0, 0.3, 0}, 0.0}, {{-0.45, 1.0, 0}, 0.1}};
        for (const auto& [normal, offset] : lines)
            for (const std::string mesh : {"quad", "tri"})
                for (const int source : {1, 2, 0}) // Bernstein of degree 1 and 2, analytic
                    for (const int q : {3, 5})
                        for (const char* part : {"phi < 0", "phi > 0", "phi = 0"})
                        {
                            Case c;
                            c.shape = Shape::line;
                            c.normal = normal;
                            c.offset = offset;
                            const SelectionTerm term = term_of(part);
                            const Part p = part_of(term);
                            const CaseLevelSet functor = {&c};
                            const AnalyticLevelSet phi_analytic = analytic_level_set(functor);
                            const int n = 7;
                            const double h = 2.0 / n;
                            double total = 0, exact_total = 0, worst = 0;
                            std::vector<double> coeffs;
                            for (int i0 = 0; i0 < n; ++i0)
                                for (int i1 = 0; i1 < n; ++i1)
                                    for (const TestCell2D& cell : grid_cells_2d(mesh, {-1 + h * i0, -1 + h * i1, 0}, h))
                                    {
                                        const Reference ref = half_plane(cell, normal, offset);
                                        const double exact_value = p == Part::negative   ? ref.area
                                                                   : p == Part::positive ? cell.area - ref.area
                                                                                         : ref.length;
                                        quadrature::QuadratureRules<double> rule;
                                        Stats stats;
                                        if (source == 0)
                                            append_cell_rules<double>(cell.type, cell.vertices, phi_analytic, term, 0, q,
                                                                      opt, 0, rule, stats);
                                        else
                                        {
                                            cell_coefficients(cell, source, [&c](const V3& x) { return phi(c, x); },
                                                              coeffs);
                                            append_cell_rules<double>(cell.type, cell.vertices, source, coeffs, term, 0,
                                                                      q, opt, 0, rule, stats);
                                        }
                                        double value = 0;
                                        for (double w : rule._weights)
                                            value += w;
                                        total += value;
                                        exact_total += exact_value;
                                        worst = std::max(worst, std::abs(value - exact_value)
                                                                    / (p == Part::interface ? h : h * h));
                                    }
                            const double rel_total = std::abs(total - exact_total) / exact_total;
                            worst_overall = std::max({worst_overall, worst, rel_total});
                            if (worst > 1e-13 || rel_total > 1e-13)
                            {
                                std::printf("FAILED: line (%g, %g) . x = %g, %s, %s, q = %d, %s: worst cell %.1e, "
                                            "total %.1e\n",
                                            normal[0], normal[1], offset, source == 0 ? "analytic" : source == 1 ? "P1" : "P2",
                                            mesh.c_str(), q, part, worst, rel_total);
                                ++failures;
                            }
                        }
        std::printf("lines: largest error %.1e (per cell relative to h^2 or h, totals relative)\n", worst_overall);
    }

    // the circle: per-cell L1 against the exact disk and circle, by q
    {
        // about three times the largest L1 of the two meshes and sources
        struct Bound
        {
            int q;
            const char* part;
            double l1;
        };
        const Bound bounds[] = {{3, "phi < 0", 3e-6},  {3, "phi = 0", 3e-5},  {5, "phi < 0", 3e-9},
                                {5, "phi = 0", 1e-7},  {8, "phi < 0", 3e-13}, {8, "phi = 0", 2e-11}};
        Case c;
        c.centre = {0.0123, -0.0371, 0};
        c.radius = 0.7;
        for (const bool analytic : {false, true})
            for (const std::string mesh : {"quad", "tri"})
                for (const Bound& b : bounds)
                {
                    const SelectionTerm term = term_of(b.part);
                    const Run r = run_case(c, mesh, 8, term, b.q, opt, analytic);
                    const double l1 = r.l1 / r.exact_total;
                    const bool ok = l1 <= b.l1 && r.fail + r.negative + r.outside + r.side == 0;
                    std::printf("circle %-9s %-4s q = %d %-8s L1 %.1e (bound %.1e), max bisections per cell %d%s\n",
                                analytic ? "analytic" : "P2", mesh.c_str(), b.q, b.part, l1, b.l1, r.max_bisections,
                                ok ? "" : " FAILED");
                    failures += !ok;
                }
    }

    // robustness cases
    std::map<std::string, double> baseline;
    for (const bool analytic : {false, true})
        for (const Case& c : all_cases())
            for (const std::string mesh : {"tri", "quad"})
                for (const char* part : {"phi < 0", "phi = 0"})
                {
                    const SelectionTerm term = term_of(part);
                    const Run r = run_case(c, mesh, 16, term, 3, opt, analytic);
                    std::vector<std::string> problems;
                    if (r.fail + r.negative + r.outside + r.side > 0)
                        problems.push_back("fail/neg/out/side " + std::to_string(r.fail) + "/" + std::to_string(r.negative)
                                           + "/" + std::to_string(r.outside) + "/" + std::to_string(r.side));
                    if (r.max_bisections > opt.max_bisections)
                        problems.push_back(std::to_string(r.max_bisections) + " bisections in one cell");
                    const double l1 = r.l1 / r.exact_total;
                    if (l1 > 1e-4)
                        problems.push_back("per-cell L1 " + std::to_string(l1));
                    const std::string key = (analytic ? "analytic " : "") + mesh + part;
                    if (c.name == "circle")
                        baseline[key] = r.total;
                    if (c.scale != 1.0 && std::abs(r.total - baseline[key]) > 1e-12 * baseline[key])
                        problems.push_back("total differs from the unscaled circle");
                    std::string message;
                    for (const std::string& p : problems)
                        message += " " + p + ";";
                    std::printf("%-9s %-20s %-4s %-8s L1 %.1e, max bisections per cell %d%s%s\n",
                                analytic ? "analytic" : "Bernstein", c.name.c_str(), mesh.c_str(), part, l1,
                                r.max_bisections, problems.empty() ? "" : " FAILED:", message.c_str());
                    failures += !problems.empty();
                }
    return failures == 0 ? 0 : 1;
}
