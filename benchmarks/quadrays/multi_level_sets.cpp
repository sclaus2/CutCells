// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// quadrays against algoim's multi-polynomial engine (Saye 2022) on parts bounded
// by several polynomial level sets: both take the same polynomials on the same
// cubes of a grid. For each case, part, mesh and points per segment q, one CSV
// line per generator: points on cut cells (total, mean per cut cell with points
// of the part, worst cell), time, per-cell L1 and worst errors against quadrays
// at q = 14 without limits, of the part's measure, of a polynomial of degree 8
// in physical coordinates, and of one in the cell's coordinates, which varies
// over the cell at any h as a product of shape functions does and so is what
// assembly sees; the error of the total, and the smallest weight.
//
// usage: quadrays_multi_level_sets [case|all] [n,n,...] [q,q,...] [generator,...]

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <cutcells/quadrays/analytic.h>
#include <cutcells/quadrays/rules.h>
#include <cutcells/selection_expr.h>

#include "quadrature_multipoly.hpp"
#include "test_mesh.h"

using namespace cutcells;
using namespace cutcells::quadrays;
using namespace cutcells::quadrays::support;

namespace
{

constexpr double radius = 0.5;
const V3 c1 = {-0.21, 0.07, 0.04}, c2 = {0.24, 0.03, -0.02};

struct LevelSet
{
    std::function<double(const V3&)> value;
    int degree = 2;
};

struct Clause
{
    int level_set = 0;
    int sign = 0; ///< -1: < 0, 1: > 0, 0: = 0
};

struct Selection
{
    std::string name;
    std::string text;
    std::vector<Clause> clauses;
};

struct Case
{
    std::string name;
    std::vector<LevelSet> level_sets;
    std::vector<Selection> parts;
};

LevelSet ball(const V3& c, double r)
{
    return {[c, r](const V3& x)
            { return (x[0] - c[0]) * (x[0] - c[0]) + (x[1] - c[1]) * (x[1] - c[1]) + (x[2] - c[2]) * (x[2] - c[2]) - r * r; },
            2};
}

LevelSet ball_distance(const V3& c, double r, int degree)
{
    return {[c, r](const V3& x)
            { return std::sqrt((x[0] - c[0]) * (x[0] - c[0]) + (x[1] - c[1]) * (x[1] - c[1]) + (x[2] - c[2]) * (x[2] - c[2])) - r; },
            degree};
}

LevelSet cylinder(const V3& c, double r, int axis)
{
    return {[c, r, axis](const V3& x)
            {
                const int i = (axis + 1) % 3, j = (axis + 2) % 3;
                return (x[i] - c[i]) * (x[i] - c[i]) + (x[j] - c[j]) * (x[j] - c[j]) - r * r;
            },
            2};
}

LevelSet plane(const V3& n, double d)
{
    return {[n, d](const V3& x) { return n[0] * x[0] + n[1] * x[1] + n[2] * x[2] - d; }, 1};
}

LevelSet torus(const V3& c, double major, double minor)
{
    return {[c, major, minor](const V3& x)
            {
                const double a = x[0] - c[0], b = x[1] - c[1], e = x[2] - c[2];
                const double s = a * a + b * b + e * e + major * major - minor * minor;
                return s * s - 4 * major * major * (a * a + b * b);
            },
            4};
}

std::vector<Case> all_cases()
{
    const Selection vol1 = {"volume", "a < 0", {{0, -1}}}, surf1 = {"interface", "a = 0", {{0, 0}}};
    const Selection vol2 = {"volume", "a < 0 and b < 0", {{0, -1}, {1, -1}}};
    const Selection surf2 = {"interface", "a = 0 and b < 0", {{0, 0}, {1, -1}}};
    const Selection vol3 = {"volume", "a < 0 and b < 0 and c < 0", {{0, -1}, {1, -1}, {2, -1}}};
    const Selection surf3 = {"interface", "a = 0 and b < 0 and c < 0", {{0, 0}, {1, -1}, {2, -1}}};
    const V3 c3 = {0.0123, -0.0371, 0.0217};
    return {
        {"ball", {ball(c1, radius)}, {vol1, surf1}},
        {"ball_in_cell", {ball({0.25, 0.25, 0.25}, 0.26)}, {vol1, surf1}},
        {"lens", {ball(c1, radius), ball(c2, radius)}, {vol2, surf2}},
        {"lens_p3", {ball_distance(c1, radius, 3), ball_distance(c2, radius, 3)}, {vol2, surf2}},
        {"ball_plane", {ball(c1, radius), plane({0.93, 0.33, -0.19}, 0.05)}, {vol2, surf2}},
        {"bicylinder", {cylinder(c3, 0.6, 2), cylinder(c3, 0.6, 1)}, {vol2, surf2}},
        {"tricylinder", {cylinder(c3, 0.6, 2), cylinder(c3, 0.6, 1), cylinder(c3, 0.6, 0)}, {vol3, surf3}},
        {"torus_ball", {torus(c3, 0.5, 0.2), ball(c1, radius)}, {vol2, surf2}},
    };
}

/// The cells of a grid with the forms of the level sets on their boxes.
struct GridCell
{
    TestCell cell;
    ClippedBox<double> box;
    std::vector<BoxBernstein<double>> forms;
    std::vector<std::vector<algoim::real>> poly_data;
    std::vector<algoim::uvector<int, 3>> poly_extents;
    double side = 0;
    bool cut = false;
};

std::vector<GridCell> make_cells(const Case& c, int n)
{
    const double h = 2.0 / n;
    std::vector<GridCell> cells;
    cells.reserve(static_cast<std::size_t>(n * n * n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
            {
                GridCell g;
                g.cell = grid_cells("hex", {-1 + i * h, -1 + j * h, -1 + k * h}, h, {0, 0, 0}).front();
                g.side = h;
                make_clipped_box<double>(g.cell.type, g.cell.vertices, 3, g.box);
                g.forms.resize(c.level_sets.size());
                std::vector<double> coeffs;
                for (std::size_t l = 0; l < c.level_sets.size(); ++l)
                {
                    const LevelSet& ls = c.level_sets[l];
                    cell_coefficients(g.cell, ls.degree, ls.value, coeffs);
                    cell_bernstein_on_box<double>(g.cell.type, ls.degree, coeffs, g.box, g.forms[l]);
                    const auto [lo, hi] = std::minmax_element(g.forms[l].coeffs.begin(), g.forms[l].coeffs.end());
                    g.cut = g.cut || (*lo < 0 && *hi > 0);
                }
                cells.push_back(std::move(g));
            }
    return cells;
}

/// The same polynomials for algoim, made only for the cells that are cut:
/// interpolated from the form, which is exact at its degree.
void make_polys(GridCell& g)
{
    g.poly_data.resize(g.forms.size());
    g.poly_extents.clear();
    for (std::size_t l = 0; l < g.forms.size(); ++l)
    {
        const Source<double> source = bernstein_source(g.forms[l]);
        const algoim::uvector<int, 3> ext(g.forms[l].degree[0] + 1, g.forms[l].degree[1] + 1, g.forms[l].degree[2] + 1);
        g.poly_data[l].assign(static_cast<std::size_t>(ext(0) * ext(1) * ext(2)), 0.0);
        algoim::xarray<algoim::real, 3> poly(g.poly_data[l].data(), ext);
        algoim::bernstein::bernsteinInterpolate<3>(
            [&source](const algoim::uvector<algoim::real, 3>& x)
            {
                const double u[3] = {x(0), x(1), x(2)};
                return evaluate(source, std::span<const double>(u, 3));
            },
            poly);
        g.poly_extents.push_back(ext);
    }
}

/// Algoim's view of polynomial @p l of a cell.
algoim::xarray<algoim::real, 3> poly_of(const GridCell& g, std::size_t l)
{
    return algoim::xarray<algoim::real, 3>(const_cast<algoim::real*>(g.poly_data[l].data()), g.poly_extents[l]);
}

/// The polynomial integrand of the comparison, in physical coordinates.
double integrand(const V3& x)
{
    return std::pow(1.0 + 0.7 * x[0] + 0.5 * x[1] - 0.3 * x[2], 8);
}

/// A polynomial in the cell's coordinates with degree 8 along each axis, the
/// squares of Legendre polynomials of degree 4: as the mass matrix of shape
/// functions of degree 4 along that axis.
double shape_integrand(double u0, double u1, double u2)
{
    const auto square = [](double u)
    {
        const double t = 2 * u - 1, t2 = t * t;
        const double p = (35 * t2 * t2 - 30 * t2 + 3) / 8;
        return p * p;
    };
    return square(u0) + square(u1) + square(u2);
}

struct CellResult
{
    int points = 0;
    double sum = 0;  ///< the measure
    double moment = 0; ///< the integral of integrand
    double shape = 0; ///< the integral of shape_integrand
    double min_weight = 0;
};

CellResult run_quadrays(const GridCell& g, const std::vector<std::vector<SelectionTerm>>& groups, int q, const Options& opt)
{
    std::vector<Source<double>> sources;
    for (const BoxBernstein<double>& f : g.forms)
        sources.push_back(bernstein_source(f));
    quadrature::QuadratureRules<double> rule;
    Stats stats;
    for (const std::vector<SelectionTerm>& terms : groups)
        append_rules<double>(g.box, std::span<const Source<double>>(sources), std::span<const SelectionTerm>(terms), q, opt,
                             0, rule, stats);
    CellResult r;
    r.points = static_cast<int>(rule._weights.size());
    r.min_weight = 1e300;
    for (std::size_t p = 0; p < rule._weights.size(); ++p)
    {
        const double w = rule._weights[p];
        r.sum += w;
        r.moment += w * integrand(physical(g.cell, {rule._points[3 * p], rule._points[3 * p + 1], rule._points[3 * p + 2]}));
        r.shape += w * shape_integrand(rule._points[3 * p], rule._points[3 * p + 1], rule._points[3 * p + 2]);
        r.min_weight = std::min(r.min_weight, w);
    }
    return r;
}

CellResult run_algoim(const GridCell& g, const Selection& part, int q, algoim::QuadStrategy strategy)
{
    algoim::ImplicitPolyQuadrature<3> ipq;
    for (std::size_t l = 0; l < g.poly_data.size(); ++l)
    {
        const algoim::xarray<algoim::real, 3> poly = poly_of(g, l);
        const auto mask = algoim::detail::nonzeroMask(poly, algoim::booluarray<3, ALGOIM_M>(true));
        if (!algoim::detail::maskEmpty(mask))
            ipq.phi.push_back(poly, mask);
    }
    ipq.build(true, false);
    const auto value = [&g](std::size_t l, const algoim::uvector<algoim::real, 3>& x)
    { return algoim::bernstein::evalBernsteinPoly(poly_of(g, l), x); };
    CellResult r;
    r.min_weight = 1e300;
    const auto add = [&r, &g](const algoim::uvector<algoim::real, 3>& x, double w)
    {
        ++r.points;
        r.sum += w;
        r.moment += w * integrand(physical(g.cell, {x(0), x(1), x(2)}));
        r.shape += w * shape_integrand(x(0), x(1), x(2));
        r.min_weight = std::min(r.min_weight, w);
    };
    const bool interface = part.clauses.front().sign == 0;
    if (interface)
    {
        const std::size_t on = static_cast<std::size_t>(part.clauses.front().level_set);
        ipq.integrate_surf(strategy, q,
                           [&](const algoim::uvector<algoim::real, 3>& x, algoim::real w,
                               const algoim::uvector<algoim::real, 3>&)
                           {
                               // the polynomial this point lies on
                               double smallest = 1e300;
                               std::size_t lying = 0;
                               for (std::size_t l = 0; l < g.poly_data.size(); ++l)
                                   if (const double v = std::abs(value(l, x)); v < smallest)
                                   {
                                       smallest = v;
                                       lying = l;
                                   }
                               if (lying != on)
                                   return;
                               for (std::size_t c = 1; c < part.clauses.size(); ++c)
                                   if (value(static_cast<std::size_t>(part.clauses[c].level_set), x) * part.clauses[c].sign <= 0)
                                       return;
                               add(x, w * g.side * g.side);
                           });
    }
    else
    {
        ipq.integrate(strategy, q,
                      [&](const algoim::uvector<algoim::real, 3>& x, algoim::real w)
                      {
                          for (const Clause& c : part.clauses)
                              if (value(static_cast<std::size_t>(c.level_set), x) * c.sign <= 0)
                                  return;
                          add(x, w * g.side * g.side * g.side);
                      });
    }
    return r;
}

std::vector<std::vector<SelectionTerm>> groups_of(const Selection& part, std::size_t n)
{
    std::vector<std::string> names;
    for (std::size_t l = 0; l < n; ++l)
        names.push_back(std::string(1, static_cast<char>('a' + l)));
    SelectionExpr expr = parse_selection_expr(part.text);
    compile_selection_expr(expr, names, {});
    return {expr.terms};
}

Options quadrays_options(const std::string& generator)
{
    Options opt;
    opt.max_points = -1;
    if (generator.find("rot12") != std::string::npos)
        opt.rotation_depth = 12;
    if (generator.find("norot") != std::string::npos)
    {
        opt.rotation_depth = 100;
        opt.diagonal_frames = false;
    }
    if (generator.find("nocheck") != std::string::npos)
        opt.tolerance = -1;
    if (generator.find("nofold") != std::string::npos)
        opt.fold_splits = false;
    if (const std::size_t m = generator.find("-m"); m != std::string::npos)
        opt.margin = std::stoi(generator.substr(m + 2)) / 100.0;
    return opt;
}

std::vector<std::string> split(const std::string& text)
{
    std::vector<std::string> out;
    std::stringstream stream(text);
    for (std::string item; std::getline(stream, item, ',');)
        out.push_back(item);
    return out;
}

} // namespace

int main(int argc, char** argv)
{
    const std::string which = argc > 1 ? argv[1] : "all";
    std::vector<int> ns, qs;
    for (const std::string& s : split(argc > 2 ? argv[2] : "3,5,9,17"))
        ns.push_back(std::stoi(s));
    for (const std::string& s : split(argc > 3 ? argv[3] : "2,3,4,6,8"))
        qs.push_back(std::stoi(s));
    const std::vector<std::string> generators
        = split(argc > 4 ? argv[4] : "algoim-gl,algoim-auto,quadrays,quadrays-nocheck,quadrays-norot-nocheck,quadrays-nocheck-m10");

    std::printf("case,part,n,q,generator,cut_cells,cells_with_points,points,mean,worst,seconds,cell_l1,cell_max,total_err,"
                "moment_l1,moment_max,shape_l1,shape_max,min_weight\n");
    for (const Case& c : all_cases())
    {
        if (which != "all" && which != c.name)
            continue;
        for (const int n : ns)
        {
            std::vector<GridCell> cells = make_cells(c, n);
            for (GridCell& g : cells)
                if (g.cut)
                    make_polys(g);
            for (const Selection& part : c.parts)
            {
                const std::vector<std::vector<SelectionTerm>> groups = groups_of(part, c.level_sets.size());
                // the reference: quadrays with many points and no limits
                std::vector<double> reference(cells.size(), 0.0), reference_moment(cells.size(), 0.0),
                    reference_shape(cells.size(), 0.0);
                double reference_total = 0;
                Options ref_opt;
                ref_opt.max_points = -1;
                int cut_cells = 0;
                for (std::size_t i = 0; i < cells.size(); ++i)
                    if (cells[i].cut)
                    {
                        const CellResult r = run_quadrays(cells[i], groups, 14, ref_opt);
                        reference[i] = r.sum;
                        reference_moment[i] = r.moment;
                        reference_shape[i] = r.shape;
                        reference_total += reference[i];
                        ++cut_cells;
                    }
                for (const int q : qs)
                    for (const std::string& generator : generators)
                    {
                        const bool algoim_engine = generator.rfind("algoim", 0) == 0;
                        const algoim::QuadStrategy strategy
                            = generator == "algoim-gl" ? algoim::AlwaysGL : algoim::AutoMixed;
                        const Options opt = quadrays_options(generator);
                        long long points = 0;
                        int worst = 0, with_points = 0;
                        double l1 = 0, largest = 0, total = 0, min_weight = 1e300, seconds = 0, moment_l1 = 0,
                               moment_max = 0, shape_l1 = 0, shape_max = 0;
                        for (std::size_t i = 0; i < cells.size(); ++i)
                        {
                            if (!cells[i].cut)
                                continue;
                            const auto t0 = std::chrono::steady_clock::now();
                            const CellResult r = algoim_engine ? run_algoim(cells[i], part, q, strategy)
                                                               : run_quadrays(cells[i], groups, q, opt);
                            seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                            points += r.points;
                            worst = std::max(worst, r.points);
                            with_points += r.points > 0;
                            total += r.sum;
                            min_weight = std::min(min_weight, r.min_weight);
                            const double e = std::abs(r.sum - reference[i]);
                            l1 += e;
                            largest = std::max(largest, e);
                            const double em = std::abs(r.moment - reference_moment[i]);
                            moment_l1 += em;
                            moment_max = std::max(moment_max, em);
                            const double es = std::abs(r.shape - reference_shape[i]);
                            shape_l1 += es;
                            shape_max = std::max(shape_max, es);
                        }
                        std::printf("%s,%s,%d,%d,%s,%d,%d,%lld,%.1f,%d,%.4f,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n", c.name.c_str(),
                                    part.name.c_str(), n, q, generator.c_str(), cut_cells, with_points, points,
                                    with_points > 0 ? double(points) / with_points : 0.0, worst, seconds, l1, largest,
                                    std::abs(total - reference_total), moment_l1, moment_max, shape_l1, shape_max, min_weight);
                        std::fflush(stdout);
                    }
            }
        }
    }
    return 0;
}
