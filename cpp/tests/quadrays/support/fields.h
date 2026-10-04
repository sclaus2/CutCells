// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Several level sets on the cells of a test mesh: balls, planes and cylinders
// as Bernstein forms or analytic level sets, the rules of a selection on every
// cell with their checks (finite, positive, inside the cell, in the part), and
// the leaves of a selection. Shared by test_level_sets and test_curves.

#pragma once

#include <cutcells/quadrays/analytic.h>
#include <cutcells/quadrays/leaves.h>
#include <cutcells/quadrays/rules.h>
#include <cutcells/quadrays/taylor.h>
#include <cutcells/selection_expr.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

#include "test_mesh.h"

namespace cutcells::quadrays::support
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

inline Field ball(const V3& c, double r)
{
    Field f;
    f.c = c;
    f.r = r;
    return f;
}

inline Field plane(const V3& n, double d)
{
    Field f;
    f.kind = Field::Kind::plane;
    f.n = n;
    f.d = d;
    return f;
}

inline Field cylinder(const V3& c, double r, int axis)
{
    Field f;
    f.kind = Field::Kind::cylinder;
    f.c = c;
    f.r = r;
    f.axis = axis;
    return f;
}

inline double value(const Field& f, const V3& x) { return f(std::array<double, 3>{x[0], x[1], x[2]}); }

inline double gradient_norm(const Field& f, const V3& x)
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

/// Does x lie in the part of one of the terms, up to @p tol h of each zero set?
inline bool in_part(const std::vector<Field>& fields, const std::vector<SelectionTerm>& terms, const V3& x, double h,
                    double tol = 1e-9)
{
    for (const SelectionTerm& t : terms)
    {
        bool ok = true;
        for (std::size_t l = 0; l < fields.size() && ok; ++l)
        {
            const std::uint64_t bit = std::uint64_t(1) << l;
            const double f = value(fields[l], x), reach = tol * h * gradient_norm(fields[l], x);
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
    int curve_lost = 0; ///< Stats::curve_lost, summed over the cells
    double total = 0, exact_total = 0, l1 = 0;
    std::vector<V3> points;      ///< every point, physical
    std::vector<double> weights; ///< its weight

    long problems() const { return fail + negative + outside + side; }
};

/// Adds the counts, totals and points of @p r to @p sum.
inline void accumulate(Run& sum, const Run& r)
{
    sum.fail += r.fail;
    sum.negative += r.negative;
    sum.outside += r.outside;
    sum.side += r.side;
    sum.max_bisections = std::max(sum.max_bisections, r.max_bisections);
    sum.curve_lost += r.curve_lost;
    sum.total += r.total;
    sum.exact_total += r.exact_total;
    sum.l1 += r.l1;
    sum.points.insert(sum.points.end(), r.points.begin(), r.points.end());
    sum.weights.insert(sum.weights.end(), r.weights.begin(), r.weights.end());
}

/// The level sets on one cell as sources: analytic, or their P2 (P1)
/// interpolants in Bernstein form, kept in @p forms.
template <typename Cell>
std::vector<Source<double>> cell_sources(const std::vector<Field>& fields, const std::vector<AnalyticLevelSet>& als,
                                         const Cell& cell, const ClippedBox<double>& box, bool analytic,
                                         std::vector<std::vector<double>>& coeffs,
                                         std::vector<BoxBernstein<double>>& forms)
{
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
        cell_bernstein_on_box<double>(cell.type, f.degree(), coeffs[l], box, forms[l]);
        sources.push_back(bernstein_source(forms[l]));
    }
    return sources;
}

/// The rules of the part @p text of the level sets @p fields (named a, b, c) on
/// every cell, their checks, and the per-cell errors if @p exact (per cell) is
/// given.
template <typename Cell>
Run run(const std::vector<Field>& fields, const char* text, const std::vector<Cell>& cells, double h, int q,
        bool analytic, const std::vector<double>& exact, const Options& opt = Options{})
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
    Run r;
    for (std::size_t i = 0; i < cells.size(); ++i)
    {
        const Cell& cell = cells[i];
        ClippedBox<double> box;
        make_clipped_box<double>(cell.type, cell.vertices, tdim, box);
        const std::vector<Source<double>> sources = cell_sources(fields, als, cell, box, analytic, coeffs, forms);
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
        r.curve_lost += stats.curve_lost;
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
            r.points.push_back(x);
            r.weights.push_back(w);
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

/// Leaves of the part @p text, of degree 2, on every cell: how many, how many
/// of the part were dropped because their nodes did not line up, how many
/// nodes lie outside the part, beyond 1e-5 h of a zero set (nodes are pulled
/// 1e-5 of their segment inside its ends), and how many of each VTK type.
struct LeafRun
{
    long leaves = 0, incomplete = 0, outside = 0;
    std::map<int, long> types;
};

template <typename Cell>
LeafRun leaves_of(const std::vector<Field>& fields, const char* text, const std::vector<Cell>& cells, double h,
                  bool analytic)
{
    constexpr int tdim = std::is_same_v<Cell, TestCell2D> ? 2 : 3;
    const std::vector<std::string> names = {"a", "b", "c"};
    SelectionExpr expr = parse_selection_expr(text);
    compile_selection_expr(expr, std::vector<std::string>(names.begin(), names.begin() + fields.size()));
    std::vector<AnalyticLevelSet> als;
    for (const Field& f : fields)
        als.push_back(analytic_level_set(f));
    std::vector<std::vector<double>> coeffs(fields.size());
    std::vector<BoxBernstein<double>> forms(fields.size());
    LeafRun r;
    for (std::size_t i = 0; i < cells.size(); ++i)
    {
        ClippedBox<double> box;
        make_clipped_box<double>(cells[i].type, cells[i].vertices, tdim, box);
        const std::vector<Source<double>> sources = cell_sources(fields, als, cells[i], box, analytic, coeffs, forms);
        LeafMesh<double> leaves;
        Stats stats;
        append_leaves<double>(box, std::span<const Source<double>>(sources),
                              std::span<const SelectionTerm>(expr.terms), 2, Options{}, static_cast<std::int32_t>(i),
                              leaves, stats);
        r.leaves += leaves.n_cells();
        r.incomplete += stats.incomplete_leaves;
        for (const std::uint8_t type : leaves.vtk_types)
            ++r.types[type];
        for (int p = 0; p < leaves.n_points(); ++p)
        {
            const double* x = leaves.points.data() + 3 * p;
            r.outside += !in_part(fields, expr.terms, V3{x[0], x[1], x[2]}, h, 1e-5);
        }
    }
    return r;
}

inline std::vector<TestCell> mesh_3d(const std::string& mesh, int n, const V3& origin)
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

inline std::vector<TestCell2D> mesh_2d(const std::string& mesh, int n)
{
    const double h = 2.0 / n;
    std::vector<TestCell2D> cells;
    for (int i0 = 0; i0 < n; ++i0)
        for (int i1 = 0; i1 < n; ++i1)
            for (const TestCell2D& cell : grid_cells_2d(mesh, {-1 + h * i0, -1 + h * i1, 0}, h))
                cells.push_back(cell);
    return cells;
}

} // namespace cutcells::quadrays::support
