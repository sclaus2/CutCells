// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// quadrays on facets: cells in a space of higher dimension than their own, as
// the facets of a mesh are, and intervals.
//  - A triangle, a parallelogram and an interval moved into 3D by a rotation
//    integrate as they do in their own dimension: a plane cutting them gives
//    the same measures up to rounding, for Bernstein and analytic level sets,
//    also when the analytic level set varies normal to the cell.
//  - The surface of [-1, 1]^3, rotated, as 16 x 16 quadrilaterals or twice as
//    many triangles per face in 3D: a sphere about a point near the centre
//    cuts a disk from each face, whose areas and circles part::cut and
//    quadrays give to 1e-9; rule points lie in the facets' reference cells
//    and on the sphere; the faces add up to 24. A plane through mesh edges of
//    that surface lies on 16 of them: rules on the edges, measured in 3D, with
//    the two triangles sharing each, on either side.
//  - The edges of the square [-1, 1]^2 as intervals in 2D and rotated into
//    3D, and the interval [-1, 1] in 1D: a circle's chords, and its crossing
//    points, each once with weight 1; a level set vanishing at mesh vertices
//    gives a rule on each, of weight 1, with the two intervals sharing it.
// Exits non-zero on failure.

#include <cutcells/bernstein.h>
#include <cutcells/level_set.h>
#include <cutcells/part/cut_result.h>
#include <cutcells/part/mesh_part.h>
#include <cutcells/part/output.h>
#include <cutcells/quadrays/analytic.h>
#include <cutcells/quadrays/rules.h>
#include <cutcells/selection_expr.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

#include "../quadrays/support/exact_reference.h"
#include "../quadrays/support/test_mesh.h"

using namespace cutcells;
using namespace cutcells::quadrays::support;

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

using M3 = std::array<V3, 3>;

/// A rotation about the axis (1, 2, 3) by 0.7, rows.
M3 rotation()
{
    const V3 a = exact::unit({1, 2, 3});
    const double c = std::cos(0.7), s = std::sin(0.7);
    M3 q;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            q[i][j] = (i == j ? c : 0) + (1 - c) * a[i] * a[j];
    q[0][1] -= s * a[2];
    q[0][2] += s * a[1];
    q[1][0] += s * a[2];
    q[1][2] -= s * a[0];
    q[2][0] -= s * a[1];
    q[2][1] += s * a[0];
    return q;
}

const M3 Q = rotation();
const V3 shift = {0.3, -0.2, 0.5};

/// Q x + shift.
V3 move(const V3& x)
{
    V3 y = shift;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            y[i] += Q[i][j] * x[j];
    return y;
}

/// n . x - d + kappa x_2 in the unmoved coordinates: on the plane z = 0 the
/// line n . x = d; kappa its slope normal to the plane.
struct TiltedPlane
{
    V3 n = {0, 0, 0};
    double d = 0, kappa = 0;
    bool moved = false;

    template <typename V>
    V operator()(const std::array<V, 3>& y) const
    {
        std::array<V, 3> x = y;
        if (moved)
        {
            for (int i = 0; i < 3; ++i)
            {
                x[i] = V(0.0);
                for (int j = 0; j < 3; ++j)
                    x[i] = x[i] + Q[j][i] * (y[j] - shift[j]);
            }
        }
        return n[0] * x[0] + n[1] * x[1] + kappa * x[2] - d;
    }
};

/// |y - c|^2 - r^2.
struct Ball
{
    V3 c = {0, 0, 0};
    double r = 0;

    template <typename V>
    V operator()(const std::array<V, 3>& y) const
    {
        const V a = y[0] - c[0], b = y[1] - c[1], e = y[2] - c[2];
        return a * a + b * b + e * e - r * r;
    }
};

SelectionTerm term_of(const char* text)
{
    SelectionExpr expr = parse_selection_expr(text);
    compile_selection_expr(expr, {"phi"});
    return expr.terms.front();
}

double sum(const std::vector<double>& w)
{
    double s = 0;
    for (const double x : w)
        s += x;
    return s;
}

// ============================================================================
// One cell moved into 3D
// ============================================================================

/// A triangle, a parallelogram and an interval cut by a line, in their own
/// dimension and moved into 3D, against each other.
void test_moved_cells()
{
    const quadrays::Options opt;
    const int q = 4;
    const V3 n = exact::unit({1, 0.4, 0});
    const double d = 0.42;
    double worst = 0;
    for (const char* kind : {"tri", "quad", "interval"})
        for (const int source : {1, 0}) // Bernstein of degree 1, analytic
            for (const char* text : {"phi < 0", "phi > 0", "phi = 0"})
            {
                // the cell in the plane z = 0 (on the x axis for the interval)
                cell::type type = cell::type::triangle;
                std::vector<V3> v = {{0.1, 0.2, 0}, {0.9, 0.35, 0}, {0.3, 0.8, 0}};
                if (std::string(kind) == "quad")
                {
                    type = cell::type::quadrilateral;
                    v = {{0.1, 0.1, 0}, {0.8, 0.2, 0}, {0.3, 0.7, 0}, {1.0, 0.8, 0}};
                }
                else if (std::string(kind) == "interval")
                {
                    type = cell::type::interval;
                    v = {{0.1, 0.0, 0}, {0.9, 0.0, 0}};
                }
                const int tdim = cell::get_tdim(type);
                std::vector<double> own, moved;
                for (const V3& x : v)
                {
                    own.insert(own.end(), x.begin(), x.begin() + tdim);
                    const V3 y = move(x);
                    moved.insert(moved.end(), y.begin(), y.end());
                }
                const SelectionTerm term = term_of(text);
                quadrature::QuadratureRules<double> a, b;
                quadrays::Stats stats;
                if (source == 0)
                {
                    // the moved level set varies normal to the cell (kappa)
                    const TiltedPlane flat = {n, d, 0.0, false}, tilted = {n, d, 0.8, true};
                    const quadrays::AnalyticLevelSet pa = quadrays::analytic_level_set(flat),
                                                     pb = quadrays::analytic_level_set(tilted);
                    quadrays::append_cell_rules<double>(type, own, pa, term, 0, q, opt, 0, a, stats);
                    quadrays::append_cell_rules<double>(type, moved, pb, term, 0, q, opt, 0, b, stats);
                }
                else
                {
                    // degree 1: the values at the vertices
                    std::vector<double> values, coeffs;
                    for (const V3& x : v)
                        values.push_back(exact::dot(n, x) - d);
                    bernstein::lagrange_to_bernstein<double>(type, 1, cell::reference_vertices<double>(type), values,
                                                             coeffs);
                    quadrays::append_cell_rules<double>(type, own, 1, coeffs, term, 0, q, opt, 0, a, stats);
                    quadrays::append_cell_rules<double>(type, moved, 1, coeffs, term, 0, q, opt, 0, b, stats);
                }
                const double sa = sum(a._weights), sb = sum(b._weights);
                const double error = std::abs(sa - sb) / std::max(1.0, std::abs(sa));
                worst = std::max(worst, error);
                check(error <= 1e-13 && sa > 0,
                      std::string("moved ") + kind + (source == 0 ? " analytic " : " P1 ") + text + ": "
                          + std::to_string(sa) + " in its own dimension, " + std::to_string(sb) + " moved");
                check(b._tdim == tdim, std::string("moved ") + kind + ": rules of the cell's dimension");
            }
    std::printf("moved cells: largest relative difference %.1e\n", worst);
}

// ============================================================================
// The surface of the cube
// ============================================================================

/// A mesh of the boundary of [-1, 1]^3, moved, with its storage.
struct Surface
{
    std::vector<double> coordinates;
    std::vector<int> connectivity, offsets = {0};
    std::vector<cell::type> types;
    MeshView<double, int> view;
};

/// The boundary of [-1, 1]^3 with n x n quadrilaterals per face, or twice as
/// many triangles, every node moved.
void make_surface(const std::string& kind, int n, Surface& s)
{
    std::map<std::array<long, 3>, int> nodes;
    auto node = [&](const V3& x)
    {
        const std::array<long, 3> key = {std::lround(x[0] * n), std::lround(x[1] * n), std::lround(x[2] * n)};
        auto [it, added] = nodes.emplace(key, static_cast<int>(nodes.size()));
        if (added)
        {
            const V3 y = move(x);
            s.coordinates.insert(s.coordinates.end(), y.begin(), y.end());
        }
        return it->second;
    };
    const double h = 2.0 / n;
    for (int a = 0; a < 3; ++a)
        for (const double side : {-1.0, 1.0})
            for (int j = 0; j < n; ++j)
                for (int i = 0; i < n; ++i)
                {
                    auto corner = [&](int di, int dj)
                    {
                        V3 x;
                        x[a] = side;
                        x[(a + 1) % 3] = -1 + h * (i + di);
                        x[(a + 2) % 3] = -1 + h * (j + dj);
                        return node(x);
                    };
                    const int v0 = corner(0, 0), v1 = corner(1, 0), v2 = corner(0, 1), v3 = corner(1, 1);
                    if (kind == "quad")
                    {
                        s.connectivity.insert(s.connectivity.end(), {v0, v1, v2, v3});
                        s.offsets.push_back(static_cast<int>(s.connectivity.size()));
                        s.types.push_back(cell::type::quadrilateral);
                        continue;
                    }
                    for (const std::array<int, 3> t : {std::array<int, 3>{v0, v1, v3}, std::array<int, 3>{v0, v3, v2}})
                    {
                        s.connectivity.insert(s.connectivity.end(), t.begin(), t.end());
                        s.offsets.push_back(static_cast<int>(s.connectivity.size()));
                        s.types.push_back(cell::type::triangle);
                    }
                }
    s.view.gdim = 3;
    s.view.tdim = 2;
    s.view.coordinates = s.coordinates;
    s.view.connectivity = s.connectivity;
    s.view.offsets = s.offsets;
    s.view.cell_types = s.types;
}

/// A Pk level set interpolating @p f on @p mesh.
LevelSetFunction<double, int> pk_level_set(const MeshView<double, int>& mesh, int degree,
                                           const std::function<double(const V3&)>& f)
{
    LevelSetMeshData<double, int> data = create_level_set_mesh_data<double, int>(mesh, degree);
    std::vector<double> values(static_cast<std::size_t>(data.num_dofs()));
    const int gdim = mesh.gdim;
    for (int k = 0; k < data.num_dofs(); ++k)
    {
        const double* x = data.dof_coordinate(k);
        V3 y = {0, 0, 0};
        for (int i = 0; i < gdim; ++i)
            y[i] = x[i];
        values[static_cast<std::size_t>(k)] = f(y);
    }
    return create_level_set_function<double, int>(std::move(data), std::span<const double>(values), "phi");
}

template <typename F>
LevelSetFunction<double, int> analytic_ls(const F& functor)
{
    auto phi = std::make_shared<const quadrays::AnalyticLevelSet>(quadrays::analytic_level_set(functor));
    return create_level_set_function<double, int>(phi, 3, "phi");
}

/// Physical points of a part's rules on a mesh of affine cells: vertex 0 plus
/// the reference coordinates along the edges to vertices 1 (and 2).
std::vector<V3> physical_points(const MeshView<double, int>& mesh, const quadrature::QuadratureRules<double>& rules)
{
    std::vector<V3> out;
    const int tdim = rules._tdim;
    std::vector<double> x;
    std::vector<int> scratch;
    for (std::size_t r = 0; r + 1 < rules._offset.size(); ++r)
    {
        cell_vertex_coords_basix(mesh, rules._parent_map[r], x, scratch);
        const int g = mesh.gdim;
        for (int p = rules._offset[r]; p < rules._offset[r + 1]; ++p)
        {
            V3 y = {0, 0, 0};
            for (int i = 0; i < g; ++i)
            {
                y[i] = x[static_cast<std::size_t>(i)];
                for (int k = 0; k < tdim; ++k)
                    y[i] += rules._points[static_cast<std::size_t>(p * tdim + k)]
                            * (x[static_cast<std::size_t>((k + 1) * g + i)] - x[static_cast<std::size_t>(i)]);
            }
            out.push_back(y);
        }
    }
    return out;
}

/// True if every rule point lies in its cell's reference cell.
bool in_reference_cells(const MeshView<double, int>& mesh, const quadrature::QuadratureRules<double>& rules)
{
    const int tdim = rules._tdim;
    for (std::size_t r = 0; r + 1 < rules._offset.size(); ++r)
    {
        const bool simplex = mesh.cell_type(rules._parent_map[r]) != cell::type::quadrilateral;
        for (int p = rules._offset[r]; p < rules._offset[r + 1]; ++p)
        {
            double s = 0;
            for (int k = 0; k < tdim; ++k)
            {
                const double xi = rules._points[static_cast<std::size_t>(p * tdim + k)];
                if (xi < -1e-12 || xi > 1 + 1e-12)
                    return false;
                s += xi;
            }
            if (simplex && s > 1 + 1e-12)
                return false;
        }
    }
    return true;
}

/// A sphere about a point near the cube's centre cuts a disk from each face.
void test_cube_surface()
{
    const V3 c = {0.1, -0.05, 0.08};
    const double r = 1.2;
    double disks = 0, circles = 0;
    for (int a = 0; a < 3; ++a)
        for (const double side : {-1.0, 1.0})
        {
            const double rho2 = r * r - (side - c[a]) * (side - c[a]);
            disks += std::numbers::pi * rho2;
            circles += 2 * std::numbers::pi * std::sqrt(rho2);
        }
    const Ball ball = {move(c), r};
    double worst = 0;
    for (const std::string kind : {"quad", "tri"})
        for (const char* source : {"P2", "analytic"})
        {
            Surface s;
            make_surface(kind, 16, s);
            const std::vector<LevelSetFunction<double, int>> level_sets
                = {std::string(source) == "P2" ? pk_level_set(s.view, 2, [&ball](const V3& y) { return ball(y); })
                                               : analytic_ls(ball)};
            const part::CutResult<double, int> result
                = part::cut<double, int>(s.view, std::span<const LevelSetFunction<double, int>>(level_sets));
            const std::string what = "cube surface (" + kind + ", " + source + ")";
            check(!result.cut_cells.empty(), what + ": cut cells");
            const auto inside = part::quadrature_rules(part::select(result, "phi < 0"), 8, true);
            const auto outside = part::quadrature_rules(part::select(result, "phi > 0"), 8, true);
            const auto circle = part::quadrature_rules(part::select(result, "phi = 0"), 8, false);
            const double e_in = std::abs(sum(inside._weights) - disks) / disks,
                         e_out = std::abs(sum(outside._weights) - (24 - disks)) / 24,
                         e_circle = std::abs(sum(circle._weights) - circles) / circles;
            worst = std::max({worst, e_in, e_out, e_circle});
            check(e_in <= 1e-9, what + ": disks " + std::to_string(sum(inside._weights)) + ", exact "
                                    + std::to_string(disks));
            check(e_out <= 1e-9, what + ": the rest of the faces");
            check(e_circle <= 1e-8, what + ": circles " + std::to_string(sum(circle._weights)) + ", exact "
                                        + std::to_string(circles));
            check(inside._tdim == 2 && circle._tdim == 2, what + ": rules of dimension 2");
            check(in_reference_cells(s.view, inside) && in_reference_cells(s.view, circle),
                  what + ": points in the reference cells");
            double off = 0;
            for (const V3& y : physical_points(s.view, circle))
                off = std::max(off, std::abs(ball(std::array<double, 3>{y[0], y[1], y[2]})));
            check(off <= 1e-10, what + ": circle points on the sphere (" + std::to_string(off) + ")");
            check(*std::min_element(circle._weights.begin(), circle._weights.end()) > 0,
                  what + ": positive weights");
        }
    std::printf("cube surface: largest relative error %.1e\n", worst);
}

/// The plane z = 0.5 of the cube, moved, on the surface with 4 x 4 x 2
/// triangles per face: its zero set, the square around the cube's sides, lies
/// on 16 mesh edges of length 0.5, each a rule on the edge with the two
/// triangles sharing it, below and above.
void test_surface_entities()
{
    Surface s;
    make_surface("tri", 4, s);
    const V3 n = {Q[0][2], Q[1][2], Q[2][2]};
    const LevelSetFunction<double, int> ls = pk_level_set(
        s.view, 1, [&n](const V3& y) { return exact::dot(n, exact::sub(y, V3{shift[0], shift[1], shift[2]})) - 0.5; });
    const std::vector<LevelSetFunction<double, int>> level_sets = {ls};
    const part::CutResult<double, int> result
        = part::cut<double, int>(s.view, std::span<const LevelSetFunction<double, int>>(level_sets));
    const part::MeshPart<double, int> zero = part::select(result, "phi = 0");
    const std::vector<part::EntityRules<double>> on = part::entity_rules(zero, 2);
    const part::EntityRules<double>& edges = on[1];
    bool sides = edges.n_rules() == 16;
    for (int r = 0; sides && r < edges.n_rules(); ++r)
    {
        const std::int32_t s0 = edges.cell_offsets[r];
        sides = edges.cell_offsets[r + 1] - s0 == 2 && (edges.negative[s0] ^ edges.negative[s0 + 1]) == 1
                && (edges.positive[s0] ^ edges.positive[s0 + 1]) == 1;
    }
    check(result.cut_cells.empty() && sides,
          "surface entities: 16 edge rules, each with its two triangles on either side ("
              + std::to_string(edges.n_rules()) + " rules)");
    check(std::abs(sum(edges.weights) - 8) <= 1e-13,
          "surface entities: the edges' length " + std::to_string(sum(edges.weights)));
    check(std::abs(sum(part::quadrature_rules(zero, 2, false)._weights) - 8) <= 1e-13,
          "surface entities: the cells' rules share the edges' length");
    check(std::abs(sum(part::quadrature_rules(part::select(result, "phi < 0"), 2, true)._weights) - 16) <= 1e-13,
          "surface entities: the surface below the plane");
    std::printf("surface entities: done\n");
}

// ============================================================================
// Intervals
// ============================================================================

/// Intervals with their storage.
struct Lines
{
    std::vector<double> coordinates;
    std::vector<int> connectivity, offsets = {0};
    std::vector<cell::type> types;
    MeshView<double, int> view;
};

/// The boundary of the square [-1, 1]^2 with n intervals per side, in 2D, or
/// moved into 3D; gdim 1: the interval [-1, 1] with n intervals.
void make_lines(int gdim, int n, Lines& l)
{
    const double h = 2.0 / n;
    auto add = [&](const V3& x)
    {
        const V3 y = gdim == 3 ? move(x) : x;
        l.coordinates.insert(l.coordinates.end(), y.begin(), y.begin() + gdim);
        return static_cast<int>(l.coordinates.size()) / gdim - 1;
    };
    if (gdim == 1)
    {
        for (int i = 0; i <= n; ++i)
            add({-1 + h * i, 0, 0});
        for (int i = 0; i < n; ++i)
            l.connectivity.insert(l.connectivity.end(), {i, i + 1});
    }
    else
    {
        // counterclockwise around the square, 4 n nodes
        for (int side = 0; side < 4; ++side)
            for (int i = 0; i < n; ++i)
            {
                const double t = -1 + h * i;
                const V3 x = side == 0 ? V3{t, -1, 0} : side == 1 ? V3{1, t, 0} : side == 2 ? V3{-t, 1, 0}
                                                                                            : V3{-1, -t, 0};
                add(x);
            }
        for (int i = 0; i < 4 * n; ++i)
            l.connectivity.insert(l.connectivity.end(), {i, (i + 1) % (4 * n)});
    }
    for (std::size_t i = 0; i < l.connectivity.size() / 2; ++i)
    {
        l.offsets.push_back(2 * static_cast<int>(i + 1));
        l.types.push_back(cell::type::interval);
    }
    l.view.gdim = gdim;
    l.view.tdim = 1;
    l.view.coordinates = l.coordinates;
    l.view.connectivity = l.connectivity;
    l.view.offsets = l.offsets;
    l.view.cell_types = l.types;
}

struct LineParts
{
    double inside = 0, outside = 0, points = 0;
    std::vector<V3> where;
};

LineParts line_parts(const MeshView<double, int>& mesh, const LevelSetFunction<double, int>& ls)
{
    const std::vector<LevelSetFunction<double, int>> level_sets = {ls};
    const part::CutResult<double, int> result
        = part::cut<double, int>(mesh, std::span<const LevelSetFunction<double, int>>(level_sets));
    LineParts out;
    out.inside = sum(part::quadrature_rules(part::select(result, "phi < 0"), 4, true)._weights);
    out.outside = sum(part::quadrature_rules(part::select(result, "phi > 0"), 4, true)._weights);
    const auto points = part::quadrature_rules(part::select(result, "phi = 0"), 4, true);
    out.points = sum(points._weights);
    out.where = physical_points(mesh, points);
    return out;
}

/// A circle's chords on the square's sides and on [-1, 1], with its crossing
/// points; level sets vanishing at mesh vertices.
void test_intervals()
{
    const V3 c = {0.13, -0.07, 0};
    const double r = 1.2;
    // the chords on the four sides, all inside them
    double chords = 0;
    for (const double d : {1 + c[0], 1 - c[0], 1 + c[1], 1 - c[1]})
        chords += 2 * std::sqrt(r * r - d * d);
    // in 3D, a sphere about a point off the square's plane with the same circle on it
    const Ball flat = {c, r}, lifted = {move({c[0], c[1], 0.5}), std::sqrt(r * r + 0.25)};
    for (const int gdim : {2, 3})
        for (const char* source : {"P2", "analytic"})
        {
            Lines l;
            make_lines(gdim, 8, l);
            const Ball& ball = gdim == 2 ? flat : lifted;
            const LevelSetFunction<double, int> ls
                = std::string(source) == "P2" ? pk_level_set(l.view, 2, [&ball](const V3& y) { return ball(y); })
                                              : analytic_ls(ball);
            const LineParts p = line_parts(l.view, ls);
            const std::string what = "square's sides in " + std::to_string(gdim) + "D (" + source + ")";
            check(std::abs(p.inside - chords) <= 1e-12 * chords,
                  what + ": chords " + std::to_string(p.inside) + ", exact " + std::to_string(chords));
            check(std::abs(p.outside - (8 - chords)) <= 1e-12 * 8, what + ": the rest of the sides");
            check(p.points == 8, what + ": 8 crossing points of weight 1, got " + std::to_string(p.points));
            double off = 0;
            for (const V3& y : p.where)
                off = std::max(off, std::abs(ball(std::array<double, 3>{y[0], y[1], y[2]})));
            check(off <= 1e-12, what + ": crossing points on the circle");
        }
    // [-1, 1] in 1D: (x - 0.13)^2 - 0.6^2 < 0 on (-0.47, 0.73)
    {
        Lines l;
        make_lines(1, 8, l);
        const Ball ball = {{0.13, 0, 0}, 0.6};
        for (const char* source : {"P2", "analytic"})
        {
            const LevelSetFunction<double, int> ls
                = std::string(source) == "P2" ? pk_level_set(l.view, 2, [&ball](const V3& y) { return ball(y); })
                                              : analytic_ls(ball);
            const LineParts p = line_parts(l.view, ls);
            const std::string what = std::string("[-1, 1] (") + source + ")";
            check(std::abs(p.inside - 1.2) <= 1e-13 && std::abs(p.outside - 0.8) <= 1e-13,
                  what + ": lengths " + std::to_string(p.inside) + " and " + std::to_string(p.outside));
            check(p.points == 2 && p.where.size() == 2, what + ": 2 points");
            if (p.where.size() == 2)
            {
                const double lo = std::min(p.where[0][0], p.where[1][0]), hi = std::max(p.where[0][0], p.where[1][0]);
                check(std::abs(lo + 0.47) <= 1e-13 && std::abs(hi - 0.73) <= 1e-13, what + ": the points");
            }
        }
        // (x - 0.25) (x + 0.5) vanishes at two vertices and changes sign across
        // them: no cell is cut, and each point is the rule of its vertex, the
        // facet of the two intervals sharing it, one on either side. By
        // default each interval takes half of it; with owner the lower takes it.
        const LevelSetFunction<double, int> ls
            = pk_level_set(l.view, 2, [](const V3& y) { return (y[0] - 0.25) * (y[0] + 0.5); });
        const LineParts p = line_parts(l.view, ls);
        check(std::abs(p.inside - 0.75) <= 1e-13 && std::abs(p.outside - 1.25) <= 1e-13,
              "zero vertices: lengths " + std::to_string(p.inside) + " and " + std::to_string(p.outside));
        const std::vector<LevelSetFunction<double, int>> level_sets = {ls};
        const part::CutResult<double, int> result
            = part::cut<double, int>(l.view, std::span<const LevelSetFunction<double, int>>(level_sets));
        const part::MeshPart<double, int> zero = part::select(result, "phi = 0");
        check(result.cut_cells.empty(), "zero vertices: no cut cell");
        const std::vector<part::EntityRules<double>> on = part::entity_rules(zero, 4);
        bool sides = on.size() == 1 && on[0].n_rules() == 2 && on[0].weights == std::vector<double>{1.0, 1.0};
        for (int r = 0; sides && r < on[0].n_rules(); ++r)
        {
            const std::int32_t s0 = on[0].cell_offsets[r], s1 = on[0].cell_offsets[r + 1];
            sides = s1 - s0 == 2 && (on[0].negative[s0] ^ on[0].negative[s0 + 1]) == 1
                    && (on[0].positive[s0] ^ on[0].positive[s0 + 1]) == 1;
        }
        check(sides, "zero vertices: two vertex rules of weight 1, each with its two intervals on either side");
        const auto split = part::quadrature_rules(zero, 4, true);
        check(split._weights == std::vector<double>(4, 0.5), "zero vertices: half a point to each sharing interval");
        const auto owner = part::quadrature_rules(zero, 4, true, quadrays::Options{}, part::EntityPolicy::owner);
        const std::vector<V3> where = physical_points(l.view, owner);
        check(owner._weights == std::vector<double>{1.0, 1.0} && where.size() == 2
                  && std::abs(std::min(where[0][0], where[1][0]) + 0.5) <= 1e-15
                  && std::abs(std::max(where[0][0], where[1][0]) - 0.25) <= 1e-15,
              "zero vertices: with owner, the points -0.5 and 0.25 of weight 1");
    }
    std::printf("intervals: done\n");
}

} // namespace

int main()
{
    test_moved_cells();
    test_cube_surface();
    test_surface_entities();
    test_intervals();
    if (failures > 0)
    {
        std::printf("%d check(s) failed\n", failures);
        return 1;
    }
    std::printf("all checks passed\n");
    return 0;
}
