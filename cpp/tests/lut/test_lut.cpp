// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// The lookup-table backend (lut/) and the front end's backend "lut":
//  - the templates' quadrilaterals and hexahedra are boxes in Basix order;
//  - lut::cut_cell on triangles, quadrilaterals, tetrahedra and hexahedra,
//    templates of order 1 to 3, against exact measures of planes (a closed
//    form on simplices, Kuhn simplices for quadrilaterals and hexahedra), whole
//    and split into simplices (classical and midpoint), in double and float;
//  - two planes in one cell: the four sides add up to the cell, the zero
//    set of one splits into its two sides of the other, and the curve where
//    both vanish has its exact length;
//  - planes through template vertices and faces: the interface is counted
//    once;
//  - parts of the analytic sphere: errors fall as (h / k)^2 with the
//    template order k, and parts of a sphere and a plane crossing in cells
//    add up.
// Exits non-zero on failure.

#include <cutcells/cell_types.h>
#include <cutcells/level_set.h>
#include <cutcells/lut/cell_pieces.h>
#include <cutcells/lut/iso_refine.h>
#include <cutcells/lut/piece_rules.h>
#include <cutcells/lut/triangulation.h>
#include <cutcells/part/cut_result.h>
#include <cutcells/part/mesh_part.h>
#include <cutcells/part/output.h>
#include <cutcells/quadrays/analytic.h>
#include <cutcells/reference_cell.h>

#include <cmath>
#include <cstdio>
#include <functional>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "../part/support/box_mesh.h"
#include "../quadrays/support/sphere_functors.h"

using namespace cutcells;
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

// ============================================================================
// Exact measures of a plane in a simplex
// ============================================================================

using P3 = std::array<double, 3>;

P3 sub(const P3& a, const P3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }

P3 cross(const P3& a, const P3& b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

double norm(const P3& a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }

double triangle_measure(const P3& a, const P3& b, const P3& c) { return 0.5 * norm(cross(sub(b, a), sub(c, a))); }

double tet_volume(const P3& a, const P3& b, const P3& c, const P3& d)
{
    const P3 n = cross(sub(b, a), sub(c, a)), e = sub(d, a);
    return std::abs(n[0] * e[0] + n[1] * e[1] + n[2] * e[2]) / 6;
}

/// The measures of {phi < 0} and {phi = 0} in a triangle or tetrahedron (2D
/// points get z = 0) for the linear phi with values v at its vertices (none
/// exactly 0), from the clipped simplex itself: no division between close
/// values, unlike the divided-difference formula.
std::pair<double, double> simplex_measures(const std::vector<P3>& x, const std::vector<double>& v)
{
    const int n = static_cast<int>(x.size()) - 1;
    std::vector<int> below, above;
    for (int i = 0; i <= n; ++i)
        (v[static_cast<std::size_t>(i)] < 0 ? below : above).push_back(i);
    auto root = [&](int a, int b)
    {
        const double t = v[static_cast<std::size_t>(a)] / (v[static_cast<std::size_t>(a)] - v[static_cast<std::size_t>(b)]);
        P3 p;
        for (int d = 0; d < 3; ++d)
            p[d] = x[static_cast<std::size_t>(a)][d] + t * (x[static_cast<std::size_t>(b)][d] - x[static_cast<std::size_t>(a)][d]);
        return p;
    };
    const double volume = n == 2 ? triangle_measure(x[0], x[1], x[2]) : tet_volume(x[0], x[1], x[2], x[3]);
    if (below.empty())
        return {0.0, 0.0};
    if (above.empty())
        return {volume, 0.0};
    // the corner cut off at a lone vertex
    auto corner = [&](int lone, const std::vector<int>& others)
    {
        std::vector<P3> r;
        for (const int o : others)
            r.push_back(root(lone, o));
        const double interface = n == 2 ? norm(sub(r[1], r[0])) : triangle_measure(r[0], r[1], r[2]);
        const double piece = n == 2 ? triangle_measure(x[static_cast<std::size_t>(lone)], r[0], r[1])
                                    : tet_volume(x[static_cast<std::size_t>(lone)], r[0], r[1], r[2]);
        return std::pair<double, double>{piece, interface};
    };
    if (below.size() == 1)
        return corner(below[0], above);
    if (above.size() == 1)
    {
        const auto [piece, interface] = corner(above[0], below);
        return {volume - piece, interface};
    }
    // two vertices on each side of a tetrahedron: a wedge between them
    const int i = below[0], j = below[1], p = above[0], q = above[1];
    const P3 ip = root(i, p), iq = root(i, q), jp = root(j, p), jq = root(j, q);
    const P3 &xi = x[static_cast<std::size_t>(i)], &xj = x[static_cast<std::size_t>(j)];
    // wedge (i, ip, iq | j, jp, jq) in three tetrahedra
    const double wedge = tet_volume(xi, iq, ip, xj) + tet_volume(ip, xj, jq, jp) + tet_volume(ip, iq, jq, xj);
    return {wedge, triangle_measure(ip, iq, jq) + triangle_measure(ip, jq, jp)};
}

/// Exact measures of {g.x + c < 0} and {g.x + c = 0} in a reference cell,
/// from its Kuhn simplices.
std::pair<double, double> exact_measures(cell::type type, const std::vector<double>& g, double c)
{
    const int tdim = cell::get_tdim(type);
    const std::vector<double> ref = cell::reference_vertices<double>(type);
    std::vector<std::vector<int>> simplices;
    if (type == cell::type::triangle || type == cell::type::tetrahedron)
    {
        simplices.emplace_back();
        for (int v = 0; v <= tdim; ++v)
            simplices.back().push_back(v);
    }
    else
    {
        int ids[8] = {0, 1, 2, 3, 4, 5, 6, 7};
        cell::triangulation(type, ids, simplices);
    }
    double below = 0, interface = 0;
    for (const std::vector<int>& s : simplices)
    {
        std::vector<P3> x;
        std::vector<double> v;
        for (const int i : s)
        {
            P3 p = {0, 0, 0};
            double value = c;
            for (int d = 0; d < tdim; ++d)
            {
                p[d] = ref[static_cast<std::size_t>(i * tdim + d)];
                value += g[static_cast<std::size_t>(d)] * p[d];
            }
            x.push_back(p);
            v.push_back(value);
        }
        const auto [b, a] = simplex_measures(x, v);
        below += b;
        interface += a;
    }
    return {below, interface};
}

double reference_volume(cell::type type)
{
    return type == cell::type::triangle ? 0.5 : (type == cell::type::tetrahedron ? 1.0 / 6.0 : 1.0);
}

/// The measure of the pieces that @p keep selects (by index), with rules of
/// degree 1 in the reference cell.
double measure(cell::type type, const lut::Pieces<double>& pieces, const std::function<bool(std::size_t)>& keep)
{
    lut::CellMap<double> map;
    map.type = type;
    map.gdim = cell::get_tdim(type);
    map.vertices = cell::reference_vertices<double>(type);
    const std::size_t tdim = static_cast<std::size_t>(pieces.tdim);
    std::vector<double> points, weights;
    for (int p = 0; p < pieces.n_pieces(); ++p)
    {
        if (!keep(static_cast<std::size_t>(p)))
            continue;
        const std::size_t first = static_cast<std::size_t>(pieces.offsets[p]) * tdim;
        const std::size_t size = static_cast<std::size_t>(pieces.offsets[p + 1] - pieces.offsets[p]) * tdim;
        lut::append_piece_rule(map, pieces.types[static_cast<std::size_t>(p)],
                               std::span<const double>(pieces.vertices).subspan(first, size), 1, points, weights);
    }
    double sum = 0;
    for (const double w : weights)
        sum += w;
    return sum;
}

double measure_below(cell::type type, const lut::Pieces<double>& pieces)
{
    return measure(type, pieces, [&](std::size_t p) { return pieces.zero[p] == 0 && (pieces.negative[p] & 1); });
}

double measure_above(cell::type type, const lut::Pieces<double>& pieces)
{
    return measure(type, pieces, [&](std::size_t p) { return pieces.zero[p] == 0 && (pieces.positive[p] & 1); });
}

double measure_zero(cell::type type, const lut::Pieces<double>& pieces)
{
    return measure(type, pieces, [&](std::size_t p) { return pieces.zero[p] == 1; });
}

/// The values of g.x + c at the template's vertices.
std::vector<double> template_values(cell::type type, int k, const std::vector<double>& g, double c)
{
    const int tdim = cell::get_tdim(type);
    const std::span<const double> tv = lut::template_vertices(type, k);
    std::vector<double> values;
    for (std::size_t v = 0; v < tv.size() / static_cast<std::size_t>(tdim); ++v)
    {
        double value = c;
        for (int d = 0; d < tdim; ++d)
            value += g[static_cast<std::size_t>(d)] * tv[v * static_cast<std::size_t>(tdim) + d];
        values.push_back(value);
    }
    return values;
}

const std::vector<cell::type> cell_types
    = {cell::type::triangle, cell::type::quadrilateral, cell::type::tetrahedron, cell::type::hexahedron};

/// Two planes in one cell: the curve where both vanish (a point in 2D, of
/// measure 1), on templates of order 1 to 3. In the tetrahedron the line
/// x = 0.3, y = 0.2 lies on x + y = 0.5, a face between sub-cells of the
/// order 2 template, and is counted once.
void test_curves()
{
    struct Case
    {
        cell::type type;
        std::vector<double> ga;
        double ca;
        std::vector<double> gb;
        double cb;
        double exact;
    };
    const std::vector<Case> cases = {
        {cell::type::hexahedron, {1, 0, 0.2}, -0.4, {-0.1, 1, 0}, -0.5, std::sqrt(1 + 0.04 + 0.0004)},
        {cell::type::tetrahedron, {1, 0, 0}, -0.3, {0, 1, 0}, -0.2, 0.5},
        {cell::type::triangle, {1, 0}, -0.3, {0, 1}, -0.2, 1.0},
        {cell::type::quadrilateral, {1, 0}, -0.3, {0, 1}, -0.6, 1.0},
    };
    for (const Case& c : cases)
    {
        double worst = 0;
        for (int k = 1; k <= 3; ++k)
        {
            std::vector<double> values = template_values(c.type, k, c.ga, c.ca);
            const std::vector<double> vb = template_values(c.type, k, c.gb, c.cb);
            values.insert(values.end(), vb.begin(), vb.end());
            lut::Pieces<double> pieces;
            lut::cut_cell<double>(c.type, k, values, 2, 3, true, cell::TriangulationStrategy::none, pieces);
            const double curve = measure(c.type, pieces, [&](std::size_t p) { return pieces.zero[p] == 3; });
            worst = std::max(worst, std::abs(curve - c.exact));
        }
        std::printf("curve of two planes in a %s: worst error %.1e\n", cell::cell_type_to_str(c.type).c_str(), worst);
        check(worst < 1e-13, "curve of two planes in a " + cell::cell_type_to_str(c.type));
    }
}

/// The sub-cells of the quadrilateral and hexahedron templates are boxes with
/// vertex 0 at their lower corner, in Basix order, as cut_cell assumes.
void test_templates()
{
    bool boxes = true;
    for (const cell::type type : {cell::type::quadrilateral, cell::type::hexahedron})
        for (int k = 1; k <= 4; ++k)
        {
            const IsoRefineTemplate& tpl = iso_p1_template(type, k);
            const int tdim = tpl.tdim, n = tpl.vertices_per_cell;
            for (int c = 0; c < tpl.n_cells; ++c)
            {
                std::vector<double> x;
                for (int j = 0; j < n; ++j)
                    for (int d = 0; d < tdim; ++d)
                        x.push_back(tpl.ref_vertex_coords[static_cast<std::size_t>(
                            tpl.cell_connectivity[static_cast<std::size_t>(c * n + j)] * tdim + d)]);
                boxes &= lut::is_parallelotope(std::span<const double>(x), tdim, tdim);
                for (int d = 0; d < tdim; ++d)
                    boxes &= x[static_cast<std::size_t>(((1 << tdim) - 1) * tdim + d)] > x[static_cast<std::size_t>(d)];
            }
        }
    check(boxes, "the templates' sub-cells are boxes in Basix order");
}

/// Random planes through each cell type on templates of order 1 to 3.
void test_planes()
{
    std::mt19937 rng(3);
    std::normal_distribution<double> normal(0.0, 1.0);
    for (const cell::type type : cell_types)
    {
        const int tdim = cell::get_tdim(type);
        const std::vector<double> ref = cell::reference_vertices<double>(type);
        const double volume = reference_volume(type);
        double worst = 0;
        for (int k = 1; k <= 3; ++k)
            for (const cell::TriangulationStrategy triangulation :
                 {cell::TriangulationStrategy::none, cell::TriangulationStrategy::classical,
                  cell::TriangulationStrategy::midpoint})
                for (int trial = 0; trial < 40; ++trial)
                {
                    std::vector<double> g(static_cast<std::size_t>(tdim));
                    double c = 0;
                    for (int d = 0; d < tdim; ++d)
                    {
                        g[static_cast<std::size_t>(d)] = normal(rng);
                        c -= g[static_cast<std::size_t>(d)] * 0.3;
                    }
                    c += 0.2 * normal(rng);
                    lut::Pieces<double> pieces;
                    lut::cut_cell<double>(type, k, template_values(type, k, g, c), 1, 1, false, triangulation, pieces);
                    const auto [below, area] = exact_measures(type, g, c);
                    worst = std::max({worst, std::abs(measure_below(type, pieces) - below),
                                      std::abs(measure_above(type, pieces) - (volume - below)),
                                      std::abs(measure_zero(type, pieces) - area)});
                }
        std::printf("planes in a %s: worst error %.1e\n", cell::cell_type_to_str(type).c_str(), worst);
        check(worst < 1e-13, "planes in a " + cell::cell_type_to_str(type));
    }
}

/// The float instantiation on one plane per cell type.
void test_float()
{
    double worst = 0;
    for (const cell::type type : cell_types)
    {
        const int tdim = cell::get_tdim(type);
        std::vector<double> g = {0.8, -0.5, 0.3};
        g.resize(static_cast<std::size_t>(tdim));
        const double c = -0.21;
        const std::vector<double> v = template_values(type, 2, g, c);
        lut::Pieces<float> pieces;
        lut::cut_cell<float>(type, 2, std::vector<float>(v.begin(), v.end()), 1, 0, false, cell::TriangulationStrategy::none,
                              pieces);
        lut::CellMap<float> map;
        map.type = type;
        map.gdim = tdim;
        const std::vector<double> ref = cell::reference_vertices<double>(type);
        map.vertices.assign(ref.begin(), ref.end());
        std::vector<float> points, weights;
        for (int p = 0; p < pieces.n_pieces(); ++p)
            if (pieces.negative[static_cast<std::size_t>(p)] & 1)
            {
                const std::size_t first = static_cast<std::size_t>(pieces.offsets[p] * tdim);
                const std::size_t size = static_cast<std::size_t>((pieces.offsets[p + 1] - pieces.offsets[p]) * tdim);
                lut::append_piece_rule<float>(map, pieces.types[static_cast<std::size_t>(p)],
                                              std::span<const float>(pieces.vertices).subspan(first, size), 1, points,
                                              weights);
            }
        double below = 0;
        for (const float w : weights)
            below += w;
        worst = std::max(worst, std::abs(below - exact_measures(type, g, c).first));
    }
    std::printf("float: worst error %.1e\n", worst);
    check(worst < 1e-6, "lut::cut_cell<float>");
}

/// Two planes in one cell: the four sides add up to the cell, and the zero set
/// of the first splits into its parts on both sides of the second.
void test_two_planes()
{
    std::mt19937 rng(5);
    std::normal_distribution<double> normal(0.0, 1.0);
    for (const cell::type type : cell_types)
    {
        const int tdim = cell::get_tdim(type);
        const double volume = reference_volume(type);
        double worst = 0;
        for (int k = 1; k <= 2; ++k)
            for (int trial = 0; trial < 60; ++trial)
            {
                std::vector<double> values;
                std::vector<double> ga(static_cast<std::size_t>(tdim));
                double ca = 0;
                for (int l = 0; l < 2; ++l)
                {
                    std::vector<double> g(static_cast<std::size_t>(tdim));
                    double c = 0.2 * normal(rng);
                    for (int d = 0; d < tdim; ++d)
                    {
                        g[static_cast<std::size_t>(d)] = normal(rng);
                        c -= g[static_cast<std::size_t>(d)] * 0.4;
                    }
                    if (l == 0)
                    {
                        ga = g;
                        ca = c;
                    }
                    const std::vector<double> v = template_values(type, k, g, c);
                    values.insert(values.end(), v.begin(), v.end());
                }
                lut::Pieces<double> pieces;
                lut::cut_cell<double>(type, k, values, 2, 1, false, cell::TriangulationStrategy::none, pieces);
                double sides = 0, below_first = 0;
                for (const std::uint64_t below_bits : {0, 1, 2, 3})
                {
                    const double w = measure(type, pieces,
                                             [&](std::size_t p)
                                             {
                                                 return pieces.zero[p] == 0 && pieces.negative[p] == below_bits
                                                        && pieces.positive[p] == (3 & ~below_bits);
                                             });
                    sides += w;
                    below_first += (below_bits & 1) ? w : 0;
                }
                // the zero set of the first plane on both sides of the second
                const double zero_split
                    = measure(type, pieces, [&](std::size_t p) { return pieces.zero[p] == 1 && (pieces.negative[p] & 2); })
                      + measure(type, pieces,
                                [&](std::size_t p) { return pieces.zero[p] == 1 && (pieces.positive[p] & 2); });
                const auto [below, area] = exact_measures(type, ga, ca);
                worst = std::max({worst, std::abs(sides - volume), std::abs(below_first - below),
                                  std::abs(zero_split - area)});
            }
        std::printf("two planes in a %s: worst error %.1e\n", cell::cell_type_to_str(type).c_str(), worst);
        check(worst < 1e-13, "two planes in a " + cell::cell_type_to_str(type));
    }
}

/// The measures of {x_0 < t} and {x_0 = t} in a reference cell.
std::pair<double, double> axis_measures(cell::type type, double t)
{
    switch (type)
    {
    case cell::type::triangle:
        return {0.5 - 0.5 * (1 - t) * (1 - t), 1 - t};
    case cell::type::tetrahedron:
        return {(1 - (1 - t) * (1 - t) * (1 - t)) / 6, 0.5 * (1 - t) * (1 - t)};
    default:
        return {t, 1.0};
    }
}

/// Planes through vertices and faces of the template's sub-cells: values of
/// exactly 0 count as positive, so the interface is counted once.
void test_planes_through_vertices()
{
    for (const cell::type type : cell_types)
    {
        const int tdim = cell::get_tdim(type);
        const double volume = reference_volume(type);
        double worst = 0;
        for (int k = 1; k <= 4; ++k)
            for (int j = 1; j < 2 * k; ++j)
            {
                // x_0 = j / (2k): through template faces for even j; and its mirror image
                for (const double sign : {1.0, -1.0})
                {
                    std::vector<double> g(static_cast<std::size_t>(tdim), 0.0);
                    g[0] = sign;
                    const double c = -sign * j / (2.0 * k);
                    lut::Pieces<double> pieces;
                    lut::cut_cell<double>(type, k, template_values(type, k, g, c), 1, 1, false, cell::TriangulationStrategy::none,
                                          pieces);
                    const auto [lower, area] = axis_measures(type, j / (2.0 * k));
                    const double below = sign > 0 ? lower : volume - lower;
                    worst = std::max({worst, std::abs(measure_below(type, pieces) - below),
                                      std::abs(measure_above(type, pieces) - (volume - below)),
                                      std::abs(measure_zero(type, pieces) - area)});
                }
            }
        std::printf("planes through template vertices in a %s: worst error %.1e\n",
                    cell::cell_type_to_str(type).c_str(), worst);
        check(worst < 1e-13, "planes through template vertices in a " + cell::cell_type_to_str(type));
    }
}

// ============================================================================
// The backend "lut" of the front end
// ============================================================================

template <typename F>
LevelSetFunction<double, int> analytic_ls(const F& functor, const std::string& name)
{
    auto phi = std::make_shared<const quadrays::AnalyticLevelSet>(quadrays::analytic_level_set(functor));
    return create_level_set_function<double, int>(phi, 3, name);
}

double total(const part::MeshPart<double, int>& p, int k, bool full)
{
    lut::Options options;
    options.template_order = k;
    double s = 0;
    for (const double w : part::quadrature_rules(p, 1, full, options)._weights)
        s += w;
    return s;
}

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

/// The analytic sphere: errors fall as (h / k)^2.
void test_template_order()
{
    const V3 centre = {0.0123, -0.0371, 0.0217};
    const double radius = 0.7;
    const double ball = 4.0 / 3.0 * M_PI * radius * radius * radius, sphere = 4.0 * M_PI * radius * radius;
    const quadrays::support::SphereDistance distance = {centre, radius};
    for (const char* kind : {"hex", "tet"})
    {
        BoxMesh mesh;
        make_box_mesh(kind, 8, centre, mesh);
        const std::vector<LevelSetFunction<double, int>> ls = {analytic_ls(distance, "phi")};
        const part::CutResult<double, int> r = part::cut<double, int>(mesh.view, ls);
        double volume_error[5], area_error[5];
        for (int k = 1; k <= 4; ++k)
        {
            volume_error[k] = std::abs(total(part::select(r, "phi < 0"), k, true) / ball - 1);
            area_error[k] = std::abs(total(part::select(r, "phi = 0"), k, false) / sphere - 1);
        }
        const double volume_rate = std::log(volume_error[1] / volume_error[4]) / std::log(4.0);
        const double area_rate = std::log(area_error[1] / area_error[4]) / std::log(4.0);
        std::printf("sphere %s, k = 1 to 4: volume %.1e to %.1e (rate %.2f), area %.1e to %.1e (rate %.2f)\n", kind,
                    volume_error[1], volume_error[4], volume_rate, area_error[1], area_error[4], area_rate);
        check(std::abs(volume_rate - 2) < 0.1 && std::abs(area_rate - 2) < 0.1,
              std::string("errors fall as (h / k)^2, ") + kind);
    }
}

/// A sphere and a plane crossing in cells: the four parts add up to the box,
/// and the plane's interface splits into its parts inside and outside the ball.
void test_crossing_level_sets()
{
    const V3 centre = {0.0123, -0.0371, 0.0217};
    const quadrays::support::SphereDistance distance = {centre, 0.7};
    const AxisPlane plane = {0, 0.1};
    for (const char* kind : {"hex", "tet"})
    {
        BoxMesh mesh;
        make_box_mesh(kind, 8, centre, mesh);
        const std::vector<LevelSetFunction<double, int>> ls = {analytic_ls(distance, "a"), analytic_ls(plane, "b")};
        const part::CutResult<double, int> r = part::cut<double, int>(mesh.view, ls);
        double sum = 0;
        for (const char* e : {"a < 0 and b < 0", "a < 0 and b > 0", "a > 0 and b < 0", "a > 0 and b > 0"})
            sum += total(part::select(r, e), 2, true);
        const double union_ = total(part::select(r, "a < 0 or b < 0"), 2, true);
        const double outside = total(part::select(r, "a > 0 and b > 0"), 2, true);
        const double plane_inside = total(part::select(r, "b = 0 and a < 0"), 2, false);
        const double plane_outside = total(part::select(r, "b = 0 and a > 0"), 2, false);
        std::printf("sphere and plane %s: four parts - 8 = %.1e, union + rest - 8 = %.1e, plane %.15f\n", kind, sum - 8,
                    union_ + outside - 8, plane_inside + plane_outside);
        check(std::abs(sum - 8) < 1e-12 && std::abs(union_ + outside - 8) < 1e-12
                  && std::abs(plane_inside + plane_outside - 4) < 1e-12,
              std::string("parts of crossing level sets add up, ") + kind);
    }
}

} // namespace

int main()
{
    test_templates();
    test_planes();
    test_float();
    test_two_planes();
    test_curves();
    test_planes_through_vertices();
    test_template_order();
    test_crossing_level_sets();
    return failures == 0 ? 0 : 1;
}
