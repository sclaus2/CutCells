// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// The robustness cases (docs/quadrays/RESULTS.md, v1.2). Batch 1: spheres and
// planes placed on the test meshes so that they pass through vertices, touch
// faces, cut tiny caps or lie in faces, and scaled level sets. Batch 2: several
// components, close roots and singular points (two touching or nearly touching
// balls, thin shells, cones, a double root, a torus). Exact per-cell values
// where they are known, totals otherwise, and the checks every rule must pass.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <string>
#include <utility>
#include <vector>

#include <cutcells/quadrature.h>
#include <cutcells/quadrays/analytic.h>
#include <cutcells/quadrays/rules.h>

#include "exact_reference.h"
#include "test_mesh.h"

namespace cutcells::quadrays::support
{

/// Placements are made on the grid with n = 16 (h = 1/8), whose planes x = 2h,
/// 4h and 6h are also grid planes for n = 8.
inline constexpr double h16 = 2.0 / 16;

enum class Shape
{
    sphere,      ///< |x - centre|^2 - radius^2
    plane,       ///< normal . x - offset
    two_spheres, ///< product of two spheres with disjoint balls (centre, radius) and (centre2, radius2)
    shell,       ///< (|x - centre|^2 - radius^2)(|x - centre|^2 - radius2^2), radius < radius2
    cone,        ///< (x - cx)^2 + (y - cy)^2 - (z - cz)^2 / 4, apex at centre
    double_root, ///< (x - cx)^2: phi >= 0, zero on a plane with zero gradient
    torus        ///< (|x - c|^2 + R^2 - r^2)^2 - 4 R^2 ((x - cx)^2 + (y - cy)^2), R = radius, r = radius2
};

/// A level set phi, multiplied by scale.
struct Case
{
    std::string name;
    std::string description;
    Shape shape = Shape::sphere;
    V3 centre = {0, 0, 0};
    double radius = 0;
    V3 centre2 = {0, 0, 0};
    double radius2 = 0;
    V3 normal = {1, 0, 0};
    double offset = 0;
    double scale = 1;
};

inline std::vector<Case> all_cases()
{
    const V3 off = {0.0123, -0.0371, 0.0217};
    std::vector<Case> cases;
    auto add = [&](const char* name, const char* description, Shape shape) -> Case&
    {
        Case c;
        c.name = name;
        c.description = description;
        c.shape = shape;
        cases.push_back(c);
        return cases.back();
    };
    auto sphere = [&](const char* name, const char* description, const V3& centre, double radius, double scale = 1.0)
    {
        Case& c = add(name, description, Shape::sphere);
        c.centre = centre;
        c.radius = radius;
        c.scale = scale;
    };
    auto plane = [&](const char* name, const char* description, const V3& normal, double offset)
    {
        Case& c = add(name, description, Shape::plane);
        c.normal = normal;
        c.offset = offset;
    };
    auto two = [&](const char* name, const char* description, double gap)
    {
        // radius 0.3 each, touching (gap 0) at the off-grid point off
        Case& c = add(name, description, Shape::two_spheres);
        c.centre = {off[0] - 0.3, off[1], off[2]};
        c.centre2 = {off[0] + 0.3 + gap, off[1], off[2]};
        c.radius = c.radius2 = 0.3;
    };
    // batch 1: placements on the grid, scaling
    sphere("sphere", "off the grid's symmetry (baseline)", off, 0.7);
    sphere("sphere-vertex-tangent", "centre on a vertex, r = 4h: through 6 vertices, tangent to grid planes there",
           {0, 0, 0}, 4 * h16);
    sphere("sphere-vertices", "centre on a vertex, r = sqrt(32) h: through 12 vertices", {0, 0, 0},
           std::sqrt(32.0) * h16);
    sphere("sphere-face-tangent", "tangent to the grid plane x = 6h inside a face", off, 6 * h16 - off[0]);
    sphere("sphere-cap-1e-6", "cap of height 1e-6 beyond x = 6h", off, 6 * h16 - off[0] + 1e-6);
    sphere("sphere-cap-1e-12", "cap of height 1e-12 beyond x = 6h", off, 6 * h16 - off[0] + 1e-12);
    sphere("sphere-scale-1e-150", "baseline sphere, phi times 1e-150", off, 0.7, 1e-150);
    sphere("sphere-scale-1e+150", "baseline sphere, phi times 1e+150", off, 0.7, 1e150);
    sphere("sphere-scale-1e-200", "baseline sphere, phi times 1e-200", off, 0.7, 1e-200);
    sphere("sphere-scale-1e+200", "baseline sphere, phi times 1e+200", off, 0.7, 1e200);
    plane("plane", "x + 0.3 y - 0.2 z = 0.01 (baseline)", {1, 0.3, -0.2}, 0.01);
    plane("plane-faces", "x = 2h: phi = 0 on whole faces", {1, 0, 0}, 2 * h16);
    plane("plane-near-faces", "x = 2h + 1e-12: slivers 1e-12 thick", {1, 0, 0}, 2 * h16 + 1e-12);
    plane("plane-vertices", "x + y + z = 2h: through vertices", {1, 1, 1}, 2 * h16);
    plane("plane-tet-faces", "x - y = h: in slanted tet faces, through hex edges", {1, -1, 0}, h16);
    // batch 2: several components, close roots, singular points (degree 4 for products)
    two("spheres-touching", "two balls of radius 0.3 touching at an off-grid point (product, degree 4)", 0.0);
    two("spheres-gap-1e-6", "the same with a gap of 1e-6", 1e-6);
    two("spheres-gap-1e-10", "the same with a gap of 1e-10", 1e-10);
    for (const auto& [name, width] : {std::pair{"shell-1e-3", 1e-3}, std::pair{"shell-1e-6", 1e-6}})
    {
        Case& c = add(name, "", Shape::shell);
        c.description = std::string("spherical shell 0.5 < |x - c| < 0.5 + ") + (width > 1e-4 ? "1e-3" : "1e-6")
                        + " (product, degree 4)";
        c.centre = off;
        c.radius = 0.5;
        c.radius2 = 0.5 + width;
    }
    add("cone-vertex", "double cone, half-opening atan(1/2), apex on a grid vertex", Shape::cone).centre = {0, 0, 0};
    add("cone-cell", "the same cone, apex inside a cell", Shape::cone).centre = off;
    add("double-root", "phi = (x - 0.0123)^2: a plane where phi and grad phi vanish", Shape::double_root).centre = off;
    Case& torus = add("torus", "torus R = 0.5, r = 0.2 around the z axis through an off-grid centre (degree 4)",
                      Shape::torus);
    torus.centre = off;
    torus.radius = 0.5;
    torus.radius2 = 0.2;
    return cases;
}

inline int degree(const Case& c)
{
    switch (c.shape)
    {
    case Shape::plane:
        return 1;
    case Shape::two_spheres:
    case Shape::shell:
    case Shape::torus:
        return 4;
    default:
        return 2;
    }
}

inline double sq_dist(const V3& x, const V3& c)
{
    return (x[0] - c[0]) * (x[0] - c[0]) + (x[1] - c[1]) * (x[1] - c[1]) + (x[2] - c[2]) * (x[2] - c[2]);
}

/// phi of a case at x, for every scalar type of quadrays/taylor.h: double for
/// values, Dual and Taylor models for the analytic interface.
template <typename V>
V case_value(const Case& c, const std::array<V, 3>& x)
{
    auto sq = [&x](const V3& p)
    {
        const V a = x[0] - p[0], b = x[1] - p[1], d = x[2] - p[2];
        return a * a + b * b + d * d;
    };
    V v;
    switch (c.shape)
    {
    case Shape::sphere:
        v = sq(c.centre) - c.radius * c.radius;
        break;
    case Shape::plane:
        v = x[0] * c.normal[0] + x[1] * c.normal[1] + x[2] * c.normal[2] - c.offset;
        break;
    case Shape::two_spheres:
        v = (sq(c.centre) - c.radius * c.radius) * (sq(c.centre2) - c.radius2 * c.radius2);
        break;
    case Shape::shell:
        v = (sq(c.centre) - c.radius * c.radius) * (sq(c.centre) - c.radius2 * c.radius2);
        break;
    case Shape::cone:
    {
        const V a = x[0] - c.centre[0], b = x[1] - c.centre[1], d = x[2] - c.centre[2];
        v = a * a + b * b - 0.25 * (d * d);
        break;
    }
    case Shape::double_root:
    {
        const V a = x[0] - c.centre[0];
        v = a * a;
        break;
    }
    case Shape::torus:
    {
        const double R = c.radius, r = c.radius2;
        const V a = sq(c.centre) + (R * R - r * r), u = x[0] - c.centre[0], w = x[1] - c.centre[1];
        v = a * a - (u * u + w * w) * (4 * R * R);
        break;
    }
    }
    return v * c.scale;
}

inline double phi(const Case& c, const V3& x) { return case_value(c, x); }

/// A case as an algoim-style functor, for quadrays/analytic.h.
struct CaseLevelSet
{
    const Case* c = nullptr;

    template <typename V>
    V operator()(const std::array<V, 3>& x) const
    {
        return case_value(*c, x);
    }
};

inline double gradient_norm(const Case& c, const V3& x)
{
    V3 g = {0, 0, 0};
    const V3 d1 = {x[0] - c.centre[0], x[1] - c.centre[1], x[2] - c.centre[2]};
    switch (c.shape)
    {
    case Shape::sphere:
        g = {2 * d1[0], 2 * d1[1], 2 * d1[2]};
        break;
    case Shape::plane:
        g = c.normal;
        break;
    case Shape::two_spheres:
    case Shape::shell:
    {
        const V3& c2 = c.shape == Shape::shell ? c.centre : c.centre2;
        const double s1 = sq_dist(x, c.centre) - c.radius * c.radius, s2 = sq_dist(x, c2) - c.radius2 * c.radius2;
        for (int d = 0; d < 3; ++d)
            g[d] = s2 * 2 * d1[d] + s1 * 2 * (x[d] - c2[d]);
        break;
    }
    case Shape::cone:
        g = {2 * d1[0], 2 * d1[1], -0.5 * d1[2]};
        break;
    case Shape::double_root:
        g = {2 * d1[0], 0, 0};
        break;
    case Shape::torus:
    {
        const double R = c.radius, r = c.radius2, a = sq_dist(x, c.centre) + R * R - r * r;
        g = {4 * a * d1[0] - 8 * R * R * d1[0], 4 * a * d1[1] - 8 * R * R * d1[1], 4 * a * d1[2]};
        break;
    }
    }
    return std::abs(c.scale) * std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
}

/// Exact part measures of one cell, where known.
struct Reference
{
    bool known = true;    ///< false: only the totals are known (exact_totals)
    double volume = 0;    ///< phi < 0
    double area = 0;      ///< phi = 0, cut through the cell
    double face_area = 0; ///< phi = 0 on a face of the cell (no owner)
};

/// The cell's faces relative to c.
inline std::vector<exact::Face> shifted(const std::vector<exact::Face>& faces, const V3& c)
{
    std::vector<exact::Face> out = faces;
    const exact::V3 cc = {c[0], c[1], c[2]};
    for (exact::Face& f : out)
    {
        for (exact::V3& v : f.loop)
            v = exact::sub(v, cc);
        f.d -= exact::dot(f.n, cc);
    }
    return out;
}

/// Ball of radius r about centre in a cell (faces in absolute coordinates).
inline Reference ball(const TestCell& cell, const V3& centre, double r)
{
    Reference out;
    const double dc = std::sqrt(sq_dist(cell.centroid, centre));
    if (dc - cell.radius >= r)
        return out;
    if (dc + cell.radius <= r)
    {
        out.volume = cell.volume;
        return out;
    }
    const std::vector<exact::Face> faces = shifted(cell.faces, centre);
    out.area = exact::sphere_area(faces, r);
    out.volume = exact::ball_volume(faces, r, out.area);
    return out;
}

/// Exact values of a cell (faces in absolute coordinates).
inline Reference reference(const Case& c, const TestCell& cell)
{
    Reference r;
    switch (c.shape)
    {
    case Shape::sphere:
        return ball(cell, c.centre, c.radius);
    case Shape::plane:
    {
        const exact::PlaneCut cut = exact::plane_cut(cell.faces, c.normal, c.offset);
        r.volume = cut.volume_below;
        r.area = cut.cut_area;
        r.face_area = cut.face_in_plane;
        return r;
    }
    case Shape::two_spheres:
    {
        const Reference a = ball(cell, c.centre, c.radius), b = ball(cell, c.centre2, c.radius2);
        r.volume = a.volume + b.volume;
        r.area = a.area + b.area;
        return r;
    }
    case Shape::shell:
    {
        const Reference a = ball(cell, c.centre, c.radius), b = ball(cell, c.centre, c.radius2);
        r.volume = b.volume - a.volume;
        r.area = a.area + b.area;
        return r;
    }
    default:
        r.known = false;
        return r;
    }
}

/// Exact totals over [-1, 1]^3 (volume of phi < 0, area of phi = 0) where
/// per-cell values are not known.
inline std::array<double, 2> exact_totals(const Case& c)
{
    switch (c.shape)
    {
    case Shape::cone:
    {
        // two nappes of heights 1 -+ cz and radius half the height, cut by z = +-1
        const double h1 = 1 - c.centre[2], h2 = 1 + c.centre[2];
        return {M_PI / 12 * (h1 * h1 * h1 + h2 * h2 * h2), M_PI * std::sqrt(5.0) / 4 * (h1 * h1 + h2 * h2)};
    }
    case Shape::double_root:
        return {0.0, 4.0}; // the area of the plane, if one counts it
    case Shape::torus:
        return {2 * M_PI * M_PI * c.radius * c.radius2 * c.radius2, 4 * M_PI * M_PI * c.radius * c.radius2};
    default:
        return {0.0, 0.0};
    }
}

/// What the checks of one rule found.
struct RuleCheck
{
    bool fail = false; ///< non-finite points or weights
    int negative = 0;  ///< negative weights
    int outside = 0;   ///< points outside the cell, beyond 1e-12 in reference coordinates
    int side = 0;      ///< volume points with phi > 0, or interface points off phi = 0, beyond 1e-9 h
    double value = 0;  ///< sum of the weights
};

/// Check the rule of one cell; points are in its reference coordinates.
inline RuleCheck check_rule(const Case& c, const TestCell& cell, const quadrature::QuadratureRules<double>& rule,
                            bool surface, double h)
{
    RuleCheck out;
    for (std::size_t p = 0; p < rule._weights.size(); ++p)
    {
        const V3 xi = {rule._points[3 * p], rule._points[3 * p + 1], rule._points[3 * p + 2]};
        const double w = rule._weights[p];
        if (!std::isfinite(w) || !std::isfinite(xi[0]) || !std::isfinite(xi[1]) || !std::isfinite(xi[2]))
        {
            out.fail = true;
            return out;
        }
        out.value += w;
        out.negative += w < 0;
        out.outside += !in_reference_cell(cell.type, xi, 1e-12);
        const V3 x = physical(cell, xi);
        const double f = phi(c, x), g = gradient_norm(c, x);
        const double reach = 1e-9 * h * g; // |phi| within 1e-9 h of the interface
        out.side += surface ? !(std::abs(f) <= reach) : !(f <= reach);
    }
    return out;
}

/// One case on a mesh of [-1, 1]^3: what the checks found, totals and
/// per-cell errors.
struct CaseRun
{
    long fail = 0, negative = 0, outside = 0, side = 0;
    int max_bisections = 0; ///< in one cell
    int two_roots = 0;      ///< boxes certified with two roots per line
    double total = 0, exact_total = 0;
    double l1 = 0; ///< over the cells with exact values, without those with phi = 0 on a face
};

/// Rules of one part of case @p c for every cell of the mesh @p mesh of
/// [-1, 1]^3 with n cells per side (grid_cells), from the case's Bernstein
/// coefficients or as an analytic level set, and their checks. Only the cells
/// from index @p first on along every axis (first = n / 2: the octant x, y, z >= 0).
inline CaseRun run_case(const Case& c, const std::string& mesh, int n, const SelectionTerm& term, int q,
                        const Options& opt, bool analytic, int first = 0)
{
    const double h = 2.0 / n;
    const bool surface = part_of(term) == Part::interface;
    std::vector<double> coeffs;
    const CaseLevelSet functor = {&c};
    const AnalyticLevelSet phi_analytic = analytic_level_set(functor);
    CaseRun r;
    for (int i0 = first; i0 < n; ++i0)
        for (int i1 = first; i1 < n; ++i1)
            for (int i2 = first; i2 < n; ++i2)
                for (const TestCell& cell : grid_cells(mesh, {-1 + h * i0, -1 + h * i1, -1 + h * i2}, h, {0, 0, 0}))
                {
                    const Reference ref = reference(c, cell);
                    const bool on_face = surface && ref.face_area > 0;
                    if (ref.known)
                        r.exact_total += surface ? ref.area + 0.5 * ref.face_area : ref.volume;
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
                    r.two_roots += stats.two_roots;
                    const RuleCheck check = check_rule(c, cell, rule, surface, h);
                    if (check.fail)
                    {
                        ++r.fail;
                        continue;
                    }
                    r.negative += check.negative;
                    r.outside += check.outside;
                    r.side += check.side;
                    r.total += check.value;
                    if (ref.known && !on_face)
                        r.l1 += std::abs(check.value - (surface ? ref.area : ref.volume));
                }
    return r;
}

} // namespace cutcells::quadrays::support
