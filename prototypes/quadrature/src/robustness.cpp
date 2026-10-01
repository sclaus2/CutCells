// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Robustness study. Batch 1: spheres and planes placed on the test meshes so that
// they pass through vertices, touch faces, cut tiny caps or lie in faces, and scaled
// level sets. Batch 2: several components, close roots and singular points (two
// touching or nearly touching balls, thin shells, cones, a double root, a torus).
// Every cell of [-1, 1]^3 (n = 16) goes to every generator, cut or not, and is
// checked against exact per-cell values; for cones, the double root and the torus
// only the totals are known, so L1, max and bad stay empty.
//
// Per run (case, mesh, q, generator, part):
//   fail   cells where the generator threw, or returned non-finite points or weights
//   neg    negative weights
//   out    points outside the cell (beyond 1e-12 in reference coordinates)
//   side   volume points on the wrong side of the level set, or interface points off
//          it, by more than 1e-9 h
//   L1     sum over cells of |error|, relative to the exact total
//   max    largest cell error, relative to h^3 (volume) or h^2 (interface)
//   bad    cells with an error above 1e-4 h^3 (volume) or 1e-4 h^2 (interface)
//   total  error of the total over the mesh, relative to the exact total
//   us     mean and largest time per cell, over all cells
// The interface of a plane lying in faces has no well-defined owner: those cells are
// left out of L1, max and bad, and the exact total counts each such face once.
//
// Usage: quadrature_robustness [--case a,b] [--mesh tet,hex] [--n 16] [--q 3] [--gen a,b]
//                              [--part "phi < 0"]... [--list]

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "certified.h"
#include "clipped_box.h"
#include "exact_reference.h"
#include "generators.h"
#include "selection_expr.h"
#include "test_mesh.h"

using namespace cutcells;
using namespace cutcells::proto;

namespace
{
/// Placements are made on the grid with n = 16 (h = 1/8), whose planes x = 2h, 4h
/// and 6h are also grid planes for n = 8. --n changes the mesh only.
constexpr double h16 = 2.0 / 16;
int n_cells = 16;           ///< cells per side of [-1, 1]^3
double h = 2.0 / n_cells;   ///< grid spacing

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
    Vec3 centre = {0, 0, 0};
    double radius = 0;
    Vec3 centre2 = {0, 0, 0};
    double radius2 = 0;
    Vec3 normal = {1, 0, 0};
    double offset = 0;
    double scale = 1;
};

std::vector<Case> all_cases()
{
    const Vec3 off = {0.0123, -0.0371, 0.0217};
    std::vector<Case> cases;
    auto add = [&](const char* name, const char* description, Shape shape)
    {
        Case c;
        c.name = name;
        c.description = description;
        c.shape = shape;
        cases.push_back(c);
        return &cases.back();
    };
    auto sphere = [&](const char* name, const char* description, const Vec3& centre, double radius, double scale = 1.0)
    {
        Case* c = add(name, description, Shape::sphere);
        c->centre = centre;
        c->radius = radius;
        c->scale = scale;
    };
    auto plane = [&](const char* name, const char* description, const Vec3& normal, double offset)
    {
        Case* c = add(name, description, Shape::plane);
        c->normal = normal;
        c->offset = offset;
    };
    auto two = [&](const char* name, const char* description, double gap)
    {
        // radius 0.3 each, touching (gap 0) at the off-grid point off
        Case* c = add(name, description, Shape::two_spheres);
        c->centre = {off[0] - 0.3, off[1], off[2]};
        c->centre2 = {off[0] + 0.3 + gap, off[1], off[2]};
        c->radius = c->radius2 = 0.3;
    };
    // batch 1: placements on the grid, scaling
    sphere("sphere", "off the grid's symmetry (baseline)", off, 0.7);
    sphere("sphere-vertex-tangent", "centre on a vertex, r = 4h: through 6 vertices, tangent to grid planes there",
           {0, 0, 0}, 4 * h16);
    sphere("sphere-vertices", "centre on a vertex, r = sqrt(32) h: through 12 vertices", {0, 0, 0}, std::sqrt(32.0) * h16);
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
        Case* c = add(name, "", Shape::shell);
        c->description = std::string("spherical shell 0.5 < |x - c| < 0.5 + ") + (width > 1e-4 ? "1e-3" : "1e-6") + " (product, degree 4)";
        c->centre = off;
        c->radius = 0.5;
        c->radius2 = 0.5 + width;
    }
    Case* cone = add("cone-vertex", "double cone, half-opening atan(1/2), apex on a grid vertex", Shape::cone);
    cone->centre = {0, 0, 0};
    cone = add("cone-cell", "the same cone, apex inside a cell", Shape::cone);
    cone->centre = off;
    Case* dr = add("double-root", "phi = (x - 0.0123)^2: a plane where phi and grad phi vanish", Shape::double_root);
    dr->centre = off;
    Case* torus = add("torus", "torus R = 0.5, r = 0.2 around the z axis through an off-grid centre (degree 4)", Shape::torus);
    torus->centre = off;
    torus->radius = 0.5;
    torus->radius2 = 0.2;
    return cases;
}

int degree(const Case& c)
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

double sq_dist(const Vec3& x, const Vec3& c)
{
    return (x[0] - c[0]) * (x[0] - c[0]) + (x[1] - c[1]) * (x[1] - c[1]) + (x[2] - c[2]) * (x[2] - c[2]);
}

double phi(const Case& c, const Vec3& x)
{
    double v = 0;
    switch (c.shape)
    {
    case Shape::sphere:
        v = sq_dist(x, c.centre) - c.radius * c.radius;
        break;
    case Shape::plane:
        v = c.normal[0] * x[0] + c.normal[1] * x[1] + c.normal[2] * x[2] - c.offset;
        break;
    case Shape::two_spheres:
        v = (sq_dist(x, c.centre) - c.radius * c.radius) * (sq_dist(x, c.centre2) - c.radius2 * c.radius2);
        break;
    case Shape::shell:
        v = (sq_dist(x, c.centre) - c.radius * c.radius) * (sq_dist(x, c.centre) - c.radius2 * c.radius2);
        break;
    case Shape::cone:
        v = (x[0] - c.centre[0]) * (x[0] - c.centre[0]) + (x[1] - c.centre[1]) * (x[1] - c.centre[1])
            - 0.25 * (x[2] - c.centre[2]) * (x[2] - c.centre[2]);
        break;
    case Shape::double_root:
        v = (x[0] - c.centre[0]) * (x[0] - c.centre[0]);
        break;
    case Shape::torus:
    {
        const double R = c.radius, r = c.radius2, a = sq_dist(x, c.centre) + R * R - r * r;
        v = a * a - 4 * R * R * ((x[0] - c.centre[0]) * (x[0] - c.centre[0]) + (x[1] - c.centre[1]) * (x[1] - c.centre[1]));
        break;
    }
    }
    return c.scale * v;
}

double gradient_norm(const Case& c, const Vec3& x)
{
    Vec3 g = {0, 0, 0};
    const Vec3 d1 = {x[0] - c.centre[0], x[1] - c.centre[1], x[2] - c.centre[2]};
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
        const Vec3& c2 = c.shape == Shape::shell ? c.centre : c.centre2;
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
    bool known = true;     ///< false: only the totals are known (exact_totals)
    double volume = 0;     ///< phi < 0
    double area = 0;       ///< phi = 0, cut through the cell
    double face_area = 0;  ///< phi = 0 on a face of the cell (no owner)
};

/// The cell's faces relative to c.
std::vector<exact::Face> shifted(const std::vector<exact::Face>& faces, const Vec3& c)
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
Reference ball(const TestCell& cell, const Vec3& centre, double r)
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

Reference reference(const Case& c, const TestCell& cell)
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

/// Exact totals over [-1, 1]^3 (volume of phi < 0, area of phi = 0) where per-cell
/// values are not known.
std::array<double, 2> exact_totals(const Case& c)
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

struct Tally
{
    long cut = 0, fail = 0, negative = 0, outside = 0, side = 0, bad = 0, points = 0;
    long bisections = 0, uncertified = 0;
    double l1 = 0, max_error = 0, total = 0, exact_total = 0, seconds = 0, max_seconds = 0;
    long cells = 0;
};

std::vector<std::string> split_list(const std::string& s)
{
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ','))
        if (!item.empty())
            out.push_back(item);
    return out;
}
} // namespace

int main(int argc, char** argv)
{
    std::vector<std::string> case_names, meshes = {"tet", "hex"}, generators
                                                                    = {"certify", "certify:0.1", "algoim-auto", "algoim-gl", "quadgen"};
    std::vector<std::string> parts;
    std::vector<int> qs = {3};
    const std::vector<Case> cases = all_cases();
    try
    {
        for (int i = 1; i < argc; ++i)
        {
            const std::string a = argv[i];
            auto next = [&]() -> std::string
            {
                if (i + 1 >= argc)
                    throw std::runtime_error("missing value after " + a);
                return argv[++i];
            };
            if (a == "--case")
                case_names = split_list(next());
            else if (a == "--mesh")
                meshes = split_list(next());
            else if (a == "--gen")
                generators = split_list(next());
            else if (a == "--part")
                parts.push_back(next());
            else if (a == "--n")
            {
                n_cells = std::stoi(next());
                h = 2.0 / n_cells;
            }
            else if (a == "--q")
            {
                qs.clear();
                for (const std::string& s : split_list(next()))
                    qs.push_back(std::stoi(s));
            }
            else if (a == "--list")
            {
                for (const Case& c : cases)
                    std::printf("%-22s %s\n", c.name.c_str(), c.description.c_str());
                return 0;
            }
            else
                throw std::runtime_error("unknown argument: " + a);
        }
    }
    catch (const std::exception& e)
    {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }
    if (parts.empty())
        parts = {"phi < 0", "phi = 0"};

    std::printf("%-22s %-4s %2s %-12s %-8s | %5s | %4s %4s %4s %5s | %8s %8s %5s | %8s | %7s | %7s %8s | %s\n", "case", "mesh",
                "q", "generator", "part", "cut", "fail", "neg", "out", "side", "L1", "max", "bad", "total", "pts/cut",
                "us/cell", "max us", "bisect/uncert");
    for (const Case& c : cases)
    {
        if (!case_names.empty() && std::find(case_names.begin(), case_names.end(), c.name) == case_names.end())
            continue;
        LevelSet ls;
        ls.degree = degree(c);
        ls.value = [&c](const Vec3& x) { return phi(c, x); };
        const Vec3 origin = {0, 0, 0}; // faces in absolute coordinates
        const std::array<double, 2> totals = exact_totals(c);
        for (const std::string& mesh : meshes)
            for (int q : qs)
                for (const std::string& gen : generators)
                {
                    const bool certify = gen.rfind("certify", 0) == 0;
                    if (gen == "quadgen" && (mesh != "hex" || c.shape != Shape::sphere || c.scale != 1.0))
                        continue; // algoim's 2015 engine here takes the sphere itself, on boxes
                    CertifyOptions copt;
                    if (certify && gen.find(':') != std::string::npos)
                        copt.margin = std::stod(gen.substr(gen.find(':') + 1));
                    const GeneratorOptions opt = gen == "quadgen" || certify ? GeneratorOptions{} : generator_preset(gen);
                    for (const std::string& part : parts)
                    {
                        SelectionExpr expr = parse_selection_expr(part);
                        compile_selection_expr(expr, {"phi"});
                        const SelectionTerm& term = expr.terms.front();
                        const PartKind kind = part_kind(term);
                        const bool surface = kind == PartKind::interface;
                        if (kind != PartKind::negative && kind != PartKind::interface)
                            throw std::runtime_error("robustness study: parts phi < 0 and phi = 0 only");
                        const double cell_scale = surface ? h * h : h * h * h;
                        auto generate = [&](const TestCell& cell, Rule& rule, CertifyStats& cs)
                        {
                            GeneratorStats stats;
                            if (gen == "quadgen")
                                algoim_quadgen_sphere(cell.box, c.centre, c.radius, term, q, rule);
                            else if (certify)
                                certified_bisection(cell.box, ls, term, q, copt, rule, cs);
                            else
                                algoim_clipped_box(cell.box, ls, term, q, opt, rule, stats);
                        };
                        {
                            // warm up (first calls set up tables), on a cell near the sphere or plane
                            Rule rule;
                            CertifyStats cs;
                            try
                            {
                                const Vec3 p = c.shape == Shape::sphere ? Vec3{c.centre[0] + c.radius, c.centre[1], c.centre[2]}
                                                                         : Vec3{0, 0, 0};
                                generate(grid_cells(mesh, {p[0] - 0.5 * h, p[1] - 0.5 * h, p[2] - 0.5 * h}, h, origin)[0], rule, cs);
                            }
                            catch (const std::exception&)
                            {
                            }
                        }
                        Tally t;
                        for (int i0 = 0; i0 < n_cells; ++i0)
                            for (int i1 = 0; i1 < n_cells; ++i1)
                                for (int i2 = 0; i2 < n_cells; ++i2)
                                {
                                    const Vec3 lo = {-1 + h * i0, -1 + h * i1, -1 + h * i2};
                                    for (const TestCell& cell : grid_cells(mesh, lo, h, origin))
                                    {
                                        ++t.cells;
                                        const Reference ref = reference(c, cell);
                                        const bool on_face = surface && ref.face_area > 0;
                                        const double exact = surface ? ref.area : ref.volume;
                                        bool cut = false;
                                        if (ref.known)
                                        {
                                            t.exact_total += surface ? ref.area + 0.5 * ref.face_area : ref.volume;
                                            cut = surface ? (ref.area > 0 || on_face) : (ref.volume > 0 && ref.volume < cell.volume);
                                        }
                                        else
                                        {
                                            // a sign change at 5^3 points of the cell's box, inside the cell
                                            bool pos = false, neg = false;
                                            for (int a = 0; a < 125; ++a)
                                            {
                                                const Vec3 u = {(a % 5) / 4.0, (a / 5 % 5) / 4.0, (a / 25) / 4.0};
                                                if (!inside_clips(cell.box, u))
                                                    continue;
                                                const double f = phi(c, physical_point(cell.box, u));
                                                pos |= f > 0;
                                                neg |= f < 0;
                                            }
                                            cut = pos && neg;
                                        }
                                        t.cut += cut;

                                        Rule rule;
                                        CertifyStats cs;
                                        bool failed = false;
                                        const auto t0 = std::chrono::steady_clock::now();
                                        try
                                        {
                                            generate(cell, rule, cs);
                                        }
                                        catch (const std::exception&)
                                        {
                                            failed = true;
                                        }
                                        const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                                        t.seconds += dt;
                                        t.max_seconds = std::max(t.max_seconds, dt);
                                        t.bisections += cs.bisections;
                                        t.uncertified += cs.uncertified;
                                        if (cut)
                                            t.points += rule.n_points();

                                        double value = 0;
                                        const bool tet = !cell.box.clips.empty();
                                        for (int p = 0; p < rule.n_points() && !failed; ++p)
                                        {
                                            const Vec3 xi = {rule.points[3 * p], rule.points[3 * p + 1], rule.points[3 * p + 2]};
                                            const double w = rule.weights[p];
                                            if (!std::isfinite(w) || !std::isfinite(xi[0]) || !std::isfinite(xi[1])
                                                || !std::isfinite(xi[2]))
                                            {
                                                failed = true;
                                                break;
                                            }
                                            value += w;
                                            t.negative += w < 0;
                                            const double tol = 1e-12;
                                            bool inside = true;
                                            for (int d = 0; d < 3; ++d)
                                                inside &= xi[d] >= -tol && (tet || xi[d] <= 1 + tol);
                                            if (tet)
                                                inside &= xi[0] + xi[1] + xi[2] <= 1 + tol;
                                            t.outside += !inside;
                                            // reference coordinates of these test cells are box coordinates
                                            const Vec3 x = physical_point(cell.box, xi);
                                            const double f = phi(c, x), g = gradient_norm(c, x);
                                            const double reach = 1e-9 * h * g; // |phi| within 1e-9 h of the interface
                                            t.side += surface ? !(std::abs(f) <= reach) : !(f <= reach);
                                        }
                                        if (failed)
                                        {
                                            ++t.fail;
                                            continue;
                                        }
                                        t.total += value;
                                        if (on_face || !ref.known)
                                            continue;
                                        const double error = std::abs(value - exact);
                                        t.l1 += error;
                                        t.max_error = std::max(t.max_error, error / cell_scale);
                                        t.bad += error > 1e-4 * cell_scale;
                                    }
                                }
                        const double exact_total = c.shape == Shape::cone || c.shape == Shape::double_root || c.shape == Shape::torus
                                                       ? totals[surface ? 1 : 0]
                                                       : t.exact_total;
                        std::printf("%-22s %-4s %2d %-12s %-8s | %5ld | %4ld %4ld %4ld %5ld | %8.1e %8.1e %5ld | %8.1e | %7.1f | "
                                    "%7.1f %8.0f | %ld/%ld\n",
                                    c.name.c_str(), mesh.c_str(), q, gen.c_str(), part.c_str(), t.cut, t.fail, t.negative,
                                    t.outside, t.side, exact_total > 0 ? t.l1 / exact_total : t.l1,
                                    t.max_error, t.bad,
                                    exact_total > 0 ? std::abs(t.total - exact_total) / exact_total : std::abs(t.total),
                                    t.cut ? double(t.points) / t.cut : 0.0, t.seconds / t.cells * 1e6, t.max_seconds * 1e6,
                                    t.bisections, t.uncertified);
                        std::fflush(stdout);
                    }
                }
    }
    return 0;
}
