// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "generators.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <stdexcept>

#include "quadrature_general.hpp"
#include "quadrature_multipoly.hpp"

namespace cutcells::proto
{

namespace
{
using algoim::real;
using algoim::SparkStack; // needed by the algoim_spark_alloc macro
using Poly3 = algoim::xarray<real, 3>;
using Mask3 = algoim::booluarray<3, ALGOIM_M>;

Vec3 to_vec(const algoim::uvector<real, 3>& u) { return {u(0), u(1), u(2)}; }

algoim::uvector<real, 3> to_uvec(const Vec3& u) { return algoim::uvector<real, 3>(u[0], u[1], u[2]); }

double dot3(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

double norm3(const Vec3& a) { return std::sqrt(dot3(a, a)); }

Vec3 cross3(const Vec3& a, const Vec3& b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

Vec3 gradient(const Poly3& p, const Vec3& u)
{
    const auto g = algoim::bernstein::evalBernsteinPolyGradient(p, to_uvec(u));
    return {g(0), g(1), g(2)};
}

/// Physical surface measure per reference surface measure, divided by |det J|:
/// |J^{-T} g| / |g| for the reference gradient g.
double surface_factor(const Mat3& inv, const Vec3& g)
{
    Vec3 y = {0, 0, 0};
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            y[i] += inv[k][i] * g[k];
    return norm3(y) / norm3(g);
}

void append_point(const ClippedBox& box, const Vec3& u, double w, Rule& rule)
{
    const Vec3 xi = reference_point(box, u);
    rule.points.insert(rule.points.end(), xi.begin(), xi.end());
    rule.weights.push_back(w);
}

/// Roots of t -> p(a + t (b - a)) in [0, 1], by sampling and bisection.
void segment_roots(const Poly3& p, const Vec3& a, const Vec3& b, std::vector<Vec3>& out)
{
    auto at = [&](double t) { return Vec3{a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])}; };
    auto f = [&](double t) { return algoim::bernstein::evalBernsteinPoly(p, to_uvec(at(t))); };
    const int ns = 32;
    double t0 = 0.0, f0 = f(0.0);
    for (int s = 1; s <= ns; ++s)
    {
        const double t1 = static_cast<double>(s) / ns, f1 = f(t1);
        if (f0 == 0.0)
            out.push_back(at(t0));
        else if (f0 * f1 < 0.0)
        {
            double lo = t0, hi = t1, flo = f0;
            for (int it = 0; it < 60; ++it)
            {
                const double mid = 0.5 * (lo + hi), fm = f(mid);
                if ((fm < 0.0) == (flo < 0.0))
                {
                    lo = mid;
                    flo = fm;
                }
                else
                    hi = mid;
            }
            out.push_back(at(0.5 * (lo + hi)));
        }
        t0 = t1;
        f0 = f1;
    }
    if (f0 == 0.0)
        out.push_back(at(1.0));
}

/// Cui et al.'s indicator for the slicing planes u_outer = const: the largest |cos|
/// between a plane's trace and the interface's trace on any face of the clipped
/// region, sampled where the interface crosses edges. Returns 1 if a face trace is
/// not a single segment.
double cui_alpha(const Poly3& p, const Polytope& region, int outer)
{
    std::vector<std::vector<Vec3>> crossings(region.edges.size());
    for (std::size_t e = 0; e < region.edges.size(); ++e)
        segment_roots(p, region.vertices[region.edges[e][0]], region.vertices[region.edges[e][1]], crossings[e]);

    Vec3 axis = {0, 0, 0};
    axis[outer] = 1.0;
    double alpha = 0.0;
    for (std::size_t f = 0; f < region.face_normals.size(); ++f)
    {
        const Vec3& n = region.face_normals[f];
        std::vector<Vec3> points;
        for (int e : region.face_edges[f])
            for (const Vec3& x : crossings[e])
            {
                bool seen = false;
                for (const Vec3& y : points)
                    seen |= std::abs(x[0] - y[0]) + std::abs(x[1] - y[1]) + std::abs(x[2] - y[2]) < 1e-10;
                if (!seen)
                    points.push_back(x);
            }
        if (!points.empty() && points.size() != 2)
            return 1.0;
        const Vec3 d = cross3(axis, n);
        const double dn = norm3(d);
        if (dn < 1e-12)
            continue; // face parallel to the slicing planes
        for (const Vec3& x : points)
        {
            const Vec3 t = cross3(n, gradient(p, x));
            const double tn = norm3(t);
            if (tn > 0.0)
                alpha = std::max(alpha, std::abs(dot3(d, t)) / (dn * tn));
        }
    }
    return alpha;
}

void integrate_box(const ClippedBox& box, const LevelSet& ls, PartKind kind, int q, const GeneratorOptions& opt,
                   int depth_left, Rule& rule, GeneratorStats& stats)
{
    Poly3 p1(nullptr, algoim::uvector<int, 3>(ls.degree + 1));
    algoim_spark_alloc(real, p1);
    algoim::bernstein::bernsteinInterpolate<3>(
        [&](const algoim::uvector<real, 3>& u) { return ls.value(physical_point(box, to_vec(u))); }, p1);

    bool pos = false, neg = false;
    for (int i = 0; i < p1.size(); ++i)
    {
        pos |= p1[i] > 0.0;
        neg |= p1[i] < 0.0;
    }
    const bool cut = pos && neg;
    if (!cut)
    {
        // uniformly signed box: all of the clipped region, or nothing
        const bool keep = kind == PartKind::whole || (kind == PartKind::negative && !pos)
                          || (kind == PartKind::positive && !neg);
        if (kind == PartKind::interface || !keep)
            return;
    }

    const int nc = static_cast<int>(box.clips.size());
    std::vector<real> clip_coeffs(8 * std::max(nc, 1));
    std::deque<Poly3> clip_polys;
    for (int j = 0; j < nc; ++j)
    {
        clip_polys.emplace_back(clip_coeffs.data() + 8 * j, algoim::uvector<int, 3>(2));
        const HalfSpace h = box.clips[j];
        algoim::bernstein::bernsteinInterpolate<3>(
            [&](const algoim::uvector<real, 3>& u) { return h.c[0] * u(0) + h.c[1] * u(1) + h.c[2] * u(2) - h.d; },
            clip_polys.back());
    }
    std::vector<const Poly3*> polys;
    if (cut)
        polys.push_back(&p1);
    for (const Poly3& c : clip_polys)
        polys.push_back(&c);

    Mask3 mask(true);
    if (opt.cell_masks && nc > 0)
        for (algoim::MultiLoop<3> i(0, ALGOIM_M); ~i; ++i)
        {
            Vec3 lo, hi;
            for (int d = 0; d < 3; ++d)
            {
                lo[d] = static_cast<double>(i(d)) / ALGOIM_M;
                hi[d] = static_cast<double>(i(d) + 1) / ALGOIM_M;
            }
            mask(i()) = may_meet_clips(box, lo, hi);
        }
    const std::vector<Mask3> masks(polys.size(), mask);

    const bool use_alpha = opt.axes == AxisChoice::alpha || (opt.split_alpha <= 1.0 && depth_left > 0);
    if (use_alpha && cut)
    {
        const Polytope region = clipped_polytope(box);
        if (region.n_vertices() < 4)
            return; // the clips leave nothing of this box
        Vec3 centre = {0, 0, 0};
        for (const Vec3& v : region.vertices)
            for (int d = 0; d < 3; ++d)
                centre[d] += v[d] / region.n_vertices();
        const Vec3 g = gradient(p1, centre);
        int k3 = 0;
        for (int d = 1; d < 3; ++d)
            if (std::abs(g[d]) > std::abs(g[k3]))
                k3 = d;
        const int a = k3 == 0 ? 1 : 0, b = k3 == 2 ? 1 : 2;
        const double alpha_a = cui_alpha(p1, region, a), alpha_b = cui_alpha(p1, region, b);
        if (std::min(alpha_a, alpha_b) >= opt.split_alpha && depth_left > 0)
        {
            ++stats.splits;
            const int axis = longest_axis(box);
            for (int half = 0; half < 2; ++half)
            {
                Vec3 lo = {0, 0, 0}, hi = {1, 1, 1};
                (half == 0 ? hi : lo)[axis] = 0.5;
                if (may_meet_clips(box, lo, hi))
                    integrate_box(sub_box(box, lo, hi), ls, kind, q, opt, depth_left - 1, rule, stats);
            }
            return;
        }
        if (opt.axes == AxisChoice::alpha)
        {
            const int outer = alpha_a <= alpha_b ? a : b;
            const int middle = outer == a ? b : a;
            algoim::force_k[3] = k3;
            algoim::force_k[2] = middle < k3 ? middle : middle - 1;
        }
    }
    ++stats.boxes;

    const double detj = std::abs(jacobian_determinant(box));
    if (polys.empty())
    {
        // uniformly signed box without clips: tensor-product Gauss-Legendre
        for (algoim::MultiLoop<3> i(0, q); ~i; ++i)
        {
            Vec3 u;
            double w = detj;
            for (int d = 0; d < 3; ++d)
            {
                u[d] = algoim::GaussQuad::x(q, i(d));
                w *= algoim::GaussQuad::w(q, i(d));
            }
            append_point(box, u, w, rule);
        }
        return;
    }

    algoim::ImplicitPolyQuadrature<3> ipq(static_cast<int>(polys.size()), polys.data(), masks.data());
    algoim::force_k[3] = algoim::force_k[2] = -1;
    const algoim::QuadStrategy strategy = opt.gauss_legendre ? algoim::AlwaysGL : algoim::AutoMixed;

    if (kind == PartKind::interface)
    {
        const Mat3 inv = inverse_jacobian(box);
        ipq.integrate_surf(strategy, q,
                           [&](const algoim::uvector<real, 3>& x, real w, const algoim::uvector<real, 3>&)
                           {
                               // keep points of the level set, not of a clip plane, inside the clips
                               const Vec3 u = to_vec(x);
                               const Vec3 g = gradient(p1, u);
                               const double d1 = std::abs(algoim::bernstein::evalBernsteinPoly(p1, x)) / norm3(g);
                               for (const HalfSpace& h : box.clips)
                                   if (std::abs(dot3(h.c, u) - h.d) / norm3(h.c) < d1)
                                       return;
                               if (inside_clips(box, u))
                                   append_point(box, u, w * detj * surface_factor(inv, g), rule);
                           });
        return;
    }
    ipq.integrate(strategy, q,
                  [&](const algoim::uvector<real, 3>& x, real w)
                  {
                      const Vec3 u = to_vec(x);
                      if (!inside_clips(box, u))
                          return;
                      if (cut)
                      {
                          const double v = algoim::bernstein::evalBernsteinPoly(p1, x);
                          if ((kind == PartKind::negative && !(v < 0.0)) || (kind == PartKind::positive && !(v > 0.0)))
                              return;
                      }
                      append_point(box, u, w * detj, rule);
                  });
}

/// |x(u) - centre|^2 - radius^2 for algoim's 2015 engine, which needs values and
/// gradients in interval arithmetic.
struct SphereInBox
{
    Vec3 offset; ///< origin - centre
    Mat3 jac;
    double radius = 0;
    double sign = 1;

    template <typename T>
    T operator()(const algoim::uvector<T, 3>& u) const
    {
        T s = T(0.0);
        for (int i = 0; i < 3; ++i)
        {
            T xi = T(offset[i]);
            for (int k = 0; k < 3; ++k)
                xi = xi + jac[i][k] * u(k);
            s = s + xi * xi;
        }
        return sign * (s - radius * radius);
    }

    template <typename T>
    algoim::uvector<T, 3> grad(const algoim::uvector<T, 3>& u) const
    {
        algoim::uvector<T, 3> g;
        for (int k = 0; k < 3; ++k)
        {
            T acc = T(0.0);
            for (int i = 0; i < 3; ++i)
            {
                T xi = T(offset[i]);
                for (int m = 0; m < 3; ++m)
                    xi = xi + jac[i][m] * u(m);
                acc = acc + (2.0 * jac[i][k]) * xi;
            }
            g(k) = sign * acc;
        }
        return g;
    }
};
} // namespace

PartKind part_kind(const SelectionTerm& term)
{
    const std::uint64_t others = ~std::uint64_t(1);
    if ((term.zero_required | term.negative_required | term.positive_required) & others)
        throw std::runtime_error("prototype: only selections on the first level set are supported");
    if (term.zero_required & 1)
        return PartKind::interface;
    if (term.negative_required & 1)
        return PartKind::negative;
    if (term.positive_required & 1)
        return PartKind::positive;
    return PartKind::whole;
}

void algoim_clipped_box(const ClippedBox& cell, const LevelSet& ls, const SelectionTerm& term, int q,
                        const GeneratorOptions& opt, Rule& rule, GeneratorStats& stats)
{
    integrate_box(cell, ls, part_kind(term), q, opt, opt.split_depth, rule, stats);
}

void algoim_quadgen_sphere(const ClippedBox& cell, const Vec3& centre, double radius, const SelectionTerm& term, int q,
                           Rule& rule)
{
    if (!cell.clips.empty())
        throw std::runtime_error("algoim_quadgen_sphere: clipped boxes are not supported");
    const PartKind kind = part_kind(term);
    SphereInBox phi;
    for (int i = 0; i < 3; ++i)
        phi.offset[i] = cell.origin[i] - centre[i];
    phi.jac = cell.jacobian;
    phi.radius = radius;
    phi.sign = kind == PartKind::positive ? -1.0 : 1.0;

    const double detj = std::abs(jacobian_determinant(cell));
    const algoim::HyperRectangle<real, 3> unit(algoim::uvector<real, 3>(0.0), algoim::uvector<real, 3>(1.0));
    if (kind == PartKind::whole)
    {
        for (algoim::MultiLoop<3> i(0, q); ~i; ++i)
        {
            Vec3 u;
            double w = detj;
            for (int d = 0; d < 3; ++d)
            {
                u[d] = algoim::GaussQuad::x(q, i(d));
                w *= algoim::GaussQuad::w(q, i(d));
            }
            append_point(cell, u, w, rule);
        }
        return;
    }
    const bool surface = kind == PartKind::interface;
    const auto qr = algoim::quadGen<3>(phi, unit, surface ? 3 : -1, -1, q);
    const Mat3 inv = inverse_jacobian(cell);
    for (const auto& node : qr.nodes)
    {
        const Vec3 u = to_vec(node.x);
        double w = node.w * detj;
        if (surface)
        {
            const auto g = phi.grad<real>(node.x);
            w *= surface_factor(inv, {g(0), g(1), g(2)});
        }
        append_point(cell, u, w, rule);
    }
}

GeneratorOptions generator_preset(const std::string& name)
{
    GeneratorOptions opt;
    opt.name = name;
    if (name == "algoim-auto")
        opt.gauss_legendre = false;
    else if (name == "algoim-gl")
    {
    }
    else if (name == "gl-cellmask")
        opt.cell_masks = true;
    else if (name == "alpha")
    {
        opt.cell_masks = true;
        opt.axes = AxisChoice::alpha;
    }
    else if (name == "split")
    {
        // algoim's axes; alpha only decides where to split
        opt.cell_masks = true;
        opt.split_alpha = 0.99;
        opt.split_depth = 3;
    }
    else if (name == "alpha-split")
    {
        opt.cell_masks = true;
        opt.axes = AxisChoice::alpha;
        opt.split_alpha = 0.99;
        opt.split_depth = 3;
    }
    else
        throw std::runtime_error("unknown generator: " + name);
    return opt;
}

} // namespace cutcells::proto

// Certify-and-bisect engine: kept in this translation unit (see certified.inl).
#include "certified.inl"
