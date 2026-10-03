// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Compression of quadrature rules: compressed rules keep every monomial
// moment of the space, use at most as many points as the space has moments,
// and have positive weights. Cases: a tensor Gauss rule on the unit cube, a
// ball cut out of a fine Gauss rule (tensor and total spaces), a rule on a
// plane inside the cube (rank deficient: fewer points), random points in a
// triangle (2D, total degree), a batch with a small rule and one with a
// negative weight, and float. Exits non-zero on failure.

#include <cutcells/compression/compress.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using namespace cutcells;
using namespace cutcells::compression;

namespace
{

int failures = 0;

void check(bool ok, const std::string& what)
{
    if (!ok)
    {
        std::printf("FAIL: %s\n", what.c_str());
        ++failures;
    }
}

/// Gauss-Legendre points and weights on [0, 1].
void gauss(int n, std::vector<double>& x, std::vector<double>& w)
{
    x.assign(n, 0.0);
    w.assign(n, 0.0);
    for (int i = 0; i < n; ++i)
    {
        double t = std::cos(M_PI * (i + 0.75) / (n + 0.5));
        for (int it = 0; it < 100; ++it)
        {
            double p0 = 1, p1 = t;
            for (int k = 1; k < n; ++k)
            {
                const double p2 = ((2 * k + 1) * t * p1 - k * p0) / (k + 1);
                p0 = p1;
                p1 = p2;
            }
            const double dp = n * (t * p1 - p0) / (t * t - 1);
            const double dt = p1 / dp;
            t -= dt;
            if (std::abs(dt) < 1e-16)
                break;
        }
        double p0 = 1, p1 = t;
        for (int k = 1; k < n; ++k)
        {
            const double p2 = ((2 * k + 1) * t * p1 - k * p0) / (k + 1);
            p0 = p1;
            p1 = p2;
        }
        const double dp = n * (t * p1 - p0) / (t * t - 1);
        x[i] = 0.5 * (t + 1);
        w[i] = 1.0 / ((1 - t * t) * dp * dp);
    }
}

/// Tensor Gauss rule with n points per direction on [0, 1]^dim.
void tensor_gauss(int dim, int n, std::vector<double>& points, std::vector<double>& weights)
{
    std::vector<double> x, w;
    gauss(n, x, w);
    points.clear();
    weights.clear();
    const int count = dim == 2 ? n * n : n * n * n;
    for (int i = 0; i < count; ++i)
    {
        int r = i;
        double wt = 1;
        for (int d = 0; d < dim; ++d)
        {
            points.push_back(x[r % n]);
            wt *= w[r % n];
            r /= n;
        }
        weights.push_back(wt);
    }
}

/// Largest difference of the monomial moments of the space between two
/// rules, relative to the first rule's total |weight|.
template <typename T>
double moment_difference(const std::vector<T>& p0, const std::vector<T>& w0, const std::vector<T>& p1,
                         const std::vector<T>& w1, int dim, int degree, MomentSpace space)
{
    double scale = 0;
    for (const T v : w0)
        scale += std::abs(static_cast<double>(v));
    double worst = 0;
    const int n = degree + 1;
    const int count = dim == 1 ? n : (dim == 2 ? n * n : n * n * n);
    for (int i = 0; i < count; ++i)
    {
        int e[3] = {0, 0, 0}, r = i, sum = 0;
        for (int d = 0; d < dim; ++d)
        {
            e[d] = r % n;
            r /= n;
            sum += e[d];
        }
        if (space == MomentSpace::total && sum > degree)
            continue;
        auto moment = [&](const std::vector<T>& p, const std::vector<T>& w)
        {
            double s = 0;
            for (std::size_t k = 0; k < w.size(); ++k)
            {
                double v = w[k];
                for (int d = 0; d < dim; ++d)
                    v *= std::pow(static_cast<double>(p[k * dim + d]), e[d]);
                s += v;
            }
            return s;
        };
        worst = std::max(worst, std::abs(moment(p0, w0) - moment(p1, w1)) / scale);
    }
    return worst;
}

/// Compress one rule and check moments, point count and positivity.
template <typename T>
void check_rule(const std::string& name, const std::vector<T>& points, const std::vector<T>& weights, int dim,
                int degree, MomentSpace space, int max_points, double tol)
{
    std::vector<T> p, w;
    double residual = 0;
    const bool compressed = compress_rule<T>(points, weights, dim, degree, space, p, w, residual);
    const double diff = moment_difference(points, weights, p, w, dim, degree, space);
    const bool positive = std::all_of(w.begin(), w.end(), [](T v) { return v > T(0); });
    std::printf("%-34s %5zu -> %4zu points (at most %4d), moment error %.1e, residual %.1e\n", name.c_str(),
                weights.size(), w.size(), max_points, diff, residual);
    check(compressed, name + ": compressed");
    check(static_cast<int>(w.size()) <= max_points, name + ": point count");
    check(diff <= tol, name + ": moments");
    check(positive, name + ": positive weights");
}

} // namespace

int main()
{
    check(n_moments(MomentSpace::tensor, 3, 4) == 125, "n_moments tensor");
    check(n_moments(MomentSpace::total, 3, 4) == 35, "n_moments total 3D");
    check(n_moments(MomentSpace::total, 2, 6) == 28, "n_moments total 2D");
    check(string_to_moment_space(moment_space_to_str(MomentSpace::total)) == MomentSpace::total, "names");

    // tensor Gauss rule on the cube: 512 points, Q_4 needs 125
    std::vector<double> points, weights;
    tensor_gauss(3, 8, points, weights);
    check_rule("cube, Gauss 8^3, Q_4", points, weights, 3, 4, MomentSpace::tensor, 125, 1e-13);
    check_rule("cube, Gauss 8^3, Q_6", points, weights, 3, 6, MomentSpace::tensor, 343, 1e-13);

    // a ball cut out of a fine Gauss rule: an irregular positive rule
    std::vector<double> gp, gw, ball_p, ball_w;
    tensor_gauss(3, 14, gp, gw);
    for (std::size_t k = 0; k < gw.size(); ++k)
    {
        const double dx = gp[3 * k] - 0.2, dy = gp[3 * k + 1] - 0.3, dz = gp[3 * k + 2] - 0.4;
        if (dx * dx + dy * dy + dz * dz < 0.6 * 0.6)
        {
            ball_p.insert(ball_p.end(), gp.begin() + 3 * k, gp.begin() + 3 * k + 3);
            ball_w.push_back(gw[k]);
        }
    }
    check_rule("ball piece, Q_4", ball_p, ball_w, 3, 4, MomentSpace::tensor, 125, 1e-13);
    check_rule("ball piece, P_4", ball_p, ball_w, 3, 4, MomentSpace::total, 35, 1e-13);
    check_rule("ball piece, Q_2", ball_p, ball_w, 3, 2, MomentSpace::tensor, 27, 1e-13);

    // a rule on the plane z = 0.3: only polynomials in x and y are resolved
    std::vector<double> sq_p, sq_w, plane_p;
    tensor_gauss(2, 12, sq_p, sq_w);
    for (std::size_t k = 0; k < sq_w.size(); ++k)
    {
        plane_p.push_back(sq_p[2 * k]);
        plane_p.push_back(sq_p[2 * k + 1]);
        plane_p.push_back(0.3);
    }
    check_rule("plane in the cube, Q_4", plane_p, sq_w, 3, 4, MomentSpace::tensor, 25, 1e-13);

    // random points with random weights in a triangle, total degree 6
    std::mt19937 gen(7);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<double> tri_p, tri_w;
    while (tri_w.size() < 400)
    {
        const double x = u(gen), y = u(gen);
        if (x + y < 1)
        {
            tri_p.push_back(x);
            tri_p.push_back(y);
            tri_w.push_back(0.1 + u(gen));
        }
    }
    check_rule("random points in a triangle, P_6", tri_p, tri_w, 2, 6, MomentSpace::total, 28, 1e-13);

    // float: the same cube rule
    std::vector<float> fp(points.begin(), points.end()), fw(weights.begin(), weights.end());
    check_rule<float>("cube, Gauss 8^3, Q_4, float", fp, fw, 3, 4, MomentSpace::tensor, 125, 1e-5);

    // a batch: the ball piece, a small rule kept as it is, a negative weight
    quadrature::QuadratureRules<double> rules;
    rules._tdim = 3;
    rules._offset = {0};
    auto append = [&](const std::vector<double>& p, const std::vector<double>& w, int parent)
    {
        rules._points.insert(rules._points.end(), p.begin(), p.end());
        rules._weights.insert(rules._weights.end(), w.begin(), w.end());
        rules._offset.push_back(static_cast<int>(rules._weights.size()));
        rules._parent_map.push_back(parent);
    };
    append(ball_p, ball_w, 4);
    const std::vector<double> small_p(points.begin(), points.begin() + 3 * 20);
    const std::vector<double> small_w(weights.begin(), weights.begin() + 20);
    append(small_p, small_w, 7);
    std::vector<double> negative_w = weights;
    negative_w[5] = -negative_w[5];
    append(points, negative_w, 9);

    quadrature::QuadratureRules<double> out;
    CompressionStats stats;
    compress_rules(rules, 4, MomentSpace::tensor, out, stats);
    const int n0 = out._offset[1], n1 = out._offset[2] - out._offset[1], n2 = out._offset[3] - out._offset[2];
    std::printf("batch: %d, %d, %d points; %d compressed, %d skipped, residual %.1e\n", n0, n1, n2,
                stats.n_compressed, stats.n_skipped, stats.max_residual);
    check(out._parent_map == rules._parent_map, "batch: parent map");
    check(out._tdim == 3 && static_cast<int>(out._points.size()) == 3 * out._offset[3], "batch: layout");
    check(n0 <= 125 && n1 == 20 && n2 == 512, "batch: point counts");
    check(stats.n_rules == 3 && stats.n_compressed == 1 && stats.n_skipped == 1, "batch: stats");
    check(stats.points_before == static_cast<std::int64_t>(rules._weights.size())
              && stats.points_after == static_cast<std::int64_t>(out._weights.size()),
          "batch: point totals");
    check(stats.max_residual < 1e-13, "batch: residual");

    if (failures > 0)
    {
        std::printf("%d failure(s)\n", failures);
        return 1;
    }
    std::printf("all compression tests passed\n");
    return 0;
}
