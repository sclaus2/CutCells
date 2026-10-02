// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Tensor Bernstein forms on boxes: the conversion of cell level sets (a check
// per degree), affine restriction, subdivision, derivatives, margins and roots.
// Exits non-zero on failure.

#include <cutcells/bernstein.h>
#include <cutcells/quadrays/box_bernstein.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <span>
#include <string>
#include <vector>

using namespace cutcells;
using namespace cutcells::quadrays;

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

/// Random point in the reference cell (inside = true) or in its box.
std::vector<double> random_point(cell::type ct, bool inside, std::mt19937& rng)
{
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const int tdim = cell::get_tdim(ct);
    std::vector<double> x(tdim);
    for (;;)
    {
        double sum = 0;
        for (int d = 0; d < tdim; ++d)
            sum += x[d] = u(rng);
        if (!inside || !bernstein::is_simplex(ct) || sum <= 1.0)
            return x;
    }
}

/// The conversion of a cell's level set reproduces it at random points, inside
/// the cell and, for simplices, in the rest of the box.
void test_conversion(cell::type ct, int max_degree)
{
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    for (int n = 0; n <= max_degree; ++n)
    {
        std::vector<double> c(bernstein::num_polynomials(ct, n));
        for (double& v : c)
            v = u(rng);
        BoxBernstein<double> box;
        cell_bernstein_on_box<double>(ct, n, c, box);
        const double scale = max_abs<double>(box.coeffs);
        double worst = 0;
        for (int k = 0; k < 200; ++k)
        {
            const std::vector<double> x = random_point(ct, k % 2 == 0, rng);
            const double ref = bernstein::evaluate<double>(ct, n, c, x);
            worst = std::max(worst, std::abs(evaluate<double>(box, x) - ref) / scale);
        }
        check(worst < 1e-13, cell::cell_type_to_str(ct) + " conversion, degree " + std::to_string(n)
                                 + ": error " + std::to_string(worst));
    }
}

/// Random tensor form of the given degrees.
BoxBernstein<double> random_form(int dim, std::array<int, 3> degree, std::mt19937& rng)
{
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    BoxBernstein<double> p;
    p.dim = dim;
    p.degree = degree;
    p.coeffs.resize(p.size());
    for (double& v : p.coeffs)
        v = u(rng);
    return p;
}

void test_restriction()
{
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<double> work;
    const BoxBernstein<double> p = random_form(3, {2, 3, 1}, rng);
    for (int m = 1; m <= 3; ++m)
        for (int trial = 0; trial < 20; ++trial)
        {
            std::vector<double> origin(3), matrix(3 * m);
            for (double& v : origin)
                v = u(rng);
            for (double& v : matrix)
                v = (rng() % 3 == 0) ? 0.0 : u(rng); // some zero entries
            BoxBernstein<double> r;
            restrict_affine<double>(p, origin, matrix, m, r, work);
            int expected_degree = 0;
            for (int i = 0; i < 3; ++i)
                expected_degree += matrix[i * m] != 0.0 ? p.degree[i] : 0;
            check(r.degree[0] == expected_degree, "restriction degree");
            double worst = 0, scale = 1e-300;
            for (int k = 0; k < 50; ++k)
            {
                std::vector<double> s(m), x(origin);
                for (double& v : s)
                    v = 0.5 * (u(rng) + 1.0);
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < m; ++j)
                        x[i] += matrix[i * m + j] * s[j];
                const double ref = evaluate<double>(p, x);
                scale = std::max(scale, std::abs(ref));
                worst = std::max(worst, std::abs(evaluate<double>(r, s) - ref));
            }
            check(worst <= 1e-13 * std::max(1.0, scale),
                  "restriction to " + std::to_string(m) + " variables: error " + std::to_string(worst));
        }

    // each path of restrict_affine: along axes (with a constant variable), a line
    // through all variables, and a general map
    const std::vector<std::vector<double>> origins = {{0.2, -0.4, 0.7}, {0.1, 0.9, -0.3}, {0.3, 0.2, 0.1}};
    const std::vector<std::vector<double>> matrices = {{0.5, 0.0, 0.0, 0.0, 0.0, 0.0}, // 3 x 2: s0 -> u0, u1 const, u2 const
                                                       {0.4, -0.7, 1.1},             // 3 x 1: a line through u0, u1, u2
                                                       {0.5, 0.0, 0.0, 0.3, 1.0, -1.0}}; // 3 x 2, mixed
    const std::vector<int> ms = {2, 1, 2};
    for (std::size_t c = 0; c < ms.size(); ++c)
    {
        BoxBernstein<double> r;
        restrict_affine<double>(p, origins[c], matrices[c], ms[c], r, work);
        const std::vector<double> s = {0.37, 0.81};
        std::vector<double> x(origins[c]);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < ms[c]; ++j)
                x[i] += matrices[c][i * ms[c] + j] * s[j];
        check(std::abs(evaluate<double>(r, std::span<const double>(s.data(), ms[c])) - evaluate<double>(p, x)) < 1e-13,
              "restriction path " + std::to_string(c));
    }

    // subdivision equals restriction to the sub-box
    BoxBernstein<double> sub;
    const std::vector<double> lo = {0.1, 0.25, 0.5}, hi = {0.6, 0.5, 0.75};
    subdivide<double>(p, lo, hi, sub, work);
    check(sub.degree == p.degree, "subdivision keeps the degree");
    const std::vector<double> s = {0.3, 0.7, 0.2};
    const std::vector<double> x = {lo[0] + s[0] * (hi[0] - lo[0]), lo[1] + s[1] * (hi[1] - lo[1]),
                                   lo[2] + s[2] * (hi[2] - lo[2])};
    check(std::abs(evaluate<double>(sub, s) - evaluate<double>(p, x)) < 1e-14, "subdivision value");
}

void test_derivatives_and_margins()
{
    std::mt19937 rng(3);
    const BoxBernstein<double> p = random_form(3, {3, 2, 4}, rng);
    const std::vector<double> s = {0.2, 0.9, 0.4};
    std::vector<double> g(3);
    gradient<double>(p, s, g);
    BoxBernstein<double> d;
    for (int k = 0; k < 3; ++k)
    {
        derivative<double>(p, k, d);
        check(std::abs(evaluate<double>(d, s) - g[k]) < 1e-12, "derivative " + std::to_string(k));
    }

    // p = 2 y0 - y1 on the box [0, 1] x [0, 0.5]: d/dy0 = 2, d/dy1 = -1
    BoxBernstein<double> lin;
    lin.dim = 2;
    lin.degree = {1, 1, 0};
    lin.coeffs = {0.0, 2.0, -0.5, 1.5}; // values at the corners (y0, y1) = (0|1, 0|0.5)
    std::vector<double> ratio(2);
    const std::vector<double> lengths = {1.0, 0.5};
    margins<double>(lin, lengths, ratio, d);
    check(std::abs(ratio[0] - 2 / std::sqrt(5.0)) < 1e-15 && std::abs(ratio[1] - 1 / std::sqrt(5.0)) < 1e-15,
          "margins of a linear function");
    check(may_vanish<double>(lin.coeffs) && !may_vanish<double>(std::vector<double>{1.0, 2.0}), "may_vanish");
    check(scaled_norm<double>(std::vector<double>{3e-200, 4e-200}) == 5e-200, "scaled_norm without underflow");
}

/// Bernstein coefficients of prod (t - r_i) on [0, 1].
std::vector<double> product_of_roots(const std::vector<double>& roots)
{
    std::vector<double> c = {1.0};
    for (double r : roots)
    {
        const int n = static_cast<int>(c.size()); // degree n - 1 -> n
        std::vector<double> next(n + 1, 0.0);
        for (int k = 0; k <= n; ++k)
        {
            // (t - r) has Bernstein coefficients (-r, 1 - r)
            if (k < n)
                next[k] += double(n - k) / n * c[k] * -r;
            if (k > 0)
                next[k] += double(k) / n * c[k - 1] * (1.0 - r);
        }
        c = next;
    }
    return c;
}

void test_roots()
{
    std::vector<double> work, roots;
    const std::vector<double> expected = {0.3, 0.31, 0.8};
    const std::vector<double> c = product_of_roots(expected);
    isolate_roots<double>(c, 0.0, 1.0, roots, work);
    std::sort(roots.begin(), roots.end());
    check(roots.size() == 3, "three roots isolated");
    for (std::size_t i = 0; i < std::min(roots.size(), expected.size()); ++i)
        check(std::abs(roots[i] - expected[i]) < 1e-14, "root " + std::to_string(i));

    // the same polynomial on [2, 4]: roots map affinely
    roots.clear();
    isolate_roots<double>(c, 2.0, 4.0, roots, work);
    std::sort(roots.begin(), roots.end());
    check(roots.size() == 3 && std::abs(roots[1] - 2.62) < 1e-13, "roots on [2, 4]");

    // a double root does not change sign: no root reported
    roots.clear();
    isolate_roots<double>(product_of_roots({0.5, 0.5}), 0.0, 1.0, roots, work);
    check(roots.empty(), "double root ignored");

    // a simple root exactly at a split point is reported once
    roots.clear();
    isolate_roots<double>(product_of_roots({0.25, 0.5, 0.75}), 0.0, 1.0, roots, work);
    std::sort(roots.begin(), roots.end());
    check(roots.size() == 3 && roots[1] == 0.5, "root at a split point");

    // bracketing a single root
    const std::vector<double> lin = product_of_roots({0.25});
    check(std::abs(bracketed_root<double>(lin, 0.0, 1.0, 0.0, 1.0, lin.front(), lin.back()) - 0.25) < 1e-15,
          "bracketed root");
}

void test_float()
{
    std::vector<float> c = {0.5f, -0.25f, 0.75f, 1.0f};
    BoxBernstein<float> box;
    cell_bernstein_on_box<float>(cell::type::tetrahedron, 1, c, box);
    const std::vector<float> x = {0.2f, 0.3f, 0.1f};
    const float ref = bernstein::evaluate<float>(cell::type::tetrahedron, 1, c, x);
    check(std::abs(evaluate<float>(box, x) - ref) < 1e-6f, "float conversion");
}

} // namespace

int main()
{
    test_conversion(cell::type::triangle, 12);
    test_conversion(cell::type::tetrahedron, 12);
    test_conversion(cell::type::quadrilateral, 8);
    test_conversion(cell::type::hexahedron, 8);
    test_restriction();
    test_derivatives_and_margins();
    test_roots();
    test_float();
    if (failures == 0)
        std::printf("test_box_bernstein: ok\n");
    return failures == 0 ? 0 : 1;
}
