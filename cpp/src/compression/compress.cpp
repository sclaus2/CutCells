// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "compress.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace cutcells::compression
{

namespace
{

/// Exponents of the basis functions, dim per function.
std::vector<int> exponents(MomentSpace space, int dim, int degree)
{
    std::vector<int> e;
    std::array<int, 3> a = {0, 0, 0};
    const int n = degree + 1;
    const int count = dim == 1 ? n : (dim == 2 ? n * n : n * n * n);
    for (int i = 0; i < count; ++i)
    {
        int r = i, sum = 0;
        for (int d = 0; d < dim; ++d)
        {
            a[d] = r % n;
            r /= n;
            sum += a[d];
        }
        if (space == MomentSpace::total && sum > degree)
            continue;
        e.insert(e.end(), a.begin(), a.begin() + dim);
    }
    return e;
}

/// Moment matrix, column-major: column i holds the basis at point i. The basis
/// is products of Legendre polynomials, orthonormal on [-1, 1], of the point's
/// coordinates scaled from the points' bounding box to [-1, 1].
template <std::floating_point T>
std::vector<double> moment_matrix(std::span<const T> points, int dim, int degree, const std::vector<int>& e)
{
    const int n_points = static_cast<int>(points.size()) / dim;
    const int m = static_cast<int>(e.size()) / dim;
    std::array<double, 3> lo, hi;
    lo.fill(std::numeric_limits<double>::infinity());
    hi.fill(-std::numeric_limits<double>::infinity());
    for (int i = 0; i < n_points; ++i)
        for (int d = 0; d < dim; ++d)
        {
            lo[d] = std::min(lo[d], static_cast<double>(points[i * dim + d]));
            hi[d] = std::max(hi[d], static_cast<double>(points[i * dim + d]));
        }
    std::vector<double> a(static_cast<std::size_t>(m) * n_points);
    std::vector<double> leg(static_cast<std::size_t>(dim) * (degree + 1));
    for (int i = 0; i < n_points; ++i)
    {
        for (int d = 0; d < dim; ++d)
        {
            const double h = hi[d] - lo[d];
            const double y = h > 0 ? 2 * (points[i * dim + d] - lo[d]) / h - 1 : 0;
            double* p = leg.data() + d * (degree + 1);
            p[0] = 1;
            if (degree > 0)
                p[1] = y;
            for (int k = 1; k < degree; ++k)
                p[k + 1] = ((2 * k + 1) * y * p[k] - k * p[k - 1]) / (k + 1);
            for (int k = 0; k <= degree; ++k)
                p[k] *= std::sqrt(2.0 * k + 1);
        }
        double* col = a.data() + static_cast<std::size_t>(i) * m;
        for (int j = 0; j < m; ++j)
        {
            double v = 1;
            for (int d = 0; d < dim; ++d)
                v *= leg[d * (degree + 1) + e[j * dim + d]];
            col[j] = v;
        }
    }
    return a;
}

/// Caratheodory reduction: x >= 0 on the k columns of g (m x k, row-major)
/// becomes x' >= 0 with g x' = g x and at most rank(g) nonzeros.
///
/// A Householder QR of g with column pivoting, g P = Q [R11 R12; 0 ~0], gives
/// the null space as P [-R11^-1 R12; I]: null vector l is nonzero only on the
/// r pivot columns and trailing column r + l. Each null vector in turn moves x
/// until one entry reaches 0, and the later null vectors are made to vanish at
/// that entry, which keeps their support within the pivots and the trailing
/// columns up to their own. The basis holds the constant function, so every
/// null vector has entries of both signs. Row-major storage keeps every inner
/// loop a contiguous update across columns.
void caratheodory(std::vector<double>& g, int m, int k, std::vector<double>& x)
{
    auto row = [&](int i) { return g.data() + static_cast<std::size_t>(i) * k; };
    std::vector<int> perm(k);
    std::iota(perm.begin(), perm.end(), 0);
    std::vector<double> sq(k, 0.0), sq_ref(k), v(m), d(k);
    for (int i = 0; i < m; ++i)
    {
        const double* gi = row(i);
        for (int j = 0; j < k; ++j)
            sq[j] += gi[j] * gi[j];
    }
    sq_ref = sq;
    auto swap_columns = [&](int a, int b)
    {
        for (int i = 0; i < m; ++i)
            std::swap(row(i)[a], row(i)[b]);
        std::swap(perm[a], perm[b]);
        std::swap(sq[a], sq[b]);
        std::swap(sq_ref[a], sq_ref[b]);
        std::swap(x[a], x[b]);
    };

    const int steps = std::min(k, m);
    const double rank_tol = 1e3 * std::numeric_limits<double>::epsilon() * std::max(k, m);
    double r00 = 0;
    int r = 0;
    for (int s = 0; s < steps; ++s)
    {
        int best = s;
        for (int j = s + 1; j < k; ++j)
            if (sq[j] > sq[best])
                best = j;
        if (best != s)
            swap_columns(s, best);
        double alpha = 0;
        for (int i = s; i < m; ++i)
            alpha += row(i)[s] * row(i)[s];
        alpha = std::sqrt(alpha);
        if (s == 0)
            r00 = alpha;
        if (!(alpha > rank_tol * r00))
            break;
        // reflector H = I - t v v^T with v = (1, g[s+1..m-1][s] / v0)
        const double beta = row(s)[s] >= 0 ? -alpha : alpha;
        const double v0 = row(s)[s] - beta;
        v[s] = 1;
        for (int i = s + 1; i < m; ++i)
            v[i] = row(i)[s] / v0;
        const double t = -v0 / beta;
        // d = t v^T G[s.., s+1..], then G[s.., s+1..] -= v d
        const int nc = k - s - 1;
        double* dd = d.data() + s + 1;
        std::copy_n(row(s) + s + 1, nc, dd);
        for (int i = s + 1; i < m; ++i)
        {
            const double vi = v[i];
            const double* gi = row(i) + s + 1;
            for (int j = 0; j < nc; ++j)
                dd[j] += vi * gi[j];
        }
        for (int j = 0; j < nc; ++j)
            dd[j] *= t;
        for (int i = s; i < m; ++i)
        {
            const double vi = v[i];
            double* gi = row(i) + s + 1;
            for (int j = 0; j < nc; ++j)
                gi[j] -= vi * dd[j];
        }
        row(s)[s] = beta;
        for (int i = s + 1; i < m; ++i)
            row(i)[s] = 0;
        // downdate the remaining column norms; recompute once cancellation sets in
        const double* gs = row(s);
        bool recompute = false;
        for (int j = s + 1; j < k; ++j)
        {
            sq[j] -= gs[j] * gs[j];
            recompute |= sq[j] < 1e-6 * sq_ref[j];
        }
        if (recompute)
        {
            std::fill(sq.begin() + s + 1, sq.end(), 0.0);
            for (int i = s + 1; i < m; ++i)
            {
                const double* gi = row(i);
                for (int j = s + 1; j < k; ++j)
                    sq[j] += gi[j] * gi[j];
            }
            std::copy(sq.begin() + s + 1, sq.end(), sq_ref.begin() + s + 1);
        }
        ++r;
    }

    // null vectors as the columns of nv (k x n_null, row-major):
    // rows 0..r-1 hold -R11^-1 R12, rows r..k-1 the identity
    const int n_null = k - r;
    std::vector<double> nv(static_cast<std::size_t>(k) * n_null, 0.0);
    auto nrow = [&](int i) { return nv.data() + static_cast<std::size_t>(i) * n_null; };
    for (int i = r - 1; i >= 0; --i)
    {
        double* y = nrow(i);
        const double* gi = row(i);
        for (int l = 0; l < n_null; ++l)
            y[l] = gi[r + l];
        for (int j = i + 1; j < r; ++j)
        {
            const double rij = gi[j];
            const double* yj = nrow(j);
            for (int l = 0; l < n_null; ++l)
                y[l] -= rij * yj[l];
        }
        const double inv = 1.0 / gi[i];
        for (int l = 0; l < n_null; ++l)
            y[l] *= inv;
    }
    for (int i = 0; i < r; ++i)
    {
        double* y = nrow(i);
        for (int l = 0; l < n_null; ++l)
            y[l] = -y[l];
    }
    for (int l = 0; l < n_null; ++l)
        nrow(r + l)[l] = 1;

    std::vector<char> alive(k, 1);
    for (int i = 0; i < k; ++i)
        if (!(x[i] > 0))
        {
            x[i] = 0;
            alive[i] = 0;
        }
    std::vector<double> f(n_null);
    for (int l = 0; l < n_null; ++l)
    {
        const int end = r + l + 1; // support of null vector l
        double vmax = 0;
        for (int i = 0; i < end; ++i)
            vmax = std::max(vmax, std::abs(nrow(i)[l]));
        if (!(vmax > 0))
            continue;
        // the entry that reaches 0 first along +v, or along -v
        int at = -1;
        double step = std::numeric_limits<double>::infinity(), sign = 1;
        for (const double sg : {1.0, -1.0})
        {
            for (int i = 0; i < end; ++i)
            {
                const double vi = sg * nrow(i)[l];
                if (alive[i] && vi > 1e-13 * vmax && x[i] / vi < step)
                {
                    step = x[i] / vi;
                    at = i;
                }
            }
            if (at >= 0)
            {
                sign = sg;
                break;
            }
        }
        if (at < 0)
            continue;
        for (int i = 0; i < end; ++i)
            if (alive[i])
                x[i] = std::max(0.0, x[i] - step * sign * nrow(i)[l]);
        x[at] = 0;
        alive[at] = 0;
        // later null vectors: w_l2 -= (w_l2[at] / v[at]) v
        const int nl = n_null - l - 1;
        if (nl == 0)
            continue;
        const double* ra = nrow(at) + l + 1;
        const double va = nrow(at)[l];
        for (int l2 = 0; l2 < nl; ++l2)
            f[l2] = ra[l2] / va;
        for (int i = 0; i < end; ++i)
        {
            double* ni = nrow(i);
            const double vi = ni[l];
            if (vi == 0)
                continue;
            double* w = ni + l + 1;
            for (int l2 = 0; l2 < nl; ++l2)
                w[l2] -= vi * f[l2];
        }
        std::fill_n(nrow(at) + l + 1, nl, 0.0);
    }

    // back to the caller's column order
    std::vector<double> xs(k);
    for (int i = 0; i < k; ++i)
        xs[perm[i]] = x[i];
    x = std::move(xs);
}

/// Least-squares weights on the columns `cols` of a (m rows) for the moments
/// mom, by Householder QR; false if the columns are numerically dependent.
bool least_squares(const std::vector<double>& a, int m, const std::vector<int>& cols,
                   const std::vector<double>& mom, std::vector<double>& x)
{
    const int n = static_cast<int>(cols.size());
    if (n > m)
        return false;
    std::vector<double> r(static_cast<std::size_t>(n) * m), rhs = mom;
    for (int j = 0; j < n; ++j)
        std::copy_n(a.begin() + static_cast<std::ptrdiff_t>(cols[j]) * m, m,
                    r.begin() + static_cast<std::ptrdiff_t>(j) * m);
    double r00 = 0;
    for (int s = 0; s < n; ++s)
    {
        double* c = r.data() + static_cast<std::size_t>(s) * m;
        double alpha = 0;
        for (int i = s; i < m; ++i)
            alpha += c[i] * c[i];
        alpha = std::sqrt(alpha);
        if (s == 0)
            r00 = alpha;
        if (!(alpha > 1e-10 * r00))
            return false;
        const double beta = c[s] >= 0 ? -alpha : alpha;
        const double v0 = c[s] - beta;
        for (int i = s + 1; i < m; ++i)
            c[i] /= v0;
        c[s] = beta;
        const double t = -v0 / beta;
        auto reflect = [&](double* y)
        {
            double d = y[s];
            for (int i = s + 1; i < m; ++i)
                d += c[i] * y[i];
            d *= t;
            y[s] -= d;
            for (int i = s + 1; i < m; ++i)
                y[i] -= d * c[i];
        };
        for (int j = s + 1; j < n; ++j)
            reflect(r.data() + static_cast<std::size_t>(j) * m);
        reflect(rhs.data());
    }
    x.assign(n, 0.0);
    for (int s = n - 1; s >= 0; --s)
    {
        double v = rhs[s];
        for (int j = s + 1; j < n; ++j)
            v -= r[static_cast<std::size_t>(j) * m + s] * x[j];
        x[s] = v / r[static_cast<std::size_t>(s) * m + s];
    }
    return true;
}

/// max_j |(a w)_j - mom_j| over the columns `cols` with weights w.
double moment_error(const std::vector<double>& a, int m, const std::vector<int>& cols,
                    const std::vector<double>& w, const std::vector<double>& mom)
{
    std::vector<double> s(m, 0.0);
    for (std::size_t c = 0; c < cols.size(); ++c)
    {
        const double* col = a.data() + static_cast<std::size_t>(cols[c]) * m;
        for (int j = 0; j < m; ++j)
            s[j] += w[c] * col[j];
    }
    double e = 0;
    for (int j = 0; j < m; ++j)
        e = std::max(e, std::abs(s[j] - mom[j]));
    return e;
}

} // namespace

std::string moment_space_to_str(MomentSpace space)
{
    return space == MomentSpace::tensor ? "tensor" : "total";
}

MomentSpace string_to_moment_space(const std::string& name)
{
    if (name == "tensor")
        return MomentSpace::tensor;
    if (name == "total")
        return MomentSpace::total;
    throw std::invalid_argument("compression: the moment space is 'tensor' or 'total', not '" + name + "'");
}

int n_moments(MomentSpace space, int dim, int degree)
{
    if (dim < 1 || dim > 3 || degree < 0)
        throw std::invalid_argument("compression: dim must be 1, 2 or 3 and degree non-negative");
    int n = 1;
    if (space == MomentSpace::tensor)
    {
        for (int d = 0; d < dim; ++d)
            n *= degree + 1;
        return n;
    }
    for (int d = 1; d <= dim; ++d)
        n = n * (degree + d) / d;
    return n;
}

template <std::floating_point T>
bool compress_rule(std::span<const T> points, std::span<const T> weights, int dim, int degree,
                   MomentSpace space, std::vector<T>& out_points, std::vector<T>& out_weights,
                   double& residual)
{
    const int n_points = static_cast<int>(weights.size());
    if (static_cast<int>(points.size()) != n_points * dim)
        throw std::invalid_argument("compression: points must hold dim coordinates per weight");
    const int m = n_moments(space, dim, degree);
    residual = 0;
    out_points.assign(points.begin(), points.end());
    out_weights.assign(weights.begin(), weights.end());
    if (n_points <= m || std::any_of(weights.begin(), weights.end(), [](T w) { return w < T(0); }))
        return false;

    const std::vector<int> e = exponents(space, dim, degree);
    const std::vector<double> a = moment_matrix(points, dim, degree, e);
    std::vector<double> w(weights.begin(), weights.end());
    double total = 0;
    for (const double v : w)
        total += v;
    if (!(total > 0))
        return false;
    std::vector<double> mom(m, 0.0);
    for (int i = 0; i < n_points; ++i)
        for (int j = 0; j < m; ++j)
            mom[j] += w[i] * a[static_cast<std::size_t>(i) * m + j];

    std::vector<int> live;
    for (int i = 0; i < n_points; ++i)
        if (w[i] > 0)
            live.push_back(i);

    // recombination: 2m groups of consecutive points, each replaced by its
    // weighted mean column (g is m x groups, row-major), until the points fit
    // one reduction
    const int groups = 2 * m;
    std::vector<double> g, x;
    while (static_cast<int>(live.size()) > groups)
    {
        const int n = static_cast<int>(live.size());
        g.assign(static_cast<std::size_t>(groups) * m, 0.0);
        x.assign(groups, 0.0);
        std::vector<int> first(groups + 1);
        for (int q = 0; q <= groups; ++q)
            first[q] = static_cast<int>(static_cast<std::int64_t>(q) * n / groups);
        std::vector<double> mean(m);
        for (int q = 0; q < groups; ++q)
        {
            std::fill(mean.begin(), mean.end(), 0.0);
            for (int p = first[q]; p < first[q + 1]; ++p)
            {
                const int i = live[p];
                x[q] += w[i];
                const double* ai = a.data() + static_cast<std::size_t>(i) * m;
                for (int j = 0; j < m; ++j)
                    mean[j] += w[i] * ai[j];
            }
            for (int j = 0; j < m; ++j)
                g[static_cast<std::size_t>(j) * groups + q] = mean[j] / x[q];
        }
        const std::vector<double> before = x;
        caratheodory(g, m, groups, x);
        std::vector<int> next;
        for (int q = 0; q < groups; ++q)
        {
            if (!(x[q] > 0))
                continue;
            const double f = x[q] / before[q];
            for (int p = first[q]; p < first[q + 1]; ++p)
            {
                w[live[p]] *= f;
                next.push_back(live[p]);
            }
        }
        if (static_cast<int>(next.size()) >= n) // no progress: reduce directly
            break;
        live = std::move(next);
    }

    // direct reduction on the remaining points
    {
        const int k = static_cast<int>(live.size());
        g.resize(static_cast<std::size_t>(k) * m);
        x.resize(k);
        for (int c = 0; c < k; ++c)
        {
            const double* ac = a.data() + static_cast<std::size_t>(live[c]) * m;
            for (int j = 0; j < m; ++j)
                g[static_cast<std::size_t>(j) * k + c] = ac[j];
            x[c] = w[live[c]];
        }
        caratheodory(g, m, k, x);
        std::vector<int> next;
        std::vector<double> xw;
        for (int c = 0; c < k; ++c)
            if (x[c] > 0)
            {
                next.push_back(live[c]);
                xw.push_back(x[c]);
            }
        live = std::move(next);
        x = std::move(xw);
    }

    // polish: exact weights on the selected points, if they stay positive
    double err = moment_error(a, m, live, x, mom);
    std::vector<double> y;
    if (least_squares(a, m, live, mom, y) && std::all_of(y.begin(), y.end(), [](double v) { return v > 0; }))
    {
        const double err_y = moment_error(a, m, live, y, mom);
        if (err_y <= err)
        {
            x = std::move(y);
            err = err_y;
        }
    }
    residual = err / total;

    const int k = static_cast<int>(live.size());
    out_points.resize(static_cast<std::size_t>(k) * dim);
    out_weights.resize(k);
    for (int c = 0; c < k; ++c)
    {
        for (int d = 0; d < dim; ++d)
            out_points[static_cast<std::size_t>(c) * dim + d] = points[static_cast<std::size_t>(live[c]) * dim + d];
        out_weights[c] = static_cast<T>(x[c]);
    }
    return true;
}

template <std::floating_point T>
void compress_rules(const quadrature::QuadratureRules<T>& rules, int degree, MomentSpace space,
                    quadrature::QuadratureRules<T>& out, CompressionStats& stats)
{
    const int dim = rules._tdim;
    const int n_rules = static_cast<int>(rules._parent_map.size());
    if (static_cast<int>(rules._offset.size()) != n_rules + 1)
        throw std::invalid_argument("compression: the rules need num_rules + 1 offsets");
    n_moments(space, dim, degree); // checks dim and degree

    std::vector<std::vector<T>> pts(n_rules), wts(n_rules);
    std::vector<double> res(n_rules, 0.0);
    std::vector<char> status(n_rules, 0); // 0 copied, 1 compressed, 2 negative weight
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
    for (int r = 0; r < n_rules; ++r)
    {
        const std::size_t b = static_cast<std::size_t>(rules._offset[r]);
        const std::size_t n = static_cast<std::size_t>(rules._offset[r + 1]) - b;
        const std::span<const T> p(rules._points.data() + b * dim, n * dim);
        const std::span<const T> w(rules._weights.data() + b, n);
        if (compress_rule(p, w, dim, degree, space, pts[r], wts[r], res[r]))
            status[r] = 1;
        else if (std::any_of(w.begin(), w.end(), [](T v) { return v < T(0); }))
            status[r] = 2;
    }

    out = quadrature::QuadratureRules<T>{};
    out._tdim = dim;
    out._parent_map = rules._parent_map;
    out._offset.reserve(n_rules + 1);
    out._offset.push_back(0);
    stats.n_rules += n_rules;
    stats.points_before += static_cast<std::int64_t>(rules._weights.size());
    for (int r = 0; r < n_rules; ++r)
    {
        out._points.insert(out._points.end(), pts[r].begin(), pts[r].end());
        out._weights.insert(out._weights.end(), wts[r].begin(), wts[r].end());
        out._offset.push_back(static_cast<std::int32_t>(out._weights.size()));
        stats.n_compressed += status[r] == 1;
        stats.n_skipped += status[r] == 2;
        stats.max_residual = std::max(stats.max_residual, res[r]);
    }
    stats.points_after += static_cast<std::int64_t>(out._weights.size());
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template bool compress_rule<float>(std::span<const float>, std::span<const float>, int, int, MomentSpace,
                                   std::vector<float>&, std::vector<float>&, double&);
template bool compress_rule<double>(std::span<const double>, std::span<const double>, int, int, MomentSpace,
                                    std::vector<double>&, std::vector<double>&, double&);
template void compress_rules<float>(const quadrature::QuadratureRules<float>&, int, MomentSpace,
                                    quadrature::QuadratureRules<float>&, CompressionStats&);
template void compress_rules<double>(const quadrature::QuadratureRules<double>&, int, MomentSpace,
                                     quadrature::QuadratureRules<double>&, CompressionStats&);

} // namespace cutcells::compression
