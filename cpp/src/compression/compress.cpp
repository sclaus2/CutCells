// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "compress.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <limits>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace cutcells::compression
{

namespace
{

/// Spaces up to this many moments try non-negative least squares first
/// (compress_rule): on rules of cut cells it reaches the moments there, and
/// beyond it stalls at a residual of 1e-10 to 1e-8.
constexpr int nnls_moments = 128;

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

/// The basis of the moments: products of Legendre polynomials, orthonormal on
/// [-1, 1], of the coordinates scaled from the points' bounding box to
/// [-1, 1]. Columns are evaluated where they are needed (basis_column), so
/// that no matrix over all points is stored.
template <std::floating_point T>
struct Basis
{
    std::span<const T> points;
    int dim = 0, degree = 0, m = 0;
    std::vector<int> e;
    std::array<double, 3> lo{}, hi{};
    std::vector<double> leg; ///< scratch: Legendre values per coordinate
};

template <std::floating_point T>
Basis<T> make_basis(std::span<const T> points, int dim, int degree, MomentSpace space)
{
    Basis<T> b;
    b.points = points;
    b.dim = dim;
    b.degree = degree;
    b.e = exponents(space, dim, degree);
    b.m = static_cast<int>(b.e.size()) / dim;
    b.lo.fill(std::numeric_limits<double>::infinity());
    b.hi.fill(-std::numeric_limits<double>::infinity());
    for (std::size_t i = 0; i < points.size() / dim; ++i)
        for (int c = 0; c < dim; ++c)
        {
            b.lo[c] = std::min(b.lo[c], static_cast<double>(points[i * dim + c]));
            b.hi[c] = std::max(b.hi[c], static_cast<double>(points[i * dim + c]));
        }
    b.leg.resize(static_cast<std::size_t>(dim) * (degree + 1));
    return b;
}

/// The m basis functions at point i, into col.
template <std::floating_point T>
void basis_column(Basis<T>& b, int i, double* col)
{
    const int dim = b.dim, degree = b.degree;
    for (int d = 0; d < dim; ++d)
    {
        const double h = b.hi[d] - b.lo[d];
        const double y = h > 0 ? 2 * (b.points[static_cast<std::size_t>(i) * dim + d] - b.lo[d]) / h - 1 : 0;
        double* p = b.leg.data() + d * (degree + 1);
        p[0] = 1;
        if (degree > 0)
            p[1] = y;
        for (int k = 1; k < degree; ++k)
            p[k + 1] = ((2 * k + 1) * y * p[k] - k * p[k - 1]) / (k + 1);
        for (int k = 0; k <= degree; ++k)
            p[k] *= std::sqrt(2.0 * k + 1);
    }
    for (int j = 0; j < b.m; ++j)
    {
        double v = 1;
        for (int d = 0; d < dim; ++d)
            v *= b.leg[d * (degree + 1) + b.e[j * dim + d]];
        col[j] = v;
    }
}

/// The columns of the points `cols`, column-major (m x cols.size()).
template <std::floating_point T>
std::vector<double> basis_matrix(Basis<T>& b, const std::vector<int>& cols)
{
    std::vector<double> a(static_cast<std::size_t>(b.m) * cols.size());
    for (std::size_t c = 0; c < cols.size(); ++c)
        basis_column(b, cols[c], a.data() + c * b.m);
    return a;
}

/// Lawson-Hanson non-negative least squares: x >= 0 minimising |a x - b| for
/// the k columns of a (m x k, column-major). The passive columns' QR, with Q
/// stored, is updated as columns enter (one Householder reflector) and leave
/// (Givens rotations), so that the whole solve costs about one QR of an
/// m x m matrix. Columns in the span of the passive ones are skipped (a
/// rank-deficient space, say points on a plane). On return @p set holds the
/// passive columns and @p x their positive weights; the result is the
/// largest entry of b - a x.
double nnls(const std::vector<double>& a, int m, int k, const std::vector<double>& b, std::vector<int>& set,
            std::vector<double>& x)
{
    const auto col = [&](int j) { return a.data() + static_cast<std::size_t>(j) * m; };
    std::vector<double> q(static_cast<std::size_t>(m) * m, 0.0), r(static_cast<std::size_t>(m) * m, 0.0);
    for (int i = 0; i < m; ++i)
        q[static_cast<std::size_t>(i) * m + i] = 1;
    const auto qcol = [&](int j) { return q.data() + static_cast<std::size_t>(j) * m; };
    const auto rcol = [&](int j) { return r.data() + static_cast<std::size_t>(j) * m; };
    std::vector<double> qb = b, v(m), res(m), z;
    std::vector<char> state(k, 1); // 1 free, 0 passive, 2 skipped
    set.clear();
    x.clear();
    double bnorm = 0, anorm = 0;
    for (const double c : b)
        bnorm = std::max(bnorm, std::abs(c));
    for (int j = 0; j < k; ++j)
        for (int i = 0; i < m; ++i)
            anorm = std::max(anorm, std::abs(col(j)[i]));
    const double eps = std::numeric_limits<double>::epsilon();

    // append column t to the passive set; false if it lies in their span
    const auto append = [&](int t) -> bool
    {
        const int p = static_cast<int>(set.size());
        if (p >= m)
            return false;
        const double* at = col(t);
        double anorm_t = 0;
        for (int i = 0; i < m; ++i)
        {
            const double* qi = qcol(i);
            double s = 0;
            for (int l = 0; l < m; ++l)
                s += qi[l] * at[l];
            v[i] = s;
            anorm_t += at[i] * at[i];
        }
        double tail = 0;
        for (int i = p; i < m; ++i)
            tail += v[i] * v[i];
        tail = std::sqrt(tail);
        if (!(tail > 1e3 * eps * std::sqrt(anorm_t) * m))
            return false;
        // Householder reflector H = I - beta u u^T on entries p.. with H v = alpha e_p
        const double alpha = v[p] > 0 ? -tail : tail;
        v[p] -= alpha;
        double uu = 0;
        for (int i = p; i < m; ++i)
            uu += v[i] * v[i];
        const double beta = 2 / uu;
        // Q <- Q H on its columns p..: Q_i -= beta u_i (Q u)
        std::fill(res.begin(), res.end(), 0.0);
        for (int i = p; i < m; ++i)
        {
            const double* qi = qcol(i);
            for (int l = 0; l < m; ++l)
                res[l] += qi[l] * v[i];
        }
        for (int i = p; i < m; ++i)
        {
            double* qi = qcol(i);
            const double f = beta * v[i];
            for (int l = 0; l < m; ++l)
                qi[l] -= f * res[l];
        }
        double s = 0;
        for (int i = p; i < m; ++i)
            s += qb[i] * v[i];
        s *= beta;
        for (int i = p; i < m; ++i)
            qb[i] -= s * v[i];
        // R gets the column (Q^T a_t)[0..p), alpha; columns 0..p of Q are unchanged
        double* rc = rcol(p);
        std::copy_n(v.begin(), p, rc);
        rc[p] = alpha;
        set.push_back(t);
        x.push_back(0);
        return true;
    };

    // remove the passive column at position i: Givens rotations restore R
    const auto remove = [&](int i)
    {
        const int p = static_cast<int>(set.size());
        for (int j = i; j + 1 < p; ++j)
            std::copy_n(rcol(j + 1), m, rcol(j));
        set.erase(set.begin() + i);
        x.erase(x.begin() + i);
        for (int j = i; j + 1 < p; ++j)
        {
            double* rj = rcol(j);
            const double f = rj[j], g = rj[j + 1];
            const double h = std::hypot(f, g);
            if (!(h > 0))
                continue;
            const double c = f / h, s = g / h;
            for (int l = j; l + 1 < p; ++l)
            {
                double* rl = rcol(l);
                const double u = rl[j], w = rl[j + 1];
                rl[j] = c * u + s * w;
                rl[j + 1] = -s * u + c * w;
            }
            const double u = qb[j], w = qb[j + 1];
            qb[j] = c * u + s * w;
            qb[j + 1] = -s * u + c * w;
            double* qj = qcol(j);
            double* qj1 = qcol(j + 1);
            for (int l = 0; l < m; ++l)
            {
                const double a0 = qj[l], a1 = qj1[l];
                qj[l] = c * a0 + s * a1;
                qj1[l] = -s * a0 + c * a1;
            }
        }
        std::fill_n(rcol(p - 1), m, 0.0);
    };

    // least-squares solution on the passive set: R z = (Q^T b)[0..p)
    const auto solve = [&]()
    {
        const int p = static_cast<int>(set.size());
        z.assign(p, 0.0);
        for (int i = p - 1; i >= 0; --i)
        {
            double s = qb[i];
            for (int j = i + 1; j < p; ++j)
                s -= rcol(j)[i] * z[j];
            z[i] = s / rcol(i)[i];
        }
    };

    const auto residual = [&]()
    {
        res = b;
        for (std::size_t c = 0; c < set.size(); ++c)
        {
            const double* ac = col(set[c]);
            for (int i = 0; i < m; ++i)
                res[i] -= x[c] * ac[i];
        }
        double e = 0;
        for (const double c : res)
            e = std::max(e, std::abs(c));
        return e;
    };

    const double stop = 1e2 * eps * bnorm;
    for (int it = 0; it < 3 * k + 3 * m; ++it)
    {
        if (residual() <= stop)
            break;
        // the free column most aligned with the residual
        int t = -1;
        double best = 0;
        for (int j = 0; j < k; ++j)
        {
            if (state[j] != 1)
                continue;
            const double* aj = col(j);
            double w = 0;
            for (int i = 0; i < m; ++i)
                w += aj[i] * res[i];
            if (w > best)
            {
                best = w;
                t = j;
            }
        }
        if (t < 0 || !(best > 1e1 * eps * anorm * bnorm * m))
            break;
        if (!append(t))
        {
            state[t] = 2;
            continue;
        }
        state[t] = 0;
        // inner loop: step towards the least-squares solution until it is positive
        for (int inner = 0; inner < m + 1; ++inner)
        {
            solve();
            bool positive = true;
            for (const double c : z)
                positive &= c > 0;
            if (positive)
            {
                x = z;
                break;
            }
            double step = 2;
            std::size_t at = 0;
            for (std::size_t c = 0; c < z.size(); ++c)
            {
                if (z[c] > 0)
                    continue;
                const double d = x[c] - z[c];
                const double ratio = d > 0 ? x[c] / d : 0.0;
                if (ratio < step)
                {
                    step = ratio;
                    at = c;
                }
            }
            for (std::size_t c = 0; c < z.size(); ++c)
                x[c] += step * (z[c] - x[c]);
            x[at] = 0;
            // the entering column leaving at once would cycle: skip it for good
            for (int c = static_cast<int>(set.size()) - 1; c >= 0; --c)
                if (!(x[c] > 0))
                {
                    state[set[c]] = set[c] == t ? 2 : 1;
                    remove(c);
                }
        }
    }
    return residual();
}

/// Householder QR with column pivoting of g (m x k, row-major), in place: on
/// return rows 0..rank-1 hold R, zero below the diagonal, of the columns in
/// the order perm gives (perm[j] = original column now at j); x, unless
/// empty, is permuted along. Only the first k_free columns are pivoted; the
/// others, a right-hand side say, are carried along. Steps stop where the
/// largest remaining column norm falls to tol times the first; the result is
/// the rank.
///
/// Reflectors reach the trailing matrix in blocks (the delayed update of
/// LAPACK's dlaqps): a step brings up to date only its pivot column and its
/// row, reads the trailing matrix once for the reflector's row of F, and a
/// block ends with A -= V F^T in one pass, so that the matrix, larger than the
/// cache from about 300 moments, is streamed once per step instead of twice.
int pivoted_qr(std::vector<double>& g, int m, int k, int k_free, double tol, std::vector<int>& perm,
               std::vector<double>& x)
{
    constexpr int nb = 32;  // reflectors per block
    constexpr int jb = 512; // columns per tile of the block update
    const auto row = [&](int i) { return g.data() + static_cast<std::size_t>(i) * k; };
    perm.resize(k);
    std::iota(perm.begin(), perm.end(), 0);
    // vb: the block's reflectors (m x nb, row-major, entry t of row s0 + t is 1);
    // ft: their rows of F (nb x k), so that the block's update is A -= vb ft
    std::vector<double> sq(k, 0.0), sq_ref(k), vb(static_cast<std::size_t>(m) * nb),
        ft(static_cast<std::size_t>(nb) * k), z(nb);
    const auto norms = [&](int s)
    {
        std::fill(sq.begin() + s, sq.begin() + k_free, 0.0);
        for (int i = s; i < m; ++i)
        {
            const double* gi = row(i);
            for (int j = s; j < k_free; ++j)
                sq[j] += gi[j] * gi[j];
        }
        std::copy(sq.begin() + s, sq.begin() + k_free, sq_ref.begin() + s);
    };
    norms(0);
    const auto swap_columns = [&](int a, int b)
    {
        for (int i = 0; i < m; ++i)
            std::swap(row(i)[a], row(i)[b]);
        for (int t = 0; t < nb; ++t)
            std::swap(ft[static_cast<std::size_t>(t) * k + a], ft[static_cast<std::size_t>(t) * k + b]);
        std::swap(perm[a], perm[b]);
        std::swap(sq[a], sq[b]);
        std::swap(sq_ref[a], sq_ref[b]);
        if (!x.empty())
            std::swap(x[a], x[b]);
    };

    const int steps = std::min(k_free, m);
    double r00 = 0;
    int s = 0;
    while (s < steps)
    {
        const int s0 = s;
        const int end = std::min(steps, s0 + nb);
        bool recompute = false;
        int jj = 0; // reflectors of this block so far
        for (; s < end && !recompute; ++s, ++jj)
        {
            int best = s;
            for (int j = s + 1; j < k_free; ++j)
                if (sq[j] > sq[best])
                    best = j;
            if (best != s)
                swap_columns(s, best);
            // the pivot column, up to date below row s - 1
            double alpha = 0;
            for (int i = s; i < m; ++i)
            {
                double c = row(i)[s];
                const double* vi = vb.data() + static_cast<std::size_t>(i) * nb;
                for (int t = 0; t < jj; ++t)
                    c -= vi[t] * ft[static_cast<std::size_t>(t) * k + s];
                row(i)[s] = c;
                alpha += c * c;
            }
            alpha = std::sqrt(alpha);
            if (s == 0)
                r00 = alpha;
            if (!(alpha > tol * r00))
                return s;
            // reflector H = I - tau v v^T with v = (1, g[s+1..m-1][s] / v0)
            const double beta = row(s)[s] >= 0 ? -alpha : alpha;
            const double v0 = row(s)[s] - beta;
            const double tau = -v0 / beta;
            vb[static_cast<std::size_t>(s) * nb + jj] = 1;
            for (int i = s + 1; i < m; ++i)
            {
                vb[static_cast<std::size_t>(i) * nb + jj] = row(i)[s] / v0;
                row(i)[s] = 0;
            }
            row(s)[s] = beta;
            // row jj of F: tau (A - V F^T)^T v on the columns after s
            double* f = ft.data() + static_cast<std::size_t>(jj) * k;
            std::fill(f, f + k, 0.0);
            std::fill(z.begin(), z.begin() + jj, 0.0);
            for (int i = s; i < m; ++i)
            {
                const double* vi = vb.data() + static_cast<std::size_t>(i) * nb;
                for (int t = 0; t < jj; ++t)
                    z[t] += vi[t] * vi[jj];
            }
            // four rows at a time, so that f is read and written once per four
            int i = s;
            for (; i + 3 < m; i += 4)
            {
                const double v0 = vb[static_cast<std::size_t>(i) * nb + jj],
                             v1 = vb[static_cast<std::size_t>(i + 1) * nb + jj],
                             v2 = vb[static_cast<std::size_t>(i + 2) * nb + jj],
                             v3 = vb[static_cast<std::size_t>(i + 3) * nb + jj];
                const double *g0 = row(i), *g1 = row(i + 1), *g2 = row(i + 2), *g3 = row(i + 3);
                for (int j = s + 1; j < k; ++j)
                    f[j] += v0 * g0[j] + v1 * g1[j] + v2 * g2[j] + v3 * g3[j];
            }
            for (; i < m; ++i)
            {
                const double vij = vb[static_cast<std::size_t>(i) * nb + jj];
                const double* gi = row(i);
                for (int j = s + 1; j < k; ++j)
                    f[j] += vij * gi[j];
            }
            for (int t = 0; t < jj; ++t)
            {
                const double zt = z[t];
                const double* ftt = ft.data() + static_cast<std::size_t>(t) * k;
                for (int j = s + 1; j < k; ++j)
                    f[j] -= zt * ftt[j];
            }
            for (int j = s + 1; j < k; ++j)
                f[j] *= tau;
            // row s of R: all reflectors of the block so far
            double* gs = row(s);
            for (int t = 0; t <= jj; ++t)
            {
                const double vt = vb[static_cast<std::size_t>(s) * nb + t];
                const double* ftt = ft.data() + static_cast<std::size_t>(t) * k;
                for (int j = s + 1; j < k; ++j)
                    gs[j] -= vt * ftt[j];
            }
            // downdate the remaining column norms; once cancellation sets in,
            // end the block and recompute them
            for (int j = s + 1; j < k_free; ++j)
            {
                sq[j] -= gs[j] * gs[j];
                recompute |= sq[j] < 1e-6 * sq_ref[j];
            }
        }
        // the block's update of the rows and columns after it
        for (int j0 = s; j0 < k; j0 += jb)
        {
            const int j1 = std::min(k, j0 + jb);
            for (int i = s; i < m; ++i)
            {
                double* gi = row(i);
                const double* vi = vb.data() + static_cast<std::size_t>(i) * nb;
                const auto frow = [&](int t) { return ft.data() + static_cast<std::size_t>(t) * k; };
                // four reflectors at a time, so that the row is read and written once per four
                int t = 0;
                for (; t + 3 < jj; t += 4)
                {
                    const double v0 = vi[t], v1 = vi[t + 1], v2 = vi[t + 2], v3 = vi[t + 3];
                    const double *f0 = frow(t), *f1 = frow(t + 1), *f2 = frow(t + 2), *f3 = frow(t + 3);
                    for (int j = j0; j < j1; ++j)
                        gi[j] -= v0 * f0[j] + v1 * f1[j] + v2 * f2[j] + v3 * f3[j];
                }
                for (; t < jj; ++t)
                {
                    const double vt = vi[t];
                    const double* ftt = frow(t);
                    for (int j = j0; j < j1; ++j)
                        gi[j] -= vt * ftt[j];
                }
            }
        }
        if (recompute)
            norms(s);
    }
    return s;
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
    const auto row = [&](int i) { return g.data() + static_cast<std::size_t>(i) * k; };
    std::vector<int> perm;
    const double rank_tol = 1e3 * std::numeric_limits<double>::epsilon() * std::max(k, m);
    const int r = pivoted_qr(g, m, k, k, rank_tol, perm, x);

    // null vectors as the columns of nv (k x n_null, row-major):
    // rows 0..r-1 hold -R11^-1 R12, rows r..k-1 the identity
    const int n_null = k - r;
    std::vector<double> nv(static_cast<std::size_t>(k) * n_null, 0.0);
    auto nrow = [&](int i) { return nv.data() + static_cast<std::size_t>(i) * n_null; };
    // back substitution in blocks of rows, so that a solved row is read once
    // per block, not once per row above it
    constexpr int ib = 16;
    for (int i1 = r; i1 > 0; i1 -= ib)
    {
        const int i0 = std::max(0, i1 - ib);
        for (int i = i0; i < i1; ++i)
            std::copy_n(row(i) + r, n_null, nrow(i));
        int j = i1;
        for (; j + 3 < r; j += 4)
        {
            const double *y0 = nrow(j), *y1 = nrow(j + 1), *y2 = nrow(j + 2), *y3 = nrow(j + 3);
            for (int i = i0; i < i1; ++i)
            {
                const double* gi = row(i);
                const double r0 = gi[j], r1 = gi[j + 1], r2 = gi[j + 2], r3 = gi[j + 3];
                double* y = nrow(i);
                for (int l = 0; l < n_null; ++l)
                    y[l] -= r0 * y0[l] + r1 * y1[l] + r2 * y2[l] + r3 * y3[l];
            }
        }
        for (; j < r; ++j)
        {
            const double* yj = nrow(j);
            for (int i = i0; i < i1; ++i)
            {
                const double rij = row(i)[j];
                double* y = nrow(i);
                for (int l = 0; l < n_null; ++l)
                    y[l] -= rij * yj[l];
            }
        }
        for (int i = i1 - 1; i >= i0; --i)
        {
            double* y = nrow(i);
            const double* gi = row(i);
            for (int j = i + 1; j < i1; ++j)
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
    // The eliminations' rank-1 updates of the later null vectors wait in
    // blocks of nb, as columns c_t and rows f_t (w_l2 -= c_t f_t[l2]), and reach
    // nv once per block; a step brings up to date only its own null vector and
    // the row it eliminates. Eliminated rows are left stale: they are zero.
    constexpr int nb = 32;
    std::vector<char> eliminated(k, 0);
    std::vector<double> c(static_cast<std::size_t>(k) * nb, 0.0), fb(static_cast<std::size_t>(nb) * n_null),
        v(k);
    const auto frow = [&](int t) { return fb.data() + static_cast<std::size_t>(t) * n_null; };
    int nt = 0, rows = 0; // waiting updates, and the rows they touch
    const auto flush = [&](int from)
    {
        for (int i = 0; i < rows; ++i)
        {
            double* ci = c.data() + static_cast<std::size_t>(i) * nb;
            if (!eliminated[i])
            {
                double* w = nrow(i);
                int t = 0;
                for (; t + 3 < nt; t += 4)
                {
                    const double c0 = ci[t], c1 = ci[t + 1], c2 = ci[t + 2], c3 = ci[t + 3];
                    const double *f0 = frow(t), *f1 = frow(t + 1), *f2 = frow(t + 2), *f3 = frow(t + 3);
                    for (int l2 = from; l2 < n_null; ++l2)
                        w[l2] -= c0 * f0[l2] + c1 * f1[l2] + c2 * f2[l2] + c3 * f3[l2];
                }
                for (; t < nt; ++t)
                {
                    const double ct = ci[t];
                    const double* ft = frow(t);
                    for (int l2 = from; l2 < n_null; ++l2)
                        w[l2] -= ct * ft[l2];
                }
            }
            std::fill_n(ci, nt, 0.0);
        }
        nt = 0;
    };
    for (int l = 0; l < n_null; ++l)
    {
        const int end = r + l + 1; // support of null vector l
        double vmax = 0;
        for (int i = 0; i < end; ++i)
        {
            double vi = 0;
            if (!eliminated[i])
            {
                vi = nrow(i)[l];
                const double* ci = c.data() + static_cast<std::size_t>(i) * nb;
                for (int t = 0; t < nt; ++t)
                    vi -= ci[t] * frow(t)[l];
            }
            v[i] = vi;
            vmax = std::max(vmax, std::abs(vi));
        }
        if (!(vmax > 0))
            continue;
        // the entry that reaches 0 first along +v, or along -v
        int at = -1;
        double step = std::numeric_limits<double>::infinity(), sign = 1;
        for (const double sg : {1.0, -1.0})
        {
            for (int i = 0; i < end; ++i)
            {
                const double vi = sg * v[i];
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
                x[i] = std::max(0.0, x[i] - step * sign * v[i]);
        x[at] = 0;
        alive[at] = 0;
        // later null vectors: w_l2 -= (w_l2[at] / v[at]) v
        if (l + 1 == n_null)
            continue;
        double* f = frow(nt);
        std::fill_n(f, l + 1, 0.0);
        {
            const double* ra = nrow(at);
            const double* ca = c.data() + static_cast<std::size_t>(at) * nb;
            for (int l2 = l + 1; l2 < n_null; ++l2)
                f[l2] = ra[l2];
            for (int t = 0; t < nt; ++t)
            {
                const double ct = ca[t];
                const double* ft = frow(t);
                for (int l2 = l + 1; l2 < n_null; ++l2)
                    f[l2] -= ct * ft[l2];
            }
            const double inv = 1.0 / v[at];
            for (int l2 = l + 1; l2 < n_null; ++l2)
                f[l2] *= inv;
        }
        eliminated[at] = 1;
        for (int i = 0; i < end; ++i)
            c[static_cast<std::size_t>(i) * nb + nt] = v[i];
        ++nt;
        rows = end;
        if (nt == nb)
            flush(l + 1);
    }

    // back to the caller's column order
    std::vector<double> xs(k);
    for (int i = 0; i < k; ++i)
        xs[perm[i]] = x[i];
    x = std::move(xs);
}

/// Least-squares weights on the columns `cols` of a (m rows) for the moments
/// mom, by pivoted_qr with the moments as a last, unpivoted column; false if
/// the columns are numerically dependent.
bool least_squares(const std::vector<double>& a, int m, const std::vector<int>& cols,
                   const std::vector<double>& mom, std::vector<double>& x)
{
    const int n = static_cast<int>(cols.size());
    if (n > m)
        return false;
    const int k = n + 1;
    std::vector<double> g(static_cast<std::size_t>(m) * k);
    for (int j = 0; j < n; ++j)
    {
        const double* c = a.data() + static_cast<std::size_t>(cols[j]) * m;
        for (int i = 0; i < m; ++i)
            g[static_cast<std::size_t>(i) * k + j] = c[i];
    }
    for (int i = 0; i < m; ++i)
        g[static_cast<std::size_t>(i) * k + n] = mom[i];
    // only exact dependence is refused here: the caller keeps the weights
    // only if they are positive and miss the moments by no more
    std::vector<int> perm;
    std::vector<double> none;
    if (pivoted_qr(g, m, k, n, 10 * std::numeric_limits<double>::epsilon(), perm, none) < n)
        return false;
    std::vector<double> y(n);
    for (int s = n - 1; s >= 0; --s)
    {
        const double* gs = g.data() + static_cast<std::size_t>(s) * k;
        double v = gs[n];
        for (int j = s + 1; j < n; ++j)
            v -= gs[j] * y[j];
        y[s] = v / gs[s];
    }
    x.assign(n, 0.0);
    for (int s = 0; s < n; ++s)
        x[perm[s]] = y[s];
    return true;
}

/// Threads for compress_rules: OpenMP's setting where the library is built with
/// it; otherwise CUTCELLS_NUM_THREADS, then OMP_NUM_THREADS, then the hardware's.
int compression_threads()
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    for (const char* name : {"CUTCELLS_NUM_THREADS", "OMP_NUM_THREADS"})
        if (const char* value = std::getenv(name))
        {
            const int n = std::atoi(value);
            if (n > 0)
                return n;
        }
    return static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
#endif
}

/// body(i) for i in [0, n), each index taken by the next free thread in order
/// (dynamic scheduling, one index at a time): OpenMP where the library has it,
/// std::thread otherwise (the macOS wheels are built without OpenMP). The first
/// exception a body throws is rethrown once all threads have finished.
template <typename F>
void parallel_for(int n, F&& body)
{
    const int threads = std::min(compression_threads(), n);
    if (threads <= 1)
    {
        for (int i = 0; i < n; ++i)
            body(i);
        return;
    }
    std::exception_ptr error;
    std::mutex error_mutex;
    const auto guarded = [&](int i)
    {
        try
        {
            body(i);
        }
        catch (...)
        {
            const std::lock_guard<std::mutex> lock(error_mutex);
            if (!error)
                error = std::current_exception();
        }
    };
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (int i = 0; i < n; ++i)
        guarded(i);
#else
    std::atomic<int> next{0};
    const auto work = [&]()
    {
        for (int i = next++; i < n; i = next++)
            guarded(i);
    };
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(threads - 1));
    for (int t = 1; t < threads; ++t)
        pool.emplace_back(work);
    work();
    for (std::thread& t : pool)
        t.join();
#endif
    if (error)
        std::rethrow_exception(error);
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

    Basis<T> basis = make_basis(points, dim, degree, space);
    std::vector<double> w(weights.begin(), weights.end());
    double total = 0;
    for (const double v : w)
        total += v;
    if (!(total > 0))
        return false;
    // the moments, in one pass over the points
    std::vector<double> mom(m, 0.0), column(m);
    std::vector<int> live;
    for (int i = 0; i < n_points; ++i)
    {
        if (!(w[i] > 0))
            continue;
        live.push_back(i);
        basis_column(basis, i, column.data());
        for (int j = 0; j < m; ++j)
            mom[j] += w[i] * column[j];
    }

    // Non-negative least squares on 2m, 4m and 8m points drawn from the rule's
    // measure (systematic sampling of the cumulative weights, so that a
    // region holds candidates in proportion to its measure, not to its points,
    // which crowd where boxes are small): a positive rule on them with the
    // moments of all, if one exists, costs about one QR of an m x m matrix.
    // Recombination below otherwise, and for larger spaces, whose bases are
    // too ill-conditioned on curved pieces for the solver to reach the
    // moments.
    std::vector<double> x, cumulative(live.size());
    {
        double c = 0;
        for (std::size_t p = 0; p < live.size(); ++p)
            cumulative[p] = c += w[live[p]];
    }
    bool found = false;
    double scale = 0;
    for (const double c : mom)
        scale = std::max(scale, std::abs(c));
    for (int k = 2 * m; m <= nnls_moments && !found && k <= 8 * m && 2 * k <= static_cast<int>(live.size()); k *= 2)
    {
        std::vector<int> candidates;
        for (int c = 0; c < k; ++c)
        {
            const double u = (c + 0.5) * cumulative.back() / k;
            const std::size_t p = static_cast<std::size_t>(
                std::lower_bound(cumulative.begin(), cumulative.end(), u) - cumulative.begin());
            const int i = live[std::min(p, live.size() - 1)];
            if (candidates.empty() || candidates.back() != i)
                candidates.push_back(i);
        }
        std::vector<int> set;
        const double err
            = nnls(basis_matrix(basis, candidates), m, static_cast<int>(candidates.size()), mom, set, x);
        found = err <= 1e-13 * scale;
        if (found)
            for (int& c : set)
                c = candidates[c];
        live = found ? std::move(set) : live;
    }
    if (!found)
    {
        // recombination: 2m groups of consecutive points, each replaced by its
        // weighted mean column (g is m x groups, row-major), until the points
        // fit one reduction
        const int groups = 2 * m;
        std::vector<double> g;
        // the columns of the points, stored where they fit in 256 MB, which
        // the rounds read again and again
        std::vector<double> stored;
        if (static_cast<std::size_t>(n_points) * m <= (std::size_t(1) << 25))
        {
            stored.resize(static_cast<std::size_t>(n_points) * m);
            for (const int i : live)
                basis_column(basis, i, stored.data() + static_cast<std::size_t>(i) * m);
        }
        const auto column_of = [&](int i) -> const double*
        {
            if (!stored.empty())
                return stored.data() + static_cast<std::size_t>(i) * m;
            basis_column(basis, i, column.data());
            return column.data();
        };
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
                    const double* ai = column_of(i);
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
            const std::vector<double> a = basis_matrix(basis, live);
            g.resize(static_cast<std::size_t>(k) * m);
            x.resize(k);
            for (int c = 0; c < k; ++c)
            {
                for (int j = 0; j < m; ++j)
                    g[static_cast<std::size_t>(j) * k + c] = a[static_cast<std::size_t>(c) * m + j];
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
    }

    // polish: exact weights on the selected points, if they stay positive
    const std::vector<double> a = basis_matrix(basis, live);
    std::vector<int> cols(live.size());
    std::iota(cols.begin(), cols.end(), 0);
    double err = moment_error(a, m, cols, x, mom);
    std::vector<double> y;
    if (least_squares(a, m, cols, mom, y) && std::all_of(y.begin(), y.end(), [](double v) { return v > 0; }))
    {
        const double err_y = moment_error(a, m, cols, y, mom);
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
    // Largest rules first, one at a time: their cost varies by orders of magnitude
    // between cells, and a large rule handed out last keeps one thread busy alone.
    std::vector<int> order(n_rules);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b)
                     { return rules._offset[a + 1] - rules._offset[a] > rules._offset[b + 1] - rules._offset[b]; });
    parallel_for(n_rules, [&](int i)
    {
        const int r = order[static_cast<std::size_t>(i)];
        const std::size_t b = static_cast<std::size_t>(rules._offset[r]);
        const std::size_t n = static_cast<std::size_t>(rules._offset[r + 1]) - b;
        const std::span<const T> p(rules._points.data() + b * dim, n * dim);
        const std::span<const T> w(rules._weights.data() + b, n);
        if (compress_rule(p, w, dim, degree, space, pts[r], wts[r], res[r]))
            status[r] = 1;
        else if (std::any_of(w.begin(), w.end(), [](T v) { return v < T(0); }))
            status[r] = 2;
    });

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
