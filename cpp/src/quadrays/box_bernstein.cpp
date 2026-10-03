// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "box_bernstein.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>

namespace cutcells::quadrays
{

namespace
{

// ============================================================================
// Univariate Bernstein bases
// ============================================================================

/// All Bernstein basis values B^n_a(s), a = 0..n, by the stable triangular
/// recursion (any s, inside or outside [0, 1]).
template <std::floating_point T>
void basis_values(int n, T s, T* out)
{
    const T t = T(1) - s;
    out[0] = T(1);
    for (int j = 1; j <= n; ++j)
    {
        T carry = T(0);
        for (int a = 0; a < j; ++a)
        {
            const T b = out[a];
            out[a] = carry + t * b;
            carry = s * b;
        }
        out[j] = carry;
    }
}

/// Values and derivatives of all Bernstein basis functions of degree n at s.
template <std::floating_point T>
void basis_values_and_derivatives(int n, T s, T* value, T* deriv)
{
    if (n == 0)
    {
        value[0] = T(1);
        deriv[0] = T(0);
        return;
    }
    basis_values(n - 1, s, value);
    for (int a = 0; a <= n; ++a)
    {
        const T left = a > 0 ? value[a - 1] : T(0);
        const T right = a < n ? value[a] : T(0);
        deriv[a] = T(n) * (left - right);
    }
    basis_values(n, s, value);
}

void check_degree(int degree)
{
    if (degree < 0 || degree > max_box_degree)
    {
        throw std::invalid_argument("quadrays: Bernstein degree "
                                    + std::to_string(degree) + " is out of range");
    }
}

// ============================================================================
// De Casteljau's algorithm with an affine argument
// ============================================================================

/// Z = (1 - l) X + l Y for tensor Bernstein forms X, Y of degree dx in up to
/// three variables and an affine function l. l depends on the variables j with
/// e[j] = 1; corner[v] is its value at the corner v (bit j of v = v_j) of the
/// unit box. Z has degree dx + e.
template <std::floating_point T>
void lerp_affine(const T* X, const T* Y, const std::array<int, 3>& dx,
                 const std::array<int, 3>& e, const std::array<T, 8>& corner, T* Z)
{
    const int nx0 = dx[0] + 1, nx1 = dx[1] + 1, nx2 = dx[2] + 1;
    const int nz0 = nx0 + e[0], nz1 = nx1 + e[1], nz2 = nx2 + e[2];
    for (int k2 = 0; k2 < nz2; ++k2)
        for (int k1 = 0; k1 < nz1; ++k1)
            for (int k0 = 0; k0 < nz0; ++k0)
            {
                T z = T(0);
                for (int v2 = 0; v2 <= e[2]; ++v2)
                {
                    const int j2 = k2 - v2;
                    if (j2 < 0 || j2 >= nx2)
                        continue;
                    const T w2 = e[2] ? (v2 ? T(k2) : T(nx2 - k2)) / T(nx2) : T(1);
                    for (int v1 = 0; v1 <= e[1]; ++v1)
                    {
                        const int j1 = k1 - v1;
                        if (j1 < 0 || j1 >= nx1)
                            continue;
                        const T w1 = e[1] ? (v1 ? T(k1) : T(nx1 - k1)) / T(nx1) : T(1);
                        for (int v0 = 0; v0 <= e[0]; ++v0)
                        {
                            const int j0 = k0 - v0;
                            if (j0 < 0 || j0 >= nx0)
                                continue;
                            const T w0 = e[0] ? (v0 ? T(k0) : T(nx0 - k0)) / T(nx0) : T(1);
                            const T l = corner[v0 + 2 * v1 + 4 * v2];
                            const int idx = j0 + nx0 * (j1 + nx1 * j2);
                            z += w0 * w1 * w2 * ((T(1) - l) * X[idx] + l * Y[idx]);
                        }
                    }
                }
                Z[k0 + nz0 * (k1 + nz1 * k2)] = z;
            }
}

int poly_size(const std::array<int, 3>& d) { return (d[0] + 1) * (d[1] + 1) * (d[2] + 1); }

/// Bernstein coefficients on s in [0, 1] of f(a + (b - a) s), f of degree n
/// with coefficients f[0..n] on [0, 1]: g_k is the blossom of f at
/// (a, ..., a, b, ..., b) with k entries b. @p f is overwritten.
template <std::floating_point T>
void compose_affine_1d(int n, T a, T b, T* f, T* g)
{
    std::array<T, max_box_degree + 1> w;
    for (int k = 0; k <= n; ++k)
    {
        // n - k de Casteljau steps with a on the array after k steps with b
        const int len = n - k + 1;
        std::copy(f, f + len, w.begin());
        for (int r = 1; r < len; ++r)
            for (int j = 0; j < len - r; ++j)
                w[j] = (T(1) - a) * w[j] + a * w[j + 1];
        g[k] = w[0];
        for (int j = 0; j + 1 < len; ++j)
            f[j] = (T(1) - b) * f[j] + b * f[j + 1];
    }
}

/// Contract the variables i of @p p with active[i] = 0 at origin[i]; the others
/// keep their order. Returns the extents of the kept variables in @p ext and
/// writes the coefficients to @p out.
template <std::floating_point T>
int contract_constants(const BoxBernstein<T>& p, std::span<const T> origin, const std::array<int, 3>& active,
                       std::array<int, 3>& kept, std::vector<T>& out, std::vector<T>& tmp)
{
    int n_kept = 0;
    out.assign(p.coeffs.begin(), p.coeffs.begin() + p.size());
    std::array<int, 3> ext = {p.degree[0] + 1, p.degree[1] + 1, p.degree[2] + 1};
    std::array<T, max_box_degree + 1> basis;
    // contract the last variables first so the strides of the earlier ones stay valid
    for (int i = p.dim - 1; i >= 0; --i)
    {
        if (active[i])
            continue;
        int inner = 1, outer = 1;
        for (int r = 0; r < i; ++r)
            inner *= ext[r];
        for (int r = i + 1; r < 3; ++r)
            outer *= ext[r];
        basis_values(p.degree[i], origin[i], basis.data());
        tmp.assign(static_cast<std::size_t>(inner) * outer, T(0));
        for (int o = 0; o < outer; ++o)
            for (int a = 0; a < ext[i]; ++a)
            {
                const T* src = out.data() + static_cast<std::size_t>(o * ext[i] + a) * inner;
                T* dst = tmp.data() + static_cast<std::size_t>(o) * inner;
                for (int in = 0; in < inner; ++in)
                    dst[in] += basis[a] * src[in];
            }
        out.swap(tmp);
        ext[i] = 1;
    }
    for (int i = 0; i < p.dim; ++i)
        if (active[i])
            kept[n_kept++] = i;
    return n_kept;
}

// ============================================================================
// Simplex to box: fixed matrices in exact rational arithmetic
// ============================================================================

__extension__ typedef __int128 int128; // exact sums of the simplex-to-box matrices

int128 gcd128(int128 a, int128 b)
{
    if (a < 0)
        a = -a;
    if (b < 0)
        b = -b;
    while (b != 0)
    {
        const int128 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

int128 factorial(int n)
{
    int128 f = 1;
    for (int k = 2; k <= n; ++k)
        f *= k;
    return f;
}

/// Simplex multi-indices in CutCells' order (bernstein.cpp): alpha_0 is the
/// power of lambda_0 = 1 - sum(xi), alpha_{d+1} the power of xi_d.
std::vector<std::array<int, 4>> simplex_multi_indices(int tdim, int n)
{
    std::vector<std::array<int, 4>> out;
    if (tdim == 2)
    {
        for (int j = 0; j <= n; ++j)
            for (int i = 0; i <= n - j; ++i)
                out.push_back({n - i - j, i, j, 0});
    }
    else
    {
        for (int k = 0; k <= n; ++k)
            for (int j = 0; j <= n - k; ++j)
                for (int i = 0; i <= n - k - j; ++i)
                    out.push_back({n - i - j - k, i, j, k});
    }
    return out;
}

/// Matrix (box index x simplex index, row-major) taking simplex Bernstein
/// coefficients of degree n to tensor Bernstein coefficients of degree n per
/// variable on the unit box. Each simplex basis function is expanded in
/// monomials, x^m = sum_i C(i, m) / C(n, m) B^n_i(x) per variable; all sums
/// are exact over the common denominator (n!)^tdim and rounded once.
std::vector<long double> simplex_to_box_exact(int tdim, int n)
{
    const std::vector<std::array<int, 4>> alphas = simplex_multi_indices(tdim, n);
    const int ns = static_cast<int>(alphas.size());
    int nb = 1;
    for (int d = 0; d < tdim; ++d)
        nb *= n + 1;

    std::vector<int128> fact(n + 1);
    for (int k = 0; k <= n; ++k)
        fact[k] = factorial(k);
    // ff[i][m] = i! / (i - m)! (n - m)!, the numerator of C(i, m) / C(n, m) over n!
    std::vector<int128> ff((n + 1) * (n + 1), 0);
    for (int i = 0; i <= n; ++i)
        for (int m = 0; m <= i; ++m)
            ff[i * (n + 1) + m] = fact[i] / fact[i - m] * fact[n - m];

    std::vector<int128> numer(static_cast<std::size_t>(nb) * ns, 0);
    for (int s = 0; s < ns; ++s)
    {
        const std::array<int, 4>& a = alphas[s];
        const int128 multinomial = fact[n] / (fact[a[0]] * fact[a[1]] * fact[a[2]] * fact[a[3]]);
        // lambda_0^a0 = (1 - x_0 - x_1 - x_2)^a0 = sum over b of a0! / b! (-1)^(b1+b2+b3) x^b
        for (int b1 = 0; b1 <= a[0]; ++b1)
            for (int b2 = 0; b2 <= a[0] - b1; ++b2)
                for (int b3 = 0; b3 <= (tdim == 3 ? a[0] - b1 - b2 : 0); ++b3)
                {
                    const int b0 = a[0] - b1 - b2 - b3;
                    int128 coef = multinomial * (fact[a[0]] / (fact[b0] * fact[b1] * fact[b2] * fact[b3]));
                    if ((b1 + b2 + b3) % 2 == 1)
                        coef = -coef;
                    const std::array<int, 3> m = {a[1] + b1, a[2] + b2, a[3] + b3};
                    for (int box = 0; box < nb; ++box)
                    {
                        int128 term = coef;
                        int rest = box;
                        for (int d = 0; d < tdim && term != 0; ++d)
                        {
                            const int i = rest % (n + 1);
                            rest /= n + 1;
                            term = i >= m[d] ? term * ff[i * (n + 1) + m[d]] : 0;
                        }
                        numer[static_cast<std::size_t>(box) * ns + s] += term;
                    }
                }
    }

    int128 denom = 1;
    for (int d = 0; d < tdim; ++d)
        denom *= fact[n];
    std::vector<long double> matrix(numer.size());
    for (std::size_t k = 0; k < numer.size(); ++k)
    {
        const int128 g = gcd128(numer[k], denom);
        const int128 num = g > 0 ? numer[k] / g : numer[k];
        const int128 den = g > 0 ? denom / g : denom;
        matrix[k] = static_cast<long double>(num) / static_cast<long double>(den);
    }
    return matrix;
}

/// Cached simplex-to-box matrix of type T.
template <std::floating_point T>
const std::vector<T>& simplex_to_box_matrix(int tdim, int n)
{
    static std::mutex mutex;
    static std::map<std::pair<int, int>, std::vector<T>> cache;
    const std::lock_guard<std::mutex> lock(mutex);
    auto it = cache.find({tdim, n});
    if (it == cache.end())
    {
        const std::vector<long double> exact = simplex_to_box_exact(tdim, n);
        it = cache.emplace(std::pair{tdim, n}, std::vector<T>(exact.begin(), exact.end())).first;
    }
    return it->second;
}

/// Matrix (box index x pyramid index, row-major) taking the coefficients of a
/// pyramid's level set phi in its basis B^m_i(s) B^m_j(t) B^n_k(z), m = n - k,
/// s = x / (1 - z), t = y / (1 - z) (bernstein.h), to the tensor Bernstein
/// coefficients of (1 - z)^n phi, of degree (n, n, 2n), on the unit box:
///   (1 - z)^n B^m_i(s) B^m_j(t) B^n_k(z)
///     = C(m, i) C(m, j) C(n, k) x^i (1 - z - x)^(m - i) y^j (1 - z - y)^(m - j) z^k (1 - z)^k,
/// expanded in powers of x and y, with x^a = sum_l C(l, a) / C(n, a) B^n_l(x)
/// and z^k (1 - z)^e = sum_l C(2n - k - e, l - k) / C(2n, l) B^2n_l(z). All
/// sums are exact over the common denominator (n!)^2 (2n)! and rounded once.
std::vector<long double> pyramid_to_box_exact(int n)
{
    const int N = 2 * n, n1 = n + 1;
    const int nb = n1 * n1 * (N + 1);
    int np = 0;
    for (int k = 0; k <= n; ++k)
        np += (n - k + 1) * (n - k + 1);
    std::vector<int128> fact(N + 1);
    for (int k = 0; k <= N; ++k)
        fact[k] = factorial(k);
    auto choose = [&fact](int a, int b) -> int128 { return b < 0 || b > a ? 0 : fact[a] / (fact[b] * fact[a - b]); };
    // ff[l][a] = l! / (l - a)! (n - a)!: C(l, a) / C(n, a) times n!
    std::vector<int128> ff(n1 * n1, 0);
    for (int l = 0; l <= n; ++l)
        for (int a = 0; a <= l; ++a)
            ff[l * n1 + a] = fact[l] / fact[l - a] * fact[n - a];

    std::vector<int128> numer(static_cast<std::size_t>(nb) * np, 0);
    int idx = 0;
    for (int k = 0; k <= n; ++k)
    {
        const int m = n - k;
        for (int i = 0; i <= m; ++i)
            for (int j = 0; j <= m; ++j, ++idx)
            {
                const int128 base = choose(m, i) * choose(m, j) * choose(n, k);
                for (int p = 0; p <= m - i; ++p)
                    for (int q = 0; q <= m - j; ++q)
                    {
                        int128 coef = base * choose(m - i, p) * choose(m - j, q);
                        if ((p + q) % 2 == 1)
                            coef = -coef;
                        const int a = i + p, b = j + q, e = (m - i - p) + (m - j - q) + k;
                        for (int l2 = k; l2 <= N - e; ++l2)
                        {
                            // C(N - k - e, l2 - k) / C(N, l2) times N!
                            const int128 zl = choose(N - k - e, l2 - k) * fact[l2] * fact[N - l2];
                            for (int l1 = b; l1 <= n; ++l1)
                                for (int l0 = a; l0 <= n; ++l0)
                                    numer[static_cast<std::size_t>(l0 + n1 * (l1 + n1 * l2)) * np + idx]
                                        += coef * ff[l0 * n1 + a] * ff[l1 * n1 + b] * zl;
                        }
                    }
            }
    }
    const int128 denom = fact[n] * fact[n] * fact[N];
    std::vector<long double> matrix(numer.size());
    for (std::size_t k = 0; k < numer.size(); ++k)
    {
        const int128 g = gcd128(numer[k], denom);
        const int128 num = g > 0 ? numer[k] / g : numer[k];
        const int128 den = g > 0 ? denom / g : denom;
        matrix[k] = static_cast<long double>(num) / static_cast<long double>(den);
    }
    return matrix;
}

/// Cached pyramid-to-box matrix of type T.
template <std::floating_point T>
const std::vector<T>& pyramid_to_box_matrix(int n)
{
    static std::mutex mutex;
    static std::map<int, std::vector<T>> cache;
    const std::lock_guard<std::mutex> lock(mutex);
    auto it = cache.find(n);
    if (it == cache.end())
    {
        const std::vector<long double> exact = pyramid_to_box_exact(n);
        it = cache.emplace(n, std::vector<T>(exact.begin(), exact.end())).first;
    }
    return it->second;
}

} // namespace

// ============================================================================
// Conversion from a cell's level set
// ============================================================================

namespace
{
/// The form on the box of the reference frame (cell_bernstein_on_box).
template <std::floating_point T>
void reference_form(cell::type cell_type, int degree, std::span<const T> coeffs, BoxBernstein<T>& out)
{
    check_degree(degree);
    const int n1 = degree + 1;
    switch (cell_type)
    {
    case cell::type::triangle:
    case cell::type::tetrahedron:
    {
        const int tdim = cell_type == cell::type::triangle ? 2 : 3;
        if (degree > 12)
        {
            throw std::invalid_argument(
                "quadrays: simplex level sets of degree above 12 are not supported");
        }
        const std::vector<T>& matrix = simplex_to_box_matrix<T>(tdim, degree);
        const int nb = tdim == 2 ? n1 * n1 : n1 * n1 * n1;
        const int ns = static_cast<int>(matrix.size()) / nb;
        if (static_cast<int>(coeffs.size()) != ns)
            throw std::invalid_argument("quadrays: wrong number of simplex Bernstein coefficients");
        out.dim = tdim;
        out.degree = {degree, degree, tdim == 3 ? degree : 0};
        out.coeffs.assign(nb, T(0));
        for (int b = 0; b < nb; ++b)
        {
            const T* row = matrix.data() + static_cast<std::size_t>(b) * ns;
            T sum = T(0);
            for (int s = 0; s < ns; ++s)
                sum += row[s] * coeffs[s];
            out.coeffs[b] = sum;
        }
        return;
    }
    case cell::type::quadrilateral:
    {
        if (static_cast<int>(coeffs.size()) != n1 * n1)
            throw std::invalid_argument("quadrays: wrong number of tensor Bernstein coefficients");
        out.dim = 2;
        out.degree = {degree, degree, 0};
        out.coeffs.resize(n1 * n1);
        // CutCells stores coeffs[i * n1 + j] with i along xi_0
        for (int a1 = 0; a1 < n1; ++a1)
            for (int a0 = 0; a0 < n1; ++a0)
                out.coeffs[a0 + n1 * a1] = coeffs[a0 * n1 + a1];
        return;
    }
    case cell::type::hexahedron:
    {
        if (static_cast<int>(coeffs.size()) != n1 * n1 * n1)
            throw std::invalid_argument("quadrays: wrong number of tensor Bernstein coefficients");
        out.dim = 3;
        out.degree = {degree, degree, degree};
        out.coeffs.resize(n1 * n1 * n1);
        // CutCells stores coeffs[(i * n1 + j) * n1 + k] with i along xi_0
        for (int a2 = 0; a2 < n1; ++a2)
            for (int a1 = 0; a1 < n1; ++a1)
                for (int a0 = 0; a0 < n1; ++a0)
                    out.coeffs[a0 + n1 * (a1 + n1 * a2)] = coeffs[(a0 * n1 + a1) * n1 + a2];
        return;
    }
    case cell::type::prism:
    {
        // each layer k along xi_2 is a triangle form: CutCells stores
        // coeffs[a * n1 + k], a the triangle's index
        if (degree > 12)
            throw std::invalid_argument("quadrays: prism level sets of degree above 12 are not supported");
        const std::vector<T>& matrix = simplex_to_box_matrix<T>(2, degree);
        const int nb = n1 * n1;
        const int ns = static_cast<int>(matrix.size()) / nb;
        if (static_cast<int>(coeffs.size()) != ns * n1)
            throw std::invalid_argument("quadrays: wrong number of prism Bernstein coefficients");
        out.dim = 3;
        out.degree = {degree, degree, degree};
        out.coeffs.assign(nb * n1, T(0));
        for (int k = 0; k < n1; ++k)
            for (int b = 0; b < nb; ++b)
            {
                const T* row = matrix.data() + static_cast<std::size_t>(b) * ns;
                T sum = T(0);
                for (int a = 0; a < ns; ++a)
                    sum += row[a] * coeffs[a * n1 + k];
                out.coeffs[b + nb * k] = sum;
            }
        return;
    }
    case cell::type::pyramid:
    {
        // the level set times (1 - z)^n, a polynomial (pyramid_to_box_exact)
        if (degree > 6)
            throw std::invalid_argument("quadrays: pyramid level sets of degree above 6 are not supported");
        const std::vector<T>& matrix = pyramid_to_box_matrix<T>(degree);
        const int nb = n1 * n1 * (2 * degree + 1);
        const int np = static_cast<int>(matrix.size()) / nb;
        if (static_cast<int>(coeffs.size()) != np)
            throw std::invalid_argument("quadrays: wrong number of pyramid Bernstein coefficients");
        out.dim = 3;
        out.degree = {degree, degree, 2 * degree};
        out.coeffs.assign(nb, T(0));
        for (int b = 0; b < nb; ++b)
        {
            const T* row = matrix.data() + static_cast<std::size_t>(b) * np;
            T sum = T(0);
            for (int a = 0; a < np; ++a)
                sum += row[a] * coeffs[a];
            out.coeffs[b] = sum;
        }
        return;
    }
    default:
        throw std::invalid_argument("quadrays: unsupported cell type "
                                    + cell::cell_type_to_str(cell_type));
    }
}

/// Inverse of the Bernstein-Vandermonde matrix of degree n at the equispaced
/// nodes i / n (row-major): coefficients from values. Computed once per
/// degree in long double.
template <std::floating_point T>
const std::vector<T>& equispaced_bernstein_inverse(int n)
{
    // per thread, the pointers into the shared cache (whose nodes stay put)
    thread_local std::array<const std::vector<T>*, max_box_degree + 1> local{};
    if (n >= 0 && n <= max_box_degree && local[static_cast<std::size_t>(n)] != nullptr)
        return *local[static_cast<std::size_t>(n)];
    static std::mutex mutex;
    static std::map<int, std::vector<T>> cache;
    const std::lock_guard<std::mutex> lock(mutex);
    auto it = cache.find(n);
    if (it == cache.end())
    {
        const int n1 = n + 1;
        std::vector<long double> a(static_cast<std::size_t>(n1 * n1)), inv(a.size(), 0.0L);
        for (int i = 0; i < n1; ++i)
        {
            const long double s = n == 0 ? 0.5L : static_cast<long double>(i) / n;
            long double binom = 1.0L;
            for (int j = 0; j < n1; ++j)
            {
                a[static_cast<std::size_t>(i * n1 + j)] = binom * std::pow(s, j) * std::pow(1.0L - s, n - j);
                binom = binom * (n - j) / (j + 1);
            }
            inv[static_cast<std::size_t>(i * n1 + i)] = 1.0L;
        }
        // Gauss-Jordan elimination with partial pivoting
        for (int c = 0; c < n1; ++c)
        {
            int pivot = c;
            for (int r = c + 1; r < n1; ++r)
                if (std::abs(a[static_cast<std::size_t>(r * n1 + c)])
                    > std::abs(a[static_cast<std::size_t>(pivot * n1 + c)]))
                    pivot = r;
            for (int k = 0; k < n1; ++k)
            {
                std::swap(a[static_cast<std::size_t>(c * n1 + k)], a[static_cast<std::size_t>(pivot * n1 + k)]);
                std::swap(inv[static_cast<std::size_t>(c * n1 + k)], inv[static_cast<std::size_t>(pivot * n1 + k)]);
            }
            const long double d = a[static_cast<std::size_t>(c * n1 + c)];
            for (int k = 0; k < n1; ++k)
            {
                a[static_cast<std::size_t>(c * n1 + k)] /= d;
                inv[static_cast<std::size_t>(c * n1 + k)] /= d;
            }
            for (int r = 0; r < n1; ++r)
            {
                if (r == c)
                    continue;
                const long double f = a[static_cast<std::size_t>(r * n1 + c)];
                for (int k = 0; k < n1; ++k)
                {
                    a[static_cast<std::size_t>(r * n1 + k)] -= f * a[static_cast<std::size_t>(c * n1 + k)];
                    inv[static_cast<std::size_t>(r * n1 + k)] -= f * inv[static_cast<std::size_t>(c * n1 + k)];
                }
            }
        }
        it = cache.emplace(n, std::vector<T>(inv.begin(), inv.end())).first;
    }
    if (n >= 0 && n <= max_box_degree)
        local[static_cast<std::size_t>(n)] = &it->second;
    return it->second;
}
} // namespace

template <std::floating_point T>
void cell_bernstein_on_box(cell::type cell_type, int degree, std::span<const T> coeffs,
                           const ClippedBox<T>& box, BoxBernstein<T>& out)
{
    if (reference_frame(box))
    {
        reference_form(cell_type, degree, coeffs, out);
        return;
    }
    if (cell_type != cell::type::triangle && cell_type != cell::type::tetrahedron && cell_type != cell::type::prism)
    {
        throw std::invalid_argument("quadrays: a " + cell::cell_type_to_str(cell_type)
                                    + " takes the box of its reference frame");
    }
    thread_local BoxBernstein<T> ref;
    reference_form(cell_type, degree, coeffs, ref);

    // degree n in each box variable (cell_bernstein_on_box)
    const int dim = box.tdim, n = degree;
    const std::array<int, 3> deg = {n, n, dim == 3 ? n : 0};
    const std::array<int, 3> ext = {deg[0] + 1, deg[1] + 1, deg[2] + 1};

    // values at the equispaced nodes of the box, then coefficients axis by axis
    out.dim = dim;
    out.degree = deg;
    out.coeffs.resize(static_cast<std::size_t>(ext[0] * ext[1] * ext[2]));
    auto node = [](int i, int d) { return d == 0 ? T(0.5) : static_cast<T>(i) / static_cast<T>(d); };
    for (int a2 = 0; a2 < ext[2]; ++a2)
        for (int a1 = 0; a1 < ext[1]; ++a1)
            for (int a0 = 0; a0 < ext[0]; ++a0)
            {
                const Vec3<T> u = {node(a0, deg[0]), node(a1, deg[1]), dim == 3 ? node(a2, deg[2]) : T(0)};
                const Vec3<T> xi = reference_point(box, u);
                out.coeffs[static_cast<std::size_t>(a0 + ext[0] * (a1 + ext[1] * a2))]
                    = evaluate(ref, std::span<const T>(xi.data(), static_cast<std::size_t>(dim)));
            }
    thread_local std::vector<T> line;
    const std::array<int, 3> stride = {1, ext[0], ext[0] * ext[1]};
    for (int axis = 0; axis < dim; ++axis)
    {
        const std::vector<T>& inv = equispaced_bernstein_inverse<T>(deg[axis]);
        const int m = ext[axis];
        line.resize(static_cast<std::size_t>(m));
        const int o1 = axis == 0 ? 1 : 0, o2 = axis == 2 ? 1 : 2; // the other two axes
        for (int i2 = 0; i2 < ext[o2]; ++i2)
            for (int i1 = 0; i1 < ext[o1]; ++i1)
            {
                const int base = i1 * stride[o1] + i2 * stride[o2];
                for (int i = 0; i < m; ++i)
                    line[static_cast<std::size_t>(i)] = out.coeffs[static_cast<std::size_t>(base + i * stride[axis])];
                for (int k = 0; k < m; ++k)
                {
                    T sum = T(0);
                    for (int i = 0; i < m; ++i)
                        sum += inv[static_cast<std::size_t>(k * m + i)] * line[static_cast<std::size_t>(i)];
                    out.coeffs[static_cast<std::size_t>(base + k * stride[axis])] = sum;
                }
            }
    }
}

// ============================================================================
// Evaluation
// ============================================================================

template <std::floating_point T>
T evaluate(const BoxBernstein<T>& p, std::span<const T> s)
{
    std::array<std::array<T, max_box_degree + 1>, 3> b;
    for (int j = 0; j < 3; ++j)
    {
        check_degree(p.degree[j]);
        if (j < p.dim)
            basis_values(p.degree[j], s[j], b[j].data());
        else
            b[j][0] = T(1);
    }
    const int n0 = p.degree[0] + 1, n1 = p.degree[1] + 1, n2 = p.degree[2] + 1;
    T result = T(0);
    for (int a2 = 0; a2 < n2; ++a2)
    {
        T s1 = T(0);
        for (int a1 = 0; a1 < n1; ++a1)
        {
            const T* c = p.coeffs.data() + n0 * (a1 + n1 * a2);
            T s0 = T(0);
            for (int a0 = 0; a0 < n0; ++a0)
                s0 += c[a0] * b[0][a0];
            s1 += s0 * b[1][a1];
        }
        result += s1 * b[2][a2];
    }
    return result;
}

template <std::floating_point T>
void gradient(const BoxBernstein<T>& p, std::span<const T> s, std::span<T> grad)
{
    std::array<std::array<T, max_box_degree + 1>, 3> b, db;
    for (int j = 0; j < 3; ++j)
    {
        check_degree(p.degree[j]);
        if (j < p.dim)
            basis_values_and_derivatives(p.degree[j], s[j], b[j].data(), db[j].data());
        else
        {
            b[j][0] = T(1);
            db[j][0] = T(0);
        }
    }
    const int n0 = p.degree[0] + 1, n1 = p.degree[1] + 1, n2 = p.degree[2] + 1;
    T g0 = T(0), g1 = T(0), g2 = T(0);
    for (int a2 = 0; a2 < n2; ++a2)
    {
        T v = T(0), d0 = T(0), d1 = T(0);
        for (int a1 = 0; a1 < n1; ++a1)
        {
            const T* c = p.coeffs.data() + n0 * (a1 + n1 * a2);
            T s0 = T(0), ds0 = T(0);
            for (int a0 = 0; a0 < n0; ++a0)
            {
                s0 += c[a0] * b[0][a0];
                ds0 += c[a0] * db[0][a0];
            }
            v += s0 * b[1][a1];
            d0 += ds0 * b[1][a1];
            d1 += s0 * db[1][a1];
        }
        g0 += d0 * b[2][a2];
        g1 += d1 * b[2][a2];
        g2 += v * db[2][a2];
    }
    const std::array<T, 3> g = {g0, g1, g2};
    for (int j = 0; j < p.dim; ++j)
        grad[j] = g[j];
}

template <std::floating_point T>
void derivative(const BoxBernstein<T>& p, int direction, BoxBernstein<T>& out)
{
    out.dim = p.dim;
    out.degree = p.degree;
    const int nk = p.degree[direction];
    if (nk == 0)
    {
        out.coeffs.assign(p.size(), T(0));
        return;
    }
    out.degree[direction] = nk - 1;
    const std::array<int, 3> ext = {p.degree[0] + 1, p.degree[1] + 1, p.degree[2] + 1};
    const std::array<int, 3> stride = {1, ext[0], ext[0] * ext[1]};
    std::array<int, 3> oext = ext;
    oext[direction] = nk;
    out.coeffs.resize(oext[0] * oext[1] * oext[2]);
    const int step = stride[direction];
    int o = 0;
    for (int a2 = 0; a2 < oext[2]; ++a2)
        for (int a1 = 0; a1 < oext[1]; ++a1)
            for (int a0 = 0; a0 < oext[0]; ++a0, ++o)
            {
                const int i = a0 + ext[0] * (a1 + ext[1] * a2);
                out.coeffs[o] = T(nk) * (p.coeffs[i + step] - p.coeffs[i]);
            }
}

// ============================================================================
// Affine restriction and subdivision
// ============================================================================

template <std::floating_point T>
void restrict_affine(const BoxBernstein<T>& p, std::span<const T> origin,
                     std::span<const T> matrix, int m, BoxBernstein<T>& out,
                     std::vector<T>& work)
{
    const int dim = p.dim;
    if (m < 0 || m > 3)
        throw std::invalid_argument("quadrays: restrict_affine needs 0 <= m <= 3");

    // variables of s each old variable depends on
    std::array<std::array<int, 3>, 3> active{};
    std::array<int, 3> n_active{};
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < m; ++j)
        {
            active[i][j] = matrix[i * m + j] != T(0) ? 1 : 0;
            n_active[i] += active[i][j];
        }

    // Fast path: every old variable depends on at most one new variable, in
    // increasing order (sub-boxes, box faces, lines along an axis). Constant
    // variables are evaluated; the others are reparametrised along their axis.
    std::array<int, 3> target = {-1, -1, -1};
    bool aligned = true;
    for (int i = 0, last = -1; i < dim && aligned; ++i)
    {
        if (n_active[i] > 1)
            aligned = false;
        for (int j = 0; j < m && aligned; ++j)
            if (active[i][j])
            {
                aligned = j > last;
                last = target[i] = j;
            }
    }
    const std::array<int, 3> is_active = {n_active[0] > 0, n_active[1] > 0, n_active[2] > 0};
    if (aligned || m == 1)
    {
        std::array<int, 3> kept{};
        std::vector<T>& coeffs = out.coeffs;
        const int n_kept = contract_constants(p, origin, is_active, kept, coeffs, work);
        std::array<T, max_box_degree + 1> fiber, result;
        if (aligned)
        {
            int inner = 1;
            for (int r = 0; r < n_kept; ++r)
            {
                const int i = kept[r], n = p.degree[i];
                int outer = 1;
                for (int rr = r + 1; rr < n_kept; ++rr)
                    outer *= p.degree[kept[rr]] + 1;
                const T a = origin[i], b = origin[i] + matrix[i * m + target[i]];
                for (int o = 0; o < outer; ++o)
                    for (int in = 0; in < inner; ++in)
                    {
                        T* base = coeffs.data() + static_cast<std::size_t>(o) * (n + 1) * inner + in;
                        for (int k = 0; k <= n; ++k)
                            fiber[k] = base[static_cast<std::size_t>(k) * inner];
                        compose_affine_1d(n, a, b, fiber.data(), result.data());
                        for (int k = 0; k <= n; ++k)
                            base[static_cast<std::size_t>(k) * inner] = result[k];
                    }
                inner *= n + 1;
            }
            out.dim = m;
            out.degree = {0, 0, 0};
            for (int r = 0; r < n_kept; ++r)
                out.degree[target[kept[r]]] = p.degree[kept[r]];
            return;
        }

        // A line through several variables: univariate polynomials in s, one
        // variable substituted at a time by de Casteljau's algorithm with the
        // argument a + (b - a) s; the fastest variable goes first.
        int count = static_cast<int>(coeffs.size()), d = 0; // polynomials of degree d
        std::vector<T>& next = work;
        for (int r = 0; r < n_kept; ++r)
        {
            const int i = kept[r], n = p.degree[i];
            const T a = origin[i], b = origin[i] + matrix[i * m];
            const int dn = d + n;
            check_degree(dn);
            const int count_next = count / (n + 1);
            next.assign(static_cast<std::size_t>(count_next) * (dn + 1), T(0));
            thread_local std::vector<std::array<T, max_box_degree + 1>> slots;
            slots.resize(n + 1);
            for (int o = 0; o < count_next; ++o)
            {
                for (int k = 0; k <= n; ++k)
                    std::copy_n(coeffs.data() + static_cast<std::size_t>(o * (n + 1) + k) * (d + 1), d + 1,
                                slots[k].begin());
                for (int t = 1; t <= n; ++t)
                {
                    const int dx = d + t - 1; // degree before this step
                    for (int k = 0; k <= n - t; ++k)
                    {
                        T* X = slots[k].data();
                        const T* Y = slots[k + 1].data();
                        for (int c = dx + 1; c >= 0; --c)
                        {
                            T z = T(0);
                            if (c <= dx)
                                z += T(dx + 1 - c) / T(dx + 1) * ((T(1) - a) * X[c] + a * Y[c]);
                            if (c >= 1)
                                z += T(c) / T(dx + 1) * ((T(1) - b) * X[c - 1] + b * Y[c - 1]);
                            X[c] = z;
                        }
                    }
                }
                std::copy_n(slots[0].begin(), dn + 1, next.begin() + static_cast<std::size_t>(o) * (dn + 1));
            }
            coeffs.swap(next);
            count = count_next;
            d = dn;
        }
        out.dim = 1;
        out.degree = {d, 0, 0};
        return;
    }

    // Substitute constant variables first and the most mixed ones last: they
    // raise the degree of every polynomial still to be combined.
    std::array<int, 3> order = {0, 1, 2};
    std::stable_sort(order.begin(), order.begin() + dim,
                     [&](int a, int b) { return n_active[a] < n_active[b]; });

    std::array<int, 3> rem = {0, 1, 2}; // variables not yet substituted, in storage order
    int n_rem = dim;
    std::array<int, 3> ext = {1, 1, 1};
    for (int i = 0; i < dim; ++i)
        ext[i] = p.degree[i] + 1;

    std::array<int, 3> d = {0, 0, 0}; // degree in s of every stored polynomial
    int ps = 1;                       // size of every stored polynomial
    int count = p.size();             // number of stored polynomials
    work.resize(static_cast<std::size_t>(count));
    std::copy(p.coeffs.begin(), p.coeffs.begin() + count, work.begin());

    std::array<T, max_box_degree + 1> basis;
    for (int step = 0; step < dim; ++step)
    {
        const int i = order[step];
        const int n = p.degree[i];
        check_degree(n);
        int pos = 0;
        while (rem[pos] != i)
            ++pos;
        int inner = 1, outer = 1;
        for (int r = 0; r < pos; ++r)
            inner *= ext[rem[r]];
        for (int r = pos + 1; r < n_rem; ++r)
            outer *= ext[rem[r]];

        const std::array<int, 3> e = {active[i][0], active[i][1], active[i][2]};
        std::array<int, 3> dn = d;
        for (int j = 0; j < 3; ++j)
            dn[j] += n * e[j];
        const int psn = poly_size(dn);

        // corner values of the affine argument u_i(s) on the unit box
        std::array<T, 8> corner{};
        for (int v = 0; v < 8; ++v)
        {
            T value = origin[i];
            for (int j = 0; j < m; ++j)
                if (e[j] && ((v >> j) & 1))
                    value += matrix[i * m + j];
            corner[v] = value;
        }

        // layout of work: [current | next | n + 1 slots | temporary]
        const std::size_t cur_size = static_cast<std::size_t>(count) * ps;
        const int count_next = count / (n + 1);
        const std::size_t next_off = cur_size;
        const std::size_t slot_off = next_off + static_cast<std::size_t>(count_next) * psn;
        const std::size_t tmp_off = slot_off + static_cast<std::size_t>(n + 1) * psn;
        work.resize(tmp_off + psn);
        T* cur = work.data();
        T* next = work.data() + next_off;
        T* slots = work.data() + slot_off;
        T* tmp = work.data() + tmp_off;

        if (n_active[i] == 0)
            basis_values(n, origin[i], basis.data());
        for (int o = 0; o < outer; ++o)
            for (int in = 0; in < inner; ++in)
            {
                T* dst = next + static_cast<std::size_t>(o * inner + in) * psn;
                auto fiber = [&](int a)
                { return cur + static_cast<std::size_t>((o * (n + 1) + a) * inner + in) * ps; };
                if (n_active[i] == 0)
                {
                    // constant argument: evaluate along this variable
                    std::fill(dst, dst + psn, T(0));
                    for (int a = 0; a <= n; ++a)
                    {
                        const T* f = fiber(a);
                        for (int k = 0; k < ps; ++k)
                            dst[k] += basis[a] * f[k];
                    }
                    continue;
                }
                for (int a = 0; a <= n; ++a)
                    std::copy(fiber(a), fiber(a) + ps, slots + static_cast<std::size_t>(a) * psn);
                std::array<int, 3> dt = d;
                for (int t = 1; t <= n; ++t)
                {
                    const int size_t1 = poly_size({dt[0] + e[0], dt[1] + e[1], dt[2] + e[2]});
                    for (int a = 0; a <= n - t; ++a)
                    {
                        lerp_affine(slots + static_cast<std::size_t>(a) * psn,
                                    slots + static_cast<std::size_t>(a + 1) * psn, dt, e, corner, tmp);
                        std::copy(tmp, tmp + size_t1, slots + static_cast<std::size_t>(a) * psn);
                    }
                    for (int j = 0; j < 3; ++j)
                        dt[j] += e[j];
                }
                std::copy(slots, slots + psn, dst);
            }

        // move next to the front
        std::copy(next, next + static_cast<std::size_t>(count_next) * psn, work.data());
        count = count_next;
        ps = psn;
        d = dn;
        for (int r = pos; r + 1 < n_rem; ++r)
            rem[r] = rem[r + 1];
        --n_rem;
    }

    out.dim = m;
    out.degree = {0, 0, 0};
    for (int j = 0; j < m; ++j)
    {
        out.degree[j] = d[j];
        check_degree(d[j]);
    }
    out.coeffs.assign(work.begin(), work.begin() + ps);
}

template <std::floating_point T>
void subdivide(const BoxBernstein<T>& p, std::span<const T> lo, std::span<const T> hi,
               BoxBernstein<T>& out, std::vector<T>& work)
{
    std::array<T, 9> matrix{};
    for (int j = 0; j < p.dim; ++j)
        matrix[j * p.dim + j] = hi[j] - lo[j];
    restrict_affine(p, lo, std::span<const T>(matrix.data(), p.dim * p.dim), p.dim, out, work);
    // keep the degree of p even where the sub-box is degenerate
    if (out.degree != p.degree)
    {
        throw std::invalid_argument("quadrays: subdivide needs a sub-box of positive size");
    }
}

// ============================================================================
// Bounds
// ============================================================================

template <std::floating_point T>
bool may_vanish(std::span<const T> coeffs)
{
    bool pos = false, neg = false;
    for (const T v : coeffs)
    {
        pos |= v > T(0);
        neg |= v < T(0);
    }
    return pos && neg;
}

template <std::floating_point T>
T max_abs(std::span<const T> coeffs)
{
    T m = T(0);
    for (const T v : coeffs)
        m = std::max(m, std::abs(v));
    return m;
}

template <std::floating_point T>
T scaled_norm(std::span<const T> v)
{
    T m = T(0);
    for (const T x : v)
        m = std::max(m, std::abs(x));
    if (m == T(0) || !std::isfinite(m))
        return m;
    T s = T(0);
    for (const T x : v)
        s += (x / m) * (x / m);
    return m * std::sqrt(s);
}

template <std::floating_point T>
void margins(const BoxBernstein<T>& p, std::span<const T> lengths, std::span<T> ratio,
             BoxBernstein<T>& work, std::span<T> upper_out)
{
    std::array<T, 3> lower{}, upper{};
    for (int k = 0; k < p.dim; ++k)
    {
        if (p.degree[k] == 0)
            continue;
        derivative(p, k, work);
        bool pos = true, neg = true;
        T amin = std::numeric_limits<T>::infinity(), amax = T(0);
        for (const T v : work.coeffs)
        {
            pos &= v > T(0);
            neg &= v < T(0);
            amin = std::min(amin, std::abs(v));
            amax = std::max(amax, std::abs(v));
        }
        const T scale = T(1) / lengths[k];
        lower[k] = (pos || neg) ? amin * scale : T(0);
        upper[k] = amax * scale;
    }
    const T norm = scaled_norm(std::span<const T>(upper.data(), p.dim));
    for (int k = 0; k < p.dim; ++k)
        ratio[k] = norm > T(0) ? lower[k] / norm : T(0);
    for (std::size_t k = 0; k < upper_out.size(); ++k)
        upper_out[k] = upper[k];
}

// ============================================================================
// Univariate roots
// ============================================================================

template <std::floating_point T>
T evaluate_1d(std::span<const T> c, T s)
{
    const int n = static_cast<int>(c.size()) - 1;
    check_degree(n);
    std::array<T, max_box_degree + 1> w;
    std::copy(c.begin(), c.end(), w.begin());
    const T t = T(1) - s;
    for (int r = 1; r <= n; ++r)
        for (int a = 0; a <= n - r; ++a)
            w[a] = t * w[a] + s * w[a + 1];
    return w[0];
}

template <std::floating_point T>
T bracketed_root(std::span<const T> c, T a, T b, T lo, T hi, T glo, T ghi)
{
    const T length = b - a;
    return illinois_root([&](T x) { return evaluate_1d(c, (x - a) / length); }, lo, hi, glo, ghi);
}

namespace
{
int sign_of(double v) { return (v > 0) - (v < 0); }

/// Recursion of isolate_roots on the piece [lo, hi] with Bernstein coefficients
/// stored at work[offset, offset + n + 1). The end coefficients of a piece are
/// the values at its ends, so neighbouring pieces agree on their shared end.
template <std::floating_point T>
void isolate_piece(std::span<const T> c, T a, T b, T lo, T hi, std::size_t offset, int depth,
                   std::vector<T>& roots, std::vector<T>& work)
{
    const int n1 = static_cast<int>(c.size());
    int changes = 0, last = 0;
    for (int i = 0; i < n1; ++i)
    {
        const int s = sign_of(work[offset + i]);
        if (s == 0)
            continue;
        if (last != 0 && s != last)
            ++changes;
        last = s;
    }
    if (changes == 0)
        return;
    // One sign change means one root inside; it is bracketed by the ends unless
    // one of them is a root itself, which needs more subdivision.
    const T glo = work[offset], ghi = work[offset + n1 - 1];
    const bool bracket = glo != T(0) && ghi != T(0);
    if ((changes == 1 && bracket) || depth >= 48)
    {
        if (bracket && (glo > T(0)) != (ghi > T(0)))
            roots.push_back(bracketed_root(c, a, b, lo, hi, glo, ghi));
        return;
    }
    // halve by de Casteljau: left at [offset + n1, ...), right at [offset + 2 n1, ...)
    const std::size_t left = offset + n1, right = offset + 2 * n1;
    if (work.size() < right + n1)
        work.resize(right + n1);
    for (int r = 0; r < n1; ++r)
    {
        work[left + r] = work[offset];
        work[right + n1 - 1 - r] = work[offset + n1 - 1 - r];
        for (int i = 0; i < n1 - 1 - r; ++i)
            work[offset + i] = T(0.5) * (work[offset + i] + work[offset + i + 1]);
    }
    const T mid = T(0.5) * (lo + hi);
    if (work[left + n1 - 1] == T(0))
    {
        // a root exactly at the split point: count it here if the sign changes
        // across it (the halves see a zero end and skip it)
        int before = 0, after = 0;
        for (int i = n1 - 2; i >= 0 && before == 0; --i)
            before = sign_of(work[left + i]);
        for (int i = 1; i < n1 && after == 0; ++i)
            after = sign_of(work[right + i]);
        if (before != 0 && after != 0 && before != after)
            roots.push_back(mid);
    }
    // the right half moves to the slot of the parent; the left half stays
    std::copy(work.begin() + right, work.begin() + right + n1, work.begin() + offset);
    isolate_piece(c, a, b, lo, mid, left, depth + 1, roots, work);
    isolate_piece(c, a, b, mid, hi, offset, depth + 1, roots, work);
}
} // namespace

template <std::floating_point T>
void isolate_roots(std::span<const T> c, T a, T b, std::vector<T>& roots, std::vector<T>& work)
{
    const int n1 = static_cast<int>(c.size());
    if (n1 < 2)
        return;
    work.resize(static_cast<std::size_t>(n1) * 3);
    std::copy(c.begin(), c.end(), work.begin());
    isolate_piece(c, a, b, a, b, 0, 0, roots, work);
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template void cell_bernstein_on_box<float>(cell::type, int, std::span<const float>, const ClippedBox<float>&,
                                          BoxBernstein<float>&);
template float evaluate<float>(const BoxBernstein<float>&, std::span<const float>);
template void gradient<float>(const BoxBernstein<float>&, std::span<const float>, std::span<float>);
template void derivative<float>(const BoxBernstein<float>&, int, BoxBernstein<float>&);
template void restrict_affine<float>(const BoxBernstein<float>&, std::span<const float>, std::span<const float>, int, BoxBernstein<float>&, std::vector<float>&);
template void subdivide<float>(const BoxBernstein<float>&, std::span<const float>, std::span<const float>, BoxBernstein<float>&, std::vector<float>&);
template bool may_vanish<float>(std::span<const float>);
template float max_abs<float>(std::span<const float>);
template float scaled_norm<float>(std::span<const float>);
template void margins<float>(const BoxBernstein<float>&, std::span<const float>, std::span<float>, BoxBernstein<float>&,
                             std::span<float>);
template float evaluate_1d<float>(std::span<const float>, float);
template float bracketed_root<float>(std::span<const float>, float, float, float, float, float, float);
template void isolate_roots<float>(std::span<const float>, float, float, std::vector<float>&, std::vector<float>&);

template void cell_bernstein_on_box<double>(cell::type, int, std::span<const double>, const ClippedBox<double>&,
                                           BoxBernstein<double>&);
template double evaluate<double>(const BoxBernstein<double>&, std::span<const double>);
template void gradient<double>(const BoxBernstein<double>&, std::span<const double>, std::span<double>);
template void derivative<double>(const BoxBernstein<double>&, int, BoxBernstein<double>&);
template void restrict_affine<double>(const BoxBernstein<double>&, std::span<const double>, std::span<const double>, int, BoxBernstein<double>&, std::vector<double>&);
template void subdivide<double>(const BoxBernstein<double>&, std::span<const double>, std::span<const double>, BoxBernstein<double>&, std::vector<double>&);
template bool may_vanish<double>(std::span<const double>);
template double max_abs<double>(std::span<const double>);
template double scaled_norm<double>(std::span<const double>);
template void margins<double>(const BoxBernstein<double>&, std::span<const double>, std::span<double>,
                              BoxBernstein<double>&, std::span<double>);
template double evaluate_1d<double>(std::span<const double>, double);
template double bracketed_root<double>(std::span<const double>, double, double, double, double, double, double);
template void isolate_roots<double>(std::span<const double>, double, double, std::vector<double>&, std::vector<double>&);

} // namespace cutcells::quadrays
