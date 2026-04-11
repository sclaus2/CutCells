// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "bernstein.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace cutcells::bernstein
{
namespace
{

// ---------------------------------------------------------------------------
// Arithmetic helpers
// ---------------------------------------------------------------------------

/// Compute base^exp for non-negative integer exponent.
/// Handles 0^0 = 1 correctly.
template <std::floating_point T>
T integer_pow(T base, int exp)
{
    T result = T(1);
    for (int i = 0; i < exp; ++i)
        result *= base;
    return result;
}

/// Binomial coefficient C(n, k) computed without overflow for moderate n.
template <std::floating_point T>
T binomial(int n, int k)
{
    if (k < 0 || k > n)
        return T(0);
    if (k == 0 || k == n)
        return T(1);
    if (k > n - k)
        k = n - k;
    T result = T(1);
    for (int i = 0; i < k; ++i)
    {
        result *= T(n - i);
        result /= T(i + 1);
    }
    return result;
}

/// Multinomial coefficient n! / (alpha[0]! * ... * alpha[d]!)
/// computed as a product of binomial coefficients.
/// alpha must satisfy sum(alpha) == n.
template <std::floating_point T>
T multinomial(int n, const int* alpha, int num_components)
{
    T result = T(1);
    int remaining = n;
    for (int c = 1; c < num_components; ++c)
    {
        result *= binomial<T>(remaining, alpha[c]);
        remaining -= alpha[c];
    }
    return result;
}

// ---------------------------------------------------------------------------
// Simplex multi-index → linear index
//
// Enumeration order (matching evaluate loops):
//   tdim=1: i = 0..n          → alpha = (n-i, i)
//   tdim=2: j = 0..n, i = 0..n-j   → alpha = (n-i-j, i, j)
//   tdim=3: k = 0..n, j = 0..n-k, i = 0..n-k-j → alpha = (n-i-j-k, i, j, k)
// ---------------------------------------------------------------------------

inline int simplex_index_1d(int i, int /*n*/)
{
    return i;
}

inline int simplex_index_2d(int i, int j, int n)
{
    return j * (n + 1) - j * (j - 1) / 2 + i;
}

inline int simplex_index_3d(int i, int j, int k, int n)
{
    // Offset for layers 0..k-1
    int offset = 0;
    for (int kk = 0; kk < k; ++kk)
    {
        int m = n - kk;
        offset += (m + 1) * (m + 2) / 2;
    }
    // Within layer k, it is a triangle of degree n-k
    return offset + simplex_index_2d(i, j, n - k);
}

// ---------------------------------------------------------------------------
// 1D Bernstein helpers (used by tensor-product evaluation)
// ---------------------------------------------------------------------------

/// Evaluate all (degree+1) 1D Bernstein basis functions at t.
template <std::floating_point T>
void bernstein_1d_all(int degree, T t, T* out)
{
    int n = degree;
    for (int i = 0; i <= n; ++i)
        out[i] = binomial<T>(n, i) * integer_pow(t, i) * integer_pow(T(1) - t, n - i);
}

/// Evaluate all degree 1D Bernstein basis functions of degree (n-1) at t,
/// and also return the derivative of each degree-n Bernstein basis function.
/// dB^n_i(t) = n * [B^{n-1}_{i-1}(t) - B^{n-1}_i(t)]
template <std::floating_point T>
void bernstein_1d_deriv(int degree, T t, T* dbasis)
{
    int n = degree;
    if (n == 0)
        return;

    // Evaluate degree-(n-1) basis
    std::vector<T> b(static_cast<std::size_t>(n)); // n = degree, n-1+1 = n entries
    bernstein_1d_all(n - 1, t, b.data());

    for (int i = 0; i <= n; ++i)
    {
        T bim1 = (i > 0) ? b[static_cast<std::size_t>(i - 1)] : T(0);
        T bi   = (i < n) ? b[static_cast<std::size_t>(i)]     : T(0);
        dbasis[i] = T(n) * (bim1 - bi);
    }
}

// ---------------------------------------------------------------------------
// Dense linear solver (Gaussian elimination with partial pivoting)
// Row-major storage: A[row * n + col]
// ---------------------------------------------------------------------------

template <std::floating_point T>
void solve_dense(int n, std::vector<T>& A, std::vector<T>& b)
{
    // Forward elimination
    for (int col = 0; col < n; ++col)
    {
        // Partial pivoting
        int max_row = col;
        T max_val = std::abs(A[static_cast<std::size_t>(col) * static_cast<std::size_t>(n)
                               + static_cast<std::size_t>(col)]);
        for (int row = col + 1; row < n; ++row)
        {
            T val = std::abs(A[static_cast<std::size_t>(row) * static_cast<std::size_t>(n)
                               + static_cast<std::size_t>(col)]);
            if (val > max_val)
            {
                max_val = val;
                max_row = row;
            }
        }
        if (max_val < T(1e-14))
            throw std::runtime_error("bernstein: singular matrix in Lagrange-to-Bernstein conversion");

        // Swap rows
        if (max_row != col)
        {
            for (int j = col; j < n; ++j)
                std::swap(A[static_cast<std::size_t>(col) * static_cast<std::size_t>(n)
                            + static_cast<std::size_t>(j)],
                          A[static_cast<std::size_t>(max_row) * static_cast<std::size_t>(n)
                            + static_cast<std::size_t>(j)]);
            std::swap(b[static_cast<std::size_t>(col)],
                      b[static_cast<std::size_t>(max_row)]);
        }

        // Eliminate below pivot
        for (int row = col + 1; row < n; ++row)
        {
            T factor = A[static_cast<std::size_t>(row) * static_cast<std::size_t>(n)
                         + static_cast<std::size_t>(col)]
                     / A[static_cast<std::size_t>(col) * static_cast<std::size_t>(n)
                         + static_cast<std::size_t>(col)];
            for (int j = col + 1; j < n; ++j)
            {
                A[static_cast<std::size_t>(row) * static_cast<std::size_t>(n)
                  + static_cast<std::size_t>(j)]
                    -= factor * A[static_cast<std::size_t>(col) * static_cast<std::size_t>(n)
                                  + static_cast<std::size_t>(j)];
            }
            b[static_cast<std::size_t>(row)] -= factor * b[static_cast<std::size_t>(col)];
        }
    }

    // Back substitution
    for (int row = n - 1; row >= 0; --row)
    {
        for (int j = row + 1; j < n; ++j)
        {
            b[static_cast<std::size_t>(row)]
                -= A[static_cast<std::size_t>(row) * static_cast<std::size_t>(n)
                     + static_cast<std::size_t>(j)]
                 * b[static_cast<std::size_t>(j)];
        }
        b[static_cast<std::size_t>(row)]
            /= A[static_cast<std::size_t>(row) * static_cast<std::size_t>(n)
                 + static_cast<std::size_t>(row)];
    }
}

// ---------------------------------------------------------------------------
// Simplex Bernstein evaluation
// ---------------------------------------------------------------------------

template <std::floating_point T>
T evaluate_simplex(int tdim, int degree, const T* coeffs, const T* xi)
{
    const int n = degree;

    // Barycentric coordinates: lambda[0] = 1 - sum(xi), lambda[k+1] = xi[k]
    T lambda[4];
    T sum = T(0);
    for (int d = 0; d < tdim; ++d)
        sum += xi[d];
    lambda[0] = T(1) - sum;
    for (int d = 0; d < tdim; ++d)
        lambda[d + 1] = xi[d];

    T result = T(0);
    int idx = 0;

    if (tdim == 1)
    {
        for (int i = 0; i <= n; ++i)
        {
            int alpha[2] = {n - i, i};
            T b = multinomial<T>(n, alpha, 2)
                * integer_pow(lambda[0], alpha[0])
                * integer_pow(lambda[1], alpha[1]);
            result += coeffs[idx++] * b;
        }
    }
    else if (tdim == 2)
    {
        for (int j = 0; j <= n; ++j)
        {
            for (int i = 0; i <= n - j; ++i)
            {
                int alpha[3] = {n - i - j, i, j};
                T b = multinomial<T>(n, alpha, 3)
                    * integer_pow(lambda[0], alpha[0])
                    * integer_pow(lambda[1], alpha[1])
                    * integer_pow(lambda[2], alpha[2]);
                result += coeffs[idx++] * b;
            }
        }
    }
    else if (tdim == 3)
    {
        for (int k = 0; k <= n; ++k)
        {
            for (int j = 0; j <= n - k; ++j)
            {
                for (int i = 0; i <= n - k - j; ++i)
                {
                    int alpha[4] = {n - i - j - k, i, j, k};
                    T b = multinomial<T>(n, alpha, 4)
                        * integer_pow(lambda[0], alpha[0])
                        * integer_pow(lambda[1], alpha[1])
                        * integer_pow(lambda[2], alpha[2])
                        * integer_pow(lambda[3], alpha[3]);
                    result += coeffs[idx++] * b;
                }
            }
        }
    }

    return result;
}

// ---------------------------------------------------------------------------
// Tensor-product Bernstein evaluation
// ---------------------------------------------------------------------------

template <std::floating_point T>
T evaluate_tensor(int tdim, int degree, const T* coeffs, const T* xi)
{
    const int n = degree;
    const int n1 = n + 1;

    if (tdim == 1)
    {
        // Same as simplex 1D
        T result = T(0);
        for (int i = 0; i <= n; ++i)
            result += coeffs[i] * binomial<T>(n, i)
                    * integer_pow(xi[0], i)
                    * integer_pow(T(1) - xi[0], n - i);
        return result;
    }
    else if (tdim == 2)
    {
        // Precompute 1D basis in each direction
        std::vector<T> Bx(static_cast<std::size_t>(n1)), By(static_cast<std::size_t>(n1));
        bernstein_1d_all(n, xi[0], Bx.data());
        bernstein_1d_all(n, xi[1], By.data());

        T result = T(0);
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n; ++j)
                result += coeffs[i * n1 + j] * Bx[static_cast<std::size_t>(i)]
                        * By[static_cast<std::size_t>(j)];
        return result;
    }
    else if (tdim == 3)
    {
        std::vector<T> Bx(static_cast<std::size_t>(n1)),
                        By(static_cast<std::size_t>(n1)),
                        Bz(static_cast<std::size_t>(n1));
        bernstein_1d_all(n, xi[0], Bx.data());
        bernstein_1d_all(n, xi[1], By.data());
        bernstein_1d_all(n, xi[2], Bz.data());

        T result = T(0);
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n; ++j)
                for (int k = 0; k <= n; ++k)
                    result += coeffs[(i * n1 + j) * n1 + k]
                            * Bx[static_cast<std::size_t>(i)]
                            * By[static_cast<std::size_t>(j)]
                            * Bz[static_cast<std::size_t>(k)];
        return result;
    }

    return T(0);
}

// ---------------------------------------------------------------------------
// Simplex Bernstein gradient
//
// Uses the derivative formula:
//   df/d(xi_k) = n * sum_{|beta|=n-1} (c[beta+e_{k+1}] - c[beta+e_0]) * B^{n-1}_beta(xi)
//
// where e_0 corresponds to lambda_0 and e_{k+1} to lambda_{k+1} = xi_k.
// ---------------------------------------------------------------------------

template <std::floating_point T>
void gradient_simplex(int tdim, int degree, const T* coeffs, const T* xi, T* grad)
{
    const int n = degree;

    for (int d = 0; d < tdim; ++d)
        grad[d] = T(0);

    if (n == 0)
        return;

    // Barycentric coordinates
    T lambda[4];
    T sum = T(0);
    for (int d = 0; d < tdim; ++d)
        sum += xi[d];
    lambda[0] = T(1) - sum;
    for (int d = 0; d < tdim; ++d)
        lambda[d + 1] = xi[d];

    const int nm1 = n - 1;

    if (tdim == 1)
    {
        // df/dxi = n * sum_{i=0}^{n-1} (c[i+1] - c[i]) * B^{n-1}_i(xi)
        T g = T(0);
        for (int i = 0; i <= nm1; ++i)
        {
            int alpha[2] = {nm1 - i, i};
            T b = multinomial<T>(nm1, alpha, 2)
                * integer_pow(lambda[0], alpha[0])
                * integer_pow(lambda[1], alpha[1]);
            T dc = coeffs[i + 1] - coeffs[i];
            g += dc * b;
        }
        grad[0] = T(n) * g;
    }
    else if (tdim == 2)
    {
        // For direction k (0-indexed reference direction):
        //   df/d(xi_k) = n * sum_{beta with |beta|=n-1}
        //                    (c[beta + e_{k+1}] - c[beta + e_0]) * B^{n-1}_beta(xi)
        //
        // beta = (nm1-i-j, i, j), enumerated by (j, i)
        // beta + e_0 has index in degree-n array: simplex_index_2d(i, j, n)
        // beta + e_1 has index: simplex_index_2d(i+1, j, n)
        // beta + e_2 has index: simplex_index_2d(i, j+1, n)

        T g0 = T(0), g1 = T(0);
        for (int j = 0; j <= nm1; ++j)
        {
            for (int i = 0; i <= nm1 - j; ++i)
            {
                int alpha[3] = {nm1 - i - j, i, j};
                T b = multinomial<T>(nm1, alpha, 3)
                    * integer_pow(lambda[0], alpha[0])
                    * integer_pow(lambda[1], alpha[1])
                    * integer_pow(lambda[2], alpha[2]);

                int idx_e0 = simplex_index_2d(i, j, n);
                int idx_e1 = simplex_index_2d(i + 1, j, n);
                int idx_e2 = simplex_index_2d(i, j + 1, n);

                g0 += (coeffs[idx_e1] - coeffs[idx_e0]) * b;
                g1 += (coeffs[idx_e2] - coeffs[idx_e0]) * b;
            }
        }
        grad[0] = T(n) * g0;
        grad[1] = T(n) * g1;
    }
    else if (tdim == 3)
    {
        T g0 = T(0), g1 = T(0), g2 = T(0);
        for (int k = 0; k <= nm1; ++k)
        {
            for (int j = 0; j <= nm1 - k; ++j)
            {
                for (int i = 0; i <= nm1 - k - j; ++i)
                {
                    int alpha[4] = {nm1 - i - j - k, i, j, k};
                    T b = multinomial<T>(nm1, alpha, 4)
                        * integer_pow(lambda[0], alpha[0])
                        * integer_pow(lambda[1], alpha[1])
                        * integer_pow(lambda[2], alpha[2])
                        * integer_pow(lambda[3], alpha[3]);

                    int idx_e0 = simplex_index_3d(i, j, k, n);
                    int idx_e1 = simplex_index_3d(i + 1, j, k, n);
                    int idx_e2 = simplex_index_3d(i, j + 1, k, n);
                    int idx_e3 = simplex_index_3d(i, j, k + 1, n);

                    g0 += (coeffs[idx_e1] - coeffs[idx_e0]) * b;
                    g1 += (coeffs[idx_e2] - coeffs[idx_e0]) * b;
                    g2 += (coeffs[idx_e3] - coeffs[idx_e0]) * b;
                }
            }
        }
        grad[0] = T(n) * g0;
        grad[1] = T(n) * g1;
        grad[2] = T(n) * g2;
    }
}

// ---------------------------------------------------------------------------
// Tensor-product Bernstein gradient
// ---------------------------------------------------------------------------

template <std::floating_point T>
void gradient_tensor(int tdim, int degree, const T* coeffs, const T* xi, T* grad)
{
    const int n = degree;
    const int n1 = n + 1;

    for (int d = 0; d < tdim; ++d)
        grad[d] = T(0);

    if (tdim == 1)
    {
        // Same as simplex 1D gradient
        if (n == 0)
            return;
        std::vector<T> bm1(static_cast<std::size_t>(n));
        bernstein_1d_all(n - 1, xi[0], bm1.data());

        T g = T(0);
        for (int i = 0; i <= n; ++i)
        {
            T bim1 = (i > 0) ? bm1[static_cast<std::size_t>(i - 1)] : T(0);
            T bi   = (i < n) ? bm1[static_cast<std::size_t>(i)]     : T(0);
            g += coeffs[i] * T(n) * (bim1 - bi);
        }
        grad[0] = g;
    }
    else if (tdim == 2)
    {
        std::vector<T> Bx(static_cast<std::size_t>(n1)),
                        By(static_cast<std::size_t>(n1));
        std::vector<T> dBx(static_cast<std::size_t>(n1)),
                        dBy(static_cast<std::size_t>(n1));
        bernstein_1d_all(n, xi[0], Bx.data());
        bernstein_1d_all(n, xi[1], By.data());
        bernstein_1d_deriv(n, xi[0], dBx.data());
        bernstein_1d_deriv(n, xi[1], dBy.data());

        T g0 = T(0), g1 = T(0);
        for (int i = 0; i <= n; ++i)
        {
            for (int j = 0; j <= n; ++j)
            {
                T c = coeffs[i * n1 + j];
                g0 += c * dBx[static_cast<std::size_t>(i)]
                        * By[static_cast<std::size_t>(j)];
                g1 += c * Bx[static_cast<std::size_t>(i)]
                        * dBy[static_cast<std::size_t>(j)];
            }
        }
        grad[0] = g0;
        grad[1] = g1;
    }
    else if (tdim == 3)
    {
        std::vector<T> Bx(static_cast<std::size_t>(n1)),
                        By(static_cast<std::size_t>(n1)),
                        Bz(static_cast<std::size_t>(n1));
        std::vector<T> dBx(static_cast<std::size_t>(n1)),
                        dBy(static_cast<std::size_t>(n1)),
                        dBz(static_cast<std::size_t>(n1));
        bernstein_1d_all(n, xi[0], Bx.data());
        bernstein_1d_all(n, xi[1], By.data());
        bernstein_1d_all(n, xi[2], Bz.data());
        bernstein_1d_deriv(n, xi[0], dBx.data());
        bernstein_1d_deriv(n, xi[1], dBy.data());
        bernstein_1d_deriv(n, xi[2], dBz.data());

        T g0 = T(0), g1 = T(0), g2 = T(0);
        for (int i = 0; i <= n; ++i)
        {
            for (int j = 0; j <= n; ++j)
            {
                for (int k = 0; k <= n; ++k)
                {
                    T c = coeffs[(i * n1 + j) * n1 + k];
                    T bx = Bx[static_cast<std::size_t>(i)];
                    T by = By[static_cast<std::size_t>(j)];
                    T bz = Bz[static_cast<std::size_t>(k)];
                    g0 += c * dBx[static_cast<std::size_t>(i)] * by * bz;
                    g1 += c * bx * dBy[static_cast<std::size_t>(j)] * bz;
                    g2 += c * bx * by * dBz[static_cast<std::size_t>(k)];
                }
            }
        }
        grad[0] = g0;
        grad[1] = g1;
        grad[2] = g2;
    }
}

// ---------------------------------------------------------------------------
// Simplex evaluate_all_basis
// ---------------------------------------------------------------------------

template <std::floating_point T>
void evaluate_all_simplex(int tdim, int degree, const T* xi, T* out)
{
    const int n = degree;

    T lambda[4];
    T sum = T(0);
    for (int d = 0; d < tdim; ++d)
        sum += xi[d];
    lambda[0] = T(1) - sum;
    for (int d = 0; d < tdim; ++d)
        lambda[d + 1] = xi[d];

    int idx = 0;

    if (tdim == 1)
    {
        for (int i = 0; i <= n; ++i)
        {
            int alpha[2] = {n - i, i};
            out[idx++] = multinomial<T>(n, alpha, 2)
                       * integer_pow(lambda[0], alpha[0])
                       * integer_pow(lambda[1], alpha[1]);
        }
    }
    else if (tdim == 2)
    {
        for (int j = 0; j <= n; ++j)
        {
            for (int i = 0; i <= n - j; ++i)
            {
                int alpha[3] = {n - i - j, i, j};
                out[idx++] = multinomial<T>(n, alpha, 3)
                           * integer_pow(lambda[0], alpha[0])
                           * integer_pow(lambda[1], alpha[1])
                           * integer_pow(lambda[2], alpha[2]);
            }
        }
    }
    else if (tdim == 3)
    {
        for (int k = 0; k <= n; ++k)
        {
            for (int j = 0; j <= n - k; ++j)
            {
                for (int i = 0; i <= n - k - j; ++i)
                {
                    int alpha[4] = {n - i - j - k, i, j, k};
                    out[idx++] = multinomial<T>(n, alpha, 4)
                               * integer_pow(lambda[0], alpha[0])
                               * integer_pow(lambda[1], alpha[1])
                               * integer_pow(lambda[2], alpha[2])
                               * integer_pow(lambda[3], alpha[3]);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Tensor-product evaluate_all_basis
// ---------------------------------------------------------------------------

template <std::floating_point T>
void evaluate_all_tensor(int tdim, int degree, const T* xi, T* out)
{
    const int n = degree;
    const int n1 = n + 1;

    if (tdim == 1)
    {
        bernstein_1d_all(n, xi[0], out);
    }
    else if (tdim == 2)
    {
        std::vector<T> Bx(static_cast<std::size_t>(n1)),
                        By(static_cast<std::size_t>(n1));
        bernstein_1d_all(n, xi[0], Bx.data());
        bernstein_1d_all(n, xi[1], By.data());

        int idx = 0;
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n; ++j)
                out[idx++] = Bx[static_cast<std::size_t>(i)]
                           * By[static_cast<std::size_t>(j)];
    }
    else if (tdim == 3)
    {
        std::vector<T> Bx(static_cast<std::size_t>(n1)),
                        By(static_cast<std::size_t>(n1)),
                        Bz(static_cast<std::size_t>(n1));
        bernstein_1d_all(n, xi[0], Bx.data());
        bernstein_1d_all(n, xi[1], By.data());
        bernstein_1d_all(n, xi[2], Bz.data());

        int idx = 0;
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n; ++j)
                for (int k = 0; k <= n; ++k)
                    out[idx++] = Bx[static_cast<std::size_t>(i)]
                               * By[static_cast<std::size_t>(j)]
                               * Bz[static_cast<std::size_t>(k)];
    }
}

} // anonymous namespace

// ===========================================================================
// Public API
// ===========================================================================

int num_polynomials(cell::type ctype, int degree)
{
    const int n = degree;
    switch (ctype)
    {
        case cell::type::interval:
            return n + 1;
        case cell::type::triangle:
            return (n + 1) * (n + 2) / 2;
        case cell::type::tetrahedron:
            return (n + 1) * (n + 2) * (n + 3) / 6;
        case cell::type::quadrilateral:
            return (n + 1) * (n + 1);
        case cell::type::hexahedron:
            return (n + 1) * (n + 1) * (n + 1);
        default:
            throw std::invalid_argument(
                "bernstein::num_polynomials: unsupported cell type");
    }
}

// --- evaluate ---

template <std::floating_point T>
T evaluate(cell::type ctype, int degree,
           std::span<const T> coeffs, std::span<const T> xi)
{
    if (is_simplex(ctype))
        return evaluate_simplex(cell::get_tdim(ctype), degree,
                                coeffs.data(), xi.data());
    else if (is_tensor_product(ctype))
        return evaluate_tensor(cell::get_tdim(ctype), degree,
                               coeffs.data(), xi.data());
    else
        throw std::invalid_argument("bernstein::evaluate: unsupported cell type");
}

template double evaluate(cell::type, int, std::span<const double>, std::span<const double>);
template float  evaluate(cell::type, int, std::span<const float>,  std::span<const float>);

// --- gradient ---

template <std::floating_point T>
void gradient(cell::type ctype, int degree,
              std::span<const T> coeffs, std::span<const T> xi,
              std::span<T> grad)
{
    if (is_simplex(ctype))
        gradient_simplex(cell::get_tdim(ctype), degree,
                         coeffs.data(), xi.data(), grad.data());
    else if (is_tensor_product(ctype))
        gradient_tensor(cell::get_tdim(ctype), degree,
                        coeffs.data(), xi.data(), grad.data());
    else
        throw std::invalid_argument("bernstein::gradient: unsupported cell type");
}

template void gradient(cell::type, int, std::span<const double>, std::span<const double>, std::span<double>);
template void gradient(cell::type, int, std::span<const float>,  std::span<const float>,  std::span<float>);

// --- evaluate_basis ---

template <std::floating_point T>
void evaluate_basis(cell::type ctype, int degree,
                    std::span<const T> xi, std::span<T> out)
{
    if (is_simplex(ctype))
        evaluate_all_simplex(cell::get_tdim(ctype), degree,
                             xi.data(), out.data());
    else if (is_tensor_product(ctype))
        evaluate_all_tensor(cell::get_tdim(ctype), degree,
                            xi.data(), out.data());
    else
        throw std::invalid_argument("bernstein::evaluate_basis: unsupported cell type");
}

template void evaluate_basis(cell::type, int, std::span<const double>, std::span<double>);
template void evaluate_basis(cell::type, int, std::span<const float>,  std::span<float>);

// --- lagrange_to_bernstein ---

template <std::floating_point T>
void lagrange_to_bernstein(cell::type ctype, int degree,
                           std::span<const T> ref_points,
                           std::span<const T> nodal_values,
                           std::vector<T>& coeffs)
{
    const int N = num_polynomials(ctype, degree);
    const int tdim = cell::get_tdim(ctype);
    const int ndofs = static_cast<int>(nodal_values.size());

    if (ndofs != N)
        throw std::invalid_argument(
            "bernstein::lagrange_to_bernstein: nodal_values size ("
            + std::to_string(ndofs)
            + ") does not match num_polynomials ("
            + std::to_string(N) + ")");

    if (static_cast<int>(ref_points.size()) != ndofs * tdim)
        throw std::invalid_argument(
            "bernstein::lagrange_to_bernstein: ref_points size mismatch");

    // Build the Bernstein evaluation matrix V:  V[i][j] = B_j(ref_point_i)
    std::vector<T> V(static_cast<std::size_t>(N) * static_cast<std::size_t>(N));
    std::vector<T> basis_row(static_cast<std::size_t>(N));

    for (int i = 0; i < N; ++i)
    {
        const T* pt = ref_points.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(tdim);
        std::span<const T> pt_span(pt, static_cast<std::size_t>(tdim));
        std::span<T> row_span(basis_row);

        evaluate_basis(ctype, degree, pt_span, row_span);

        for (int j = 0; j < N; ++j)
            V[static_cast<std::size_t>(i) * static_cast<std::size_t>(N)
              + static_cast<std::size_t>(j)] = basis_row[static_cast<std::size_t>(j)];
    }

    // Solve V * coeffs = nodal_values
    coeffs.assign(nodal_values.begin(), nodal_values.end());
    solve_dense(N, V, coeffs);
}

template void lagrange_to_bernstein(cell::type, int,
                                    std::span<const double>, std::span<const double>,
                                    std::vector<double>&);
template void lagrange_to_bernstein(cell::type, int,
                                    std::span<const float>, std::span<const float>,
                                    std::vector<float>&);

// --- derivative_coefficients ---

template <std::floating_point T>
void derivative_coefficients(cell::type ctype, int degree,
                             std::span<const T> coeffs,
                             int direction,
                             std::vector<T>& deriv_coeffs)
{
    const int n = degree;
    if (n <= 0)
    {
        deriv_coeffs.clear();
        return;
    }

    const int tdim = cell::get_tdim(ctype);
    if (direction < 0 || direction >= tdim)
        throw std::invalid_argument(
            "derivative_coefficients: direction out of range");

    if (is_simplex(ctype))
    {
        // Derivative of degree-n simplex Bernstein → degree-(n-1) expansion.
        // c'_beta = n * (c_{beta + e_{dir+1}} - c_{beta + e_0})
        // where beta ranges over multi-indices with |beta| = n-1.
        const int nm1 = n - 1;
        const int n_out = num_polynomials(ctype, nm1);
        deriv_coeffs.resize(static_cast<std::size_t>(n_out));

        if (tdim == 1)
        {
            // df/dx: c'_i = n * (c_{i+1} - c_i), i = 0..n-1
            for (int i = 0; i <= nm1; ++i)
                deriv_coeffs[static_cast<std::size_t>(i)] =
                    T(n) * (coeffs[static_cast<std::size_t>(i + 1)]
                          - coeffs[static_cast<std::size_t>(i)]);
        }
        else if (tdim == 2)
        {
            int idx = 0;
            for (int j = 0; j <= nm1; ++j)
            {
                for (int i = 0; i <= nm1 - j; ++i, ++idx)
                {
                    int idx_e0 = simplex_index_2d(i, j, n);
                    int idx_ek;
                    if (direction == 0)
                        idx_ek = simplex_index_2d(i + 1, j, n);
                    else
                        idx_ek = simplex_index_2d(i, j + 1, n);

                    deriv_coeffs[static_cast<std::size_t>(idx)] =
                        T(n) * (coeffs[static_cast<std::size_t>(idx_ek)]
                              - coeffs[static_cast<std::size_t>(idx_e0)]);
                }
            }
        }
        else if (tdim == 3)
        {
            int idx = 0;
            for (int k = 0; k <= nm1; ++k)
            {
                for (int j = 0; j <= nm1 - k; ++j)
                {
                    for (int i = 0; i <= nm1 - k - j; ++i, ++idx)
                    {
                        int idx_e0 = simplex_index_3d(i, j, k, n);
                        int idx_ek;
                        if (direction == 0)
                            idx_ek = simplex_index_3d(i + 1, j, k, n);
                        else if (direction == 1)
                            idx_ek = simplex_index_3d(i, j + 1, k, n);
                        else
                            idx_ek = simplex_index_3d(i, j, k + 1, n);

                        deriv_coeffs[static_cast<std::size_t>(idx)] =
                            T(n) * (coeffs[static_cast<std::size_t>(idx_ek)]
                                  - coeffs[static_cast<std::size_t>(idx_e0)]);
                    }
                }
            }
        }
    }
    else if (is_tensor_product(ctype))
    {
        // Tensor-product: derivative in direction k reduces degree from n to n-1
        // in that direction; degree stays n in other directions.
        const int nm1 = n - 1;
        const int n1 = n + 1;

        if (tdim == 1)
        {
            deriv_coeffs.resize(static_cast<std::size_t>(n));
            for (int i = 0; i < n; ++i)
                deriv_coeffs[static_cast<std::size_t>(i)] =
                    T(n) * (coeffs[static_cast<std::size_t>(i + 1)]
                          - coeffs[static_cast<std::size_t>(i)]);
        }
        else if (tdim == 2)
        {
            // Index: coeffs[ix * n1 + iy]
            if (direction == 0)
            {
                // d/dx: degree (n-1) in x, n in y → n * (n+1) entries
                deriv_coeffs.resize(static_cast<std::size_t>(n * n1));
                for (int ix = 0; ix < n; ++ix)
                    for (int iy = 0; iy <= n; ++iy)
                        deriv_coeffs[static_cast<std::size_t>(ix * n1 + iy)] =
                            T(n) * (coeffs[static_cast<std::size_t>((ix + 1) * n1 + iy)]
                                  - coeffs[static_cast<std::size_t>(ix * n1 + iy)]);
            }
            else
            {
                // d/dy: degree n in x, (n-1) in y → (n+1) * n entries
                deriv_coeffs.resize(static_cast<std::size_t>(n1 * n));
                for (int ix = 0; ix <= n; ++ix)
                    for (int iy = 0; iy < n; ++iy)
                        deriv_coeffs[static_cast<std::size_t>(ix * n + iy)] =
                            T(n) * (coeffs[static_cast<std::size_t>(ix * n1 + iy + 1)]
                                  - coeffs[static_cast<std::size_t>(ix * n1 + iy)]);
            }
        }
        else if (tdim == 3)
        {
            // Index: coeffs[(ix * n1 + iy) * n1 + iz]
            if (direction == 0)
            {
                deriv_coeffs.resize(static_cast<std::size_t>(n * n1 * n1));
                for (int ix = 0; ix < n; ++ix)
                    for (int iy = 0; iy <= n; ++iy)
                        for (int iz = 0; iz <= n; ++iz)
                            deriv_coeffs[static_cast<std::size_t>((ix * n1 + iy) * n1 + iz)] =
                                T(n) * (coeffs[static_cast<std::size_t>(((ix + 1) * n1 + iy) * n1 + iz)]
                                      - coeffs[static_cast<std::size_t>((ix * n1 + iy) * n1 + iz)]);
            }
            else if (direction == 1)
            {
                deriv_coeffs.resize(static_cast<std::size_t>(n1 * n * n1));
                for (int ix = 0; ix <= n; ++ix)
                    for (int iy = 0; iy < n; ++iy)
                        for (int iz = 0; iz <= n; ++iz)
                            deriv_coeffs[static_cast<std::size_t>((ix * n + iy) * n1 + iz)] =
                                T(n) * (coeffs[static_cast<std::size_t>((ix * n1 + iy + 1) * n1 + iz)]
                                      - coeffs[static_cast<std::size_t>((ix * n1 + iy) * n1 + iz)]);
            }
            else
            {
                deriv_coeffs.resize(static_cast<std::size_t>(n1 * n1 * n));
                for (int ix = 0; ix <= n; ++ix)
                    for (int iy = 0; iy <= n; ++iy)
                        for (int iz = 0; iz < n; ++iz)
                            deriv_coeffs[static_cast<std::size_t>((ix * n1 + iy) * n + iz)] =
                                T(n) * (coeffs[static_cast<std::size_t>((ix * n1 + iy) * n1 + iz + 1)]
                                      - coeffs[static_cast<std::size_t>((ix * n1 + iy) * n1 + iz)]);
            }
        }
    }
    else
    {
        throw std::invalid_argument(
            "derivative_coefficients: unsupported cell type");
    }
}

template void derivative_coefficients(cell::type, int,
                                      std::span<const double>, int,
                                      std::vector<double>&);
template void derivative_coefficients(cell::type, int,
                                      std::span<const float>, int,
                                      std::vector<float>&);

} // namespace cutcells::bernstein
