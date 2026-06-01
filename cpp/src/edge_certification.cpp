// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "edge_certification.h"
#include "bernstein.h"
#include "cell_topology.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cutcells
{
namespace
{

template <std::floating_point T>
T integer_pow(T base, int exp)
{
    T result = T(1);
    for (int i = 0; i < exp; ++i)
        result *= base;
    return result;
}

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

template <std::floating_point T>
std::vector<T> linear_power_poly(T a, T b, int n)
{
    std::vector<T> out(static_cast<std::size_t>(n + 1), T(0));
    for (int i = 0; i <= n; ++i)
    {
        out[static_cast<std::size_t>(i)] =
            binomial<T>(n, i) * integer_pow(a, n - i) * integer_pow(b, i);
    }
    return out;
}

template <std::floating_point T>
std::vector<T> multiply_poly(std::span<const T> a, std::span<const T> b)
{
    std::vector<T> out(a.size() + b.size() - 1, T(0));
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j)
            out[i + j] += a[i] * b[j];
    return out;
}

template <std::floating_point T>
void add_scaled_poly(std::vector<T>& dst, std::span<const T> src, T scale)
{
    if (dst.size() < src.size())
        dst.resize(src.size(), T(0));
    for (std::size_t i = 0; i < src.size(); ++i)
        dst[i] += scale * src[i];
}

template <std::floating_point T>
void power_to_bernstein_1d(std::span<const T> power_coeffs,
                           std::vector<T>& bernstein_coeffs)
{
    const int degree = static_cast<int>(power_coeffs.size()) - 1;
    bernstein_coeffs.assign(static_cast<std::size_t>(degree + 1), T(0));

    for (int i = 0; i <= degree; ++i)
    {
        T coeff = T(0);
        for (int m = 0; m <= i; ++m)
        {
            coeff += power_coeffs[static_cast<std::size_t>(m)]
                   * binomial<T>(i, m)
                   / binomial<T>(degree, m);
        }
        bernstein_coeffs[static_cast<std::size_t>(i)] = coeff;
    }
}

// =====================================================================
// Simplex multi-index helpers (same ordering as bernstein.cpp)
// =====================================================================

inline int simplex_index_1d(int i, int /*n*/) { return i; }

inline int simplex_index_2d(int i, int j, int n)
{
    return j * (n + 1) - j * (j - 1) / 2 + i;
}

inline int simplex_index_3d(int i, int j, int k, int n)
{
    int offset = 0;
    for (int kk = 0; kk < k; ++kk)
    {
        int m = n - kk;
        offset += (m + 1) * (m + 2) / 2;
    }
    return offset + simplex_index_2d(i, j, n - k);
}

// =====================================================================
// Tensor-product index helpers
// =====================================================================

inline int tp_index_2d(int ix, int iy, int n)
{
    return iy * (n + 1) + ix;
}

inline int tp_index_3d(int ix, int iy, int iz, int n)
{
    return iz * (n + 1) * (n + 1) + iy * (n + 1) + ix;
}

// =====================================================================
// Edge-to-multi-index mapping for simplex parent edge extraction
// =====================================================================

/// For a triangle edge with local vertices (v0, v1) in the canonical numbering
/// {0, 1, 2}, the Bernstein multi-index alpha = (alpha_0, alpha_1, alpha_2)
/// with sum = degree, varies only along the edge.
/// We parameterize: alpha[v0] = n-k, alpha[v1] = k, all others = 0.
///
/// Returns the linear Bernstein index for the k-th point on that edge.
inline int triangle_edge_coeff_index(int v0, int v1, int k, int n)
{
    // multi-index: alpha[v0] = n-k, alpha[v1] = k, alpha[other] = 0
    int alpha[3] = {0, 0, 0};
    alpha[v0] = n - k;
    alpha[v1] = k;
    // simplex_index_2d expects (i, j, n) where alpha = (n-i-j, i, j)
    return simplex_index_2d(alpha[1], alpha[2], n);
}

/// Same logic for tetrahedron.
inline int tetrahedron_edge_coeff_index(int v0, int v1, int k, int n)
{
    int alpha[4] = {0, 0, 0, 0};
    alpha[v0] = n - k;
    alpha[v1] = k;
    // simplex_index_3d expects (i, j, k_idx, n) where alpha = (n-i-j-k_idx, i, j, k_idx)
    return simplex_index_3d(alpha[1], alpha[2], alpha[3], n);
}

} // anonymous namespace

// =====================================================================
// subdivide_bernstein_1d
// =====================================================================

template <std::floating_point T>
void subdivide_bernstein_1d(std::span<const T> coeffs,
                            T t_split,
                            std::vector<T>& left,
                            std::vector<T>& right)
{
    const int p = static_cast<int>(coeffs.size()) - 1;
    assert(p >= 0);

    left.resize(static_cast<std::size_t>(p + 1));
    right.resize(static_cast<std::size_t>(p + 1));

    // Work array: de Casteljau triangle columns.
    // Starting column = coeffs.
    std::vector<T> work(coeffs.begin(), coeffs.end());

    // The left child's i-th coeff is work[0] after i reduction steps.
    // The right child's i-th coeff is work[p-i] after (p-i) reduction steps.
    left[0] = work[0];
    right[static_cast<std::size_t>(p)] = work[static_cast<std::size_t>(p)];

    for (int level = 1; level <= p; ++level)
    {
        // Reduce: work[j] = (1-t)*work[j] + t*work[j+1]  for j = 0..p-level
        for (int j = 0; j <= p - level; ++j)
        {
            work[static_cast<std::size_t>(j)] =
                (T(1) - t_split) * work[static_cast<std::size_t>(j)]
                + t_split * work[static_cast<std::size_t>(j + 1)];
        }
        left[static_cast<std::size_t>(level)] = work[0];
        right[static_cast<std::size_t>(p - level)] =
            work[static_cast<std::size_t>(p - level)];
    }
}

// =====================================================================
// Sign-hull helpers
// =====================================================================

template <std::floating_point T>
bool bernstein_all_positive(std::span<const T> coeffs, T tol)
{
    for (const auto& c : coeffs)
        if (c <= tol) return false;
    return true;
}

template <std::floating_point T>
bool bernstein_all_negative(std::span<const T> coeffs, T tol)
{
    for (const auto& c : coeffs)
        if (c >= -tol) return false;
    return true;
}

template <std::floating_point T>
bool bernstein_all_zero(std::span<const T> coeffs, T tol)
{
    for (const auto& c : coeffs)
        if (std::fabs(c) > tol) return false;
    return true;
}

// =====================================================================
// find_root_intervals_1d
// =====================================================================

template <std::floating_point T>
void find_root_intervals_1d(std::span<const T> coeffs,
                            T t0, T t1,
                            T zero_tol, T sign_tol,
                            int depth, int max_depth,
                            std::vector<EdgeRootInterval<T>>& intervals,
                            bool& has_zero_segment)
{
    // If convex hull excludes zero → no root in this interval.
    if (bernstein_all_positive(coeffs, sign_tol)
        || bernstein_all_negative(coeffs, sign_tol))
        return;

    // If all coefficients are zero → zero segment.
    if (bernstein_all_zero(coeffs, zero_tol))
    {
        has_zero_segment = true;
        return;
    }

    // If max depth reached, record this interval as containing a root.
    if (depth >= max_depth)
    {
        intervals.push_back({t0, t1});
        return;
    }

    // Subdivide at midpoint.
    T t_mid = (t0 + t1) * T(0.5);
    std::vector<T> left_c, right_c;
    subdivide_bernstein_1d(coeffs, T(0.5), left_c, right_c);

    find_root_intervals_1d<T>(
        std::span<const T>(left_c), t0, t_mid,
        zero_tol, sign_tol, depth + 1, max_depth, intervals, has_zero_segment);

    find_root_intervals_1d<T>(
        std::span<const T>(right_c), t_mid, t1,
        zero_tol, sign_tol, depth + 1, max_depth, intervals, has_zero_segment);
}

// =====================================================================
// Merge overlapping root intervals
// =====================================================================

namespace
{

template <std::floating_point T>
void merge_root_intervals(std::vector<EdgeRootInterval<T>>& intervals,
                          T merge_tol)
{
    if (intervals.size() <= 1) return;

    // Sort by t0.
    std::sort(intervals.begin(), intervals.end(),
              [](const EdgeRootInterval<T>& a, const EdgeRootInterval<T>& b)
              { return a.t0 < b.t0; });

    std::vector<EdgeRootInterval<T>> merged;
    merged.push_back(intervals[0]);

    for (std::size_t i = 1; i < intervals.size(); ++i)
    {
        auto& back = merged.back();
        if (intervals[i].t0 <= back.t1 + merge_tol)
        {
            // Overlapping or adjacent — extend.
            back.t1 = std::max(back.t1, intervals[i].t1);
        }
        else
        {
            merged.push_back(intervals[i]);
        }
    }

    intervals = std::move(merged);
}

template <std::floating_point T>
bool linear_one_root_parameter_from_endpoint_values(T phi0,
                                                    T phi1,
                                                    T zero_tol,
                                                    T& root_t)
{
    const bool left_zero = std::fabs(phi0) <= zero_tol;
    const bool right_zero = std::fabs(phi1) <= zero_tol;

    root_t = T(0);
    if (left_zero && right_zero)
        return false;
    if (left_zero)
        return true;
    if (right_zero)
    {
        root_t = T(1);
        return true;
    }

    const bool brackets_zero = (phi0 < T(0) && phi1 > T(0))
                            || (phi0 > T(0) && phi1 < T(0));
    if (!brackets_zero)
        return false;

    const T denom = phi0 - phi1;
    const T denom_tol =
        std::max(zero_tol, T(64) * std::numeric_limits<T>::epsilon());
    if (std::fabs(denom) <= denom_tol)
        return false;

    root_t = std::clamp(phi0 / denom, T(0), T(1));
    return true;
}

} // anonymous namespace

// =====================================================================
// extract_parent_edge_bernstein
// =====================================================================

template <std::floating_point T>
void extract_parent_edge_bernstein(cell::type parent_cell_type,
                                   int degree,
                                   std::span<const T> parent_coeffs,
                                   int parent_local_edge_id,
                                   std::vector<T>& edge_coeffs)
{
    const int p = degree;
    edge_coeffs.resize(static_cast<std::size_t>(p + 1));

    auto edge_verts = cell::edges(parent_cell_type);
    const int v0 = edge_verts[static_cast<std::size_t>(parent_local_edge_id)][0];
    const int v1 = edge_verts[static_cast<std::size_t>(parent_local_edge_id)][1];

    switch (parent_cell_type)
    {
    case cell::type::interval:
        // Trivial: the edge IS the cell.
        for (int k = 0; k <= p; ++k)
            edge_coeffs[static_cast<std::size_t>(k)] =
                parent_coeffs[static_cast<std::size_t>(k)];
        break;

    case cell::type::triangle:
        for (int k = 0; k <= p; ++k)
            edge_coeffs[static_cast<std::size_t>(k)] =
                parent_coeffs[static_cast<std::size_t>(
                    triangle_edge_coeff_index(v0, v1, k, p))];
        break;

    case cell::type::tetrahedron:
        for (int k = 0; k <= p; ++k)
            edge_coeffs[static_cast<std::size_t>(k)] =
                parent_coeffs[static_cast<std::size_t>(
                    tetrahedron_edge_coeff_index(v0, v1, k, p))];
        break;

    case cell::type::quadrilateral:
    {
        // Basix quad vertices: 0=(0,0),1=(1,0),2=(0,1),3=(1,1)
        // Basix quad edges: 0={0,1}, 1={0,2}, 2={1,3}, 3={2,3}
        // Bernstein tensor-product index: tp_index_2d(ix,iy,p) = iy*(p+1)+ix
        for (int k = 0; k <= p; ++k)
        {
            int ix = 0, iy = 0;
            switch (parent_local_edge_id)
            {
            case 0: ix = k;     iy = 0;     break;  // {0,1}: ix varies, iy=0
            case 1: ix = 0;     iy = k;     break;  // {0,2}: iy varies, ix=0
            case 2: ix = p;     iy = k;     break;  // {1,3}: iy varies, ix=p
            case 3: ix = k;     iy = p;     break;  // {2,3}: ix varies, iy=p
            }
            edge_coeffs[static_cast<std::size_t>(k)] =
                parent_coeffs[static_cast<std::size_t>(tp_index_2d(ix, iy, p))];
        }
        break;
    }

    case cell::type::hexahedron:
    {
        // Basix hex vertex numbering:
        //   0=(0,0,0), 1=(1,0,0), 2=(0,1,0), 3=(1,1,0),
        //   4=(0,0,1), 5=(1,0,1), 6=(0,1,1), 7=(1,1,1)
        // For Basix ordering: v = bx + 2*by + 4*bz
        // cell_topology.h now returns Basix vertex IDs directly.
        for (int k = 0; k <= p; ++k)
        {
            const int x0 = (v0 & 1) * p, y0 = ((v0 >> 1) & 1) * p, z0 = ((v0 >> 2) & 1) * p;
            const int x1 = (v1 & 1) * p, y1 = ((v1 >> 1) & 1) * p, z1 = ((v1 >> 2) & 1) * p;

            const int ix = x0 + (x1 - x0) * k / p;
            const int iy = y0 + (y1 - y0) * k / p;
            const int iz = z0 + (z1 - z0) * k / p;

            edge_coeffs[static_cast<std::size_t>(k)] =
                parent_coeffs[static_cast<std::size_t>(tp_index_3d(ix, iy, iz, p))];
        }
        break;
    }

    default:
        throw std::runtime_error(
            "extract_parent_edge_bernstein: unsupported cell type");
    }
}

// =====================================================================
// restrict_edge_bernstein_exact
// =====================================================================

template <std::floating_point T>
void restrict_edge_bernstein_exact(cell::type parent_cell_type,
                                   int degree,
                                   std::span<const T> parent_coeffs,
                                   std::span<const T> xi_a,
                                   std::span<const T> xi_b,
                                   std::vector<T>& edge_coeffs)
{
    const int p = degree;
    const int tdim = cell::get_tdim(parent_cell_type);
    if (static_cast<int>(xi_a.size()) != tdim
        || static_cast<int>(xi_b.size()) != tdim)
    {
        throw std::runtime_error(
            "restrict_edge_bernstein_exact: endpoint dimension mismatch");
    }

    std::vector<T> power_coeffs(static_cast<std::size_t>(p + 1), T(0));

    if (bernstein::is_simplex(parent_cell_type))
    {
        std::vector<T> lambda_a(static_cast<std::size_t>(tdim + 1), T(0));
        std::vector<T> lambda_b(static_cast<std::size_t>(tdim + 1), T(0));

        lambda_a[0] = T(1);
        lambda_b[0] = T(1);
        for (int d = 0; d < tdim; ++d)
        {
            lambda_a[0] -= xi_a[static_cast<std::size_t>(d)];
            lambda_b[0] -= xi_b[static_cast<std::size_t>(d)];
            lambda_a[static_cast<std::size_t>(d + 1)] = xi_a[static_cast<std::size_t>(d)];
            lambda_b[static_cast<std::size_t>(d + 1)] = xi_b[static_cast<std::size_t>(d)];
        }

        int coeff_index = 0;
        if (tdim == 1)
        {
            for (int i = 0; i <= p; ++i, ++coeff_index)
            {
                const int alpha[2] = {p - i, i};
                std::vector<T> term = {T(1)};
                for (int c = 0; c < 2; ++c)
                {
                    const std::vector<T> factor = linear_power_poly(
                        lambda_a[static_cast<std::size_t>(c)],
                        lambda_b[static_cast<std::size_t>(c)] - lambda_a[static_cast<std::size_t>(c)],
                        alpha[c]);
                    term = multiply_poly<T>(term, factor);
                }
                add_scaled_poly(power_coeffs, std::span<const T>(term),
                                parent_coeffs[static_cast<std::size_t>(coeff_index)]
                                * multinomial<T>(p, alpha, 2));
            }
        }
        else if (tdim == 2)
        {
            for (int j = 0; j <= p; ++j)
            {
                for (int i = 0; i <= p - j; ++i, ++coeff_index)
                {
                    const int alpha[3] = {p - i - j, i, j};
                    std::vector<T> term = {T(1)};
                    for (int c = 0; c < 3; ++c)
                    {
                        const std::vector<T> factor = linear_power_poly(
                            lambda_a[static_cast<std::size_t>(c)],
                            lambda_b[static_cast<std::size_t>(c)] - lambda_a[static_cast<std::size_t>(c)],
                            alpha[c]);
                        term = multiply_poly<T>(term, factor);
                    }
                    add_scaled_poly(power_coeffs, std::span<const T>(term),
                                    parent_coeffs[static_cast<std::size_t>(coeff_index)]
                                    * multinomial<T>(p, alpha, 3));
                }
            }
        }
        else if (tdim == 3)
        {
            for (int k = 0; k <= p; ++k)
            {
                for (int j = 0; j <= p - k; ++j)
                {
                    for (int i = 0; i <= p - k - j; ++i, ++coeff_index)
                    {
                        const int alpha[4] = {p - i - j - k, i, j, k};
                        std::vector<T> term = {T(1)};
                        for (int c = 0; c < 4; ++c)
                        {
                            const std::vector<T> factor = linear_power_poly(
                                lambda_a[static_cast<std::size_t>(c)],
                                lambda_b[static_cast<std::size_t>(c)] - lambda_a[static_cast<std::size_t>(c)],
                                alpha[c]);
                            term = multiply_poly<T>(term, factor);
                        }
                        add_scaled_poly(power_coeffs, std::span<const T>(term),
                                        parent_coeffs[static_cast<std::size_t>(coeff_index)]
                                        * multinomial<T>(p, alpha, 4));
                    }
                }
            }
        }
        else
        {
            throw std::runtime_error(
                "restrict_edge_bernstein_exact: unsupported simplex tdim");
        }
    }
    else if (bernstein::is_tensor_product(parent_cell_type))
    {
        if (tdim == 2)
        {
            int coeff_index = 0;
            for (int iy = 0; iy <= p; ++iy)
            {
                for (int ix = 0; ix <= p; ++ix, ++coeff_index)
                {
                    std::vector<T> term = {T(1)};
                    const int ids[2] = {ix, iy};
                    for (int d = 0; d < 2; ++d)
                    {
                        const T a = xi_a[static_cast<std::size_t>(d)];
                        const T b = xi_b[static_cast<std::size_t>(d)] - a;
                        std::vector<T> dim_poly = linear_power_poly(a, b, ids[d]);
                        const std::vector<T> one_minus_poly =
                            linear_power_poly(T(1) - a, -b, p - ids[d]);
                        dim_poly = multiply_poly<T>(dim_poly, one_minus_poly);
                        const T scale = binomial<T>(p, ids[d]);
                        for (T& x : dim_poly)
                            x *= scale;
                        term = multiply_poly<T>(term, dim_poly);
                    }
                    add_scaled_poly(power_coeffs, std::span<const T>(term),
                                    parent_coeffs[static_cast<std::size_t>(coeff_index)]);
                }
            }
        }
        else if (tdim == 3)
        {
            int coeff_index = 0;
            for (int iz = 0; iz <= p; ++iz)
            {
                for (int iy = 0; iy <= p; ++iy)
                {
                    for (int ix = 0; ix <= p; ++ix, ++coeff_index)
                    {
                        std::vector<T> term = {T(1)};
                        const int ids[3] = {ix, iy, iz};
                        for (int d = 0; d < 3; ++d)
                        {
                            const T a = xi_a[static_cast<std::size_t>(d)];
                            const T b = xi_b[static_cast<std::size_t>(d)] - a;
                            std::vector<T> dim_poly = linear_power_poly(a, b, ids[d]);
                            const std::vector<T> one_minus_poly =
                                linear_power_poly(T(1) - a, -b, p - ids[d]);
                            dim_poly = multiply_poly<T>(dim_poly, one_minus_poly);
                            const T scale = binomial<T>(p, ids[d]);
                            for (T& x : dim_poly)
                                x *= scale;
                            term = multiply_poly<T>(term, dim_poly);
                        }
                        add_scaled_poly(power_coeffs, std::span<const T>(term),
                                        parent_coeffs[static_cast<std::size_t>(coeff_index)]);
                    }
                }
            }
        }
        else
        {
            throw std::runtime_error(
                "restrict_edge_bernstein_exact: unsupported tensor-product tdim");
        }
    }
    else
    {
        throw std::runtime_error(
            "restrict_edge_bernstein_exact: unsupported parent cell type");
    }

    power_to_bernstein_1d(std::span<const T>(power_coeffs), edge_coeffs);
}

// =====================================================================
// edge_is_on_single_parent_edge
// =====================================================================

template <std::floating_point T>
bool edge_is_on_single_parent_edge(const AdaptCell<T>& adapt_cell,
                                   int edge_id,
                                   int& parent_edge_id)
{
    // Get the two vertex indices of this edge.
    auto verts = adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(edge_id)];
    assert(verts.size() == 2);

    const int v0 = verts[0];
    const int v1 = verts[1];

    // Both must be on parent dimension 0 or 1 (vertex or edge).
    const int dim0 = adapt_cell.vertex_parent_dim[static_cast<std::size_t>(v0)];
    const int dim1 = adapt_cell.vertex_parent_dim[static_cast<std::size_t>(v1)];

    // Collect the set of parent edges each vertex lies on.
    // A vertex with parent_dim == 0 (parent vertex) lies on all parent edges
    // incident to that vertex. A vertex with parent_dim == 1 lies on exactly
    // one parent edge (parent_id).
    // For simplicity, and because the initial AdaptCell has only parent-vertex
    // provenance, we handle the case where both vertices are parent vertices:
    // check if they form a parent edge.

    if (dim0 == 1 && dim1 == 1)
    {
        // Both on parent edges; must be the same one.
        int pe0 = adapt_cell.vertex_parent_id[static_cast<std::size_t>(v0)];
        int pe1 = adapt_cell.vertex_parent_id[static_cast<std::size_t>(v1)];
        if (pe0 == pe1)
        {
            parent_edge_id = pe0;
            return true;
        }
        return false;
    }

    if (dim0 == 0 && dim1 == 0)
    {
        // Both are initial parent vertices.  In make_adapt_cell the vertices
        // are stored in Basix order, so AdaptCell vertex index v == Basix
        // local vertex index v.  Compare local indices directly against the
        // parent edge table (which also uses Basix vertex indices).
        auto parent_edges = cell::edges(adapt_cell.parent_cell_type);
        for (int e = 0; e < static_cast<int>(parent_edges.size()); ++e)
        {
            const auto& pe = parent_edges[static_cast<std::size_t>(e)];
            if ((pe[0] == v0 && pe[1] == v1)
                || (pe[0] == v1 && pe[1] == v0))
            {
                parent_edge_id = e;
                return true;
            }
        }
        return false;
    }

    if ((dim0 == 0 && dim1 == 1) || (dim0 == 1 && dim1 == 0))
    {
        // One is a parent vertex, the other is on a parent edge.
        // Check that the parent vertex is an endpoint of that parent edge.
        int vert_v   = (dim0 == 0) ? v0 : v1;
        int edge_v   = (dim0 == 0) ? v1 : v0;
        int gv = adapt_cell.vertex_parent_id[static_cast<std::size_t>(vert_v)];
        int pe = adapt_cell.vertex_parent_id[static_cast<std::size_t>(edge_v)];

        auto parent_edges = cell::edges(adapt_cell.parent_cell_type);
        const auto& edge_pair = parent_edges[static_cast<std::size_t>(pe)];
        if (edge_pair[0] == gv || edge_pair[1] == gv)
        {
            parent_edge_id = pe;
            return true;
        }
        return false;
    }

    // Other combinations (face, interior) → not on a single parent edge.
    return false;
}

// =====================================================================
// classify_edge_roots
// =====================================================================

template <std::floating_point T>
EdgeRootTag classify_edge_roots(std::span<const T> edge_coeffs,
                                T zero_tol, T sign_tol,
                                int max_depth,
                                T& green_split_t,
                                bool& has_green_split_t)
{
    has_green_split_t = false;
    green_split_t = T(0);

    if (edge_coeffs.empty())
        return EdgeRootTag::not_classified;

    // 1. All zero?
    if (bernstein_all_zero(edge_coeffs, zero_tol))
        return EdgeRootTag::zero;

    // 2. All positive or all negative?
    if (bernstein_all_positive(edge_coeffs, sign_tol))
        return EdgeRootTag::no_root;
    if (bernstein_all_negative(edge_coeffs, sign_tol))
        return EdgeRootTag::no_root;

    const bool left_zero = std::fabs(edge_coeffs.front()) <= zero_tol;
    const bool right_zero = std::fabs(edge_coeffs.back()) <= zero_tol;

    // 3. Find root intervals.
    std::vector<EdgeRootInterval<T>> intervals;
    bool has_zero_segment = false;
    find_root_intervals_1d<T>(edge_coeffs, T(0), T(1),
                              zero_tol, sign_tol, 0, max_depth,
                              intervals, has_zero_segment);

    if (has_zero_segment)
        return EdgeRootTag::zero;

    // ---------------------------------------------------------------
    // Check: "one vertex zero AND an interior root" must be multiple_roots.
    //
    // If one endpoint has phi≈0 (left_zero / right_zero) AND the root
    // finder also found a strictly interior interval (t0>0 or t1<1),
    // merging them into a single one_root would place the cut at the
    // existing vertex (wrong position) and produce degenerate geometry.
    // Classify as multiple_roots here, *before* any merge step, so
    // green refinement splits the edge between the vertex-zero event
    // and the interior sign change.
    // ---------------------------------------------------------------
    if (left_zero && !intervals.empty())
    {
        T ep_zero_end     = T(0);     // rightmost bound of zero-endpoint intervals
        T first_interior  = T(1);     // leftmost start of strictly interior intervals
        bool found_ep     = false;
        bool found_int    = false;
        for (const auto& iv : intervals)
        {
            if (iv.t0 == T(0))
            {
                ep_zero_end = std::max(ep_zero_end, iv.t1);
                found_ep = true;
            }
            else // iv.t0 > 0 → strictly interior
            {
                first_interior = std::min(first_interior, iv.t0);
                found_int = true;
            }
        }
        if (found_int)
        {
            T gap_lo = found_ep ? ep_zero_end : T(0);
            green_split_t = (gap_lo + first_interior) * T(0.5);
            green_split_t = std::clamp(green_split_t, T(0), T(1));
            has_green_split_t = true;
            return EdgeRootTag::multiple_roots;
        }
    }
    if (right_zero && !intervals.empty())
    {
        T ep_zero_start   = T(1);     // leftmost bound of zero-endpoint intervals
        T last_interior   = T(0);     // rightmost end of strictly interior intervals
        bool found_ep     = false;
        bool found_int    = false;
        for (const auto& iv : intervals)
        {
            if (iv.t1 == T(1))
            {
                ep_zero_start = std::min(ep_zero_start, iv.t0);
                found_ep = true;
            }
            else // iv.t1 < 1 → strictly interior
            {
                last_interior = std::max(last_interior, iv.t1);
                found_int = true;
            }
        }
        if (found_int)
        {
            T gap_hi = found_ep ? ep_zero_start : T(1);
            green_split_t = (last_interior + gap_hi) * T(0.5);
            green_split_t = std::clamp(green_split_t, T(0), T(1));
            has_green_split_t = true;
            return EdgeRootTag::multiple_roots;
        }
    }

    // Merge overlapping intervals.
    T merge_tol = T(4) * std::numeric_limits<T>::epsilon();
    if (max_depth > 0)
    {
        T width = T(1);
        for (int i = 0; i < max_depth; ++i)
            width *= T(0.5);
        merge_tol = std::max(merge_tol, T(2) * width);
    }
    merge_root_intervals(intervals, merge_tol);

    std::vector<EdgeRootInterval<T>> clusters = intervals;

    bool left_covered = false;
    bool right_covered = false;
    for (const auto& interval : clusters)
    {
        if (interval.t0 <= merge_tol)
            left_covered = true;
        if (interval.t1 >= T(1) - merge_tol)
            right_covered = true;
    }

    if (left_zero && !left_covered)
        clusters.push_back({T(0), T(0)});
    if (right_zero && !right_covered)
        clusters.push_back({T(1), T(1)});

    if (clusters.empty())
        return EdgeRootTag::no_root;

    std::sort(clusters.begin(), clusters.end(),
              [](const EdgeRootInterval<T>& a, const EdgeRootInterval<T>& b)
              { return a.t0 < b.t0; });
    merge_root_intervals(clusters, merge_tol);

    if (clusters.size() == 1)
        return EdgeRootTag::one_root;

    green_split_t = (clusters[0].t1 + clusters[1].t0) * T(0.5);
    green_split_t = std::clamp(green_split_t, T(0), T(1));
    has_green_split_t = true;
    return EdgeRootTag::multiple_roots;
}

// =====================================================================
// classify_new_edges
// =====================================================================

template <std::floating_point T, std::integral I>
void classify_new_edges(AdaptCell<T>& adapt_cell,
                        const LevelSetCell<T, I>& ls_cell,
                        int level_set_id,
                        T zero_tol, T sign_tol,
                        int max_depth)
{
    const int n_edges = adapt_cell.n_entities(1);

    // Ensure tag storage is large enough.
    if (adapt_cell.edge_root_tag_num_level_sets <= level_set_id)
        adapt_cell.resize_edge_root_tags(level_set_id + 1);
    if (adapt_cell.edge_green_split_has_value.size()
        < static_cast<std::size_t>((level_set_id + 1) * n_edges))
        adapt_cell.resize_green_split_data(level_set_id + 1);
    if (adapt_cell.edge_one_root_has_value.size()
        < static_cast<std::size_t>((level_set_id + 1) * n_edges))
        adapt_cell.resize_one_root_data(level_set_id + 1);

    std::vector<T> edge_coeffs;

    for (int e = 0; e < n_edges; ++e)
    {
        // Skip already classified edges.
        if (adapt_cell.get_edge_root_tag(level_set_id, e)
            != EdgeRootTag::not_classified)
            continue;

        // Extract edge Bernstein coefficients on the actual adaptive edge segment.
        //
        // Even if an adaptive edge lies on a parent edge, using the parent-edge
        // fast-path can be incorrect for subsegments. Restrict by endpoints so
        // the endpoint coefficients and green split parameter use this edge's
        // local parameterization.
        auto verts = adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(e)];
        const int tdim = adapt_cell.tdim;
        std::span<const T> xi_a(
            adapt_cell.vertex_coords.data()
                + static_cast<std::size_t>(verts[0]) * static_cast<std::size_t>(tdim),
            static_cast<std::size_t>(tdim));
        std::span<const T> xi_b(
            adapt_cell.vertex_coords.data()
                + static_cast<std::size_t>(verts[1]) * static_cast<std::size_t>(tdim),
            static_cast<std::size_t>(tdim));

        restrict_edge_bernstein_exact(
            ls_cell.cell_type, ls_cell.bernstein_order,
            std::span<const T>(ls_cell.bernstein_coeffs),
            xi_a, xi_b, edge_coeffs);

        // Classify.
        T green_split_t = T(0);
        bool has_green = false;
        EdgeRootTag tag = classify_edge_roots<T>(
            std::span<const T>(edge_coeffs),
            zero_tol, sign_tol, max_depth,
            green_split_t, has_green);

        T linear_root_t = T(0);
        const bool has_linear_root =
            tag == EdgeRootTag::one_root
            && edge_coeffs.size() >= 2
            && linear_one_root_parameter_from_endpoint_values<T>(
                edge_coeffs.front(), edge_coeffs.back(), zero_tol, linear_root_t);
        if (tag == EdgeRootTag::one_root && !has_linear_root)
            tag = EdgeRootTag::no_root;

        adapt_cell.set_edge_root_tag(level_set_id, e, tag);

        const auto idx = static_cast<std::size_t>(level_set_id * n_edges + e);
        adapt_cell.edge_green_split_param[idx] = green_split_t;
        adapt_cell.edge_green_split_has_value[idx] = has_green ? 1 : 0;
        if (tag == EdgeRootTag::one_root)
        {
            adapt_cell.edge_one_root_param[idx] = linear_root_t;
            adapt_cell.edge_one_root_vertex_id[idx] = -1;
            adapt_cell.edge_one_root_has_value[idx] = 1;
        }
        else
        {
            adapt_cell.edge_one_root_param[idx] = T(0);
            adapt_cell.edge_one_root_vertex_id[idx] = -1;
            adapt_cell.edge_one_root_has_value[idx] = 0;
        }
    }
}

// =====================================================================
// Explicit template instantiations
// =====================================================================

template void subdivide_bernstein_1d(std::span<const double>, double,
                                     std::vector<double>&, std::vector<double>&);
template void subdivide_bernstein_1d(std::span<const float>, float,
                                     std::vector<float>&, std::vector<float>&);

template bool bernstein_all_positive(std::span<const double>, double);
template bool bernstein_all_positive(std::span<const float>, float);
template bool bernstein_all_negative(std::span<const double>, double);
template bool bernstein_all_negative(std::span<const float>, float);
template bool bernstein_all_zero(std::span<const double>, double);
template bool bernstein_all_zero(std::span<const float>, float);

template void find_root_intervals_1d(std::span<const double>, double, double,
                                     double, double, int, int,
                                     std::vector<EdgeRootInterval<double>>&, bool&);
template void find_root_intervals_1d(std::span<const float>, float, float,
                                     float, float, int, int,
                                     std::vector<EdgeRootInterval<float>>&, bool&);

template void extract_parent_edge_bernstein(cell::type, int,
                                            std::span<const double>, int,
                                            std::vector<double>&);
template void extract_parent_edge_bernstein(cell::type, int,
                                            std::span<const float>, int,
                                            std::vector<float>&);

template void restrict_edge_bernstein_exact(cell::type, int,
                                            std::span<const double>,
                                            std::span<const double>,
                                            std::span<const double>,
                                            std::vector<double>&);
template void restrict_edge_bernstein_exact(cell::type, int,
                                            std::span<const float>,
                                            std::span<const float>,
                                            std::span<const float>,
                                            std::vector<float>&);

template bool edge_is_on_single_parent_edge(const AdaptCell<double>&, int, int&);
template bool edge_is_on_single_parent_edge(const AdaptCell<float>&, int, int&);

template EdgeRootTag classify_edge_roots(std::span<const double>, double, double,
                                         int, double&, bool&);
template EdgeRootTag classify_edge_roots(std::span<const float>, float, float,
                                         int, float&, bool&);

template void classify_new_edges(AdaptCell<double>&, const LevelSetCell<double, int>&,
                                 int, double, double, int);
template void classify_new_edges(AdaptCell<float>&, const LevelSetCell<float, int>&,
                                 int, float, float, int);
template void classify_new_edges(AdaptCell<double>&, const LevelSetCell<double, long>&,
                                 int, double, double, int);
template void classify_new_edges(AdaptCell<float>&, const LevelSetCell<float, long>&,
                                 int, float, float, int);

} // namespace cutcells
