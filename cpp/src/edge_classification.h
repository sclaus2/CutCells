// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include "local_level_set.h"
#include "local_mesh.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace cutcells
{

template <std::floating_point T, std::integral I>
bool classify_edges_and_mark_refine(
    LocalMesh<T>&                  mesh,
    const LevelSetFunction<T, I>&  level_set,
    int                            level_set_id = 0,
    T                              tol = static_cast<T>(1e-14),
    bool                           refine_on_uncertain = false);

template <std::floating_point T>
inline int sign_with_tol(T v, T tol)
{
    if (std::abs(v) <= tol)
        return 0;
    return (v < T(0)) ? -1 : 1;
}

template <std::floating_point T>
inline EdgeState classify_edge_from_nodal_samples(
    T                   v0,
    T                   v1,
    std::span<const T>  interior_values,
    T                   tol)
{
    const int s0 = sign_with_tol(v0, tol);
    const int s1 = sign_with_tol(v1, tol);

    if (s0 == 0 && s1 == 0)
        return EdgeState::zero_edge;

    for (const T vi : interior_values)
    {
        if (sign_with_tol(vi, tol) == 0)
            return EdgeState::uncertain;
    }

    if (s0 != s1)
        return EdgeState::single_cross;

    for (const T vi : interior_values)
    {
        const int si = sign_with_tol(vi, tol);
        if (si != 0 && si != s0)
            return EdgeState::multi_cross;
    }

    return EdgeState::uncut;
}

template <std::floating_point T>
struct EdgeZeroTopology
{
    uint8_t endpoint_zero_count = 0;
    uint8_t interior_root_count = 0;
    EdgeState state = EdgeState::uncertain;
};

template <std::floating_point T>
inline uint8_t count_distinct_interior_roots_from_samples(
    std::span<const std::pair<T, T>> samples,
    T                                tol,
    T                                root_endpoint_tol)
{
    if (samples.empty())
        return 0;

    uint8_t count = 0;
    int prev_sign = 0;
    bool prev_kept = false;

    for (const auto& [t, value] : samples)
    {
        const int sv = sign_with_tol(value, tol);
        if (sv == 0)
        {
            if (t > root_endpoint_tol && t < T(1) - root_endpoint_tol)
                ++count;
            prev_sign = 0;
            prev_kept = true;
            continue;
        }

        if (prev_kept && prev_sign != 0 && sv != prev_sign)
            ++count;

        prev_sign = sv;
        prev_kept = true;
    }

    return count;
}

template <std::floating_point T>
inline EdgeZeroTopology<T> classify_edge_zero_topology(
    T                               v0,
    T                               v1,
    std::span<const std::pair<T, T>> interior_samples,
    EdgeState                       raw_state,
    T                               tol,
    T                               root_endpoint_tol)
{
    EdgeZeroTopology<T> out;
    out.endpoint_zero_count = static_cast<uint8_t>(
        (sign_with_tol(v0, tol) == 0 ? 1 : 0)
      + (sign_with_tol(v1, tol) == 0 ? 1 : 0));

    out.interior_root_count = count_distinct_interior_roots_from_samples(
        interior_samples, tol, root_endpoint_tol);

    if (out.endpoint_zero_count == 2)
    {
        out.state = (out.interior_root_count == 0)
                  ? EdgeState::zero_edge
                  : EdgeState::multi_cross;
        return out;
    }

    if (out.endpoint_zero_count == 1)
    {
        if (out.interior_root_count == 0)
            out.state = EdgeState::uncut;
        else if (out.interior_root_count == 1)
            out.state = EdgeState::single_cross;
        else
            out.state = EdgeState::multi_cross;
        return out;
    }

    if (out.interior_root_count > 1)
    {
        out.state = EdgeState::multi_cross;
        return out;
    }

    out.state = raw_state;
    return out;
}

template <std::floating_point T>
inline std::vector<T> bernstein_coefficients_from_equispaced_nodal(
    std::span<const T> nodal_values)
{
    const int p = static_cast<int>(nodal_values.size()) - 1;
    if (p < 0)
        return {};

    // Build interpolation matrix A_{j,i} = B_i^p(j/p), solve A b = f.
    std::vector<T> A(static_cast<std::size_t>((p + 1) * (p + 1)), T(0));
    auto a = [&](int r, int c) -> T& { return A[static_cast<std::size_t>(r * (p + 1) + c)]; };

    for (int j = 0; j <= p; ++j)
    {
        const T x = (p == 0) ? T(0) : static_cast<T>(j) / static_cast<T>(p);
        for (int i = 0; i <= p; ++i)
        {
            const T bi = binomial_coeff<T>(p, i)
                         * std::pow(x, static_cast<T>(i))
                         * std::pow(T(1) - x, static_cast<T>(p - i));
            a(j, i) = bi;
        }
    }

    std::vector<T> b(nodal_values.begin(), nodal_values.end());

    // Gaussian elimination with partial pivoting.
    for (int k = 0; k <= p; ++k)
    {
        int piv = k;
        T maxv = std::abs(a(k, k));
        for (int r = k + 1; r <= p; ++r)
        {
            const T v = std::abs(a(r, k));
            if (v > maxv)
            {
                maxv = v;
                piv = r;
            }
        }
        if (maxv <= std::numeric_limits<T>::epsilon())
            return {};

        if (piv != k)
        {
            for (int c = k; c <= p; ++c)
                std::swap(a(k, c), a(piv, c));
            std::swap(b[static_cast<std::size_t>(k)], b[static_cast<std::size_t>(piv)]);
        }

        const T diag = a(k, k);
        for (int c = k; c <= p; ++c)
            a(k, c) /= diag;
        b[static_cast<std::size_t>(k)] /= diag;

        for (int r = k + 1; r <= p; ++r)
        {
            const T f = a(r, k);
            if (std::abs(f) <= std::numeric_limits<T>::epsilon())
                continue;
            for (int c = k; c <= p; ++c)
                a(r, c) -= f * a(k, c);
            b[static_cast<std::size_t>(r)] -= f * b[static_cast<std::size_t>(k)];
        }
    }

    std::vector<T> x(static_cast<std::size_t>(p + 1), T(0));
    for (int r = p; r >= 0; --r)
    {
        T s = b[static_cast<std::size_t>(r)];
        for (int c = r + 1; c <= p; ++c)
            s -= a(r, c) * x[static_cast<std::size_t>(c)];
        x[static_cast<std::size_t>(r)] = s;
    }
    return x;
}

template <std::floating_point T>
inline T evaluate_bernstein_poly_1d(std::span<const T> coeffs, T t)
{
    if (coeffs.empty())
        return T(0);

    std::vector<T> work(coeffs.begin(), coeffs.end());
    const int p = static_cast<int>(work.size()) - 1;
    for (int r = 1; r <= p; ++r)
    {
        for (int i = 0; i <= p - r; ++i)
            work[static_cast<std::size_t>(i)]
                = (T(1) - t) * work[static_cast<std::size_t>(i)]
                + t * work[static_cast<std::size_t>(i + 1)];
    }
    return work[0];
}

template <std::floating_point T>
inline EdgeState classify_edge_from_bernstein_coeffs(
    std::span<const T> coeffs,
    bool endpoint_sign_change,
    T tol)
{
    if (coeffs.empty())
        return EdgeState::uncertain;

    if (bernstein_all_strict_one_sign(coeffs, tol))
        return EdgeState::uncut;
    bool all_near_zero = true;
    for (const T c : coeffs)
        all_near_zero = all_near_zero && (std::abs(c) <= tol);
    if (all_near_zero)
        return EdgeState::zero_edge;

    if (bernstein_derivative_strict_one_sign(coeffs, tol))
        return endpoint_sign_change ? EdgeState::single_cross : EdgeState::uncut;

    int sign_changes = 0;
    int prev_sign = 0;
    for (const T c : coeffs)
    {
        const int sc = sign_with_tol(c, tol);
        if (sc == 0)
        {
            continue;
        }
        if (prev_sign != 0 && sc != prev_sign)
            ++sign_changes;
        prev_sign = sc;
    }

    if (sign_changes == 0)
        return EdgeState::uncut;
    if (sign_changes == 1 && endpoint_sign_change)
        return EdgeState::single_cross;
    if (sign_changes >= 2)
        return EdgeState::multi_cross;
    return EdgeState::uncertain;
}

template <std::floating_point T>
inline bool bernstein_all_near_zero(std::span<const T> coeffs, T tol)
{
    for (const T c : coeffs)
    {
        if (std::abs(c) > tol)
            return false;
    }
    return true;
}

template <std::floating_point T>
struct EdgeBernsteinResolution
{
    EdgeState state = EdgeState::uncertain;
    EdgeCertification cert = EdgeCertification::unresolved;
};

template <std::floating_point T>
struct BernsteinSubdivisionStats
{
    int n_cross = 0;
    int n_zero_intervals = 0;
    bool unresolved = false;
};

template <std::floating_point T>
inline void analyze_bernstein_edge_subdivision(
    std::span<const T>          coeffs,
    int                         depth,
    int                         max_depth,
    T                           tol,
    BernsteinSubdivisionStats<T>& stats)
{
    if (coeffs.empty())
    {
        stats.unresolved = true;
        return;
    }

    if (bernstein_all_strict_one_sign(coeffs, tol))
        return;

    if (bernstein_all_near_zero(coeffs, tol))
    {
        ++stats.n_zero_intervals;
        return;
    }

    const bool monotone = bernstein_derivative_strict_one_sign(coeffs, tol);
    const int s0 = sign_with_tol(coeffs.front(), tol);
    const int s1 = sign_with_tol(coeffs.back(), tol);
    if (monotone)
    {
        if (s0 != 0 && s1 != 0 && s0 != s1)
            ++stats.n_cross;
        return;
    }

    if (depth >= max_depth)
    {
        stats.unresolved = true;
        return;
    }

    std::vector<T> left;
    std::vector<T> right;
    de_casteljau_split_1d<T>(coeffs, left, right, T(0.5));
    analyze_bernstein_edge_subdivision(
        std::span<const T>(left.data(), left.size()), depth + 1, max_depth, tol, stats);
    if (stats.unresolved)
        return;
    analyze_bernstein_edge_subdivision(
        std::span<const T>(right.data(), right.size()), depth + 1, max_depth, tol, stats);
}

template <std::floating_point T>
inline EdgeBernsteinResolution<T> classify_edge_from_bernstein_subdivision(
    std::span<const T> coeffs,
    T                  tol,
    int                max_depth = 10)
{
    BernsteinSubdivisionStats<T> stats;
    analyze_bernstein_edge_subdivision(coeffs, 0, max_depth, tol, stats);

    EdgeBernsteinResolution<T> out;
    out.cert = stats.unresolved
                   ? EdgeCertification::unresolved
                   : EdgeCertification::resolved_by_subdivision;
    if (stats.unresolved)
    {
        out.state = EdgeState::uncertain;
        return out;
    }

    if (stats.n_zero_intervals > 0 && stats.n_cross == 0)
    {
        out.state = (stats.n_zero_intervals == 1)
                        ? EdgeState::zero_edge
                        : EdgeState::multi_cross;
        return out;
    }

    const int complexity = stats.n_cross + stats.n_zero_intervals;
    if (complexity == 0)
    {
        out.state = EdgeState::uncut;
        return out;
    }
    if (complexity == 1)
    {
        out.state = (stats.n_cross == 1)
                        ? EdgeState::single_cross
                        : EdgeState::multi_cross;
        return out;
    }
    out.state = EdgeState::multi_cross;
    return out;
}

template <std::floating_point T>
inline EdgeBernsteinResolution<T> classify_edge_from_full_cell_bernstein(
    std::span<const T> edge_coeffs,
    T                  tol)
{
    EdgeBernsteinResolution<T> out;

    if (edge_coeffs.empty())
        return out;

    if (bernstein_all_strict_one_sign(edge_coeffs, tol))
    {
        out.state = EdgeState::uncut;
        out.cert = EdgeCertification::certified;
        return out;
    }

    if (bernstein_all_near_zero(edge_coeffs, tol))
    {
        out.state = EdgeState::zero_edge;
        out.cert = EdgeCertification::certified;
        return out;
    }

    const int sign_changes = bernstein_sign_variation_count(edge_coeffs, tol);
    const bool monotone = bernstein_derivative_strict_one_sign(edge_coeffs, tol);
    const int s0 = sign_with_tol(edge_coeffs.front(), tol);
    const int s1 = sign_with_tol(edge_coeffs.back(), tol);
    const bool endpoint_sign_change = (s0 != 0 && s1 != 0 && s0 != s1);

    if (monotone && endpoint_sign_change)
    {
        out.state = EdgeState::single_cross;
        out.cert = EdgeCertification::certified;
        return out;
    }
    if (monotone && s0 != 0 && s1 != 0 && s0 == s1)
    {
        out.state = EdgeState::uncut;
        out.cert = EdgeCertification::certified;
        return out;
    }
    if (sign_changes >= 2)
    {
        out.state = EdgeState::multi_cross;
        out.cert = EdgeCertification::resolved_by_sampling;
        return out;
    }

    return classify_edge_from_bernstein_subdivision<T>(edge_coeffs, tol);
}

inline EdgeState combine_nodal_and_bernstein(EdgeState nodal, EdgeState bernstein)
{
    if (nodal == EdgeState::multi_cross || bernstein == EdgeState::multi_cross)
        return EdgeState::multi_cross;
    if (nodal == EdgeState::uncertain)
        return bernstein;
    if (bernstein == EdgeState::uncertain)
        return nodal;

    if (nodal == EdgeState::zero_edge || bernstein == EdgeState::zero_edge)
    {
        if (nodal == EdgeState::single_cross || bernstein == EdgeState::single_cross)
            return EdgeState::multi_cross;
        return EdgeState::zero_edge;
    }

    if (nodal == bernstein)
        return nodal;
    return EdgeState::uncertain;
}

template <std::floating_point T>
inline std::vector<std::pair<T, T>> sample_interior_bernstein_edge(
    std::span<const T> edge_coeffs)
{
    std::vector<std::pair<T, T>> samples;
    const int p = static_cast<int>(edge_coeffs.size()) - 1;
    if (p <= 1)
        return samples;

    samples.reserve(static_cast<std::size_t>(p - 1));
    for (int i = 1; i < p; ++i)
    {
        const T t = static_cast<T>(i) / static_cast<T>(p);
        samples.push_back({t, evaluate_bernstein_poly_1d(edge_coeffs, t)});
    }
    return samples;
}

template <std::floating_point T, std::integral I = int>
inline void populate_local_level_set_values(
    LocalMesh<T>&                 mesh,
    const LevelSetFunction<T, I>& level_set,
    int                           level_set_id,
    T                             tol,
    std::vector<uint8_t>&         value_available,
    int                           nodal_prefix_count = -1)
{
    const int nv = mesh.n_vertices();
    const bool use_nodal = level_set.has_nodal_values()
                           && nodal_prefix_count > 0
                           && static_cast<int>(level_set.nodal_values.size()) >= nodal_prefix_count;
    const bool use_value = level_set.has_value();

    value_available.assign(static_cast<std::size_t>(nv), 0);

    if (mesh.vertex_phi.size() != static_cast<std::size_t>(nv * mesh.n_level_sets))
        mesh.vertex_phi.assign(static_cast<std::size_t>(nv * mesh.n_level_sets), T(0));
    if (mesh.vertex_zero_mask.size() != static_cast<std::size_t>(nv))
        mesh.vertex_zero_mask.assign(static_cast<std::size_t>(nv), 0);
    if (mesh.vertex_inside_mask.size() != static_cast<std::size_t>(nv))
        mesh.vertex_inside_mask.assign(static_cast<std::size_t>(nv), 0);

    for (int i = 0; i < nv; ++i)
    {
        T v = T(0);
        bool ok = false;
        if (use_nodal && i < nodal_prefix_count)
        {
            v = level_set.nodal_values[static_cast<std::size_t>(i)];
            ok = true;
        }
        else if (use_value)
        {
            const T* x = mesh.vertex_x.data() + static_cast<std::size_t>(i * mesh.gdim);
            v = level_set.value(x, static_cast<I>(mesh.parent_cell_id));
            ok = true;
        }

        if (!ok)
            continue;

        value_available[static_cast<std::size_t>(i)] = 1;
        mesh.vertex_phi[static_cast<std::size_t>(i * mesh.n_level_sets + level_set_id)] = v;

        const uint64_t mask = (uint64_t(1) << level_set_id);
        mesh.vertex_zero_mask[static_cast<std::size_t>(i)] &= ~mask;
        mesh.vertex_inside_mask[static_cast<std::size_t>(i)] &= ~mask;
        const int s = sign_with_tol(v, tol);
        if (s == 0)
            mesh.vertex_zero_mask[static_cast<std::size_t>(i)] |= mask;
        else if (s < 0)
            mesh.vertex_inside_mask[static_cast<std::size_t>(i)] |= mask;
    }
}

template <std::floating_point T, std::integral I = int>
inline void populate_local_level_set_values(
    LocalMesh<T>&                      mesh,
    const LocalLevelSetFunction<T, I>& level_set,
    int                                level_set_id,
    T                                  tol,
    std::vector<uint8_t>&              value_available)
{
    const int nv = mesh.n_vertices();
    value_available.assign(static_cast<std::size_t>(nv), 0);

    if (mesh.vertex_phi.size() != static_cast<std::size_t>(nv * mesh.n_level_sets))
        mesh.vertex_phi.assign(static_cast<std::size_t>(nv * mesh.n_level_sets), T(0));
    if (mesh.vertex_zero_mask.size() != static_cast<std::size_t>(nv))
        mesh.vertex_zero_mask.assign(static_cast<std::size_t>(nv), 0);
    if (mesh.vertex_inside_mask.size() != static_cast<std::size_t>(nv))
        mesh.vertex_inside_mask.assign(static_cast<std::size_t>(nv), 0);

    if (!level_set.has_value())
        return;

    const bool use_ref = !mesh.vertex_ref_x.empty()
                         && mesh.vertex_ref_x.size()
                                == static_cast<std::size_t>(mesh.n_vertices() * mesh.tdim);
    const int coord_dim = use_ref ? mesh.tdim : mesh.gdim;
    for (int i = 0; i < nv; ++i)
    {
        const T* x = use_ref
                         ? mesh.vertex_ref_x.data() + static_cast<std::size_t>(i * mesh.tdim)
                         : mesh.vertex_x.data() + static_cast<std::size_t>(i * mesh.gdim);
        const T v = level_set.value(x);
        value_available[static_cast<std::size_t>(i)] = 1;
        mesh.vertex_phi[static_cast<std::size_t>(i * mesh.n_level_sets + level_set_id)] = v;

        const uint64_t mask = (uint64_t(1) << level_set_id);
        mesh.vertex_zero_mask[static_cast<std::size_t>(i)] &= ~mask;
        mesh.vertex_inside_mask[static_cast<std::size_t>(i)] &= ~mask;
        const int s = sign_with_tol(v, tol);
        if (s == 0)
            mesh.vertex_zero_mask[static_cast<std::size_t>(i)] |= mask;
        else if (s < 0)
            mesh.vertex_inside_mask[static_cast<std::size_t>(i)] |= mask;
        (void)coord_dim;
    }
}

template <std::floating_point T, std::integral I = int>
inline int bernstein_level_set_degree(
    const LocalMesh<T>&                 mesh,
    const LevelSetFunction<T, I>&       level_set)
{
    if (level_set.degree >= 1)
        return level_set.degree;
    return infer_lagrange_order_from_num_nodes(
        mesh.parent_cell_type, static_cast<int>(level_set.nodal_values.size()));
}

template <std::floating_point T>
inline bool interior_node_refinement_test(
    const LocalMesh<T>&           mesh,
    const std::vector<uint8_t>&   value_available,
    int                           level_set_id,
    T                             tol)
{
    const int nv = mesh.n_vertices();
    for (int c = 0; c < mesh.n_cells(); ++c)
    {
        const int c0 = mesh.cell_offsets[static_cast<std::size_t>(c)];
        const int c1 = mesh.cell_offsets[static_cast<std::size_t>(c + 1)];

        int vertex_sign = 0;
        bool mixed = false;
        for (int j = c0; j < c1; ++j)
        {
            const int vi = mesh.cell_vertices[static_cast<std::size_t>(j)];
            if (!value_available[static_cast<std::size_t>(vi)])
                continue;
            const T vv = mesh.vertex_phi[static_cast<std::size_t>(vi * mesh.n_level_sets + level_set_id)];
            const int sv = sign_with_tol(vv, tol);
            if (sv == 0)
                continue;
            if (vertex_sign == 0)
                vertex_sign = sv;
            else if (vertex_sign != sv)
            {
                mixed = true;
                break;
            }
        }
        if (mixed || vertex_sign == 0)
            continue;

        for (int vi = 0; vi < nv; ++vi)
        {
            if (mesh.vertex_parent_dim[static_cast<std::size_t>(vi)] != mesh.tdim)
                continue;
            if (!value_available[static_cast<std::size_t>(vi)])
                continue;
            const T v = mesh.vertex_phi[static_cast<std::size_t>(vi * mesh.n_level_sets + level_set_id)];
            const int si = sign_with_tol(v, tol);
            if (si != 0 && si != vertex_sign)
                return true;
        }
    }

    return false;
}

template <std::floating_point T, std::integral I = int>
inline void populate_local_level_set_values_from_parent_polynomial(
    LocalMesh<T>&                       mesh,
    const ParentPolynomialContext<T>&   parent_poly,
    const LevelSetFunction<T, I>&       level_set,
    int                                 level_set_id,
    T                                   tol,
    std::vector<uint8_t>&               value_available)
{
    const int nv = mesh.n_vertices();
    const int nodal_prefix_count = lagrange_node_count(parent_poly.cell_type, parent_poly.degree);
    value_available.assign(static_cast<std::size_t>(nv), 0);

    if (mesh.vertex_phi.size() != static_cast<std::size_t>(nv * mesh.n_level_sets))
        mesh.vertex_phi.assign(static_cast<std::size_t>(nv * mesh.n_level_sets), T(0));
    if (mesh.vertex_zero_mask.size() != static_cast<std::size_t>(nv))
        mesh.vertex_zero_mask.assign(static_cast<std::size_t>(nv), 0);
    if (mesh.vertex_inside_mask.size() != static_cast<std::size_t>(nv))
        mesh.vertex_inside_mask.assign(static_cast<std::size_t>(nv), 0);

    for (int i = 0; i < nv; ++i)
    {
        T v = T(0);
        bool ok = false;
        if (i < nodal_prefix_count)
        {
            v = level_set.nodal_values[static_cast<std::size_t>(i)];
            ok = true;
        }
        else if (!mesh.vertex_ref_x.empty())
        {
            const std::span<const T> x_ref(
                mesh.vertex_ref_x.data() + static_cast<std::size_t>(i * mesh.tdim),
                static_cast<std::size_t>(mesh.tdim));
            v = evaluate_bernstein_cell<T>(parent_poly.bernstein, x_ref);
            ok = true;
        }
        else if (level_set.has_value())
        {
            const T* x = mesh.vertex_x.data() + static_cast<std::size_t>(i * mesh.gdim);
            v = level_set.value(x, static_cast<I>(mesh.parent_cell_id));
            ok = true;
        }

        if (!ok)
            continue;

        value_available[static_cast<std::size_t>(i)] = 1;
        mesh.vertex_phi[static_cast<std::size_t>(i * mesh.n_level_sets + level_set_id)] = v;

        const uint64_t mask = (uint64_t(1) << level_set_id);
        mesh.vertex_zero_mask[static_cast<std::size_t>(i)] &= ~mask;
        mesh.vertex_inside_mask[static_cast<std::size_t>(i)] &= ~mask;
        const int s = sign_with_tol(v, tol);
        if (s == 0)
            mesh.vertex_zero_mask[static_cast<std::size_t>(i)] |= mask;
        else if (s < 0)
            mesh.vertex_inside_mask[static_cast<std::size_t>(i)] |= mask;
    }
}

template <std::floating_point T>
inline bool point_in_local_cell_ref(
    const LocalMesh<T>& mesh,
    int                 cell_id,
    std::span<const T>  x_ref,
    T                   tol)
{
    const int c0 = mesh.cell_offsets[static_cast<std::size_t>(cell_id)];
    const auto ct = mesh.cell_types[static_cast<std::size_t>(cell_id)];
    const auto xv = [&](int local_v, int d) -> T
    {
        const int gv = mesh.cell_vertices[static_cast<std::size_t>(c0 + local_v)];
        return mesh.vertex_ref_x[static_cast<std::size_t>(gv * mesh.tdim + d)];
    };

    switch (ct)
    {
    case cell::type::interval:
    {
        const T xmin = std::min(xv(0, 0), xv(1, 0)) - tol;
        const T xmax = std::max(xv(0, 0), xv(1, 0)) + tol;
        return x_ref[0] >= xmin && x_ref[0] <= xmax;
    }
    case cell::type::triangle:
    {
        const T x0 = xv(0, 0), y0 = xv(0, 1);
        const T x1 = xv(1, 0), y1 = xv(1, 1);
        const T x2 = xv(2, 0), y2 = xv(2, 1);
        const T det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
        if (std::abs(det) <= tol)
            return false;
        const T l1 = ((x_ref[0] - x0) * (y2 - y0) - (x2 - x0) * (x_ref[1] - y0)) / det;
        const T l2 = ((x1 - x0) * (x_ref[1] - y0) - (x_ref[0] - x0) * (y1 - y0)) / det;
        const T l0 = T(1) - l1 - l2;
        return l0 >= -tol && l1 >= -tol && l2 >= -tol;
    }
    case cell::type::quadrilateral:
    case cell::type::hexahedron:
    {
        for (int d = 0; d < mesh.tdim; ++d)
        {
            T xmin = xv(0, d);
            T xmax = xv(0, d);
            for (int j = 1; j < mesh.cell_num_vertices(cell_id); ++j)
            {
                xmin = std::min(xmin, xv(j, d));
                xmax = std::max(xmax, xv(j, d));
            }
            if (x_ref[static_cast<std::size_t>(d)] < xmin - tol
                || x_ref[static_cast<std::size_t>(d)] > xmax + tol)
                return false;
        }
        return true;
    }
    case cell::type::tetrahedron:
    {
        const T x0 = xv(0, 0), y0 = xv(0, 1), z0 = xv(0, 2);
        const T x1 = xv(1, 0), y1 = xv(1, 1), z1 = xv(1, 2);
        const T x2 = xv(2, 0), y2 = xv(2, 1), z2 = xv(2, 2);
        const T x3 = xv(3, 0), y3 = xv(3, 1), z3 = xv(3, 2);

        const T j00 = x1 - x0, j01 = x2 - x0, j02 = x3 - x0;
        const T j10 = y1 - y0, j11 = y2 - y0, j12 = y3 - y0;
        const T j20 = z1 - z0, j21 = z2 - z0, j22 = z3 - z0;
        const T det = j00 * (j11 * j22 - j12 * j21)
                    - j01 * (j10 * j22 - j12 * j20)
                    + j02 * (j10 * j21 - j11 * j20);
        if (std::abs(det) <= tol)
            return false;

        const T rx = x_ref[0] - x0;
        const T ry = x_ref[1] - y0;
        const T rz = x_ref[2] - z0;
        const T l1 = (rx * (j11 * j22 - j12 * j21)
                    - j01 * (ry * j22 - j12 * rz)
                    + j02 * (ry * j21 - j11 * rz)) / det;
        const T l2 = (j00 * (ry * j22 - j12 * rz)
                    - rx * (j10 * j22 - j12 * j20)
                    + j02 * (j10 * rz - ry * j20)) / det;
        const T l3 = (j00 * (j11 * rz - ry * j21)
                    - j01 * (j10 * rz - ry * j20)
                    + rx * (j10 * j21 - j11 * j20)) / det;
        const T l0 = T(1) - l1 - l2 - l3;
        return l0 >= -tol && l1 >= -tol && l2 >= -tol && l3 >= -tol;
    }
    default:
        return false;
    }
}

template <std::floating_point T>
inline bool mark_cells_for_polynomial_backend(
    const LocalMesh<T>&         mesh,
    const std::vector<uint8_t>& value_available,
    int                         level_set_id,
    T                           tol,
    bool                        mark_uncertain,
    std::vector<uint8_t>&       marked_cells)
{
    const int nc = mesh.n_cells();
    marked_cells.assign(static_cast<std::size_t>(nc), 0);
    // Refinement policy:
    // 1) Default to edge/topology-driven refinement.
    // 2) Do not refine solely from broad interior-sign contradictions by default,
    //    because that over-refines neighboring same-sign cells.
    constexpr bool use_interior_contradiction_refinement = false;
    constexpr bool debug_mark_reasons = false;
    enum class CellMarkReason : uint8_t
    {
        not_marked = 0,
        multi_cross = 1,
        unresolved_zero_topology = 2,
        zero_edge_plus_extra_root = 3,
        zero_vertex_plus_same_edge_root = 4,
        interior_sign_contradiction = 5
    };
    std::vector<uint8_t> mark_reason;
    if constexpr (debug_mark_reasons)
        mark_reason.assign(
            static_cast<std::size_t>(nc),
            static_cast<uint8_t>(CellMarkReason::not_marked));

    for (int c = 0; c < nc; ++c)
    {
        uint8_t reason = static_cast<uint8_t>(CellMarkReason::not_marked);
        bool has_uncertain_edge = false;
        int n_single_cross_edges = 0;
        int n_zero_edge_edges = 0;
        bool one_root_shares_zero_endpoint = false;
        auto emit_debug = [&](int vertex_sign_dbg, bool mixed_dbg, bool has_zero_dbg)
        {
            if constexpr (debug_mark_reasons)
            {
                mark_reason[static_cast<std::size_t>(c)] = reason;
                std::printf(
                    "[poly-mark] c=%d marked=%u reason=%u v_sign=%d mixed=%d has_zero=%d single=%d zero_edge=%d uncertain=%d\n",
                    c,
                    static_cast<unsigned>(marked_cells[static_cast<std::size_t>(c)]),
                    static_cast<unsigned>(mark_reason[static_cast<std::size_t>(c)]),
                    vertex_sign_dbg,
                    mixed_dbg ? 1 : 0,
                    has_zero_dbg ? 1 : 0,
                    n_single_cross_edges,
                    n_zero_edge_edges,
                    has_uncertain_edge ? 1 : 0);
            }
        };
        const int e0 = mesh.cell_edge_offsets[static_cast<std::size_t>(c)];
        const int e1 = mesh.cell_edge_offsets[static_cast<std::size_t>(c + 1)];
        for (int j = e0; j < e1; ++j)
        {
            const int eid = mesh.cell_edges_flat[static_cast<std::size_t>(j)];
            const auto st = static_cast<EdgeState>(
                mesh.edge_state_for(eid, level_set_id));
            if (st == EdgeState::multi_cross)
            {
                marked_cells[static_cast<std::size_t>(c)] = 1;
                reason = static_cast<uint8_t>(CellMarkReason::multi_cross);
                break;
            }
            if (mark_uncertain && st == EdgeState::uncertain)
                has_uncertain_edge = true;
            if (st == EdgeState::single_cross)
            {
                ++n_single_cross_edges;
                const int v0 = mesh.edge_vertices[static_cast<std::size_t>(2 * eid)];
                const int v1 = mesh.edge_vertices[static_cast<std::size_t>(2 * eid + 1)];
                const uint64_t mask = (uint64_t(1) << level_set_id);
                const bool zero0 = (mesh.vertex_zero_mask[static_cast<std::size_t>(v0)] & mask) != 0;
                const bool zero1 = (mesh.vertex_zero_mask[static_cast<std::size_t>(v1)] & mask) != 0;
                if (zero0 || zero1)
                    one_root_shares_zero_endpoint = true;
            }
            else if (st == EdgeState::zero_edge)
            {
                ++n_zero_edge_edges;
            }
        }

        if (marked_cells[static_cast<std::size_t>(c)])
        {
            emit_debug(0, false, false);
            continue;
        }

        int vertex_sign = 0;
        bool mixed = false;
        bool has_zero = false;
        const int c0 = mesh.cell_offsets[static_cast<std::size_t>(c)];
        const int c1 = mesh.cell_offsets[static_cast<std::size_t>(c + 1)];
        for (int j = c0; j < c1; ++j)
        {
            const int vi = mesh.cell_vertices[static_cast<std::size_t>(j)];
            if (!value_available[static_cast<std::size_t>(vi)])
                continue;
            const T vv = mesh.vertex_phi[static_cast<std::size_t>(vi * mesh.n_level_sets + level_set_id)];
            const int sv = sign_with_tol(vv, tol);
            if (sv == 0)
            {
                has_zero = true;
                continue;
            }
            if (vertex_sign == 0)
                vertex_sign = sv;
            else if (vertex_sign != sv)
            {
                mixed = true;
                break;
            }
        }

        const bool pure_zero_touch
            = (n_zero_edge_edges > 0 && n_single_cross_edges == 0);
        const bool zero_vertex_plus_single_root_ok
            = has_zero
              && n_zero_edge_edges == 0
              && n_single_cross_edges == 1
              && !one_root_shares_zero_endpoint;
        const bool zero_edge_plus_extra_root
            = (n_zero_edge_edges > 0 && n_single_cross_edges > 0);
        const bool zero_vertex_plus_same_edge_root
            = one_root_shares_zero_endpoint;
        const bool unresolved_without_interface_entity
            = has_uncertain_edge
              && n_single_cross_edges == 0
              && n_zero_edge_edges == 0;

        if (zero_edge_plus_extra_root || zero_vertex_plus_same_edge_root)
        {
            marked_cells[static_cast<std::size_t>(c)] = 1;
            reason = zero_edge_plus_extra_root
                         ? static_cast<uint8_t>(CellMarkReason::zero_edge_plus_extra_root)
                         : static_cast<uint8_t>(CellMarkReason::zero_vertex_plus_same_edge_root);
            emit_debug(vertex_sign, mixed, has_zero);
            continue;
        }

        // Zero-touch cases are acceptable without further refinement:
        // - a pure zero edge/face with no additional root
        // - a zero vertex plus one separate root edge
        if (unresolved_without_interface_entity
            && !pure_zero_touch
            && !zero_vertex_plus_single_root_ok
            && (mixed || vertex_sign == 0))
        {
            marked_cells[static_cast<std::size_t>(c)] = 1;
            reason = static_cast<uint8_t>(CellMarkReason::unresolved_zero_topology);
            emit_debug(vertex_sign, mixed, has_zero);
            continue;
        }

        if (use_interior_contradiction_refinement && !mixed && vertex_sign != 0)
        {
            const int nv = mesh.n_vertices();
            for (int vi = 0; vi < nv; ++vi)
            {
                if (mesh.vertex_parent_dim[static_cast<std::size_t>(vi)] != mesh.tdim)
                    continue;
                if (!value_available[static_cast<std::size_t>(vi)])
                    continue;
                const std::span<const T> x_ref(
                    mesh.vertex_ref_x.data() + static_cast<std::size_t>(vi * mesh.tdim),
                    static_cast<std::size_t>(mesh.tdim));
                if (!point_in_local_cell_ref(mesh, c, x_ref, tol))
                    continue;
                const T v = mesh.vertex_phi[static_cast<std::size_t>(vi * mesh.n_level_sets + level_set_id)];
                const int si = sign_with_tol(v, tol);
                if (si != 0 && si != vertex_sign)
                {
                    marked_cells[static_cast<std::size_t>(c)] = 1;
                    reason = static_cast<uint8_t>(CellMarkReason::interior_sign_contradiction);
                    break;
                }
            }
        }

        emit_debug(vertex_sign, mixed, has_zero);
    }

    for (const uint8_t m : marked_cells)
    {
        if (m != 0)
            return true;
    }
    return false;
}

template <std::floating_point T, std::integral I = int>
inline bool classify_edges_from_parent_polynomial(
    LocalMesh<T>&                       mesh,
    const LevelSetFunction<T, I>&       level_set,
    const ParentPolynomialContext<T>&   parent_poly,
    int                                 level_set_id,
    T                                   tol,
    bool                                refine_on_uncertain,
    std::vector<uint8_t>*               marked_cells = nullptr)
{
    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("classify_edges_from_parent_polynomial: invalid level_set_id");
    if (mesh.vertex_ref_x.size() != static_cast<std::size_t>(mesh.n_vertices() * mesh.tdim))
        throw std::invalid_argument("classify_edges_from_parent_polynomial: Bernstein backend requires reference coordinates");

    std::vector<uint8_t> value_available;
    populate_local_level_set_values_from_parent_polynomial(
        mesh, parent_poly, level_set, level_set_id, tol, value_available);
    const T root_endpoint_tol = std::max(T(32) * tol, T(1e-10));

    const int ne = mesh.n_edges();
    for (int e = 0; e < ne; ++e)
    {
        const int v0 = mesh.edge_vertices[static_cast<std::size_t>(2 * e)];
        const int v1 = mesh.edge_vertices[static_cast<std::size_t>(2 * e + 1)];
        if (!value_available[static_cast<std::size_t>(v0)]
            || !value_available[static_cast<std::size_t>(v1)])
        {
            mesh.edge_state_for(e, level_set_id) = static_cast<uint8_t>(EdgeState::uncertain);
            mesh.edge_cert_for(e, level_set_id)
                = static_cast<uint8_t>(EdgeCertification::unresolved);
            continue;
        }

        const std::span<const T> x0_ref(
            mesh.vertex_ref_x.data() + static_cast<std::size_t>(v0 * mesh.tdim),
            static_cast<std::size_t>(mesh.tdim));
        const std::span<const T> x1_ref(
            mesh.vertex_ref_x.data() + static_cast<std::size_t>(v1 * mesh.tdim),
            static_cast<std::size_t>(mesh.tdim));

        std::vector<T> edge_coeffs;
        restrict_bernstein_to_segment_1d<T>(parent_poly.bernstein, x0_ref, x1_ref, edge_coeffs);
        const auto raw = classify_edge_from_full_cell_bernstein<T>(
            std::span<const T>(edge_coeffs.data(), edge_coeffs.size()), tol);
        const auto interior_samples = sample_interior_bernstein_edge<T>(
            std::span<const T>(edge_coeffs.data(), edge_coeffs.size()));
        const T ev0 = mesh.vertex_phi[static_cast<std::size_t>(v0 * mesh.n_level_sets + level_set_id)];
        const T ev1 = mesh.vertex_phi[static_cast<std::size_t>(v1 * mesh.n_level_sets + level_set_id)];
        const auto topo = classify_edge_zero_topology<T>(
            ev0, ev1,
            std::span<const std::pair<T, T>>(interior_samples.data(), interior_samples.size()),
            raw.state, tol, root_endpoint_tol);
        mesh.edge_state_for(e, level_set_id) = static_cast<uint8_t>(topo.state);
        mesh.edge_cert_for(e, level_set_id) = static_cast<uint8_t>(raw.cert);
    }

    std::vector<uint8_t> local_marks;
    const bool refine = mark_cells_for_polynomial_backend(
        mesh, value_available, level_set_id, tol, refine_on_uncertain, local_marks);
    if (marked_cells)
        *marked_cells = std::move(local_marks);
    return refine;
}

template <std::floating_point T, std::integral I = int>
inline bool classify_edges_from_local_level_set(
    LocalMesh<T>&                       mesh,
    const LocalLevelSetFunction<T, I>&  level_set,
    int                                 level_set_id,
    T                                   tol,
    bool                                refine_on_uncertain,
    std::vector<uint8_t>*               marked_cells = nullptr)
{
    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("classify_edges_from_local_level_set: invalid level_set_id");
    if (mesh.vertex_ref_x.size() != static_cast<std::size_t>(mesh.n_vertices() * mesh.tdim))
        throw std::invalid_argument(
            "classify_edges_from_local_level_set: Bernstein backend requires reference coordinates");

    std::vector<uint8_t> value_available;
    populate_local_level_set_values(mesh, level_set, level_set_id, tol, value_available);
    const T root_endpoint_tol = std::max(T(32) * tol, T(1e-10));

    const int ne = mesh.n_edges();
    for (int e = 0; e < ne; ++e)
    {
        const int v0 = mesh.edge_vertices[static_cast<std::size_t>(2 * e)];
        const int v1 = mesh.edge_vertices[static_cast<std::size_t>(2 * e + 1)];
        if (!value_available[static_cast<std::size_t>(v0)]
            || !value_available[static_cast<std::size_t>(v1)])
        {
            mesh.edge_state_for(e, level_set_id) = static_cast<uint8_t>(EdgeState::uncertain);
            mesh.edge_cert_for(e, level_set_id)
                = static_cast<uint8_t>(EdgeCertification::unresolved);
            continue;
        }

        std::vector<T> edge_coeffs;
        if (level_set.has_segment_restriction())
        {
            const std::span<const T> x0_ref(
                mesh.vertex_ref_x.data() + static_cast<std::size_t>(v0 * mesh.tdim),
                static_cast<std::size_t>(mesh.tdim));
            const std::span<const T> x1_ref(
                mesh.vertex_ref_x.data() + static_cast<std::size_t>(v1 * mesh.tdim),
                static_cast<std::size_t>(mesh.tdim));
            level_set.segment_restriction(x0_ref, x1_ref, edge_coeffs);
        }
        else if (level_set.has_edge_restriction()
                 && mesh.edge_parent_dim[static_cast<std::size_t>(e)] == 1
                 && mesh.edge_parent_id[static_cast<std::size_t>(e)] >= 0)
        {
            edge_coeffs = level_set.edge_restriction(
                mesh.edge_parent_id[static_cast<std::size_t>(e)]).coeffs;
        }
        else
        {
            mesh.edge_state_for(e, level_set_id) = static_cast<uint8_t>(EdgeState::uncertain);
            mesh.edge_cert_for(e, level_set_id)
                = static_cast<uint8_t>(EdgeCertification::unresolved);
            continue;
        }

        const auto raw = classify_edge_from_full_cell_bernstein<T>(
            std::span<const T>(edge_coeffs.data(), edge_coeffs.size()), tol);
        const auto interior_samples = sample_interior_bernstein_edge<T>(
            std::span<const T>(edge_coeffs.data(), edge_coeffs.size()));
        const T ev0 = mesh.vertex_phi[static_cast<std::size_t>(v0 * mesh.n_level_sets + level_set_id)];
        const T ev1 = mesh.vertex_phi[static_cast<std::size_t>(v1 * mesh.n_level_sets + level_set_id)];
        const auto topo = classify_edge_zero_topology<T>(
            ev0, ev1,
            std::span<const std::pair<T, T>>(interior_samples.data(), interior_samples.size()),
            raw.state, tol, root_endpoint_tol);
        mesh.edge_state_for(e, level_set_id) = static_cast<uint8_t>(topo.state);
        mesh.edge_cert_for(e, level_set_id) = static_cast<uint8_t>(raw.cert);
    }

    std::vector<uint8_t> local_marks;
    const bool refine = mark_cells_for_polynomial_backend(
        mesh, value_available, level_set_id, tol, refine_on_uncertain, local_marks);
    if (marked_cells)
        *marked_cells = std::move(local_marks);
    return refine;
}

template <std::floating_point T, std::integral I = int>
bool classify_edges_with_backend(
    LocalMesh<T>&                 mesh,
    const LevelSetFunction<T, I>& level_set,
    LocalLevelSetBackend          backend,
    int                           level_set_id = 0,
    T                             tol = static_cast<T>(1e-14),
    bool                          refine_on_uncertain = false)
{
    if (backend != LocalLevelSetBackend::bernstein)
        return classify_edges_and_mark_refine(mesh, level_set, level_set_id, tol, refine_on_uncertain);

    const auto local_level_set
        = level_set.mesh != nullptr
              ? make_local_level_set_function_bernstein<T, I>(
                    *level_set.mesh, level_set, static_cast<I>(mesh.parent_cell_id))
              : make_local_level_set_function_bernstein<T, I>(
                    mesh.parent_cell_type, mesh.gdim, level_set,
                    static_cast<I>(mesh.parent_cell_id));
    return classify_edges_from_local_level_set(
        mesh, local_level_set, level_set_id, tol, refine_on_uncertain, nullptr);
}

/// Classify all local-mesh edges from endpoint/interior interpolation-node signs.
///
/// Behavior:
/// - Uses `level_set.nodal_values` if available and its size matches mesh vertices.
/// - Otherwise falls back to `level_set.value(...)`.
/// - If neither source can provide values, the corresponding edge is `uncertain`.
/// - If any edge is `multi_cross`, `refine` is set true.
/// - If a cell-interior interpolation node has a sign that conflicts with a uniform
///   corner-vertex sign, `refine` is set true.
///
/// Returns whether the cell should be refined.
template <std::floating_point T, std::integral I = int>
bool classify_edges_and_mark_refine(
    LocalMesh<T>&                  mesh,
    const LevelSetFunction<T, I>&  level_set,
    int                            level_set_id,
    T                              tol,
    bool                           refine_on_uncertain)
{
    const int nv = mesh.n_vertices();
    const int ne = mesh.n_edges();

    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("classify_edges_and_mark_refine: invalid level_set_id");

    const bool use_nodal = level_set.has_nodal_values()
                           && static_cast<int>(level_set.nodal_values.size()) == nv;
    const bool use_value = level_set.has_value();

    std::vector<uint8_t> value_available(static_cast<std::size_t>(nv), 0);

    if (mesh.vertex_phi.size() != static_cast<std::size_t>(nv * mesh.n_level_sets))
        mesh.vertex_phi.assign(static_cast<std::size_t>(nv * mesh.n_level_sets), T(0));
    if (mesh.vertex_zero_mask.size() != static_cast<std::size_t>(nv))
        mesh.vertex_zero_mask.assign(static_cast<std::size_t>(nv), 0);
    if (mesh.vertex_inside_mask.size() != static_cast<std::size_t>(nv))
        mesh.vertex_inside_mask.assign(static_cast<std::size_t>(nv), 0);

    // Populate nodal values and masks for this level set.
    for (int i = 0; i < nv; ++i)
    {
        T v = T(0);
        bool ok = false;
        if (use_nodal)
        {
            v = level_set.nodal_values[static_cast<std::size_t>(i)];
            ok = true;
        }
        else if (use_value)
        {
            const T* x = mesh.vertex_x.data() + static_cast<std::size_t>(i * mesh.gdim);
            v = level_set.value(x, static_cast<I>(mesh.parent_cell_id));
            ok = true;
        }

        if (!ok)
            continue;

        value_available[static_cast<std::size_t>(i)] = 1;
        mesh.vertex_phi[static_cast<std::size_t>(i * mesh.n_level_sets + level_set_id)] = v;

        const uint64_t mask = (uint64_t(1) << level_set_id);
        mesh.vertex_zero_mask[static_cast<std::size_t>(i)] &= ~mask;
        mesh.vertex_inside_mask[static_cast<std::size_t>(i)] &= ~mask;
        const int s = sign_with_tol(v, tol);
        if (s == 0)
            mesh.vertex_zero_mask[static_cast<std::size_t>(i)] |= mask;
        else if (s < 0)
            mesh.vertex_inside_mask[static_cast<std::size_t>(i)] |= mask;
    }

    bool refine = false;
    const T root_endpoint_tol = std::max(T(32) * tol, T(1e-10));

    // Edge classification.
    for (int e = 0; e < ne; ++e)
    {
        const int v0 = mesh.edge_vertices[static_cast<std::size_t>(2 * e)];
        const int v1 = mesh.edge_vertices[static_cast<std::size_t>(2 * e + 1)];

        if (!value_available[static_cast<std::size_t>(v0)]
            || !value_available[static_cast<std::size_t>(v1)])
        {
            mesh.edge_state_for(e, level_set_id) = static_cast<uint8_t>(EdgeState::uncertain);
            mesh.edge_cert_for(e, level_set_id)
                = static_cast<uint8_t>(EdgeCertification::unresolved);
            if (refine_on_uncertain)
                refine = true;
            continue;
        }

        std::vector<std::pair<T, T>> interior_samples;
        interior_samples.reserve(8);
        std::vector<T> interior_vals;
        interior_vals.reserve(8);

        const bool use_ref = (mesh.tdim > 0)
                             && (mesh.vertex_ref_x.size() == static_cast<std::size_t>(nv * mesh.tdim));
        const int edim = use_ref ? mesh.tdim : mesh.gdim;
        const auto coord_ptr = [&](int vi) -> const T*
        {
            if (use_ref)
                return mesh.vertex_ref_x.data() + static_cast<std::size_t>(vi * mesh.tdim);
            return mesh.vertex_x.data() + static_cast<std::size_t>(vi * mesh.gdim);
        };
        const T* x0 = coord_ptr(v0);
        const T* x1 = coord_ptr(v1);
        T den = T(0);
        for (int d = 0; d < edim; ++d)
        {
            const T dx = x1[d] - x0[d];
            den += dx * dx;
        }

        for (int vi = 0; vi < nv; ++vi)
        {
            if (vi == v0 || vi == v1)
                continue;
            if (!value_available[static_cast<std::size_t>(vi)])
                continue;
            // Skip root vertices — they lie on edges by construction and their
            // near-zero phi values would corrupt the Bernstein classification.
            if (!mesh.vertex_root_edge_id.empty()
                && mesh.vertex_root_edge_id[static_cast<std::size_t>(vi)] >= 0)
                continue;

            const T* xi = coord_ptr(vi);
            T t = T(0.5);
            bool on_edge = false;
            if (den > std::numeric_limits<T>::epsilon())
            {
                T num = T(0);
                for (int d = 0; d < edim; ++d)
                    num += (xi[d] - x0[d]) * (x1[d] - x0[d]);
                t = num / den;

                T dist2 = T(0);
                for (int d = 0; d < edim; ++d)
                {
                    const T diff = xi[d] - (x0[d] + t * (x1[d] - x0[d]));
                    dist2 += diff * diff;
                }
                on_edge = (t > tol && t < T(1) - tol && dist2 <= T(64) * tol * tol);
            }

            if (on_edge)
            {
                const T val = mesh.vertex_phi[static_cast<std::size_t>(vi * mesh.n_level_sets + level_set_id)];
                interior_vals.push_back(val);
                interior_samples.push_back({t, val});
            }
        }

        const T ev0 = mesh.vertex_phi[static_cast<std::size_t>(v0 * mesh.n_level_sets + level_set_id)];
        const T ev1 = mesh.vertex_phi[static_cast<std::size_t>(v1 * mesh.n_level_sets + level_set_id)];
        const EdgeState nodal_st = classify_edge_from_nodal_samples(
            ev0, ev1, std::span<const T>(interior_vals.data(), interior_vals.size()), tol);

        EdgeState bernstein_st = EdgeState::uncertain;
        if (!interior_samples.empty())
        {
            std::sort(interior_samples.begin(), interior_samples.end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });
            std::vector<T> edge_samples;
            edge_samples.reserve(interior_samples.size() + 2);
            edge_samples.push_back(ev0);
            for (const auto& tv : interior_samples)
                edge_samples.push_back(tv.second);
            edge_samples.push_back(ev1);

            const auto bernstein_coeffs = bernstein_coefficients_from_equispaced_nodal<T>(
                std::span<const T>(edge_samples.data(), edge_samples.size()));
            const bool endpoint_sign_change = sign_with_tol(ev0, tol) != sign_with_tol(ev1, tol);
            bernstein_st = classify_edge_from_bernstein_coeffs<T>(
                std::span<const T>(bernstein_coeffs.data(), bernstein_coeffs.size()),
                endpoint_sign_change,
                tol);
        }
        const EdgeState raw_st = interior_samples.empty()
            ? nodal_st
            : combine_nodal_and_bernstein(nodal_st, bernstein_st);
        const auto topo = classify_edge_zero_topology<T>(
            ev0, ev1,
            std::span<const std::pair<T, T>>(interior_samples.data(), interior_samples.size()),
            raw_st, tol, root_endpoint_tol);
        const EdgeState st = topo.state;
        mesh.edge_state_for(e, level_set_id) = static_cast<uint8_t>(st);
        mesh.edge_cert_for(e, level_set_id)
            = static_cast<uint8_t>(EdgeCertification::resolved_by_sampling);

        if (st == EdgeState::multi_cross)
            refine = true;
        if (refine_on_uncertain && st == EdgeState::uncertain)
            refine = true;
    }

    // Interior interpolation-node test:
    // if corner vertices of a local cell have a uniform sign, interior nodes must
    // keep the same sign; otherwise we mark for refinement.
    for (int c = 0; c < mesh.n_cells(); ++c)
    {
        const int c0 = mesh.cell_offsets[static_cast<std::size_t>(c)];
        const int c1 = mesh.cell_offsets[static_cast<std::size_t>(c + 1)];

        int vertex_sign = 0;
        bool mixed = false;
        for (int j = c0; j < c1; ++j)
        {
            const int vi = mesh.cell_vertices[static_cast<std::size_t>(j)];
            if (!value_available[static_cast<std::size_t>(vi)])
                continue;
            const T vv = mesh.vertex_phi[static_cast<std::size_t>(vi * mesh.n_level_sets + level_set_id)];
            const int sv = sign_with_tol(vv, tol);
            if (sv == 0)
                continue;
            if (vertex_sign == 0)
                vertex_sign = sv;
            else if (vertex_sign != sv)
            {
                mixed = true;
                break;
            }
        }
        if (mixed || vertex_sign == 0)
            continue;

        for (int vi = 0; vi < nv; ++vi)
        {
            if (mesh.vertex_parent_dim[static_cast<std::size_t>(vi)] != mesh.tdim)
                continue;
            if (!value_available[static_cast<std::size_t>(vi)])
                continue;
            const T v = mesh.vertex_phi[static_cast<std::size_t>(vi * mesh.n_level_sets + level_set_id)];
            const int si = sign_with_tol(v, tol);
            if (si != 0 && si != vertex_sign)
            {
                refine = true;
                break;
            }
        }
        if (refine)
            break;
    }

    return refine;
}

/// Classify edges and, if required, refine the local mesh using a refinement template.
///
/// Returns true if refinement was triggered.
template <std::floating_point T, std::integral I = int>
bool classify_edges_and_refine_cell(
    LocalMesh<T>&                  mesh,
    const LevelSetFunction<T, I>&  level_set,
    const RefinementTemplate&      refine_template,
    int                            level_set_id = 0,
    T                              tol = static_cast<T>(1e-14),
    bool                           refine_on_uncertain = false)
{
    const bool do_refine = classify_edges_and_mark_refine(
        mesh, level_set, level_set_id, tol, refine_on_uncertain);
    if (do_refine)
    {
        refine_local_mesh_from_template(mesh, refine_template);
    }
    return do_refine;
}

template <std::floating_point T, std::integral I = int>
struct EdgeRefinementResult
{
    bool converged = false;
    bool has_ambiguous_edges = true;
    int iterations = 0;
};

/// Repeatedly classify edges and refine local mesh until ambiguous edge states
/// are eliminated or max_iterations is reached.
///
/// "Converged" means no edge is in {multi_cross, uncertain}.
template <std::floating_point T, std::integral I = int>
EdgeRefinementResult<T, I> refine_until_edges_single_root(
    LocalMesh<T>&                  mesh,
    const LevelSetFunction<T, I>&  level_set,
    const RefinementTemplate&      refine_template,
    int                            level_set_id = 0,
    T                              tol = static_cast<T>(1e-14),
    int                            max_iterations = 8,
    bool                           refine_on_uncertain = false)
{
    if (max_iterations <= 0)
        throw std::invalid_argument("refine_until_edges_single_root: max_iterations must be > 0");

    EdgeRefinementResult<T, I> out;
    for (int it = 0; it < max_iterations; ++it)
    {
        const bool do_refine = classify_edges_and_mark_refine(
            mesh, level_set, level_set_id, tol, refine_on_uncertain);
        out.iterations = it + 1;

        bool has_ambiguous = false;
        for (const uint8_t st_u8 : mesh.edge_state)
        {
            const auto st = static_cast<EdgeState>(st_u8);
            if (st == EdgeState::multi_cross
                || st == EdgeState::uncertain)
            {
                has_ambiguous = true;
                break;
            }
        }

        if (!has_ambiguous)
        {
            out.converged = true;
            out.has_ambiguous_edges = false;
            return out;
        }
        if (!do_refine)
        {
            out.converged = false;
            out.has_ambiguous_edges = true;
            return out;
        }
        refine_local_mesh_from_template(mesh, refine_template);
    }

    out.converged = false;
    out.has_ambiguous_edges = true;
    return out;
}

template <std::floating_point T, std::integral I = int>
struct PolynomialCertificationResult
{
    bool certified = false;
    bool has_marked_cells = true;
    int iterations = 0;
};

template <std::floating_point T, std::integral I = int>
PolynomialCertificationResult<T, I> certify_local_mesh_polynomial(
    LocalMesh<T>&                  mesh,
    const LevelSetFunction<T, I>&  level_set,
    int                            level_set_id = 0,
    int                            max_refine_levels = 8,
    T                              tol = static_cast<T>(1e-14),
    bool                           refine_on_uncertain = true)
{
    if (max_refine_levels <= 0)
        throw std::invalid_argument("certify_local_mesh_polynomial: max_refine_levels must be > 0");

    const auto local_level_set
        = level_set.mesh != nullptr
              ? make_local_level_set_function_bernstein<T, I>(
                    *level_set.mesh, level_set, static_cast<I>(mesh.parent_cell_id))
              : make_local_level_set_function_bernstein<T, I>(
                    mesh.parent_cell_type, mesh.gdim, level_set,
                    static_cast<I>(mesh.parent_cell_id));

    PolynomialCertificationResult<T, I> out;
    for (int it = 0; it < max_refine_levels; ++it)
    {
        std::vector<uint8_t> marked_cells;
        const bool need_refine = classify_edges_from_local_level_set(
            mesh, local_level_set, level_set_id, tol, refine_on_uncertain, &marked_cells);
        out.iterations = it + 1;

        if (!need_refine)
        {
            out.certified = true;
            out.has_marked_cells = false;
            return out;
        }

        red_refine_marked_cells(mesh, std::span<const uint8_t>(marked_cells.data(), marked_cells.size()));
    }

    out.certified = false;
    out.has_marked_cells = true;
    return out;
}

} // namespace cutcells
