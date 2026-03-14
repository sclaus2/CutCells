// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include "local_mesh.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace cutcells
{

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

    if (s0 == 0 || s1 == 0)
        return EdgeState::near_tangency;

    for (const T vi : interior_values)
    {
        if (sign_with_tol(vi, tol) == 0)
            return EdgeState::near_tangency;
    }

    if (s0 != s1)
        return EdgeState::one_root;

    for (const T vi : interior_values)
    {
        const int si = sign_with_tol(vi, tol);
        if (si != 0 && si != s0)
            return EdgeState::multiple_roots;
    }

    return EdgeState::no_root;
}

template <std::floating_point T>
inline T binomial_coeff(int n, int k)
{
    if (k < 0 || k > n)
        return T(0);
    if (k == 0 || k == n)
        return T(1);
    if (k > n - k)
        k = n - k;
    T c = T(1);
    for (int i = 1; i <= k; ++i)
        c = c * static_cast<T>(n - (k - i)) / static_cast<T>(i);
    return c;
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
inline EdgeState classify_edge_from_bernstein_coeffs(
    std::span<const T> coeffs,
    bool endpoint_sign_change,
    T tol)
{
    if (coeffs.empty())
        return EdgeState::uncertain;

    int sign_changes = 0;
    int prev_sign = 0;
    for (const T c : coeffs)
    {
        const int sc = sign_with_tol(c, tol);
        if (sc == 0)
            return EdgeState::near_tangency;
        if (prev_sign != 0 && sc != prev_sign)
            ++sign_changes;
        prev_sign = sc;
    }

    if (sign_changes == 0)
        return EdgeState::no_root;
    if (sign_changes == 1 && endpoint_sign_change)
        return EdgeState::one_root;
    if (sign_changes >= 2)
        return EdgeState::multiple_roots;
    return EdgeState::uncertain;
}

inline EdgeState combine_nodal_and_bernstein(EdgeState nodal, EdgeState bernstein)
{
    if (nodal == EdgeState::multiple_roots || nodal == EdgeState::near_tangency
        || nodal == EdgeState::uncertain)
        return nodal;
    if (bernstein == EdgeState::multiple_roots || bernstein == EdgeState::near_tangency
        || bernstein == EdgeState::uncertain)
        return bernstein;
    if (nodal == bernstein)
        return nodal;
    return EdgeState::uncertain;
}

/// Classify all local-mesh edges from endpoint/interior interpolation-node signs.
///
/// Behavior:
/// - Uses `level_set.nodal_values` if available and its size matches mesh vertices.
/// - Otherwise falls back to `level_set.value(...)`.
/// - If neither source can provide values, the corresponding edge is `uncertain`.
/// - If any edge is `multiple_roots` or `near_tangency`, `refine` is set true.
/// - If a cell-interior interpolation node has a sign that conflicts with a uniform
///   corner-vertex sign, `refine` is set true.
///
/// Returns whether the cell should be refined.
template <std::floating_point T, std::integral I = int>
bool classify_edges_and_mark_refine(
    LocalMesh<T>&                  mesh,
    const LevelSetFunction<T, I>&  level_set,
    int                            level_set_id = 0,
    T                              tol = static_cast<T>(1e-14),
    bool                           refine_on_uncertain = false)
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

    // Edge classification.
    for (int e = 0; e < ne; ++e)
    {
        const int v0 = mesh.edge_vertices[static_cast<std::size_t>(2 * e)];
        const int v1 = mesh.edge_vertices[static_cast<std::size_t>(2 * e + 1)];

        if (!value_available[static_cast<std::size_t>(v0)]
            || !value_available[static_cast<std::size_t>(v1)])
        {
            mesh.edge_state[static_cast<std::size_t>(e)] = static_cast<uint8_t>(EdgeState::uncertain);
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
            if (mesh.vertex_parent_dim[static_cast<std::size_t>(vi)] == 1
                && mesh.vertex_parent_id[static_cast<std::size_t>(vi)] == e
                && value_available[static_cast<std::size_t>(vi)])
            {
                const T val = mesh.vertex_phi[static_cast<std::size_t>(vi * mesh.n_level_sets + level_set_id)];
                interior_vals.push_back(val);
                T t = T(0.5);
                if (den > std::numeric_limits<T>::epsilon())
                {
                    const T* xi = coord_ptr(vi);
                    T num = T(0);
                    for (int d = 0; d < edim; ++d)
                        num += (xi[d] - x0[d]) * (x1[d] - x0[d]);
                    t = num / den;
                }
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
        const EdgeState st = interior_samples.empty() ? nodal_st : combine_nodal_and_bernstein(nodal_st, bernstein_st);
        mesh.edge_state[static_cast<std::size_t>(e)] = static_cast<uint8_t>(st);

        if (st == EdgeState::multiple_roots || st == EdgeState::near_tangency)
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
            if (si == 0 || si != vertex_sign)
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
/// "Converged" means no edge is in {multiple_roots, near_tangency, uncertain}.
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
            if (st == EdgeState::multiple_roots
                || st == EdgeState::near_tangency
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

} // namespace cutcells
