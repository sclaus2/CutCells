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
            if (si == 0 || si != vertex_sign)
                return true;
        }
    }

    return false;
}

template <std::floating_point T>
inline EdgeState classify_edge_from_full_cell_bernstein(
    std::span<const T> edge_coeffs,
    T                  tol)
{
    if (edge_coeffs.empty())
        return EdgeState::uncertain;

    if (bernstein_all_strict_one_sign(edge_coeffs, tol))
        return EdgeState::no_root;
    if (bernstein_has_near_zero(edge_coeffs, tol))
        return EdgeState::near_tangency;

    const int sign_changes = bernstein_sign_variation_count(edge_coeffs, tol);
    const bool monotone = bernstein_derivative_strict_one_sign(edge_coeffs, tol);
    const int s0 = sign_with_tol(edge_coeffs.front(), tol);
    const int s1 = sign_with_tol(edge_coeffs.back(), tol);
    const bool endpoint_sign_change = (s0 != 0 && s1 != 0 && s0 != s1);

    if (monotone && endpoint_sign_change)
        return EdgeState::one_root;
    if (sign_changes >= 2)
        return EdgeState::multiple_roots;
    if (sign_changes == 0)
        return EdgeState::no_root;
    return EdgeState::uncertain;
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
    const int nv = mesh.n_vertices();
    marked_cells.assign(static_cast<std::size_t>(nc), 0);

    for (int c = 0; c < nc; ++c)
    {
        const int e0 = mesh.cell_edge_offsets[static_cast<std::size_t>(c)];
        const int e1 = mesh.cell_edge_offsets[static_cast<std::size_t>(c + 1)];
        for (int j = e0; j < e1; ++j)
        {
            const auto st = static_cast<EdgeState>(
                mesh.edge_state[static_cast<std::size_t>(mesh.cell_edges_flat[static_cast<std::size_t>(j)])]);
            if (st == EdgeState::multiple_roots || st == EdgeState::near_tangency
                || (mark_uncertain && st == EdgeState::uncertain))
            {
                marked_cells[static_cast<std::size_t>(c)] = 1;
                break;
            }
        }

        if (marked_cells[static_cast<std::size_t>(c)])
            continue;

        int vertex_sign = 0;
        bool mixed = false;
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
            const std::span<const T> x_ref(
                mesh.vertex_ref_x.data() + static_cast<std::size_t>(vi * mesh.tdim),
                static_cast<std::size_t>(mesh.tdim));
            if (!point_in_local_cell_ref(mesh, c, x_ref, tol))
                continue;
            const T v = mesh.vertex_phi[static_cast<std::size_t>(vi * mesh.n_level_sets + level_set_id)];
            const int si = sign_with_tol(v, tol);
            if (si == 0 || si != vertex_sign)
            {
                marked_cells[static_cast<std::size_t>(c)] = 1;
                break;
            }
        }
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

    const int ne = mesh.n_edges();
    for (int e = 0; e < ne; ++e)
    {
        const int v0 = mesh.edge_vertices[static_cast<std::size_t>(2 * e)];
        const int v1 = mesh.edge_vertices[static_cast<std::size_t>(2 * e + 1)];
        if (!value_available[static_cast<std::size_t>(v0)]
            || !value_available[static_cast<std::size_t>(v1)])
        {
            mesh.edge_state[static_cast<std::size_t>(e)] = static_cast<uint8_t>(EdgeState::uncertain);
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
        mesh.edge_state[static_cast<std::size_t>(e)] = static_cast<uint8_t>(
            classify_edge_from_full_cell_bernstein<T>(
                std::span<const T>(edge_coeffs.data(), edge_coeffs.size()), tol));
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

    const auto parent_poly = make_parent_polynomial_context<T, I>(
        mesh.parent_cell_type, level_set);
    return classify_edges_from_parent_polynomial(
        mesh, level_set, parent_poly, level_set_id, tol, refine_on_uncertain, nullptr);
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
            if (!value_available[static_cast<std::size_t>(vi)])
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

    const auto parent_poly = make_parent_polynomial_context<T, I>(
        mesh.parent_cell_type, level_set);

    PolynomialCertificationResult<T, I> out;
    for (int it = 0; it < max_refine_levels; ++it)
    {
        std::vector<uint8_t> marked_cells;
        const bool need_refine = classify_edges_from_parent_polynomial(
            mesh, level_set, parent_poly, level_set_id, tol, refine_on_uncertain, &marked_cells);
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
