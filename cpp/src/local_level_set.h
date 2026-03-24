// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include "bernstein_backend.h"
#include "cell_topology.h"
#include "mesh_view.h"

#include <concepts>
#include <cstdint>
#include <functional>
#include <memory>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cutcells
{

enum class LocalLevelSetKind : uint8_t
{
    callable = 0,
    bernstein = 1,
    taylor = 2
};

template <std::floating_point T>
struct LocalLevelSetEdgeRestriction
{
    int degree = -1;
    std::vector<T> coeffs;
};

template <std::floating_point T>
struct LocalBernsteinData
{
    BernsteinCell<T> cell_poly;
    cell::type parent_cell_type = cell::type::point;
    int gdim = 0;
    int tdim = 0;
    int degree = -1;
};

template <std::floating_point T, std::integral I = int>
struct LocalLevelSetFunction
{
    std::function<T(const T*)> value_fn;
    std::function<void(const T*, T*)> grad_fn;
    std::function<void(int, std::vector<T>&)> edge_restrict_fn;
    std::function<void(std::span<const T>, std::span<const T>, std::vector<T>&)>
        segment_restrict_fn;

    I parent_cell_id = static_cast<I>(-1);
    int gdim = 0;
    int tdim = 0;
    int degree = -1;
    int level_set_id = -1;
    LocalLevelSetKind backend = LocalLevelSetKind::callable;
    std::shared_ptr<void> owner;

    bool has_value() const
    {
        return static_cast<bool>(value_fn);
    }

    bool has_gradient() const
    {
        return static_cast<bool>(grad_fn);
    }

    bool has_edge_restriction() const
    {
        return static_cast<bool>(edge_restrict_fn);
    }

    bool has_segment_restriction() const
    {
        return static_cast<bool>(segment_restrict_fn);
    }

    T value(const T* x_ref) const
    {
        if (!value_fn)
            throw std::runtime_error("LocalLevelSetFunction::value not available");
        return value_fn(x_ref);
    }

    void grad(const T* x_ref, T* g_ref) const
    {
        if (!grad_fn)
            throw std::runtime_error("LocalLevelSetFunction::grad not available");
        grad_fn(x_ref, g_ref);
    }

    LocalLevelSetEdgeRestriction<T> edge_restriction(int edge_id) const
    {
        if (!edge_restrict_fn)
            throw std::runtime_error("LocalLevelSetFunction::edge_restriction not available");
        LocalLevelSetEdgeRestriction<T> out;
        out.degree = degree;
        edge_restrict_fn(edge_id, out.coeffs);
        return out;
    }

    void segment_restriction(
        std::span<const T> x0_ref,
        std::span<const T> x1_ref,
        std::vector<T>&    coeffs) const
    {
        if (!segment_restrict_fn)
            throw std::runtime_error("LocalLevelSetFunction::segment_restriction not available");
        segment_restrict_fn(x0_ref, x1_ref, coeffs);
    }
};

template <std::floating_point T, std::integral I>
inline std::vector<T> extract_cell_nodal_values(
    const MeshView<T, I>& mesh,
    std::span<const T>    global_nodal_values,
    I                     cell_id)
{
    const auto cell_nodes = mesh.cell_nodes(cell_id);
    std::vector<T> local_values(static_cast<std::size_t>(cell_nodes.size()), T(0));
    for (std::size_t i = 0; i < cell_nodes.size(); ++i)
    {
        const I node_id = cell_nodes[i];
        if (node_id < 0 || static_cast<std::size_t>(node_id) >= global_nodal_values.size())
            throw std::invalid_argument("extract_cell_nodal_values: cell node id out of range");
        local_values[i] = global_nodal_values[static_cast<std::size_t>(node_id)];
    }
    return local_values;
}

template <std::floating_point T>
inline T bernstein_basis_interval_derivative(int p, int i, T x)
{
    if (p <= 0)
        return T(0);

    T value = T(0);
    if (i > 0)
    {
        value += binomial_coeff<T>(p, i) * static_cast<T>(i)
                 * powi(x, i - 1) * powi(T(1) - x, p - i);
    }
    if (p - i > 0)
    {
        value -= binomial_coeff<T>(p, i) * static_cast<T>(p - i)
                 * powi(x, i) * powi(T(1) - x, p - i - 1);
    }
    return value;
}

template <std::floating_point T>
inline void evaluate_bernstein_cell_gradient(
    const BernsteinCell<T>& cell_poly,
    std::span<const T>      x_ref,
    std::span<T>            grad_ref)
{
    const int p = cell_poly.degree;
    if (cell_poly.tdim != static_cast<int>(x_ref.size()))
        throw std::invalid_argument(
            "evaluate_bernstein_cell_gradient: point dimension mismatch");
    if (cell_poly.tdim != static_cast<int>(grad_ref.size()))
        throw std::invalid_argument(
            "evaluate_bernstein_cell_gradient: gradient dimension mismatch");

    std::fill(grad_ref.begin(), grad_ref.end(), T(0));

    switch (cell_poly.cell_type)
    {
    case cell::type::interval:
        for (int i = 0; i <= p; ++i)
        {
            const int idx = interval_i_to_flat(p, i);
            grad_ref[0] += cell_poly.coeffs[static_cast<std::size_t>(idx)]
                           * bernstein_basis_interval_derivative<T>(p, i, x_ref[0]);
        }
        return;
    case cell::type::triangle:
    {
        const T l0 = T(1) - x_ref[0] - x_ref[1];
        const T l1 = x_ref[0];
        const T l2 = x_ref[1];
        for (int i = 0; i <= p; ++i)
        {
            for (int j = 0; j <= p - i; ++j)
            {
                const int k = p - i - j;
                const int idx = triangle_ijk_to_flat(p, i, j, k);
                const T coeff = cell_poly.coeffs[static_cast<std::size_t>(idx)]
                                * multinomial_coeff_triangle<T>(p, i, j, k);
                if (i > 0)
                {
                    const T term = coeff * powi(l0, i - 1) * powi(l1, j) * powi(l2, k);
                    grad_ref[0] -= static_cast<T>(i) * term;
                    grad_ref[1] -= static_cast<T>(i) * term;
                }
                if (j > 0)
                {
                    grad_ref[0] += coeff * static_cast<T>(j)
                                   * powi(l0, i) * powi(l1, j - 1) * powi(l2, k);
                }
                if (k > 0)
                {
                    grad_ref[1] += coeff * static_cast<T>(k)
                                   * powi(l0, i) * powi(l1, j) * powi(l2, k - 1);
                }
            }
        }
        return;
    }
    case cell::type::quadrilateral:
        for (int j = 0; j <= p; ++j)
        {
            for (int i = 0; i <= p; ++i)
            {
                const int idx = quad_ij_to_flat(p, i, j);
                const T coeff = cell_poly.coeffs[static_cast<std::size_t>(idx)];
                grad_ref[0] += coeff
                               * bernstein_basis_interval_derivative<T>(p, i, x_ref[0])
                               * bernstein_basis_interval<T>(p, j, x_ref[1]);
                grad_ref[1] += coeff
                               * bernstein_basis_interval<T>(p, i, x_ref[0])
                               * bernstein_basis_interval_derivative<T>(p, j, x_ref[1]);
            }
        }
        return;
    case cell::type::tetrahedron:
    {
        const T l0 = T(1) - x_ref[0] - x_ref[1] - x_ref[2];
        const T l1 = x_ref[0];
        const T l2 = x_ref[1];
        const T l3 = x_ref[2];
        for (int i = 0; i <= p; ++i)
        {
            for (int j = 0; j <= p - i; ++j)
            {
                for (int k = 0; k <= p - i - j; ++k)
                {
                    const int l = p - i - j - k;
                    const int idx = tetrahedron_ijkl_to_flat(p, i, j, k, l);
                    const T coeff = cell_poly.coeffs[static_cast<std::size_t>(idx)]
                                    * multinomial_coeff_tetrahedron<T>(p, i, j, k, l);
                    if (i > 0)
                    {
                        const T term = coeff * powi(l0, i - 1)
                                       * powi(l1, j) * powi(l2, k) * powi(l3, l);
                        grad_ref[0] -= static_cast<T>(i) * term;
                        grad_ref[1] -= static_cast<T>(i) * term;
                        grad_ref[2] -= static_cast<T>(i) * term;
                    }
                    if (j > 0)
                    {
                        grad_ref[0] += coeff * static_cast<T>(j)
                                       * powi(l0, i) * powi(l1, j - 1) * powi(l2, k) * powi(l3, l);
                    }
                    if (k > 0)
                    {
                        grad_ref[1] += coeff * static_cast<T>(k)
                                       * powi(l0, i) * powi(l1, j) * powi(l2, k - 1) * powi(l3, l);
                    }
                    if (l > 0)
                    {
                        grad_ref[2] += coeff * static_cast<T>(l)
                                       * powi(l0, i) * powi(l1, j) * powi(l2, k) * powi(l3, l - 1);
                    }
                }
            }
        }
        return;
    }
    case cell::type::hexahedron:
        for (int k = 0; k <= p; ++k)
        {
            for (int j = 0; j <= p; ++j)
            {
                for (int i = 0; i <= p; ++i)
                {
                    const int idx = hex_ijk_to_flat(p, i, j, k);
                    const T coeff = cell_poly.coeffs[static_cast<std::size_t>(idx)];
                    grad_ref[0] += coeff
                                   * bernstein_basis_interval_derivative<T>(p, i, x_ref[0])
                                   * bernstein_basis_interval<T>(p, j, x_ref[1])
                                   * bernstein_basis_interval<T>(p, k, x_ref[2]);
                    grad_ref[1] += coeff
                                   * bernstein_basis_interval<T>(p, i, x_ref[0])
                                   * bernstein_basis_interval_derivative<T>(p, j, x_ref[1])
                                   * bernstein_basis_interval<T>(p, k, x_ref[2]);
                    grad_ref[2] += coeff
                                   * bernstein_basis_interval<T>(p, i, x_ref[0])
                                   * bernstein_basis_interval<T>(p, j, x_ref[1])
                                   * bernstein_basis_interval_derivative<T>(p, k, x_ref[2]);
                }
            }
        }
        return;
    default:
        throw std::invalid_argument(
            "evaluate_bernstein_cell_gradient: unsupported cell type");
    }
}

template <std::floating_point T, std::integral I>
inline LocalLevelSetFunction<T, I> make_local_level_set_function_bernstein(
    cell::type         cell_type,
    int                gdim,
    std::span<const T> local_values,
    int                degree,
    I                  cell_id)
{
    if (degree < 1)
        degree = infer_lagrange_order_from_num_nodes(
            cell_type, static_cast<int>(local_values.size()));

    auto data = std::make_shared<LocalBernsteinData<T>>();
    data->cell_poly = make_bernstein_cell<T>(
        cell_type, degree, std::span<const T>(local_values.data(), local_values.size()));
    data->parent_cell_type = cell_type;
    data->gdim = gdim;
    data->tdim = cell::get_tdim(cell_type);
    data->degree = degree;

    LocalLevelSetFunction<T, I> local_phi;
    local_phi.parent_cell_id = cell_id;
    local_phi.gdim = gdim;
    local_phi.tdim = data->tdim;
    local_phi.degree = degree;
    local_phi.backend = LocalLevelSetKind::bernstein;
    local_phi.owner = data;
    local_phi.value_fn = [data](const T* x_ref) -> T
    {
        return evaluate_bernstein_cell<T>(
            data->cell_poly,
            std::span<const T>(x_ref, static_cast<std::size_t>(data->tdim)));
    };
    local_phi.grad_fn = [data](const T* x_ref, T* g_ref)
    {
        evaluate_bernstein_cell_gradient<T>(
            data->cell_poly,
            std::span<const T>(x_ref, static_cast<std::size_t>(data->tdim)),
            std::span<T>(g_ref, static_cast<std::size_t>(data->tdim)));
    };
    local_phi.edge_restrict_fn = [data](int edge_id, std::vector<T>& coeffs)
    {
        const auto ref_vertices = p1_ref_coords(data->parent_cell_type);
        const auto edge_vertices = cell::edges(data->parent_cell_type);
        if (edge_id < 0 || edge_id >= static_cast<int>(edge_vertices.size()))
            throw std::invalid_argument(
                "LocalLevelSetFunction::edge_restriction: invalid edge id");

        const int tdim = data->tdim;
        std::vector<T> x0(static_cast<std::size_t>(tdim), T(0));
        std::vector<T> x1(static_cast<std::size_t>(tdim), T(0));
        const int v0 = edge_vertices[static_cast<std::size_t>(edge_id)][0];
        const int v1 = edge_vertices[static_cast<std::size_t>(edge_id)][1];
        for (int d = 0; d < tdim; ++d)
        {
            x0[static_cast<std::size_t>(d)] = static_cast<T>(
                ref_vertices[static_cast<std::size_t>(v0 * tdim + d)]);
            x1[static_cast<std::size_t>(d)] = static_cast<T>(
                ref_vertices[static_cast<std::size_t>(v1 * tdim + d)]);
        }
        restrict_bernstein_to_segment_1d<T>(
            data->cell_poly,
            std::span<const T>(x0.data(), x0.size()),
            std::span<const T>(x1.data(), x1.size()),
            coeffs);
    };
    local_phi.segment_restrict_fn =
        [data](std::span<const T> x0_ref, std::span<const T> x1_ref, std::vector<T>& coeffs)
    {
        restrict_bernstein_to_segment_1d<T>(data->cell_poly, x0_ref, x1_ref, coeffs);
    };
    return local_phi;
}

template <std::floating_point T, std::integral I>
inline LocalLevelSetFunction<T, I> make_local_level_set_function_bernstein(
    cell::type                    cell_type,
    int                           gdim,
    const LevelSetFunction<T, I>& level_set,
    I                             cell_id)
{
    if (!level_set.has_nodal_values())
        throw std::invalid_argument(
            "make_local_level_set_function_bernstein: nodal values are required");

    return make_local_level_set_function_bernstein<T, I>(
        cell_type,
        gdim,
        std::span<const T>(level_set.nodal_values.data(), level_set.nodal_values.size()),
        level_set.degree,
        cell_id);
}

template <std::floating_point T, std::integral I>
inline LocalLevelSetFunction<T, I> make_local_level_set_function_bernstein(
    const MeshView<T, I>& mesh,
    const LevelSetFunction<T, I>& level_set,
    I cell_id)
{
    const auto cell_type = static_cast<cell::type>(mesh.cell_type(cell_id));
    const auto local_values = extract_cell_nodal_values(mesh, level_set.nodal_values, cell_id);
    return make_local_level_set_function_bernstein<T, I>(
        cell_type,
        mesh.gdim,
        std::span<const T>(local_values.data(), local_values.size()),
        level_set.degree,
        cell_id);
}

template <std::floating_point T, std::integral I>
inline LocalLevelSetFunction<T, I> local_level_set_from_callable(
    const LevelSetFunction<T, I>& level_set,
    I                             cell_id)
{
    LocalLevelSetFunction<T, I> local_phi;
    local_phi.parent_cell_id = cell_id;
    local_phi.gdim = level_set.gdim;
    local_phi.tdim = level_set.gdim;
    local_phi.degree = level_set.degree;
    local_phi.backend = LocalLevelSetKind::callable;
    local_phi.owner = level_set.owner;
    if (level_set.has_value())
    {
        auto value_fn = level_set.value_fn;
        local_phi.value_fn = [value_fn, cell_id](const T* x_ref) -> T
        {
            return value_fn(x_ref, cell_id);
        };
    }
    if (level_set.has_gradient())
    {
        auto grad_fn = level_set.grad_fn;
        const int gdim = level_set.gdim;
        local_phi.grad_fn = [grad_fn, cell_id, gdim](const T* x_ref, T* g_ref)
        {
            if (gdim <= 0)
                throw std::runtime_error("LocalLevelSetFunction callable grad: invalid gdim");
            grad_fn(x_ref, cell_id, g_ref);
        };
    }
    return local_phi;
}

template <std::floating_point T, std::integral I>
inline LocalLevelSetEdgeRestriction<T> local_edge_restriction(
    const LocalLevelSetFunction<T, I>& phi,
    int                                edge_id)
{
    return phi.edge_restriction(edge_id);
}

template <std::floating_point T, std::integral I>
inline LocalLevelSetFunction<T, I> LevelSetFunction<T, I>::local(I cell_id) const
{
    if (kind == Kind::fem_nodal)
    {
        if (mesh == nullptr)
            throw std::runtime_error("LevelSetFunction::local requires mesh in FEM mode");
        return make_local_level_set_function_bernstein(*mesh, *this, cell_id);
    }
    return local_level_set_from_callable(*this, cell_id);
}

} // namespace cutcells
