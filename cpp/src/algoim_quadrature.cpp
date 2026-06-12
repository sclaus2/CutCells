// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "algoim_quadrature.h"

#include "cell_topology.h"
#include "edge_root.h"
#include "ho_mesh_part_output.h"
#include "mapping.h"
#include "reference_cell.h"
#include "selection_expr.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <limits>
#include <span>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef CUTCELLS_HAS_ALGOIM
#include <algoim/quadrature_general.hpp>
#endif

#ifdef CUTCELLS_HAS_ALGOIM_MULTIPOLY
#include <algoim/quadrature_multipoly.hpp>
#endif

namespace cutcells::output
{
namespace
{

#ifdef CUTCELLS_HAS_ALGOIM

constexpr int max_algoim_bernstein_degree = 16;
constexpr int derivative_bound_cache_base = 4;

template <typename S>
using BasisArray = std::array<S, max_algoim_bernstein_degree + 1>;

inline void validate_algoim_bernstein_degree(int degree)
{
    if (degree < 0 || degree > max_algoim_bernstein_degree)
    {
        throw std::runtime_error(
            "Algoim quadrature currently supports Bernstein level-set degree <= "
            + std::to_string(max_algoim_bernstein_degree));
    }
}

template <typename S>
S integer_pow(S base, int exp)
{
    S result = S(1);
    for (int i = 0; i < exp; ++i)
        result *= base;
    return result;
}

template <typename S>
S binomial(int n, int k)
{
    if (k < 0 || k > n)
        return S(0);
    if (k == 0 || k == n)
        return S(1);
    if (k > n - k)
        k = n - k;

    S result = S(1);
    for (int i = 0; i < k; ++i)
    {
        result *= S(n - i);
        result /= S(i + 1);
    }
    return result;
}

template <typename S>
void bernstein_1d_all(int degree, S t, BasisArray<S>& out)
{
    validate_algoim_bernstein_degree(degree);
    const int n = degree;
    for (int i = 0; i <= n; ++i)
        out[static_cast<std::size_t>(i)] =
            binomial<S>(n, i) * integer_pow(t, i)
            * integer_pow(S(1) - t, n - i);
}

template <typename S>
void bernstein_1d_deriv(int degree, S t, BasisArray<S>& dbasis)
{
    validate_algoim_bernstein_degree(degree);
    const int n = degree;
    if (n == 0)
    {
        dbasis[0] = S(0);
        return;
    }

    BasisArray<S> b;
    bernstein_1d_all(n - 1, t, b);
    for (int i = 0; i <= n; ++i)
    {
        const S bim1 = (i > 0) ? b[static_cast<std::size_t>(i - 1)] : S(0);
        const S bi = (i < n) ? b[static_cast<std::size_t>(i)] : S(0);
        dbasis[static_cast<std::size_t>(i)] = S(n) * (bim1 - bi);
    }
}

template <int N, typename F>
void for_each_multi_index(const std::array<int, N>& upper, F&& f)
{
    if constexpr (N == 2)
    {
        for (int i = 0; i <= upper[0]; ++i)
            for (int j = 0; j <= upper[1]; ++j)
                f(std::array<int, N>{{i, j}});
    }
    else if constexpr (N == 3)
    {
        for (int i = 0; i <= upper[0]; ++i)
            for (int j = 0; j <= upper[1]; ++j)
                for (int k = 0; k <= upper[2]; ++k)
                    f(std::array<int, N>{{i, j, k}});
    }
}

template <std::floating_point T>
T falling_factorial(int n, int k)
{
    T value = T(1);
    for (int i = 0; i < k; ++i)
        value *= T(n - i);
    return value;
}

template <int N>
std::size_t tensor_coeff_index(const std::array<int, N>& index, int n1)
{
    if constexpr (N == 2)
    {
        return static_cast<std::size_t>(index[0] * n1 + index[1]);
    }
    else
    {
        return static_cast<std::size_t>((index[0] * n1 + index[1]) * n1
                                        + index[2]);
    }
}

template <int N, std::floating_point T>
algoim::real derivative_coefficient(int degree,
                                    std::span<const T> coeffs,
                                    const std::array<int, N>& orders,
                                    const std::array<int, N>& index)
{
    for (int d = 0; d < N; ++d)
        if (orders[static_cast<std::size_t>(d)] > degree)
            return algoim::real(0);

    algoim::real scale = algoim::real(1);
    for (int d = 0; d < N; ++d)
    {
        scale *= falling_factorial<algoim::real>(
            degree, orders[static_cast<std::size_t>(d)]);
    }

    algoim::real value = algoim::real(0);
    for_each_multi_index<N>(orders, [&](const std::array<int, N>& shift) {
        std::array<int, N> coeff_index = index;
        algoim::real sign = algoim::real(1);
        algoim::real weight = algoim::real(1);
        for (int d = 0; d < N; ++d)
        {
            const int od = orders[static_cast<std::size_t>(d)];
            const int sd = shift[static_cast<std::size_t>(d)];
            coeff_index[static_cast<std::size_t>(d)] += sd;
            if ((od - sd) % 2 != 0)
                sign = -sign;
            weight *= binomial<algoim::real>(od, sd);
        }
        value += sign * weight
               * static_cast<algoim::real>(
                   coeffs[tensor_coeff_index<N>(coeff_index, degree + 1)]);
    });

    return scale * value;
}

template <int N, std::floating_point T>
algoim::real evaluate_tensor_bernstein_derivative(
    int degree,
    std::span<const T> coeffs,
    const std::array<int, N>& orders,
    const std::array<algoim::real, N>& x)
{
    std::array<int, N> upper = {};
    for (int d = 0; d < N; ++d)
    {
        if (orders[static_cast<std::size_t>(d)] > degree)
            return algoim::real(0);
        upper[static_cast<std::size_t>(d)] =
            degree - orders[static_cast<std::size_t>(d)];
    }

    std::array<BasisArray<algoim::real>, N> basis;
    for (int d = 0; d < N; ++d)
    {
        bernstein_1d_all(
            upper[static_cast<std::size_t>(d)],
            x[static_cast<std::size_t>(d)],
            basis[static_cast<std::size_t>(d)]);
    }

    algoim::real value = algoim::real(0);
    for_each_multi_index<N>(upper, [&](const std::array<int, N>& index) {
        algoim::real term = derivative_coefficient<N>(
            degree, coeffs, orders, index);
        for (int d = 0; d < N; ++d)
        {
            term *= basis[static_cast<std::size_t>(d)]
                         [static_cast<std::size_t>(
                             index[static_cast<std::size_t>(d)])];
        }
        value += term;
    });
    return value;
}

template <int N, std::floating_point T>
algoim::real derivative_abs_bound(int degree,
                                  std::span<const T> coeffs,
                                  const std::array<int, N>& orders)
{
    std::array<int, N> upper = {};
    for (int d = 0; d < N; ++d)
    {
        if (orders[static_cast<std::size_t>(d)] > degree)
            return algoim::real(0);
        upper[static_cast<std::size_t>(d)] =
            degree - orders[static_cast<std::size_t>(d)];
    }

    algoim::real bound = algoim::real(0);
    for_each_multi_index<N>(upper, [&](const std::array<int, N>& index) {
        using std::abs;
        bound = std::max(
            bound, abs(derivative_coefficient<N>(
                       degree, coeffs, orders, index)));
    });
    return bound;
}

template <int N>
constexpr int derivative_bound_cache_size()
{
    int size = 1;
    for (int i = 0; i < N; ++i)
        size *= derivative_bound_cache_base;
    return size;
}

template <int N>
int derivative_bound_cache_index(const std::array<int, N>& orders)
{
    int index = 0;
    int stride = 1;
    for (int d = 0; d < N; ++d)
    {
        const int order = orders[static_cast<std::size_t>(d)];
        if (order < 0 || order >= derivative_bound_cache_base)
            return -1;
        index += order * stride;
        stride *= derivative_bound_cache_base;
    }
    return index;
}

template <int N, std::floating_point T>
std::array<algoim::real, derivative_bound_cache_size<N>()>
make_derivative_bound_cache(int degree, std::span<const T> coeffs)
{
    std::array<algoim::real, derivative_bound_cache_size<N>()> bounds = {};
    std::array<int, N> upper = {};
    upper.fill(derivative_bound_cache_base - 1);
    for_each_multi_index<N>(upper, [&](const std::array<int, N>& orders) {
        bounds[static_cast<std::size_t>(derivative_bound_cache_index<N>(orders))] =
            derivative_abs_bound<N>(degree, coeffs, orders);
    });
    return bounds;
}

template <int N, std::floating_point T>
algoim::real cached_derivative_abs_bound(
    int degree,
    std::span<const T> coeffs,
    const std::array<int, N>& orders,
    std::span<const algoim::real> derivative_bounds)
{
    const int index = derivative_bound_cache_index<N>(orders);
    if (index >= 0
        && static_cast<std::size_t>(index) < derivative_bounds.size())
    {
        return derivative_bounds[static_cast<std::size_t>(index)];
    }
    return derivative_abs_bound<N>(degree, coeffs, orders);
}

template <int N>
std::array<int, N> derivative_order_with_increment(
    const std::array<int, N>& base,
    int first,
    int second = -1)
{
    std::array<int, N> out = base;
    ++out[static_cast<std::size_t>(first)];
    if (second >= 0)
        ++out[static_cast<std::size_t>(second)];
    return out;
}

template <int N, std::floating_point T>
algoim::real second_order_remainder_bound(
    int degree,
    std::span<const T> coeffs,
    const std::array<int, N>& base_orders,
    const std::array<algoim::real, N>& half_width,
    std::span<const algoim::real> derivative_bounds)
{
    algoim::real eps = algoim::real(0);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            const auto orders =
                derivative_order_with_increment<N>(base_orders, i, j);
            eps += cached_derivative_abs_bound<N>(
                       degree, coeffs, orders, derivative_bounds)
                 * half_width[static_cast<std::size_t>(i)]
                 * half_width[static_cast<std::size_t>(j)];
        }
    }
    return algoim::real(0.5) * eps;
}

template <int N, std::floating_point T>
algoim::Interval<N> evaluate_tensor_bernstein_interval(
    int tdim,
    int degree,
    std::span<const T> coeffs,
    std::span<const algoim::real> derivative_bounds,
    const algoim::uvector<algoim::Interval<N>, N>& x)
{
    if (tdim != N)
        throw std::runtime_error("Algoim Bernstein evaluator dimension mismatch");

    std::array<algoim::real, N> center = {};
    std::array<algoim::real, N> half_width = {};
    for (int d = 0; d < N; ++d)
    {
        center[static_cast<std::size_t>(d)] = x(d).alpha;
        half_width[static_cast<std::size_t>(d)] = x(d).maxDeviation();
    }

    const std::array<int, N> value_orders = {};
    algoim::Interval<N> out(evaluate_tensor_bernstein_derivative<N>(
        degree, coeffs, value_orders, center));

    for (int k = 0; k < N; ++k)
    {
        algoim::real beta = algoim::real(0);
        for (int d = 0; d < N; ++d)
        {
            const auto grad_orders =
                derivative_order_with_increment<N>(value_orders, d);
            beta += evaluate_tensor_bernstein_derivative<N>(
                        degree, coeffs, grad_orders, center)
                  * x(d).beta(k);
        }
        out.beta(k) = beta;
    }

    out.eps = second_order_remainder_bound<N>(
        degree, coeffs, value_orders, half_width, derivative_bounds);
    return out;
}

template <int N, std::floating_point T>
algoim::uvector<algoim::Interval<N>, N> gradient_tensor_bernstein_interval(
    int tdim,
    int degree,
    std::span<const T> coeffs,
    std::span<const algoim::real> derivative_bounds,
    const algoim::uvector<algoim::Interval<N>, N>& x)
{
    if (tdim != N)
        throw std::runtime_error("Algoim Bernstein gradient dimension mismatch");

    std::array<algoim::real, N> center = {};
    std::array<algoim::real, N> half_width = {};
    for (int d = 0; d < N; ++d)
    {
        center[static_cast<std::size_t>(d)] = x(d).alpha;
        half_width[static_cast<std::size_t>(d)] = x(d).maxDeviation();
    }

    algoim::uvector<algoim::Interval<N>, N> out;
    const std::array<int, N> value_orders = {};
    for (int component = 0; component < N; ++component)
    {
        const auto grad_orders =
            derivative_order_with_increment<N>(value_orders, component);
        out(component) = algoim::Interval<N>(
            evaluate_tensor_bernstein_derivative<N>(
                degree, coeffs, grad_orders, center));

        for (int k = 0; k < N; ++k)
        {
            algoim::real beta = algoim::real(0);
            for (int d = 0; d < N; ++d)
            {
                const auto hess_orders =
                    derivative_order_with_increment<N>(grad_orders, d);
                beta += evaluate_tensor_bernstein_derivative<N>(
                            degree, coeffs, hess_orders, center)
                      * x(d).beta(k);
            }
            out(component).beta(k) = beta;
        }

        out(component).eps = second_order_remainder_bound<N>(
            degree, coeffs, grad_orders, half_width, derivative_bounds);
    }
    return out;
}

template <typename S, std::floating_point T>
S evaluate_tensor_bernstein(int tdim,
                            int degree,
                            std::span<const T> coeffs,
                            std::span<const algoim::real> derivative_bounds,
                            const algoim::uvector<S, 2>& x)
{
    if constexpr (std::is_same_v<S, algoim::Interval<2>>)
    {
        return evaluate_tensor_bernstein_interval<2>(
            tdim, degree, coeffs, derivative_bounds, x);
    }
    else
    {
    const int n = degree;
    const int n1 = n + 1;

    BasisArray<S> bx;
    BasisArray<S> by;
    bernstein_1d_all(n, x(0), bx);
    bernstein_1d_all(n, x(1), by);

    S result = S(0);
    if (tdim == 2)
    {
        for (int i = 0; i <= n; ++i)
        {
            for (int j = 0; j <= n; ++j)
            {
                result += S(coeffs[static_cast<std::size_t>(i * n1 + j)])
                        * bx[static_cast<std::size_t>(i)]
                        * by[static_cast<std::size_t>(j)];
            }
        }
        return result;
    }

    throw std::runtime_error("Algoim Bernstein evaluator: expected tdim=2");
    }
}

template <typename S, std::floating_point T>
S evaluate_tensor_bernstein(int tdim,
                            int degree,
                            std::span<const T> coeffs,
                            std::span<const algoim::real> derivative_bounds,
                            const algoim::uvector<S, 3>& x)
{
    if constexpr (std::is_same_v<S, algoim::Interval<3>>)
    {
        return evaluate_tensor_bernstein_interval<3>(
            tdim, degree, coeffs, derivative_bounds, x);
    }
    else
    {
    const int n = degree;
    const int n1 = n + 1;

    BasisArray<S> bx;
    BasisArray<S> by;
    BasisArray<S> bz;
    bernstein_1d_all(n, x(0), bx);
    bernstein_1d_all(n, x(1), by);
    bernstein_1d_all(n, x(2), bz);

    S result = S(0);
    if (tdim == 3)
    {
        for (int i = 0; i <= n; ++i)
        {
            for (int j = 0; j <= n; ++j)
            {
                for (int k = 0; k <= n; ++k)
                {
                    result += S(coeffs[static_cast<std::size_t>(
                                  (i * n1 + j) * n1 + k)])
                            * bx[static_cast<std::size_t>(i)]
                            * by[static_cast<std::size_t>(j)]
                            * bz[static_cast<std::size_t>(k)];
                }
            }
        }
        return result;
    }

    throw std::runtime_error("Algoim Bernstein evaluator: expected tdim=3");
    }
}

template <typename S, std::floating_point T>
algoim::uvector<S, 2> gradient_tensor_bernstein(int tdim,
                                                int degree,
                                                std::span<const T> coeffs,
                                                std::span<const algoim::real> derivative_bounds,
                                                const algoim::uvector<S, 2>& x)
{
    if (tdim != 2)
        throw std::runtime_error("Algoim Bernstein gradient: expected tdim=2");

    if constexpr (std::is_same_v<S, algoim::Interval<2>>)
    {
        return gradient_tensor_bernstein_interval<2>(
            tdim, degree, coeffs, derivative_bounds, x);
    }
    else
    {
    const int n = degree;
    const int n1 = n + 1;

    BasisArray<S> bx;
    BasisArray<S> by;
    BasisArray<S> dbx;
    BasisArray<S> dby;
    bernstein_1d_all(n, x(0), bx);
    bernstein_1d_all(n, x(1), by);
    bernstein_1d_deriv(n, x(0), dbx);
    bernstein_1d_deriv(n, x(1), dby);

    S g0 = S(0);
    S g1 = S(0);
    for (int i = 0; i <= n; ++i)
    {
        for (int j = 0; j <= n; ++j)
        {
            const S c = S(coeffs[static_cast<std::size_t>(i * n1 + j)]);
            g0 += c * dbx[static_cast<std::size_t>(i)]
                    * by[static_cast<std::size_t>(j)];
            g1 += c * bx[static_cast<std::size_t>(i)]
                    * dby[static_cast<std::size_t>(j)];
        }
    }

    return algoim::uvector<S, 2>(g0, g1);
    }
}

template <typename S, std::floating_point T>
algoim::uvector<S, 3> gradient_tensor_bernstein(int tdim,
                                                int degree,
                                                std::span<const T> coeffs,
                                                std::span<const algoim::real> derivative_bounds,
                                                const algoim::uvector<S, 3>& x)
{
    if (tdim != 3)
        throw std::runtime_error("Algoim Bernstein gradient: expected tdim=3");

    if constexpr (std::is_same_v<S, algoim::Interval<3>>)
    {
        return gradient_tensor_bernstein_interval<3>(
            tdim, degree, coeffs, derivative_bounds, x);
    }
    else
    {
    const int n = degree;
    const int n1 = n + 1;

    BasisArray<S> bx;
    BasisArray<S> by;
    BasisArray<S> bz;
    BasisArray<S> dbx;
    BasisArray<S> dby;
    BasisArray<S> dbz;
    bernstein_1d_all(n, x(0), bx);
    bernstein_1d_all(n, x(1), by);
    bernstein_1d_all(n, x(2), bz);
    bernstein_1d_deriv(n, x(0), dbx);
    bernstein_1d_deriv(n, x(1), dby);
    bernstein_1d_deriv(n, x(2), dbz);

    S g0 = S(0);
    S g1 = S(0);
    S g2 = S(0);
    for (int i = 0; i <= n; ++i)
    {
        for (int j = 0; j <= n; ++j)
        {
            for (int k = 0; k <= n; ++k)
            {
                const S c = S(coeffs[static_cast<std::size_t>(
                              (i * n1 + j) * n1 + k)]);
                const S bxi = bx[static_cast<std::size_t>(i)];
                const S byj = by[static_cast<std::size_t>(j)];
                const S bzk = bz[static_cast<std::size_t>(k)];
                g0 += c * dbx[static_cast<std::size_t>(i)] * byj * bzk;
                g1 += c * bxi * dby[static_cast<std::size_t>(j)] * bzk;
                g2 += c * bxi * byj * dbz[static_cast<std::size_t>(k)];
            }
        }
    }

    return algoim::uvector<S, 3>(g0, g1, g2);
    }
}

template <std::floating_point T, int N>
struct BernsteinLevelSet
{
    cell::type cell_type = cell::type::point;
    int tdim = 0;
    int degree = 0;
    std::span<const T> coeffs;
    std::span<const algoim::real> derivative_bounds;
    T sign = T(1);

    template <typename S>
    S operator()(const algoim::uvector<S, N>& x) const
    {
        return S(sign)
             * evaluate_tensor_bernstein(
                   tdim, degree, coeffs, derivative_bounds, x);
    }

    template <typename S>
    algoim::uvector<S, N> grad(const algoim::uvector<S, N>& x) const
    {
        auto g = gradient_tensor_bernstein(
            tdim, degree, coeffs, derivative_bounds, x);
        g *= S(sign);
        return g;
    }
};

template <std::floating_point T, std::integral I>
std::vector<T> parent_cell_vertex_coords_basix(const MeshView<T, I>& mesh,
                                               I cell_id)
{
    const auto ctype = mesh.cell_type(cell_id);
    const int nv = cell::get_num_vertices(ctype);
    std::vector<T> coords(static_cast<std::size_t>(nv * mesh.gdim), T(0));
    std::vector<I> cell_node_scratch;
    const auto parent_nodes = mesh.cell_nodes(cell_id, cell_node_scratch);

    for (int basix_v = 0; basix_v < nv; ++basix_v)
    {
        const int local_v = mesh.vtk_vertex_order
                                ? cell::basix_to_vtk_vertex(ctype, basix_v)
                                : basix_v;
        const I node_id = parent_nodes[static_cast<std::size_t>(local_v)];
        const T* x = mesh.node(node_id);
        for (int d = 0; d < mesh.gdim; ++d)
            coords[static_cast<std::size_t>(basix_v * mesh.gdim + d)] = x[d];
    }
    return coords;
}

inline bool supported_tensor_cell(cell::type type)
{
    return type == cell::type::quadrilateral
        || type == cell::type::hexahedron;
}

inline int single_clause_level_set_index(const SelectionExpr& expr,
                                         Relation& relation)
{
    if (expr.terms.size() != 1 || expr.terms.front().clauses.size() != 1)
    {
        throw std::runtime_error(
            "Algoim quadrature currently supports a single level-set clause");
    }

    const auto& clause = expr.terms.front().clauses.front();
    relation = clause.relation;
    if (clause.level_set_index < 0)
    {
        throw std::runtime_error(
            "Algoim quadrature requires a compiled level-set selection");
    }
    return clause.level_set_index;
}

template <std::floating_point T>
T sign_for_relation(Relation relation)
{
    switch (relation)
    {
    case Relation::LessThan:
    case Relation::LessEqual:
    case Relation::EqualTo:
        return T(1);
    case Relation::GreaterThan:
    case Relation::GreaterEqual:
        return T(-1);
    }
    throw std::runtime_error("Algoim quadrature: unsupported relation");
}

inline bool relation_selects_interface(Relation relation)
{
    return relation == Relation::EqualTo;
}

template <std::floating_point T>
T interval_root_xtol()
{
    if constexpr (std::is_same_v<T, float>)
        return T(1e-6);
    else
        return T(1e-12);
}

inline bool relation_selects_negative_volume(Relation relation)
{
    return relation == Relation::LessThan || relation == Relation::LessEqual;
}

inline bool relation_selects_positive_volume(Relation relation)
{
    return relation == Relation::GreaterThan || relation == Relation::GreaterEqual;
}

template <std::floating_point T>
bool relation_selects_value(Relation relation, T value)
{
    if (relation_selects_negative_volume(relation))
        return value < T(0);
    if (relation_selects_positive_volume(relation))
        return value > T(0);
    throw std::runtime_error("Algoim volume quadrature requires a sign relation");
}

template <std::floating_point T, std::integral I>
const LevelSetCell<T, I>& active_level_set_cell(const HOCutCells<T, I>& cut_cells,
                                                int cut_id,
                                                int level_set_index)
{
    const int begin = cut_cells.ls_offsets.at(static_cast<std::size_t>(cut_id));
    const int end = cut_cells.ls_offsets.at(static_cast<std::size_t>(cut_id + 1));
    for (int i = begin; i < end; ++i)
    {
        const auto& ls_cell =
            cut_cells.level_set_cells[static_cast<std::size_t>(i)];
        if (ls_cell.level_set_id == level_set_index)
            return ls_cell;
    }

    throw std::runtime_error(
        "Algoim quadrature could not find the selected active level set");
}

template <std::floating_point T>
void initialise_rules(quadrature::QuadratureRules<T>& rules, int tdim)
{
    if (rules._tdim == 0)
        rules._tdim = tdim;
    if (rules._offset.empty())
        rules._offset.push_back(0);
}

template <std::floating_point T, std::integral I>
T evaluate_interval_level_set(const LevelSetCell<T, I>& ls_cell, T t)
{
    t = std::clamp(t, T(0), T(1));
    std::array<T, 1> xi = {t};
    return ls_cell.value(
        std::span<const T>(xi.data(), xi.size()));
}

template <std::floating_point T, std::integral I>
bool solve_interval_root_between(const LevelSetCell<T, I>& ls_cell,
                                 T a,
                                 T b,
                                 T& root)
{
    constexpr int max_iter = 64;
    const T xtol = interval_root_xtol<T>();
    const T ftol = cell::edge_root::default_value_tolerance<T>();

    a = std::clamp(a, T(0), T(1));
    b = std::clamp(b, T(0), T(1));
    if (b < a)
        std::swap(a, b);

    const T fa = evaluate_interval_level_set(ls_cell, a);
    if (std::abs(fa) <= ftol)
    {
        root = a;
        return true;
    }
    const T fb = evaluate_interval_level_set(ls_cell, b);
    if (std::abs(fb) <= ftol)
    {
        root = b;
        return true;
    }
    if (fa * fb > T(0))
        return false;

    auto phi = [&](std::span<const T> xi) -> T
    {
        return ls_cell.value(xi);
    };
    std::array<T, 1> p0 = {a};
    std::array<T, 1> p1 = {b};
    const auto info = cell::edge_root::find_root_parameter_info<T>(
        std::span<const T>(p0.data(), p0.size()),
        std::span<const T>(p1.data(), p1.size()),
        phi,
        cell::edge_root::method::brent,
        T(0),
        max_iter,
        xtol,
        ftol);
    root = a + info.t * (b - a);
    return true;
}

template <std::floating_point T, std::integral I>
T polish_interval_interface_root(const LevelSetCell<T, I>& ls_cell,
                                 T guess)
{
    const T xtol = interval_root_xtol<T>();
    const T ftol = cell::edge_root::default_value_tolerance<T>();
    guess = std::clamp(guess, T(0), T(1));

    if (std::abs(evaluate_interval_level_set(ls_cell, guess)) <= ftol)
        return guess;

    std::vector<T> roots;
    auto add_root = [&](T root)
    {
        root = std::clamp(root, T(0), T(1));
        for (T existing : roots)
        {
            if (std::abs(existing - root) <= T(16) * xtol)
                return;
        }
        roots.push_back(root);
    };

    const int degree = std::max(1, ls_cell.bernstein_order);
    const int intervals = std::max(16, 4 * (degree + 1));
    std::vector<T> grid;
    grid.reserve(static_cast<std::size_t>(intervals + 2));
    for (int i = 0; i <= intervals; ++i)
        grid.push_back(T(i) / T(intervals));
    grid.push_back(guess);
    std::ranges::sort(grid);
    grid.erase(std::unique(
                   grid.begin(), grid.end(),
                   [xtol](T a, T b)
                   { return std::abs(a - b) <= T(16) * xtol; }),
               grid.end());

    for (std::size_t i = 1; i < grid.size(); ++i)
    {
        T root = guess;
        if (solve_interval_root_between(
                ls_cell, grid[i - 1], grid[i], root))
        {
            add_root(root);
        }
    }

    if (!roots.empty())
    {
        return *std::min_element(
            roots.begin(), roots.end(),
            [guess](T a, T b)
            {
                return std::abs(a - guess) < std::abs(b - guess);
            });
    }

    T root = guess;
    if (solve_interval_root_between(ls_cell, T(0), T(1), root))
        return root;

    return guess;
}

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> interval_interface_quadrature_rules(
    const HOMeshPart<T, I>& part,
    int level_set_index)
{
    if (!part.mesh || !part.cut_cells || !part.parent_cells)
        throw std::runtime_error("HOMeshPart is not attached to cut-cell storage");
    if (part.mesh->tdim != 1 || part.dim != 0)
    {
        throw std::runtime_error(
            "Interval interface quadrature expects a zero-dimensional part "
            "on an interval host mesh");
    }

    quadrature::QuadratureRules<T> rules;
    initialise_rules(rules, part.mesh->tdim);

    const auto infos = selected_zero_entity_infos(part);
    for (const auto& info : infos)
    {
        if (info.dimension != 0)
        {
            throw std::runtime_error(
                "Interval interface quadrature found a non-point zero entity");
        }

        const auto& adapt_cell =
            part.cut_cells->adapt_cells[static_cast<std::size_t>(
                info.cut_cell_id)];
        const int zero_id = static_cast<int>(info.local_zero_entity_id);
        if (zero_id < 0
            || static_cast<std::size_t>(zero_id)
                   >= adapt_cell.zero_entity_dim.size()
            || adapt_cell.zero_entity_dim[static_cast<std::size_t>(zero_id)] != 0)
        {
            throw std::runtime_error(
                "Interval interface quadrature found invalid zero-entity data");
        }

        const int vertex_id =
            adapt_cell.zero_entity_id[static_cast<std::size_t>(zero_id)];
        if (vertex_id < 0
            || static_cast<std::size_t>(vertex_id * adapt_cell.tdim)
                   >= adapt_cell.vertex_coords.size())
        {
            throw std::runtime_error(
                "Interval interface quadrature found invalid root vertex data");
        }

        const auto& ls_cell = active_level_set_cell(
            *part.cut_cells, info.cut_cell_id, level_set_index);
        if (ls_cell.cell_type != cell::type::interval || ls_cell.tdim != 1)
        {
            throw std::runtime_error(
                "Interval interface quadrature requires interval level-set cells");
        }

        const T guess =
            adapt_cell.vertex_coords[static_cast<std::size_t>(
                vertex_id * adapt_cell.tdim)];
        const T root = polish_interval_interface_root(ls_cell, guess);
        rules._points.push_back(root);
        rules._weights.push_back(T(1));
        rules._parent_map.push_back(info.parent_cell_id);
        rules._offset.push_back(
            static_cast<std::int32_t>(rules._weights.size()));
    }

    return rules;
}

template <std::floating_point T>
void merge_rules(quadrature::QuadratureRules<T>& rules,
                 const quadrature::QuadratureRules<T>& extra)
{
    if (extra._weights.empty())
        return;

    initialise_rules(rules, extra._tdim);
    if (rules._tdim != extra._tdim)
        throw std::runtime_error("quadrature rule dimension mismatch");

    const auto offset_shift = static_cast<std::int32_t>(rules._weights.size());
    rules._points.insert(rules._points.end(),
                         extra._points.begin(),
                         extra._points.end());
    rules._weights.insert(rules._weights.end(),
                          extra._weights.begin(),
                          extra._weights.end());
    rules._parent_map.insert(rules._parent_map.end(),
                             extra._parent_map.begin(),
                             extra._parent_map.end());
    for (std::size_t i = 1; i < extra._offset.size(); ++i)
        rules._offset.push_back(extra._offset[i] + offset_shift);
}

template <std::floating_point T>
struct AffineCellGeometry
{
    T volume_scale = T(0);
    std::array<T, 9> metric_inverse = {};
    int tdim = 0;
    int gdim = 0;
};

template <std::floating_point T>
AffineCellGeometry<T> affine_cell_geometry(cell::type cell_type,
                                           std::span<const T> parent_vertices,
                                           int tdim,
                                           int gdim)
{
    AffineCellGeometry<T> geometry;
    geometry.tdim = tdim;
    geometry.gdim = gdim;
    if (tdim < 1 || tdim > 3 || gdim < tdim || gdim > 3)
    {
        throw std::runtime_error(
            "Algoim affine geometry requires 1 <= tdim <= gdim <= 3");
    }

    std::array<T, 9> J = {};
    const auto cols = cell::jacobian_col_indices(cell_type);
    for (int c = 0; c < tdim; ++c)
    {
        const int vertex_id = cols[static_cast<std::size_t>(c)];
        for (int r = 0; r < gdim; ++r)
        {
            J[static_cast<std::size_t>(c * gdim + r)] =
                parent_vertices[static_cast<std::size_t>(vertex_id * gdim + r)]
              - parent_vertices[static_cast<std::size_t>(r)];
        }
    }

    std::array<T, 9> metric = {};
    for (int i = 0; i < tdim; ++i)
    {
        for (int j = 0; j < tdim; ++j)
        {
            T value = T(0);
            for (int r = 0; r < gdim; ++r)
                value += J[static_cast<std::size_t>(i * gdim + r)]
                       * J[static_cast<std::size_t>(j * gdim + r)];
            metric[static_cast<std::size_t>(i * tdim + j)] = value;
        }
    }

    if (tdim == 1)
    {
        const T det = metric[0];
        if (det <= T(0))
            throw std::runtime_error("Algoim quadrature found a degenerate cell");
        geometry.volume_scale = std::sqrt(det);
        geometry.metric_inverse[0] = T(1) / det;
        return geometry;
    }

    if (tdim == 2)
    {
        const T a = metric[0];
        const T b = metric[1];
        const T c = metric[2];
        const T d = metric[3];
        const T det = a * d - b * c;
        if (det <= T(0))
            throw std::runtime_error("Algoim quadrature found a degenerate cell");
        geometry.volume_scale = std::sqrt(det);
        geometry.metric_inverse[0] = d / det;
        geometry.metric_inverse[1] = -b / det;
        geometry.metric_inverse[2] = -c / det;
        geometry.metric_inverse[3] = a / det;
        return geometry;
    }

    if (tdim == 3)
    {
        const T a00 = metric[0];
        const T a01 = metric[1];
        const T a02 = metric[2];
        const T a10 = metric[3];
        const T a11 = metric[4];
        const T a12 = metric[5];
        const T a20 = metric[6];
        const T a21 = metric[7];
        const T a22 = metric[8];
        const T det = a00 * (a11 * a22 - a12 * a21)
                    - a01 * (a10 * a22 - a12 * a20)
                    + a02 * (a10 * a21 - a11 * a20);
        if (det <= T(0))
            throw std::runtime_error("Algoim quadrature found a degenerate cell");
        geometry.volume_scale = std::sqrt(det);
        geometry.metric_inverse[0] = (a11 * a22 - a12 * a21) / det;
        geometry.metric_inverse[1] = (a02 * a21 - a01 * a22) / det;
        geometry.metric_inverse[2] = (a01 * a12 - a02 * a11) / det;
        geometry.metric_inverse[3] = (a12 * a20 - a10 * a22) / det;
        geometry.metric_inverse[4] = (a00 * a22 - a02 * a20) / det;
        geometry.metric_inverse[5] = (a02 * a10 - a00 * a12) / det;
        geometry.metric_inverse[6] = (a10 * a21 - a11 * a20) / det;
        geometry.metric_inverse[7] = (a01 * a20 - a00 * a21) / det;
        geometry.metric_inverse[8] = (a00 * a11 - a01 * a10) / det;
        return geometry;
    }

    throw std::runtime_error(
        "Algoim affine geometry supports only dimensions 1, 2, and 3");
}

template <std::floating_point T>
T norm2(const std::array<T, 3>& v, int n)
{
    T sum = T(0);
    for (int i = 0; i < n; ++i)
        sum += v[static_cast<std::size_t>(i)]
             * v[static_cast<std::size_t>(i)];
    return std::sqrt(sum);
}

template <std::floating_point T>
T interface_scale(const AffineCellGeometry<T>& geometry,
                  std::span<const T> grad_ref)
{
    std::array<T, 3> grad = {};
    for (int i = 0; i < geometry.tdim; ++i)
        grad[static_cast<std::size_t>(i)] = grad_ref[static_cast<std::size_t>(i)];

    const T grad_norm = norm2(grad, geometry.tdim);
    if (grad_norm <= T(0))
        throw std::runtime_error("Algoim interface quadrature found zero gradient");

    T metric_grad_norm_sq = T(0);
    for (int i = 0; i < geometry.tdim; ++i)
    {
        for (int j = 0; j < geometry.tdim; ++j)
        {
            metric_grad_norm_sq += grad[static_cast<std::size_t>(i)]
                * geometry.metric_inverse[static_cast<std::size_t>(
                    i * geometry.tdim + j)]
                * grad[static_cast<std::size_t>(j)];
        }
    }
    if (metric_grad_norm_sq <= T(0))
        throw std::runtime_error("Algoim interface quadrature found zero metric gradient");
    return geometry.volume_scale * std::sqrt(metric_grad_norm_sq) / grad_norm;
}

template <std::floating_point T, int N>
std::array<T, 3> reference_gradient(const BernsteinLevelSet<T, N>& phi,
                                    const algoim::uvector<algoim::real, N>& x)
{
    algoim::uvector<T, N> xt;
    for (int i = 0; i < N; ++i)
        xt(i) = static_cast<T>(x(i));

    const auto g = phi.template grad<T>(xt);
    std::array<T, 3> out = {};
    for (int i = 0; i < N; ++i)
        out[static_cast<std::size_t>(i)] = g(i);
    return out;
}

template <std::floating_point T, std::integral I, int N>
void append_algoim_general_cell(quadrature::QuadratureRules<T>& rules,
                                const LevelSetCell<T, I>& ls_cell,
                                std::span<const T> parent_vertices,
                                I parent_cell_id,
                                Relation relation,
                                int order)
{
    const bool interface = relation_selects_interface(relation);
    BernsteinLevelSet<T, N> phi;
    phi.cell_type = ls_cell.cell_type;
    phi.tdim = ls_cell.tdim;
    phi.degree = ls_cell.bernstein_order;
    phi.coeffs = std::span<const T>(
        ls_cell.bernstein_coeffs.data(), ls_cell.bernstein_coeffs.size());
    const auto derivative_bounds = make_derivative_bound_cache<N>(
        phi.degree, phi.coeffs);
    phi.derivative_bounds = std::span<const algoim::real>(
        derivative_bounds.data(), derivative_bounds.size());
    phi.sign = sign_for_relation<T>(relation);

    auto range = algoim::HyperRectangle<algoim::real, N>(
        algoim::uvector<algoim::real, N>(algoim::real(0)),
        algoim::uvector<algoim::real, N>(algoim::real(1)));
    const auto q = algoim::quadGen<N>(
        phi, range, interface ? N : -1, -1, order);

    initialise_rules(rules, ls_cell.tdim);
    const auto geometry = affine_cell_geometry<T>(
        ls_cell.cell_type, parent_vertices, ls_cell.tdim, ls_cell.gdim);
    const auto start = rules._weights.size();
    rules._points.reserve(rules._points.size() + q.nodes.size() * N);
    rules._weights.reserve(rules._weights.size() + q.nodes.size());

    const T det_scale = interface ? T(0) : geometry.volume_scale;

    for (const auto& node : q.nodes)
    {
        for (int d = 0; d < N; ++d)
            rules._points.push_back(static_cast<T>(node.x(d)));

        T weight_scale = det_scale;
        if (interface)
        {
            const auto grad = reference_gradient(phi, node.x);
            weight_scale = interface_scale<T>(
                geometry,
                std::span<const T>(grad.data(), static_cast<std::size_t>(N)));
        }

        rules._weights.push_back(static_cast<T>(node.w) * weight_scale);
    }

    if (rules._weights.size() != start)
    {
        rules._parent_map.push_back(static_cast<std::int32_t>(parent_cell_id));
        rules._offset.push_back(static_cast<std::int32_t>(rules._weights.size()));
    }
}

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> algoim_general_quadrature_rules_impl(
    const HOMeshPart<T, I>& part,
    int order,
    bool include_uncut_cells)
{
    if (!part.mesh)
        throw std::runtime_error("HOMeshPart is not attached to a mesh");
    if (!part.cut_cells || !part.parent_cells)
        throw std::runtime_error("Algoim quadrature requires cut-cell storage");
    if (order < 1 || order > 10)
        throw std::runtime_error("Algoim quadrature order must be in [1, 10]");

    Relation relation = Relation::LessThan;
    const int level_set_index =
        single_clause_level_set_index(part.expr, relation);
    const bool interface = relation_selects_interface(relation);

    const int mesh_tdim = part.mesh->tdim;
    if (interface)
    {
        if (part.dim != mesh_tdim - 1)
        {
            throw std::runtime_error(
                "Algoim interface quadrature expects a codimension-one mesh part");
        }
    }
    else if (part.dim != mesh_tdim)
    {
        throw std::runtime_error(
            "Algoim volume quadrature expects a volume mesh part");
    }

    if (mesh_tdim == 1 && interface)
        return interval_interface_quadrature_rules(part, level_set_index);

    quadrature::QuadratureRules<T> rules;
    initialise_rules(rules, mesh_tdim);

    for (const std::int32_t cut_id : part.cut_cell_ids)
    {
        const I parent_cell_id =
            part.cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_id)];
        const auto& ls_cell = active_level_set_cell(
            *part.cut_cells, cut_id, level_set_index);

        if (!supported_tensor_cell(ls_cell.cell_type))
        {
            throw std::runtime_error(
                "Algoim quadrature currently supports quadrilateral and "
                "hexahedron parent cells only");
        }
        if (ls_cell.bernstein_coeffs.empty())
        {
            throw std::runtime_error(
                "Algoim quadrature currently requires Bernstein level-set data");
        }
        if (ls_cell.tdim != mesh_tdim)
        {
            throw std::runtime_error(
                "Algoim quadrature found inconsistent mesh dimensions");
        }
        std::vector<T> parent_vertex_storage;
        std::span<const T> parent_vertices;
        const int nv = cell::get_num_vertices(ls_cell.cell_type);
        if (ls_cell.parent_vertex_coords.size()
            == static_cast<std::size_t>(nv * ls_cell.gdim))
        {
            parent_vertices = std::span<const T>(
                ls_cell.parent_vertex_coords.data(),
                ls_cell.parent_vertex_coords.size());
        }
        else
        {
            parent_vertex_storage =
                parent_cell_vertex_coords_basix(*part.mesh, parent_cell_id);
            parent_vertices = std::span<const T>(
                parent_vertex_storage.data(), parent_vertex_storage.size());
        }

        if (mesh_tdim == 2)
        {
            append_algoim_general_cell<T, I, 2>(
                rules, ls_cell, parent_vertices,
                parent_cell_id, relation, order);
        }
        else if (mesh_tdim == 3)
        {
            append_algoim_general_cell<T, I, 3>(
                rules, ls_cell, parent_vertices,
                parent_cell_id, relation, order);
        }
        else
        {
            throw std::runtime_error(
                "Algoim quadrature currently supports 2D and 3D only");
        }
    }

    if (include_uncut_cells && !interface && !part.uncut_cell_ids.empty())
    {
        HOMeshPart<T, I> uncut_part = part;
        uncut_part.cut_cell_ids.clear();
        const auto uncut_rules = quadrature_rules(
            uncut_part, order, /*include_uncut_cells=*/true);
        merge_rules(rules, uncut_rules);
    }

    return rules;
}

#ifdef CUTCELLS_HAS_ALGOIM_MULTIPOLY

template <std::floating_point T, std::integral I>
struct PairedRequest
{
    std::string name;
    const HOMeshPart<T, I>* part = nullptr;
    Relation relation = Relation::LessThan;
    int level_set_index = -1;
    bool interface = false;
};

template <int N>
algoim::uvector<int, N> tensor_extent(int n1)
{
    if constexpr (N == 2)
        return algoim::uvector<int, N>(n1, n1);
    else
        return algoim::uvector<int, N>(n1, n1, n1);
}

template <int N, std::floating_point T, std::integral I>
std::vector<algoim::real> algoim_coefficients(const LevelSetCell<T, I>& ls_cell)
{
    const int n1 = ls_cell.bernstein_order + 1;
    std::size_t expected = 1;
    for (int d = 0; d < N; ++d)
        expected *= static_cast<std::size_t>(n1);
    if (ls_cell.bernstein_coeffs.size() != expected)
    {
        throw std::runtime_error(
            "Algoim Bernstein quadrature found inconsistent coefficient storage");
    }

    std::vector<algoim::real> coeffs(expected);
    for (std::size_t i = 0; i < expected; ++i)
        coeffs[i] = static_cast<algoim::real>(ls_cell.bernstein_coeffs[i]);
    return coeffs;
}

template <std::floating_point T, int N>
void append_multipoly_point(quadrature::QuadratureRules<T>& rules,
                            const algoim::uvector<algoim::real, N>& x,
                            T weight)
{
    for (int d = 0; d < N; ++d)
        rules._points.push_back(static_cast<T>(x(d)));
    rules._weights.push_back(weight);
}

template <std::floating_point T, int N>
std::array<T, 3> multipoly_reference_gradient(
    const algoim::xarray<algoim::real, N>& poly,
    const algoim::uvector<algoim::real, N>& x)
{
    const auto g = algoim::bernstein::evalBernsteinPolyGradient(poly, x);
    std::array<T, 3> out = {};
    for (int i = 0; i < N; ++i)
        out[static_cast<std::size_t>(i)] = static_cast<T>(g(i));
    return out;
}

template <std::floating_point T, std::integral I, int N>
void append_algoim_bernstein_cell(
    std::vector<std::pair<std::string, quadrature::QuadratureRules<T>>>& out,
    const std::vector<PairedRequest<T, I>>& requests,
    const LevelSetCell<T, I>& ls_cell,
    std::span<const T> parent_vertices,
    I parent_cell_id,
    int order)
{
    auto coeffs = algoim_coefficients<N>(ls_cell);
    algoim::xarray<algoim::real, N> poly(
        coeffs.data(), tensor_extent<N>(ls_cell.bernstein_order + 1));
    algoim::ImplicitPolyQuadrature<N> ipquad(poly);

    const auto geometry = affine_cell_geometry<T>(
        ls_cell.cell_type, parent_vertices, ls_cell.tdim, ls_cell.gdim);

    std::vector<std::size_t> starts(out.size());
    for (std::size_t i = 0; i < out.size(); ++i)
        starts[i] = out[i].second._weights.size();

    bool has_volume_request = false;
    bool has_interface_request = false;
    for (const auto& request : requests)
    {
        has_volume_request = has_volume_request || !request.interface;
        has_interface_request = has_interface_request || request.interface;
    }

    if (has_volume_request)
    {
        ipquad.integrate(algoim::AutoMixed, order,
                         [&](const algoim::uvector<algoim::real, N>& x,
                             algoim::real w)
        {
            const T value = static_cast<T>(
                algoim::bernstein::evalBernsteinPoly(poly, x));
            const T scaled_weight = static_cast<T>(w) * geometry.volume_scale;
            for (std::size_t i = 0; i < requests.size(); ++i)
            {
                const auto& request = requests[i];
                if (request.interface)
                    continue;
                if (!relation_selects_value(request.relation, value))
                    continue;
                append_multipoly_point(out[i].second, x, scaled_weight);
            }
        });
    }

    if (has_interface_request)
    {
        ipquad.integrate_surf(
            algoim::AutoMixed, order,
            [&](const algoim::uvector<algoim::real, N>& x,
                algoim::real w,
                const algoim::uvector<algoim::real, N>&)
        {
            const auto grad = multipoly_reference_gradient<T>(poly, x);
            const T scale = interface_scale<T>(
                geometry,
                std::span<const T>(grad.data(), static_cast<std::size_t>(N)));
            const T scaled_weight = static_cast<T>(w) * scale;
            for (std::size_t i = 0; i < requests.size(); ++i)
            {
                if (!requests[i].interface)
                    continue;
                append_multipoly_point(out[i].second, x, scaled_weight);
            }
        });
    }

    for (std::size_t i = 0; i < out.size(); ++i)
    {
        auto& rules = out[i].second;
        if (rules._weights.size() == starts[i])
            continue;

        rules._parent_map.push_back(static_cast<std::int32_t>(parent_cell_id));
        rules._offset.push_back(static_cast<std::int32_t>(rules._weights.size()));
    }
}

template <std::floating_point T, std::integral I>
std::vector<PairedRequest<T, I>> prepare_paired_requests(
    const std::vector<std::pair<std::string, HOMeshPart<T, I>>>& parts,
    int mesh_tdim)
{
    std::vector<PairedRequest<T, I>> requests;
    requests.reserve(parts.size());
    int level_set_index = -1;

    for (const auto& [name, part] : parts)
    {
        if (!part.mesh)
            throw std::runtime_error("HOMeshPart is not attached to a mesh");
        if (!part.cut_cells || !part.parent_cells)
            throw std::runtime_error("Algoim quadrature requires cut-cell storage");

        Relation relation = Relation::LessThan;
        const int part_level_set =
            single_clause_level_set_index(part.expr, relation);
        if (level_set_index < 0)
            level_set_index = part_level_set;
        else if (part_level_set != level_set_index)
        {
            throw std::runtime_error(
                "Algoim paired quadrature currently requires all selectors to "
                "refer to the same level set");
        }

        const bool interface = relation_selects_interface(relation);
        if (interface)
        {
            if (part.dim != mesh_tdim - 1)
            {
                throw std::runtime_error(
                    "Algoim interface quadrature expects a codimension-one mesh part");
            }
        }
        else if (part.dim != mesh_tdim)
        {
            throw std::runtime_error(
                "Algoim volume quadrature expects a volume mesh part");
        }

        requests.push_back(
            PairedRequest<T, I>{name, &part, relation, part_level_set, interface});
    }

    return requests;
}

template <std::floating_point T, std::integral I>
std::vector<std::int32_t> paired_cut_cell_ids(
    const std::vector<PairedRequest<T, I>>& requests)
{
    std::vector<std::int32_t> cut_ids;
    for (const auto& request : requests)
    {
        cut_ids.insert(cut_ids.end(),
                       request.part->cut_cell_ids.begin(),
                       request.part->cut_cell_ids.end());
    }
    std::ranges::sort(cut_ids);
    cut_ids.erase(std::unique(cut_ids.begin(), cut_ids.end()), cut_ids.end());
    return cut_ids;
}

template <std::floating_point T, std::integral I>
std::span<const T> parent_vertices_for_algoim_cell(
    const MeshView<T, I>& mesh,
    const LevelSetCell<T, I>& ls_cell,
    I parent_cell_id,
    std::vector<T>& storage)
{
    const int nv = cell::get_num_vertices(ls_cell.cell_type);
    if (ls_cell.parent_vertex_coords.size()
        == static_cast<std::size_t>(nv * ls_cell.gdim))
    {
        return std::span<const T>(
            ls_cell.parent_vertex_coords.data(),
            ls_cell.parent_vertex_coords.size());
    }

    storage = parent_cell_vertex_coords_basix(mesh, parent_cell_id);
    return std::span<const T>(storage.data(), storage.size());
}

template <std::floating_point T, std::integral I>
std::vector<std::pair<std::string, quadrature::QuadratureRules<T>>>
algoim_paired_quadrature_rules_impl(
    const std::vector<std::pair<std::string, HOMeshPart<T, I>>>& parts,
    int order,
    bool include_uncut_cells)
{
    if (parts.empty())
        return {};
    if (order < 1 || order > 10)
        throw std::runtime_error("Algoim quadrature order must be in [1, 10]");

    const auto& reference_part = parts.front().second;
    if (!reference_part.mesh)
        throw std::runtime_error("HOMeshPart is not attached to a mesh");
    if (!reference_part.cut_cells || !reference_part.parent_cells)
        throw std::runtime_error("Algoim quadrature requires cut-cell storage");

    const int mesh_tdim = reference_part.mesh->tdim;
    const auto* mesh = reference_part.mesh;
    const auto* cut_cells = reference_part.cut_cells;
    const auto* parent_cells = reference_part.parent_cells;

    for (const auto& [name, part] : parts)
    {
        (void)name;
        if (part.mesh != mesh || part.cut_cells != cut_cells
            || part.parent_cells != parent_cells)
        {
            throw std::runtime_error(
                "Algoim paired quadrature requires parts from the same cut data");
        }
    }

    const auto requests = prepare_paired_requests(parts, mesh_tdim);
    std::vector<std::pair<std::string, quadrature::QuadratureRules<T>>> out;
    out.reserve(parts.size());
    for (const auto& request : requests)
    {
        out.emplace_back(request.name, quadrature::QuadratureRules<T>{});
        initialise_rules(out.back().second, mesh_tdim);
    }

    if (mesh_tdim == 1)
    {
        bool all_interface = true;
        for (const auto& request : requests)
            all_interface = all_interface && request.interface;
        if (all_interface)
        {
            for (std::size_t i = 0; i < requests.size(); ++i)
            {
                out[i].second = interval_interface_quadrature_rules(
                    *requests[i].part, requests[i].level_set_index);
            }
            return out;
        }
    }

    const auto cut_ids = paired_cut_cell_ids(requests);
    for (const std::int32_t cut_id : cut_ids)
    {
        const I parent_cell_id =
            cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_id)];
        const auto& ls_cell = active_level_set_cell(
            *cut_cells, cut_id, requests.front().level_set_index);

        if (!supported_tensor_cell(ls_cell.cell_type))
        {
            throw std::runtime_error(
                "Algoim quadrature currently supports quadrilateral and "
                "hexahedron parent cells only");
        }
        if (ls_cell.bernstein_coeffs.empty())
        {
            throw std::runtime_error(
                "Algoim Bernstein quadrature requires Bernstein level-set data");
        }
        if (ls_cell.tdim != mesh_tdim)
        {
            throw std::runtime_error(
                "Algoim quadrature found inconsistent mesh dimensions");
        }
        std::vector<T> parent_vertex_storage;
        const auto parent_vertices = parent_vertices_for_algoim_cell(
            *mesh, ls_cell, parent_cell_id, parent_vertex_storage);

        if (mesh_tdim == 2)
        {
            append_algoim_bernstein_cell<T, I, 2>(
                out, requests, ls_cell, parent_vertices, parent_cell_id, order);
        }
        else if (mesh_tdim == 3)
        {
            append_algoim_bernstein_cell<T, I, 3>(
                out, requests, ls_cell, parent_vertices, parent_cell_id, order);
        }
        else
        {
            throw std::runtime_error(
                "Algoim quadrature currently supports 2D and 3D only");
        }
    }

    if (include_uncut_cells)
    {
        for (std::size_t i = 0; i < requests.size(); ++i)
        {
            const auto& request = requests[i];
            if (request.interface || request.part->uncut_cell_ids.empty())
                continue;

            HOMeshPart<T, I> uncut_part = *request.part;
            uncut_part.cut_cell_ids.clear();
            const auto uncut_rules = quadrature_rules(
                uncut_part, order, /*include_uncut_cells=*/true);
            merge_rules(out[i].second, uncut_rules);
        }
    }

    return out;
}

#endif // CUTCELLS_HAS_ALGOIM_MULTIPOLY

#endif // CUTCELLS_HAS_ALGOIM

template <std::floating_point T>
T root_only_interval_xtol()
{
    if constexpr (std::is_same_v<T, float>)
        return T(1e-6);
    else
        return T(1e-12);
}

template <std::floating_point T>
void initialise_root_only_rules(quadrature::QuadratureRules<T>& rules, int tdim)
{
    if (rules._tdim == 0)
        rules._tdim = tdim;
    if (rules._offset.empty())
        rules._offset.push_back(0);
}

inline bool root_only_interval_interface_level_set_index(
    const SelectionExpr& expr,
    int& level_set_index)
{
    if (expr.terms.size() != 1 || expr.terms.front().clauses.size() != 1)
        return false;

    const auto& clause = expr.terms.front().clauses.front();
    if (clause.relation != Relation::EqualTo || clause.level_set_index < 0)
        return false;

    level_set_index = clause.level_set_index;
    return true;
}

template <std::floating_point T, std::integral I>
const LevelSetCell<T, I>& root_only_active_level_set_cell(
    const HOCutCells<T, I>& cut_cells,
    int cut_id,
    int level_set_index)
{
    const int begin = cut_cells.ls_offsets.at(static_cast<std::size_t>(cut_id));
    const int end = cut_cells.ls_offsets.at(static_cast<std::size_t>(cut_id + 1));
    for (int i = begin; i < end; ++i)
    {
        const auto& ls_cell =
            cut_cells.level_set_cells[static_cast<std::size_t>(i)];
        if (ls_cell.level_set_id == level_set_index)
            return ls_cell;
    }

    throw std::runtime_error(
        "Interval interface quadrature could not find the selected level set");
}

template <std::floating_point T, std::integral I>
T root_only_interval_value(const LevelSetCell<T, I>& ls_cell, T t)
{
    t = std::clamp(t, T(0), T(1));
    std::array<T, 1> xi = {t};
    return ls_cell.value(std::span<const T>(xi.data(), xi.size()));
}

template <std::floating_point T, std::integral I>
bool root_only_solve_interval_between(const LevelSetCell<T, I>& ls_cell,
                                      T a,
                                      T b,
                                      T& root)
{
    constexpr int max_iter = 64;
    const T xtol = root_only_interval_xtol<T>();
    const T ftol = cell::edge_root::default_value_tolerance<T>();

    a = std::clamp(a, T(0), T(1));
    b = std::clamp(b, T(0), T(1));
    if (b < a)
        std::swap(a, b);

    const T fa = root_only_interval_value(ls_cell, a);
    if (std::abs(fa) <= ftol)
    {
        root = a;
        return true;
    }
    const T fb = root_only_interval_value(ls_cell, b);
    if (std::abs(fb) <= ftol)
    {
        root = b;
        return true;
    }
    if (fa * fb > T(0))
        return false;

    auto phi = [&](std::span<const T> xi) -> T
    {
        return ls_cell.value(xi);
    };
    std::array<T, 1> p0 = {a};
    std::array<T, 1> p1 = {b};
    const auto info = cell::edge_root::find_root_parameter_info<T>(
        std::span<const T>(p0.data(), p0.size()),
        std::span<const T>(p1.data(), p1.size()),
        phi,
        cell::edge_root::method::brent,
        T(0),
        max_iter,
        xtol,
        ftol);
    root = a + info.t * (b - a);
    return true;
}

template <std::floating_point T, std::integral I>
T root_only_polish_interval_root(const LevelSetCell<T, I>& ls_cell, T guess)
{
    const T xtol = root_only_interval_xtol<T>();
    const T ftol = cell::edge_root::default_value_tolerance<T>();
    guess = std::clamp(guess, T(0), T(1));

    if (std::abs(root_only_interval_value(ls_cell, guess)) <= ftol)
        return guess;

    std::vector<T> roots;
    auto add_root = [&](T root)
    {
        root = std::clamp(root, T(0), T(1));
        for (T existing : roots)
        {
            if (std::abs(existing - root) <= T(16) * xtol)
                return;
        }
        roots.push_back(root);
    };

    const int degree = std::max(1, ls_cell.bernstein_order);
    const int intervals = std::max(16, 4 * (degree + 1));
    std::vector<T> grid;
    grid.reserve(static_cast<std::size_t>(intervals + 2));
    for (int i = 0; i <= intervals; ++i)
        grid.push_back(T(i) / T(intervals));
    grid.push_back(guess);
    std::ranges::sort(grid);
    grid.erase(std::unique(
                   grid.begin(), grid.end(),
                   [xtol](T a, T b)
                   { return std::abs(a - b) <= T(16) * xtol; }),
               grid.end());

    for (std::size_t i = 1; i < grid.size(); ++i)
    {
        T root = guess;
        if (root_only_solve_interval_between(
                ls_cell, grid[i - 1], grid[i], root))
        {
            add_root(root);
        }
    }

    if (!roots.empty())
    {
        return *std::min_element(
            roots.begin(), roots.end(),
            [guess](T a, T b)
            {
                return std::abs(a - guess) < std::abs(b - guess);
            });
    }

    T root = guess;
    if (root_only_solve_interval_between(ls_cell, T(0), T(1), root))
        return root;

    return guess;
}

template <std::floating_point T, std::integral I>
bool is_root_only_interval_interface_part(const HOMeshPart<T, I>& part,
                                          int& level_set_index)
{
    if (!part.mesh || !part.cut_cells || !part.parent_cells)
        return false;
    if (part.mesh->tdim != 1 || part.dim != 0)
        return false;
    return root_only_interval_interface_level_set_index(
        part.expr, level_set_index);
}

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> root_only_interval_interface_quadrature_rules(
    const HOMeshPart<T, I>& part,
    int level_set_index)
{
    quadrature::QuadratureRules<T> rules;
    initialise_root_only_rules(rules, part.mesh->tdim);

    const auto infos = selected_zero_entity_infos(part);
    for (const auto& info : infos)
    {
        if (info.dimension != 0)
        {
            throw std::runtime_error(
                "Interval interface quadrature found a non-point zero entity");
        }

        const auto& adapt_cell =
            part.cut_cells->adapt_cells[static_cast<std::size_t>(
                info.cut_cell_id)];
        const int zero_id = static_cast<int>(info.local_zero_entity_id);
        if (zero_id < 0
            || static_cast<std::size_t>(zero_id)
                   >= adapt_cell.zero_entity_dim.size()
            || adapt_cell.zero_entity_dim[static_cast<std::size_t>(zero_id)] != 0)
        {
            throw std::runtime_error(
                "Interval interface quadrature found invalid zero-entity data");
        }

        const int vertex_id =
            adapt_cell.zero_entity_id[static_cast<std::size_t>(zero_id)];
        if (vertex_id < 0
            || static_cast<std::size_t>(vertex_id * adapt_cell.tdim)
                   >= adapt_cell.vertex_coords.size())
        {
            throw std::runtime_error(
                "Interval interface quadrature found invalid root vertex data");
        }

        const auto& ls_cell = root_only_active_level_set_cell(
            *part.cut_cells, info.cut_cell_id, level_set_index);
        if (ls_cell.cell_type != cell::type::interval || ls_cell.tdim != 1)
        {
            throw std::runtime_error(
                "Interval interface quadrature requires interval level-set cells");
        }

        const T guess =
            adapt_cell.vertex_coords[static_cast<std::size_t>(
                vertex_id * adapt_cell.tdim)];
        const T root = root_only_polish_interval_root(ls_cell, guess);
        rules._points.push_back(root);
        rules._weights.push_back(T(1));
        rules._parent_map.push_back(info.parent_cell_id);
        rules._offset.push_back(
            static_cast<std::int32_t>(rules._weights.size()));
    }

    return rules;
}

} // namespace

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> algoim_quadrature_rules(
    const HOMeshPart<T, I>& part,
    int order,
    bool include_uncut_cells)
{
    int level_set_index = -1;
    if (is_root_only_interval_interface_part(part, level_set_index))
        return root_only_interval_interface_quadrature_rules(
            part, level_set_index);

#ifdef CUTCELLS_HAS_ALGOIM_MULTIPOLY
    std::vector<std::pair<std::string, HOMeshPart<T, I>>> parts;
    parts.emplace_back(std::string{}, part);
    auto rules = algoim_paired_quadrature_rules_impl(
        parts, order, include_uncut_cells);
    if (rules.empty())
        return {};
    return std::move(rules.front().second);
#else
    (void)part;
    (void)order;
    (void)include_uncut_cells;
    throw std::runtime_error(
        "CutCells was built without Algoim Bernstein multipoly support. "
        "Rebuild with CUTCELLS_WITH_ALGOIM=ON and a LAPACK library exporting "
        "LAPACKE symbols to use backend='algoim'.");
#endif
}

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> algoim_general_quadrature_rules(
    const HOMeshPart<T, I>& part,
    int order,
    bool include_uncut_cells)
{
    int level_set_index = -1;
    if (is_root_only_interval_interface_part(part, level_set_index))
        return root_only_interval_interface_quadrature_rules(
            part, level_set_index);

#ifdef CUTCELLS_HAS_ALGOIM
    return algoim_general_quadrature_rules_impl(
        part, order, include_uncut_cells);
#else
    (void)part;
    (void)order;
    (void)include_uncut_cells;
    throw std::runtime_error(
        "CutCells was built without Algoim support. Rebuild with "
        "CUTCELLS_WITH_ALGOIM=ON to use backend='algoim_general'.");
#endif
}

template <std::floating_point T, std::integral I>
std::vector<std::pair<std::string, quadrature::QuadratureRules<T>>>
algoim_paired_quadrature_rules(
    const std::vector<std::pair<std::string, HOMeshPart<T, I>>>& parts,
    int order,
    bool include_uncut_cells)
{
    if (!parts.empty())
    {
        bool all_root_only_interval = true;
        std::vector<int> level_set_indices;
        level_set_indices.reserve(parts.size());
        for (const auto& [name, part] : parts)
        {
            (void)name;
            int level_set_index = -1;
            const bool matches =
                is_root_only_interval_interface_part(part, level_set_index);
            all_root_only_interval = all_root_only_interval && matches;
            level_set_indices.push_back(level_set_index);
        }
        if (all_root_only_interval)
        {
            std::vector<std::pair<std::string, quadrature::QuadratureRules<T>>> out;
            out.reserve(parts.size());
            for (std::size_t i = 0; i < parts.size(); ++i)
            {
                out.emplace_back(
                    parts[i].first,
                    root_only_interval_interface_quadrature_rules(
                        parts[i].second, level_set_indices[i]));
            }
            return out;
        }
    }

#ifdef CUTCELLS_HAS_ALGOIM_MULTIPOLY
    return algoim_paired_quadrature_rules_impl(
        parts, order, include_uncut_cells);
#else
    (void)parts;
    (void)order;
    (void)include_uncut_cells;
    throw std::runtime_error(
        "CutCells was built without Algoim Bernstein multipoly support. "
        "Rebuild with CUTCELLS_WITH_ALGOIM=ON and a LAPACK library exporting "
        "LAPACKE symbols to use backend='algoim'.");
#endif
}

template quadrature::QuadratureRules<double> algoim_quadrature_rules(
    const HOMeshPart<double, int>&, int, bool);
template quadrature::QuadratureRules<float> algoim_quadrature_rules(
    const HOMeshPart<float, int>&, int, bool);
template quadrature::QuadratureRules<double> algoim_quadrature_rules(
    const HOMeshPart<double, long>&, int, bool);
template quadrature::QuadratureRules<float> algoim_quadrature_rules(
    const HOMeshPart<float, long>&, int, bool);

template quadrature::QuadratureRules<double> algoim_general_quadrature_rules(
    const HOMeshPart<double, int>&, int, bool);
template quadrature::QuadratureRules<float> algoim_general_quadrature_rules(
    const HOMeshPart<float, int>&, int, bool);
template quadrature::QuadratureRules<double> algoim_general_quadrature_rules(
    const HOMeshPart<double, long>&, int, bool);
template quadrature::QuadratureRules<float> algoim_general_quadrature_rules(
    const HOMeshPart<float, long>&, int, bool);

template std::vector<std::pair<std::string, quadrature::QuadratureRules<double>>>
algoim_paired_quadrature_rules(
    const std::vector<std::pair<std::string, HOMeshPart<double, int>>>&,
    int, bool);
template std::vector<std::pair<std::string, quadrature::QuadratureRules<float>>>
algoim_paired_quadrature_rules(
    const std::vector<std::pair<std::string, HOMeshPart<float, int>>>&,
    int, bool);
template std::vector<std::pair<std::string, quadrature::QuadratureRules<double>>>
algoim_paired_quadrature_rules(
    const std::vector<std::pair<std::string, HOMeshPart<double, long>>>&,
    int, bool);
template std::vector<std::pair<std::string, quadrature::QuadratureRules<float>>>
algoim_paired_quadrature_rules(
    const std::vector<std::pair<std::string, HOMeshPart<float, long>>>&,
    int, bool);

} // namespace cutcells::output
