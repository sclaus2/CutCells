// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "ho_mesh_part_output.h"

#include "mapping.h"
#include "quadrature_tables.h"
#include "reference_cell.h"
#include "cell_topology.h"
#include "quad_midpoint_split.h"
#include "prism_midpoint_split.h"
#include "triangulation.h"
#include "write_vtk.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cutcells::output
{
namespace
{

struct SelectedEntity
{
    cell::type type = cell::type::point;
    std::vector<int> vertices;
    int zero_entity_index = -1;
};

constexpr int vtk_lagrange_curve = 68;
constexpr int vtk_lagrange_triangle = 69;
constexpr int vtk_lagrange_quadrilateral = 70;
constexpr int vtk_lagrange_tetrahedron = 71;
constexpr int vtk_lagrange_wedge = 73;

template <std::floating_point T>
struct CurvedBoundaryEdge
{
    std::array<int, 2> local_vertices = {-1, -1};
    std::array<int, 2> local_global_vertices = {-1, -1};
    std::array<int, 2> zero_vertices = {-1, -1};
    const curving::CurvedZeroEntityState<T>* state = nullptr;
};

template <std::floating_point T>
struct CurvedBoundaryFace
{
    std::array<int, 4> local_vertices = {-1, -1, -1, -1};
    int num_local_vertices = 0;
    std::array<int, 4> zero_vertices = {-1, -1, -1, -1};
    int num_zero_vertices = 0;
    cell::type zero_face_type = cell::type::triangle;
    const curving::CurvedZeroEntityState<T>* state = nullptr;
};

template <std::floating_point T>
struct LocalCurvedSimplexMap
{
    cell::type simplex_type = cell::type::point;
    int dim = 0;
    int parent_tdim = 0;
    int gdim = 0;
    cell::type parent_cell_type = cell::type::point;
    int geometry_order = 1;
    curving::NodeFamily node_family = curving::NodeFamily::lagrange;
    std::vector<int> vertices;
    std::vector<T> ref_vertex_coords;
    std::vector<T> parent_physical_coords;
    std::vector<CurvedBoundaryEdge<T>> curved_edges;
    std::vector<CurvedBoundaryFace<T>> curved_faces;
};

inline bool is_simplex(cell::type cell_type)
{
    return cell_type == cell::type::interval
        || cell_type == cell::type::triangle
        || cell_type == cell::type::tetrahedron;
}

inline cell::type simplex_type_for_dim(int dim)
{
    switch (dim)
    {
        case 1:
            return cell::type::interval;
        case 2:
            return cell::type::triangle;
        case 3:
            return cell::type::tetrahedron;
        default:
            throw std::runtime_error("Unsupported simplex dimension");
    }
}

inline bool supports_curved_lagrange_cell(cell::type cell_type)
{
    return is_simplex(cell_type)
        || cell_type == cell::type::quadrilateral
        || cell_type == cell::type::prism;
}

template <std::floating_point T>
std::pair<T, T> legendre_value_and_derivative(int order, T x)
{
    if (order == 0)
        return {T(1), T(0)};

    T pm2 = T(1);
    T pm1 = x;
    for (int n = 2; n <= order; ++n)
    {
        const T p = ((T(2 * n - 1) * x * pm1) - T(n - 1) * pm2) / T(n);
        pm2 = pm1;
        pm1 = p;
    }

    const T denom = T(1) - x * x;
    if (std::abs(denom) <= T(64) * std::numeric_limits<T>::epsilon())
        return {pm1, T(0)};
    const T derivative = T(order) * (pm2 - x * pm1) / denom;
    return {pm1, derivative};
}

template <std::floating_point T>
std::vector<T> gll_parameters(int order)
{
    order = std::max(order, 1);
    std::vector<T> params(static_cast<std::size_t>(order + 1), T(0));
    params.front() = T(0);
    params.back() = T(1);
    if (order == 1)
        return params;

    const T pi = std::acos(T(-1));
    const T eps = T(128) * std::numeric_limits<T>::epsilon();
    for (int i = 1; i < order; ++i)
    {
        T x = -std::cos(pi * T(i) / T(order));
        for (int iter = 0; iter < 32; ++iter)
        {
            const auto [p, dp] = legendre_value_and_derivative<T>(order, x);
            const T denom = T(1) - x * x;
            if (std::abs(denom) <= eps)
                break;
            const T d2p = (T(2) * x * dp - T(order * (order + 1)) * p) / denom;
            if (std::abs(d2p) <= eps)
                break;
            const T step = dp / d2p;
            x -= step;
            x = std::clamp(x, -T(1) + eps, T(1) - eps);
            if (std::abs(step) <= eps)
                break;
        }
        params[static_cast<std::size_t>(i)] = T(0.5) * (x + T(1));
    }
    return params;
}

template <std::floating_point T>
std::vector<T> interpolation_parameters(int order, curving::NodeFamily family)
{
    order = std::max(order, 1);
    std::vector<T> params(static_cast<std::size_t>(order + 1), T(0));
    if (family == curving::NodeFamily::gll)
        return gll_parameters<T>(order);

    for (int i = 0; i <= order; ++i)
        params[static_cast<std::size_t>(i)] = T(i) / T(order);
    return params;
}

inline curving::NodeFamily construction_node_family(int geometry_order,
                                                    curving::NodeFamily output_family)
{
    (void)geometry_order;
    return output_family;
}

template <std::floating_point T>
std::vector<std::vector<T>> lagrange_simplex_nodes_basix(cell::type simplex_type, int order)
{
    order = std::max(order, 1);
    std::vector<std::vector<T>> nodes;
    if (simplex_type == cell::type::interval)
    {
        nodes.push_back({T(0)});
        nodes.push_back({T(1)});
        for (int i = 1; i < order; ++i)
            nodes.push_back({T(i) / T(order)});
        return nodes;
    }

    if (simplex_type == cell::type::triangle)
    {
        nodes.push_back({T(0), T(0)});
        nodes.push_back({T(1), T(0)});
        nodes.push_back({T(0), T(1)});
        const auto edges = cell::edges(cell::type::triangle);
        const auto verts = cell::reference_vertices<T>(cell::type::triangle);
        for (const auto& edge : edges)
        {
            for (int i = 1; i < order; ++i)
            {
                const T s = T(i) / T(order);
                std::vector<T> x(2, T(0));
                for (int d = 0; d < 2; ++d)
                {
                    const T x0 = verts[static_cast<std::size_t>(edge[0] * 2 + d)];
                    const T x1 = verts[static_cast<std::size_t>(edge[1] * 2 + d)];
                    x[static_cast<std::size_t>(d)] = (T(1) - s) * x0 + s * x1;
                }
                nodes.push_back(std::move(x));
            }
        }
        for (int j = 1; j < order; ++j)
            for (int i = 1; i < order - j; ++i)
                nodes.push_back({T(i) / T(order), T(j) / T(order)});
        return nodes;
    }

    if (simplex_type == cell::type::tetrahedron)
    {
        nodes.push_back({T(0), T(0), T(0)});
        nodes.push_back({T(1), T(0), T(0)});
        nodes.push_back({T(0), T(1), T(0)});
        nodes.push_back({T(0), T(0), T(1)});
        const auto edges = cell::edges(cell::type::tetrahedron);
        const auto verts = cell::reference_vertices<T>(cell::type::tetrahedron);
        for (const auto& edge : edges)
        {
            for (int i = 1; i < order; ++i)
            {
                const T s = T(i) / T(order);
                std::vector<T> x(3, T(0));
                for (int d = 0; d < 3; ++d)
                {
                    const T x0 = verts[static_cast<std::size_t>(edge[0] * 3 + d)];
                    const T x1 = verts[static_cast<std::size_t>(edge[1] * 3 + d)];
                    x[static_cast<std::size_t>(d)] = (T(1) - s) * x0 + s * x1;
                }
                nodes.push_back(std::move(x));
            }
        }
        for (int f = 0; f < 4; ++f)
        {
            const auto fv = cell::face_vertices(cell::type::tetrahedron, f);
            for (int j = 1; j < order; ++j)
            {
                for (int i = 1; i < order - j; ++i)
                {
                    const T w0 = T(1) - T(i + j) / T(order);
                    const T w1 = T(i) / T(order);
                    const T w2 = T(j) / T(order);
                    std::vector<T> x(3, T(0));
                    for (int d = 0; d < 3; ++d)
                    {
                        x[static_cast<std::size_t>(d)] =
                            w0 * verts[static_cast<std::size_t>(fv[0] * 3 + d)]
                          + w1 * verts[static_cast<std::size_t>(fv[1] * 3 + d)]
                          + w2 * verts[static_cast<std::size_t>(fv[2] * 3 + d)];
                    }
                    nodes.push_back(std::move(x));
                }
            }
        }
        for (int k = 1; k < order; ++k)
            for (int j = 1; j < order - k; ++j)
                for (int i = 1; i < order - j - k; ++i)
                    nodes.push_back({T(i) / T(order), T(j) / T(order), T(k) / T(order)});
        return nodes;
    }

    throw std::runtime_error("lagrange_simplex_nodes_basix: unsupported cell type");
}

template <std::floating_point T>
std::vector<std::vector<T>> lagrange_quadrilateral_nodes_vtk(int order)
{
    order = std::max(order, 1);
    std::vector<std::vector<T>> nodes;
    nodes.reserve(static_cast<std::size_t>((order + 1) * (order + 1)));

    nodes.push_back({T(0), T(0)});
    nodes.push_back({T(1), T(0)});
    nodes.push_back({T(1), T(1)});
    nodes.push_back({T(0), T(1)});

    for (int i = 1; i < order; ++i)
        nodes.push_back({T(i) / T(order), T(0)});
    for (int j = 1; j < order; ++j)
        nodes.push_back({T(1), T(j) / T(order)});
    for (int i = 1; i < order; ++i)
        nodes.push_back({T(i) / T(order), T(1)});
    for (int j = 1; j < order; ++j)
        nodes.push_back({T(0), T(j) / T(order)});

    for (int j = 1; j < order; ++j)
        for (int i = 1; i < order; ++i)
            nodes.push_back({T(i) / T(order), T(j) / T(order)});

    return nodes;
}

template <std::floating_point T>
std::vector<std::vector<T>> lagrange_prism_nodes_vtk(int order)
{
    order = std::max(order, 1);
    const int tri_nodes = (order + 1) * (order + 2) / 2;
    std::vector<std::vector<T>> nodes;
    nodes.reserve(static_cast<std::size_t>(tri_nodes * (order + 1)));

    nodes.push_back({T(0), T(0), T(0)});
    nodes.push_back({T(1), T(0), T(0)});
    nodes.push_back({T(0), T(1), T(0)});
    nodes.push_back({T(0), T(0), T(1)});
    nodes.push_back({T(1), T(0), T(1)});
    nodes.push_back({T(0), T(1), T(1)});

    auto append_triangle_edges = [&](T z)
    {
        for (int i = 1; i < order; ++i)
            nodes.push_back({T(i) / T(order), T(0), z});
        for (int r = 1; r < order; ++r)
            nodes.push_back({T(order - r) / T(order), T(r) / T(order), z});
        for (int r = 1; r < order; ++r)
            nodes.push_back({T(0), T(order - r) / T(order), z});
    };

    append_triangle_edges(T(0));
    append_triangle_edges(T(1));

    for (const auto base : {std::array<T, 2>{T(0), T(0)},
                            std::array<T, 2>{T(1), T(0)},
                            std::array<T, 2>{T(0), T(1)}})
    {
        for (int k = 1; k < order; ++k)
            nodes.push_back({base[0], base[1], T(k) / T(order)});
    }

    auto append_triangle_interior = [&](T z)
    {
        for (int j = 1; j < order; ++j)
            for (int i = 1; i < order - j; ++i)
                nodes.push_back({T(i) / T(order), T(j) / T(order), z});
    };

    append_triangle_interior(T(0));
    append_triangle_interior(T(1));

    for (int k = 1; k < order; ++k)
        for (int i = 1; i < order; ++i)
            nodes.push_back({T(i) / T(order), T(0), T(k) / T(order)});
    for (int k = 1; k < order; ++k)
        for (int r = 1; r < order; ++r)
            nodes.push_back({T(order - r) / T(order), T(r) / T(order), T(k) / T(order)});
    for (int k = 1; k < order; ++k)
        for (int r = 1; r < order; ++r)
            nodes.push_back({T(0), T(order - r) / T(order), T(k) / T(order)});

    for (int k = 1; k < order; ++k)
        for (int j = 1; j < order; ++j)
            for (int i = 1; i < order - j; ++i)
                nodes.push_back({T(i) / T(order), T(j) / T(order), T(k) / T(order)});

    return nodes;
}

template <std::floating_point T>
std::vector<std::vector<T>> lagrange_cell_nodes(cell::type cell_type, int order)
{
    if (is_simplex(cell_type))
        return lagrange_simplex_nodes_basix<T>(cell_type, order);
    if (cell_type == cell::type::quadrilateral)
        return lagrange_quadrilateral_nodes_vtk<T>(order);
    if (cell_type == cell::type::prism)
        return lagrange_prism_nodes_vtk<T>(order);
    throw std::runtime_error("lagrange_cell_nodes: unsupported cell type");
}

template <std::floating_point T, std::integral I>
std::vector<T> parent_cell_vertex_coords_vtk(const MeshView<T, I>& mesh, I cell_id)
{
    const auto ctype = mesh.cell_type(cell_id);
    const int nv = cell::get_num_vertices(ctype);
    std::vector<T> coords(static_cast<std::size_t>(nv * mesh.gdim), T(0));

    for (int vtk_v = 0; vtk_v < nv; ++vtk_v)
    {
        const int local_v = mesh.vtk_vertex_order
                                ? vtk_v
                                : cell::vtk_to_basix_vertex(ctype, vtk_v);
        const I node_id = mesh.cell_node(cell_id, static_cast<I>(local_v));
        const T* x = mesh.node(node_id);
        for (int d = 0; d < mesh.gdim; ++d)
            coords[static_cast<std::size_t>(vtk_v * mesh.gdim + d)] = x[d];
    }

    return coords;
}

template <std::floating_point T>
std::vector<T> barycentric_from_simplex_point(cell::type simplex_type,
                                              std::span<const T> xi)
{
    if (simplex_type == cell::type::interval)
        return {T(1) - xi[0], xi[0]};
    if (simplex_type == cell::type::triangle)
        return {T(1) - xi[0] - xi[1], xi[0], xi[1]};
    if (simplex_type == cell::type::tetrahedron)
        return {T(1) - xi[0] - xi[1] - xi[2], xi[0], xi[1], xi[2]};
    throw std::runtime_error("barycentric_from_simplex_point: unsupported simplex");
}

template <std::floating_point T>
std::vector<T> cell_vertex_shape_weights(cell::type cell_type,
                                         std::span<const T> xi)
{
    if (cell_type == cell::type::interval)
        return {T(1) - xi[0], xi[0]};
    if (cell_type == cell::type::triangle)
        return {T(1) - xi[0] - xi[1], xi[0], xi[1]};
    if (cell_type == cell::type::quadrilateral)
    {
        const T u = xi[0];
        const T v = xi[1];
        return {(T(1) - u) * (T(1) - v),
                u * (T(1) - v),
                (T(1) - u) * v,
                u * v};
    }
    if (cell_type == cell::type::tetrahedron)
        return {T(1) - xi[0] - xi[1] - xi[2], xi[0], xi[1], xi[2]};
    if (cell_type == cell::type::prism)
    {
        const T u = xi[0];
        const T v = xi[1];
        const T z = xi[2];
        const T w0 = T(1) - u - v;
        return {w0 * (T(1) - z),
                u * (T(1) - z),
                v * (T(1) - z),
                w0 * z,
                u * z,
                v * z};
    }
    throw std::runtime_error("cell_vertex_shape_weights: unsupported cell type");
}

template <std::floating_point T>
std::vector<T> affine_ref_from_bary(const LocalCurvedSimplexMap<T>& map,
                                    std::span<const T> bary)
{
    std::vector<T> x(static_cast<std::size_t>(map.parent_tdim), T(0));
    for (std::size_t v = 0; v < bary.size(); ++v)
    {
        for (int d = 0; d < map.parent_tdim; ++d)
        {
            x[static_cast<std::size_t>(d)] += bary[v] * map.ref_vertex_coords[
                v * static_cast<std::size_t>(map.parent_tdim) + static_cast<std::size_t>(d)];
        }
    }
    return x;
}

template <std::floating_point T>
std::vector<T> straight_ref_from_cell_point(const LocalCurvedSimplexMap<T>& map,
                                            std::span<const T> xi)
{
    const auto weights = cell_vertex_shape_weights<T>(map.simplex_type, xi);
    std::vector<T> x(static_cast<std::size_t>(map.parent_tdim), T(0));
    for (std::size_t v = 0; v < weights.size(); ++v)
    {
        for (int d = 0; d < map.parent_tdim; ++d)
        {
            x[static_cast<std::size_t>(d)] += weights[v] * map.ref_vertex_coords[
                v * static_cast<std::size_t>(map.parent_tdim) + static_cast<std::size_t>(d)];
        }
    }
    return x;
}

template <std::floating_point T>
std::vector<T> push_parent_ref_to_physical(const LocalCurvedSimplexMap<T>& map,
                                           std::span<const T> ref_point)
{
    const auto phys = cell::push_forward_affine_map<T>(
        map.parent_cell_type,
        map.parent_physical_coords,
        map.gdim,
        ref_point);
    return phys;
}

inline bool same_unordered(std::span<const int> a, std::span<const int> b)
{
    if (a.size() != b.size())
        return false;
    for (const int av : a)
    {
        bool found = false;
        for (const int bv : b)
            found = found || (av == bv);
        if (!found)
            return false;
    }
    return true;
}

template <std::floating_point T>
std::vector<int> zero_entity_vertex_ids(const AdaptCell<T>& ac, int local_zero_entity_id)
{
    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    const int zid = ac.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim == 0)
        return {zid};
    auto verts = ac.entity_to_vertex[zdim][static_cast<std::int32_t>(zid)];
    std::vector<int> out;
    out.reserve(verts.size());
    for (const auto v : verts)
        out.push_back(static_cast<int>(v));
    return out;
}

template <std::floating_point T>
T lagrange_basis_1d(int i, std::span<const T> params, T x)
{
    T value = T(1);
    const T xi = params[static_cast<std::size_t>(i)];
    for (int j = 0; j < static_cast<int>(params.size()); ++j)
    {
        if (j == i)
            continue;
        value *= (x - params[static_cast<std::size_t>(j)])
               / (xi - params[static_cast<std::size_t>(j)]);
    }
    return value;
}

template <std::floating_point T>
T warp_factor(int order, T r)
{
    if (order <= 1)
        return T(0);

    std::vector<T> equispaced(static_cast<std::size_t>(order + 1), T(0));
    std::vector<T> gll = gll_parameters<T>(order);
    for (int i = 0; i <= order; ++i)
    {
        equispaced[static_cast<std::size_t>(i)] = -T(1) + T(2 * i) / T(order);
        gll[static_cast<std::size_t>(i)] = T(2) * gll[static_cast<std::size_t>(i)] - T(1);
    }

    T warp = T(0);
    for (int i = 0; i <= order; ++i)
    {
        const T Li = lagrange_basis_1d<T>(
            i, std::span<const T>(equispaced.data(), equispaced.size()), r);
        warp += Li * (gll[static_cast<std::size_t>(i)]
                    - equispaced[static_cast<std::size_t>(i)]);
    }

    const T edge_factor = T(1) - r * r;
    if (std::abs(edge_factor) > T(64) * std::numeric_limits<T>::epsilon())
        warp /= edge_factor;
    return warp;
}

template <std::floating_point T>
std::array<T, 3> equilateral_to_reference_barycentric(T x, T y)
{
    const T sqrt3 = std::sqrt(T(3));
    std::array<T, 3> w = {};
    w[2] = (sqrt3 * y + T(1)) / T(3);
    w[1] = (x + T(1) - w[2]) / T(2);
    w[0] = T(1) - w[1] - w[2];

    T sum = T(0);
    for (T& wi : w)
    {
        if (std::abs(wi) < T(256) * std::numeric_limits<T>::epsilon())
            wi = T(0);
        if (std::abs(wi - T(1)) < T(256) * std::numeric_limits<T>::epsilon())
            wi = T(1);
        wi = std::clamp(wi, T(0), T(1));
        sum += wi;
    }
    if (sum > T(0))
        for (T& wi : w)
            wi /= sum;
    return w;
}

template <std::floating_point T>
std::vector<std::array<T, 3>> triangle_interpolation_barycentric_nodes(
    int order,
    curving::NodeFamily family)
{
    order = std::max(order, 1);
    std::vector<std::array<T, 3>> nodes;
    nodes.reserve(static_cast<std::size_t>((order + 1) * (order + 2) / 2));

    if (family != curving::NodeFamily::gll)
    {
        for (int j = 0; j <= order; ++j)
        {
            for (int i = 0; i <= order - j; ++i)
            {
                const T u = T(i) / T(order);
                const T v = T(j) / T(order);
                nodes.push_back({T(1) - u - v, u, v});
            }
        }
        return nodes;
    }

    constexpr std::array<T, 16> alpha_opt = {
        T(0.0), T(0.0), T(1.4152), T(0.1001),
        T(0.2751), T(0.9800), T(1.0999), T(1.2832),
        T(1.3648), T(1.4773), T(1.4959), T(1.5743),
        T(1.5770), T(1.6223), T(1.6258), T(1.6530)};
    const T alpha = (order < static_cast<int>(alpha_opt.size()))
                      ? alpha_opt[static_cast<std::size_t>(order)]
                      : T(5) / T(3);
    const T sqrt3 = std::sqrt(T(3));
    const T cos120 = -T(0.5);
    const T sin120 = sqrt3 / T(2);
    const T cos240 = -T(0.5);
    const T sin240 = -sqrt3 / T(2);

    for (int j = 0; j <= order; ++j)
    {
        for (int i = 0; i <= order - j; ++i)
        {
            const T u = T(i) / T(order);
            const T v = T(j) / T(order);
            const T lambda0 = T(1) - u - v;
            const T lambda1 = u;
            const T lambda2 = v;

            T x = -lambda0 + lambda1;
            T y = (-lambda0 - lambda1 + T(2) * lambda2) / sqrt3;

            const T L1 = lambda2;
            const T L2 = lambda0;
            const T L3 = lambda1;
            const T warp1 = T(4) * L2 * L3 * warp_factor<T>(order, L3 - L2)
                          * (T(1) + (alpha * L1) * (alpha * L1));
            const T warp2 = T(4) * L1 * L3 * warp_factor<T>(order, L1 - L3)
                          * (T(1) + (alpha * L2) * (alpha * L2));
            const T warp3 = T(4) * L1 * L2 * warp_factor<T>(order, L2 - L1)
                          * (T(1) + (alpha * L3) * (alpha * L3));

            x += warp1 + cos120 * warp2 + cos240 * warp3;
            y += sin120 * warp2 + sin240 * warp3;
            nodes.push_back(equilateral_to_reference_barycentric<T>(x, y));
        }
    }
    return nodes;
}

template <std::floating_point T>
std::vector<T> triangle_monomials(int order, std::span<const T> bary)
{
    const T u = bary[1];
    const T v = bary[2];
    std::vector<T> values;
    values.reserve(static_cast<std::size_t>((order + 1) * (order + 2) / 2));
    for (int total = 0; total <= order; ++total)
    {
        for (int j = 0; j <= total; ++j)
        {
            const int i = total - j;
            values.push_back(std::pow(u, i) * std::pow(v, j));
        }
    }
    return values;
}

template <std::floating_point T>
bool solve_dense(std::vector<T> A, std::vector<T> b, int n, std::vector<T>& x)
{
    x.assign(static_cast<std::size_t>(n), T(0));
    for (int k = 0; k < n; ++k)
    {
        int pivot = k;
        T best = std::abs(A[static_cast<std::size_t>(k * n + k)]);
        for (int r = k + 1; r < n; ++r)
        {
            const T value = std::abs(A[static_cast<std::size_t>(r * n + k)]);
            if (value > best)
            {
                best = value;
                pivot = r;
            }
        }
        if (best <= T(256) * std::numeric_limits<T>::epsilon())
            return false;
        if (pivot != k)
        {
            for (int c = k; c < n; ++c)
                std::swap(A[static_cast<std::size_t>(k * n + c)],
                          A[static_cast<std::size_t>(pivot * n + c)]);
            std::swap(b[static_cast<std::size_t>(k)], b[static_cast<std::size_t>(pivot)]);
        }

        const T diag = A[static_cast<std::size_t>(k * n + k)];
        for (int c = k; c < n; ++c)
            A[static_cast<std::size_t>(k * n + c)] /= diag;
        b[static_cast<std::size_t>(k)] /= diag;

        for (int r = 0; r < n; ++r)
        {
            if (r == k)
                continue;
            const T factor = A[static_cast<std::size_t>(r * n + k)];
            if (factor == T(0))
                continue;
            for (int c = k; c < n; ++c)
                A[static_cast<std::size_t>(r * n + c)] -=
                    factor * A[static_cast<std::size_t>(k * n + c)];
            b[static_cast<std::size_t>(r)] -= factor * b[static_cast<std::size_t>(k)];
        }
    }
    x = std::move(b);
    return true;
}

template <std::floating_point T>
std::vector<T> triangle_lagrange_basis(int order,
                                       curving::NodeFamily node_family,
                                       std::span<const T> bary)
{
    const auto nodes =
        triangle_interpolation_barycentric_nodes<T>(order, node_family);
    const int n = static_cast<int>(nodes.size());
    std::vector<T> matrix(static_cast<std::size_t>(n * n), T(0));
    for (int row = 0; row < n; ++row)
    {
        const auto mono = triangle_monomials<T>(
            order,
            std::span<const T>(nodes[static_cast<std::size_t>(row)].data(), 3));
        for (int col = 0; col < n; ++col)
            matrix[static_cast<std::size_t>(col * n + row)] =
                mono[static_cast<std::size_t>(col)];
    }

    const auto rhs = triangle_monomials<T>(order, bary);
    std::vector<T> basis;
    if (!solve_dense<T>(std::move(matrix), rhs, n, basis))
        throw std::runtime_error("triangle_lagrange_basis: singular interpolation matrix");
    return basis;
}

template <std::floating_point T>
T simplex_lagrange_factor(int alpha, T lambda, int order)
{
    T value = T(1);
    for (int r = 0; r < alpha; ++r)
        value *= (T(order) * lambda - T(r)) / T(r + 1);
    return value;
}

template <std::floating_point T>
std::vector<T> eval_curved_edge_ref(const CurvedBoundaryEdge<T>& edge,
                                    int order,
                                    curving::NodeFamily family,
                                    T t_local)
{
    const bool same_orientation =
        edge.zero_vertices[0] == edge.local_global_vertices[0]
        && edge.zero_vertices[1] == edge.local_global_vertices[1];
    const bool reverse_orientation =
        edge.zero_vertices[0] == edge.local_global_vertices[1]
        && edge.zero_vertices[1] == edge.local_global_vertices[0];
    if (!same_orientation && !reverse_orientation)
        throw std::runtime_error("eval_curved_edge_ref: edge orientation mismatch");
    const T t = same_orientation ? t_local : T(1) - t_local;
    const auto params = interpolation_parameters<T>(order, family);
    const int tdim = static_cast<int>(edge.state->ref_nodes.size()) / (order + 1);
    std::vector<T> x(static_cast<std::size_t>(tdim), T(0));
    for (int i = 0; i <= order; ++i)
    {
        const T Li = lagrange_basis_1d<T>(i, std::span<const T>(params.data(), params.size()), t);
        for (int d = 0; d < tdim; ++d)
            x[static_cast<std::size_t>(d)] += Li * edge.state->ref_nodes[
                static_cast<std::size_t>(i * tdim + d)];
    }
    return x;
}

template <std::floating_point T>
std::vector<T> eval_curved_face_ref(const CurvedBoundaryFace<T>& face,
                                    int order,
                                    curving::NodeFamily node_family,
                                    std::span<const T> coordinates)
{
    const int nodes_per_face =
        (face.zero_face_type == cell::type::quadrilateral)
            ? (order + 1) * (order + 1)
            : (order + 1) * (order + 2) / 2;
    const int tdim = static_cast<int>(face.state->ref_nodes.size()) / nodes_per_face;
    std::vector<T> x(static_cast<std::size_t>(tdim), T(0));

    if (face.zero_face_type == cell::type::quadrilateral)
    {
        const T u = coordinates[0];
        const T v = coordinates[1];
        const auto params = interpolation_parameters<T>(order, node_family);
        int node = 0;
        for (int j = 0; j <= order; ++j)
        {
            const T Lj = lagrange_basis_1d<T>(
                j, std::span<const T>(params.data(), params.size()), v);
            for (int i = 0; i <= order; ++i)
            {
                const T Li = lagrange_basis_1d<T>(
                    i, std::span<const T>(params.data(), params.size()), u);
                const T L = Li * Lj;
                for (int d = 0; d < tdim; ++d)
                    x[static_cast<std::size_t>(d)] += L * face.state->ref_nodes[
                        static_cast<std::size_t>(node * tdim + d)];
                ++node;
            }
        }
        return x;
    }

    const auto basis = triangle_lagrange_basis<T>(order, node_family, coordinates);
    for (int node = 0; node < static_cast<int>(basis.size()); ++node)
    {
        const T L = basis[static_cast<std::size_t>(node)];
        for (int d = 0; d < tdim; ++d)
            x[static_cast<std::size_t>(d)] += L * face.state->ref_nodes[
                static_cast<std::size_t>(node * tdim + d)];
    }
    return x;
}

template <std::floating_point T>
bool edge_is_in_any_curved_face(const CurvedBoundaryEdge<T>& edge,
                                std::span<const CurvedBoundaryFace<T>> faces)
{
    for (const auto& face : faces)
    {
        bool has0 = false;
        bool has1 = false;
        for (int i = 0; i < face.num_local_vertices; ++i)
        {
            const int fv = face.local_vertices[static_cast<std::size_t>(i)];
            has0 = has0 || (fv == edge.local_vertices[0]);
            has1 = has1 || (fv == edge.local_vertices[1]);
        }
        if (has0 && has1)
            return true;
    }
    return false;
}

template <std::floating_point T>
std::vector<T> curved_map_ref_point(const LocalCurvedSimplexMap<T>& map,
                                    std::span<const T> xi)
{
    const auto vertex_weights = cell_vertex_shape_weights<T>(map.simplex_type, xi);
    std::vector<T> ref = straight_ref_from_cell_point<T>(map, xi);
    constexpr T eps = T(64) * std::numeric_limits<T>::epsilon();

    if (map.dim == 1)
    {
        for (const auto& edge : map.curved_edges)
        {
            const auto curved = eval_curved_edge_ref<T>(
                edge, map.geometry_order, map.node_family, xi[0]);
            for (int d = 0; d < map.parent_tdim; ++d)
                ref[static_cast<std::size_t>(d)] = curved[static_cast<std::size_t>(d)];
        }
        return ref;
    }

    for (const auto& face : map.curved_faces)
    {
        T s = T(0);
        std::array<T, 4> local_w = {};
        for (int i = 0; i < face.num_local_vertices; ++i)
        {
            local_w[static_cast<std::size_t>(i)] =
                vertex_weights[static_cast<std::size_t>(
                    face.local_vertices[static_cast<std::size_t>(i)])];
            s += local_w[static_cast<std::size_t>(i)];
        }
        if (s <= eps)
            continue;

        std::array<T, 3> zero_w = {};
        std::array<T, 2> zero_uv = {};
        for (int zi = 0; zi < face.num_zero_vertices; ++zi)
        {
            for (int li = 0; li < face.num_local_vertices; ++li)
            {
                if (face.zero_vertices[static_cast<std::size_t>(zi)]
                    == map.vertices[static_cast<std::size_t>(
                        face.local_vertices[static_cast<std::size_t>(li)])])
                {
                    const T w = local_w[static_cast<std::size_t>(li)] / s;
                    if (face.zero_face_type == cell::type::quadrilateral)
                    {
                        const T u = (zi == 1 || zi == 3) ? T(1) : T(0);
                        const T v = (zi == 2 || zi == 3) ? T(1) : T(0);
                        zero_uv[0] += w * u;
                        zero_uv[1] += w * v;
                    }
                    else
                    {
                        zero_w[static_cast<std::size_t>(zi)] = w;
                    }
                }
            }
        }

        std::span<const T> face_coordinates =
            (face.zero_face_type == cell::type::quadrilateral)
                ? std::span<const T>(zero_uv.data(), zero_uv.size())
                : std::span<const T>(zero_w.data(), zero_w.size());
        const auto curved = eval_curved_face_ref<T>(
            face, map.geometry_order, map.node_family, face_coordinates);
        std::vector<T> straight(static_cast<std::size_t>(map.parent_tdim), T(0));
        for (int li = 0; li < face.num_local_vertices; ++li)
        {
            const int lv = face.local_vertices[static_cast<std::size_t>(li)];
            const T w = local_w[static_cast<std::size_t>(li)] / s;
            for (int d = 0; d < map.parent_tdim; ++d)
                straight[static_cast<std::size_t>(d)] += w * map.ref_vertex_coords[
                    static_cast<std::size_t>(lv * map.parent_tdim + d)];
        }
        for (int d = 0; d < map.parent_tdim; ++d)
            ref[static_cast<std::size_t>(d)] += s * (curved[static_cast<std::size_t>(d)]
                                                   - straight[static_cast<std::size_t>(d)]);
    }

    for (const auto& edge : map.curved_edges)
    {
        if (edge_is_in_any_curved_face<T>(
                edge, std::span<const CurvedBoundaryFace<T>>(map.curved_faces.data(), map.curved_faces.size())))
        {
            continue;
        }
        const T a = vertex_weights[static_cast<std::size_t>(edge.local_vertices[0])];
        const T b = vertex_weights[static_cast<std::size_t>(edge.local_vertices[1])];
        const T s = a + b;
        if (s <= eps)
            continue;
        const T t = b / s;
        const auto curved = eval_curved_edge_ref<T>(
            edge, map.geometry_order, map.node_family, t);
        std::vector<T> straight(static_cast<std::size_t>(map.parent_tdim), T(0));
        for (int d = 0; d < map.parent_tdim; ++d)
        {
            const T x0 = map.ref_vertex_coords[
                static_cast<std::size_t>(edge.local_vertices[0] * map.parent_tdim + d)];
            const T x1 = map.ref_vertex_coords[
                static_cast<std::size_t>(edge.local_vertices[1] * map.parent_tdim + d)];
            straight[static_cast<std::size_t>(d)] = (T(1) - t) * x0 + t * x1;
            ref[static_cast<std::size_t>(d)] += s * (curved[static_cast<std::size_t>(d)]
                                                   - straight[static_cast<std::size_t>(d)]);
        }
    }

    return ref;
}

template <std::floating_point T>
std::vector<T> curved_map_physical_point(const LocalCurvedSimplexMap<T>& map,
                                         std::span<const T> xi)
{
    const auto ref = curved_map_ref_point<T>(map, xi);
    return push_parent_ref_to_physical<T>(map, std::span<const T>(ref.data(), ref.size()));
}

template <std::floating_point T>
std::vector<T> canonical_simplex_vertices(cell::type simplex_type)
{
    return cell::reference_vertices<T>(simplex_type);
}

template <std::floating_point T>
std::vector<T> map_child_xi_to_parent_xi(cell::type simplex_type,
                                         std::span<const T> child_vertices,
                                         std::span<const T> xi)
{
    const int dim = cell::get_tdim(simplex_type);
    const auto weights = cell_vertex_shape_weights<T>(simplex_type, xi);
    std::vector<T> out(static_cast<std::size_t>(dim), T(0));
    for (std::size_t v = 0; v < weights.size(); ++v)
    {
        for (int d = 0; d < dim; ++d)
        {
            out[static_cast<std::size_t>(d)] += weights[v] * child_vertices[
                v * static_cast<std::size_t>(dim) + static_cast<std::size_t>(d)];
        }
    }
    return out;
}

template <std::floating_point T>
std::vector<T> child_point(cell::type cell_type,
                           std::span<const T> child_vertices,
                           std::initializer_list<T> xi)
{
    std::vector<T> coords(xi);
    return map_child_xi_to_parent_xi<T>(
        cell_type,
        child_vertices,
        std::span<const T>(coords.data(), coords.size()));
}

template <std::floating_point T>
std::vector<std::vector<T>> subdivide_child_simplex(cell::type simplex_type,
                                                    std::span<const T> vertices)
{
    const int dim = cell::get_tdim(simplex_type);
    auto vertex = [&](int i)
    {
        return std::vector<T>(
            vertices.begin() + static_cast<std::ptrdiff_t>(i * dim),
            vertices.begin() + static_cast<std::ptrdiff_t>((i + 1) * dim));
    };
    auto midpoint = [&](const std::vector<T>& a, const std::vector<T>& b)
    {
        std::vector<T> m(static_cast<std::size_t>(dim), T(0));
        for (int d = 0; d < dim; ++d)
            m[static_cast<std::size_t>(d)] = T(0.5) * (a[static_cast<std::size_t>(d)]
                                                     + b[static_cast<std::size_t>(d)]);
        return m;
    };
    auto pack = [&](std::initializer_list<std::vector<T>> pts)
    {
        std::vector<T> child;
        child.reserve(pts.size() * static_cast<std::size_t>(dim));
        for (const auto& p : pts)
            child.insert(child.end(), p.begin(), p.end());
        return child;
    };

    if (simplex_type == cell::type::interval)
    {
        const auto v0 = vertex(0);
        const auto v1 = vertex(1);
        const auto m01 = midpoint(v0, v1);
        return {pack({v0, m01}), pack({m01, v1})};
    }

    if (simplex_type == cell::type::triangle)
    {
        const auto v0 = vertex(0);
        const auto v1 = vertex(1);
        const auto v2 = vertex(2);
        const auto m01 = midpoint(v0, v1);
        const auto m12 = midpoint(v1, v2);
        const auto m20 = midpoint(v2, v0);
        return {
            pack({v0, m01, m20}),
            pack({m01, v1, m12}),
            pack({m20, m12, v2}),
            pack({m01, m12, m20})
        };
    }

    if (simplex_type == cell::type::tetrahedron)
    {
        const auto v0 = vertex(0);
        const auto v1 = vertex(1);
        const auto v2 = vertex(2);
        const auto v3 = vertex(3);
        const auto m01 = midpoint(v0, v1);
        const auto m02 = midpoint(v0, v2);
        const auto m03 = midpoint(v0, v3);
        const auto m12 = midpoint(v1, v2);
        const auto m13 = midpoint(v1, v3);
        const auto m23 = midpoint(v2, v3);
        return {
            pack({v0, m01, m02, m03}),
            pack({m01, v1, m12, m13}),
            pack({m02, m12, v2, m23}),
            pack({m03, m13, m23, v3}),
            pack({m01, m02, m03, m13}),
            pack({m01, m02, m12, m13}),
            pack({m02, m12, m13, m23}),
            pack({m02, m03, m13, m23})
        };
    }

    if (simplex_type == cell::type::quadrilateral)
    {
        return {
            pack({child_point<T>(simplex_type, vertices, {T(0), T(0)}),
                  child_point<T>(simplex_type, vertices, {T(0.5), T(0)}),
                  child_point<T>(simplex_type, vertices, {T(0), T(0.5)}),
                  child_point<T>(simplex_type, vertices, {T(0.5), T(0.5)})}),
            pack({child_point<T>(simplex_type, vertices, {T(0.5), T(0)}),
                  child_point<T>(simplex_type, vertices, {T(1), T(0)}),
                  child_point<T>(simplex_type, vertices, {T(0.5), T(0.5)}),
                  child_point<T>(simplex_type, vertices, {T(1), T(0.5)})}),
            pack({child_point<T>(simplex_type, vertices, {T(0), T(0.5)}),
                  child_point<T>(simplex_type, vertices, {T(0.5), T(0.5)}),
                  child_point<T>(simplex_type, vertices, {T(0), T(1)}),
                  child_point<T>(simplex_type, vertices, {T(0.5), T(1)})}),
            pack({child_point<T>(simplex_type, vertices, {T(0.5), T(0.5)}),
                  child_point<T>(simplex_type, vertices, {T(1), T(0.5)}),
                  child_point<T>(simplex_type, vertices, {T(0.5), T(1)}),
                  child_point<T>(simplex_type, vertices, {T(1), T(1)})})
        };
    }

    if (simplex_type == cell::type::prism)
    {
        const std::array<std::array<T, 2>, 6> tri = {{
            {T(0), T(0)}, {T(1), T(0)}, {T(0), T(1)},
            {T(0.5), T(0)}, {T(0.5), T(0.5)}, {T(0), T(0.5)}
        }};
        const std::array<std::array<int, 3>, 4> sub_tri = {{
            {0, 3, 5},
            {3, 1, 4},
            {5, 4, 2},
            {3, 4, 5}
        }};
        std::vector<std::vector<T>> children;
        children.reserve(8);
        for (int slab = 0; slab < 2; ++slab)
        {
            const T z0 = T(slab) * T(0.5);
            const T z1 = T(slab + 1) * T(0.5);
            for (const auto& st : sub_tri)
            {
                children.push_back(pack({
                    child_point<T>(simplex_type, vertices, {tri[st[0]][0], tri[st[0]][1], z0}),
                    child_point<T>(simplex_type, vertices, {tri[st[1]][0], tri[st[1]][1], z0}),
                    child_point<T>(simplex_type, vertices, {tri[st[2]][0], tri[st[2]][1], z0}),
                    child_point<T>(simplex_type, vertices, {tri[st[0]][0], tri[st[0]][1], z1}),
                    child_point<T>(simplex_type, vertices, {tri[st[1]][0], tri[st[1]][1], z1}),
                    child_point<T>(simplex_type, vertices, {tri[st[2]][0], tri[st[2]][1], z1})
                }));
            }
        }
        return children;
    }

    throw std::runtime_error("subdivide_child_simplex: unsupported simplex type");
}

template <std::floating_point T>
T finite_difference_measure(const LocalCurvedSimplexMap<T>& map,
                            std::span<const T> xi)
{
    const int d = map.dim;
    constexpr T h = T(1e-6);
    std::vector<T> J(static_cast<std::size_t>(map.gdim * d), T(0));
    for (int c = 0; c < d; ++c)
    {
        std::vector<T> xp(xi.begin(), xi.end());
        std::vector<T> xm(xi.begin(), xi.end());
        xp[static_cast<std::size_t>(c)] += h;
        xm[static_cast<std::size_t>(c)] -= h;
        const auto pp = curved_map_physical_point<T>(map, std::span<const T>(xp.data(), xp.size()));
        const auto pm = curved_map_physical_point<T>(map, std::span<const T>(xm.data(), xm.size()));
        for (int r = 0; r < map.gdim; ++r)
            J[static_cast<std::size_t>(c * map.gdim + r)] =
                (pp[static_cast<std::size_t>(r)] - pm[static_cast<std::size_t>(r)]) / (T(2) * h);
    }

    if (d == 1)
    {
        T n2 = T(0);
        for (int r = 0; r < map.gdim; ++r)
            n2 += J[static_cast<std::size_t>(r)] * J[static_cast<std::size_t>(r)];
        return std::sqrt(n2);
    }
    if (d == map.gdim)
    {
        if (d == 2)
            return J[0] * J[3] - J[2] * J[1];
        return J[0] * (J[4] * J[8] - J[7] * J[5])
             - J[3] * (J[1] * J[8] - J[7] * J[2])
             + J[6] * (J[1] * J[5] - J[4] * J[2]);
    }

    T G[9] = {};
    for (int i = 0; i < d; ++i)
        for (int j = 0; j < d; ++j)
            for (int r = 0; r < map.gdim; ++r)
                G[i * d + j] += J[static_cast<std::size_t>(i * map.gdim + r)]
                              * J[static_cast<std::size_t>(j * map.gdim + r)];
    if (d == 2)
        return std::sqrt(std::max(T(0), G[0] * G[3] - G[1] * G[2]));
    return T(0);
}

template <std::floating_point T>
T finite_difference_straight_measure(const LocalCurvedSimplexMap<T>& map,
                                     std::span<const T> xi)
{
    const int d = map.dim;
    constexpr T h = T(1e-6);
    std::vector<T> J(static_cast<std::size_t>(map.gdim * d), T(0));
    auto straight_physical = [&](std::span<const T> xref)
    {
        const auto ref = straight_ref_from_cell_point<T>(map, xref);
        return push_parent_ref_to_physical<T>(
            map, std::span<const T>(ref.data(), ref.size()));
    };

    for (int c = 0; c < d; ++c)
    {
        std::vector<T> xp(xi.begin(), xi.end());
        std::vector<T> xm(xi.begin(), xi.end());
        xp[static_cast<std::size_t>(c)] += h;
        xm[static_cast<std::size_t>(c)] -= h;
        const auto pp = straight_physical(std::span<const T>(xp.data(), xp.size()));
        const auto pm = straight_physical(std::span<const T>(xm.data(), xm.size()));
        for (int r = 0; r < map.gdim; ++r)
            J[static_cast<std::size_t>(c * map.gdim + r)] =
                (pp[static_cast<std::size_t>(r)] - pm[static_cast<std::size_t>(r)]) / (T(2) * h);
    }

    if (d == 1)
    {
        T n2 = T(0);
        for (int r = 0; r < map.gdim; ++r)
            n2 += J[static_cast<std::size_t>(r)] * J[static_cast<std::size_t>(r)];
        return std::sqrt(n2);
    }
    if (d == map.gdim)
    {
        if (d == 2)
            return std::abs(J[0] * J[3] - J[2] * J[1]);
        return std::abs(J[0] * (J[4] * J[8] - J[7] * J[5])
                      - J[3] * (J[1] * J[8] - J[7] * J[2])
                      + J[6] * (J[1] * J[5] - J[4] * J[2]));
    }

    T G[9] = {};
    for (int i = 0; i < d; ++i)
        for (int j = 0; j < d; ++j)
            for (int r = 0; r < map.gdim; ++r)
                G[i * d + j] += J[static_cast<std::size_t>(i * map.gdim + r)]
                              * J[static_cast<std::size_t>(j * map.gdim + r)];
    if (d == 2)
        return std::sqrt(std::max(T(0), G[0] * G[3] - G[1] * G[2]));
    return T(0);
}

template <std::floating_point T>
T finite_difference_child_measure(cell::type cell_type,
                                  std::span<const T> child_vertices,
                                  std::span<const T> xi)
{
    const int d = cell::get_tdim(cell_type);
    constexpr T h = T(1e-6);
    std::vector<T> J(static_cast<std::size_t>(d * d), T(0));
    for (int c = 0; c < d; ++c)
    {
        std::vector<T> xp(xi.begin(), xi.end());
        std::vector<T> xm(xi.begin(), xi.end());
        xp[static_cast<std::size_t>(c)] += h;
        xm[static_cast<std::size_t>(c)] -= h;
        const auto pp = map_child_xi_to_parent_xi<T>(
            cell_type, child_vertices, std::span<const T>(xp.data(), xp.size()));
        const auto pm = map_child_xi_to_parent_xi<T>(
            cell_type, child_vertices, std::span<const T>(xm.data(), xm.size()));
        for (int r = 0; r < d; ++r)
            J[static_cast<std::size_t>(c * d + r)] =
                (pp[static_cast<std::size_t>(r)] - pm[static_cast<std::size_t>(r)]) / (T(2) * h);
    }

    if (d == 1)
        return std::abs(J[0]);
    if (d == 2)
        return std::abs(J[0] * J[3] - J[2] * J[1]);
    return std::abs(J[0] * (J[4] * J[8] - J[7] * J[5])
                  - J[3] * (J[1] * J[8] - J[7] * J[2])
                  + J[6] * (J[1] * J[5] - J[4] * J[2]));
}

template <std::floating_point T>
bool curved_map_valid(const LocalCurvedSimplexMap<T>& map)
{
    std::vector<std::vector<T>> samples;
    if (map.simplex_type == cell::type::interval)
        samples = {{T(0.25)}, {T(0.5)}, {T(0.75)}};
    else if (map.simplex_type == cell::type::triangle)
        samples = {{T(1) / T(3), T(1) / T(3)},
                   {T(0.5), T(0.25)},
                   {T(0.25), T(0.5)},
                   {T(0.25), T(0.25)}};
    else if (map.simplex_type == cell::type::quadrilateral)
        samples = {{T(0.5), T(0.5)},
                   {T(0.25), T(0.25)},
                   {T(0.75), T(0.25)},
                   {T(0.25), T(0.75)},
                   {T(0.75), T(0.75)}};
    else if (map.simplex_type == cell::type::tetrahedron)
        samples = {{T(0.25), T(0.25), T(0.25)},
                   {T(0.5), T(1) / T(6), T(1) / T(6)},
                   {T(1) / T(6), T(0.5), T(1) / T(6)},
                   {T(1) / T(6), T(1) / T(6), T(0.5)}};
    else if (map.simplex_type == cell::type::prism)
        samples = {{T(1) / T(3), T(1) / T(3), T(0.5)},
                   {T(0.2), T(0.2), T(0.25)},
                   {T(0.6), T(0.2), T(0.25)},
                   {T(0.2), T(0.6), T(0.75)},
                   {T(0.2), T(0.2), T(0.75)}};
    else
        return true;

    constexpr T tol = T(1e-12);
    for (const auto& sample : samples)
    {
        const T measure = finite_difference_measure<T>(
            map, std::span<const T>(sample.data(), sample.size()));
        if (!(std::abs(measure) > tol))
            return false;
    }
    return true;
}

template <std::floating_point T>
bool curved_child_map_valid(const LocalCurvedSimplexMap<T>& map,
                            std::span<const T> child_vertices)
{
    if (map.curved_edges.empty() && map.curved_faces.empty())
        return true;

    std::vector<std::vector<T>> samples;
    if (map.simplex_type == cell::type::interval)
        samples = {{T(0.25)}, {T(0.5)}, {T(0.75)}};
    else if (map.simplex_type == cell::type::triangle)
        samples = {{T(1) / T(3), T(1) / T(3)},
                   {T(0.5), T(0.25)},
                   {T(0.25), T(0.5)},
                   {T(0.25), T(0.25)}};
    else if (map.simplex_type == cell::type::quadrilateral)
        samples = {{T(0.5), T(0.5)},
                   {T(0.25), T(0.25)},
                   {T(0.75), T(0.25)},
                   {T(0.25), T(0.75)},
                   {T(0.75), T(0.75)}};
    else if (map.simplex_type == cell::type::tetrahedron)
        samples = {{T(0.25), T(0.25), T(0.25)},
                   {T(0.5), T(1) / T(6), T(1) / T(6)},
                   {T(1) / T(6), T(0.5), T(1) / T(6)},
                   {T(1) / T(6), T(1) / T(6), T(0.5)}};
    else if (map.simplex_type == cell::type::prism)
        samples = {{T(1) / T(3), T(1) / T(3), T(0.5)},
                   {T(0.2), T(0.2), T(0.25)},
                   {T(0.6), T(0.2), T(0.25)},
                   {T(0.2), T(0.6), T(0.75)},
                   {T(0.2), T(0.2), T(0.75)}};
    else
        return true;

    constexpr T tol = T(1e-12);
    for (const auto& sample : samples)
    {
        const auto xi_parent = map_child_xi_to_parent_xi<T>(
            map.simplex_type,
            child_vertices,
            std::span<const T>(sample.data(), sample.size()));
        const T measure = finite_difference_measure<T>(
            map, std::span<const T>(xi_parent.data(), xi_parent.size()));
        if (!(std::abs(measure) > tol))
            return false;
    }
    return true;
}


template <std::floating_point T>
bool vertex_is_zero_for_level_set(const AdaptCell<T>& adapt_cell,
                                  int vertex_id,
                                  int level_set_id)
{
    const std::uint64_t bit = std::uint64_t(1) << level_set_id;
    return (adapt_cell.zero_mask_per_vertex[static_cast<std::size_t>(vertex_id)] & bit) != 0;
}

template <std::floating_point T>
bool cell_contains_all_vertices(std::span<const std::int32_t> cell_verts,
                                std::span<const int> entity_verts)
{
    for (const int v : entity_verts)
    {
        bool found = false;
        for (const auto cv : cell_verts)
        {
            if (cv == v)
            {
                found = true;
                break;
            }
        }
        if (!found)
            return false;
    }
    return true;
}

template <std::floating_point T>
void vertex_state_for_level_set(const AdaptCell<T>& ac,
                                int vertex_id,
                                int level_set_id,
                                bool& is_negative,
                                bool& is_positive,
                                bool& is_zero)
{
    const std::uint64_t bit = std::uint64_t(1) << level_set_id;
    const auto zm = ac.zero_mask_per_vertex[static_cast<std::size_t>(vertex_id)];
    const auto nm = ac.negative_mask_per_vertex[static_cast<std::size_t>(vertex_id)];
    is_zero = (zm & bit) != 0;
    is_negative = !is_zero && ((nm & bit) != 0);
    is_positive = !is_zero && !is_negative;
}

template <std::floating_point T, std::integral I>
bool graph_checks_allow_curving(const HOMeshPart<T, I>& part, int cut_cell_id)
{
    if (!part.cut_cells)
        return true;

    const auto& diagnostics = part.cut_cells->graph_diagnostics;
    if (cut_cell_id < 0
        || cut_cell_id >= static_cast<int>(diagnostics.size()))
    {
        return true;
    }
    return diagnostics[static_cast<std::size_t>(cut_cell_id)].accepted;
}

template <std::floating_point T, std::integral I>
bool graph_checks_allow_zero_entity_curving(const HOMeshPart<T, I>& part,
                                            int cut_cell_id,
                                            int local_zero_entity_id)
{
    if (!part.cut_cells)
        return true;

    const auto& diagnostics = part.cut_cells->graph_diagnostics;
    if (cut_cell_id < 0
        || cut_cell_id >= static_cast<int>(diagnostics.size()))
    {
        return true;
    }

    const auto& cell_diag = diagnostics[static_cast<std::size_t>(cut_cell_id)];
    for (const auto& record : cell_diag.zero_entities)
    {
        if (record.local_zero_entity_id == local_zero_entity_id)
            return record.accepted;
    }

    return true;
}

template <std::floating_point T, std::integral I>
const curving::CurvedZeroEntityState<T>* accepted_curved_state(
    const HOMeshPart<T, I>& part,
    int cut_cell_id,
    int local_zero_entity_id,
    const curving::CurvingOptions<T>& options)
{
    if (!graph_checks_allow_zero_entity_curving<T, I>(
            part, cut_cell_id, local_zero_entity_id))
        return nullptr;

    const auto& state = curving::ensure_curved<T, I>(
        part.cut_cells->curving,
        std::span<const I>(part.cut_cells->parent_cell_ids),
        std::span<const AdaptCell<T>>(part.cut_cells->adapt_cells),
        std::span<const LevelSetCell<T, I>>(part.cut_cells->level_set_cells),
        std::span<const int>(part.cut_cells->ls_offsets),
        cut_cell_id,
        local_zero_entity_id,
        options);
    if (state.status != curving::CurvingStatus::curved)
        return nullptr;
    return &state;
}

template <std::floating_point T>
bool zero_edge_is_same_mask_boundary_of_zero_face(const AdaptCell<T>& ac,
                                                  const LocalCurvedSimplexMap<T>& map,
                                                  int zero_edge_id,
                                                  std::span<const int> edge_vertices)
{
    const auto edge_mask = ac.zero_entity_zero_mask[static_cast<std::size_t>(zero_edge_id)];
    for (int z = 0; z < ac.n_zero_entities(); ++z)
    {
        if (ac.zero_entity_dim[static_cast<std::size_t>(z)] != 2)
            continue;
        if (ac.zero_entity_zero_mask[static_cast<std::size_t>(z)] != edge_mask)
            continue;

        const auto face_vertices = zero_entity_vertex_ids<T>(ac, z);
        bool edge_in_face = true;
        for (const int ev : edge_vertices)
        {
            bool found = false;
            for (const int fv : face_vertices)
                found = found || (ev == fv);
            edge_in_face = edge_in_face && found;
        }
        if (!edge_in_face)
            continue;

        bool face_in_map = true;
        for (const int fv : face_vertices)
        {
            bool found = false;
            for (const int mv : map.vertices)
                found = found || (fv == mv);
            face_in_map = face_in_map && found;
        }
        if (face_in_map)
            return true;
    }
    return false;
}

template <std::floating_point T, std::integral I>
void attach_curved_boundaries(LocalCurvedSimplexMap<T>& map,
                              const HOMeshPart<T, I>& part,
                              const AdaptCell<T>& ac,
                              int cut_cell_id,
                              const curving::CurvingOptions<T>& options,
                              int required_zero_entity_index = -1)
{
    const int nzero = ac.n_zero_entities();
    for (int z = 0; z < nzero; ++z)
    {
        if (required_zero_entity_index >= 0 && z != required_zero_entity_index)
            continue;

        const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(z)];
        if (zdim != 1 && zdim != 2)
            continue;

        const auto zverts = zero_entity_vertex_ids<T>(ac, z);

        if (zdim == 1)
        {
            if (map.parent_tdim == 3 && map.dim >= 2
                && zero_edge_is_same_mask_boundary_of_zero_face<T>(
                    ac, map, z, std::span<const int>(zverts.data(), zverts.size())))
            {
                continue;
            }

            for (const auto& edge : cell::edges(map.simplex_type))
            {
                std::array<int, 2> local_pair = {
                    map.vertices[static_cast<std::size_t>(edge[0])],
                    map.vertices[static_cast<std::size_t>(edge[1])]
                };
                if (!same_unordered(
                        std::span<const int>(local_pair.data(), local_pair.size()),
                        std::span<const int>(zverts.data(), zverts.size())))
                {
                    continue;
                }

                const auto* state = accepted_curved_state<T, I>(part, cut_cell_id, z, options);
                if (state == nullptr)
                    continue;

                CurvedBoundaryEdge<T> ce;
                ce.local_vertices = edge;
                ce.local_global_vertices = local_pair;
                ce.zero_vertices = {zverts[0], zverts[1]};
                ce.state = state;
                map.curved_edges.push_back(ce);
            }
        }
        else if (zdim == 2
                 && (map.simplex_type == cell::type::tetrahedron
                     || map.simplex_type == cell::type::triangle
                     || map.simplex_type == cell::type::quadrilateral
                     || map.simplex_type == cell::type::prism)
                 && (zverts.size() == 3 || zverts.size() == 4))
        {
            auto local_face_is_on_zero_face = [&](std::span<const int> local_face)
            {
                for (const int v : local_face)
                {
                    bool found = false;
                    for (const int zv : zverts)
                        found = found || (v == zv);
                    if (!found)
                        return false;
                }
                return true;
            };

            auto append_face = [&](std::array<int, 4> local_ids, int nlocal)
            {
                std::array<int, 4> local_face = {-1, -1, -1, -1};
                for (int i = 0; i < nlocal; ++i)
                    local_face[static_cast<std::size_t>(i)] =
                        map.vertices[static_cast<std::size_t>(local_ids[static_cast<std::size_t>(i)])];
                if (!local_face_is_on_zero_face(
                        std::span<const int>(local_face.data(), static_cast<std::size_t>(nlocal))))
                {
                    return;
                }

                const auto* state = accepted_curved_state<T, I>(part, cut_cell_id, z, options);
                if (state == nullptr)
                    return;

                CurvedBoundaryFace<T> cf;
                cf.local_vertices = local_ids;
                cf.num_local_vertices = nlocal;
                cf.num_zero_vertices = static_cast<int>(zverts.size());
                cf.zero_face_type = ac.entity_types[2][
                    static_cast<std::size_t>(ac.zero_entity_id[static_cast<std::size_t>(z)])];
                for (std::size_t i = 0; i < zverts.size(); ++i)
                    cf.zero_vertices[i] = zverts[i];
                cf.state = state;
                map.curved_faces.push_back(cf);
            };

            if (map.simplex_type == cell::type::triangle)
            {
                append_face({0, 1, 2, -1}, 3);
            }
            else if (map.simplex_type == cell::type::quadrilateral)
            {
                append_face({0, 1, 2, 3}, 4);
            }
            else if (map.simplex_type == cell::type::prism)
            {
                for (int f = 0; f < cell::num_faces(cell::type::prism); ++f)
                {
                    const auto fv = cell::face_vertices(cell::type::prism, f);
                    std::array<int, 4> local_ids = {-1, -1, -1, -1};
                    for (std::size_t i = 0; i < fv.size(); ++i)
                        local_ids[i] = fv[i];
                    append_face(local_ids, static_cast<int>(fv.size()));
                }
            }
            else
            {
                for (int f = 0; f < cell::num_faces(cell::type::tetrahedron); ++f)
                {
                    const auto fv = cell::face_vertices(cell::type::tetrahedron, f);
                    append_face({fv[0], fv[1], fv[2], -1}, 3);
                }
            }
        }
    }
}

template <std::floating_point T, std::integral I>
bool leaf_cell_matches_sign_requirements(
    const AdaptCell<T>& ac,
    std::span<const std::int32_t> cell_verts,
    const SelectionExpr& expr,
    std::uint64_t cut_cell_active_mask,
    const BackgroundMeshData<T, I>& bg,
    I parent_cell_id)
{
    const int nls = std::min(bg.num_level_sets, 64);
    for (int li = 0; li < nls; ++li)
    {
        const std::uint64_t bit = std::uint64_t(1) << li;
        const bool require_neg = (expr.negative_required & bit) != 0;
        const bool require_pos = (expr.positive_required & bit) != 0;
        if (!require_neg && !require_pos)
            continue;

        const bool ls_is_active = (cut_cell_active_mask & bit) != 0;
        if (!ls_is_active)
        {
            const auto dom = bg.domain(li, parent_cell_id);
            if (require_neg && dom != cell::domain::inside)
                return false;
            if (require_pos && dom != cell::domain::outside)
                return false;
            continue;
        }

        bool has_neg = false;
        bool has_pos = false;
        for (const auto cv : cell_verts)
        {
            bool is_neg = false;
            bool is_pos = false;
            bool is_zero = false;
            vertex_state_for_level_set(ac, static_cast<int>(cv), li, is_neg, is_pos, is_zero);
            has_neg = has_neg || is_neg;
            has_pos = has_pos || is_pos;
        }

        if (require_neg)
        {
            if (has_pos || !has_neg)
                return false;
        }
        if (require_pos)
        {
            if (has_neg || !has_pos)
                return false;
        }
    }

    return true;
}

template <std::floating_point T, std::integral I>
bool zero_entity_matches(const HOMeshPart<T, I>& part,
                         const AdaptCell<T>& ac,
                         int zero_entity_index,
                         std::uint64_t cut_cell_active_mask,
                         I parent_cell_id)
{
    if (ac.zero_entity_dim[static_cast<std::size_t>(zero_entity_index)] != part.dim)
        return false;

    const auto zero_mask = ac.zero_entity_zero_mask[static_cast<std::size_t>(zero_entity_index)];
    if ((zero_mask & part.expr.zero_required) != part.expr.zero_required)
        return false;

    if (part.expr.negative_required == 0 && part.expr.positive_required == 0)
        return true;

    std::vector<int> zero_verts;
    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(zero_entity_index)];
    const int zid = ac.zero_entity_id[static_cast<std::size_t>(zero_entity_index)];
    if (zdim == 0)
    {
        zero_verts.push_back(zid);
    }
    else
    {
        auto verts = ac.entity_to_vertex[zdim][static_cast<std::int32_t>(zid)];
        zero_verts.reserve(verts.size());
        for (const auto v : verts)
            zero_verts.push_back(static_cast<int>(v));
    }

    const int tdim = ac.tdim;
    const int n_cells = ac.n_entities(tdim);
    for (int c = 0; c < n_cells; ++c)
    {
        auto cell_verts = ac.entity_to_vertex[tdim][static_cast<std::int32_t>(c)];
        if (!cell_contains_all_vertices<T>(cell_verts, std::span<const int>(zero_verts)))
            continue;

        if (leaf_cell_matches_sign_requirements(
                ac, cell_verts, part.expr, cut_cell_active_mask, *part.bg, parent_cell_id))
        {
            return true;
        }
    }

    return false;
}

template <std::floating_point T, std::integral I>
std::vector<SelectedEntity> selected_entities(const HOMeshPart<T, I>& part,
                                              const AdaptCell<T>& adapt_cell,
                                              int cut_cell_id)
{
    std::vector<SelectedEntity> entities;
    const I parent_cell_id =
        part.cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_cell_id)];
    const std::uint64_t cut_active_mask =
        part.cut_cells->active_level_set_mask[static_cast<std::size_t>(cut_cell_id)];

    if (part.dim == adapt_cell.tdim)
    {
        const int n_cells = adapt_cell.n_entities(adapt_cell.tdim);
        entities.reserve(static_cast<std::size_t>(n_cells));
        for (int c = 0; c < n_cells; ++c)
        {
            auto verts = adapt_cell.entity_to_vertex[adapt_cell.tdim][static_cast<std::int32_t>(c)];
            if (!leaf_cell_matches_sign_requirements(
                    adapt_cell, verts, part.expr, cut_active_mask, *part.bg, parent_cell_id))
            {
                continue;
            }

            SelectedEntity entity;
            entity.type = adapt_cell.entity_types[adapt_cell.tdim][static_cast<std::size_t>(c)];
            entity.vertices.assign(verts.begin(), verts.end());
            entities.push_back(std::move(entity));
        }
        return entities;
    }

    if ((cut_active_mask & part.expr.zero_required) != part.expr.zero_required)
        return entities;

    const int n_zero = adapt_cell.n_zero_entities();
    entities.reserve(static_cast<std::size_t>(n_zero));
    for (int z = 0; z < n_zero; ++z)
    {
        if (!zero_entity_matches(part, adapt_cell, z, cut_active_mask, parent_cell_id))
            continue;

        const int zdim = adapt_cell.zero_entity_dim[static_cast<std::size_t>(z)];
        const int zid = adapt_cell.zero_entity_id[static_cast<std::size_t>(z)];
        SelectedEntity entity;
        entity.type = (zdim == 0)
                          ? cell::type::point
                          : adapt_cell.entity_types[zdim][static_cast<std::size_t>(zid)];
        entity.zero_entity_index = z;
        if (zdim == 0)
        {
            entity.vertices.push_back(zid);
        }
        else
        {
            auto verts = adapt_cell.entity_to_vertex[zdim][static_cast<std::int32_t>(zid)];
            entity.vertices.assign(verts.begin(), verts.end());
        }

        entities.push_back(std::move(entity));
    }

    return entities;
}

} // namespace

template <std::floating_point T, std::integral I>
std::vector<SelectedZeroEntityInfo> selected_zero_entity_infos(
    const HOMeshPart<T, I>& part)
{
    if (!part.cut_cells)
        throw std::runtime_error("HOMeshPart is not attached to cut-cell storage");

    std::vector<SelectedZeroEntityInfo> out;
    for (std::int32_t cut_id : part.cut_cell_ids)
    {
        if (cut_id < 0
            || cut_id >= static_cast<std::int32_t>(part.cut_cells->adapt_cells.size()))
        {
            continue;
        }

        const auto& adapt_cell =
            part.cut_cells->adapt_cells[static_cast<std::size_t>(cut_id)];
        const auto entities = selected_entities(part, adapt_cell, cut_id);
        const auto parent_cell_id =
            part.cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_id)];
        for (const auto& entity : entities)
        {
            if (entity.zero_entity_index < 0)
                continue;

            SelectedZeroEntityInfo info;
            info.cut_cell_id = cut_id;
            info.parent_cell_id = static_cast<std::int32_t>(parent_cell_id);
            info.local_zero_entity_id =
                static_cast<std::int32_t>(entity.zero_entity_index);
            info.dimension = static_cast<std::int32_t>(
                adapt_cell.zero_entity_dim[
                    static_cast<std::size_t>(entity.zero_entity_index)]);
            out.push_back(info);
        }
    }
    return out;
}

namespace
{

template <std::floating_point T>
std::vector<T> entity_reference_coords(const AdaptCell<T>& adapt_cell,
                                       std::span<const int> entity_vertices)
{
    std::vector<T> coords(
        static_cast<std::size_t>(entity_vertices.size() * adapt_cell.tdim), T(0));

    for (std::size_t j = 0; j < entity_vertices.size(); ++j)
    {
        const int gv = entity_vertices[j];
        for (int d = 0; d < adapt_cell.tdim; ++d)
        {
            coords[static_cast<std::size_t>(j * adapt_cell.tdim + d)] =
                adapt_cell.vertex_coords[static_cast<std::size_t>(gv * adapt_cell.tdim + d)];
        }
    }

    return coords;
}

template <std::floating_point T>
void gather_subcell_vertices(std::span<const T> coords,
                             int coord_dim,
                             std::span<const int> vertex_ids,
                             std::vector<T>& out)
{
    out.resize(static_cast<std::size_t>(vertex_ids.size() * coord_dim));
    for (std::size_t j = 0; j < vertex_ids.size(); ++j)
    {
        const int local_v = vertex_ids[j];
        for (int d = 0; d < coord_dim; ++d)
        {
            out[static_cast<std::size_t>(j * coord_dim + d)] =
                coords[static_cast<std::size_t>(local_v * coord_dim + d)];
        }
    }
}

template <std::floating_point T>
void map_canonical_to_subcell_points(const T* canonical_points,
                                     int num_points,
                                     int simplex_dim,
                                     const T* subcell_vertices,
                                     int parent_tdim,
                                     T* out_points)
{
    const T* v0 = subcell_vertices;

    for (int q = 0; q < num_points; ++q)
    {
        const T* X = canonical_points + q * simplex_dim;
        T* x = out_points + q * parent_tdim;

        for (int d = 0; d < parent_tdim; ++d)
            x[d] = v0[d];

        for (int i = 1; i <= simplex_dim; ++i)
        {
            const T* vi = subcell_vertices + i * parent_tdim;
            for (int d = 0; d < parent_tdim; ++d)
                x[d] += X[i - 1] * (vi[d] - v0[d]);
        }
    }
}

template <std::floating_point T>
T simplex_physical_measure(const T* vertices,
                           int simplex_dim,
                           int gdim)
{
    T J[9] = {};
    const T* v0 = vertices;

    for (int col = 0; col < simplex_dim; ++col)
    {
        const T* vi = vertices + (col + 1) * gdim;
        for (int row = 0; row < gdim; ++row)
            J[col * gdim + row] = vi[row] - v0[row];
    }

    if (simplex_dim == gdim)
    {
        if (simplex_dim == 1)
            return std::abs(J[0]);
        if (simplex_dim == 2)
            return std::abs(J[0] * J[3] - J[2] * J[1]);

        const T det =
            J[0] * (J[4] * J[8] - J[7] * J[5])
          - J[3] * (J[1] * J[8] - J[7] * J[2])
          + J[6] * (J[1] * J[5] - J[4] * J[2]);
        return std::abs(det);
    }

    T G[9] = {};
    for (int i = 0; i < simplex_dim; ++i)
    {
        for (int j = 0; j < simplex_dim; ++j)
        {
            T sum = 0;
            for (int k = 0; k < gdim; ++k)
                sum += J[i * gdim + k] * J[j * gdim + k];
            G[i * simplex_dim + j] = sum;
        }
    }

    if (simplex_dim == 1)
        return std::sqrt(G[0]);
    if (simplex_dim == 2)
        return std::sqrt(G[0] * G[3] - G[1] * G[2]);

    const T det =
        G[0] * (G[4] * G[8] - G[7] * G[5])
      - G[3] * (G[1] * G[8] - G[7] * G[2])
      + G[6] * (G[1] * G[5] - G[4] * G[2]);
    return std::sqrt(det);
}

template <std::floating_point T>
std::vector<T> reorder_vertex_coords_to_vtk(cell::type cell_type,
                                            std::span<const T> coords,
                                            int coord_dim)
{
    const auto perm = cell::basix_to_vtk_vertex_permutation(cell_type);
    return cell::permute_vertex_data(coords, coord_dim, perm);
}

inline std::vector<int> vtk_local_ids_from_basix(cell::type cell_type)
{
    const auto perm = cell::vtk_to_basix_vertex_permutation(cell_type);
    return std::vector<int>(perm.begin(), perm.end());
}

template <std::floating_point T>
void append_mesh_entity(mesh::CutMesh<T>& out,
                        std::span<const T> physical_coords,
                        int gdim,
                        cell::type cell_type,
                        int parent_cell_id,
                        bool triangulate,
                        bool input_is_basix,
                        cell::type source_parent_type,
                        std::span<const int> root_vertex_flags)
{
    if (out._gdim == 0)
        out._gdim = gdim;
    if (out._tdim == 0)
        out._tdim = cell::get_tdim(cell_type);

    // CutMesh stores vertex coordinates in basix ordering internally.
    // Basix-ordered inputs are stored as-is; VTK-ordered inputs are
    // permuted to basix ordering before storage.
    std::vector<T> reordered_coords;
    std::span<const T> output_coords = physical_coords;
    if (!input_is_basix && !is_simplex(cell_type))
    {
        const auto perm = cell::vtk_to_basix_vertex_permutation(cell_type);
        reordered_coords = cell::permute_vertex_data(physical_coords, gdim, perm);
        output_coords = std::span<const T>(reordered_coords.data(), reordered_coords.size());
    }

    const int nv = static_cast<int>(output_coords.size()) / gdim;
    const int vertex_base = out._num_vertices;

    if (triangulate
        && cell_type == cell::type::quadrilateral
        && source_parent_type == cell::type::triangle
        && root_vertex_flags.size() == static_cast<std::size_t>(nv))
    {
        std::array<int, 4> quad_tokens = {};
        for (int i = 0; i < nv; ++i)
        {
            quad_tokens[static_cast<std::size_t>(i)]
                = root_vertex_flags[static_cast<std::size_t>(i)] ? i : 100 + i;
        }

        int next_token_base = 200;
        auto split = cell::quad_midpoint::split_triangle_derived_quadrilateral<T>(
            output_coords,
            gdim,
            std::span<const int>(quad_tokens.data(), quad_tokens.size()),
            next_token_base);

        out._vertex_coords.insert(
            out._vertex_coords.end(), output_coords.begin(), output_coords.end());
        out._num_vertices += nv;

        const int midpoint_vertex_base = out._num_vertices;
        out._vertex_coords.insert(
            out._vertex_coords.end(),
            split.added_vertex_coords.begin(),
            split.added_vertex_coords.end());
        out._num_vertices += static_cast<int>(split.added_vertex_tokens.size());

        std::array<int, 256> token_to_vertex;
        token_to_vertex.fill(-1);
        for (int i = 0; i < nv; ++i)
            token_to_vertex[quad_tokens[static_cast<std::size_t>(i)]] = vertex_base + i;
        for (std::size_t i = 0; i < split.added_vertex_tokens.size(); ++i)
        {
            token_to_vertex[split.added_vertex_tokens[i]]
                = midpoint_vertex_base + static_cast<int>(i);
        }

        for (const auto& tri_tokens : split.triangles)
        {
            for (int k = 0; k < 3; ++k)
                out._connectivity.push_back(token_to_vertex[tri_tokens[static_cast<std::size_t>(k)]]);
            out._offset.push_back(static_cast<int>(out._connectivity.size()));
            out._types.push_back(cell::type::triangle);
            out._parent_map.push_back(parent_cell_id);
            out._num_cells += 1;
        }
        return;
    }

    if (triangulate
        && cell_type == cell::type::prism
        && source_parent_type == cell::type::tetrahedron
        && root_vertex_flags.size() == static_cast<std::size_t>(nv))
    {
        std::array<int, 6> prism_tokens = {};
        for (int i = 0; i < nv; ++i)
            prism_tokens[static_cast<std::size_t>(i)] = root_vertex_flags[static_cast<std::size_t>(i)] ? i : 100 + i;

        int next_token_base = 200;
        auto split = cell::prism_midpoint::split_tetra_derived_prism<T>(
            output_coords,
            gdim,
            std::span<const int>(prism_tokens.data(), prism_tokens.size()),
            next_token_base);

        out._vertex_coords.insert(
            out._vertex_coords.end(), output_coords.begin(), output_coords.end());
        out._num_vertices += nv;

        const int midpoint_vertex_base = out._num_vertices;
        out._vertex_coords.insert(
            out._vertex_coords.end(),
            split.added_vertex_coords.begin(),
            split.added_vertex_coords.end());
        out._num_vertices += static_cast<int>(split.added_vertex_tokens.size());

        std::array<int, 256> token_to_vertex;
        token_to_vertex.fill(-1);
        for (int i = 0; i < nv; ++i)
            token_to_vertex[prism_tokens[static_cast<std::size_t>(i)]] = vertex_base + i;
        for (std::size_t i = 0; i < split.added_vertex_tokens.size(); ++i)
            token_to_vertex[split.added_vertex_tokens[i]] = midpoint_vertex_base + static_cast<int>(i);

        for (const auto& tet_tokens : split.tets)
        {
            for (int k = 0; k < 4; ++k)
                out._connectivity.push_back(token_to_vertex[tet_tokens[static_cast<std::size_t>(k)]]);
            out._offset.push_back(static_cast<int>(out._connectivity.size()));
            out._types.push_back(cell::type::tetrahedron);
            out._parent_map.push_back(parent_cell_id);
            out._num_cells += 1;
        }
        return;
    }

    out._vertex_coords.insert(
        out._vertex_coords.end(), output_coords.begin(), output_coords.end());
    out._num_vertices += nv;

    if (triangulate && !is_simplex(cell_type) && cell::get_tdim(cell_type) >= 2)
    {
        std::vector<int> local_ids(static_cast<std::size_t>(nv));
        std::iota(local_ids.begin(), local_ids.end(), 0);

        std::vector<std::vector<int>> simplices;
        cell::triangulation(cell_type, local_ids.data(), simplices);
        const auto simplex_type = simplex_type_for_dim(cell::get_tdim(cell_type));

        for (const auto& simplex : simplices)
        {
            for (int lv : simplex)
                out._connectivity.push_back(vertex_base + lv);
            out._offset.push_back(static_cast<int>(out._connectivity.size()));
            out._types.push_back(simplex_type);
            out._parent_map.push_back(parent_cell_id);
            out._num_cells += 1;
        }
        return;
    }

    for (int lv = 0; lv < nv; ++lv)
        out._connectivity.push_back(vertex_base + lv);
    out._offset.push_back(static_cast<int>(out._connectivity.size()));
    out._types.push_back(cell_type);
    out._parent_map.push_back(parent_cell_id);
    out._num_cells += 1;
}

template <std::floating_point T>
void append_simplex_quadrature(quadrature::QuadratureRules<T>& rules,
                               cell::type simplex_type,
                               std::span<const T> ref_vertices,
                               std::span<const T> physical_vertices,
                               int parent_tdim,
                               int gdim,
                               int order)
{
    const auto ref_rule = quadrature::get_reference_rule<T>(simplex_type, order);
    const int num_points = ref_rule._num_points;
    const int simplex_dim = ref_rule._tdim;

    std::vector<T> mapped_ref_points(
        static_cast<std::size_t>(num_points * parent_tdim), T(0));
    map_canonical_to_subcell_points(
        ref_rule._points.data(),
        num_points,
        simplex_dim,
        ref_vertices.data(),
        parent_tdim,
        mapped_ref_points.data());

    rules._points.insert(
        rules._points.end(), mapped_ref_points.begin(), mapped_ref_points.end());

    const T measure = simplex_physical_measure(
        physical_vertices.data(), simplex_dim, gdim);
    for (int q = 0; q < num_points; ++q)
        rules._weights.push_back(ref_rule._weights[q] * measure);
}

template <std::floating_point T>
void append_entity_quadrature(quadrature::QuadratureRules<T>& rules,
                              cell::type cell_type,
                              std::span<const T> ref_vertices,
                              std::span<const T> physical_vertices,
                              int parent_tdim,
                              int gdim,
                              int parent_cell_id,
                              int order,
                              bool triangulate,
                              bool input_is_basix,
                              cell::type source_parent_type,
                              std::span<const int> root_vertex_flags)
{
    if (rules._tdim == 0)
        rules._tdim = parent_tdim;
    if (rules._offset.empty())
        rules._offset.push_back(0);

    std::span<const T> ref_use = ref_vertices;
    std::span<const T> phys_use = physical_vertices;
    std::vector<T> ref_vertices_basix;
    std::vector<T> phys_vertices_basix;
    if (triangulate && !input_is_basix && !is_simplex(cell_type))
    {
        const auto perm = cell::vtk_to_basix_vertex_permutation(cell_type);
        ref_vertices_basix = cell::permute_vertex_data(ref_vertices, parent_tdim, perm);
        phys_vertices_basix = cell::permute_vertex_data(physical_vertices, gdim, perm);
        ref_use = std::span<const T>(ref_vertices_basix.data(), ref_vertices_basix.size());
        phys_use = std::span<const T>(phys_vertices_basix.data(), phys_vertices_basix.size());
    }

    const int entity_dim = cell::get_tdim(cell_type);
    const int nv = static_cast<int>(ref_use.size()) / parent_tdim;

    if (triangulate
        && cell_type == cell::type::quadrilateral
        && source_parent_type == cell::type::triangle
        && root_vertex_flags.size() == static_cast<std::size_t>(nv))
    {
        std::array<int, 4> quad_tokens = {};
        for (int i = 0; i < nv; ++i)
        {
            quad_tokens[static_cast<std::size_t>(i)]
                = root_vertex_flags[static_cast<std::size_t>(i)] ? i : 100 + i;
        }

        int next_ref_token = 200;
        auto ref_split = cell::quad_midpoint::split_triangle_derived_quadrilateral<T>(
            ref_use,
            parent_tdim,
            std::span<const int>(quad_tokens.data(), quad_tokens.size()),
            next_ref_token);

        int next_phys_token = 200;
        auto phys_split = cell::quad_midpoint::split_triangle_derived_quadrilateral<T>(
            phys_use,
            gdim,
            std::span<const int>(quad_tokens.data(), quad_tokens.size()),
            next_phys_token);

        std::array<int, 256> token_to_local;
        token_to_local.fill(-1);
        for (int i = 0; i < nv; ++i)
            token_to_local[quad_tokens[static_cast<std::size_t>(i)]] = i;
        for (std::size_t i = 0; i < ref_split.added_vertex_tokens.size(); ++i)
            token_to_local[ref_split.added_vertex_tokens[i]] = nv + static_cast<int>(i);

        std::vector<T> ref_all(ref_use.begin(), ref_use.end());
        ref_all.insert(
            ref_all.end(),
            ref_split.added_vertex_coords.begin(),
            ref_split.added_vertex_coords.end());

        std::vector<T> phys_all(phys_use.begin(), phys_use.end());
        phys_all.insert(
            phys_all.end(),
            phys_split.added_vertex_coords.begin(),
            phys_split.added_vertex_coords.end());

        std::vector<T> ref_simplex;
        std::vector<T> phys_simplex;
        for (const auto& tri_tokens : ref_split.triangles)
        {
            const std::array<int, 3> local_ids = {
                token_to_local[tri_tokens[0]],
                token_to_local[tri_tokens[1]],
                token_to_local[tri_tokens[2]],
            };
            gather_subcell_vertices(
                std::span<const T>(ref_all.data(), ref_all.size()),
                parent_tdim,
                std::span<const int>(local_ids.data(), local_ids.size()),
                ref_simplex);
            gather_subcell_vertices(
                std::span<const T>(phys_all.data(), phys_all.size()),
                gdim,
                std::span<const int>(local_ids.data(), local_ids.size()),
                phys_simplex);
            append_simplex_quadrature(
                rules,
                cell::type::triangle,
                std::span<const T>(ref_simplex.data(), ref_simplex.size()),
                std::span<const T>(phys_simplex.data(), phys_simplex.size()),
                parent_tdim,
                gdim,
                order);
        }

        rules._parent_map.push_back(parent_cell_id);
        rules._offset.push_back(static_cast<int32_t>(rules._weights.size()));
        return;
    }

    if (triangulate
        && cell_type == cell::type::prism
        && source_parent_type == cell::type::tetrahedron
        && root_vertex_flags.size() == static_cast<std::size_t>(nv))
    {
        std::array<int, 6> prism_tokens = {};
        for (int i = 0; i < nv; ++i)
        {
            prism_tokens[static_cast<std::size_t>(i)]
                = root_vertex_flags[static_cast<std::size_t>(i)] ? i : 100 + i;
        }

        int next_ref_token = 200;
        auto ref_split = cell::prism_midpoint::split_tetra_derived_prism<T>(
            ref_use,
            parent_tdim,
            std::span<const int>(prism_tokens.data(), prism_tokens.size()),
            next_ref_token);

        int next_phys_token = 200;
        auto phys_split = cell::prism_midpoint::split_tetra_derived_prism<T>(
            phys_use,
            gdim,
            std::span<const int>(prism_tokens.data(), prism_tokens.size()),
            next_phys_token);

        std::array<int, 256> token_to_local;
        token_to_local.fill(-1);
        for (int i = 0; i < nv; ++i)
            token_to_local[prism_tokens[static_cast<std::size_t>(i)]] = i;
        for (std::size_t i = 0; i < ref_split.added_vertex_tokens.size(); ++i)
            token_to_local[ref_split.added_vertex_tokens[i]] = nv + static_cast<int>(i);

        std::vector<T> ref_all(ref_use.begin(), ref_use.end());
        ref_all.insert(
            ref_all.end(),
            ref_split.added_vertex_coords.begin(),
            ref_split.added_vertex_coords.end());

        std::vector<T> phys_all(phys_use.begin(), phys_use.end());
        phys_all.insert(
            phys_all.end(),
            phys_split.added_vertex_coords.begin(),
            phys_split.added_vertex_coords.end());

        std::vector<T> ref_simplex;
        std::vector<T> phys_simplex;
        for (const auto& tet_tokens : ref_split.tets)
        {
            const std::array<int, 4> local_ids = {
                token_to_local[tet_tokens[0]],
                token_to_local[tet_tokens[1]],
                token_to_local[tet_tokens[2]],
                token_to_local[tet_tokens[3]],
            };
            gather_subcell_vertices(
                std::span<const T>(ref_all.data(), ref_all.size()),
                parent_tdim,
                std::span<const int>(local_ids.data(), local_ids.size()),
                ref_simplex);
            gather_subcell_vertices(
                std::span<const T>(phys_all.data(), phys_all.size()),
                gdim,
                std::span<const int>(local_ids.data(), local_ids.size()),
                phys_simplex);
            append_simplex_quadrature(
                rules,
                cell::type::tetrahedron,
                std::span<const T>(ref_simplex.data(), ref_simplex.size()),
                std::span<const T>(phys_simplex.data(), phys_simplex.size()),
                parent_tdim,
                gdim,
                order);
        }

        rules._parent_map.push_back(parent_cell_id);
        rules._offset.push_back(static_cast<int32_t>(rules._weights.size()));
        return;
    }

    if (triangulate && !is_simplex(cell_type) && entity_dim >= 2)
    {
        std::vector<int> local_ids(static_cast<std::size_t>(nv));
        std::iota(local_ids.begin(), local_ids.end(), 0);

        std::vector<std::vector<int>> simplices;
        cell::triangulation(cell_type, local_ids.data(), simplices);
        const auto simplex_type = simplex_type_for_dim(entity_dim);

        std::vector<T> ref_simplex;
        std::vector<T> phys_simplex;
        for (const auto& simplex : simplices)
        {
            gather_subcell_vertices(
                ref_use, parent_tdim,
                std::span<const int>(simplex.data(), simplex.size()),
                ref_simplex);
            gather_subcell_vertices(
                phys_use, gdim,
                std::span<const int>(simplex.data(), simplex.size()),
                phys_simplex);
            append_simplex_quadrature(
                rules,
                simplex_type,
                std::span<const T>(ref_simplex.data(), ref_simplex.size()),
                std::span<const T>(phys_simplex.data(), phys_simplex.size()),
                parent_tdim,
                gdim,
                order);
        }
    }
    else if (is_simplex(cell_type))
    {
        append_simplex_quadrature(
            rules,
            cell_type,
            ref_vertices,
            physical_vertices,
            parent_tdim,
            gdim,
            order);
    }
    else
    {
        const auto ref_rule = quadrature::get_reference_rule<T>(cell_type, order);
        const T measure = cell::affine_volume_factor<T>(
            cell_type, phys_use.data(), gdim);
        rules._points.insert(
            rules._points.end(), ref_rule._points.begin(), ref_rule._points.end());
        for (int q = 0; q < ref_rule._num_points; ++q)
            rules._weights.push_back(ref_rule._weights[q] * measure);
    }

    rules._parent_map.push_back(parent_cell_id);
    rules._offset.push_back(static_cast<int32_t>(rules._weights.size()));
}

template <std::floating_point T>
void append_lagrange_child_simplex(CurvedVTUGrid<T>& out,
                                   const LocalCurvedSimplexMap<T>& map,
                                   int parent_cell_id,
                                   int curving_status,
                                   std::span<const T> child_vertices,
                                   int subdivision_depth,
                                   int max_subdivision_depth)
{
    const int order = std::max(map.geometry_order, 1);
    const bool curving_failed =
        curving_status == static_cast<int>(curving::CurvingStatus::failed);
    const bool map_valid = curved_child_map_valid<T>(map, child_vertices);
    const bool valid = !curving_failed && map_valid;
    if (!map_valid && !curving_failed && subdivision_depth < max_subdivision_depth)
    {
        const auto children = subdivide_child_simplex<T>(map.simplex_type, child_vertices);
        for (const auto& child : children)
        {
            append_lagrange_child_simplex<T>(
                out,
                map,
                parent_cell_id,
                curving_status,
                std::span<const T>(child.data(), child.size()),
                subdivision_depth + 1,
                max_subdivision_depth);
        }
        return;
    }

    const auto local_nodes = lagrange_cell_nodes<T>(map.simplex_type, order);
    const int local_dofs = static_cast<int>(local_nodes.size());
    const int point_base = static_cast<int>(out.points.size()) / out.gdim;

    for (const auto& xi : local_nodes)
    {
        const auto xi_parent = map_child_xi_to_parent_xi<T>(
            map.simplex_type,
            child_vertices,
            std::span<const T>(xi.data(), xi.size()));
        std::vector<T> x;
        if (valid)
        {
            x = curved_map_physical_point<T>(
                map, std::span<const T>(xi_parent.data(), xi_parent.size()));
        }
        else
        {
            const auto ref = straight_ref_from_cell_point<T>(
                map, std::span<const T>(xi_parent.data(), xi_parent.size()));
            x = push_parent_ref_to_physical<T>(map, std::span<const T>(ref.data(), ref.size()));
        }
        out.points.insert(out.points.end(), x.begin(), x.end());
    }

    std::vector<int> perm;
    int vtk_type = 0;
    if (map.simplex_type == cell::type::interval)
    {
        vtk_type = vtk_lagrange_curve;
        perm.resize(static_cast<std::size_t>(local_dofs));
        std::iota(perm.begin(), perm.end(), 0);
    }
    else if (map.simplex_type == cell::type::triangle)
    {
        vtk_type = vtk_lagrange_triangle;
        perm = io::basix_to_vtk_lagrange_permutation(
            cell::type::triangle, local_dofs, order);
    }
    else if (map.simplex_type == cell::type::quadrilateral)
    {
        vtk_type = vtk_lagrange_quadrilateral;
        perm.resize(static_cast<std::size_t>(local_dofs));
        std::iota(perm.begin(), perm.end(), 0);
    }
    else if (map.simplex_type == cell::type::tetrahedron)
    {
        vtk_type = vtk_lagrange_tetrahedron;
        perm = io::basix_to_vtk_lagrange_permutation(
            cell::type::tetrahedron, local_dofs, order);
    }
    else if (map.simplex_type == cell::type::prism)
    {
        vtk_type = vtk_lagrange_wedge;
        perm.resize(static_cast<std::size_t>(local_dofs));
        std::iota(perm.begin(), perm.end(), 0);
    }
    else
    {
        throw std::runtime_error("append_lagrange_simplex: unsupported simplex type");
    }

    for (const int p : perm)
        out.connectivity.push_back(point_base + p);
    out.offsets.push_back(static_cast<int>(out.connectivity.size()));
    out.vtk_types.push_back(vtk_type);
    out.parent_map.push_back(static_cast<std::int32_t>(parent_cell_id));
    out.curved_valid.push_back(valid ? 1 : 0);
    out.subdivision_depth.push_back(static_cast<std::int32_t>(subdivision_depth));
    out.curving_status.push_back(static_cast<std::int32_t>(curving_status));
}

template <std::floating_point T>
void append_lagrange_simplex(CurvedVTUGrid<T>& out,
                             const LocalCurvedSimplexMap<T>& map,
                             int parent_cell_id,
                             int curving_status,
                             int max_subdivision_depth)
{
    const auto child = canonical_simplex_vertices<T>(map.simplex_type);
    append_lagrange_child_simplex<T>(
        out,
        map,
        parent_cell_id,
        curving_status,
        std::span<const T>(child.data(), child.size()),
        0,
        max_subdivision_depth);
}

template <std::floating_point T>
void append_linear_entity_to_curved_grid(CurvedVTUGrid<T>& out,
                                         const AdaptCell<T>& adapt_cell,
                                         const SelectedEntity& entity,
                                         std::span<const T> parent_vertex_coords,
                                         int parent_cell_id)
{
    const auto ref_coords = entity_reference_coords<T>(
        adapt_cell,
        std::span<const int>(entity.vertices.data(), entity.vertices.size()));
    const std::vector<T> parent_coords(parent_vertex_coords.begin(),
                                       parent_vertex_coords.end());
    const auto phys_coords = cell::push_forward_affine_map<T>(
        adapt_cell.parent_cell_type,
        parent_coords,
        out.gdim,
        std::span<const T>(ref_coords.data(), ref_coords.size()));

    const int point_base = static_cast<int>(out.points.size()) / out.gdim;
    out.points.insert(out.points.end(), phys_coords.begin(), phys_coords.end());

    if (is_simplex(entity.type))
    {
        for (int i = 0; i < static_cast<int>(entity.vertices.size()); ++i)
            out.connectivity.push_back(point_base + i);
    }
    else
    {
        const auto perm = cell::basix_to_vtk_vertex_permutation(entity.type);
        for (const int p : perm)
            out.connectivity.push_back(point_base + p);
    }

    out.offsets.push_back(static_cast<int>(out.connectivity.size()));
    out.vtk_types.push_back(static_cast<int>(cell::map_cell_type_to_vtk(entity.type)));
    out.parent_map.push_back(static_cast<std::int32_t>(parent_cell_id));
    out.curved_valid.push_back(0);
    out.subdivision_depth.push_back(0);
    out.curving_status.push_back(0);
}

template <std::floating_point T>
void append_curved_child_simplex_quadrature(quadrature::QuadratureRules<T>& rules,
                                            const LocalCurvedSimplexMap<T>& map,
                                            int parent_cell_id,
                                            int order,
                                            std::span<const T> child_vertices,
                                            int subdivision_depth,
                                            int max_subdivision_depth)
{
    if (rules._tdim == 0)
        rules._tdim = map.parent_tdim;
    if (rules._offset.empty())
        rules._offset.push_back(0);

    const bool valid = curved_child_map_valid<T>(map, child_vertices);
    if (!valid && subdivision_depth < max_subdivision_depth)
    {
        const auto children = subdivide_child_simplex<T>(map.simplex_type, child_vertices);
        for (const auto& child : children)
        {
            append_curved_child_simplex_quadrature<T>(
                rules,
                map,
                parent_cell_id,
                order,
                std::span<const T>(child.data(), child.size()),
                subdivision_depth + 1,
                max_subdivision_depth);
        }
        return;
    }

    const auto ref_rule = quadrature::get_reference_rule<T>(map.simplex_type, order);
    for (int q = 0; q < ref_rule._num_points; ++q)
    {
        std::span<const T> xi_child(
            ref_rule._points.data() + static_cast<std::size_t>(q * ref_rule._tdim),
            static_cast<std::size_t>(ref_rule._tdim));
        const auto xi = map_child_xi_to_parent_xi<T>(
            map.simplex_type,
            child_vertices,
            xi_child);
        std::vector<T> ref;
        if (valid)
            ref = curved_map_ref_point<T>(map, std::span<const T>(xi.data(), xi.size()));
        else
        {
            ref = straight_ref_from_cell_point<T>(map, std::span<const T>(xi.data(), xi.size()));
        }
        rules._points.insert(rules._points.end(), ref.begin(), ref.end());

        const T child_measure = finite_difference_child_measure<T>(
            map.simplex_type, child_vertices, xi_child);
        T measure = T(0);
        if (valid)
        {
            measure = std::abs(finite_difference_measure<T>(
                map, std::span<const T>(xi.data(), xi.size()))) * child_measure;
        }
        else
        {
            measure = finite_difference_straight_measure<T>(
                map, std::span<const T>(xi.data(), xi.size())) * child_measure;
        }
        rules._weights.push_back(ref_rule._weights[static_cast<std::size_t>(q)] * measure);
    }
    rules._parent_map.push_back(static_cast<std::int32_t>(parent_cell_id));
    rules._offset.push_back(static_cast<int32_t>(rules._weights.size()));
}

template <std::floating_point T>
void append_curved_simplex_quadrature(quadrature::QuadratureRules<T>& rules,
                                      const LocalCurvedSimplexMap<T>& map,
                                      int parent_cell_id,
                                      int order,
                                      int max_subdivision_depth)
{
    const auto child = canonical_simplex_vertices<T>(map.simplex_type);
    append_curved_child_simplex_quadrature<T>(
        rules,
        map,
        parent_cell_id,
        order,
        std::span<const T>(child.data(), child.size()),
        0,
        max_subdivision_depth);
}

template <std::floating_point T, std::integral I>
LocalCurvedSimplexMap<T> make_local_map(const HOMeshPart<T, I>& part,
                                        const AdaptCell<T>& adapt_cell,
                                        int cut_cell_id,
                                        std::span<const int> simplex_vertices,
                                        cell::type simplex_type,
                                        std::span<const T> parent_vertex_coords,
                                        int geometry_order,
                                        curving::NodeFamily node_family,
                                        const curving::CurvingOptions<T>& options,
                                        int required_zero_entity_index = -1)
{
    const auto& mesh = *part.bg->mesh;
    LocalCurvedSimplexMap<T> map;
    map.simplex_type = simplex_type;
    map.dim = cell::get_tdim(simplex_type);
    map.parent_tdim = mesh.tdim;
    map.gdim = mesh.gdim;
    map.parent_cell_type = adapt_cell.parent_cell_type;
    map.geometry_order = geometry_order;
    map.node_family = node_family;
    map.vertices.assign(simplex_vertices.begin(), simplex_vertices.end());
    map.parent_physical_coords.assign(parent_vertex_coords.begin(), parent_vertex_coords.end());
    map.ref_vertex_coords = entity_reference_coords<T>(adapt_cell, simplex_vertices);

    attach_curved_boundaries<T, I>(
        map, part, adapt_cell, cut_cell_id, options, required_zero_entity_index);
    return map;
}

template <std::floating_point T, std::integral I>
void append_curved_entity(CurvedVTUGrid<T>& out,
                          const HOMeshPart<T, I>& part,
                          const AdaptCell<T>& adapt_cell,
                          int cut_cell_id,
                          const SelectedEntity& entity,
                          std::span<const T> parent_vertex_coords,
                          int geometry_order,
                          curving::NodeFamily node_family,
                          const curving::CurvingOptions<T>& options)
{
    const int entity_dim = cell::get_tdim(entity.type);
    if (entity_dim == 0)
        return;

    std::vector<std::vector<int>> simplices;
    cell::type simplex_type = entity.type;
    if (is_simplex(entity.type))
    {
        simplices.push_back(entity.vertices);
    }
    else if (supports_curved_lagrange_cell(entity.type))
    {
        simplices.push_back(entity.vertices);
    }
    else
    {
        append_linear_entity_to_curved_grid<T>(
            out,
            adapt_cell,
            entity,
            parent_vertex_coords,
            static_cast<int>(
                part.cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_cell_id)]));
        return;
    }

    for (const auto& simplex_vertices : simplices)
    {
        auto map = make_local_map<T, I>(
            part,
            adapt_cell,
            cut_cell_id,
            std::span<const int>(simplex_vertices.data(), simplex_vertices.size()),
            simplex_type,
            parent_vertex_coords,
            geometry_order,
            node_family,
            options,
            entity.zero_entity_index);
        int status = 0;
        if (entity.zero_entity_index >= 0)
        {
            if (graph_checks_allow_zero_entity_curving<T, I>(
                    part, cut_cell_id, entity.zero_entity_index))
            {
                const auto& state = curving::ensure_curved<T, I>(
                    part.cut_cells->curving,
                    std::span<const I>(part.cut_cells->parent_cell_ids),
                    std::span<const AdaptCell<T>>(part.cut_cells->adapt_cells),
                    std::span<const LevelSetCell<T, I>>(part.cut_cells->level_set_cells),
                    std::span<const int>(part.cut_cells->ls_offsets),
                    cut_cell_id,
                    entity.zero_entity_index,
                    options);
                status = static_cast<int>(state.status);
            }
            else
            {
                status = static_cast<int>(curving::CurvingStatus::failed);
            }
        }
        else if (!map.curved_faces.empty() || !map.curved_edges.empty())
        {
            status = static_cast<int>(curving::CurvingStatus::curved);
        }
        append_lagrange_simplex<T>(
            out, map,
            static_cast<int>(part.cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_cell_id)]),
            status,
            options.max_subdivision_depth);
    }
}

template <std::floating_point T, std::integral I>
void append_curved_entity_quadrature(quadrature::QuadratureRules<T>& rules,
                                     const HOMeshPart<T, I>& part,
                                     const AdaptCell<T>& adapt_cell,
                                     int cut_cell_id,
                                     const SelectedEntity& entity,
                                     std::span<const T> parent_vertex_coords,
                                     int geometry_order,
                                     curving::NodeFamily node_family,
                                     const curving::CurvingOptions<T>& options,
                                     int quadrature_order)
{
    const int entity_dim = cell::get_tdim(entity.type);
    if (entity_dim == 0)
        return;

    std::vector<std::vector<int>> simplices;
    cell::type simplex_type = entity.type;
    if (is_simplex(entity.type))
    {
        simplices.push_back(entity.vertices);
    }
    else if (supports_curved_lagrange_cell(entity.type))
    {
        simplices.push_back(entity.vertices);
    }
    else
    {
        return;
    }

    const int parent_cell_id =
        static_cast<int>(part.cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_cell_id)]);
    for (const auto& simplex_vertices : simplices)
    {
        auto map = make_local_map<T, I>(
            part,
            adapt_cell,
            cut_cell_id,
            std::span<const int>(simplex_vertices.data(), simplex_vertices.size()),
            simplex_type,
            parent_vertex_coords,
            geometry_order,
            node_family,
            options,
            entity.zero_entity_index);
        append_curved_simplex_quadrature<T>(
            rules, map, parent_cell_id, quadrature_order,
            options.max_subdivision_depth);
    }
}

template <std::floating_point T, std::integral I>
void append_curved_cut_quadrature(quadrature::QuadratureRules<T>& rules,
                                  const HOMeshPart<T, I>& part,
                                  int quadrature_order,
                                  int geometry_order,
                                  curving::NodeFamily node_family,
                                  const curving::CurvingOptions<T>& options)
{
    const auto& mesh = *part.bg->mesh;
    for (std::int32_t cut_id : part.cut_cell_ids)
    {
        const auto& adapt_cell = part.cut_cells->adapt_cells[static_cast<std::size_t>(cut_id)];
        const auto entities = selected_entities(part, adapt_cell, cut_id);
        if (entities.empty())
            continue;

        const I parent_cell_id = part.cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_id)];
        const auto cut_active_mask =
            part.cut_cells->active_level_set_mask[static_cast<std::size_t>(cut_id)];
        const auto parent_vertex_coords = parent_cell_vertex_coords_vtk(mesh, parent_cell_id);
        for (const auto& entity : entities)
        {
            append_curved_entity_quadrature<T, I>(
                rules,
                part,
                adapt_cell,
                cut_id,
                entity,
                std::span<const T>(parent_vertex_coords.data(), parent_vertex_coords.size()),
                geometry_order,
                node_family,
                options,
                quadrature_order);
        }
    }
}

template <std::floating_point T, std::integral I>
void append_cut_entities(mesh::CutMesh<T>& out,
                         quadrature::QuadratureRules<T>* rules,
                         const HOMeshPart<T, I>& part,
                         int quadrature_order)
{
    const auto& mesh = *part.bg->mesh;
    for (std::int32_t cut_id : part.cut_cell_ids)
    {
        const auto& adapt_cell = part.cut_cells->adapt_cells[static_cast<std::size_t>(cut_id)];
        const auto entities = selected_entities(part, adapt_cell, cut_id);
        if (entities.empty())
            continue;

        const I parent_cell_id = part.cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_id)];
        const auto parent_vertex_coords = parent_cell_vertex_coords_vtk(mesh, parent_cell_id);

        for (const auto& entity : entities)
        {
            std::vector<int> root_vertex_flags;
            if ((entity.type == cell::type::prism
                 && adapt_cell.parent_cell_type == cell::type::tetrahedron)
                || (entity.type == cell::type::quadrilateral
                    && adapt_cell.parent_cell_type == cell::type::triangle))
            {
                root_vertex_flags.resize(entity.vertices.size(), 0);
                for (std::size_t j = 0; j < entity.vertices.size(); ++j)
                {
                    const auto zm = adapt_cell.zero_mask_per_vertex[
                        static_cast<std::size_t>(entity.vertices[j])];
                    bool is_root = false;
                    if (part.expr.zero_required != 0)
                    {
                        is_root = (zm & part.expr.zero_required) == part.expr.zero_required;
                    }
                    else
                    {
                        is_root = zm != 0;
                    }
                    root_vertex_flags[j] = is_root ? 1 : 0;
                }
            }

            const auto ref_coords = entity_reference_coords(
                adapt_cell, std::span<const int>(entity.vertices.data(), entity.vertices.size()));
            const auto phys_coords = cell::push_forward_affine_map<T>(
                adapt_cell.parent_cell_type,
                parent_vertex_coords,
                mesh.gdim,
                std::span<const T>(ref_coords.data(), ref_coords.size()));

            append_mesh_entity(
                out,
                std::span<const T>(phys_coords.data(), phys_coords.size()),
                mesh.gdim,
                entity.type,
                static_cast<int>(parent_cell_id),
                /*triangulate=*/false,
                /*input_is_basix=*/true,
                adapt_cell.parent_cell_type,
                std::span<const int>(root_vertex_flags.data(), root_vertex_flags.size()));

            if (rules != nullptr)
            {
                append_entity_quadrature(
                    *rules,
                    entity.type,
                    std::span<const T>(ref_coords.data(), ref_coords.size()),
                    std::span<const T>(phys_coords.data(), phys_coords.size()),
                    mesh.tdim,
                    mesh.gdim,
                    static_cast<int>(parent_cell_id),
                    quadrature_order,
                    /*triangulate=*/false,
                    /*input_is_basix=*/true,
                    adapt_cell.parent_cell_type,
                    std::span<const int>(root_vertex_flags.data(), root_vertex_flags.size()));
            }
        }
    }
}

template <std::floating_point T, std::integral I>
void append_uncut_volume_cells(mesh::CutMesh<T>& out,
                               quadrature::QuadratureRules<T>* rules,
                               const HOMeshPart<T, I>& part,
                               int quadrature_order)
{
    const auto& mesh = *part.bg->mesh;
    if (part.dim != mesh.tdim)
        return;

    for (I cell_id : part.uncut_cell_ids)
    {
        const auto ctype = mesh.cell_type(cell_id);
        const auto phys_coords = parent_cell_vertex_coords_vtk(mesh, cell_id);
        append_mesh_entity(
            out,
            std::span<const T>(phys_coords.data(), phys_coords.size()),
            mesh.gdim,
            ctype,
            static_cast<int>(cell_id),
            /*triangulate=*/false,
            /*input_is_basix=*/false,
            cell::type::point,
            std::span<const int>());

        if (rules != nullptr)
        {
            const auto ref_coords = cell::canonical_vertices<T>(ctype);
            append_entity_quadrature(
                *rules,
                ctype,
                std::span<const T>(ref_coords.data(), ref_coords.size()),
                std::span<const T>(phys_coords.data(), phys_coords.size()),
                mesh.tdim,
                mesh.gdim,
                static_cast<int>(cell_id),
                quadrature_order,
                /*triangulate=*/false,
                /*input_is_basix=*/false,
                cell::type::point,
                std::span<const int>());
        }
    }
}

template <std::floating_point T, std::integral I>
void append_uncut_curved_cells(CurvedVTUGrid<T>& out,
                               const HOMeshPart<T, I>& part,
                               int geometry_order,
                               curving::NodeFamily node_family)
{
    const auto& mesh = *part.bg->mesh;
    if (part.dim != mesh.tdim)
        return;

    for (I cell_id : part.uncut_cell_ids)
    {
        const auto ctype = mesh.cell_type(cell_id);
        if (!supports_curved_lagrange_cell(ctype))
            continue;

        const auto parent_coords = parent_cell_vertex_coords_vtk(mesh, cell_id);
        const auto ref_vertices = cell::reference_vertices<T>(ctype);
        const int nv = cell::get_num_vertices(ctype);
        LocalCurvedSimplexMap<T> map;
        map.simplex_type = ctype;
        map.dim = cell::get_tdim(ctype);
        map.parent_tdim = mesh.tdim;
        map.gdim = mesh.gdim;
        map.parent_cell_type = ctype;
        map.geometry_order = geometry_order;
        map.node_family = node_family;
        map.vertices.resize(static_cast<std::size_t>(nv));
        std::iota(map.vertices.begin(), map.vertices.end(), 0);
        map.ref_vertex_coords = ref_vertices;
        map.parent_physical_coords = parent_coords;
        append_lagrange_simplex<T>(out, map, static_cast<int>(cell_id), 0, 0);
    }
}

} // namespace

template <std::floating_point T, std::integral I>
mesh::CutMesh<T> visualization_mesh(const HOMeshPart<T, I>& part,
                                    bool include_uncut_cells,
                                    int geometry_order,
                                    curving::NodeFamily node_family,
                                    curving::CurvingDirectionMode direction_mode)
{
    if (!part.cut_cells || !part.bg || !part.bg->mesh)
        throw std::runtime_error("HOMeshPart is not attached to cut-cell storage");

    mesh::CutMesh<T> out;
    out._gdim = part.bg->mesh->gdim;
    out._tdim = part.dim;
    out._offset.push_back(0);

    append_cut_entities(
        out,
        static_cast<quadrature::QuadratureRules<T>*>(nullptr),
        part,
        /*quadrature_order=*/0);
    if (include_uncut_cells)
    {
        append_uncut_volume_cells(
            out,
            static_cast<quadrature::QuadratureRules<T>*>(nullptr),
            part,
            /*quadrature_order=*/0);
    }

    return out;
}

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> quadrature_rules(const HOMeshPart<T, I>& part,
                                                int order,
                                                bool include_uncut_cells,
                                                int geometry_order,
                                                curving::NodeFamily node_family,
                                                curving::CurvingDirectionMode direction_mode)
{
    if (!part.cut_cells || !part.bg || !part.bg->mesh)
        throw std::runtime_error("HOMeshPart is not attached to cut-cell storage");
 
    curving::CurvingOptions<T> options;
    if (geometry_order > 1)
    {
        options.geometry_order = geometry_order;
        options.node_family = construction_node_family(geometry_order, node_family);
        options.direction_mode = direction_mode;
    }

    mesh::CutMesh<T> unused_mesh;
    quadrature::QuadratureRules<T> rules;
    rules._offset.push_back(0);

    if (geometry_order > 1)
    {
        append_curved_cut_quadrature<T, I>(
            rules, part, order, geometry_order, options.node_family, options);
    }
    else
    {
        append_cut_entities(unused_mesh, &rules, part, order);
    }
    if (include_uncut_cells)
        append_uncut_volume_cells(unused_mesh, &rules, part, order);

    return rules;
}

template <std::floating_point T, std::integral I>
CurvedVTUGrid<T> curved_lagrange_grid(const HOMeshPart<T, I>& part,
                                      bool include_uncut_cells,
                                      int geometry_order,
                                      curving::NodeFamily node_family,
                                      curving::CurvingDirectionMode direction_mode)
{
    if (!part.cut_cells || !part.bg || !part.bg->mesh)
        throw std::runtime_error("HOMeshPart is not attached to cut-cell storage");
    if (geometry_order <= 1)
        throw std::runtime_error("curved_lagrange_grid requires geometry_order > 1");

    curving::CurvingOptions<T> options;
    options.geometry_order = geometry_order;
    options.node_family = construction_node_family(geometry_order, node_family);
    options.direction_mode = direction_mode;

    CurvedVTUGrid<T> out;
    out.gdim = part.bg->mesh->gdim;
    out.tdim = part.dim;
    out.geometry_order = geometry_order;
    out.offsets.push_back(0);

    const auto& mesh = *part.bg->mesh;
    for (std::int32_t cut_id : part.cut_cell_ids)
    {
        const auto& adapt_cell = part.cut_cells->adapt_cells[static_cast<std::size_t>(cut_id)];
        const auto entities = selected_entities(part, adapt_cell, cut_id);
        if (entities.empty())
            continue;

        const I parent_cell_id = part.cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_id)];
        const auto parent_vertex_coords = parent_cell_vertex_coords_vtk(mesh, parent_cell_id);
        for (const auto& entity : entities)
        {
            append_curved_entity<T, I>(
                out,
                part,
                adapt_cell,
                cut_id,
                entity,
                std::span<const T>(parent_vertex_coords.data(), parent_vertex_coords.size()),
                geometry_order,
                options.node_family,
                options);
        }
    }

    if (include_uncut_cells)
        append_uncut_curved_cells<T, I>(out, part, geometry_order, node_family);

    return out;
}

template mesh::CutMesh<double> visualization_mesh(
    const HOMeshPart<double, int>&, bool, int, curving::NodeFamily,
    curving::CurvingDirectionMode);
template mesh::CutMesh<float> visualization_mesh(
    const HOMeshPart<float, int>&, bool, int, curving::NodeFamily,
    curving::CurvingDirectionMode);
template mesh::CutMesh<double> visualization_mesh(
    const HOMeshPart<double, long>&, bool, int, curving::NodeFamily,
    curving::CurvingDirectionMode);
template mesh::CutMesh<float> visualization_mesh(
    const HOMeshPart<float, long>&, bool, int, curving::NodeFamily,
    curving::CurvingDirectionMode);

template std::vector<SelectedZeroEntityInfo> selected_zero_entity_infos(
    const HOMeshPart<double, int>&);
template std::vector<SelectedZeroEntityInfo> selected_zero_entity_infos(
    const HOMeshPart<float, int>&);
template std::vector<SelectedZeroEntityInfo> selected_zero_entity_infos(
    const HOMeshPart<double, long>&);
template std::vector<SelectedZeroEntityInfo> selected_zero_entity_infos(
    const HOMeshPart<float, long>&);

template quadrature::QuadratureRules<double> quadrature_rules(
    const HOMeshPart<double, int>&, int, bool, int, curving::NodeFamily,
    curving::CurvingDirectionMode);
template quadrature::QuadratureRules<float> quadrature_rules(
    const HOMeshPart<float, int>&, int, bool, int, curving::NodeFamily,
    curving::CurvingDirectionMode);
template quadrature::QuadratureRules<double> quadrature_rules(
    const HOMeshPart<double, long>&, int, bool, int, curving::NodeFamily,
    curving::CurvingDirectionMode);
template quadrature::QuadratureRules<float> quadrature_rules(
    const HOMeshPart<float, long>&, int, bool, int, curving::NodeFamily,
    curving::CurvingDirectionMode);

template CurvedVTUGrid<double> curved_lagrange_grid(
    const HOMeshPart<double, int>&, bool, int, curving::NodeFamily,
    curving::CurvingDirectionMode);
template CurvedVTUGrid<float> curved_lagrange_grid(
    const HOMeshPart<float, int>&, bool, int, curving::NodeFamily,
    curving::CurvingDirectionMode);
template CurvedVTUGrid<double> curved_lagrange_grid(
    const HOMeshPart<double, long>&, bool, int, curving::NodeFamily,
    curving::CurvingDirectionMode);
template CurvedVTUGrid<float> curved_lagrange_grid(
    const HOMeshPart<float, long>&, bool, int, curving::NodeFamily,
    curving::CurvingDirectionMode);

} // namespace cutcells::output
