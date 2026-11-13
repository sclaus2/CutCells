// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "level_set_cell.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "bernstein.h"
#include "cell_topology.h"
#include "mapping.h"
#include "reference_cell.h"

namespace cutcells
{
namespace
{

template <std::floating_point T>
struct RefPointKey
{
    cell::type ctype = cell::type::point;
    int degree = 0;

    bool operator==(const RefPointKey&) const = default;
};

template <std::floating_point T>
struct RefPointKeyHash
{
    std::size_t operator()(const RefPointKey<T>& key) const noexcept
    {
        std::size_t seed = std::hash<int>{}(static_cast<int>(key.ctype));
        seed ^= std::hash<int>{}(key.degree) + 0x9e3779b97f4a7c15ULL
              + (seed << 6) + (seed >> 2);
        return seed;
    }
};

template <typename T>
struct FaceDef
{
    cell::type face_type = cell::type::point;
    int nverts = 0;
    std::array<int, 4> verts = {-1, -1, -1, -1};
};

inline std::span<const FaceDef<int>> basix_faces(cell::type ctype)
{
    static constexpr std::array<FaceDef<int>, 4> tetrahedron = {{
        {cell::type::triangle, 3, {1, 2, 3, -1}},
        {cell::type::triangle, 3, {0, 2, 3, -1}},
        {cell::type::triangle, 3, {0, 1, 3, -1}},
        {cell::type::triangle, 3, {0, 1, 2, -1}},
    }};
    static constexpr std::array<FaceDef<int>, 6> hexahedron = {{
        {cell::type::quadrilateral, 4, {0, 1, 2, 3}},
        {cell::type::quadrilateral, 4, {0, 1, 4, 5}},
        {cell::type::quadrilateral, 4, {0, 2, 4, 6}},
        {cell::type::quadrilateral, 4, {1, 3, 5, 7}},
        {cell::type::quadrilateral, 4, {2, 3, 6, 7}},
        {cell::type::quadrilateral, 4, {4, 5, 6, 7}},
    }};
    static constexpr std::array<FaceDef<int>, 5> prism = {{
        {cell::type::triangle, 3, {0, 1, 2, -1}},
        {cell::type::triangle, 3, {3, 4, 5, -1}},
        {cell::type::quadrilateral, 4, {0, 1, 3, 4}},
        {cell::type::quadrilateral, 4, {0, 2, 3, 5}},
        {cell::type::quadrilateral, 4, {1, 2, 4, 5}},
    }};
    static constexpr std::array<FaceDef<int>, 5> pyramid = {{
        {cell::type::quadrilateral, 4, {0, 1, 2, 3}},
        {cell::type::triangle, 3, {0, 1, 4, -1}},
        {cell::type::triangle, 3, {0, 2, 4, -1}},
        {cell::type::triangle, 3, {1, 3, 4, -1}},
        {cell::type::triangle, 3, {2, 3, 4, -1}},
    }};

    switch (ctype)
    {
        case cell::type::tetrahedron: return std::span(tetrahedron);
        case cell::type::hexahedron: return std::span(hexahedron);
        case cell::type::prism: return std::span(prism);
        case cell::type::pyramid: return std::span(pyramid);
        default: return {};
    }
}

template <std::floating_point T>
void append_point(std::vector<T>& out, std::span<const T> x)
{
    out.insert(out.end(), x.begin(), x.end());
}

template <std::floating_point T>
void append_affine_edge_point(std::vector<T>& out,
                              std::span<const T> x0,
                              std::span<const T> x1,
                              T t)
{
    for (std::size_t d = 0; d < x0.size(); ++d)
        out.push_back((T(1) - t) * x0[d] + t * x1[d]);
}

template <std::floating_point T>
void append_triangle_face_point(std::vector<T>& out,
                                std::span<const T> x0,
                                std::span<const T> x1,
                                std::span<const T> x2,
                                T u, T v)
{
    const T w0 = T(1) - u - v;
    for (std::size_t d = 0; d < x0.size(); ++d)
        out.push_back(w0 * x0[d] + u * x1[d] + v * x2[d]);
}

template <std::floating_point T>
void append_quad_face_point(std::vector<T>& out,
                            std::span<const T> x0,
                            std::span<const T> x1,
                            std::span<const T> x2,
                            std::span<const T> x3,
                            T u, T v)
{
    const T w0 = (T(1) - u) * (T(1) - v);
    const T w1 = u * (T(1) - v);
    const T w2 = (T(1) - u) * v;
    const T w3 = u * v;
    for (std::size_t d = 0; d < x0.size(); ++d)
        out.push_back(w0 * x0[d] + w1 * x1[d] + w2 * x2[d] + w3 * x3[d]);
}

template <std::floating_point T>
void append_tetra_cell_point(std::vector<T>& out,
                             std::span<const T> x0,
                             std::span<const T> x1,
                             std::span<const T> x2,
                             std::span<const T> x3,
                             T u, T v, T w)
{
    const T w0 = T(1) - u - v - w;
    for (std::size_t d = 0; d < x0.size(); ++d)
        out.push_back(w0 * x0[d] + u * x1[d] + v * x2[d] + w * x3[d]);
}

template <std::floating_point T>
void append_hex_cell_point(std::vector<T>& out,
                           const std::vector<T>& verts, int tdim,
                           T u, T v, T w)
{
    for (int d = 0; d < tdim; ++d)
    {
        out.push_back(
            (T(1) - u) * (T(1) - v) * (T(1) - w) * verts[static_cast<std::size_t>(d)] +
            u * (T(1) - v) * (T(1) - w) * verts[static_cast<std::size_t>(tdim + d)] +
            (T(1) - u) * v * (T(1) - w) * verts[static_cast<std::size_t>(2 * tdim + d)] +
            u * v * (T(1) - w) * verts[static_cast<std::size_t>(3 * tdim + d)] +
            (T(1) - u) * (T(1) - v) * w * verts[static_cast<std::size_t>(4 * tdim + d)] +
            u * (T(1) - v) * w * verts[static_cast<std::size_t>(5 * tdim + d)] +
            (T(1) - u) * v * w * verts[static_cast<std::size_t>(6 * tdim + d)] +
            u * v * w * verts[static_cast<std::size_t>(7 * tdim + d)]);
    }
}

template <std::floating_point T>
std::vector<T> build_reference_lagrange_points(cell::type ctype, int degree)
{
    const int tdim = cell::get_tdim(ctype);
    const int nv = cell::get_num_vertices(ctype);
    const std::vector<T> vertices = cell::reference_vertices<T>(ctype);
    std::vector<T> ref_points;
    ref_points.reserve(static_cast<std::size_t>(
        bernstein::num_polynomials(ctype, degree) * tdim));

    for (int v = 0; v < nv; ++v)
    {
        append_point(
            ref_points,
            std::span<const T>(
                vertices.data() + static_cast<std::size_t>(v * tdim),
                static_cast<std::size_t>(tdim)));
    }

    if (degree == 1)
        return ref_points;

    const auto edges = cell::edges(ctype);
    for (const auto& edge : edges)
    {
        const T* x0 = vertices.data() + static_cast<std::size_t>(edge[0] * tdim);
        const T* x1 = vertices.data() + static_cast<std::size_t>(edge[1] * tdim);
        for (int k = 1; k < degree; ++k)
        {
            const T t = static_cast<T>(k) / static_cast<T>(degree);
            append_affine_edge_point(
                ref_points,
                std::span<const T>(x0, static_cast<std::size_t>(tdim)),
                std::span<const T>(x1, static_cast<std::size_t>(tdim)),
                t);
        }
    }

    if (tdim == 3)
    {
        const auto faces = basix_faces(ctype);
        for (const auto& face : faces)
        {
            if (face.face_type == cell::type::triangle)
            {
                if (degree <= 2)
                    continue;

                const T* x0 = vertices.data() + static_cast<std::size_t>(face.verts[0] * tdim);
                const T* x1 = vertices.data() + static_cast<std::size_t>(face.verts[1] * tdim);
                const T* x2 = vertices.data() + static_cast<std::size_t>(face.verts[2] * tdim);
                for (int j = 1; j <= degree - 2; ++j)
                {
                    for (int i = 1; i <= degree - j - 1; ++i)
                    {
                        const T u = static_cast<T>(i) / static_cast<T>(degree);
                        const T v = static_cast<T>(j) / static_cast<T>(degree);
                        append_triangle_face_point(
                            ref_points,
                            std::span<const T>(x0, static_cast<std::size_t>(tdim)),
                            std::span<const T>(x1, static_cast<std::size_t>(tdim)),
                            std::span<const T>(x2, static_cast<std::size_t>(tdim)),
                            u, v);
                    }
                }
            }
            else if (face.face_type == cell::type::quadrilateral)
            {
                const T* x0 = vertices.data() + static_cast<std::size_t>(face.verts[0] * tdim);
                const T* x1 = vertices.data() + static_cast<std::size_t>(face.verts[1] * tdim);
                const T* x2 = vertices.data() + static_cast<std::size_t>(face.verts[2] * tdim);
                const T* x3 = vertices.data() + static_cast<std::size_t>(face.verts[3] * tdim);
                for (int j = 1; j < degree; ++j)
                {
                    for (int i = 1; i < degree; ++i)
                    {
                        const T u = static_cast<T>(i) / static_cast<T>(degree);
                        const T v = static_cast<T>(j) / static_cast<T>(degree);
                        append_quad_face_point(
                            ref_points,
                            std::span<const T>(x0, static_cast<std::size_t>(tdim)),
                            std::span<const T>(x1, static_cast<std::size_t>(tdim)),
                            std::span<const T>(x2, static_cast<std::size_t>(tdim)),
                            std::span<const T>(x3, static_cast<std::size_t>(tdim)),
                            u, v);
                    }
                }
            }
        }
    }

    switch (ctype)
    {
        case cell::type::triangle:
            if (degree > 2)
            {
                const T* x0 = vertices.data();
                const T* x1 = vertices.data() + tdim;
                const T* x2 = vertices.data() + 2 * tdim;
                for (int j = 1; j <= degree - 2; ++j)
                {
                    for (int i = 1; i <= degree - j - 1; ++i)
                    {
                        const T u = static_cast<T>(i) / static_cast<T>(degree);
                        const T v = static_cast<T>(j) / static_cast<T>(degree);
                        append_triangle_face_point(
                            ref_points,
                            std::span<const T>(x0, static_cast<std::size_t>(tdim)),
                            std::span<const T>(x1, static_cast<std::size_t>(tdim)),
                            std::span<const T>(x2, static_cast<std::size_t>(tdim)),
                            u, v);
                    }
                }
            }
            break;
        case cell::type::quadrilateral:
            {
                const T* x0 = vertices.data();
                const T* x1 = vertices.data() + tdim;
                const T* x2 = vertices.data() + 2 * tdim;
                const T* x3 = vertices.data() + 3 * tdim;
                for (int j = 1; j < degree; ++j)
                {
                    for (int i = 1; i < degree; ++i)
                    {
                        const T u = static_cast<T>(i) / static_cast<T>(degree);
                        const T v = static_cast<T>(j) / static_cast<T>(degree);
                        append_quad_face_point(
                            ref_points,
                            std::span<const T>(x0, static_cast<std::size_t>(tdim)),
                            std::span<const T>(x1, static_cast<std::size_t>(tdim)),
                            std::span<const T>(x2, static_cast<std::size_t>(tdim)),
                            std::span<const T>(x3, static_cast<std::size_t>(tdim)),
                            u, v);
                    }
                }
            }
            break;
        case cell::type::tetrahedron:
            if (degree > 3)
            {
                const T* x0 = vertices.data();
                const T* x1 = vertices.data() + tdim;
                const T* x2 = vertices.data() + 2 * tdim;
                const T* x3 = vertices.data() + 3 * tdim;
                for (int k = 1; k <= degree - 3; ++k)
                {
                    for (int j = 1; j <= degree - k - 2; ++j)
                    {
                        for (int i = 1; i <= degree - j - k - 1; ++i)
                        {
                            const T u = static_cast<T>(i) / static_cast<T>(degree);
                            const T v = static_cast<T>(j) / static_cast<T>(degree);
                            const T w = static_cast<T>(k) / static_cast<T>(degree);
                            append_tetra_cell_point(
                                ref_points,
                                std::span<const T>(x0, static_cast<std::size_t>(tdim)),
                                std::span<const T>(x1, static_cast<std::size_t>(tdim)),
                                std::span<const T>(x2, static_cast<std::size_t>(tdim)),
                                std::span<const T>(x3, static_cast<std::size_t>(tdim)),
                                u, v, w);
                        }
                    }
                }
            }
            break;
        case cell::type::hexahedron:
            for (int k = 1; k < degree; ++k)
            {
                for (int j = 1; j < degree; ++j)
                {
                    for (int i = 1; i < degree; ++i)
                    {
                        append_hex_cell_point(
                            ref_points, vertices, tdim,
                            static_cast<T>(i) / static_cast<T>(degree),
                            static_cast<T>(j) / static_cast<T>(degree),
                            static_cast<T>(k) / static_cast<T>(degree));
                    }
                }
            }
            break;
        case cell::type::prism:
            // the triangle's interior points on each interior layer, as
            // create_level_set_mesh_data places them
            for (int k = 1; k < degree; ++k)
            {
                for (int j = 1; j <= degree - 2; ++j)
                {
                    for (int i = 1; i <= degree - j - 1; ++i)
                    {
                        const std::array<T, 3> x = {static_cast<T>(i) / static_cast<T>(degree),
                                                    static_cast<T>(j) / static_cast<T>(degree),
                                                    static_cast<T>(k) / static_cast<T>(degree)};
                        append_point(ref_points, std::span<const T>(x));
                    }
                }
            }
            break;
        case cell::type::pyramid:
            if (degree > 2)
                throw std::invalid_argument("make_cell_level_set: pyramids of degree above 2 are not supported");
            break;
        default:
            break;
    }

    return ref_points;
}

template <std::floating_point T>
const std::vector<T>& cached_reference_lagrange_points(cell::type ctype, int degree)
{
    thread_local std::unordered_map<RefPointKey<T>, std::vector<T>, RefPointKeyHash<T>> cache;

    RefPointKey<T> key;
    key.ctype = ctype;
    key.degree = degree;

    auto it = cache.find(key);
    if (it != cache.end())
        return it->second;

    auto points = build_reference_lagrange_points<T>(ctype, degree);
    auto [inserted_it, inserted] = cache.emplace(key, std::move(points));
    (void)inserted;
    return inserted_it->second;
}

/// Reference coordinates of a cell's dofs: the dof coordinates pulled back
/// through the affine map of the cell (axes from the push-forward, so the
/// convention is CutCells' own; least squares, so entities of lower dimension
/// than the space, such as facets in 3D, work too).
template <std::floating_point T, std::integral I>
std::vector<T> pull_back_dof_coordinates(const LevelSetMeshData<T, I>& md,
                                         cell::type ctype,
                                         const std::vector<T>& vertex_coords,
                                         int gdim, int tdim,
                                         std::span<const I> cell_dofs)
{
    std::vector<T> ref(static_cast<std::size_t>((tdim + 1) * tdim), T(0));
    for (int k = 0; k < tdim; ++k)
        ref[static_cast<std::size_t>((k + 1) * tdim + k)] = T(1);
    const std::vector<T> x = cell::push_forward_affine_map<T>(
        ctype, vertex_coords, gdim, std::span<const T>(ref.data(), ref.size()));
    // J (gdim x tdim): column k = F(e_k) - F(0); G = J^T J
    std::vector<T> J(static_cast<std::size_t>(gdim * tdim));
    for (int k = 0; k < tdim; ++k)
        for (int d = 0; d < gdim; ++d)
            J[static_cast<std::size_t>(d * tdim + k)]
                = x[static_cast<std::size_t>((k + 1) * gdim + d)] - x[static_cast<std::size_t>(d)];
    std::array<T, 9> G{};
    for (int a = 0; a < tdim; ++a)
        for (int b = 0; b < tdim; ++b)
            for (int d = 0; d < gdim; ++d)
                G[static_cast<std::size_t>(a * 3 + b)]
                    += J[static_cast<std::size_t>(d * tdim + a)] * J[static_cast<std::size_t>(d * tdim + b)];
    // inverse of G (tdim <= 3)
    std::array<T, 9> Gi{};
    if (tdim == 1)
        Gi[0] = T(1) / G[0];
    else if (tdim == 2)
    {
        const T det = G[0] * G[4] - G[1] * G[3];
        Gi[0] = G[4] / det; Gi[1] = -G[1] / det; Gi[3] = -G[3] / det; Gi[4] = G[0] / det;
    }
    else
    {
        const T det = G[0] * (G[4] * G[8] - G[5] * G[7]) - G[1] * (G[3] * G[8] - G[5] * G[6])
                      + G[2] * (G[3] * G[7] - G[4] * G[6]);
        Gi[0] = (G[4] * G[8] - G[5] * G[7]) / det; Gi[1] = (G[2] * G[7] - G[1] * G[8]) / det;
        Gi[2] = (G[1] * G[5] - G[2] * G[4]) / det; Gi[3] = (G[5] * G[6] - G[3] * G[8]) / det;
        Gi[4] = (G[0] * G[8] - G[2] * G[6]) / det; Gi[5] = (G[2] * G[3] - G[0] * G[5]) / det;
        Gi[6] = (G[3] * G[7] - G[4] * G[6]) / det; Gi[7] = (G[1] * G[6] - G[0] * G[7]) / det;
        Gi[8] = (G[0] * G[4] - G[1] * G[3]) / det;
    }
    std::vector<T> X(cell_dofs.size() * static_cast<std::size_t>(tdim));
    for (std::size_t j = 0; j < cell_dofs.size(); ++j)
    {
        const T* xj = md.dof_coordinate(cell_dofs[j]);
        std::array<T, 3> rhs{};
        for (int a = 0; a < tdim; ++a)
            for (int d = 0; d < gdim; ++d)
                rhs[static_cast<std::size_t>(a)]
                    += J[static_cast<std::size_t>(d * tdim + a)] * (xj[d] - x[static_cast<std::size_t>(d)]);
        for (int a = 0; a < tdim; ++a)
        {
            T v = T(0);
            for (int b = 0; b < tdim; ++b)
                v += Gi[static_cast<std::size_t>(a * 3 + b)] * rhs[static_cast<std::size_t>(b)];
            X[j * static_cast<std::size_t>(tdim) + static_cast<std::size_t>(a)] = v;
        }
    }
    return X;
}

/// Determine the cell::type for a cell in the level-set mesh data.
/// Uses cell_types array if available, otherwise infers from tdim and
/// number of DOFs per cell (assuming equispaced Lagrange).
template <std::floating_point T, std::integral I>
cell::type infer_ls_cell_type(const LevelSetMeshData<T, I>& md, I cell_id)
{
    if (md.has_cell_types())
        return md.cell_type(cell_id);

    // Infer from tdim and number of DOFs per cell
    const I ndofs = md.cell_num_dofs(cell_id);
    const int n = md.degree;

    switch (md.tdim)
    {
        case 1:
            return cell::type::interval;
        case 2:
        {
            // triangle: C(n+2,2) = (n+1)(n+2)/2
            // quad: (n+1)^2
            int tri_dofs = (n + 1) * (n + 2) / 2;
            if (ndofs == static_cast<I>(tri_dofs))
                return cell::type::triangle;
            else
                return cell::type::quadrilateral;
        }
        case 3:
        {
            // tet: C(n+3,3) = (n+1)(n+2)(n+3)/6
            // hex: (n+1)^3
            int tet_dofs = (n + 1) * (n + 2) * (n + 3) / 6;
            int hex_dofs = (n + 1) * (n + 1) * (n + 1);
            if (ndofs == static_cast<I>(tet_dofs))
                return cell::type::tetrahedron;
            else if (ndofs == static_cast<I>(hex_dofs))
                return cell::type::hexahedron;
            else
            {
                throw std::runtime_error(
                    "make_cell_level_set: cannot infer cell type from "
                    + std::to_string(ndofs) + " DOFs in 3D");
            }
        }
        default:
            throw std::runtime_error(
                "make_cell_level_set: unsupported tdim = "
                + std::to_string(md.tdim));
    }
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// make_cell_level_set
// ---------------------------------------------------------------------------

template <std::floating_point T, std::integral I>
LevelSetCell<T, I>
make_cell_level_set(const LevelSetFunction<T, I>& global_ls,
                    I cell_id)
{
    LevelSetCell<T, I> cell_ls;
    make_cell_level_set(global_ls, cell_id, cell_ls);
    return cell_ls;
}

template <std::floating_point T, std::integral I>
void make_cell_level_set(const LevelSetFunction<T, I>& global_ls,
                         I cell_id, LevelSetCell<T, I>& cell_ls)
{
    // Reset to a default-constructed state, keeping vector capacity.
    cell_ls.level_set_id = 0;
    cell_ls.cell_type = cell::type::point;
    cell_ls.gdim = 0;
    cell_ls.tdim = 0;
    cell_ls.parent_vertex_coords.clear();
    cell_ls.bernstein_coeffs.clear();
    cell_ls.bernstein_order = 0;
    cell_ls.nodal_values.clear();
    cell_ls.nodal_order = 0;
    cell_ls.global_level_set = &global_ls;
    cell_ls.cell_id = cell_id;

    // --- Polynomial path: extract nodal values, convert to Bernstein ---
    if (global_ls.type == LevelSetType::Polynomial
        && global_ls.has_mesh_data()
        && global_ls.has_dof_values())
    {
        const auto& md = global_ls.mesh_data;
        const int degree = md.degree;
        const int gdim = md.gdim;

        // Determine cell type
        cell::type ctype = infer_ls_cell_type(md, cell_id);
        cell_ls.cell_type = ctype;
        cell_ls.gdim = gdim;
        cell_ls.tdim = cell::get_tdim(ctype);

        // Get the cell's DOF indices. Direct layouts return a view; callback
        // layouts use the scratch buffer.
        std::vector<I> cell_dof_scratch;
        auto cell_dofs = md.cell_dofs_span(cell_id, cell_dof_scratch);
        const int ndofs = static_cast<int>(cell_dofs.size());
        const int nv = cell::get_num_vertices(ctype);
        if (ndofs < nv)
        {
            throw std::runtime_error(
                "make_cell_level_set: cell dof map is missing vertex dofs");
        }
        cell_ls.parent_vertex_coords.resize(
            static_cast<std::size_t>(nv * gdim), T(0));
        for (int i = 0; i < nv; ++i)
        {
            const I dof = cell_dofs[static_cast<std::size_t>(i)];
            const T* x = md.dof_coordinate(dof);
            for (int d = 0; d < gdim; ++d)
            {
                cell_ls.parent_vertex_coords[
                    static_cast<std::size_t>(i * gdim + d)] = x[d];
            }
        }

        // Extract nodal values from the global DOF array
        cell_ls.nodal_values.resize(static_cast<std::size_t>(ndofs));
        for (int i = 0; i < ndofs; ++i)
            cell_ls.nodal_values[static_cast<std::size_t>(i)]
                = global_ls.dof_values[static_cast<std::size_t>(cell_dofs[static_cast<std::size_t>(i)])];
        cell_ls.nodal_order = degree;
        std::span<const T> dof_ref;
        std::vector<T> pulled_back;
        if (!md.cell_reference_points.empty())
        {
            dof_ref = std::span<const T>(
                md.cell_reference_points.data(), md.cell_reference_points.size());
        }
        else
        {
            const auto& cached_dof_ref =
                cached_reference_lagrange_points<T>(ctype, degree);
            dof_ref = std::span<const T>(
                cached_dof_ref.data(), cached_dof_ref.size());
            // Without reference points, read the node positions off the dof
            // coordinates: pull them back through the cell's affine map (the
            // map the cut uses). Equispaced data reproduce the cached lattice
            // and keep it; data at other nodes (from degree 3 on, e.g. the
            // GLL-warped default of Basix/DOLFINx) get their own positions
            // instead of being read as equispaced values.
            if (static_cast<int>(dof_ref.size()) == ndofs * cell_ls.tdim)
            {
                pulled_back = pull_back_dof_coordinates<T, I>(
                    md, ctype, cell_ls.parent_vertex_coords, gdim, cell_ls.tdim,
                    cell_dofs);
                bool lattice = true;
                for (std::size_t i = 0; i < pulled_back.size(); ++i)
                {
                    if (std::abs(pulled_back[i] - dof_ref[i]) > T(1e-10))
                    {
                        lattice = false;
                        break;
                    }
                }
                if (!lattice)
                    dof_ref = std::span<const T>(pulled_back.data(), pulled_back.size());
            }
        }
        if (static_cast<int>(dof_ref.size()) != ndofs * cell_ls.tdim)
        {
            throw std::runtime_error(
                "make_cell_level_set: cached reference lattice size mismatch");
        }

        // Convert Lagrange nodal values to Bernstein coefficients
        cell_ls.bernstein_order = degree;
        bernstein::lagrange_to_bernstein(
            ctype, degree,
            dof_ref,
            std::span<const T>(cell_ls.nodal_values),
            cell_ls.bernstein_coeffs);
    }
}

// ---------------------------------------------------------------------------
// LevelSetCell::value
// ---------------------------------------------------------------------------

template <std::floating_point T, std::integral I>
T LevelSetCell<T, I>::value(std::span<const T> xi) const
{
    // Polynomial backend: evaluate the Bernstein expansion
    if (!bernstein_coeffs.empty())
    {
        return bernstein::evaluate(
            cell_type, bernstein_order,
            std::span<const T>(bernstein_coeffs), xi);
    }

    // Analytical fallback: delegate to the global level-set value in physical
    // coordinates after mapping xi to x.
    if (global_level_set != nullptr
        && global_level_set->has_value()
        && !parent_vertex_coords.empty())
    {
        const auto x = cell::push_forward_affine_map<T>(
            cell_type, parent_vertex_coords, gdim, xi);
        return global_level_set->value(x.data(), cell_id);
    }

    throw std::runtime_error(
        "LevelSetCell::value: no Bernstein coefficients available "
        "and no analytical fallback available");
}

// ---------------------------------------------------------------------------
// LevelSetCell::grad
// ---------------------------------------------------------------------------

template <std::floating_point T, std::integral I>
void LevelSetCell<T, I>::grad(std::span<const T> xi, std::span<T> g) const
{
    if (static_cast<int>(g.size()) != tdim)
    {
        throw std::runtime_error(
            "LevelSetCell::grad: output gradient has wrong dimension");
    }

    // Polynomial backend: gradient of the Bernstein expansion in reference
    // coordinates. This is the gradient of the local approximation phi_h.
    if (!bernstein_coeffs.empty())
    {
        bernstein::gradient(
            cell_type, bernstein_order,
            std::span<const T>(bernstein_coeffs), xi, g);
        return;
    }

    // Analytical fallback: LevelSetFunction::grad is physical dphi/dx.
    // Pull it back to the parent reference element as dphi/dxi = J^T dphi/dx.
    if (global_level_set != nullptr
        && global_level_set->has_gradient()
        && !parent_vertex_coords.empty())
    {
        const auto x = cell::push_forward_affine_map<T>(
            cell_type, parent_vertex_coords, gdim, xi);
        std::vector<T> grad_phys(static_cast<std::size_t>(gdim), T(0));
        global_level_set->grad(x.data(), cell_id, grad_phys.data());

        const auto cols = cell::jacobian_col_indices(cell_type);
        for (int a = 0; a < tdim; ++a)
        {
            const int va = cols[static_cast<std::size_t>(a)];
            if (va < 0)
            {
                throw std::runtime_error(
                    "LevelSetCell::grad: invalid affine Jacobian column");
            }
            T value = T(0);
            for (int r = 0; r < gdim; ++r)
            {
                const T j_ra =
                    parent_vertex_coords[
                        static_cast<std::size_t>(va * gdim + r)]
                  - parent_vertex_coords[static_cast<std::size_t>(r)];
                value += j_ra * grad_phys[static_cast<std::size_t>(r)];
            }
            g[static_cast<std::size_t>(a)] = value;
        }
        return;
    }

    throw std::runtime_error(
        "LevelSetCell::grad: no Bernstein coefficients available "
        "and no analytical fallback available");
}

// ---------------------------------------------------------------------------
// Explicit template instantiations
// ---------------------------------------------------------------------------

template LevelSetCell<double, int>
make_cell_level_set(const LevelSetFunction<double, int>&, int);

template LevelSetCell<float, int>
make_cell_level_set(const LevelSetFunction<float, int>&, int);

template LevelSetCell<double, long>
make_cell_level_set(const LevelSetFunction<double, long>&, long);

template LevelSetCell<float, long>
make_cell_level_set(const LevelSetFunction<float, long>&, long);
template void make_cell_level_set(const LevelSetFunction<double, int>&, int, LevelSetCell<double, int>&);
template void make_cell_level_set(const LevelSetFunction<float, int>&, int, LevelSetCell<float, int>&);
template void make_cell_level_set(const LevelSetFunction<double, long>&, long, LevelSetCell<double, long>&);
template void make_cell_level_set(const LevelSetFunction<float, long>&, long, LevelSetCell<float, long>&);

template double LevelSetCell<double, int>::value(std::span<const double>) const;
template float  LevelSetCell<float,  int>::value(std::span<const float>)  const;
template double LevelSetCell<double, long>::value(std::span<const double>) const;
template float  LevelSetCell<float,  long>::value(std::span<const float>)  const;

template void LevelSetCell<double, int>::grad(std::span<const double>, std::span<double>) const;
template void LevelSetCell<float,  int>::grad(std::span<const float>,  std::span<float>)  const;
template void LevelSetCell<double, long>::grad(std::span<const double>, std::span<double>) const;
template void LevelSetCell<float,  long>::grad(std::span<const float>,  std::span<float>)  const;

} // namespace cutcells
