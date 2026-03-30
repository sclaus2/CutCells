// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "quadrature_curved.h"
#include "quadrature.h"
#include "quadrature_tables.h"
#include "local_mesh.h"
#include "bernstein_backend.h"
#include "mapping.h"
#include "mapping_curved.h"
#include "cell_flags.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cutcells::quadrature
{

// Use mapping helpers from the shared header
using mapping::eval_pk_map;
using mapping::eval_pk_deriv;
using mapping::eval_collapsed_triangle;
using mapping::eval_gh_triangle_one_curved_edge;
using mapping::eval_gh_quad_one_curved_edge;
using mapping::eval_gh_tet_one_curved_face;
using mapping::eval_gh_tet_one_curved_tri_face;
using mapping::eval_gh_face_triangle;
using mapping::eval_gh_face_quad;
using mapping::face_surface_metric;
using mapping::CurvedMappingBackend;

namespace
{

template <std::floating_point T>
inline void append_physical_point(
    std::vector<T>* physical_points,
    const T*        x_phys,
    int             gdim)
{
    if (physical_points == nullptr)
        return;
    for (int d = 0; d < gdim; ++d)
        physical_points->push_back(x_phys[d]);
}

template <std::floating_point T>
inline void map_parent_ref_to_phys(
    const LocalMesh<T>& mesh,
    const T*            x_ref,
    T*                  x_phys)
{
    cell::ref_to_phys_affine<T>(
        mesh.parent_cell_type,
        std::span<const T>(mesh.parent_cell_coords_p1.data(),
                           mesh.parent_cell_coords_p1.size()),
        mesh.gdim,
        mesh.tdim,
        x_ref,
        x_phys);
}

template <std::floating_point T>
inline cell::type zero_face_type(const LocalMesh<T>& mesh, int zero_entity_id)
{
    const int z0 = mesh.zero_entity_offsets[static_cast<std::size_t>(zero_entity_id)];
    const int z1 = mesh.zero_entity_offsets[static_cast<std::size_t>(zero_entity_id + 1)];
    const int n_face_verts = z1 - z0;
    if (n_face_verts == 3)
        return cell::type::triangle;
    if (n_face_verts == 4)
        return cell::type::quadrilateral;
    return cell::type::point;
}

template <std::floating_point T>
inline void build_zero_face_coordinate_polys(
    const LocalMesh<T>&          mesh,
    int                          zero_entity_id,
    int                          geom_order,
    std::vector<BernsteinCell<T>>& coord_polys)
{
    const cell::type face_type = zero_face_type(mesh, zero_entity_id);
    if (face_type != cell::type::triangle
        && face_type != cell::type::quadrilateral)
        throw std::invalid_argument(
            "build_zero_face_coordinate_polys: unsupported zero-face type");

    const int tdim = mesh.tdim;
    const int z0 = mesh.zero_entity_offsets[static_cast<std::size_t>(zero_entity_id)];
    const int z1 = mesh.zero_entity_offsets[static_cast<std::size_t>(zero_entity_id + 1)];
    const int n_face_verts = z1 - z0;
    if ((face_type == cell::type::triangle && n_face_verts != 3)
        || (face_type == cell::type::quadrilateral && n_face_verts != 4))
        throw std::invalid_argument(
            "build_zero_face_coordinate_polys: inconsistent face vertex count");

    const bool has_curved =
        geom_order >= 2
        && !mesh.curved_zero_offsets.empty()
        && static_cast<int>(mesh.curved_zero_offsets.size()) == mesh.n_zero_entities() + 1;
    const int node_start = has_curved
        ? mesh.curved_zero_offsets[static_cast<std::size_t>(zero_entity_id)] : 0;
    const int node_end = has_curved
        ? mesh.curved_zero_offsets[static_cast<std::size_t>(zero_entity_id + 1)] : 0;
    const int n_face_interior = node_end - node_start;
    const int expected_face_interior = (face_type == cell::type::triangle)
        ? (geom_order - 1) * (geom_order - 2) / 2
        : (geom_order - 1) * (geom_order - 1);
    const bool use_curved = has_curved && n_face_interior == expected_face_interior;

    const int degree = use_curved ? geom_order : 1;
    const RefinementTemplate& tpl = (degree == 1)
        ? p1_template(face_type)
        : iso_p1_template(face_type, degree);

    std::vector<T> face_nodes_ref(
        static_cast<std::size_t>(tpl.n_vertices * tdim), T(0));

    std::array<const T*, 4> face_corners = {nullptr, nullptr, nullptr, nullptr};
    for (int i = 0; i < n_face_verts; ++i)
    {
        const int lv = mesh.zero_entity_vertices[static_cast<std::size_t>(z0 + i)];
        face_corners[static_cast<std::size_t>(i)]
            = &mesh.vertex_ref_x[static_cast<std::size_t>(lv * tdim)];
    }

    int interior_id = 0;
    for (int i = 0; i < tpl.n_vertices; ++i)
    {
        const T* cached = nullptr;
        if (use_curved && tpl.vertex_parent_dim[static_cast<std::size_t>(i)] == 2)
        {
            cached = &mesh.curved_zero_ref_nodes[
                static_cast<std::size_t>((node_start + interior_id) * tdim)];
            ++interior_id;
        }

        T* x_node = &face_nodes_ref[static_cast<std::size_t>(i * tdim)];
        if (cached != nullptr)
        {
            for (int d = 0; d < tdim; ++d)
                x_node[d] = cached[d];
            continue;
        }

        const T u = static_cast<T>(
            tpl.ref_vertex_coords[static_cast<std::size_t>(i * tpl.tdim + 0)]);
        const T v = static_cast<T>(
            tpl.ref_vertex_coords[static_cast<std::size_t>(i * tpl.tdim + 1)]);

        if (face_type == cell::type::triangle)
        {
            for (int d = 0; d < tdim; ++d)
                x_node[d] = (T(1) - u - v) * face_corners[0][d]
                          + u * face_corners[1][d]
                          + v * face_corners[2][d];
        }
        else
        {
            for (int d = 0; d < tdim; ++d)
                x_node[d] = (T(1) - u) * (T(1) - v) * face_corners[0][d]
                          + u * (T(1) - v) * face_corners[1][d]
                          + u * v * face_corners[2][d]
                          + (T(1) - u) * v * face_corners[3][d];
        }
    }

    coord_polys.assign(static_cast<std::size_t>(tdim), BernsteinCell<T>{});
    std::vector<T> lagrange_values(static_cast<std::size_t>(tpl.n_vertices), T(0));
    for (int d = 0; d < tdim; ++d)
    {
        for (int i = 0; i < tpl.n_vertices; ++i)
            lagrange_values[static_cast<std::size_t>(i)]
                = face_nodes_ref[static_cast<std::size_t>(i * tdim + d)];
        coord_polys[static_cast<std::size_t>(d)] = make_bernstein_cell<T>(
            face_type, degree,
            std::span<const T>(lagrange_values.data(), lagrange_values.size()));
    }
}

} // anonymous namespace

// ============================================================================
// append_interface_quadrature_curved
// ============================================================================

template <std::floating_point T>
void append_interface_quadrature_curved(
    const LocalMesh<T>& mesh,
    int zero_entity_id,
    int order,
    QuadratureRules<T>& rules)
{
    std::vector<T> ignored_physical_points;
    append_interface_quadrature_curved(
        mesh, zero_entity_id, order, rules, ignored_physical_points);
}

template <std::floating_point T>
void append_interface_quadrature_curved(
    const LocalMesh<T>& mesh,
    int zero_entity_id,
    int order,
    QuadratureRules<T>& rules,
    std::vector<T>& physical_points)
{
    const int gdim = mesh.gdim;
    const int tdim = mesh.tdim;
    const int codim1_dim = std::max(0, tdim - 1);

    // Validate entity id
    const int n_zero = mesh.n_zero_entities();
    if (zero_entity_id < 0 || zero_entity_id >= n_zero)
        return;
    if (mesh.zero_entity_dim[static_cast<std::size_t>(zero_entity_id)] != codim1_dim)
        return;

    // Geo order from cache, default to 1 (straight) if not built
    const int geom_order = (mesh.curved_geometry_order >= 1)
        ? mesh.curved_geometry_order : 1;

    const int z0 = mesh.zero_entity_offsets[static_cast<std::size_t>(zero_entity_id)];
    const int z1 = mesh.zero_entity_offsets[static_cast<std::size_t>(zero_entity_id + 1)];
    const int n_ent_verts = z1 - z0;
    if (n_ent_verts < 2)
        return;

    // 3D interface face quadrature.
    //
    // Phase D currently stores curved interior face nodes in the zero-entity
    // cache, while the boundary edges remain straight. For triangular zero
    // faces, build a regular degree-p lattice on the face, overwrite the
    // interior lattice nodes with the cached curved positions, and integrate
    // over the resulting piecewise-planar surface patch.
    if (tdim == 3 && (n_ent_verts == 3 || n_ent_verts == 4))
    {
        const cell::type face_type = zero_face_type(mesh, zero_entity_id);
        std::vector<BernsteinCell<T>> coord_polys;
        build_zero_face_coordinate_polys(
            mesh, zero_entity_id, geom_order, coord_polys);

        if (rules._tdim == 0)
            rules._tdim = tdim;
        if (rules._offset.empty())
            rules._offset.push_back(0);

        const auto& face_rule = get_reference_rule<T>(face_type, order);
        const int nq = face_rule._num_points;
        std::vector<T> J(static_cast<std::size_t>(gdim * tdim), T(0));

        for (int q = 0; q < nq; ++q)
        {
            const T* can_pt = face_rule._points.data()
                + static_cast<std::size_t>(q) * face_rule._tdim;
            T x_ref_q[3] = {T(0), T(0), T(0)};
            T dF_du_ref[3] = {T(0), T(0), T(0)};
            T dF_dv_ref[3] = {T(0), T(0), T(0)};
            if (face_type == cell::type::triangle)
            {
                eval_gh_face_triangle<T>(
                    std::span<const BernsteinCell<T>>(coord_polys.data(), coord_polys.size()),
                    can_pt[0], can_pt[1],
                    x_ref_q, dF_du_ref, dF_dv_ref);
            }
            else
            {
                eval_gh_face_quad<T>(
                    std::span<const BernsteinCell<T>>(coord_polys.data(), coord_polys.size()),
                    can_pt[0], can_pt[1],
                    x_ref_q, dF_du_ref, dF_dv_ref);
            }

            cell::compute_jacobian<T>(
                mesh.parent_cell_type,
                std::span<const T>(mesh.parent_cell_coords_p1.data(),
                                   mesh.parent_cell_coords_p1.size()),
                gdim,
                tdim,
                x_ref_q,
                J.data());

            const T jac = face_surface_metric<T>(
                std::span<const T>(J.data(), J.size()),
                dF_du_ref, dF_dv_ref, gdim, tdim);

            for (int d = 0; d < tdim; ++d)
                rules._points.push_back(x_ref_q[d]);
            T x_phys_q[3] = {};
            map_parent_ref_to_phys(mesh, x_ref_q, x_phys_q);
            append_physical_point(&physical_points, x_phys_q, gdim);
            rules._weights.push_back(face_rule._weights[q] * jac);
        }

        rules._parent_map.push_back(static_cast<int32_t>(mesh.parent_cell_id));
        rules._offset.push_back(static_cast<int32_t>(rules._weights.size()));
        return;
    }

    int lv0 = mesh.zero_entity_endpoint_v0[static_cast<std::size_t>(zero_entity_id)];
    int lv1 = mesh.zero_entity_endpoint_v1[static_cast<std::size_t>(zero_entity_id)];
    if (lv0 < 0 || lv1 < 0)
    {
        lv0 = mesh.zero_entity_vertices[static_cast<std::size_t>(z0)];
        lv1 = mesh.zero_entity_vertices[static_cast<std::size_t>(z0 + 1)];
    }

    const T* r_v0 = &mesh.vertex_ref_x[static_cast<std::size_t>(lv0 * tdim)];
    const T* r_v1 = &mesh.vertex_ref_x[static_cast<std::size_t>(lv1 * tdim)];

    // Interior curved nodes in parent reference coordinates.
    const bool has_interior =
        !mesh.curved_zero_offsets.empty()
        && static_cast<int>(mesh.curved_zero_offsets.size()) == n_zero + 1;

    const int node_start = has_interior
        ? mesh.curved_zero_offsets[static_cast<std::size_t>(zero_entity_id)] : 0;
    const int node_end = has_interior
        ? mesh.curved_zero_offsets[static_cast<std::size_t>(zero_entity_id + 1)] : 0;
    const int n_interior = node_end - node_start;

    // Expected n_interior for curved case
    const int expected_interior = geom_order - 1;
    const bool is_curved = has_interior && n_interior == expected_interior
                           && geom_order >= 2;

    std::vector<T> x_int_ref;
    if (is_curved)
    {
        x_int_ref.resize(static_cast<std::size_t>(n_interior * tdim));
        for (int ni = 0; ni < n_interior; ++ni)
        {
            const int ci = node_start + ni;
            const T* x_ref = &mesh.curved_zero_ref_nodes[
                static_cast<std::size_t>(ci * tdim)];
            T* x_ref_i = &x_int_ref[static_cast<std::size_t>(ni * tdim)];
            for (int d = 0; d < tdim; ++d)
                x_ref_i[d] = x_ref[d];
        }
    }

    // Gauss-Lobatto interior points for this geom order
    const auto gl_pts_vec =
        is_curved ? ::cutcells::detail::gauss_lobatto_interior_points_1d<T>(geom_order)
                  : std::vector<T>{};
    const T* gl_pts = is_curved ? gl_pts_vec.data() : nullptr;

    // ---- quadrature on the interval [0,1] ----
    const auto& quad_rule = get_reference_rule<T>(cell::type::interval, order);
    const int nq = quad_rule._num_points;

    // Quadrature points are on [0,1] (interval reference)
    const T* q_pts = quad_rule._points.data();
    const T* q_wts = quad_rule._weights.data();

    // Initialise output if needed
    if (rules._tdim == 0)
        rules._tdim = tdim;
    if (rules._offset.empty())
        rules._offset.push_back(0);

    for (int q = 0; q < nq; ++q)
    {
        const T s = q_pts[q];   // parameter in [0,1]

        // Parent-reference position and tangent of the curved segment.
        T x_ref_q[3] = {};
        T dx_ref_ds[3] = {};

        if (is_curved)
        {
            eval_pk_map(s,
                        r_v0, r_v1,
                        x_int_ref.data(),
                        gl_pts, geom_order, tdim,
                        x_ref_q);
            eval_pk_deriv(s,
                          r_v0, r_v1,
                          x_int_ref.data(),
                          gl_pts, geom_order, tdim,
                          dx_ref_ds);
        }
        else
        {
            for (int d = 0; d < tdim; ++d)
            {
                x_ref_q[d] = (T(1) - s) * r_v0[d] + s * r_v1[d];
                dx_ref_ds[d] = r_v1[d] - r_v0[d];
            }
        }

        // Surface Jacobian uses parent Jacobian at x_ref_q.
        std::vector<T> J(static_cast<std::size_t>(gdim * tdim), T(0));
        cell::compute_jacobian<T>(
            mesh.parent_cell_type,
            std::span<const T>(mesh.parent_cell_coords_p1.data(),
                               mesh.parent_cell_coords_p1.size()),
            gdim,
            tdim,
            x_ref_q,
            J.data());

        std::array<T, 3> dx_phys_ds = {T(0), T(0), T(0)};
        for (int d = 0; d < gdim; ++d)
        {
            T val = T(0);
            for (int k = 0; k < tdim; ++k)
                val += J[static_cast<std::size_t>(d * tdim + k)]
                     * dx_ref_ds[static_cast<std::size_t>(k)];
            dx_phys_ds[static_cast<std::size_t>(d)] = val;
        }

        T jac = T(0);
        for (int d = 0; d < gdim; ++d)
            jac += dx_phys_ds[static_cast<std::size_t>(d)]
                 * dx_phys_ds[static_cast<std::size_t>(d)];
        jac = std::sqrt(jac);

        for (int d = 0; d < tdim; ++d)
            rules._points.push_back(x_ref_q[d]);
        T x_phys_q[3] = {};
        map_parent_ref_to_phys(mesh, x_ref_q, x_phys_q);
        append_physical_point(&physical_points, x_phys_q, gdim);

        rules._weights.push_back(q_wts[q] * jac);
    }

    rules._parent_map.push_back(static_cast<int32_t>(mesh.parent_cell_id));
    rules._offset.push_back(static_cast<int32_t>(rules._weights.size()));
}

// ============================================================================
// append_volume_quadrature_curved
// ============================================================================

template <std::floating_point T>
void append_volume_quadrature_curved(
    const LocalMesh<T>& mesh,
    int cell_id,
    int level_set_id,
    int order,
    QuadratureRules<T>& rules,
    int& fallback_count)
{
    std::vector<T> ignored_physical_points;
    append_volume_quadrature_curved(
        mesh, cell_id, level_set_id, order, rules, fallback_count,
        ignored_physical_points, CurvedMappingBackend::collapsed);
}

template <std::floating_point T>
void append_volume_quadrature_curved(
    const LocalMesh<T>& mesh,
    int cell_id,
    int level_set_id,
    int order,
    QuadratureRules<T>& rules,
    int& fallback_count,
    std::vector<T>& physical_points)
{
    append_volume_quadrature_curved(
        mesh, cell_id, level_set_id, order, rules, fallback_count,
        physical_points, CurvedMappingBackend::collapsed);
}

template <std::floating_point T>
void append_volume_quadrature_curved(
    const LocalMesh<T>& mesh,
    int cell_id,
    int level_set_id,
    int order,
    QuadratureRules<T>& rules,
    int& fallback_count,
    CurvedMappingBackend vol_backend)
{
    std::vector<T> ignored_physical_points;
    append_volume_quadrature_curved(
        mesh, cell_id, level_set_id, order, rules, fallback_count,
        ignored_physical_points, vol_backend);
}

template <std::floating_point T>
void append_volume_quadrature_curved(
    const LocalMesh<T>& mesh,
    int cell_id,
    int level_set_id,
    int order,
    QuadratureRules<T>& rules,
    int& fallback_count,
    std::vector<T>& physical_points,
    CurvedMappingBackend vol_backend)
{
    const int gdim  = mesh.gdim;
    const int tdim  = mesh.tdim;
    const int nc    = mesh.n_cells();

    if (cell_id < 0 || cell_id >= nc)
        return;

    // ---- check domain ----
    const uint8_t domain_val = mesh.cell_domain[static_cast<std::size_t>(cell_id)];
    if (domain_val != static_cast<uint8_t>(cell::domain::inside))
        return;

    // Initialise output if needed
    if (rules._tdim == 0)
        rules._tdim = tdim;
    if (rules._offset.empty())
        rules._offset.push_back(0);

    const cell::type ctype = mesh.cell_types[static_cast<std::size_t>(cell_id)];
    const int c0 = mesh.cell_offsets[static_cast<std::size_t>(cell_id)];
    const int c1 = mesh.cell_offsets[static_cast<std::size_t>(cell_id + 1)];
    const int nv = c1 - c0;

    // Gather physical vertex coordinates for this sub-cell
    std::vector<T> phys_verts(static_cast<std::size_t>(nv * gdim));
    std::vector<T> ref_verts(static_cast<std::size_t>(nv * tdim));
    for (int j = 0; j < nv; ++j)
    {
        const int lv = mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];
        for (int d = 0; d < gdim; ++d)
            phys_verts[static_cast<std::size_t>(j * gdim + d)] =
                mesh.vertex_x[static_cast<std::size_t>(lv * gdim + d)];
        for (int d = 0; d < tdim; ++d)
            ref_verts[static_cast<std::size_t>(j * tdim + d)] =
                mesh.vertex_ref_x[static_cast<std::size_t>(lv * tdim + d)];
    }

    auto append_affine_subcell = [&]()
    {
        const auto& ref_rule = get_reference_rule<T>(ctype, order);
        const int nq = ref_rule._num_points;

        if (tdim == 2 && ctype == cell::type::triangle && nv == 3)
        {
            const T* r0 = &ref_verts[0];
            const T* r1 = &ref_verts[static_cast<std::size_t>(tdim)];
            const T* r2 = &ref_verts[static_cast<std::size_t>(2 * tdim)];
            const T a0 = r1[0] - r0[0];
            const T a1 = r1[1] - r0[1];
            const T b0 = r2[0] - r0[0];
            const T b1 = r2[1] - r0[1];

            for (int q = 0; q < nq; ++q)
            {
                const T* can_pt = ref_rule._points.data()
                    + static_cast<std::size_t>(q) * ref_rule._tdim;
                const T xi = can_pt[0];
                const T eta = can_pt[1];

                T x_ref_q[3] = {T(0), T(0), T(0)};
                x_ref_q[0] = r0[0] + xi * a0 + eta * b0;
                x_ref_q[1] = r0[1] + xi * a1 + eta * b1;

                std::vector<T> J(static_cast<std::size_t>(gdim * tdim), T(0));
                cell::compute_jacobian<T>(
                    mesh.parent_cell_type,
                    std::span<const T>(mesh.parent_cell_coords_p1.data(),
                                       mesh.parent_cell_coords_p1.size()),
                    gdim,
                    tdim,
                    x_ref_q,
                    J.data());

                std::array<T, 3> col0 = {T(0), T(0), T(0)};
                std::array<T, 3> col1 = {T(0), T(0), T(0)};
                for (int d = 0; d < gdim; ++d)
                {
                    col0[static_cast<std::size_t>(d)] =
                        J[static_cast<std::size_t>(d * tdim + 0)] * a0
                        + J[static_cast<std::size_t>(d * tdim + 1)] * a1;
                    col1[static_cast<std::size_t>(d)] =
                        J[static_cast<std::size_t>(d * tdim + 0)] * b0
                        + J[static_cast<std::size_t>(d * tdim + 1)] * b1;
                }

                T det = T(0);
                if (gdim == 2)
                {
                    det = std::abs(col0[0] * col1[1] - col0[1] * col1[0]);
                }
                else if (gdim == 3)
                {
                    const T cx = col0[1] * col1[2] - col0[2] * col1[1];
                    const T cy = col0[2] * col1[0] - col0[0] * col1[2];
                    const T cz = col0[0] * col1[1] - col0[1] * col1[0];
                    det = std::sqrt(cx * cx + cy * cy + cz * cz);
                }
                else
                {
                    det = T(0);
                }

                rules._points.push_back(x_ref_q[0]);
                rules._points.push_back(x_ref_q[1]);
                T x_phys_q[3] = {};
                map_parent_ref_to_phys(mesh, x_ref_q, x_phys_q);
                append_physical_point(&physical_points, x_phys_q, gdim);
                rules._weights.push_back(ref_rule._weights[q] * det);
            }
        }
        else
        {
            const T det_J = cell::affine_volume_factor<T>(
                ctype, phys_verts.data(), gdim);
            for (int q = 0; q < nq; ++q)
            {
                const T* can_pt = ref_rule._points.data()
                    + static_cast<std::size_t>(q) * ref_rule._tdim;

                T x_phys_q[3] = {};
                const int sdim = ref_rule._tdim;
                for (int d = 0; d < gdim; ++d)
                    x_phys_q[d] = phys_verts[d];
                for (int i = 1; i <= sdim; ++i)
                {
                    for (int d = 0; d < gdim; ++d)
                        x_phys_q[d] += can_pt[i - 1]
                            * (phys_verts[static_cast<std::size_t>(i * gdim + d)]
                               - phys_verts[static_cast<std::size_t>(d)]);
                }

                T x_ref_q[3] = {};
                cell::pull_back_affine<T>(
                    mesh.parent_cell_type,
                    mesh.parent_cell_coords_p1,
                    gdim,
                    std::span<const T>(x_phys_q, static_cast<std::size_t>(gdim)),
                    std::span<T>(x_ref_q, static_cast<std::size_t>(tdim)));
                for (int d = 0; d < tdim; ++d)
                    rules._points.push_back(x_ref_q[d]);
                append_physical_point(&physical_points, x_phys_q, gdim);
                rules._weights.push_back(ref_rule._weights[q] * det_J);
            }
        }

        rules._parent_map.push_back(static_cast<int32_t>(mesh.parent_cell_id));
        rules._offset.push_back(static_cast<int32_t>(rules._weights.size()));
    };

    const bool has_curved = mesh.curved_geometry_order >= 2
        && !mesh.curved_zero_offsets.empty();

    int touching_ze = -1;
    int cell_if_a = -1;
    int cell_if_b = -1;
    bool iface_reversed_for_cell = false;
    // For 3D: face zero entity vertex IDs
    std::vector<int32_t> touching_face_verts;

    if (has_curved && tdim == 2 && ctype == cell::type::triangle)
    {
        const uint64_t zero_mask = uint64_t(1) << level_set_id;
        const int n_zero = mesh.n_zero_entities();
        const int codim1_dim = std::max(0, mesh.tdim - 1);
        for (int ze = 0; ze < n_zero && touching_ze < 0; ++ze)
        {
            if (mesh.zero_entity_dim[static_cast<std::size_t>(ze)] != codim1_dim)
                continue;
            if ((mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] & zero_mask) == 0)
                continue;

            const int if_a = mesh.zero_entity_endpoint_v0[static_cast<std::size_t>(ze)];
            const int if_b = mesh.zero_entity_endpoint_v1[static_cast<std::size_t>(ze)];
            if (if_a < 0 || if_b < 0)
                continue;

            // Check each edge of the triangle
            for (int e = 0; e < nv; ++e)
            {
                const int ea = mesh.cell_vertices[static_cast<std::size_t>(c0 + e)];
                const int eb = mesh.cell_vertices[static_cast<std::size_t>(c0 + (e + 1) % nv)];

                if ((ea == if_a && eb == if_b) || (ea == if_b && eb == if_a))
                {
                    touching_ze = ze;
                    cell_if_a = ea;
                    cell_if_b = eb;
                    iface_reversed_for_cell = (ea == if_b && eb == if_a);
                    break;
                }
            }
        }
    }
    else if (has_curved && tdim == 3 && ctype == cell::type::tetrahedron)
    {
        // 3D: check if any face of this tet matches a face zero entity
        const uint64_t zero_mask = uint64_t(1) << level_set_id;
        const int n_zero = mesh.n_zero_entities();

        // Gather this cell's vertices as a set for face matching
        std::vector<int32_t> cell_verts(static_cast<std::size_t>(nv));
        for (int j = 0; j < nv; ++j)
            cell_verts[static_cast<std::size_t>(j)]
                = mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];

        for (int ze = 0; ze < n_zero && touching_ze < 0; ++ze)
        {
            if (mesh.zero_entity_dim[static_cast<std::size_t>(ze)] != 2)
                continue;
            if ((mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] & zero_mask) == 0)
                continue;

            // Get face vertices
            const int fe0 = mesh.zero_entity_offsets[static_cast<std::size_t>(ze)];
            const int fe1 = mesh.zero_entity_offsets[static_cast<std::size_t>(ze + 1)];
            const int n_fv = fe1 - fe0;

            // Check if all face vertices belong to this cell
            bool all_in_cell = true;
            std::vector<int32_t> fv(static_cast<std::size_t>(n_fv));
            for (int i = 0; i < n_fv; ++i)
            {
                fv[static_cast<std::size_t>(i)]
                    = mesh.zero_entity_vertices[static_cast<std::size_t>(fe0 + i)];
                bool found = false;
                for (int j = 0; j < nv; ++j)
                {
                    if (fv[static_cast<std::size_t>(i)] == cell_verts[static_cast<std::size_t>(j)])
                    {
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    all_in_cell = false;
                    break;
                }
            }
            if (all_in_cell && n_fv >= 3)
            {
                touching_ze = ze;
                touching_face_verts = fv;
            }
        }
    }

    // ---- no interface contact or unsupported cell type: affine path ----
    const bool can_do_curved_2d = (tdim == 2 && ctype == cell::type::triangle
                                   && touching_ze >= 0);
    const bool can_do_curved_3d = (tdim == 3 && ctype == cell::type::tetrahedron
                                   && touching_ze >= 0
                                   && !touching_face_verts.empty());
    if (!has_curved || (!can_do_curved_2d && !can_do_curved_3d))
    {
        append_affine_subcell();
        return;
    }

    if (can_do_curved_3d)
    {
        // ---- 3D: GH tet quadrature for tet touching one curved triangular face ----
        const int ze = touching_ze;

        // Identify the face vertices and the opposite (apex) vertex
        const int n_face_verts = static_cast<int>(touching_face_verts.size());
        if (n_face_verts != 3)
        {
            // Only triangular faces supported for now
            append_affine_subcell();
            return;
        }

        // Find the apex vertex (not on the face)
        int apex_vid = -1;
        for (int j = 0; j < nv; ++j)
        {
            const int lv = mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];
            bool on_face = false;
            for (int i = 0; i < n_face_verts; ++i)
            {
                if (lv == touching_face_verts[static_cast<std::size_t>(i)])
                {
                    on_face = true;
                    break;
                }
            }
            if (!on_face)
            {
                apex_vid = lv;
                break;
            }
        }
        if (apex_vid < 0)
        {
            append_affine_subcell();
            return;
        }

        const T* r_apex = &mesh.vertex_ref_x[static_cast<std::size_t>(apex_vid * tdim)];
        std::vector<BernsteinCell<T>> face_coord_polys;
        build_zero_face_coordinate_polys(mesh, ze, mesh.curved_geometry_order, face_coord_polys);

        const auto& tet_rule = get_reference_rule<T>(cell::type::tetrahedron, order);
        const int nq = tet_rule._num_points;

        bool jac_valid = true;
        int det_sign_3d = 0;
        std::vector<T> J(static_cast<std::size_t>(gdim * tdim), T(0));
        const int n_pts_total = nq;
        std::vector<T> pts_out, wts_out;
        pts_out.reserve(static_cast<std::size_t>(n_pts_total * tdim));
        wts_out.reserve(static_cast<std::size_t>(n_pts_total));

        for (int q = 0; q < nq && jac_valid; ++q)
        {
            const T* can_pt = tet_rule._points.data()
                + static_cast<std::size_t>(q) * tet_rule._tdim;
            const T xi_q = can_pt[0];
            const T eta_q = can_pt[1];
            const T zeta_q = can_pt[2];

            T F_ref[3] = {};
            T dF_dxi[3] = {};
            T dF_deta[3] = {};
            T dF_dzeta[3] = {};

            eval_gh_tet_one_curved_tri_face<T>(
                xi_q, eta_q, zeta_q,
                std::span<const BernsteinCell<T>>(face_coord_polys.data(), face_coord_polys.size()),
                r_apex,
                F_ref, dF_dxi, dF_deta, dF_dzeta);

            cell::compute_jacobian<T>(
                mesh.parent_cell_type,
                std::span<const T>(mesh.parent_cell_coords_p1.data(),
                                   mesh.parent_cell_coords_p1.size()),
                gdim, tdim, F_ref, J.data());

            std::array<T, 3> col_xi = {T(0), T(0), T(0)};
            std::array<T, 3> col_eta = {T(0), T(0), T(0)};
            std::array<T, 3> col_zeta = {T(0), T(0), T(0)};
            for (int d = 0; d < gdim; ++d)
            {
                for (int k = 0; k < tdim; ++k)
                {
                    const T Jdk = J[static_cast<std::size_t>(d * tdim + k)];
                    col_xi[static_cast<std::size_t>(d)]   += Jdk * dF_dxi[k];
                    col_eta[static_cast<std::size_t>(d)]  += Jdk * dF_deta[k];
                    col_zeta[static_cast<std::size_t>(d)] += Jdk * dF_dzeta[k];
                }
            }

            const T det_J = mapping::det_3d(
                col_xi.data(), col_eta.data(), col_zeta.data());
            const T det_abs = std::abs(det_J);
            const T det_tol = T(64) * std::numeric_limits<T>::epsilon();

            if (det_abs > det_tol)
            {
                const int sgn = (det_J > T(0)) ? 1 : -1;
                if (det_sign_3d == 0)
                    det_sign_3d = sgn;
                else if (sgn != det_sign_3d)
                {
                    jac_valid = false;
                    break;
                }
            }

            for (int d = 0; d < tdim; ++d)
                pts_out.push_back(F_ref[d]);
            T x_phys_q[3] = {};
            map_parent_ref_to_phys(mesh, F_ref, x_phys_q);
            append_physical_point(&physical_points, x_phys_q, gdim);
            wts_out.push_back(tet_rule._weights[q] * det_abs);
        }

        if (!jac_valid)
        {
            ++fallback_count;
            append_affine_subcell();
            return;
        }

        rules._points.insert(rules._points.end(), pts_out.begin(), pts_out.end());
        rules._weights.insert(rules._weights.end(), wts_out.begin(), wts_out.end());
        rules._parent_map.push_back(static_cast<int32_t>(mesh.parent_cell_id));
        rules._offset.push_back(static_cast<int32_t>(rules._weights.size()));
        return;
    }

    // ---- blended Pk quadrature for triangle touching the interface ----
    //
    // The triangle has 3 vertices. One edge [v_a, v_b] lies on the interface.
    // The opposite vertex is v_opp.
    //
    // Blended map F(s, t) = (1-t) * gamma(s) + t * v_opp
    //   where gamma(s) is the Pk curved interface map,
    //         t ∈ [0,1] controls blending from interface (t=0) to opposite vertex (t=1).
    //   s ∈ [0,1] parametrises the interface edge.
    //
    // The Jacobian is:
    //   J = [d F/d s, d F/d t] = [(1-t)*dgamma/ds, v_opp - gamma(s)]
    //
    // Integration weight = |det J| * w_q (2D case).
    //
    // The (s,t) quadrature is built as a product of 1D Gauss rules on [0,1]^2,
    // NOT a triangle rule, because the collapsed triangle is degenerate at t=1.
    // We use `order` points in each direction.

    const int ze = touching_ze;

    // Identify the opposite vertex (not on the zero entity segment)
    const int if_a = (cell_if_a >= 0) ? cell_if_a
                                      : mesh.zero_entity_endpoint_v0[static_cast<std::size_t>(ze)];
    const int if_b = (cell_if_b >= 0) ? cell_if_b
                                      : mesh.zero_entity_endpoint_v1[static_cast<std::size_t>(ze)];

    // Interface edge endpoints in parent reference coordinates.
    const T* r_va = &mesh.vertex_ref_x[static_cast<std::size_t>(if_a * tdim)];
    const T* r_vb = &mesh.vertex_ref_x[static_cast<std::size_t>(if_b * tdim)];

    // Find v_opp (the triangle vertex that is NOT if_a or if_b)
    int opp_idx = -1;
    for (int j = 0; j < nv; ++j)
    {
        const int lv = mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];
        if (lv != if_a && lv != if_b)
        {
            opp_idx = j;
            break;
        }
    }

    if (opp_idx < 0)
    {
        append_affine_subcell();
        return;
    }

    const int lv_opp = mesh.cell_vertices[static_cast<std::size_t>(c0 + opp_idx)];
    const T* r_opp = &mesh.vertex_ref_x[static_cast<std::size_t>(lv_opp * tdim)];

    // Interior curved nodes in parent reference coordinates.
    const int node_start = mesh.curved_zero_offsets[static_cast<std::size_t>(ze)];
    const int node_end   = mesh.curved_zero_offsets[static_cast<std::size_t>(ze + 1)];
    const int n_int = node_end - node_start;
    const int geom_order = mesh.curved_geometry_order;

    std::vector<T> x_int_ref(static_cast<std::size_t>(n_int * tdim));
    for (int ni = 0; ni < n_int; ++ni)
    {
        const int src_ni = iface_reversed_for_cell ? (n_int - 1 - ni) : ni;
        const int ci = node_start + src_ni;
        const T* x_ref_i = &mesh.curved_zero_ref_nodes[
            static_cast<std::size_t>(ci * tdim)];
        T* rp = &x_int_ref[static_cast<std::size_t>(ni * tdim)];
        for (int d = 0; d < tdim; ++d)
            rp[d] = x_ref_i[d];
    }

    const auto gl_pts_vec = ::cutcells::detail::gauss_lobatto_interior_points_1d<T>(geom_order);
    const T* gl_pts = gl_pts_vec.data();

    // Product Gauss rule on [0,1]^2
    const auto& quad1d = get_reference_rule<T>(cell::type::interval, order);
    const int nq1 = quad1d._num_points;

    // Check if iface edge vertices are ordered consistently with gamma(s) = v0+s*(v1-v0).
    // gamma(0) = x_va, gamma(1) = x_vb regardless of Lagrange basis endpoint ordering
    // (L_0(0) = 1, L_p(1) = 1).

    bool jac_valid = true;
    int det_sign = 0; // 0 unset, +1 positive orientation, -1 negative orientation

    const int n_pts_total = nq1 * nq1;
    std::vector<T> pts_out, wts_out;
    pts_out.reserve(static_cast<std::size_t>(n_pts_total * tdim));
    wts_out.reserve(static_cast<std::size_t>(n_pts_total));

    for (int qi = 0; qi < nq1; ++qi)
    {
        const T s = quad1d._points[qi];
        const T ws = quad1d._weights[qi];

        for (int qj = 0; qj < nq1; ++qj)
        {
            const T t = quad1d._points[qj];
            const T wt = quad1d._weights[qj];

            // Evaluate mapping in parent reference space.
            // Both backends produce identical results for triangles with
            // one curved edge. The GH formulation generalises to quads
            // and multiple curved edges.
            T F_ref[3] = {};
            T col_s_ref[3] = {};
            T col_t_ref[3] = {};

            if (vol_backend == CurvedMappingBackend::gordon_hall)
            {
                eval_gh_triangle_one_curved_edge(
                    s, t,
                    r_va, r_vb,
                    x_int_ref.data(), gl_pts, geom_order,
                    r_opp, tdim,
                    F_ref, col_s_ref, col_t_ref);
            }
            else
            {
                eval_collapsed_triangle(
                    s, t,
                    r_va, r_vb,
                    x_int_ref.data(), gl_pts, geom_order,
                    r_opp, tdim,
                    F_ref, col_s_ref, col_t_ref);
            }

            // Parent Jacobian J_parent(F_ref), then map columns to physical:
            //   col_*^phys = J_parent * col_*^ref.
            std::vector<T> J(static_cast<std::size_t>(gdim * tdim), T(0));
            cell::compute_jacobian<T>(
                mesh.parent_cell_type,
                std::span<const T>(mesh.parent_cell_coords_p1.data(),
                                   mesh.parent_cell_coords_p1.size()),
                gdim,
                tdim,
                F_ref,
                J.data());

            std::array<T, 3> col_s_phys = {T(0), T(0), T(0)};
            std::array<T, 3> col_t_phys = {T(0), T(0), T(0)};
            for (int d = 0; d < gdim; ++d)
            {
                T vs = T(0);
                T vt = T(0);
                for (int k = 0; k < tdim; ++k)
                {
                    const T Jdk = J[static_cast<std::size_t>(d * tdim + k)];
                    vs += Jdk * col_s_ref[k];
                    vt += Jdk * col_t_ref[k];
                }
                col_s_phys[static_cast<std::size_t>(d)] = vs;
                col_t_phys[static_cast<std::size_t>(d)] = vt;
            }

            T det_J = T(0);
            if (gdim == 2)
            {
                const T det_oriented = col_s_phys[0] * col_t_phys[1]
                                     - col_s_phys[1] * col_t_phys[0];
                const T det_abs = std::abs(det_oriented);
                const T det_tol = T(64) * std::numeric_limits<T>::epsilon();

                // Accept either global orientation (CW or CCW), but reject
                // sign flips across quadrature points (fold-over).
                if (det_abs > det_tol)
                {
                    const int sgn = (det_oriented > T(0)) ? 1 : -1;
                    if (det_sign == 0)
                        det_sign = sgn;
                    else if (sgn != det_sign)
                    {
                        jac_valid = false;
                        break;
                    }
                }
                det_J = det_abs;
            }
            else if (gdim == 3)
            {
                const T cx = col_s_phys[1] * col_t_phys[2] - col_s_phys[2] * col_t_phys[1];
                const T cy = col_s_phys[2] * col_t_phys[0] - col_s_phys[0] * col_t_phys[2];
                const T cz = col_s_phys[0] * col_t_phys[1] - col_s_phys[1] * col_t_phys[0];
                det_J = std::sqrt(cx * cx + cy * cy + cz * cz);
            }

            for (int d = 0; d < tdim; ++d)
                pts_out.push_back(F_ref[d]);
            T x_phys_q[3] = {};
            map_parent_ref_to_phys(mesh, F_ref, x_phys_q);
            append_physical_point(&physical_points, x_phys_q, gdim);
            wts_out.push_back(ws * wt * det_J);
        }
        if (!jac_valid)
            break;
    }

    if (!jac_valid)
    {
        ++fallback_count;
        append_affine_subcell();
        return;
    }

    rules._points.insert(rules._points.end(), pts_out.begin(), pts_out.end());
    rules._weights.insert(rules._weights.end(), wts_out.begin(), wts_out.end());
    rules._parent_map.push_back(static_cast<int32_t>(mesh.parent_cell_id));
    rules._offset.push_back(static_cast<int32_t>(rules._weights.size()));
}

// ============================================================================
// make_quadrature_curved
// ============================================================================

template <std::floating_point T>
void make_quadrature_curved(
    const LocalMesh<T>& mesh,
    int level_set_id,
    int order,
    QuadratureRules<T>& volume_rules,
    QuadratureRules<T>& interface_rules)
{
    std::vector<T> ignored_volume_points;
    std::vector<T> ignored_interface_points;
    make_quadrature_curved(
        mesh,
        level_set_id,
        order,
        volume_rules,
        ignored_volume_points,
        interface_rules,
        ignored_interface_points,
        CurvedMappingBackend::collapsed);
}

template <std::floating_point T>
void make_quadrature_curved(
    const LocalMesh<T>& mesh,
    int level_set_id,
    int order,
    QuadratureRules<T>& volume_rules,
    std::vector<T>& volume_physical_points,
    QuadratureRules<T>& interface_rules,
    std::vector<T>& interface_physical_points)
{
    make_quadrature_curved(
        mesh, level_set_id, order,
        volume_rules, volume_physical_points,
        interface_rules, interface_physical_points,
        CurvedMappingBackend::collapsed);
}

template <std::floating_point T>
void make_quadrature_curved(
    const LocalMesh<T>& mesh,
    int level_set_id,
    int order,
    QuadratureRules<T>& volume_rules,
    QuadratureRules<T>& interface_rules,
    CurvedMappingBackend volume_backend)
{
    std::vector<T> ignored_volume_points;
    std::vector<T> ignored_interface_points;
    make_quadrature_curved(
        mesh,
        level_set_id,
        order,
        volume_rules,
        ignored_volume_points,
        interface_rules,
        ignored_interface_points,
        volume_backend);
}

template <std::floating_point T>
void make_quadrature_curved(
    const LocalMesh<T>& mesh,
    int level_set_id,
    int order,
    QuadratureRules<T>& volume_rules,
    std::vector<T>& volume_physical_points,
    QuadratureRules<T>& interface_rules,
    std::vector<T>& interface_physical_points,
    CurvedMappingBackend volume_backend)
{
    // Reset outputs
    volume_rules    = {};
    interface_rules = {};
    volume_physical_points.clear();
    interface_physical_points.clear();

    const int nc = mesh.n_cells();
    const int n_zero = mesh.n_zero_entities();
    const uint64_t zero_mask = uint64_t(1) << level_set_id;
    const int codim1_dim = std::max(0, mesh.tdim - 1);

    // Interface quadrature (independent of volume mapping backend)
    for (int ze = 0; ze < n_zero; ++ze)
    {
        if (mesh.zero_entity_dim[static_cast<std::size_t>(ze)] != codim1_dim)
            continue;
        if ((mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] & zero_mask) == 0)
            continue;
        if (static_cast<std::size_t>(ze) < mesh.zero_entity_is_owned.size()
            && mesh.zero_entity_is_owned[static_cast<std::size_t>(ze)] == uint8_t(0))
            continue;

        append_interface_quadrature_curved(
            mesh, ze, order, interface_rules, interface_physical_points);
    }

    // Volume quadrature over inside sub-cells
    int fallback_count = 0;
    for (int c = 0; c < nc; ++c)
    {
        const uint8_t dom = mesh.cell_domain[static_cast<std::size_t>(c)];
        if (dom != static_cast<uint8_t>(cell::domain::inside))
            continue;
        append_volume_quadrature_curved(
            mesh, c, level_set_id, order, volume_rules, fallback_count,
            volume_physical_points, volume_backend);
    }
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template void append_interface_quadrature_curved<float>(
    const LocalMesh<float>&, int, int, QuadratureRules<float>&);
template void append_interface_quadrature_curved<double>(
    const LocalMesh<double>&, int, int, QuadratureRules<double>&);
template void append_interface_quadrature_curved<float>(
    const LocalMesh<float>&, int, int, QuadratureRules<float>&, std::vector<float>&);
template void append_interface_quadrature_curved<double>(
    const LocalMesh<double>&, int, int, QuadratureRules<double>&, std::vector<double>&);

template void append_volume_quadrature_curved<float>(
    const LocalMesh<float>&, int, int, int, QuadratureRules<float>&, int&);
template void append_volume_quadrature_curved<double>(
    const LocalMesh<double>&, int, int, int, QuadratureRules<double>&, int&);
template void append_volume_quadrature_curved<float>(
    const LocalMesh<float>&, int, int, int, QuadratureRules<float>&, int&, std::vector<float>&);
template void append_volume_quadrature_curved<double>(
    const LocalMesh<double>&, int, int, int, QuadratureRules<double>&, int&, std::vector<double>&);
template void append_volume_quadrature_curved<float>(
    const LocalMesh<float>&, int, int, int, QuadratureRules<float>&, int&, CurvedMappingBackend);
template void append_volume_quadrature_curved<double>(
    const LocalMesh<double>&, int, int, int, QuadratureRules<double>&, int&, CurvedMappingBackend);
template void append_volume_quadrature_curved<float>(
    const LocalMesh<float>&, int, int, int, QuadratureRules<float>&, int&, std::vector<float>&, CurvedMappingBackend);
template void append_volume_quadrature_curved<double>(
    const LocalMesh<double>&, int, int, int, QuadratureRules<double>&, int&, std::vector<double>&, CurvedMappingBackend);

template void make_quadrature_curved<float>(
    const LocalMesh<float>&, int, int,
    QuadratureRules<float>&, QuadratureRules<float>&);
template void make_quadrature_curved<double>(
    const LocalMesh<double>&, int, int,
    QuadratureRules<double>&, QuadratureRules<double>&);
template void make_quadrature_curved<float>(
    const LocalMesh<float>&, int, int,
    QuadratureRules<float>&, std::vector<float>&,
    QuadratureRules<float>&, std::vector<float>&);
template void make_quadrature_curved<double>(
    const LocalMesh<double>&, int, int,
    QuadratureRules<double>&, std::vector<double>&,
    QuadratureRules<double>&, std::vector<double>&);
template void make_quadrature_curved<float>(
    const LocalMesh<float>&, int, int,
    QuadratureRules<float>&, QuadratureRules<float>&,
    CurvedMappingBackend);
template void make_quadrature_curved<double>(
    const LocalMesh<double>&, int, int,
    QuadratureRules<double>&, QuadratureRules<double>&,
    CurvedMappingBackend);
template void make_quadrature_curved<float>(
    const LocalMesh<float>&, int, int,
    QuadratureRules<float>&, std::vector<float>&,
    QuadratureRules<float>&, std::vector<float>&,
    CurvedMappingBackend);
template void make_quadrature_curved<double>(
    const LocalMesh<double>&, int, int,
    QuadratureRules<double>&, std::vector<double>&,
    QuadratureRules<double>&, std::vector<double>&,
    CurvedMappingBackend);

} // namespace cutcells::quadrature
