// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "green_split.h"
#include "bernstein_backend.h"

#include <cmath>
#include <stdexcept>

namespace cutcells
{

template <std::floating_point T>
std::vector<int> cells_incident_on_edge(const LocalMesh<T>& mesh, int edge_id)
{
    std::vector<int> result;
    const int nc = mesh.n_cells();
    for (int c = 0; c < nc; ++c)
    {
        const int e0 = mesh.cell_edge_offsets[static_cast<std::size_t>(c)];
        const int e1 = mesh.cell_edge_offsets[static_cast<std::size_t>(c + 1)];
        for (int j = e0; j < e1; ++j)
        {
            if (mesh.cell_edges_flat[static_cast<std::size_t>(j)] == edge_id)
            {
                result.push_back(c);
                break;
            }
        }
    }
    return result;
}

template <std::floating_point T>
bool green_split_one_edge(LocalMesh<T>& mesh, int edge_id, T t_sep)
{
    const int32_t ev0 = mesh.edge_vertices[static_cast<std::size_t>(2 * edge_id)];
    const int32_t ev1 = mesh.edge_vertices[static_cast<std::size_t>(2 * edge_id + 1)];

    // Compute separator vertex coordinates (reference and physical)
    std::array<T, 3> x_ref  = {T(0), T(0), T(0)};
    std::array<T, 3> x_phys = {T(0), T(0), T(0)};
    for (int d = 0; d < mesh.tdim; ++d)
    {
        const T r0 = mesh.vertex_ref_x[static_cast<std::size_t>(ev0 * mesh.tdim + d)];
        const T r1 = mesh.vertex_ref_x[static_cast<std::size_t>(ev1 * mesh.tdim + d)];
        x_ref[static_cast<std::size_t>(d)] = r0 + t_sep * (r1 - r0);
    }
    for (int d = 0; d < mesh.gdim; ++d)
    {
        const T p0 = mesh.vertex_x[static_cast<std::size_t>(ev0 * mesh.gdim + d)];
        const T p1 = mesh.vertex_x[static_cast<std::size_t>(ev1 * mesh.gdim + d)];
        x_phys[static_cast<std::size_t>(d)] = p0 + t_sep * (p1 - p0);
    }

    // Infer parent-entity tracking for the new vertex
    const auto [pdim, pid] = infer_background_edge_parent(mesh, ev0, ev1);

    // Identify incident cells before the mesh is modified
    const auto incident_cells = cells_incident_on_edge(mesh, edge_id);
    if (incident_cells.empty())
        return true;

    // Prism and pyramid splits are not yet implemented — fall back to red
    // refinement for these cell types (TODO: add later).
    for (const int c : incident_cells)
    {
        const cell::type ct = mesh.cell_types[static_cast<std::size_t>(c)];
        if (ct == cell::type::prism || ct == cell::type::pyramid)
            return false;
    }

    // Append the separator vertex
    const int32_t mid_v = static_cast<int32_t>(mesh.vertex_parent_dim.size());
    for (int d = 0; d < mesh.gdim; ++d)
        mesh.vertex_x.push_back(x_phys[static_cast<std::size_t>(d)]);
    for (int d = 0; d < mesh.tdim; ++d)
        mesh.vertex_ref_x.push_back(x_ref[static_cast<std::size_t>(d)]);
    mesh.vertex_parent_dim.push_back(pdim);
    mesh.vertex_parent_id.push_back(pid);
    mesh.vertex_root_edge_id.push_back(-1);

    const int nv_new = mesh.n_vertices();
    mesh.vertex_zero_mask.resize(static_cast<std::size_t>(nv_new), 0);
    mesh.vertex_inside_mask.resize(static_cast<std::size_t>(nv_new), 0);
    mesh.vertex_phi.resize(static_cast<std::size_t>(nv_new * mesh.n_level_sets), T(0));

    // Rebuild cell arrays: keep non-incident cells, green-split incident ones
    std::vector<int32_t>    new_cell_vertices;
    std::vector<int32_t>    new_cell_offsets = {0};
    std::vector<cell::type> new_cell_types;

    const int old_nc = mesh.n_cells();
    std::vector<uint8_t> is_incident(static_cast<std::size_t>(old_nc), 0);
    for (const int c : incident_cells)
        is_incident[static_cast<std::size_t>(c)] = 1;

    for (int c = 0; c < old_nc; ++c)
    {
        const int        c0 = mesh.cell_offsets[static_cast<std::size_t>(c)];
        const int        c1 = mesh.cell_offsets[static_cast<std::size_t>(c + 1)];
        const cell::type ct = mesh.cell_types[static_cast<std::size_t>(c)];

        if (!is_incident[static_cast<std::size_t>(c)])
        {
            for (int j = c0; j < c1; ++j)
                new_cell_vertices.push_back(
                    mesh.cell_vertices[static_cast<std::size_t>(j)]);
            new_cell_offsets.push_back(static_cast<int32_t>(new_cell_vertices.size()));
            new_cell_types.push_back(ct);
            continue;
        }

        int n_children = 0;
        switch (ct)
        {
        case cell::type::interval:
        {
            n_children = green_split_interval(
                mesh.cell_vertices[static_cast<std::size_t>(c0)],
                mesh.cell_vertices[static_cast<std::size_t>(c0 + 1)],
                ev0, ev1, mid_v,
                new_cell_vertices, new_cell_offsets, new_cell_types);
            break;
        }
        case cell::type::triangle:
        {
            n_children = green_split_triangle(
                mesh.cell_vertices[static_cast<std::size_t>(c0)],
                mesh.cell_vertices[static_cast<std::size_t>(c0 + 1)],
                mesh.cell_vertices[static_cast<std::size_t>(c0 + 2)],
                ev0, ev1, mid_v,
                new_cell_vertices, new_cell_offsets, new_cell_types);
            break;
        }
        case cell::type::quadrilateral:
        {
            n_children = green_split_quadrilateral(
                mesh.cell_vertices[static_cast<std::size_t>(c0)],
                mesh.cell_vertices[static_cast<std::size_t>(c0 + 1)],
                mesh.cell_vertices[static_cast<std::size_t>(c0 + 2)],
                mesh.cell_vertices[static_cast<std::size_t>(c0 + 3)],
                ev0, ev1, mid_v,
                new_cell_vertices, new_cell_offsets, new_cell_types);
            break;
        }
        case cell::type::tetrahedron:
        {
            n_children = green_split_tetrahedron(
                mesh.cell_vertices[static_cast<std::size_t>(c0)],
                mesh.cell_vertices[static_cast<std::size_t>(c0 + 1)],
                mesh.cell_vertices[static_cast<std::size_t>(c0 + 2)],
                mesh.cell_vertices[static_cast<std::size_t>(c0 + 3)],
                ev0, ev1, mid_v,
                new_cell_vertices, new_cell_offsets, new_cell_types);
            break;
        }
        case cell::type::hexahedron:
        {
            const std::span<const int32_t> hex_verts(
                mesh.cell_vertices.data() + static_cast<std::size_t>(c0),
                static_cast<std::size_t>(c1 - c0));
            n_children = green_split_hexahedron(
                hex_verts, ev0, ev1, mid_v,
                new_cell_vertices, new_cell_offsets, new_cell_types);
            break;
        }
        default:
            // Prism / pyramid — should not reach here, checked above
            return false;
        }

        if (n_children == 0)
            return false;
    }

    // Commit
    mesh.cell_vertices = std::move(new_cell_vertices);
    mesh.cell_offsets  = std::move(new_cell_offsets);
    mesh.cell_types    = std::move(new_cell_types);
    mesh.cell_domain.assign(
        mesh.cell_types.size(), static_cast<uint8_t>(cell::domain::unset));

    // Rebuild connectivity
    build_local_edges(mesh);
    build_local_faces(mesh);
    rebuild_parent_entity_maps(mesh);

    return true;
}

template <std::floating_point T>
T find_separator_on_edge_bernstein(std::span<const T> edge_coeffs, T tol)
{
    const int p = static_cast<int>(edge_coeffs.size()) - 1;
    if (p < 1)
        return T(0.5);

    struct SubInterval
    {
        std::vector<T> coeffs;
        T t_lo;
        T t_hi;
    };

    std::vector<SubInterval> work;
    work.push_back(
        {std::vector<T>(edge_coeffs.begin(), edge_coeffs.end()), T(0), T(1)});

    for (int depth = 0; depth < 20; ++depth)
    {
        std::vector<SubInterval> next;
        next.reserve(work.size() * 2);

        for (auto& seg : work)
        {
            const T           t_mid = (seg.t_lo + seg.t_hi) / T(2);
            std::vector<T>    left, right;
            de_casteljau_split_1d<T>(
                std::span<const T>(seg.coeffs.data(), seg.coeffs.size()),
                left, right, T(0.5));
            next.push_back({std::move(left),  seg.t_lo, t_mid});
            next.push_back({std::move(right), t_mid,    seg.t_hi});
        }

        struct Tagged
        {
            int  index;
            bool has_root;
            bool strict_one_sign;
        };
        std::vector<Tagged> tags(next.size());
        for (std::size_t i = 0; i < next.size(); ++i)
        {
            const auto& c  = next[i].coeffs;
            const bool  os = bernstein_all_strict_one_sign<T>(
                std::span<const T>(c.data(), c.size()), tol);
            const int sv = bernstein_sign_variation_count<T>(
                std::span<const T>(c.data(), c.size()), tol);
            tags[i] = {static_cast<int>(i), sv > 0, os};
        }

        // A strict-one-sign interval between two root intervals is the separator
        for (std::size_t i = 0; i + 2 < tags.size(); ++i)
        {
            if (tags[i].has_root && tags[i + 1].strict_one_sign && tags[i + 2].has_root)
                return (next[i + 1].t_lo + next[i + 1].t_hi) / T(2);
        }

        work = std::move(next);
    }

    return T(0.5); // fallback
}

template <std::floating_point T>
T find_separator_on_edge_midpoint(
    const LocalMesh<T>& mesh, int edge_id, int level_set_id, T tol)
{
    const int32_t v0 = mesh.edge_vertices[static_cast<std::size_t>(2 * edge_id)];
    const int32_t v1 = mesh.edge_vertices[static_cast<std::size_t>(2 * edge_id + 1)];
    const T phi0 = mesh.vertex_phi[static_cast<std::size_t>(
        v0 * mesh.n_level_sets + level_set_id)];
    const T phi1 = mesh.vertex_phi[static_cast<std::size_t>(
        v1 * mesh.n_level_sets + level_set_id)];

    const T phi_mid = phi0 + T(0.5) * (phi1 - phi0);
    if (std::abs(phi_mid) > tol)
        return T(0.5);

    for (const T t : {T(0.25), T(0.75), T(0.375), T(0.625)})
    {
        const T phi_t = phi0 + t * (phi1 - phi0);
        if (std::abs(phi_t) > tol)
            return t;
    }

    return T(0.5);
}

// Explicit instantiations
template std::vector<int> cells_incident_on_edge<double>(
    const LocalMesh<double>&, int);
template std::vector<int> cells_incident_on_edge<float>(
    const LocalMesh<float>&, int);

template bool green_split_one_edge<double>(LocalMesh<double>&, int, double);
template bool green_split_one_edge<float>(LocalMesh<float>&, int, float);

template double find_separator_on_edge_bernstein<double>(
    std::span<const double>, double);
template float find_separator_on_edge_bernstein<float>(
    std::span<const float>, float);

template double find_separator_on_edge_midpoint<double>(
    const LocalMesh<double>&, int, int, double);
template float find_separator_on_edge_midpoint<float>(
    const LocalMesh<float>&, int, int, float);

} // namespace cutcells
