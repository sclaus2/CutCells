// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "interface_split.h"
#include "edge_classification.h"
#include "edge_root.h"
#include "macro_split.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <unordered_map>

namespace cutcells
{

// ============================================================================
// Internal helpers
// ============================================================================

namespace
{

/// Pack two vertex IDs into a single 64-bit key (order-independent).
uint64_t iface_edge_key(int32_t a, int32_t b)
{
    const auto lo = static_cast<uint64_t>(std::min(a, b));
    const auto hi = static_cast<uint64_t>(std::max(a, b));
    return (lo << 32) | hi;
}

template <std::floating_point T>
int count_all_zero_tets_in_range(
    const LocalMesh<T>& mesh,
    int                 cell_begin,
    int                 cell_end,
    int                 level_set_id,
    T                   tol)
{
    int count = 0;
    const int n_ls = mesh.n_level_sets;
    for (int ci = std::max(0, cell_begin); ci < std::min(cell_end, mesh.n_cells()); ++ci)
    {
        if (mesh.cell_types[static_cast<std::size_t>(ci)] != cell::type::tetrahedron)
            continue;
        const int c0 = mesh.cell_offsets[static_cast<std::size_t>(ci)];
        const int c1 = mesh.cell_offsets[static_cast<std::size_t>(ci + 1)];
        if (c1 - c0 != 4)
            continue;

        bool all_zero = true;
        for (int j = c0; j < c1; ++j)
        {
            const int32_t vid = mesh.cell_vertices[static_cast<std::size_t>(j)];
            const T phi = mesh.vertex_phi[static_cast<std::size_t>(vid * n_ls + level_set_id)];
            if (std::abs(phi) > tol)
            {
                all_zero = false;
                break;
            }
        }
        if (all_zero)
            ++count;
    }
    return count;
}

template <std::floating_point T>
int find_local_edge_id(const LocalMesh<T>& mesh, int32_t va, int32_t vb)
{
    for (int e = 0; e < mesh.n_edges(); ++e)
    {
        const int32_t ev0 = mesh.edge_vertices[static_cast<std::size_t>(2 * e)];
        const int32_t ev1 = mesh.edge_vertices[static_cast<std::size_t>(2 * e + 1)];
        if ((ev0 == va && ev1 == vb) || (ev0 == vb && ev1 == va))
            return e;
    }
    return -1;
}

template <std::floating_point T>
std::pair<int32_t, int32_t> infer_edge_root_parent(
    const LocalMesh<T>& mesh,
    int32_t             va,
    int32_t             vb)
{
    const int edge_id = find_local_edge_id(mesh, va, vb);
    if (edge_id < 0)
        return {-1, -1};
    return {
        mesh.edge_parent_dim[static_cast<std::size_t>(edge_id)],
        mesh.edge_parent_id[static_cast<std::size_t>(edge_id)]};
}

template <std::floating_point T>
bool vertex_lies_on_parent_tet_face(
    const LocalMesh<T>& mesh,
    int32_t             lv,
    int                 bg_face_id)
{
    if (mesh.parent_cell_type != cell::type::tetrahedron)
        return false;
    if (bg_face_id < 0 || bg_face_id >= cell::num_faces(cell::type::tetrahedron))
        return false;

    const int pdim = mesh.vertex_parent_dim[static_cast<std::size_t>(lv)];
    const int pid = mesh.vertex_parent_id[static_cast<std::size_t>(lv)];
    if (pdim == 2 && pid == bg_face_id)
        return true;

    const auto bg_faces = cell::faces(cell::type::tetrahedron);
    const auto bg_edges = cell::edges(cell::type::tetrahedron);
    const auto face = bg_faces[static_cast<std::size_t>(bg_face_id)];

    if (pdim == 1 && pid >= 0 && pid < static_cast<int>(bg_edges.size()))
    {
        const auto edge = bg_edges[static_cast<std::size_t>(pid)];
        bool v0_on_face = false;
        bool v1_on_face = false;
        for (const int fv : face)
        {
            if (edge[0] == fv)
                v0_on_face = true;
            if (edge[1] == fv)
                v1_on_face = true;
        }
        return v0_on_face && v1_on_face;
    }

    if (pdim == 0)
    {
        for (const int fv : face)
        {
            if (pid == fv)
                return true;
        }
    }

    return false;
}

template <std::floating_point T>
int infer_parent_tet_face_id(
    const LocalMesh<T>&      mesh,
    std::span<const int32_t> face_vertices)
{
    if (mesh.parent_cell_type != cell::type::tetrahedron)
        return -1;

    for (int f = 0; f < cell::num_faces(cell::type::tetrahedron); ++f)
    {
        bool all_on_face = true;
        for (const int32_t lv : face_vertices)
        {
            if (!vertex_lies_on_parent_tet_face(mesh, lv, f))
            {
                all_on_face = false;
                break;
            }
        }
        if (all_on_face)
            return f;
    }
    return -1;
}

/// Create a new vertex with explicit parent ancestry and phi value.
template <std::floating_point T>
int32_t create_vertex(
    LocalMesh<T>& mesh,
    std::span<const T> phys_coords,
    std::span<const T> ref_coords,
    int parent_dim,
    int parent_id,
    int root_edge_id,
    int level_set_id,
    T phi_value)
{
    const int n_ls = mesh.n_level_sets;
    const int gdim = mesh.gdim;
    const int tdim = mesh.tdim;

    const auto new_vid = static_cast<int32_t>(mesh.vertex_parent_dim.size());
    for (int d = 0; d < gdim; ++d)
        mesh.vertex_x.push_back(phys_coords[static_cast<std::size_t>(d)]);
    for (int d = 0; d < tdim; ++d)
        mesh.vertex_ref_x.push_back(ref_coords[static_cast<std::size_t>(d)]);

    mesh.vertex_parent_dim.push_back(parent_dim);
    mesh.vertex_parent_id.push_back(parent_id);
    mesh.vertex_root_edge_id.push_back(root_edge_id);

    const auto nv = static_cast<std::size_t>(new_vid + 1);
    mesh.vertex_zero_mask.resize(nv, 0);
    mesh.vertex_inside_mask.resize(nv, 0);
    mesh.vertex_phi.resize(nv * static_cast<std::size_t>(n_ls), T(0));
    mesh.vertex_phi[static_cast<std::size_t>(new_vid * n_ls + level_set_id)] = phi_value;

    const uint64_t mask = (uint64_t(1) << static_cast<unsigned>(level_set_id));
    if (phi_value == T(0))
        mesh.vertex_zero_mask[static_cast<std::size_t>(new_vid)] |= mask;
    else if (phi_value < T(0))
        mesh.vertex_inside_mask[static_cast<std::size_t>(new_vid)] |= mask;

    return new_vid;
}

/// Create a new vertex on the interface (phi = 0) with explicit ancestry.
template <std::floating_point T>
int32_t create_root_vertex(
    LocalMesh<T>& mesh,
    std::span<const T> phys_coords,
    std::span<const T> ref_coords,
    int parent_dim,
    int parent_id,
    int level_set_id)
{
    return create_vertex(
        mesh, phys_coords, ref_coords, parent_dim, parent_id, -1, level_set_id, T(0));
}

/// Find root on edge (va, vb) via ITP solver and create a new vertex at
/// the root location.  If root_cache is non-null, deduplicate across cells.
///
/// @return new vertex ID, or -1 if root-finding failed
template <std::floating_point T>
int32_t find_and_create_edge_root(
    LocalMesh<T>& mesh,
    int32_t va,
    int32_t vb,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol,
    std::unordered_map<uint64_t, int32_t>* root_cache)
{
    // Check cache first
    if (root_cache)
    {
        const uint64_t key = iface_edge_key(va, vb);
        auto it = root_cache->find(key);
        if (it != root_cache->end())
            return it->second;
    }

    const int gdim = mesh.gdim;
    const int tdim = mesh.tdim;
    const int n_ls = mesh.n_level_sets;

    // Get phi values at endpoints
    const T phi_a = mesh.vertex_phi[static_cast<std::size_t>(
        va * n_ls + level_set_id)];
    const T phi_b = mesh.vertex_phi[static_cast<std::size_t>(
        vb * n_ls + level_set_id)];

    // Must have a sign change
    if (phi_a * phi_b > T(0))
        return -1;

    // Endpoint coordinates
    std::array<T, 3> pa = {T(0), T(0), T(0)};
    std::array<T, 3> pb = {T(0), T(0), T(0)};
    std::array<T, 3> ra = {T(0), T(0), T(0)};
    std::array<T, 3> rb = {T(0), T(0), T(0)};

    for (int d = 0; d < gdim; ++d)
    {
        pa[static_cast<std::size_t>(d)]
            = mesh.vertex_x[static_cast<std::size_t>(va * gdim + d)];
        pb[static_cast<std::size_t>(d)]
            = mesh.vertex_x[static_cast<std::size_t>(vb * gdim + d)];
    }
    for (int d = 0; d < tdim; ++d)
    {
        ra[static_cast<std::size_t>(d)]
            = mesh.vertex_ref_x[static_cast<std::size_t>(va * tdim + d)];
        rb[static_cast<std::size_t>(d)]
            = mesh.vertex_ref_x[static_cast<std::size_t>(vb * tdim + d)];
    }

    // Parametric evaluation along edge [0,1]
    std::array<T, 3> x_eval = {T(0), T(0), T(0)};
    auto eval_edge = [&](T t) -> T
    {
        for (int d = 0; d < gdim; ++d)
            x_eval[static_cast<std::size_t>(d)]
                = pa[static_cast<std::size_t>(d)]
                  + t * (pb[static_cast<std::size_t>(d)]
                         - pa[static_cast<std::size_t>(d)]);
        return ls.value(x_eval.data(), mesh.parent_cell_id);
    };

    const T t_linear
        = cell::edge_root::linear_root_parameter(phi_a, phi_b);
    int iterations = 0;
    bool converged = false;
    const T t_root = cell::edge_root::itp_parameter<T>(
        eval_edge, T(0), T(1), phi_a, phi_b,
        t_linear, /*max_iter=*/100, /*xtol=*/tol, /*ftol=*/tol,
        &iterations, &converged);

    const T residual = std::abs(eval_edge(t_root));
    if (!converged && residual > T(1e-8))
        return -1;
    if (!converged)
    {
        std::cerr << "Debug: interface edge-root solve accepted without solver convergence"
                  << " on parent cell " << mesh.parent_cell_id
                  << ", residual=" << residual
                  << ", edge=(" << va << "," << vb << ")\n";
    }

    // Compute physical and reference coords at root
    std::array<T, 3> p_root = {T(0), T(0), T(0)};
    std::array<T, 3> r_root = {T(0), T(0), T(0)};
    for (int d = 0; d < gdim; ++d)
        p_root[static_cast<std::size_t>(d)]
            = pa[static_cast<std::size_t>(d)]
              + t_root * (pb[static_cast<std::size_t>(d)]
                          - pa[static_cast<std::size_t>(d)]);
    for (int d = 0; d < tdim; ++d)
        r_root[static_cast<std::size_t>(d)]
            = ra[static_cast<std::size_t>(d)]
              + t_root * (rb[static_cast<std::size_t>(d)]
                          - ra[static_cast<std::size_t>(d)]);

    const auto [parent_dim, parent_id] = infer_edge_root_parent(mesh, va, vb);

    // Create vertex on the interface
    const int32_t vid = create_root_vertex(
        mesh,
        std::span<const T>(p_root.data(), static_cast<std::size_t>(gdim)),
        std::span<const T>(r_root.data(), static_cast<std::size_t>(tdim)),
        parent_dim,
        parent_id,
        level_set_id);

    // Store in cache
    if (root_cache)
        (*root_cache)[iface_edge_key(va, vb)] = vid;

    return vid;
}

/// Create an interior vertex on the interface via root-finding on the
/// ray from the apex through the centroid of the 3 edge root vertices.
///
/// The centroid of 3 edge roots (which are on the interface) is used as
/// a target direction from the apex.  The root on this ray is found by
/// ITP solving, placing the new vertex exactly on the interface.
///
/// @return new vertex ID, or -1 if root-finding fails
template <std::floating_point T>
int32_t create_root_centroid(
    LocalMesh<T>& mesh,
    int32_t apex_vid,
    int32_t r0,
    int32_t r1,
    int32_t r2,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol)
{
    const int gdim = mesh.gdim;
    const int tdim = mesh.tdim;
    const int n_ls = mesh.n_level_sets;
    const T third = T(1) / T(3);

    // Apex coords
    std::array<T, 3> apex_x = {T(0), T(0), T(0)};
    std::array<T, 3> apex_ref = {T(0), T(0), T(0)};
    for (int d = 0; d < gdim; ++d)
        apex_x[static_cast<std::size_t>(d)]
            = mesh.vertex_x[static_cast<std::size_t>(apex_vid * gdim + d)];
    for (int d = 0; d < tdim; ++d)
        apex_ref[static_cast<std::size_t>(d)]
            = mesh.vertex_ref_x[static_cast<std::size_t>(apex_vid * tdim + d)];

    // Root centroid (target point on/near the interface)
    std::array<T, 3> rc_x = {T(0), T(0), T(0)};
    std::array<T, 3> rc_ref = {T(0), T(0), T(0)};
    for (int d = 0; d < gdim; ++d)
    {
        rc_x[static_cast<std::size_t>(d)] = third * (
            mesh.vertex_x[static_cast<std::size_t>(r0 * gdim + d)]
          + mesh.vertex_x[static_cast<std::size_t>(r1 * gdim + d)]
          + mesh.vertex_x[static_cast<std::size_t>(r2 * gdim + d)]);
    }
    for (int d = 0; d < tdim; ++d)
    {
        rc_ref[static_cast<std::size_t>(d)] = third * (
            mesh.vertex_ref_x[static_cast<std::size_t>(r0 * tdim + d)]
          + mesh.vertex_ref_x[static_cast<std::size_t>(r1 * tdim + d)]
          + mesh.vertex_ref_x[static_cast<std::size_t>(r2 * tdim + d)]);
    }

    // Evaluate phi at apex and at root centroid
    const T phi_apex = mesh.vertex_phi[static_cast<std::size_t>(
        apex_vid * n_ls + level_set_id)];

    std::array<T, 3> x_tmp = rc_x;
    const T phi_rc = ls.value(x_tmp.data(), mesh.parent_cell_id);

    // Check for sign change on ray [apex, root_centroid]
    // If same sign: the root centroid is on the same side as apex.
    // Extend the ray further beyond the root centroid to find a bracket.
    T t_bracket_end = T(1);
    T phi_bracket_end = phi_rc;

    if (phi_apex * phi_rc >= T(0))
    {
        // Extend the ray: try t = 2, 4, 8 along the ray from apex
        // through root centroid
        bool found = false;
        for (T ext = T(2); ext <= T(8); ext *= T(2))
        {
            std::array<T, 3> x_ext = {T(0), T(0), T(0)};
            for (int d = 0; d < gdim; ++d)
                x_ext[static_cast<std::size_t>(d)]
                    = apex_x[static_cast<std::size_t>(d)]
                      + ext * (rc_x[static_cast<std::size_t>(d)]
                               - apex_x[static_cast<std::size_t>(d)]);
            const T phi_ext = ls.value(x_ext.data(), mesh.parent_cell_id);
            if (phi_apex * phi_ext < T(0))
            {
                t_bracket_end = ext;
                phi_bracket_end = phi_ext;
                found = true;
                break;
            }
        }
        if (!found)
            return -1; // cannot bracket
    }

    // ITP root solving on [0, t_bracket_end]
    std::array<T, 3> x_eval = {T(0), T(0), T(0)};
    auto eval_ray = [&](T t) -> T
    {
        for (int d = 0; d < gdim; ++d)
            x_eval[static_cast<std::size_t>(d)]
                = apex_x[static_cast<std::size_t>(d)]
                  + t * (rc_x[static_cast<std::size_t>(d)]
                         - apex_x[static_cast<std::size_t>(d)]);
        return ls.value(x_eval.data(), mesh.parent_cell_id);
    };

    const T t_linear = cell::edge_root::linear_root_parameter(
        phi_apex, phi_bracket_end);
    const T t_guess = t_linear * t_bracket_end;

    int iterations = 0;
    bool converged = false;
    const T t_root = cell::edge_root::brent_solve<T>(
        eval_ray, T(0), t_bracket_end, phi_apex, phi_bracket_end,
        /*max_iter=*/100, /*xtol=*/tol, /*ftol=*/tol,
        &iterations, &converged);

    const T residual = std::abs(eval_ray(t_root));
    if (!converged && residual > T(1e-8))
        return -1;
    if (!converged)
    {
        std::cerr << "Debug: interface centroid ray solve accepted without solver convergence"
                  << " on parent cell " << mesh.parent_cell_id
                  << ", residual=" << residual
                  << ", apex_vid=" << apex_vid << "\n";
    }

    // Compute the root point coordinates
    std::array<T, 3> p_root = {T(0), T(0), T(0)};
    std::array<T, 3> r_root = {T(0), T(0), T(0)};
    for (int d = 0; d < gdim; ++d)
        p_root[static_cast<std::size_t>(d)]
            = apex_x[static_cast<std::size_t>(d)]
              + t_root * (rc_x[static_cast<std::size_t>(d)]
                          - apex_x[static_cast<std::size_t>(d)]);
    for (int d = 0; d < tdim; ++d)
        r_root[static_cast<std::size_t>(d)]
            = apex_ref[static_cast<std::size_t>(d)]
              + t_root * (rc_ref[static_cast<std::size_t>(d)]
                          - apex_ref[static_cast<std::size_t>(d)]);

    // Create vertex exactly on the interface
    return create_root_vertex(
        mesh,
        std::span<const T>(p_root.data(), static_cast<std::size_t>(gdim)),
        std::span<const T>(r_root.data(), static_cast<std::size_t>(tdim)),
        3, // parent_dim = 3 (cell interior)
        mesh.parent_cell_id,
        level_set_id);
}

template <std::floating_point T>
T signed_tet_volume_from_vertices(
    const LocalMesh<T>&           mesh,
    const std::array<int32_t, 4>& tet)
{
    if (mesh.gdim < 3)
        return T(0);

    auto load = [&](int32_t vid, int axis) -> T
    {
        return mesh.vertex_x[static_cast<std::size_t>(vid * mesh.gdim + axis)];
    };

    const T ax = load(tet[0], 0);
    const T ay = load(tet[0], 1);
    const T az = load(tet[0], 2);
    const T bx = load(tet[1], 0) - ax;
    const T by = load(tet[1], 1) - ay;
    const T bz = load(tet[1], 2) - az;
    const T cx = load(tet[2], 0) - ax;
    const T cy = load(tet[2], 1) - ay;
    const T cz = load(tet[2], 2) - az;
    const T dx = load(tet[3], 0) - ax;
    const T dy = load(tet[3], 1) - ay;
    const T dz = load(tet[3], 2) - az;

    return (bx * (cy * dz - cz * dy)
          - by * (cx * dz - cz * dx)
          + bz * (cx * dy - cy * dx))
           / T(6);
}

template <std::floating_point T>
T signed_tet_volume_from_cell(
    const LocalMesh<T>& mesh,
    int                 cell_id)
{
    const int c0 = mesh.cell_offsets[static_cast<std::size_t>(cell_id)];
    const int c1 = mesh.cell_offsets[static_cast<std::size_t>(cell_id + 1)];
    if (c1 - c0 != 4)
        return T(0);

    const std::array<int32_t, 4> tet = {
        mesh.cell_vertices[static_cast<std::size_t>(c0)],
        mesh.cell_vertices[static_cast<std::size_t>(c0 + 1)],
        mesh.cell_vertices[static_cast<std::size_t>(c0 + 2)],
        mesh.cell_vertices[static_cast<std::size_t>(c0 + 3)]};
    return signed_tet_volume_from_vertices(mesh, tet);
}

template <std::floating_point T>
void orient_child_tets_like_parent(
    const LocalMesh<T>&                  mesh,
    int                                  parent_cell_id,
    std::span<std::array<int32_t, 4>>    children)
{
    if (mesh.gdim < 3)
        return;

    const T parent_vol = signed_tet_volume_from_cell(mesh, parent_cell_id);
    if (parent_vol == T(0))
        return;

    for (auto& child : children)
    {
        const T child_vol = signed_tet_volume_from_vertices(mesh, child);
        if (child_vol == T(0))
            continue;
        if ((child_vol > T(0)) != (parent_vol > T(0)))
            std::swap(child[2], child[3]);
    }
}

/// Core implementation of the interface-adapted split.
template <std::floating_point T>
bool interface_split_impl(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol,
    std::unordered_map<uint64_t, int32_t>* root_cache);

/// Replace a single cell with N child tets in the mesh.
///
/// This is shared by the 1-vs-3 (10 children) and 2-vs-2 (12 children)
/// implementations.
template <std::floating_point T>
void replace_cell_with_children(
    LocalMesh<T>&                           mesh,
    int                                     cell_id,
    std::span<const std::array<int32_t, 4>> children);

template <std::floating_point T, std::size_t N>
void replace_cell_with_children(
    LocalMesh<T>& mesh,
    int cell_id,
    const std::array<std::array<int32_t, 4>, N>& children)
{
    replace_cell_with_children(
        mesh,
        cell_id,
        std::span<const std::array<int32_t, 4>>(children.data(), children.size()));
}

/// Replace a single cell with a dynamic list of child tets in the mesh.
template <std::floating_point T>
void replace_cell_with_children(
    LocalMesh<T>&                              mesh,
    int                                        cell_id,
    std::span<const std::array<int32_t, 4>>    children)
{
    const int old_nc = mesh.n_cells();
    const std::size_t n_children = children.size();

    std::vector<int32_t>    new_cell_vertices;
    std::vector<int32_t>    new_cell_offsets = {0};
    std::vector<cell::type> new_cell_types;

    new_cell_vertices.reserve(
        mesh.cell_vertices.size() + n_children * 4);
    new_cell_offsets.reserve(
        static_cast<std::size_t>(old_nc) + n_children);
    new_cell_types.reserve(
        static_cast<std::size_t>(old_nc) + n_children - 1);

    for (int ci = 0; ci < old_nc; ++ci)
    {
        const int s0 = mesh.cell_offsets[static_cast<std::size_t>(ci)];
        const int s1 = mesh.cell_offsets[static_cast<std::size_t>(ci + 1)];

        if (ci != cell_id)
        {
            for (int j = s0; j < s1; ++j)
                new_cell_vertices.push_back(
                    mesh.cell_vertices[static_cast<std::size_t>(j)]);
            new_cell_offsets.push_back(
                static_cast<int32_t>(new_cell_vertices.size()));
            new_cell_types.push_back(
                mesh.cell_types[static_cast<std::size_t>(ci)]);
        }
        else
        {
            for (const auto& child : children)
            {
                for (int32_t vid : child)
                    new_cell_vertices.push_back(vid);
                new_cell_offsets.push_back(
                    static_cast<int32_t>(new_cell_vertices.size()));
                new_cell_types.push_back(cell::type::tetrahedron);
            }
        }
    }

    mesh.cell_vertices = std::move(new_cell_vertices);
    mesh.cell_offsets = std::move(new_cell_offsets);
    mesh.cell_types = std::move(new_cell_types);
    mesh.cell_domain.assign(
        mesh.cell_types.size(),
        static_cast<uint8_t>(cell::domain::unset));

    build_local_edges(mesh);
    build_local_faces(mesh);
    rebuild_parent_entity_maps(mesh);
}

template <std::floating_point T>
bool certify_split_children(
    LocalMesh<T>&                mesh,
    const LevelSetFunction<T>&   ls,
    int                          first_child,
    int                          child_count,
    int                          level_set_id,
    T                            tol)
{
    if (!ls.has_value())
        return false;

    classify_edges_and_mark_refine(
        mesh, ls, level_set_id, tol, false);

    const int last_child = first_child + child_count;
    for (int c = first_child; c < last_child; ++c)
    {
        const int ce0 = mesh.cell_edge_offsets[static_cast<std::size_t>(c)];
        const int ce1 = mesh.cell_edge_offsets[static_cast<std::size_t>(c + 1)];
        for (int idx = ce0; idx < ce1; ++idx)
        {
            const int edge_id = mesh.cell_edges_flat[static_cast<std::size_t>(idx)];
            const auto state = static_cast<EdgeState>(
                mesh.edge_state_for(edge_id, level_set_id));
            if (state == EdgeState::single_cross
                || state == EdgeState::multi_cross
                || state == EdgeState::uncertain)
            {
                return false;
            }
        }
    }

    return true;
}

/// Create an interior vertex on the interface via root-finding on a
/// ray from a source vertex through the centroid of N edge root vertices.
///
/// @return new vertex ID, or -1 if root-finding fails
template <std::floating_point T>
int32_t create_root_centroid_n(
    LocalMesh<T>& mesh,
    int32_t source_vid,
    std::span<const int32_t> root_vids,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol)
{
    const int gdim = mesh.gdim;
    const int tdim = mesh.tdim;
    const int n_ls = mesh.n_level_sets;
    const int n_roots = static_cast<int>(root_vids.size());
    const T inv_n = T(1) / static_cast<T>(n_roots);

    // Source vertex coords
    std::array<T, 3> src_x = {T(0), T(0), T(0)};
    std::array<T, 3> src_ref = {T(0), T(0), T(0)};
    for (int d = 0; d < gdim; ++d)
        src_x[static_cast<std::size_t>(d)]
            = mesh.vertex_x[static_cast<std::size_t>(source_vid * gdim + d)];
    for (int d = 0; d < tdim; ++d)
        src_ref[static_cast<std::size_t>(d)]
            = mesh.vertex_ref_x[static_cast<std::size_t>(source_vid * tdim + d)];

    // Root centroid
    std::array<T, 3> rc_x = {T(0), T(0), T(0)};
    std::array<T, 3> rc_ref = {T(0), T(0), T(0)};
    for (int ri = 0; ri < n_roots; ++ri)
    {
        const int32_t rv = root_vids[static_cast<std::size_t>(ri)];
        for (int d = 0; d < gdim; ++d)
            rc_x[static_cast<std::size_t>(d)]
                += mesh.vertex_x[static_cast<std::size_t>(rv * gdim + d)];
        for (int d = 0; d < tdim; ++d)
            rc_ref[static_cast<std::size_t>(d)]
                += mesh.vertex_ref_x[static_cast<std::size_t>(rv * tdim + d)];
    }
    for (int d = 0; d < gdim; ++d)
        rc_x[static_cast<std::size_t>(d)] *= inv_n;
    for (int d = 0; d < tdim; ++d)
        rc_ref[static_cast<std::size_t>(d)] *= inv_n;

    // Evaluate phi at source and root centroid
    const T phi_src = mesh.vertex_phi[static_cast<std::size_t>(
        source_vid * n_ls + level_set_id)];

    std::array<T, 3> x_tmp = rc_x;
    const T phi_rc = ls.value(x_tmp.data(), mesh.parent_cell_id);

    // Check for sign change on ray [source, root_centroid]
    T t_bracket_end = T(1);
    T phi_bracket_end = phi_rc;

    if (phi_src * phi_rc >= T(0))
    {
        // Extend the ray beyond the root centroid
        bool found = false;
        for (T ext = T(2); ext <= T(8); ext *= T(2))
        {
            std::array<T, 3> x_ext = {T(0), T(0), T(0)};
            for (int d = 0; d < gdim; ++d)
                x_ext[static_cast<std::size_t>(d)]
                    = src_x[static_cast<std::size_t>(d)]
                      + ext * (rc_x[static_cast<std::size_t>(d)]
                               - src_x[static_cast<std::size_t>(d)]);
            const T phi_ext = ls.value(x_ext.data(), mesh.parent_cell_id);
            if (phi_src * phi_ext < T(0))
            {
                t_bracket_end = ext;
                phi_bracket_end = phi_ext;
                found = true;
                break;
            }
        }
        if (!found)
            return -1;
    }

    // Brent solver on [0, t_bracket_end]
    std::array<T, 3> x_eval = {T(0), T(0), T(0)};
    auto eval_ray = [&](T t) -> T
    {
        for (int d = 0; d < gdim; ++d)
            x_eval[static_cast<std::size_t>(d)]
                = src_x[static_cast<std::size_t>(d)]
                  + t * (rc_x[static_cast<std::size_t>(d)]
                         - src_x[static_cast<std::size_t>(d)]);
        return ls.value(x_eval.data(), mesh.parent_cell_id);
    };

    int iterations = 0;
    bool converged = false;
    const T t_root = cell::edge_root::brent_solve<T>(
        eval_ray, T(0), t_bracket_end, phi_src, phi_bracket_end,
        /*max_iter=*/100, /*xtol=*/tol, /*ftol=*/tol,
        &iterations, &converged);

    const T residual = std::abs(eval_ray(t_root));
    if (!converged && residual > T(1e-8))
        return -1;
    if (!converged)
    {
        std::cerr << "Debug: interface centroid-n ray solve accepted without solver convergence"
                  << " on parent cell " << mesh.parent_cell_id
                  << ", residual=" << residual
                  << ", source_vid=" << source_vid << "\n";
    }

    // Compute root point coordinates
    std::array<T, 3> p_root = {T(0), T(0), T(0)};
    std::array<T, 3> r_root = {T(0), T(0), T(0)};
    for (int d = 0; d < gdim; ++d)
        p_root[static_cast<std::size_t>(d)]
            = src_x[static_cast<std::size_t>(d)]
              + t_root * (rc_x[static_cast<std::size_t>(d)]
                          - src_x[static_cast<std::size_t>(d)]);
    for (int d = 0; d < tdim; ++d)
        r_root[static_cast<std::size_t>(d)]
            = src_ref[static_cast<std::size_t>(d)]
              + t_root * (rc_ref[static_cast<std::size_t>(d)]
                          - src_ref[static_cast<std::size_t>(d)]);

    return create_root_vertex(
        mesh,
        std::span<const T>(p_root.data(), static_cast<std::size_t>(gdim)),
        std::span<const T>(r_root.data(), static_cast<std::size_t>(tdim)),
        3, // parent_dim = 3 (cell interior)
        mesh.parent_cell_id,
        level_set_id);
}

/// Pick the positive/negative pairs for a 2-vs-2 tetrahedron.
///
/// @param phi  level-set values at the 4 tet vertices
/// @return {p0, p1, n0, n1} as local indices, or {-1,-1,-1,-1} if not 2-vs-2
template <std::floating_point T>
std::array<int, 4> pick_topology2_pairs(std::span<const T> phi)
{
    assert(phi.size() == 4);
    int pos[2] = {-1, -1};
    int neg[2] = {-1, -1};
    int np = 0, nn = 0;
    for (int i = 0; i < 4; ++i)
    {
        if (phi[static_cast<std::size_t>(i)] > T(0))
        {
            if (np >= 2) return {-1, -1, -1, -1};
            pos[np++] = i;
        }
        else if (phi[static_cast<std::size_t>(i)] < T(0))
        {
            if (nn >= 2) return {-1, -1, -1, -1};
            neg[nn++] = i;
        }
    }
    if (np != 2 || nn != 2)
        return {-1, -1, -1, -1};
    return {pos[0], pos[1], neg[0], neg[1]};
}

/// Core implementation of the 2-vs-2 interface-adapted split.
template <std::floating_point T>
bool interface_split_impl2(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol,
    std::unordered_map<uint64_t, int32_t>* root_cache)
{
    if (!ls.has_value())
        return false;

    const LocalMesh<T> mesh_before = mesh;
    std::unordered_map<uint64_t, int32_t> root_cache_before;
    if (root_cache)
        root_cache_before = *root_cache;
    const auto rollback = [&]() -> bool
    {
        mesh = mesh_before;
        if (root_cache)
            *root_cache = root_cache_before;
        return false;
    };

    const int n_ls = mesh.n_level_sets;

    const int c0 = mesh.cell_offsets[static_cast<std::size_t>(cell_id)];
    const int c1 = mesh.cell_offsets[static_cast<std::size_t>(cell_id + 1)];
    if (c1 - c0 != 4)
        return rollback();

    // Read tet vertex IDs and phi values
    std::array<int32_t, 4> vids;
    std::array<T, 4> phi_vals;
    for (int j = 0; j < 4; ++j)
    {
        vids[static_cast<std::size_t>(j)]
            = mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];
        phi_vals[static_cast<std::size_t>(j)]
            = mesh.vertex_phi[static_cast<std::size_t>(
                  vids[static_cast<std::size_t>(j)] * n_ls + level_set_id)];
    }

    // Identify the 2-vs-2 pairs
    const auto pairs = pick_topology2_pairs<T>(
        std::span<const T>(phi_vals.data(), 4));
    if (pairs[0] < 0)
        return rollback();

    const int32_t p0 = vids[static_cast<std::size_t>(pairs[0])];
    const int32_t p1 = vids[static_cast<std::size_t>(pairs[1])];
    const int32_t n0 = vids[static_cast<std::size_t>(pairs[2])];
    const int32_t n1 = vids[static_cast<std::size_t>(pairs[3])];

    // -----------------------------------------------------------------
    // Find 4 edge roots on crossing edges
    // -----------------------------------------------------------------
    const int32_t r_p0n0 = find_and_create_edge_root(
        mesh, p0, n0, ls, level_set_id, tol, root_cache);
    if (r_p0n0 < 0) return rollback();

    const int32_t r_p0n1 = find_and_create_edge_root(
        mesh, p0, n1, ls, level_set_id, tol, root_cache);
    if (r_p0n1 < 0) return rollback();

    const int32_t r_p1n0 = find_and_create_edge_root(
        mesh, p1, n0, ls, level_set_id, tol, root_cache);
    if (r_p1n0 < 0) return rollback();

    const int32_t r_p1n1 = find_and_create_edge_root(
        mesh, p1, n1, ls, level_set_id, tol, root_cache);
    if (r_p1n1 < 0) return rollback();

    // -----------------------------------------------------------------
    // Create interior vertex on interface via ray through centroid
    // of the 4 edge roots
    // -----------------------------------------------------------------
    // Pick source vertex: vertex with largest |phi|
    int source_local = 0;
    T max_abs_phi = std::abs(phi_vals[0]);
    for (int i = 1; i < 4; ++i)
    {
        const T a = std::abs(phi_vals[static_cast<std::size_t>(i)]);
        if (a > max_abs_phi)
        {
            source_local = i;
            max_abs_phi = a;
        }
    }
    const int32_t source_vid = vids[static_cast<std::size_t>(source_local)];

    std::array<int32_t, 4> root_arr = {r_p0n0, r_p0n1, r_p1n0, r_p1n1};
    const int32_t c = create_root_centroid_n(
        mesh, source_vid,
        std::span<const int32_t>(root_arr.data(), 4),
        ls, level_set_id, tol);
    if (c < 0)
        return rollback();

    // -----------------------------------------------------------------
    // Build 12 child tetrahedra
    //
    // All 4 faces are intersected (1-vs-2 pattern). Each face has 2
    // crossing edges with roots, giving 1 lone-side triangle and 1
    // quad (triangulated into 2) on the pair side.
    //
    // Cone from interior vertex c to each of 12 boundary sub-triangles:
    //
    //   Face (p0, p1, n0), roots r_p0n0 on (p0,n0), r_p1n0 on (p1,n0):
    //     C0  = (c, n0, r_p0n0, r_p1n0)     — negative lone
    //     C1  = (c, p0, p1, r_p1n0)          — positive quad part 1
    //     C2  = (c, p0, r_p1n0, r_p0n0)      — positive quad part 2
    //
    //   Face (p0, p1, n1), roots r_p0n1 on (p0,n1), r_p1n1 on (p1,n1):
    //     C3  = (c, n1, r_p0n1, r_p1n1)      — negative lone
    //     C4  = (c, p0, p1, r_p1n1)           — positive quad part 1
    //     C5  = (c, p0, r_p1n1, r_p0n1)       — positive quad part 2
    //
    //   Face (p0, n0, n1), roots r_p0n0 on (p0,n0), r_p0n1 on (p0,n1):
    //     C6  = (c, p0, r_p0n0, r_p0n1)       — positive lone
    //     C7  = (c, n0, n1, r_p0n1)            — negative quad part 1
    //     C8  = (c, n0, r_p0n1, r_p0n0)        — negative quad part 2
    //
    //   Face (p1, n0, n1), roots r_p1n0 on (p1,n0), r_p1n1 on (p1,n1):
    //     C9  = (c, p1, r_p1n0, r_p1n1)        — positive lone
    //     C10 = (c, n0, n1, r_p1n1)             — negative quad part 1
    //     C11 = (c, n0, r_p1n1, r_p1n0)         — negative quad part 2
    //
    // Interface faces (4 triangles, fan from c through the quad):
    //     (c, r_p0n0, r_p1n0) — between C0 and C2
    //     (c, r_p0n1, r_p1n1) — between C3 and C5
    //     (c, r_p0n0, r_p0n1) — between C6 and C8
    //     (c, r_p1n0, r_p1n1) — between C9 and C11
    //
    // Sign patterns: every child has at most 1-2 non-zero vertices
    // on the same side → 0 crossing edges.
    // -----------------------------------------------------------------
    auto children = std::array<std::array<int32_t, 4>, 12>{{
        {c, n0, r_p0n0, r_p1n0},     // C0:  neg lone
        {c, p0, p1, r_p1n0},          // C1:  pos quad 1
        {c, p0, r_p1n0, r_p0n0},      // C2:  pos quad 2
        {c, n1, r_p0n1, r_p1n1},      // C3:  neg lone
        {c, p0, p1, r_p1n1},          // C4:  pos quad 1
        {c, p0, r_p1n1, r_p0n1},      // C5:  pos quad 2
        {c, p0, r_p0n0, r_p0n1},      // C6:  pos lone
        {c, n0, n1, r_p0n1},          // C7:  neg quad 1
        {c, n0, r_p0n1, r_p0n0},      // C8:  neg quad 2
        {c, p1, r_p1n0, r_p1n1},      // C9:  pos lone
        {c, n0, n1, r_p1n1},          // C10: neg quad 1
        {c, n0, r_p1n1, r_p1n0},      // C11: neg quad 2
    }};

    orient_child_tets_like_parent(
        mesh,
        cell_id,
        std::span<std::array<int32_t, 4>>(children.data(), children.size()));
    replace_cell_with_children(mesh, cell_id, children);
    if (!certify_split_children(mesh, ls, cell_id, 12, level_set_id, tol))
        return rollback();
    return true;
}

/// Core implementation of the interface-adapted split.
template <std::floating_point T>
bool interface_split_impl(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol,
    std::unordered_map<uint64_t, int32_t>* root_cache)
{
    if (!ls.has_value())
        return false;

    const LocalMesh<T> mesh_before = mesh;
    std::unordered_map<uint64_t, int32_t> root_cache_before;
    if (root_cache)
        root_cache_before = *root_cache;
    const auto rollback = [&]() -> bool
    {
        mesh = mesh_before;
        if (root_cache)
            *root_cache = root_cache_before;
        return false;
    };

    const int n_ls = mesh.n_level_sets;

    const int c0 = mesh.cell_offsets[static_cast<std::size_t>(cell_id)];
    const int c1 = mesh.cell_offsets[static_cast<std::size_t>(cell_id + 1)];
    if (c1 - c0 != 4)
        return rollback();

    // Read tet vertex IDs and phi values
    std::array<int32_t, 4> vids;
    std::array<T, 4> phi_vals;
    for (int j = 0; j < 4; ++j)
    {
        vids[static_cast<std::size_t>(j)]
            = mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];
        phi_vals[static_cast<std::size_t>(j)]
            = mesh.vertex_phi[static_cast<std::size_t>(
                  vids[static_cast<std::size_t>(j)] * n_ls + level_set_id)];
    }

    // Pick apex: the lone-sign vertex
    const int apex_local = pick_topology1_apex<T>(
        std::span<const T>(phi_vals.data(), 4));
    if (apex_local < 0)
        return rollback();

    // Face vertex local indices (the 3 that are not apex)
    std::array<int, 3> fl;
    int fi = 0;
    for (int i = 0; i < 4; ++i)
        if (i != apex_local)
            fl[static_cast<std::size_t>(fi++)] = i;

    const int32_t v0  = vids[static_cast<std::size_t>(apex_local)];
    const int32_t fv0 = vids[static_cast<std::size_t>(fl[0])];
    const int32_t fv1 = vids[static_cast<std::size_t>(fl[1])];
    const int32_t fv2 = vids[static_cast<std::size_t>(fl[2])];

    // -----------------------------------------------------------------
    // Find 3 edge roots on crossing edges (v0-fv0), (v0-fv1), (v0-fv2)
    // -----------------------------------------------------------------
    const int32_t r0 = find_and_create_edge_root(
        mesh, v0, fv0, ls, level_set_id, tol, root_cache);
    if (r0 < 0)
        return rollback();

    const int32_t r1 = find_and_create_edge_root(
        mesh, v0, fv1, ls, level_set_id, tol, root_cache);
    if (r1 < 0)
        return rollback();

    const int32_t r2 = find_and_create_edge_root(
        mesh, v0, fv2, ls, level_set_id, tol, root_cache);
    if (r2 < 0)
        return rollback();

    // -----------------------------------------------------------------
    // Create interior vertex on interface via ray from apex through
    // centroid of the 3 edge roots
    // -----------------------------------------------------------------
    const int32_t c = create_root_centroid(
        mesh, v0, r0, r1, r2, ls, level_set_id, tol);
    if (c < 0)
        return rollback();

    // -----------------------------------------------------------------
    // Build 10 child tetrahedra
    //
    // The boundary of the original tet is decomposed as:
    //
    //   3 intersected lateral faces (each containing v0 and 2 base verts):
    //     Face (v0,fv0,fv1), roots r0 on (v0,fv0), r1 on (v0,fv1):
    //       sub-tris: (v0,r0,r1), (r0,fv0,fv1), (r0,r1,fv1)
    //
    //     Face (v0,fv1,fv2), roots r1 on (v0,fv1), r2 on (v0,fv2):
    //       sub-tris: (v0,r1,r2), (r1,fv1,fv2), (r1,r2,fv2)
    //
    //     Face (v0,fv2,fv0), roots r2 on (v0,fv2), r0 on (v0,fv0):
    //       sub-tris: (v0,r2,r0), (r2,fv2,fv0), (r2,r0,fv0)
    //
    //   1 base face: (fv0,fv1,fv2)
    //
    //   Total: 10 boundary sub-triangles
    //
    // Cone from interior root c to each sub-triangle → 10 tets:
    //
    //   Apex-side (3 tets, all contain v0 — single layer):
    //     T0 = (c, v0, r0, r1)
    //     T1 = (c, v0, r1, r2)
    //     T2 = (c, v0, r2, r0)
    //
    //   Base-side transition (6 tets, all vertices on base side or
    //   interface — clean):
    //     T3 = (c, r0, fv0, fv1)
    //     T4 = (c, r0, r1,  fv1)
    //     T5 = (c, r1, fv1, fv2)
    //     T6 = (c, r1, r2,  fv2)
    //     T7 = (c, r2, fv2, fv0)
    //     T8 = (c, r2, r0,  fv0)
    //
    //   Base (1 tet):
    //     T9 = (c, fv0, fv1, fv2)
    // -----------------------------------------------------------------
    auto children = std::array<std::array<int32_t, 4>, 10>{{
        {c, v0, r0, r1},     // T0: apex-side
        {c, v0, r1, r2},     // T1: apex-side
        {c, v0, r2, r0},     // T2: apex-side
        {c, r0, fv0, fv1},   // T3: base-side transition
        {c, r0, r1,  fv1},   // T4: base-side transition
        {c, r1, fv1, fv2},   // T5: base-side transition
        {c, r1, r2,  fv2},   // T6: base-side transition
        {c, r2, fv2, fv0},   // T7: base-side transition
        {c, r2, r0,  fv0},   // T8: base-side transition
        {c, fv0, fv1, fv2}   // T9: base
    }};

    orient_child_tets_like_parent(
        mesh,
        cell_id,
        std::span<std::array<int32_t, 4>>(children.data(), children.size()));
    replace_cell_with_children(mesh, cell_id, children);
    if (!certify_split_children(mesh, ls, cell_id, 10, level_set_id, tol))
        return rollback();
    return true;
}

// ============================================================================
// Face-enriched helpers
// ============================================================================

/// Find a root inside a triangular face on the interface using two existing
/// edge roots as the centroid target.  The source vertex must have phi != 0.
/// Implemented as a thin wrapper over create_root_centroid_n with n = 2.
///
/// @return new vertex ID, or -1 on failure
template <std::floating_point T>
int32_t find_and_create_face_root(
    LocalMesh<T>& mesh,
    int32_t source_vid,
    int32_t ra,
    int32_t rb,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol)
{
    const std::array<int32_t, 2> roots = {ra, rb};
    return create_root_centroid_n(
        mesh, source_vid,
        std::span<const int32_t>(roots.data(), 2),
        ls, level_set_id, tol);
}

/// Core implementation of the face-enriched 1-vs-3 interface split.
/// Produces 16 child tetrahedra (5 per intersected lateral face + 1 base).
///
/// Compared to interface_split_impl (10 tets), each of the 3 intersected
/// lateral faces receives one additional interface vertex in its interior
/// (f01, f12, f20).  This subdivides each face into 5 sub-triangles instead
/// of 3, providing better anchoring of the interface for non-linear level sets.
template <std::floating_point T>
bool interface_split_impl_face(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol,
    std::unordered_map<uint64_t, int32_t>* root_cache)
{
    const int n_ls = mesh.n_level_sets;

    const int c0 = mesh.cell_offsets[static_cast<std::size_t>(cell_id)];
    const int c1 = mesh.cell_offsets[static_cast<std::size_t>(cell_id + 1)];
    if (c1 - c0 != 4)
        return false;

    std::array<int32_t, 4> vids;
    std::array<T, 4> phi_vals;
    for (int j = 0; j < 4; ++j)
    {
        vids[static_cast<std::size_t>(j)]
            = mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];
        phi_vals[static_cast<std::size_t>(j)]
            = mesh.vertex_phi[static_cast<std::size_t>(
                  vids[static_cast<std::size_t>(j)] * n_ls + level_set_id)];
    }

    const int apex_local = pick_topology1_apex<T>(
        std::span<const T>(phi_vals.data(), 4));
    if (apex_local < 0)
        return false;

    std::array<int, 3> fl;
    int fi = 0;
    for (int i = 0; i < 4; ++i)
        if (i != apex_local)
            fl[static_cast<std::size_t>(fi++)] = i;

    const int32_t v0  = vids[static_cast<std::size_t>(apex_local)];
    const int32_t fv0 = vids[static_cast<std::size_t>(fl[0])];
    const int32_t fv1 = vids[static_cast<std::size_t>(fl[1])];
    const int32_t fv2 = vids[static_cast<std::size_t>(fl[2])];

    // -----------------------------------------------------------------------
    // Edge roots on the 3 crossing edges (v0-fvk)
    // -----------------------------------------------------------------------
    const int32_t r0 = find_and_create_edge_root(
        mesh, v0, fv0, ls, level_set_id, tol, root_cache);
    if (r0 < 0) return false;

    const int32_t r1 = find_and_create_edge_root(
        mesh, v0, fv1, ls, level_set_id, tol, root_cache);
    if (r1 < 0) return false;

    const int32_t r2 = find_and_create_edge_root(
        mesh, v0, fv2, ls, level_set_id, tol, root_cache);
    if (r2 < 0) return false;

    // -----------------------------------------------------------------------
    // Interior cell vertex on interface
    // -----------------------------------------------------------------------
    const int32_t c = create_root_centroid(
        mesh, v0, r0, r1, r2, ls, level_set_id, tol);
    if (c < 0) return false;

    // -----------------------------------------------------------------------
    // Face-interior interface vertices — one per intersected lateral face.
    // Use the lone-sign apex v0 as ray source: the ray from v0 through the
    // midpoint of (ra, rb) reliably brackets the interface near that midpoint.
    // -----------------------------------------------------------------------
    const int32_t f01 = find_and_create_face_root(
        mesh, v0, r0, r1, ls, level_set_id, tol);  // inside face (v0,fv0,fv1)
    if (f01 < 0) return false;

    const int32_t f12 = find_and_create_face_root(
        mesh, v0, r1, r2, ls, level_set_id, tol);  // inside face (v0,fv1,fv2)
    if (f12 < 0) return false;

    const int32_t f20 = find_and_create_face_root(
        mesh, v0, r2, r0, ls, level_set_id, tol);  // inside face (v0,fv2,fv0)
    if (f20 < 0) return false;

    // -----------------------------------------------------------------------
    // Build 16 child tetrahedra — cone from interior vertex c to
    // 5 sub-triangles per lateral face + 1 base triangle.
    //
    // Lateral face (v0, fv0, fv1): roots r0, r1; face vertex f01
    //   L0  = (c, v0, r0,  f01)    apex sub-tri 1
    //   L1  = (c, v0, f01, r1)     apex sub-tri 2
    //   L2  = (c, r0, fv0, f01)    base sub-tri 1
    //   L3  = (c, fv0, fv1, f01)   base sub-tri 2
    //   L4  = (c, fv1, r1,  f01)   base sub-tri 3
    //
    // Lateral face (v0, fv1, fv2): roots r1, r2; face vertex f12
    //   L5  = (c, v0, r1,  f12)
    //   L6  = (c, v0, f12, r2)
    //   L7  = (c, r1, fv1, f12)
    //   L8  = (c, fv1, fv2, f12)
    //   L9  = (c, fv2, r2,  f12)
    //
    // Lateral face (v0, fv2, fv0): roots r2, r0; face vertex f20
    //   L10 = (c, v0, r2,  f20)
    //   L11 = (c, v0, f20, r0)
    //   L12 = (c, r2, fv2, f20)
    //   L13 = (c, fv2, fv0, f20)
    //   L14 = (c, fv0, r0,  f20)
    //
    // Base face (fv0, fv1, fv2) — not intersected:
    //   L15 = (c, fv0, fv1, fv2)
    // -----------------------------------------------------------------------
    auto children = std::array<std::array<int32_t, 4>, 16>{{
        {c, v0,  r0,  f01},   // L0
        {c, v0,  f01, r1},    // L1
        {c, r0,  fv0, f01},   // L2
        {c, fv0, fv1, f01},   // L3
        {c, fv1, r1,  f01},   // L4
        {c, v0,  r1,  f12},   // L5
        {c, v0,  f12, r2},    // L6
        {c, r1,  fv1, f12},   // L7
        {c, fv1, fv2, f12},   // L8
        {c, fv2, r2,  f12},   // L9
        {c, v0,  r2,  f20},   // L10
        {c, v0,  f20, r0},    // L11
        {c, r2,  fv2, f20},   // L12
        {c, fv2, fv0, f20},   // L13
        {c, fv0, r0,  f20},   // L14
        {c, fv0, fv1, fv2},   // L15 base
    }};

    orient_child_tets_like_parent(
        mesh,
        cell_id,
        std::span<std::array<int32_t, 4>>(children.data(), children.size()));
    replace_cell_with_children(mesh, cell_id, children);
    return true;
}

/// Core implementation of the face-enriched 2-vs-2 interface split.
/// Produces 20 child tetrahedra (5 per intersected face × 4 faces).
///
/// Compared to interface_split_impl2 (12 tets), each face receives one
/// additional interface vertex in its interior (fpn0, fpn1, fn0p0, fn0p1).
template <std::floating_point T>
bool interface_split_impl2_face(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol,
    std::unordered_map<uint64_t, int32_t>* root_cache)
{
    const int n_ls = mesh.n_level_sets;

    const int c0 = mesh.cell_offsets[static_cast<std::size_t>(cell_id)];
    const int c1 = mesh.cell_offsets[static_cast<std::size_t>(cell_id + 1)];
    if (c1 - c0 != 4)
        return false;

    std::array<int32_t, 4> vids;
    std::array<T, 4> phi_vals;
    for (int j = 0; j < 4; ++j)
    {
        vids[static_cast<std::size_t>(j)]
            = mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];
        phi_vals[static_cast<std::size_t>(j)]
            = mesh.vertex_phi[static_cast<std::size_t>(
                  vids[static_cast<std::size_t>(j)] * n_ls + level_set_id)];
    }

    const auto pairs = pick_topology2_pairs<T>(
        std::span<const T>(phi_vals.data(), 4));
    if (pairs[0] < 0)
        return false;

    const int32_t p0 = vids[static_cast<std::size_t>(pairs[0])];
    const int32_t p1 = vids[static_cast<std::size_t>(pairs[1])];
    const int32_t n0 = vids[static_cast<std::size_t>(pairs[2])];
    const int32_t n1 = vids[static_cast<std::size_t>(pairs[3])];

    // -----------------------------------------------------------------------
    // Edge roots on all 4 crossing edges
    // -----------------------------------------------------------------------
    const int32_t r00 = find_and_create_edge_root(
        mesh, p0, n0, ls, level_set_id, tol, root_cache);
    if (r00 < 0) return false;

    const int32_t r01 = find_and_create_edge_root(
        mesh, p0, n1, ls, level_set_id, tol, root_cache);
    if (r01 < 0) return false;

    const int32_t r10 = find_and_create_edge_root(
        mesh, p1, n0, ls, level_set_id, tol, root_cache);
    if (r10 < 0) return false;

    const int32_t r11 = find_and_create_edge_root(
        mesh, p1, n1, ls, level_set_id, tol, root_cache);
    if (r11 < 0) return false;

    // -----------------------------------------------------------------------
    // Interior cell vertex: use vertex with largest |phi| as ray source
    // -----------------------------------------------------------------------
    int source_local = 0;
    T max_abs_phi = std::abs(phi_vals[0]);
    for (int i = 1; i < 4; ++i)
    {
        const T a = std::abs(phi_vals[static_cast<std::size_t>(i)]);
        if (a > max_abs_phi)
        {
            source_local = i;
            max_abs_phi = a;
        }
    }
    const int32_t source_vid = vids[static_cast<std::size_t>(source_local)];

    std::array<int32_t, 4> root_arr = {r00, r01, r10, r11};
    const int32_t c = create_root_centroid_n(
        mesh, source_vid,
        std::span<const int32_t>(root_arr.data(), 4),
        ls, level_set_id, tol);
    if (c < 0) return false;

    // -----------------------------------------------------------------------
    // Face-interior interface vertices — one per face, using the lone-sign
    // vertex of that face as ray source.
    // -----------------------------------------------------------------------
    // Face (p0, p1, n0) — lone-negative n0; roots r00, r10
    const int32_t fpn0 = find_and_create_face_root(
        mesh, n0, r00, r10, ls, level_set_id, tol);
    if (fpn0 < 0) return false;

    // Face (p0, p1, n1) — lone-negative n1; roots r01, r11
    const int32_t fpn1 = find_and_create_face_root(
        mesh, n1, r01, r11, ls, level_set_id, tol);
    if (fpn1 < 0) return false;

    // Face (p0, n0, n1) — lone-positive p0; roots r00, r01
    const int32_t fn0p0 = find_and_create_face_root(
        mesh, p0, r00, r01, ls, level_set_id, tol);
    if (fn0p0 < 0) return false;

    // Face (p1, n0, n1) — lone-positive p1; roots r10, r11
    const int32_t fn0p1 = find_and_create_face_root(
        mesh, p1, r10, r11, ls, level_set_id, tol);
    if (fn0p1 < 0) return false;

    // -----------------------------------------------------------------------
    // Build 20 child tetrahedra — cone from interior vertex c to
    // 5 sub-triangles per face × 4 faces.
    //
    // Face (p0, p1, n0) — lone=n0; face vertex fpn0; roots r00, r10
    //   F0  = (c, n0, r00,  fpn0)    lone sub-tri 1
    //   F1  = (c, n0, fpn0, r10)     lone sub-tri 2
    //   F2  = (c, p0, r00,  fpn0)    pair sub-tri 1
    //   F3  = (c, p0, p1,   fpn0)    pair sub-tri 2
    //   F4  = (c, p1, r10,  fpn0)    pair sub-tri 3
    //
    // Face (p0, p1, n1) — lone=n1; face vertex fpn1; roots r01, r11
    //   F5  = (c, n1, r01,  fpn1)
    //   F6  = (c, n1, fpn1, r11)
    //   F7  = (c, p0, r01,  fpn1)
    //   F8  = (c, p0, p1,   fpn1)
    //   F9  = (c, p1, r11,  fpn1)
    //
    // Face (p0, n0, n1) — lone=p0; face vertex fn0p0; roots r00, r01
    //   F10 = (c, p0, r00,  fn0p0)   lone sub-tri 1
    //   F11 = (c, p0, fn0p0, r01)    lone sub-tri 2
    //   F12 = (c, n0, r00,  fn0p0)   pair sub-tri 1
    //   F13 = (c, n0, n1,   fn0p0)   pair sub-tri 2
    //   F14 = (c, n1, r01,  fn0p0)   pair sub-tri 3
    //
    // Face (p1, n0, n1) — lone=p1; face vertex fn0p1; roots r10, r11
    //   F15 = (c, p1, r10,  fn0p1)   lone sub-tri 1
    //   F16 = (c, p1, fn0p1, r11)    lone sub-tri 2
    //   F17 = (c, n0, r10,  fn0p1)   pair sub-tri 1
    //   F18 = (c, n0, n1,   fn0p1)   pair sub-tri 2
    //   F19 = (c, n1, r11,  fn0p1)   pair sub-tri 3
    // -----------------------------------------------------------------------
    auto children = std::array<std::array<int32_t, 4>, 20>{{
        {c, n0, r00,   fpn0},   // F0
        {c, n0, fpn0,  r10},    // F1
        {c, p0, r00,   fpn0},   // F2
        {c, p0, p1,    fpn0},   // F3
        {c, p1, r10,   fpn0},   // F4
        {c, n1, r01,   fpn1},   // F5
        {c, n1, fpn1,  r11},    // F6
        {c, p0, r01,   fpn1},   // F7
        {c, p0, p1,    fpn1},   // F8
        {c, p1, r11,   fpn1},   // F9
        {c, p0, r00,   fn0p0},  // F10
        {c, p0, fn0p0, r01},    // F11
        {c, n0, r00,   fn0p0},  // F12
        {c, n0, n1,    fn0p0},  // F13
        {c, n1, r01,   fn0p0},  // F14
        {c, p1, r10,   fn0p1},  // F15
        {c, p1, fn0p1, r11},    // F16
        {c, n0, r10,   fn0p1},  // F17
        {c, n0, n1,    fn0p1},  // F18
        {c, n1, r11,   fn0p1},  // F19
    }};

    orient_child_tets_like_parent(
        mesh,
        cell_id,
        std::span<std::array<int32_t, 4>>(children.data(), children.size()));
    replace_cell_with_children(mesh, cell_id, children);
    return true;
}

template <std::floating_point T>
int32_t create_face_midpoint_vertex(
    LocalMesh<T>& mesh,
    int32_t v0,
    int32_t v1,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    int parent_face_id)
{
    if (!ls.has_value())
        return -1;

    const int gdim = mesh.gdim;
    const int tdim = mesh.tdim;
    std::array<T, 3> p_mid = {T(0), T(0), T(0)};
    std::array<T, 3> r_mid = {T(0), T(0), T(0)};
    for (int d = 0; d < gdim; ++d)
    {
        p_mid[static_cast<std::size_t>(d)] = T(0.5) * (
            mesh.vertex_x[static_cast<std::size_t>(v0 * gdim + d)]
          + mesh.vertex_x[static_cast<std::size_t>(v1 * gdim + d)]);
    }
    for (int d = 0; d < tdim; ++d)
    {
        r_mid[static_cast<std::size_t>(d)] = T(0.5) * (
            mesh.vertex_ref_x[static_cast<std::size_t>(v0 * tdim + d)]
          + mesh.vertex_ref_x[static_cast<std::size_t>(v1 * tdim + d)]);
    }

    const T phi_mid = ls.value(p_mid.data(), mesh.parent_cell_id);
    return create_vertex(
        mesh,
        std::span<const T>(p_mid.data(), static_cast<std::size_t>(gdim)),
        std::span<const T>(r_mid.data(), static_cast<std::size_t>(tdim)),
        parent_face_id >= 0 ? 2 : -1,
        parent_face_id,
        -1,
        level_set_id,
        phi_mid);
}

template <std::floating_point T>
int32_t create_face_root_from_vertex_to_midpoint(
    LocalMesh<T>& mesh,
    int32_t source_vid,
    int32_t midpoint_vid,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    int parent_face_id,
    T tol)
{
    if (!ls.has_value())
        return -1;

    const int gdim = mesh.gdim;
    const int tdim = mesh.tdim;
    const int n_ls = mesh.n_level_sets;

    std::array<T, 3> src_x = {T(0), T(0), T(0)};
    std::array<T, 3> src_ref = {T(0), T(0), T(0)};
    std::array<T, 3> mid_x = {T(0), T(0), T(0)};
    std::array<T, 3> mid_ref = {T(0), T(0), T(0)};
    for (int d = 0; d < gdim; ++d)
    {
        src_x[static_cast<std::size_t>(d)]
            = mesh.vertex_x[static_cast<std::size_t>(source_vid * gdim + d)];
        mid_x[static_cast<std::size_t>(d)]
            = mesh.vertex_x[static_cast<std::size_t>(midpoint_vid * gdim + d)];
    }
    for (int d = 0; d < tdim; ++d)
    {
        src_ref[static_cast<std::size_t>(d)]
            = mesh.vertex_ref_x[static_cast<std::size_t>(source_vid * tdim + d)];
        mid_ref[static_cast<std::size_t>(d)]
            = mesh.vertex_ref_x[static_cast<std::size_t>(midpoint_vid * tdim + d)];
    }

    const T phi_src = mesh.vertex_phi[static_cast<std::size_t>(
        source_vid * n_ls + level_set_id)];
    const T phi_mid = mesh.vertex_phi[static_cast<std::size_t>(
        midpoint_vid * n_ls + level_set_id)];

    T t_bracket_end = T(1);
    T phi_bracket_end = phi_mid;
    if (phi_src * phi_mid >= T(0))
    {
        bool found = false;
        for (T ext = T(2); ext <= T(8); ext *= T(2))
        {
            std::array<T, 3> x_ext = {T(0), T(0), T(0)};
            for (int d = 0; d < gdim; ++d)
            {
                x_ext[static_cast<std::size_t>(d)]
                    = src_x[static_cast<std::size_t>(d)]
                      + ext * (mid_x[static_cast<std::size_t>(d)]
                               - src_x[static_cast<std::size_t>(d)]);
            }
            const T phi_ext = ls.value(x_ext.data(), mesh.parent_cell_id);
            if (phi_src * phi_ext < T(0))
            {
                t_bracket_end = ext;
                phi_bracket_end = phi_ext;
                found = true;
                break;
            }
        }
        if (!found)
            return -1;
    }

    std::array<T, 3> x_eval = {T(0), T(0), T(0)};
    auto eval_ray = [&](T t) -> T
    {
        for (int d = 0; d < gdim; ++d)
        {
            x_eval[static_cast<std::size_t>(d)]
                = src_x[static_cast<std::size_t>(d)]
                  + t * (mid_x[static_cast<std::size_t>(d)]
                         - src_x[static_cast<std::size_t>(d)]);
        }
        return ls.value(x_eval.data(), mesh.parent_cell_id);
    };

    int iterations = 0;
    bool converged = false;
    const T t_root = cell::edge_root::brent_solve<T>(
        eval_ray, T(0), t_bracket_end, phi_src, phi_bracket_end,
        100, tol, tol, &iterations, &converged);
    const T residual = std::abs(eval_ray(t_root));
    if (!converged && residual > T(1e-8))
        return -1;
    if (!converged)
    {
        std::cerr << "Debug: interface face-root solve accepted without solver convergence"
                  << " on parent cell " << mesh.parent_cell_id
                  << ", residual=" << residual
                  << ", source_vid=" << source_vid
                  << ", midpoint_vid=" << midpoint_vid << "\n";
    }

    std::array<T, 3> p_root = {T(0), T(0), T(0)};
    std::array<T, 3> r_root = {T(0), T(0), T(0)};
    for (int d = 0; d < gdim; ++d)
    {
        p_root[static_cast<std::size_t>(d)]
            = src_x[static_cast<std::size_t>(d)]
              + t_root * (mid_x[static_cast<std::size_t>(d)]
                          - src_x[static_cast<std::size_t>(d)]);
    }
    for (int d = 0; d < tdim; ++d)
    {
        r_root[static_cast<std::size_t>(d)]
            = src_ref[static_cast<std::size_t>(d)]
              + t_root * (mid_ref[static_cast<std::size_t>(d)]
                          - src_ref[static_cast<std::size_t>(d)]);
    }

    return create_root_vertex(
        mesh,
        std::span<const T>(p_root.data(), static_cast<std::size_t>(gdim)),
        std::span<const T>(r_root.data(), static_cast<std::size_t>(tdim)),
        parent_face_id >= 0 ? 2 : -1,
        parent_face_id,
        level_set_id);
}

inline void append_face_enriched_half_face(
    std::vector<std::array<int32_t, 3>>& boundary_tris,
    int32_t lone_vertex,
    int32_t pair_vertex,
    int32_t midpoint_vertex,
    int32_t edge_root,
    int32_t face_root)
{
    boundary_tris.push_back({lone_vertex, edge_root, face_root});
    boundary_tris.push_back({edge_root, pair_vertex, midpoint_vertex});
    boundary_tris.push_back({edge_root, midpoint_vertex, face_root});
}

template <std::floating_point T>
void build_topology1_face_enriched_boundary_tris(
    std::vector<std::array<int32_t, 3>>& boundary_tris,
    int32_t v0,
    int32_t fv0,
    int32_t fv1,
    int32_t fv2,
    int32_t r0,
    int32_t r1,
    int32_t r2,
    int32_t m01,
    int32_t m12,
    int32_t m20,
    int32_t f01,
    int32_t f12,
    int32_t f20)
{
    boundary_tris.clear();
    boundary_tris.reserve(19);

    append_face_enriched_half_face(boundary_tris, v0, fv0, m01, r0, f01);
    boundary_tris.push_back({v0, f01, r1});
    boundary_tris.push_back({f01, m01, r1});
    boundary_tris.push_back({f01, fv1, r1});

    append_face_enriched_half_face(boundary_tris, v0, fv1, m12, r1, f12);
    boundary_tris.push_back({v0, f12, r2});
    boundary_tris.push_back({f12, m12, r2});
    boundary_tris.push_back({f12, fv2, r2});

    append_face_enriched_half_face(boundary_tris, v0, fv2, m20, r2, f20);
    boundary_tris.push_back({v0, f20, r0});
    boundary_tris.push_back({f20, m20, r0});
    boundary_tris.push_back({f20, fv0, r0});

    boundary_tris.push_back({fv0, fv1, fv2});
}

template <std::floating_point T>
void build_topology2_face_enriched_boundary_tris(
    std::vector<std::array<int32_t, 3>>& boundary_tris,
    int32_t p0,
    int32_t p1,
    int32_t n0,
    int32_t n1,
    int32_t r00,
    int32_t r01,
    int32_t r10,
    int32_t r11,
    int32_t m_ppn0,
    int32_t m_ppn1,
    int32_t m_pnn0,
    int32_t m_pnn1,
    int32_t f_ppn0,
    int32_t f_ppn1,
    int32_t f_pnn0,
    int32_t f_pnn1)
{
    boundary_tris.clear();
    boundary_tris.reserve(24);

    append_face_enriched_half_face(boundary_tris, n0, p0, m_ppn0, r00, f_ppn0);
    boundary_tris.push_back({n0, f_ppn0, r10});
    boundary_tris.push_back({f_ppn0, m_ppn0, r10});
    boundary_tris.push_back({f_ppn0, p1, r10});

    append_face_enriched_half_face(boundary_tris, n1, p0, m_ppn1, r01, f_ppn1);
    boundary_tris.push_back({n1, f_ppn1, r11});
    boundary_tris.push_back({f_ppn1, m_ppn1, r11});
    boundary_tris.push_back({f_ppn1, p1, r11});

    append_face_enriched_half_face(boundary_tris, p0, n0, m_pnn0, r00, f_pnn0);
    boundary_tris.push_back({p0, f_pnn0, r01});
    boundary_tris.push_back({f_pnn0, m_pnn0, r01});
    boundary_tris.push_back({f_pnn0, n1, r01});

    append_face_enriched_half_face(boundary_tris, p1, n0, m_pnn1, r10, f_pnn1);
    boundary_tris.push_back({p1, f_pnn1, r11});
    boundary_tris.push_back({f_pnn1, m_pnn1, r11});
    boundary_tris.push_back({f_pnn1, n1, r11});
}

template <std::floating_point T>
bool interface_split_impl_faces(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol,
    std::unordered_map<uint64_t, int32_t>* root_cache)
{
    if (!ls.has_value())
        return false;

    const LocalMesh<T> mesh_before = mesh;
    std::unordered_map<uint64_t, int32_t> root_cache_before;
    if (root_cache)
        root_cache_before = *root_cache;
    const auto rollback = [&]() -> bool
    {
        mesh = mesh_before;
        if (root_cache)
            *root_cache = root_cache_before;
        return false;
    };

    const int n_ls = mesh.n_level_sets;
    const int c0 = mesh.cell_offsets[static_cast<std::size_t>(cell_id)];
    const int c1 = mesh.cell_offsets[static_cast<std::size_t>(cell_id + 1)];
    if (c1 - c0 != 4)
        return rollback();

    std::array<int32_t, 4> vids;
    std::array<T, 4> phi_vals;
    for (int j = 0; j < 4; ++j)
    {
        vids[static_cast<std::size_t>(j)]
            = mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];
        phi_vals[static_cast<std::size_t>(j)]
            = mesh.vertex_phi[static_cast<std::size_t>(
                vids[static_cast<std::size_t>(j)] * n_ls + level_set_id)];
    }

    const int apex_local = pick_topology1_apex<T>(
        std::span<const T>(phi_vals.data(), 4));
    if (apex_local < 0)
        return rollback();

    std::array<int, 3> fl;
    int fi = 0;
    for (int i = 0; i < 4; ++i)
    {
        if (i != apex_local)
            fl[static_cast<std::size_t>(fi++)] = i;
    }

    const int32_t v0 = vids[static_cast<std::size_t>(apex_local)];
    const int32_t fv0 = vids[static_cast<std::size_t>(fl[0])];
    const int32_t fv1 = vids[static_cast<std::size_t>(fl[1])];
    const int32_t fv2 = vids[static_cast<std::size_t>(fl[2])];

    const int32_t r0 = find_and_create_edge_root(
        mesh, v0, fv0, ls, level_set_id, tol, root_cache);
    const int32_t r1 = find_and_create_edge_root(
        mesh, v0, fv1, ls, level_set_id, tol, root_cache);
    const int32_t r2 = find_and_create_edge_root(
        mesh, v0, fv2, ls, level_set_id, tol, root_cache);
    if (r0 < 0 || r1 < 0 || r2 < 0)
        return rollback();

    const int32_t c = create_root_centroid(
        mesh, v0, r0, r1, r2, ls, level_set_id, tol);
    if (c < 0)
        return rollback();

    const std::array<int32_t, 3> face01 = {v0, fv0, fv1};
    const std::array<int32_t, 3> face12 = {v0, fv1, fv2};
    const std::array<int32_t, 3> face20 = {v0, fv2, fv0};
    const int bg_face01 = infer_parent_tet_face_id(
        mesh, std::span<const int32_t>(face01.data(), face01.size()));
    const int bg_face12 = infer_parent_tet_face_id(
        mesh, std::span<const int32_t>(face12.data(), face12.size()));
    const int bg_face20 = infer_parent_tet_face_id(
        mesh, std::span<const int32_t>(face20.data(), face20.size()));

    const int32_t m01 = create_face_midpoint_vertex(
        mesh, fv0, fv1, ls, level_set_id, bg_face01);
    const int32_t m12 = create_face_midpoint_vertex(
        mesh, fv1, fv2, ls, level_set_id, bg_face12);
    const int32_t m20 = create_face_midpoint_vertex(
        mesh, fv2, fv0, ls, level_set_id, bg_face20);
    if (m01 < 0 || m12 < 0 || m20 < 0)
        return rollback();

    const int32_t f01 = create_face_root_from_vertex_to_midpoint(
        mesh, v0, m01, ls, level_set_id, bg_face01, tol);
    const int32_t f12 = create_face_root_from_vertex_to_midpoint(
        mesh, v0, m12, ls, level_set_id, bg_face12, tol);
    const int32_t f20 = create_face_root_from_vertex_to_midpoint(
        mesh, v0, m20, ls, level_set_id, bg_face20, tol);
    if (f01 < 0 || f12 < 0 || f20 < 0)
        return rollback();

    std::vector<std::array<int32_t, 3>> boundary_tris;
    build_topology1_face_enriched_boundary_tris<T>(
        boundary_tris,
        v0, fv0, fv1, fv2,
        r0, r1, r2,
        m01, m12, m20,
        f01, f12, f20);

    std::vector<std::array<int32_t, 4>> children;
    children.reserve(boundary_tris.size());
    for (const auto& tri : boundary_tris)
        children.push_back({c, tri[0], tri[1], tri[2]});

    orient_child_tets_like_parent(
        mesh,
        cell_id,
        std::span<std::array<int32_t, 4>>(children.data(), children.size()));
    replace_cell_with_children(
        mesh,
        cell_id,
        std::span<const std::array<int32_t, 4>>(children.data(), children.size()));
    if (!certify_split_children(mesh, ls, cell_id, 19, level_set_id, tol))
        return rollback();
    return true;
}

template <std::floating_point T>
bool interface_split_impl2_faces(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol,
    std::unordered_map<uint64_t, int32_t>* root_cache)
{
    if (!ls.has_value())
        return false;

    const LocalMesh<T> mesh_before = mesh;
    std::unordered_map<uint64_t, int32_t> root_cache_before;
    if (root_cache)
        root_cache_before = *root_cache;
    const auto rollback = [&]() -> bool
    {
        mesh = mesh_before;
        if (root_cache)
            *root_cache = root_cache_before;
        return false;
    };

    const int n_ls = mesh.n_level_sets;
    const int c0 = mesh.cell_offsets[static_cast<std::size_t>(cell_id)];
    const int c1 = mesh.cell_offsets[static_cast<std::size_t>(cell_id + 1)];
    if (c1 - c0 != 4)
        return rollback();

    std::array<int32_t, 4> vids;
    std::array<T, 4> phi_vals;
    for (int j = 0; j < 4; ++j)
    {
        vids[static_cast<std::size_t>(j)]
            = mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];
        phi_vals[static_cast<std::size_t>(j)]
            = mesh.vertex_phi[static_cast<std::size_t>(
                vids[static_cast<std::size_t>(j)] * n_ls + level_set_id)];
    }

    const auto pairs = pick_topology2_pairs<T>(
        std::span<const T>(phi_vals.data(), 4));
    if (pairs[0] < 0)
        return rollback();

    const int32_t p0 = vids[static_cast<std::size_t>(pairs[0])];
    const int32_t p1 = vids[static_cast<std::size_t>(pairs[1])];
    const int32_t n0 = vids[static_cast<std::size_t>(pairs[2])];
    const int32_t n1 = vids[static_cast<std::size_t>(pairs[3])];

    const int32_t r00 = find_and_create_edge_root(
        mesh, p0, n0, ls, level_set_id, tol, root_cache);
    const int32_t r01 = find_and_create_edge_root(
        mesh, p0, n1, ls, level_set_id, tol, root_cache);
    const int32_t r10 = find_and_create_edge_root(
        mesh, p1, n0, ls, level_set_id, tol, root_cache);
    const int32_t r11 = find_and_create_edge_root(
        mesh, p1, n1, ls, level_set_id, tol, root_cache);
    if (r00 < 0 || r01 < 0 || r10 < 0 || r11 < 0)
        return rollback();

    int source_local = 0;
    T max_abs_phi = std::abs(phi_vals[0]);
    for (int i = 1; i < 4; ++i)
    {
        const T a = std::abs(phi_vals[static_cast<std::size_t>(i)]);
        if (a > max_abs_phi)
        {
            source_local = i;
            max_abs_phi = a;
        }
    }
    const int32_t source_vid = vids[static_cast<std::size_t>(source_local)];
    const std::array<int32_t, 4> roots = {r00, r01, r10, r11};
    const int32_t c = create_root_centroid_n(
        mesh,
        source_vid,
        std::span<const int32_t>(roots.data(), roots.size()),
        ls,
        level_set_id,
        tol);
    if (c < 0)
        return rollback();

    const std::array<int32_t, 3> face_ppn0 = {p0, p1, n0};
    const std::array<int32_t, 3> face_ppn1 = {p0, p1, n1};
    const std::array<int32_t, 3> face_pnn0 = {p0, n0, n1};
    const std::array<int32_t, 3> face_pnn1 = {p1, n0, n1};
    const int bg_face_ppn0 = infer_parent_tet_face_id(
        mesh, std::span<const int32_t>(face_ppn0.data(), face_ppn0.size()));
    const int bg_face_ppn1 = infer_parent_tet_face_id(
        mesh, std::span<const int32_t>(face_ppn1.data(), face_ppn1.size()));
    const int bg_face_pnn0 = infer_parent_tet_face_id(
        mesh, std::span<const int32_t>(face_pnn0.data(), face_pnn0.size()));
    const int bg_face_pnn1 = infer_parent_tet_face_id(
        mesh, std::span<const int32_t>(face_pnn1.data(), face_pnn1.size()));

    const int32_t m_ppn0 = create_face_midpoint_vertex(
        mesh, p0, p1, ls, level_set_id, bg_face_ppn0);
    const int32_t m_ppn1 = create_face_midpoint_vertex(
        mesh, p0, p1, ls, level_set_id, bg_face_ppn1);
    const int32_t m_pnn0 = create_face_midpoint_vertex(
        mesh, n0, n1, ls, level_set_id, bg_face_pnn0);
    const int32_t m_pnn1 = create_face_midpoint_vertex(
        mesh, n0, n1, ls, level_set_id, bg_face_pnn1);
    if (m_ppn0 < 0 || m_ppn1 < 0 || m_pnn0 < 0 || m_pnn1 < 0)
        return rollback();

    const int32_t f_ppn0 = create_face_root_from_vertex_to_midpoint(
        mesh, n0, m_ppn0, ls, level_set_id, bg_face_ppn0, tol);
    const int32_t f_ppn1 = create_face_root_from_vertex_to_midpoint(
        mesh, n1, m_ppn1, ls, level_set_id, bg_face_ppn1, tol);
    const int32_t f_pnn0 = create_face_root_from_vertex_to_midpoint(
        mesh, p0, m_pnn0, ls, level_set_id, bg_face_pnn0, tol);
    const int32_t f_pnn1 = create_face_root_from_vertex_to_midpoint(
        mesh, p1, m_pnn1, ls, level_set_id, bg_face_pnn1, tol);
    if (f_ppn0 < 0 || f_ppn1 < 0 || f_pnn0 < 0 || f_pnn1 < 0)
        return rollback();

    std::vector<std::array<int32_t, 3>> boundary_tris;
    build_topology2_face_enriched_boundary_tris<T>(
        boundary_tris,
        p0, p1, n0, n1,
        r00, r01, r10, r11,
        m_ppn0, m_ppn1, m_pnn0, m_pnn1,
        f_ppn0, f_ppn1, f_pnn0, f_pnn1);

    std::vector<std::array<int32_t, 4>> children;
    children.reserve(boundary_tris.size());
    for (const auto& tri : boundary_tris)
        children.push_back({c, tri[0], tri[1], tri[2]});

    orient_child_tets_like_parent(
        mesh,
        cell_id,
        std::span<std::array<int32_t, 4>>(children.data(), children.size()));
    replace_cell_with_children(
        mesh,
        cell_id,
        std::span<const std::array<int32_t, 4>>(children.data(), children.size()));
    if (!certify_split_children(mesh, ls, cell_id, 24, level_set_id, tol))
        return rollback();
    return true;
}

} // anonymous namespace

// ============================================================================
// interface_split_topology1_tet (public)
// ============================================================================

template <std::floating_point T>
bool interface_split_topology1_tet(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol)
{
    return interface_split_impl(
        mesh, cell_id, ls, level_set_id, tol, nullptr);
}

// ============================================================================
// interface_refine_topology1 (public)
// ============================================================================

template <std::floating_point T>
int interface_refine_topology1(
    LocalMesh<T>& mesh,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    int max_depth,
    T tol)
{
    const int n_ls = mesh.n_level_sets;
    int total_splits = 0;

    // Root cache persists across depth levels so edge-root vertices
    // shared between sibling cells are reused.
    std::unordered_map<uint64_t, int32_t> root_cache;

    for (int depth = 0; depth < max_depth; ++depth)
    {
        // Collect all 1-vs-3 intersected tets at this level
        std::vector<int> candidates;
        const int nc = mesh.n_cells();
        for (int ci = 0; ci < nc; ++ci)
        {
            if (mesh.cell_types[static_cast<std::size_t>(ci)]
                != cell::type::tetrahedron)
                continue;

            const int s0 = mesh.cell_offsets[static_cast<std::size_t>(ci)];
            const int s1 = mesh.cell_offsets[static_cast<std::size_t>(ci + 1)];
            if (s1 - s0 != 4)
                continue;

            // Count sign pattern (ignoring vertices with phi = 0)
            int n_pos = 0;
            int n_neg = 0;
            for (int j = s0; j < s1; ++j)
            {
                const int32_t vid
                    = mesh.cell_vertices[static_cast<std::size_t>(j)];
                const T phi
                    = mesh.vertex_phi[static_cast<std::size_t>(
                          vid * n_ls + level_set_id)];
                if (phi > tol)
                    ++n_pos;
                else if (phi < -tol)
                    ++n_neg;
            }

            // Only split cells with exactly 1-vs-3 sign pattern
            // among non-zero vertices.  After an interface split,
            // the apex-side children have 1 non-zero vertex (v0)
            // and 3 zero vertices (c, ri, rj), which gives n_pos=1
            // n_neg=0 or n_pos=0 n_neg=1.  These are NOT candidates
            // since they don't have the opposite-sign count of 3.
            if ((n_pos == 1 && n_neg == 3) || (n_pos == 3 && n_neg == 1))
                candidates.push_back(ci);
        }

        if (candidates.empty())
            break;

        // Process in reverse order so that cell replacement doesn't
        // invalidate earlier indices.
        std::sort(candidates.rbegin(), candidates.rend());

        for (int ci : candidates)
        {
            if (interface_split_impl(
                    mesh, ci, ls, level_set_id, tol, &root_cache))
            {
                ++total_splits;
            }
        }
    }

    return total_splits;
}

// ============================================================================
// Explicit template instantiations
// ============================================================================

template bool interface_split_topology1_tet<float>(
    LocalMesh<float>&, int, const LevelSetFunction<float>&, int, float);
template bool interface_split_topology1_tet<double>(
    LocalMesh<double>&, int, const LevelSetFunction<double>&, int, double);

template int interface_refine_topology1<float>(
    LocalMesh<float>&, const LevelSetFunction<float>&, int, int, float);
template int interface_refine_topology1<double>(
    LocalMesh<double>&, const LevelSetFunction<double>&, int, int, double);

// ============================================================================
// interface_split_topology2_tet (public)
// ============================================================================

template <std::floating_point T>
bool interface_split_topology2_tet(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol)
{
    return interface_split_impl2(
        mesh, cell_id, ls, level_set_id, tol, nullptr);
}

// ============================================================================
// interface_refine (public) — combined 1-vs-3 + 2-vs-2
// ============================================================================

template <std::floating_point T>
int interface_refine(
    LocalMesh<T>& mesh,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    int max_depth,
    T tol)
{
    const int n_ls = mesh.n_level_sets;
    int total_splits = 0;

    std::unordered_map<uint64_t, int32_t> root_cache;

    for (int depth = 0; depth < max_depth; ++depth)
    {
        // Collect candidates: 1-vs-3 and 2-vs-2 intersected tets
        // topo_type: 1 = 1-vs-3, 2 = 2-vs-2
        std::vector<std::pair<int, int>> candidates; // (cell_id, topo_type)
        const int nc = mesh.n_cells();
        for (int ci = 0; ci < nc; ++ci)
        {
            if (mesh.cell_types[static_cast<std::size_t>(ci)]
                != cell::type::tetrahedron)
                continue;

            const int s0 = mesh.cell_offsets[static_cast<std::size_t>(ci)];
            const int s1 = mesh.cell_offsets[static_cast<std::size_t>(ci + 1)];
            if (s1 - s0 != 4)
                continue;

            int n_pos = 0;
            int n_neg = 0;
            for (int j = s0; j < s1; ++j)
            {
                const int32_t vid
                    = mesh.cell_vertices[static_cast<std::size_t>(j)];
                const T phi
                    = mesh.vertex_phi[static_cast<std::size_t>(
                          vid * n_ls + level_set_id)];
                if (phi > tol)
                    ++n_pos;
                else if (phi < -tol)
                    ++n_neg;
            }

            if ((n_pos == 1 && n_neg == 3) || (n_pos == 3 && n_neg == 1))
                candidates.emplace_back(ci, 1);
            else if (n_pos == 2 && n_neg == 2)
                candidates.emplace_back(ci, 2);
        }

        if (candidates.empty())
            break;

        // Process in reverse order (by cell_id) so that cell replacement
        // doesn't invalidate earlier indices.
        std::sort(candidates.begin(), candidates.end(),
                  [](const auto& a, const auto& b) {
                      return a.first > b.first;
                  });

        for (const auto& [ci, topo] : candidates)
        {
            bool ok = false;
            if (topo == 1)
                ok = interface_split_impl(
                    mesh, ci, ls, level_set_id, tol, &root_cache);
            else
                ok = interface_split_impl2(
                    mesh, ci, ls, level_set_id, tol, &root_cache);

            if (ok)
                ++total_splits;
        }
    }

    return total_splits;
}

// ============================================================================
// Explicit template instantiations
// ============================================================================

template bool interface_split_topology2_tet<float>(
    LocalMesh<float>&, int, const LevelSetFunction<float>&, int, float);
template bool interface_split_topology2_tet<double>(
    LocalMesh<double>&, int, const LevelSetFunction<double>&, int, double);

template int interface_refine<float>(
    LocalMesh<float>&, const LevelSetFunction<float>&, int, int, float);
template int interface_refine<double>(
    LocalMesh<double>&, const LevelSetFunction<double>&, int, int, double);

// ============================================================================
// interface_split_topology1_tet_face (public) — 16-tet with face vertices
// ============================================================================

template <std::floating_point T>
bool interface_split_topology1_tet_face(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol)
{
    return interface_split_impl_face(
        mesh, cell_id, ls, level_set_id, tol, nullptr);
}

// ============================================================================
// interface_split_topology2_tet_face (public) — 20-tet with face vertices
// ============================================================================

template <std::floating_point T>
bool interface_split_topology2_tet_face(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol)
{
    return interface_split_impl2_face(
        mesh, cell_id, ls, level_set_id, tol, nullptr);
}

// ============================================================================
// interface_refine_face (public) — face-enriched 1-vs-3 + 2-vs-2
// ============================================================================

template <std::floating_point T>
int interface_refine_face(
    LocalMesh<T>& mesh,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    int max_depth,
    T tol)
{
    const int n_ls = mesh.n_level_sets;
    int total_splits = 0;

    std::unordered_map<uint64_t, int32_t> root_cache;

    for (int depth = 0; depth < max_depth; ++depth)
    {
        std::vector<std::pair<int, int>> candidates;
        const int nc = mesh.n_cells();
        for (int ci = 0; ci < nc; ++ci)
        {
            if (mesh.cell_types[static_cast<std::size_t>(ci)]
                != cell::type::tetrahedron)
                continue;

            const int s0 = mesh.cell_offsets[static_cast<std::size_t>(ci)];
            const int s1 = mesh.cell_offsets[static_cast<std::size_t>(ci + 1)];
            if (s1 - s0 != 4)
                continue;

            int n_pos = 0;
            int n_neg = 0;
            for (int j = s0; j < s1; ++j)
            {
                const int32_t vid
                    = mesh.cell_vertices[static_cast<std::size_t>(j)];
                const T phi
                    = mesh.vertex_phi[static_cast<std::size_t>(
                          vid * n_ls + level_set_id)];
                if (phi > tol)
                    ++n_pos;
                else if (phi < -tol)
                    ++n_neg;
            }

            if ((n_pos == 1 && n_neg == 3) || (n_pos == 3 && n_neg == 1))
                candidates.emplace_back(ci, 1);
            else if (n_pos == 2 && n_neg == 2)
                candidates.emplace_back(ci, 2);
        }

        if (candidates.empty())
            break;

        std::sort(candidates.begin(), candidates.end(),
                  [](const auto& a, const auto& b) {
                      return a.first > b.first;
                  });

        for (const auto& [ci, topo] : candidates)
        {
            bool ok = false;
            if (topo == 1)
                ok = interface_split_impl_face(
                    mesh, ci, ls, level_set_id, tol, &root_cache);
            else
                ok = interface_split_impl2_face(
                    mesh, ci, ls, level_set_id, tol, &root_cache);

            if (ok)
                ++total_splits;
        }
    }

    return total_splits;
}

// ============================================================================
// interface_split_topology1_tet_faces (public) — 19-tet median face split
// ============================================================================

template <std::floating_point T>
bool interface_split_topology1_tet_faces(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol)
{
    return interface_split_impl_faces(
        mesh, cell_id, ls, level_set_id, tol, nullptr);
}

// ============================================================================
// interface_split_topology2_tet_faces (public) — 24-tet median face split
// ============================================================================

template <std::floating_point T>
bool interface_split_topology2_tet_faces(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol)
{
    return interface_split_impl2_faces(
        mesh, cell_id, ls, level_set_id, tol, nullptr);
}

// ============================================================================
// interface_refine_faces (public) — median face-enriched 1-vs-3 + 2-vs-2
// ============================================================================

template <std::floating_point T>
int interface_refine_faces(
    LocalMesh<T>& mesh,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    int max_depth,
    T tol)
{
    const int n_ls = mesh.n_level_sets;
    int total_splits = 0;

    std::unordered_map<uint64_t, int32_t> root_cache;

    for (int depth = 0; depth < max_depth; ++depth)
    {
        std::vector<std::pair<int, int>> candidates;
        const int nc = mesh.n_cells();
        for (int ci = 0; ci < nc; ++ci)
        {
            if (mesh.cell_types[static_cast<std::size_t>(ci)]
                != cell::type::tetrahedron)
                continue;

            const int s0 = mesh.cell_offsets[static_cast<std::size_t>(ci)];
            const int s1 = mesh.cell_offsets[static_cast<std::size_t>(ci + 1)];
            if (s1 - s0 != 4)
                continue;

            int n_pos = 0;
            int n_neg = 0;
            for (int j = s0; j < s1; ++j)
            {
                const int32_t vid
                    = mesh.cell_vertices[static_cast<std::size_t>(j)];
                const T phi
                    = mesh.vertex_phi[static_cast<std::size_t>(
                          vid * n_ls + level_set_id)];
                if (phi > tol)
                    ++n_pos;
                else if (phi < -tol)
                    ++n_neg;
            }

            if ((n_pos == 1 && n_neg == 3) || (n_pos == 3 && n_neg == 1))
                candidates.emplace_back(ci, 1);
            else if (n_pos == 2 && n_neg == 2)
                candidates.emplace_back(ci, 2);
        }

        if (candidates.empty())
            break;

        std::sort(candidates.begin(), candidates.end(),
                  [](const auto& a, const auto& b) {
                      return a.first > b.first;
                  });

        for (const auto& [ci, topo] : candidates)
        {
            bool ok = false;
            if (topo == 1)
                ok = interface_split_impl_faces(
                    mesh, ci, ls, level_set_id, tol, &root_cache);
            else
                ok = interface_split_impl2_faces(
                    mesh, ci, ls, level_set_id, tol, &root_cache);

            if (ok)
                ++total_splits;
        }
    }

    return total_splits;
}

// ============================================================================
// Explicit template instantiations — face-enriched variants
// ============================================================================

template bool interface_split_topology1_tet_face<float>(
    LocalMesh<float>&, int, const LevelSetFunction<float>&, int, float);
template bool interface_split_topology1_tet_face<double>(
    LocalMesh<double>&, int, const LevelSetFunction<double>&, int, double);

template bool interface_split_topology2_tet_face<float>(
    LocalMesh<float>&, int, const LevelSetFunction<float>&, int, float);
template bool interface_split_topology2_tet_face<double>(
    LocalMesh<double>&, int, const LevelSetFunction<double>&, int, double);

template int interface_refine_face<float>(
    LocalMesh<float>&, const LevelSetFunction<float>&, int, int, float);
template int interface_refine_face<double>(
    LocalMesh<double>&, const LevelSetFunction<double>&, int, int, double);

template bool interface_split_topology1_tet_faces<float>(
    LocalMesh<float>&, int, const LevelSetFunction<float>&, int, float);
template bool interface_split_topology1_tet_faces<double>(
    LocalMesh<double>&, int, const LevelSetFunction<double>&, int, double);

template bool interface_split_topology2_tet_faces<float>(
    LocalMesh<float>&, int, const LevelSetFunction<float>&, int, float);
template bool interface_split_topology2_tet_faces<double>(
    LocalMesh<double>&, int, const LevelSetFunction<double>&, int, double);

template int interface_refine_faces<float>(
    LocalMesh<float>&, const LevelSetFunction<float>&, int, int, float);
template int interface_refine_faces<double>(
    LocalMesh<double>&, const LevelSetFunction<double>&, int, int, double);

} // namespace cutcells
