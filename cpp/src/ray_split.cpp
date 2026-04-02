// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "ray_split.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace cutcells
{

// ============================================================================
// pick_ray_apex
// ============================================================================

template <std::floating_point T>
int pick_ray_apex(std::span<const T> phi)
{
    assert(phi.size() == 4);

    // Count positive and negative vertices
    int n_pos = 0;
    int n_neg = 0;
    for (int i = 0; i < 4; ++i)
    {
        if (phi[static_cast<std::size_t>(i)] > T(0))
            ++n_pos;
        else if (phi[static_cast<std::size_t>(i)] < T(0))
            ++n_neg;
    }

    // 1-vs-3 case: pick the lone vertex
    if (n_pos == 1)
    {
        for (int i = 0; i < 4; ++i)
            if (phi[static_cast<std::size_t>(i)] > T(0))
                return i;
    }
    if (n_neg == 1)
    {
        for (int i = 0; i < 4; ++i)
            if (phi[static_cast<std::size_t>(i)] < T(0))
                return i;
    }

    // 2-vs-2 case: pick the vertex with largest |phi|.
    // This gives the longest ray, which tends to produce a root
    // closer to the face centroid (better-shaped children).
    int best = 0;
    T best_abs = std::abs(phi[0]);
    for (int i = 1; i < 4; ++i)
    {
        const T a = std::abs(phi[static_cast<std::size_t>(i)]);
        if (a > best_abs)
        {
            best = i;
            best_abs = a;
        }
    }
    return best;
}

// ============================================================================
// stellar_split_tetrahedron
// ============================================================================

int stellar_split_tetrahedron(
    std::span<const int32_t> cv,
    int32_t p,
    std::vector<int32_t>& out_cv,
    std::vector<int32_t>& out_off,
    std::vector<cell::type>& out_ct)
{
    assert(cv.size() == 4);

    // 4 children: replace each vertex with p in turn
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            out_cv.push_back(j == i ? p : cv[static_cast<std::size_t>(j)]);
        }
        out_off.push_back(static_cast<int32_t>(out_cv.size()));
        out_ct.push_back(cell::type::tetrahedron);
    }
    return 4;
}

// ============================================================================
// ray_split_one_tet
// ============================================================================

template <std::floating_point T>
bool ray_split_one_tet(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol)
{
    const int n_ls = mesh.n_level_sets;
    const int c0 = mesh.cell_offsets[static_cast<std::size_t>(cell_id)];
    const int c1 = mesh.cell_offsets[static_cast<std::size_t>(cell_id + 1)];
    if (c1 - c0 != 4)
        return false; // not a tet

    // Read tet vertex indices
    std::array<int32_t, 4> vids;
    for (int j = 0; j < 4; ++j)
        vids[static_cast<std::size_t>(j)]
            = mesh.cell_vertices[static_cast<std::size_t>(c0 + j)];

    // Read phi values at vertices
    std::array<T, 4> phi_vals;
    for (int j = 0; j < 4; ++j)
        phi_vals[static_cast<std::size_t>(j)]
            = mesh.vertex_phi[static_cast<std::size_t>(
                  vids[static_cast<std::size_t>(j)] * n_ls + level_set_id)];

    // Pick the apex vertex (local index 0..3)
    // First try the preferred apex from pick_ray_apex, then fall back
    // to others if the ray doesn't give a sign change bracket.
    const int preferred_apex = pick_ray_apex<T>(
        std::span<const T>(phi_vals.data(), 4));

    int apex_local = -1;
    std::array<T, 3> face_x = {T(0), T(0), T(0)};
    std::array<T, 3> face_ref = {T(0), T(0), T(0)};
    T phi_apex = T(0);
    T phi_face = T(0);

    // Order of vertices to try: preferred first, then the rest
    std::array<int, 4> try_order;
    try_order[0] = preferred_apex;
    int idx = 1;
    for (int i = 0; i < 4; ++i)
    {
        if (i != preferred_apex)
            try_order[static_cast<std::size_t>(idx++)] = i;
    }

    for (int attempt = 0; attempt < 4; ++attempt)
    {
        const int candidate = try_order[static_cast<std::size_t>(attempt)];
        const T phi_cand = phi_vals[static_cast<std::size_t>(candidate)];

        // Compute centroid of the opposite face (the 3 other vertices)
        std::array<T, 3> fx = {T(0), T(0), T(0)};
        std::array<T, 3> fr = {T(0), T(0), T(0)};
        for (int j = 0; j < 4; ++j)
        {
            if (j == candidate)
                continue;
            const int32_t vid = vids[static_cast<std::size_t>(j)];
            for (int d = 0; d < mesh.gdim; ++d)
                fx[static_cast<std::size_t>(d)]
                    += mesh.vertex_x[static_cast<std::size_t>(vid * mesh.gdim + d)];
            for (int d = 0; d < mesh.tdim; ++d)
                fr[static_cast<std::size_t>(d)]
                    += mesh.vertex_ref_x[static_cast<std::size_t>(vid * mesh.tdim + d)];
        }
        for (int d = 0; d < mesh.gdim; ++d)
            fx[static_cast<std::size_t>(d)] /= T(3);
        for (int d = 0; d < mesh.tdim; ++d)
            fr[static_cast<std::size_t>(d)] /= T(3);

        // Evaluate phi at face centroid
        std::array<T, 3> x_tmp = fx;
        const T pf = ls.value(x_tmp.data(), mesh.parent_cell_id);

        if (phi_cand * pf < T(0))
        {
            apex_local = candidate;
            face_x = fx;
            face_ref = fr;
            phi_apex = phi_cand;
            phi_face = pf;
            break;
        }
    }

    if (apex_local < 0)
        return false; // no vertex gives a sign-change bracket

    // Apex vertex coordinates (physical and reference)
    const int32_t apex_vid = vids[static_cast<std::size_t>(apex_local)];
    std::array<T, 3> apex_x = {T(0), T(0), T(0)};
    std::array<T, 3> apex_ref = {T(0), T(0), T(0)};
    for (int d = 0; d < mesh.gdim; ++d)
        apex_x[static_cast<std::size_t>(d)]
            = mesh.vertex_x[static_cast<std::size_t>(apex_vid * mesh.gdim + d)];
    for (int d = 0; d < mesh.tdim; ++d)
        apex_ref[static_cast<std::size_t>(d)]
            = mesh.vertex_ref_x[static_cast<std::size_t>(apex_vid * mesh.tdim + d)];

    // Root finding using ITP on [0, 1]
    std::array<T, 3> x_eval = {T(0), T(0), T(0)};
    auto eval_ray = [&](T t) -> T
    {
        for (int d = 0; d < mesh.gdim; ++d)
            x_eval[static_cast<std::size_t>(d)]
                = apex_x[static_cast<std::size_t>(d)]
                  + t * (face_x[static_cast<std::size_t>(d)]
                         - apex_x[static_cast<std::size_t>(d)]);
        return ls.value(x_eval.data(), mesh.parent_cell_id);
    };

    const T t_linear = cell::edge_root::linear_root_parameter(phi_apex, phi_face);
    int iterations = 0;
    bool converged = false;
    const T t_root = cell::edge_root::itp_parameter<T>(
        eval_ray, T(0), T(1), phi_apex, phi_face,
        t_linear, /*max_iter=*/100, /*xtol=*/tol, /*ftol=*/tol,
        &iterations, &converged);

    if (!converged && std::abs(eval_ray(t_root)) > T(1e-8))
        return false; // root finding failed

    // Compute coordinates of the new interior vertex
    std::array<T, 3> p_x = {T(0), T(0), T(0)};
    std::array<T, 3> p_ref = {T(0), T(0), T(0)};
    for (int d = 0; d < mesh.gdim; ++d)
        p_x[static_cast<std::size_t>(d)]
            = apex_x[static_cast<std::size_t>(d)]
              + t_root * (face_x[static_cast<std::size_t>(d)]
                          - apex_x[static_cast<std::size_t>(d)]);
    for (int d = 0; d < mesh.tdim; ++d)
        p_ref[static_cast<std::size_t>(d)]
            = apex_ref[static_cast<std::size_t>(d)]
              + t_root * (face_ref[static_cast<std::size_t>(d)]
                          - apex_ref[static_cast<std::size_t>(d)]);

    // Append the new interior vertex
    const int32_t p_vid = static_cast<int32_t>(mesh.vertex_parent_dim.size());
    for (int d = 0; d < mesh.gdim; ++d)
        mesh.vertex_x.push_back(p_x[static_cast<std::size_t>(d)]);
    for (int d = 0; d < mesh.tdim; ++d)
        mesh.vertex_ref_x.push_back(p_ref[static_cast<std::size_t>(d)]);

    // Parent tracking: interior point → parent_dim = 3 (cell interior)
    mesh.vertex_parent_dim.push_back(3);
    mesh.vertex_parent_id.push_back(mesh.parent_cell_id);
    mesh.vertex_root_edge_id.push_back(-1);

    const int nv = mesh.n_vertices();
    mesh.vertex_zero_mask.resize(static_cast<std::size_t>(nv), 0);
    mesh.vertex_inside_mask.resize(static_cast<std::size_t>(nv), 0);
    mesh.vertex_phi.resize(static_cast<std::size_t>(nv * n_ls), T(0));

    // Set phi = 0 at the root vertex (it is on the interface)
    mesh.vertex_phi[static_cast<std::size_t>(p_vid * n_ls + level_set_id)] = T(0);
    // Set zero_mask
    mesh.vertex_zero_mask[static_cast<std::size_t>(p_vid)]
        |= (uint64_t(1) << static_cast<unsigned>(level_set_id));
    // Set inside_mask based on the apex sign (same side as apex is NOT inside
    // unless apex is inside)
    // The root is on the interface, so inside_mask for this level set should
    // be 0 (not strictly inside).

    // Rebuild cell arrays: replace cell_id with 4 children, keep others
    const int old_nc = mesh.n_cells();

    std::vector<int32_t>    new_cell_vertices;
    std::vector<int32_t>    new_cell_offsets = {0};
    std::vector<cell::type> new_cell_types;

    new_cell_vertices.reserve(mesh.cell_vertices.size() + 12); // 4 new - 1 old = +12
    new_cell_offsets.reserve(static_cast<std::size_t>(old_nc + 3 + 1));
    new_cell_types.reserve(static_cast<std::size_t>(old_nc + 3));

    for (int c = 0; c < old_nc; ++c)
    {
        const int s0 = mesh.cell_offsets[static_cast<std::size_t>(c)];
        const int s1 = mesh.cell_offsets[static_cast<std::size_t>(c + 1)];

        if (c != cell_id)
        {
            for (int j = s0; j < s1; ++j)
                new_cell_vertices.push_back(
                    mesh.cell_vertices[static_cast<std::size_t>(j)]);
            new_cell_offsets.push_back(
                static_cast<int32_t>(new_cell_vertices.size()));
            new_cell_types.push_back(mesh.cell_types[static_cast<std::size_t>(c)]);
        }
        else
        {
            // Replace with 4 stellar children
            stellar_split_tetrahedron(
                std::span<const int32_t>(
                    mesh.cell_vertices.data() + static_cast<std::size_t>(s0),
                    static_cast<std::size_t>(s1 - s0)),
                p_vid,
                new_cell_vertices, new_cell_offsets, new_cell_types);
        }
    }

    // Commit
    mesh.cell_vertices = std::move(new_cell_vertices);
    mesh.cell_offsets = std::move(new_cell_offsets);
    mesh.cell_types = std::move(new_cell_types);
    mesh.cell_domain.assign(
        mesh.cell_types.size(), static_cast<uint8_t>(cell::domain::unset));

    // Rebuild connectivity
    build_local_edges(mesh);
    build_local_faces(mesh);
    rebuild_parent_entity_maps(mesh);

    return true;
}

// ============================================================================
// ray_refine_local_mesh
// ============================================================================

template <std::floating_point T>
int ray_refine_local_mesh(
    LocalMesh<T>& mesh,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    int max_depth,
    T tol)
{
    const int n_ls = mesh.n_level_sets;
    int total_splits = 0;

    for (int depth = 0; depth < max_depth; ++depth)
    {
        // Find all intersected tets at this level
        std::vector<int> intersected;
        const int nc = mesh.n_cells();
        for (int c = 0; c < nc; ++c)
        {
            if (mesh.cell_types[static_cast<std::size_t>(c)]
                != cell::type::tetrahedron)
                continue;

            const int s0 = mesh.cell_offsets[static_cast<std::size_t>(c)];
            const int s1 = mesh.cell_offsets[static_cast<std::size_t>(c + 1)];
            if (s1 - s0 != 4)
                continue;

            bool has_pos = false;
            bool has_neg = false;
            for (int j = s0; j < s1; ++j)
            {
                const int32_t vid
                    = mesh.cell_vertices[static_cast<std::size_t>(j)];
                const T phi
                    = mesh.vertex_phi[static_cast<std::size_t>(
                          vid * n_ls + level_set_id)];
                if (phi > T(0))
                    has_pos = true;
                else if (phi < T(0))
                    has_neg = true;
            }

            // A cell is intersected if it has both positive and negative vertices
            // (vertices with phi=0 alone don't count as intersected for
            //  ray refinement since we need a sign change on the ray)
            if (has_pos && has_neg)
                intersected.push_back(c);
        }

        if (intersected.empty())
            break;

        // Process in reverse order so indices remain valid
        // (splitting cell c doesn't change cells 0..c-1)
        std::sort(intersected.rbegin(), intersected.rend());

        int splits_this_level = 0;
        for (const int c : intersected)
        {
            if (ray_split_one_tet(mesh, c, ls, level_set_id, tol))
                ++splits_this_level;
        }

        total_splits += splits_this_level;
        if (splits_this_level == 0)
            break;

        // Evaluate phi at any new vertices that have phi = 0 already set
        // (from the root computation).  For vertices that were inherited
        // from the parent, their phi is already correct.  We don't need
        // to re-evaluate because ray_split_one_tet already handles this.
    }

    return total_splits;
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template int pick_ray_apex<float>(std::span<const float>);
template int pick_ray_apex<double>(std::span<const double>);

template bool ray_split_one_tet<float>(
    LocalMesh<float>&, int, const LevelSetFunction<float>&, int, float);
template bool ray_split_one_tet<double>(
    LocalMesh<double>&, int, const LevelSetFunction<double>&, int, double);

template int ray_refine_local_mesh<float>(
    LocalMesh<float>&, const LevelSetFunction<float>&, int, int, float);
template int ray_refine_local_mesh<double>(
    LocalMesh<double>&, const LevelSetFunction<double>&, int, int, double);

} // namespace cutcells
