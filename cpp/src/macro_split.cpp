// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "macro_split.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <unordered_map>

namespace cutcells
{

// ============================================================================
// pick_topology1_apex
// ============================================================================

template <std::floating_point T>
int pick_topology1_apex(std::span<const T> phi)
{
    assert(phi.size() == 4);

    int n_pos = 0;
    int n_neg = 0;
    for (int i = 0; i < 4; ++i)
    {
        if (phi[static_cast<std::size_t>(i)] > T(0))
            ++n_pos;
        else if (phi[static_cast<std::size_t>(i)] < T(0))
            ++n_neg;
    }

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

    return -1; // not a 1-vs-3 pattern
}

// ============================================================================
// Internal helpers
// ============================================================================

namespace
{

/// Pack two vertex IDs into a single 64-bit key (order-independent).
uint64_t edge_key(int32_t a, int32_t b)
{
    const auto lo = static_cast<uint64_t>(std::min(a, b));
    const auto hi = static_cast<uint64_t>(std::max(a, b));
    return (lo << 32) | hi;
}

/// Implementation of the topology-1 tet macro split.
///
/// When midpoint_cache is non-null, midpoints on shared face edges
/// are deduplicated across sibling cells at the same depth level.
template <std::floating_point T>
bool macro_split_impl(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol,
    std::unordered_map<uint64_t, int32_t>* midpoint_cache)
{
    const int n_ls = mesh.n_level_sets;
    const int gdim = mesh.gdim;
    const int tdim = mesh.tdim;

    const int c0 = mesh.cell_offsets[static_cast<std::size_t>(cell_id)];
    const int c1 = mesh.cell_offsets[static_cast<std::size_t>(cell_id + 1)];
    if (c1 - c0 != 4)
        return false;

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

    const int apex_local = pick_topology1_apex<T>(
        std::span<const T>(phi_vals.data(), 4));
    if (apex_local < 0)
        return false;

    // Face vertex local indices (the 3 that are not apex), preserving order
    std::array<int, 3> fl;
    int fi = 0;
    for (int i = 0; i < 4; ++i)
        if (i != apex_local)
            fl[static_cast<std::size_t>(fi++)] = i;

    const int32_t v0  = vids[static_cast<std::size_t>(apex_local)];
    const int32_t fv0 = vids[static_cast<std::size_t>(fl[0])];
    const int32_t fv1 = vids[static_cast<std::size_t>(fl[1])];
    const int32_t fv2 = vids[static_cast<std::size_t>(fl[2])];
    const T phi_apex = phi_vals[static_cast<std::size_t>(apex_local)];

    // -----------------------------------------------------------------
    // Helper: get or create a midpoint vertex between two existing verts
    // -----------------------------------------------------------------
    auto get_or_create_midpoint = [&](int32_t va, int32_t vb,
                                      int parent_dim) -> int32_t
    {
        if (midpoint_cache)
        {
            const uint64_t key = edge_key(va, vb);
            auto it = midpoint_cache->find(key);
            if (it != midpoint_cache->end())
                return it->second;
        }

        const auto new_vid
            = static_cast<int32_t>(mesh.vertex_parent_dim.size());

        for (int d = 0; d < gdim; ++d)
        {
            const T xa
                = mesh.vertex_x[static_cast<std::size_t>(va * gdim + d)];
            const T xb
                = mesh.vertex_x[static_cast<std::size_t>(vb * gdim + d)];
            mesh.vertex_x.push_back(T(0.5) * (xa + xb));
        }
        for (int d = 0; d < tdim; ++d)
        {
            const T ra
                = mesh.vertex_ref_x[static_cast<std::size_t>(va * tdim + d)];
            const T rb
                = mesh.vertex_ref_x[static_cast<std::size_t>(vb * tdim + d)];
            mesh.vertex_ref_x.push_back(T(0.5) * (ra + rb));
        }

        mesh.vertex_parent_dim.push_back(parent_dim);
        mesh.vertex_parent_id.push_back(mesh.parent_cell_id);
        mesh.vertex_root_edge_id.push_back(-1);

        // Evaluate phi at the new vertex
        std::array<T, 3> x_eval = {T(0), T(0), T(0)};
        for (int d = 0; d < gdim; ++d)
            x_eval[static_cast<std::size_t>(d)]
                = mesh.vertex_x[static_cast<std::size_t>(new_vid * gdim + d)];

        const T phi_val = ls.value(x_eval.data(), mesh.parent_cell_id);

        const auto nv = static_cast<std::size_t>(new_vid + 1);
        mesh.vertex_zero_mask.resize(nv, 0);
        mesh.vertex_inside_mask.resize(nv, 0);
        mesh.vertex_phi.resize(nv * static_cast<std::size_t>(n_ls), T(0));
        mesh.vertex_phi[static_cast<std::size_t>(
            new_vid * n_ls + level_set_id)] = phi_val;

        if (std::abs(phi_val) <= tol)
            mesh.vertex_zero_mask[static_cast<std::size_t>(new_vid)]
                |= (uint64_t(1) << static_cast<unsigned>(level_set_id));
        if (phi_val < -tol)
            mesh.vertex_inside_mask[static_cast<std::size_t>(new_vid)]
                |= (uint64_t(1) << static_cast<unsigned>(level_set_id));

        if (midpoint_cache)
            (*midpoint_cache)[edge_key(va, vb)] = new_vid;

        return new_vid;
    };

    // -----------------------------------------------------------------
    // Helper: create a face-centroid vertex (never shared)
    // -----------------------------------------------------------------
    auto create_centroid = [&](int32_t va, int32_t vb, int32_t vc) -> int32_t
    {
        const auto new_vid
            = static_cast<int32_t>(mesh.vertex_parent_dim.size());
        const T third = T(1) / T(3);

        for (int d = 0; d < gdim; ++d)
        {
            const T xa
                = mesh.vertex_x[static_cast<std::size_t>(va * gdim + d)];
            const T xb
                = mesh.vertex_x[static_cast<std::size_t>(vb * gdim + d)];
            const T xc
                = mesh.vertex_x[static_cast<std::size_t>(vc * gdim + d)];
            mesh.vertex_x.push_back(third * (xa + xb + xc));
        }
        for (int d = 0; d < tdim; ++d)
        {
            const T ra
                = mesh.vertex_ref_x[static_cast<std::size_t>(va * tdim + d)];
            const T rb
                = mesh.vertex_ref_x[static_cast<std::size_t>(vb * tdim + d)];
            const T rc
                = mesh.vertex_ref_x[static_cast<std::size_t>(vc * tdim + d)];
            mesh.vertex_ref_x.push_back(third * (ra + rb + rc));
        }

        mesh.vertex_parent_dim.push_back(2); // face interior
        mesh.vertex_parent_id.push_back(mesh.parent_cell_id);
        mesh.vertex_root_edge_id.push_back(-1);

        std::array<T, 3> x_eval = {T(0), T(0), T(0)};
        for (int d = 0; d < gdim; ++d)
            x_eval[static_cast<std::size_t>(d)]
                = mesh.vertex_x[static_cast<std::size_t>(new_vid * gdim + d)];

        const T phi_val = ls.value(x_eval.data(), mesh.parent_cell_id);

        const auto nv = static_cast<std::size_t>(new_vid + 1);
        mesh.vertex_zero_mask.resize(nv, 0);
        mesh.vertex_inside_mask.resize(nv, 0);
        mesh.vertex_phi.resize(nv * static_cast<std::size_t>(n_ls), T(0));
        mesh.vertex_phi[static_cast<std::size_t>(
            new_vid * n_ls + level_set_id)] = phi_val;

        if (std::abs(phi_val) <= tol)
            mesh.vertex_zero_mask[static_cast<std::size_t>(new_vid)]
                |= (uint64_t(1) << static_cast<unsigned>(level_set_id));
        if (phi_val < -tol)
            mesh.vertex_inside_mask[static_cast<std::size_t>(new_vid)]
                |= (uint64_t(1) << static_cast<unsigned>(level_set_id));

        return new_vid;
    };

    // -----------------------------------------------------------------
    // Create the 4 new vertices: m01, m12, m20, cf
    // -----------------------------------------------------------------
    const int32_t m01_id = get_or_create_midpoint(fv0, fv1, 1);
    const int32_t m12_id = get_or_create_midpoint(fv1, fv2, 1);
    const int32_t m20_id = get_or_create_midpoint(fv2, fv0, 1);
    const int32_t cf_id  = create_centroid(fv0, fv1, fv2);

    // -----------------------------------------------------------------
    // Acceptance: all 4 spokes from apex must have a sign change
    // -----------------------------------------------------------------
    auto spoke_ok = [&](int32_t vid) -> bool
    {
        const T phi_new = mesh.vertex_phi[static_cast<std::size_t>(
            vid * n_ls + level_set_id)];
        return phi_apex * phi_new < T(0);
    };

    if (!spoke_ok(m01_id) || !spoke_ok(m12_id)
        || !spoke_ok(m20_id) || !spoke_ok(cf_id))
    {
        // Vertices already created but cell not split.  Orphaned vertices
        // are harmless — they don't affect later processing.
        return false;
    }

    // -----------------------------------------------------------------
    // Build 6 child tets (fan from apex through subdivided face)
    //
    //   T0 = (v0, fv0, m01, cf)
    //   T1 = (v0, m01, fv1, cf)
    //   T2 = (v0, fv1, m12, cf)
    //   T3 = (v0, m12, fv2, cf)
    //   T4 = (v0, fv2, m20, cf)
    //   T5 = (v0, m20, fv0, cf)
    // -----------------------------------------------------------------
    const std::array<std::array<int32_t, 4>, 6> children = {{
        {v0, fv0,   m01_id, cf_id},
        {v0, m01_id, fv1,   cf_id},
        {v0, fv1,   m12_id, cf_id},
        {v0, m12_id, fv2,   cf_id},
        {v0, fv2,   m20_id, cf_id},
        {v0, m20_id, fv0,   cf_id}
    }};

    // -----------------------------------------------------------------
    // Replace cell arrays
    // -----------------------------------------------------------------
    const int old_nc = mesh.n_cells();

    std::vector<int32_t>    new_cell_vertices;
    std::vector<int32_t>    new_cell_offsets = {0};
    std::vector<cell::type> new_cell_types;

    new_cell_vertices.reserve(mesh.cell_vertices.size() + 20);
    new_cell_offsets.reserve(static_cast<std::size_t>(old_nc + 5 + 1));
    new_cell_types.reserve(static_cast<std::size_t>(old_nc + 5));

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
            new_cell_types.push_back(
                mesh.cell_types[static_cast<std::size_t>(c)]);
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

    return true;
}

} // anonymous namespace

// ============================================================================
// macro_split_topology1_tet (public — no deduplication)
// ============================================================================

template <std::floating_point T>
bool macro_split_topology1_tet(
    LocalMesh<T>& mesh,
    int cell_id,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    T tol)
{
    return macro_split_impl(mesh, cell_id, ls, level_set_id, tol, nullptr);
}

// ============================================================================
// macro_refine_topology1 (public — with deduplication)
// ============================================================================

template <std::floating_point T>
int macro_refine_topology1(
    LocalMesh<T>& mesh,
    const LevelSetFunction<T>& ls,
    int level_set_id,
    int max_depth,
    T tol)
{
    const int n_ls = mesh.n_level_sets;
    int total_splits = 0;

    // Cache persists across depth levels so that vertices already
    // created on shared face edges are reused.
    std::unordered_map<uint64_t, int32_t> midpoint_cache;

    for (int depth = 0; depth < max_depth; ++depth)
    {
        // Collect all 1-vs-3 intersected tets at this level
        std::vector<int> candidates;
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

            int n_pos = 0;
            int n_neg = 0;
            for (int j = s0; j < s1; ++j)
            {
                const int32_t vid
                    = mesh.cell_vertices[static_cast<std::size_t>(j)];
                const T phi
                    = mesh.vertex_phi[static_cast<std::size_t>(
                          vid * n_ls + level_set_id)];
                if (phi > T(0))
                    ++n_pos;
                else if (phi < T(0))
                    ++n_neg;
            }

            if ((n_pos == 1 && n_neg == 3) || (n_pos == 3 && n_neg == 1))
                candidates.push_back(c);
        }

        if (candidates.empty())
            break;

        // Process in reverse order so that replacing a cell does
        // not invalidate the indices of earlier cells.
        std::sort(candidates.rbegin(), candidates.rend());

        int splits_this_level = 0;
        for (const int c : candidates)
        {
            if (macro_split_impl(
                    mesh, c, ls, level_set_id, tol, &midpoint_cache))
                ++splits_this_level;
        }

        total_splits += splits_this_level;
        if (splits_this_level == 0)
            break;
    }

    return total_splits;
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template int pick_topology1_apex<float>(std::span<const float>);
template int pick_topology1_apex<double>(std::span<const double>);

template bool macro_split_topology1_tet<float>(
    LocalMesh<float>&, int, const LevelSetFunction<float>&, int, float);
template bool macro_split_topology1_tet<double>(
    LocalMesh<double>&, int, const LevelSetFunction<double>&, int, double);

template int macro_refine_topology1<float>(
    LocalMesh<float>&, const LevelSetFunction<float>&, int, int, float);
template int macro_refine_topology1<double>(
    LocalMesh<double>&, const LevelSetFunction<double>&, int, int, double);

} // namespace cutcells
