// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include "cell_flags.h"
#include "cell_topology.h"
#include "cell_types.h"
#include "local_mesh.h"

#include <array>
#include <concepts>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace cutcells
{

// ============================================================================
// Parent-entity helpers
// ============================================================================

/// Test whether a local vertex lies on a specific background-cell edge.
///
/// A vertex lies on background edge bg_edge_id when:
///   - its parent dimension is 1 and its parent id equals bg_edge_id, OR
///   - its parent dimension is 0 and its parent id is one of the two
///     background-cell vertices that bound bg_edge_id.
template <std::floating_point T>
inline bool vertex_lies_on_background_edge(
    const LocalMesh<T>& mesh, int lv, int bg_edge_id)
{
    if (bg_edge_id < 0 || bg_edge_id >= cell::num_edges(mesh.parent_cell_type))
        return false;

    const int pdim = mesh.vertex_parent_dim[static_cast<std::size_t>(lv)];
    const int pid  = mesh.vertex_parent_id[static_cast<std::size_t>(lv)];
    if (pdim == 1 && pid == bg_edge_id)
        return true;
    if (pdim == 0)
    {
        const auto edge
            = cell::edges(mesh.parent_cell_type)[static_cast<std::size_t>(bg_edge_id)];
        return pid == edge[0] || pid == edge[1];
    }
    return false;
}

/// Infer the background-cell edge that both local vertices lv0 and lv1 lie on.
///
/// Returns {1, edge_id} when a shared background edge is found, or {-1, -1}
/// if no such edge exists.
template <std::floating_point T>
inline std::pair<int32_t, int32_t> infer_background_edge_parent(
    const LocalMesh<T>& mesh, int lv0, int lv1)
{
    const int n_parent_edges = cell::num_edges(mesh.parent_cell_type);
    for (int pe = 0; pe < n_parent_edges; ++pe)
    {
        if (vertex_lies_on_background_edge(mesh, lv0, pe)
            && vertex_lies_on_background_edge(mesh, lv1, pe))
        {
            return {1, pe};
        }
    }
    return {-1, -1};
}

// ============================================================================
// Per-cell-type green split helpers (topology only, no template dependency)
// ============================================================================

/// Split an interval cell (cv0, cv1) at a separator vertex mid_v on edge
/// (ev0, ev1).  Produces 2 child intervals.
/// Returns the number of children (2), or 0 if the edge is not part of this cell.
inline int green_split_interval(
    int32_t cv0, int32_t cv1,
    int32_t ev0, int32_t ev1, int32_t mid_v,
    std::vector<int32_t>&   out_verts,
    std::vector<int32_t>&   out_offsets,
    std::vector<cell::type>& out_types)
{
    if (!((cv0 == ev0 && cv1 == ev1) || (cv0 == ev1 && cv1 == ev0)))
        return 0;

    out_verts.push_back(ev0);
    out_verts.push_back(mid_v);
    out_offsets.push_back(static_cast<int32_t>(out_verts.size()));
    out_types.push_back(cell::type::interval);

    out_verts.push_back(mid_v);
    out_verts.push_back(ev1);
    out_offsets.push_back(static_cast<int32_t>(out_verts.size()));
    out_types.push_back(cell::type::interval);

    return 2;
}

/// Split a triangle cell (cv0, cv1, cv2) at a separator vertex mid_v on
/// edge (ev0, ev1).
///
/// Identifies the opposite vertex and produces 2 child triangles:
///   child 1: (ev0, mid_v, v_opp)
///   child 2: (mid_v, ev1, v_opp)
///
/// Returns 2 on success, 0 if the edge is not part of this triangle.
inline int green_split_triangle(
    int32_t cv0, int32_t cv1, int32_t cv2,
    int32_t ev0, int32_t ev1, int32_t mid_v,
    std::vector<int32_t>&   out_verts,
    std::vector<int32_t>&   out_offsets,
    std::vector<cell::type>& out_types)
{
    std::array<int32_t, 3> verts = {cv0, cv1, cv2};
    int idx_a = -1, idx_b = -1;
    for (int i = 0; i < 3; ++i)
    {
        if (verts[static_cast<std::size_t>(i)] == ev0) idx_a = i;
        if (verts[static_cast<std::size_t>(i)] == ev1) idx_b = i;
    }
    if (idx_a < 0 || idx_b < 0)
        return 0;

    const int     idx_opp = 3 - idx_a - idx_b; // 0+1+2 = 3
    const int32_t v_opp   = verts[static_cast<std::size_t>(idx_opp)];

    out_verts.push_back(ev0);
    out_verts.push_back(mid_v);
    out_verts.push_back(v_opp);
    out_offsets.push_back(static_cast<int32_t>(out_verts.size()));
    out_types.push_back(cell::type::triangle);

    out_verts.push_back(mid_v);
    out_verts.push_back(ev1);
    out_verts.push_back(v_opp);
    out_offsets.push_back(static_cast<int32_t>(out_verts.size()));
    out_types.push_back(cell::type::triangle);

    return 2;
}

/// Split a quadrilateral cell (cv0..cv3) at a separator vertex mid_v on
/// edge (ev0, ev1).
///
/// Basix quad ordering: v0(0,0), v1(1,0), v2(0,1), v3(1,1).
/// Perimeter order: 0 → 1 → 3 → 2.
///
/// Inserts mid_v into the perimeter ring to create a pentagon, then
/// fan-triangulates from ring[0] into 3 child triangles.
///
/// Returns 3 on success, 0 on failure.
inline int green_split_quadrilateral(
    int32_t cv0, int32_t cv1, int32_t cv2, int32_t cv3,
    int32_t ev0, int32_t ev1, int32_t mid_v,
    std::vector<int32_t>&   out_verts,
    std::vector<int32_t>&   out_offsets,
    std::vector<cell::type>& out_types)
{
    std::array<int32_t, 4> verts = {cv0, cv1, cv2, cv3};
    int idx_a = -1, idx_b = -1;
    for (int i = 0; i < 4; ++i)
    {
        if (verts[static_cast<std::size_t>(i)] == ev0) idx_a = i;
        if (verts[static_cast<std::size_t>(i)] == ev1) idx_b = i;
    }
    if (idx_a < 0 || idx_b < 0)
        return 0;

    // Perimeter order for Basix quad: 0 → 1 → 3 → 2
    static constexpr std::array<int, 4> perimeter = {0, 1, 3, 2};

    std::array<int32_t, 5> ring;
    int ring_count = 0;
    for (int pi = 0; pi < 4; ++pi)
    {
        const int32_t pv
            = verts[static_cast<std::size_t>(perimeter[static_cast<std::size_t>(pi)])];
        ring[static_cast<std::size_t>(ring_count++)] = pv;
        const int32_t next_pv
            = verts[static_cast<std::size_t>(
                perimeter[static_cast<std::size_t>((pi + 1) % 4)])];
        if ((pv == ev0 && next_pv == ev1) || (pv == ev1 && next_pv == ev0))
            ring[static_cast<std::size_t>(ring_count++)] = mid_v;
    }
    if (ring_count != 5)
        return 0;

    // Fan triangulation from ring[0]
    for (int i = 1; i < 4; ++i)
    {
        out_verts.push_back(ring[0]);
        out_verts.push_back(ring[static_cast<std::size_t>(i)]);
        out_verts.push_back(ring[static_cast<std::size_t>(i + 1)]);
        out_offsets.push_back(static_cast<int32_t>(out_verts.size()));
        out_types.push_back(cell::type::triangle);
    }

    return 3;
}

/// Split a tetrahedron cell (cv0..cv3) at a separator vertex mid_v on
/// edge (ev0, ev1).
///
/// Produces 2 child tets by replacing each split-edge endpoint in turn
/// with mid_v (preserving the parent vertex ordering and orientation):
///   child 1: parent tet with the occurrence of ev1 replaced by mid_v
///   child 2: parent tet with the occurrence of ev0 replaced by mid_v
///
/// Returns 2 on success, 0 if the edge is not part of this tet.
inline int green_split_tetrahedron(
    int32_t cv0, int32_t cv1, int32_t cv2, int32_t cv3,
    int32_t ev0, int32_t ev1, int32_t mid_v,
    std::vector<int32_t>&   out_verts,
    std::vector<int32_t>&   out_offsets,
    std::vector<cell::type>& out_types)
{
    std::array<int32_t, 4> verts = {cv0, cv1, cv2, cv3};
    int idx_a = -1, idx_b = -1;
    for (int i = 0; i < 4; ++i)
    {
        if (verts[static_cast<std::size_t>(i)] == ev0) idx_a = i;
        if (verts[static_cast<std::size_t>(i)] == ev1) idx_b = i;
    }
    if (idx_a < 0 || idx_b < 0)
        return 0;

    // Child 1: replace ev1's slot with mid_v
    for (int i = 0; i < 4; ++i)
        out_verts.push_back(i == idx_b ? mid_v : verts[static_cast<std::size_t>(i)]);
    out_offsets.push_back(static_cast<int32_t>(out_verts.size()));
    out_types.push_back(cell::type::tetrahedron);

    // Child 2: replace ev0's slot with mid_v
    for (int i = 0; i < 4; ++i)
        out_verts.push_back(i == idx_a ? mid_v : verts[static_cast<std::size_t>(i)]);
    out_offsets.push_back(static_cast<int32_t>(out_verts.size()));
    out_types.push_back(cell::type::tetrahedron);

    return 2;
}

/// Split a hexahedron cell at a separator vertex mid_v on edge (ev0, ev1).
///
/// Decomposes the hex into 5 tetrahedra using the Kuhn-type decomposition:
///
///   Basix hex ordering: v0(0,0,0) v1(1,0,0) v2(0,1,0) v3(1,1,0)
///                       v4(0,0,1) v5(1,0,1) v6(0,1,1) v7(1,1,1)
///
///   T0: (v0, v1, v3, v5)  — hex edges 0-1, 1-3, 1-5
///   T1: (v0, v3, v2, v6)  — hex edges 0-2, 2-3, 2-6
///   T2: (v0, v4, v5, v6)  — hex edges 0-4, 4-5, 4-6
///   T3: (v3, v6, v5, v7)  — hex edges 3-7, 5-7, 6-7
///   T4: (v0, v5, v3, v6)  — face diagonals only (never contains a hex edge)
///
/// Every hex edge is in exactly one of T0..T3.  The tet containing the
/// split edge is green-split (→2 children); the rest pass through unchanged.
///
/// Returns the total number of child tets (always 6 for a single-edge split),
/// or 0 if cell_verts does not have exactly 8 vertices.
inline int green_split_hexahedron(
    std::span<const int32_t> cell_verts,
    int32_t                  ev0, int32_t ev1, int32_t mid_v,
    std::vector<int32_t>&   out_verts,
    std::vector<int32_t>&   out_offsets,
    std::vector<cell::type>& out_types)
{
    if (cell_verts.size() != 8)
        return 0;

    auto v = [&](int i) -> int32_t
    { return cell_verts[static_cast<std::size_t>(i)]; };

    static constexpr std::array<std::array<int, 4>, 5> tet_decomp = {{
        {0, 1, 3, 5}, // T0
        {0, 3, 2, 6}, // T1
        {0, 4, 5, 6}, // T2
        {3, 6, 5, 7}, // T3
        {0, 5, 3, 6}, // T4 (central, no hex edges)
    }};

    int n_children = 0;
    for (int t = 0; t < 5; ++t)
    {
        const auto&           td = tet_decomp[static_cast<std::size_t>(t)];
        std::array<int32_t, 4> tv = {v(td[0]), v(td[1]), v(td[2]), v(td[3])};

        int idx_a = -1, idx_b = -1;
        for (int i = 0; i < 4; ++i)
        {
            if (tv[static_cast<std::size_t>(i)] == ev0) idx_a = i;
            if (tv[static_cast<std::size_t>(i)] == ev1) idx_b = i;
        }

        if (idx_a >= 0 && idx_b >= 0)
        {
            // Green-split this tet
            for (int i = 0; i < 4; ++i)
                out_verts.push_back(i == idx_b ? mid_v : tv[static_cast<std::size_t>(i)]);
            out_offsets.push_back(static_cast<int32_t>(out_verts.size()));
            out_types.push_back(cell::type::tetrahedron);

            for (int i = 0; i < 4; ++i)
                out_verts.push_back(i == idx_a ? mid_v : tv[static_cast<std::size_t>(i)]);
            out_offsets.push_back(static_cast<int32_t>(out_verts.size()));
            out_types.push_back(cell::type::tetrahedron);

            n_children += 2;
        }
        else
        {
            for (int i = 0; i < 4; ++i)
                out_verts.push_back(tv[static_cast<std::size_t>(i)]);
            out_offsets.push_back(static_cast<int32_t>(out_verts.size()));
            out_types.push_back(cell::type::tetrahedron);

            n_children += 1;
        }
    }

    return n_children;
}

// ============================================================================
// Template function declarations
// ============================================================================

/// Find all cells in mesh that are incident on the given edge.
/// @param mesh      local mesh
/// @param edge_id   index into mesh.edge_vertices
/// @return unsorted list of cell indices
template <std::floating_point T>
std::vector<int> cells_incident_on_edge(const LocalMesh<T>& mesh, int edge_id);

/// Insert a separator vertex at parameter t_sep on edge edge_id and
/// green-split all cells incident on that edge.
///
/// Steps:
///   1. Interpolate separator coordinates (reference + physical).
///   2. Infer parent-entity tracking for the new vertex.
///   3. Replace each incident cell with its shape-specific children.
///   4. Rebuild edge and face layers and parent-entity maps.
///
/// Returns true if all incident cells were split successfully.
/// Returns false when an incident cell is a prism or pyramid (not yet
/// implemented); in that case the mesh is left unchanged and the caller
/// should fall back to red refinement.
template <std::floating_point T>
bool green_split_one_edge(LocalMesh<T>& mesh, int edge_id, T t_sep);

/// Find a separator parameter t ∈ (0, 1) strictly between two consecutive
/// roots on a multi-cross edge, using recursive Bernstein subdivision.
///
/// @param edge_coeffs  1D Bernstein coefficients of the level-set restricted
///                     to the edge (size p+1 for degree p).
/// @param tol          zero tolerance for sign classification.
/// @return t where the level set is strictly nonzero between two root
///         intervals.  Falls back to 0.5 when no separator is found within
///         20 subdivision levels.
template <std::floating_point T>
T find_separator_on_edge_bernstein(std::span<const T> edge_coeffs, T tol);

/// Find a separator parameter for a multi-cross edge using level-set values
/// already stored in the local mesh (analytic / callable backend).
///
/// Evaluates the linear interpolant at t = 0.5, 0.25, 0.75, 0.375, 0.625
/// and returns the first parameter where |phi| > tol.  Falls back to 0.5.
template <std::floating_point T>
T find_separator_on_edge_midpoint(
    const LocalMesh<T>& mesh, int edge_id, int level_set_id, T tol);

} // namespace cutcells
