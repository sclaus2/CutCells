// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <concepts>
#include <span>
#include <vector>

#include "cell_types.h"
#include "edge_root.h"
#include "iso_refine.h"
#include "level_set.h"

namespace cutcells
{

// ============================================================================
// EdgeState
// ============================================================================

/// Classification of a single edge with respect to one level set.
enum class EdgeState : uint8_t
{
    no_root = 0,         ///< Certified absence of a zero on the edge.
    one_root = 1,        ///< Certified existence of exactly one zero on the edge.
    multiple_roots = 2,  ///< Evidence that more than one root may exist.
    near_tangency = 3,   ///< Edge likely contains a tangential/near tangential contact.
    uncertain = 4        ///< Current tests cannot certify this edge.
};

// ============================================================================
// LocalMesh  — flat SoA local mesh for one background cell
// ============================================================================

/// Flat structure-of-arrays mesh local to a single background cell.
/// Holds reference and physical coordinates, edge and cell connectivity,
/// level-set evaluation results, edge classification, and root vertices.
///
/// All free functions operate on this struct.  No virtual dispatch.
/// Designed for easy Python / nanobind export (contiguous flat arrays).
template <std::floating_point T>
struct LocalMesh
{
    // ---- geometry meta ----
    int gdim = 0;   ///< geometric dimension
    int tdim = 0;   ///< topological dimension of cells in this mesh

    // ---- parent information ----
    int32_t parent_cell_id  = -1;   ///< which background cell this mesh belongs to
    cell::type parent_cell_type = cell::type::point; ///< parent background cell type
    int     n_level_sets    = 1;    ///< number of active level sets
    std::vector<T> parent_cell_coords_p1; ///< parent-cell P1 corner coordinates (size = n_p1_vertices * gdim)
    // (cell types are stored per cell in cell_types; see below)
    std::vector<int32_t> parent_vertex_to_local_vertex; ///< map: parent corner vertex id -> local vertex id
    std::vector<int32_t> parent_edge_to_local_edge;     ///< map: parent edge id -> local edge id

    // ---- vertices ----
    /// Physical coordinates, size = n_vertices * gdim.
    std::vector<T> vertex_x;

    /// Reference coordinates in the parent background cell, size = n_vertices * tdim.
    std::vector<T> vertex_ref_x;

    /// Parent background entity dimension per vertex:
    ///   0 = background corner vertex
    ///   1 = background edge (midpoint or root)
    ///   2 = background face midpoint or root
    ///   3 = interior of background cell
    ///  -1 = unknown
    std::vector<int32_t> vertex_parent_dim;

    /// Parent background entity id.
    /// For dim=0: local interpolation-node index in the parent cell.
    /// For dim=1: local edge index within the background cell.
    std::vector<int32_t> vertex_parent_id;

    /// If this vertex is a computed root on a local-mesh edge, store that local
    /// edge id; otherwise -1.
    std::vector<int32_t> vertex_root_edge_id;

    /// Bit mask (up to 64 level sets): bit i set if phi_i ≈ 0 at this vertex.
    std::vector<uint64_t> vertex_zero_mask;

    /// Bit mask (up to 64 level sets): bit i set if phi_i < 0 at this vertex.
    std::vector<uint64_t> vertex_inside_mask;

    /// Evaluated level-set values, size = n_vertices * n_level_sets.
    /// Layout: [phi_0(v0), phi_1(v0), ..., phi_0(v1), phi_1(v1), ...]
    std::vector<T> vertex_phi;

    // ---- edges ----
    /// Vertex connectivity (local indices), size = n_edges * 2.
    std::vector<int32_t> edge_vertices;

    /// Parent entity dimension per edge:
    ///   1  →  lies on a background boundary edge
    ///  -1  →  interior refined edge (no background parent)
    std::vector<int32_t> edge_parent_dim;

    /// Local background edge index for edges with parent_dim == 1; -1 otherwise.
    std::vector<int32_t> edge_parent_id;

    /// EdgeState per edge, stored as uint8_t.
    std::vector<uint8_t> edge_state;

    /// Index of the root vertex in this mesh for each edge (-1 if not computed).
    std::vector<int32_t> edge_root_vertex;

    /// Parameter t in [0, 1] of the computed edge root (-1 if unset).
    std::vector<T> edge_root_parameter;

    /// Root-solver diagnostics per edge.
    std::vector<int32_t> edge_root_iterations;
    std::vector<int32_t> edge_root_evaluations;
    std::vector<uint8_t> edge_root_converged;
    std::vector<T> edge_root_residual;

    // ---- cells ----
    /// Flattened cell vertex connectivity in CSR layout.
    /// Use cell_offsets to index into cell_vertices.
    std::vector<int32_t> cell_vertices;

    /// CSR offsets, size = n_cells + 1.
    /// Cell i uses cell_vertices[cell_offsets[i] .. cell_offsets[i+1] - 1].
    std::vector<int32_t> cell_offsets;

    /// Cell type per cell — supports any mix of types.
    std::vector<cell::type> cell_types;

    /// For each cell: starting index of that cell's edges in cell_edges_flat.
    /// Size = n_cells + 1.
    std::vector<int32_t> cell_edge_offsets;

    /// Flattened per-cell local edge index lists (index into edge_vertices).
    std::vector<int32_t> cell_edges_flat;

    /// Domain classification per cell (cast from cell::domain).
    std::vector<uint8_t> cell_domain;

    // ---- convenience accessors ----
    int n_vertices() const
    {
        if (gdim <= 0) return 0;
        return static_cast<int>(vertex_x.size()) / gdim;
    }

    int n_edges() const
    {
        return static_cast<int>(edge_vertices.size()) / 2;
    }

    int n_cells() const
    {
        return cell_offsets.empty() ? 0 : static_cast<int>(cell_offsets.size()) - 1;
    }

    /// Number of vertices for cell i.
    int cell_num_vertices(int i) const
    {
        return cell_offsets[i + 1] - cell_offsets[i];
    }
};

// ============================================================================
// EdgeCache  — optional, cell-local cache of already-computed root vertices
// ============================================================================

/// Maps local edge indices to already-computed root vertex indices.
/// Purely an optimisation; not required for correctness.
struct EdgeCache
{
    std::vector<int32_t> key_edge_id;  ///< local edge index in LocalMesh
    std::vector<int32_t> vertex_id;    ///< corresponding root vertex index
};

// ============================================================================
// Free function declarations
// ============================================================================

/// Initialize a LocalMesh from a RefinementTemplate and one parent background cell.
///
/// After this call:
///   - vertex_x is filled from parent-cell interpolation-node coordinates.
///   - vertex_parent_dim / vertex_parent_id are copied from template parent entities.
///   - edge_vertices, edge_parent_dim/id are derived from cell connectivity.
///   - cell_vertices, cell_offsets, cell_types, cell_edge_offsets, cell_edges_flat are filled.
///   - All classification arrays are zero-initialised (unknown / unset).
///
/// Level-set values are NOT evaluated here.
///
/// @param mesh           output LocalMesh (will be fully reset)
/// @param tpl            refinement template
/// @param parent_cell_coords parent-cell interpolation-node coordinates
///                          (size = tpl.n_vertices * gdim, Basix ordering)
/// @param parent_cell_type parent background cell type
/// @param parent_cell_id user-provided parent-cell id
/// @param n_level_sets   number of active level sets to reserve space for
template <std::floating_point T>
void init_local_mesh_from_template(
    LocalMesh<T>&             mesh,
    const RefinementTemplate& tpl,
    std::span<const T>        parent_cell_coords,
    cell::type                parent_cell_type,
    int                       parent_cell_id,
    int                       n_level_sets = 1);

/// Initialize a LocalMesh from one parent background cell (no refinement yet).
///
/// `parent_cell_coords` may contain:
/// - P1 corner nodes only, or
/// - higher-order interpolation nodes in Basix ordering (p=2,3,4 currently).
///
/// The local mesh topology starts with a single cell (the parent cell).
/// If higher-order nodes are provided, they are still stored as local vertices
/// with parent-entity tags so edge/cell classification can use edge-interior
/// and cell-interior interpolation nodes.
template <std::floating_point T>
void init_local_mesh_from_cell(
    LocalMesh<T>&      mesh,
    std::span<const T> parent_cell_coords,
    cell::type         parent_cell_type,
    int                parent_cell_id,
    int                n_level_sets = 1);

/// Reinitialize a local mesh from a refinement template.
///
/// Reinitialize `mesh` from `tpl` using parent-cell metadata already stored in
/// `mesh` by `init_local_mesh_from_cell`.
template <std::floating_point T>
void refine_local_mesh_from_template(
    LocalMesh<T>&             mesh,
    const RefinementTemplate& tpl);

/// Evaluate all active level sets at every vertex and fill vertex_phi,
/// vertex_zero_mask, and vertex_inside_mask.
///
/// For all vertices, value_fn is called at physical coordinates.
/// (No global-node mapping is assumed in the local-mesh API.)
///
/// @param mesh  LocalMesh (vertex_x must already be filled)
/// @param phi   level-set functions, size >= n_ls
/// @param n_ls  number of active level sets to evaluate
/// @param tol   tolerance for the zero-level detection
template <std::floating_point T>
void evaluate_levelsets_on_vertices(
    LocalMesh<T>&                           mesh,
    const std::vector<LevelSetFunction<T>>& phi,
    int                                     n_ls,
    T                                       tol = static_cast<T>(1e-14));

/// Classify every edge of the local mesh using the sign information already stored
/// in vertex_inside_mask and vertex_zero_mask.
///
/// Fills edge_state.  Only point-sign tests are performed in this version;
/// Bernstein or interval-arithmetic upgrades can be added later.
///
/// @param mesh          LocalMesh (vertex masks must already be filled)
/// @param level_set_id  which bit in the masks to test (default 0)
template <std::floating_point T>
void classify_local_edges(
    LocalMesh<T>& mesh,
    int           level_set_id = 0);

/// Compute the root position on one edge by linear interpolation and append
/// the result as a new vertex in the LocalMesh.
///
/// Requires edge_state[edge_id] == EdgeState::one_root and vertex_phi filled.
/// Updates edge_root_vertex[edge_id] with the index of the new vertex.
///
/// @returns local vertex index of the newly appended root vertex
template <std::floating_point T>
int compute_edge_root_linear(
    LocalMesh<T>& mesh,
    int           edge_id,
    int           level_set_id = 0);

/// Compute one edge root using selected root finder.
///
/// - `linear`: uses endpoint nodal values (fast path).
/// - `brent` / `itp` / `newton`: require `level_set.value(...)` and start from the linear guess.
template <std::floating_point T, std::integral I = int>
int compute_edge_root(
    LocalMesh<T>&                      mesh,
    const LevelSetFunction<T, I>&      level_set,
    int                                edge_id,
    int                                level_set_id = 0,
    cell::edge_root::method            root_method = cell::edge_root::method::linear);

/// Call compute_edge_root_linear for every edge with state == one_root
/// that does not yet have a root vertex computed.
template <std::floating_point T>
void compute_all_roots_linear(
    LocalMesh<T>& mesh,
    int           level_set_id = 0);

/// Compute roots on all `one_root` edges with the selected root finder.
template <std::floating_point T, std::integral I = int>
void compute_all_roots(
    LocalMesh<T>&                      mesh,
    const LevelSetFunction<T, I>&      level_set,
    int                                level_set_id = 0,
    cell::edge_root::method            root_method = cell::edge_root::method::linear);

/// Decompose all intersected local cells into both `phi<0` and `phi>0` fragments
/// using the linear cutting path.
///
/// This routine computes missing edge roots once (via edge_root_vertex cache),
/// then reuses those root vertices when inserting cut fragments back into the
/// local mesh.
template <std::floating_point T>
void decompose_local_mesh_linear(
    LocalMesh<T>& mesh,
    int           level_set_id = 0,
    bool          triangulate = true);

/// Decompose using cached roots after computing roots with selected method.
template <std::floating_point T, std::integral I = int>
void decompose_local_mesh(
    LocalMesh<T>&                      mesh,
    const LevelSetFunction<T, I>&      level_set,
    int                                level_set_id = 0,
    cell::edge_root::method            root_method = cell::edge_root::method::linear,
    bool                               triangulate = true);

} // namespace cutcells
