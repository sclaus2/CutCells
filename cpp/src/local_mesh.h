// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <concepts>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

#include "bernstein_backend.h"
#include "cell_types.h"
#include "cut_cell.h"
#include "edge_root.h"
#include "iso_refine.h"
#include "level_set.h"
#include "mapping.h"

namespace cutcells
{

enum class CurveStatus : uint8_t
{
    accepted = 0,
    fallback_straight = 1,
    ambiguous = 2,
    failed = 3
};

template <std::floating_point T, std::integral I>
struct LocalLevelSetFunction;

// ============================================================================
// EdgeState
// ============================================================================

/// Classification of a single edge with respect to one level set.
enum class EdgeState : uint8_t
{
    uncut = 0,         ///< Certified absence of a crossing (touch-only allowed).
    single_cross = 1,  ///< Certified exactly one crossing root.
    multi_cross = 2,   ///< More than one crossing/zero-interval may exist.
    zero_edge = 3,     ///< Edge lies on the interface.
    uncertain = 4,     ///< Current tests cannot certify this edge yet.

    // Legacy aliases kept for compatibility in existing call sites.
    no_root = uncut,
    one_root = single_cross,
    multiple_roots = multi_cross,
    near_tangency = multi_cross,
    vertex_zero_only = 5,
    on_interface = zero_edge
};

/// Numerical certification provenance for an edge state.
enum class EdgeCertification : uint8_t
{
    certified = 0,
    resolved_by_sampling = 1,
    resolved_by_subdivision = 2,
    unresolved = 3
};

/// Provenance of a cut-subcell edge during LUT certification.
enum class EdgeOrigin : uint8_t
{
    original = 0,
    root_split = 1,
    lut_new = 2
};

template <std::floating_point T, std::integral I = int>
struct LUTCertificationResult
{
    bool certified = true;
    bool hit_max_depth = false;
    int refine_iterations = 0;
    int invalid_cells = 0;
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
    /// Edge certification status per edge, stored as uint8_t (EdgeCertification).
    std::vector<uint8_t> edge_cert;

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

    // ---- faces (populated only when tdim >= 3) ----
    /// Flat vertex ids for all faces (CSR).
    std::vector<int32_t>    face_vertices;
    /// CSR offsets, size = n_faces + 1.
    std::vector<int32_t>    face_offsets;
    /// Cell type per face (triangle or quadrilateral).
    std::vector<cell::type> face_types;
    /// Parent entity dimension per face: 2 = lies on a background boundary face, -1 = interior.
    std::vector<int32_t>    face_parent_dim;
    /// Local background face id for faces with parent_dim == 2; -1 otherwise.
    std::vector<int32_t>    face_parent_id;
    /// Cell-to-face CSR offsets, size = n_cells + 1.
    std::vector<int32_t>    cell_face_offsets;
    /// Cell-to-face connectivity (indices into face arrays).
    std::vector<int32_t>    cell_faces_flat;
    /// Map: background face id -> local face id (-1 if not present).
    std::vector<int32_t>    parent_face_to_local_face;

    // ---- straight zero entities ----
    // Mixed-dimensional zero-level entities extracted from the decomposed local
    // mesh. Entities are filtered by (dimension, zero mask) when assembling
    // curves, surfaces, and intersection topology.
    std::vector<uint8_t>  zero_entity_dim;          ///< topological dim per entity
    std::vector<uint64_t> zero_entity_zero_mask;    ///< active zero level sets
    std::vector<uint64_t> zero_entity_sign_mask;    ///< reserved for side classification
    std::vector<int32_t>  zero_entity_vertices;     ///< CSR vertex ids
    std::vector<int32_t>  zero_entity_offsets;      ///< size = n_zero_entities + 1
    std::vector<int32_t>  zero_entity_parent_cell;  ///< owner sub-cell
    std::vector<int32_t>  zero_entity_parent_dim;   ///< shared background carrier dim or -1
    std::vector<int32_t>  zero_entity_parent_id;    ///< shared background carrier id or -1
    std::vector<uint8_t>  zero_entity_is_owned;     ///< 1 if this LocalMesh owns the entity for quadrature

    // 1D entities: explicit endpoint storage for chain construction.
    std::vector<int32_t>  zero_entity_endpoint_v0;
    std::vector<int32_t>  zero_entity_endpoint_v1;

    // 2D entities: boundary-edge topology stored as flat endpoint pairs.
    std::vector<int32_t>  zero_face_edge_vertices;  ///< flat pairs
    std::vector<int32_t>  zero_face_edge_offsets;   ///< size = n_zero_entities + 1

    // ---- curved zero-entity cache ----
    // Stores only interior high-order nodes in canonical entity-local ordering.
    std::vector<T>        curved_zero_ref_nodes;
    std::vector<int32_t>  curved_zero_offsets;      ///< size = n_zero_entities + 1
    std::vector<uint8_t>  curved_zero_converged;
    std::vector<uint8_t>  curved_zero_status;

    // ---- topology caches ----
    std::vector<int32_t>  zero_chain_offsets;       ///< size = n_zero_chains + 1
    std::vector<int32_t>  zero_chain_entity_ids;
    std::vector<uint8_t>  zero_chain_entity_reversed;
    std::vector<uint8_t>  zero_chain_is_closed;
    std::vector<uint64_t> zero_chain_zero_mask;

    std::vector<int32_t>  zero_patch_offsets;       ///< size = n_zero_patches + 1
    std::vector<int32_t>  zero_patch_entity_ids;
    std::vector<uint8_t>  zero_patch_entity_oriented;
    std::vector<uint8_t>  zero_patch_is_closed;
    std::vector<uint64_t> zero_patch_zero_mask;

    // ---- compatibility interface views ----
    // Legacy single-level-set codim-1 interface arrays. These mirror the subset
    // of zero_entity_* with dim = tdim - 1 and a single-bit zero mask.
    std::vector<int32_t>  iface_vertices;
    std::vector<int32_t>  iface_offsets;
    std::vector<int32_t>  iface_level_set_id;
    std::vector<int32_t>  iface_parent_cell;

    std::vector<T>        curved_iface_ref_nodes;
    std::vector<int32_t>  curved_iface_offsets;
    std::vector<uint8_t>  curved_iface_converged;
    std::vector<uint8_t>  curved_iface_status;

    int curved_geometry_order = -1;  ///< -1 = cache not built
    int curved_fallback_count = 0;   ///< number of nodes that fell back to straight
    int curved_ambiguous_count = 0;  ///< number of nodes with two-sided valid roots
    int curved_failed_count = 0;     ///< number of nodes where no valid projected root was accepted
    uint32_t zero_entity_version = 0;       ///< incremented whenever zero entities are rebuilt
    uint32_t curved_cache_version = 0;      ///< incremented whenever curved_zero_* cache is rebuilt
    uint32_t curved_cache_zero_version = 0; ///< zero_entity_version captured by latest curved cache

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

    /// Number of faces.
    int n_faces() const
    {
        return face_offsets.empty() ? 0 : static_cast<int>(face_offsets.size()) - 1;
    }

    /// Number of stored zero-level entities across all dimensions.
    int n_zero_entities() const
    {
        return zero_entity_offsets.empty()
            ? 0
            : static_cast<int>(zero_entity_offsets.size()) - 1;
    }

    int n_zero_chains() const
    {
        return zero_chain_offsets.empty()
            ? 0
            : static_cast<int>(zero_chain_offsets.size()) - 1;
    }

    int n_zero_patches() const
    {
        return zero_patch_offsets.empty()
            ? 0
            : static_cast<int>(zero_patch_offsets.size()) - 1;
    }

    /// Number of straight interface entities (legacy codim-1 compatibility view).
    int n_iface_entities() const
    {
        return iface_offsets.empty() ? 0 : static_cast<int>(iface_offsets.size()) - 1;
    }

    // ---- per-level-set edge array accessors ----
    /// Access edge_state for edge e and level set ls.
    uint8_t& edge_state_for(int e, int ls)
        { return edge_state[static_cast<std::size_t>(e * n_level_sets + ls)]; }
    const uint8_t& edge_state_for(int e, int ls) const
        { return edge_state[static_cast<std::size_t>(e * n_level_sets + ls)]; }

    /// Access edge_cert for edge e and level set ls.
    uint8_t& edge_cert_for(int e, int ls)
        { return edge_cert[static_cast<std::size_t>(e * n_level_sets + ls)]; }
    const uint8_t& edge_cert_for(int e, int ls) const
        { return edge_cert[static_cast<std::size_t>(e * n_level_sets + ls)]; }

    /// Access edge_root_vertex for edge e and level set ls.
    int32_t& edge_root_vertex_for(int e, int ls)
        { return edge_root_vertex[static_cast<std::size_t>(e * n_level_sets + ls)]; }
    const int32_t& edge_root_vertex_for(int e, int ls) const
        { return edge_root_vertex[static_cast<std::size_t>(e * n_level_sets + ls)]; }

    /// Access edge_root_parameter for edge e and level set ls.
    T& edge_root_parameter_for(int e, int ls)
        { return edge_root_parameter[static_cast<std::size_t>(e * n_level_sets + ls)]; }
    const T& edge_root_parameter_for(int e, int ls) const
        { return edge_root_parameter[static_cast<std::size_t>(e * n_level_sets + ls)]; }

    /// Access edge_root_iterations for edge e and level set ls.
    int32_t& edge_root_iterations_for(int e, int ls)
        { return edge_root_iterations[static_cast<std::size_t>(e * n_level_sets + ls)]; }

    /// Access edge_root_evaluations for edge e and level set ls.
    int32_t& edge_root_evaluations_for(int e, int ls)
        { return edge_root_evaluations[static_cast<std::size_t>(e * n_level_sets + ls)]; }

    /// Access edge_root_converged for edge e and level set ls.
    uint8_t& edge_root_converged_for(int e, int ls)
        { return edge_root_converged[static_cast<std::size_t>(e * n_level_sets + ls)]; }

    /// Access edge_root_residual for edge e and level set ls.
    T& edge_root_residual_for(int e, int ls)
        { return edge_root_residual[static_cast<std::size_t>(e * n_level_sets + ls)]; }
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

/// Red-refine the currently marked local cells in-place.
///
/// This refines the current local subcells rather than restarting from the
/// parent-cell iso-P1 template. New vertices inherit stable parent-entity
/// ancestry when they lie on original parent edges/faces; otherwise they are
/// marked as local interior vertices.
template <std::floating_point T>
void red_refine_marked_cells(
    LocalMesh<T>&             mesh,
    std::span<const uint8_t>  marked_cells);

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

/// Compute one edge root using a full-cell Bernstein representation of the
/// local FEM level-set polynomial restricted to the current local edge.
template <std::floating_point T, std::integral I = int>
int compute_edge_root_bernstein(
    LocalMesh<T>&                      mesh,
    const LocalLevelSetFunction<T, I>& level_set,
    int                                edge_id,
    int                                level_set_id = 0,
    cell::edge_root::method            root_method = cell::edge_root::method::itp,
    T                                  tol = static_cast<T>(1e-12));

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

/// Compute roots on all `one_root` edges using the selected local backend.
template <std::floating_point T, std::integral I = int>
void compute_all_roots_with_backend(
    LocalMesh<T>&                      mesh,
    const LevelSetFunction<T, I>&      level_set,
    LocalLevelSetBackend               backend,
    int                                level_set_id = 0,
    cell::edge_root::method            root_method = cell::edge_root::method::linear,
    T                                  tol = static_cast<T>(1e-12));

template <std::floating_point T, std::integral I = int>
bool certify_cut_subcells(
    const LevelSetFunction<T, I>&      level_set,
    const cell::CutCell<T>&            cut_cell,
    cell::type                         parent_cell_type,
    std::span<const T>                 parent_ls_values,
    int                                parent_cell_id,
    T                                  tol = static_cast<T>(1e-14),
    bool                               debug = false);

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

/// Build the face layer of a LocalMesh for tdim >= 3.
///
/// Enumerates all faces of all cells, deduplicates them, fills
/// face_vertices, face_offsets, face_types, face_parent_dim, face_parent_id,
/// cell_face_offsets, cell_faces_flat, and parent_face_to_local_face.
///
/// This is a no-op if mesh.tdim < 3.
template <std::floating_point T>
void build_local_faces(LocalMesh<T>& mesh);

/// Build or rebuild the edge layer of a local mesh (cell_edge_offsets,
/// cell_edges_flat, edge_vertices, n_edges_*).
template <std::floating_point T>
void build_local_edges(LocalMesh<T>& mesh);

/// Build or rebuild the face layer of a local mesh (face_vertices,
/// face_offsets, face_types, face_parent_*, cell_face_offsets,
/// cell_faces_flat).  No-op for tdim < 3.
template <std::floating_point T>
void build_local_faces(LocalMesh<T>& mesh);

/// Rebuild vertex_parent_*, edge_parent_*, and face_parent_* maps for
/// all vertices and edges (and faces when tdim == 3) of the local mesh.
template <std::floating_point T>
void rebuild_parent_entity_maps(LocalMesh<T>& mesh);

/// Build mixed-dimensional straight zero entities after decomposition.
///
/// The current extraction path populates codim-1 entities directly from the
/// decomposed sub-mesh and keeps the data shape-general and level-set-aware so
/// additional dim-1 and dim-0 intersection entities can be inserted later
/// without redesigning LocalMesh storage.
template <std::floating_point T>
void build_zero_entities(LocalMesh<T>& mesh, int level_set_id = 0);

/// Assign deterministic quadrature ownership for zero entities that lie on a
/// shared background carrier. Entities without a shared carrier remain owned by
/// their local mesh.
template <std::floating_point T>
void assign_zero_entity_ownership(
    std::span<LocalMesh<T>*> local_meshes,
    uint64_t                 zero_mask);

/// Legacy compatibility wrapper. Mirrors the codim-1 subset of zero entities
/// into iface_* arrays for existing callers.
template <std::floating_point T>
void build_interface_entities(LocalMesh<T>& mesh, int level_set_id = 0);

/// Build connected chains of dim-1 zero entities filtered by zero_mask.
template <std::floating_point T>
void build_zero_chains(LocalMesh<T>& mesh, uint64_t zero_mask);

/// Build connected patches of dim-2 zero entities filtered by zero_mask.
template <std::floating_point T>
void build_zero_patches(LocalMesh<T>& mesh, uint64_t zero_mask);

template <std::floating_point T, std::integral I = int>
LUTCertificationResult<T, I> decompose_local_mesh_with_backend(
    LocalMesh<T>&                      mesh,
    const LevelSetFunction<T, I>&      level_set,
    LocalLevelSetBackend               backend,
    int                                level_set_id = 0,
    cell::edge_root::method            root_method = cell::edge_root::method::linear,
    bool                               triangulate = true,
    int                                max_refine_depth = 6,
    T                                  tol = static_cast<T>(1e-14),
    bool                               debug = false);

/// Resolve all multi-root edges on the local mesh by recursively applying
/// green splits before LUT cutting.
///
/// For each edge classified as multi_cross or uncertain, this routine:
///   1. Identifies a separator parameter strictly between two consecutive
///      roots (using Bernstein subdivision for the polynomial backend,
///      or midpoint evaluation for the analytic backend).
///   2. Inserts a separator vertex on that edge.
///   3. Applies a shape-specific green refinement split to all cells
///      incident on that edge.
///   4. Reclassifies the affected child edges.
///   5. Repeats until every local edge has at most one root.
///
/// Returns true if all edges are resolved (at most one_root).
/// Returns false if max_iterations was hit with unresolved edges remaining.
template <std::floating_point T, std::integral I = int>
bool resolve_multi_root_edges(
    LocalMesh<T>&                      mesh,
    const LevelSetFunction<T, I>&      level_set,
    LocalLevelSetBackend               backend,
    int                                level_set_id = 0,
    T                                  tol = static_cast<T>(1e-14),
    int                                max_iterations = 20);

/// Curve all interface entities with backend dispatch.
///
/// - Bernstein backend: converts to BernsteinCell, evaluates in reference space.
/// - Analytical callbacks: uses LevelSetFunction value_fn/grad_fn with
///   reference-to-physical coordinate mapping and J^T gradient transform.
///
/// Must be called after build_interface_entities(mesh, level_set_id).
template <std::floating_point T, std::integral I = int>
void curve_interface_entities_with_backend(
    LocalMesh<T>&                      mesh,
    const LevelSetFunction<T, I>&      level_set,
    LocalLevelSetBackend               backend,
    int                                level_set_id = 0,
    int                                geom_order = 2,
    T                                  tol = static_cast<T>(1e-12));

template <std::floating_point T, std::integral I = int>
void curve_zero_entities_with_backend(
    LocalMesh<T>&                      mesh,
    const LevelSetFunction<T, I>&      level_set,
    LocalLevelSetBackend               backend,
    int                                level_set_id = 0,
    int                                geom_order = 2,
    T                                  tol = static_cast<T>(1e-12));

// ============================================================================
// Detail helpers — Gauss–Lobatto node computation
// ============================================================================

namespace detail
{

/// Evaluate Legendre polynomial P_n(x) and its first two derivatives
/// via the standard three-term recurrence.
template <std::floating_point T>
struct LegendreResult { T P; T dP; T d2P; };

template <std::floating_point T>
inline LegendreResult<T> legendre_eval(int n, T x)
{
    if (n == 0) return {T(1), T(0), T(0)};
    if (n == 1) return {x, T(1), T(0)};

    T P0 = T(1), P1 = x;
    T dP0 = T(0), dP1 = T(1);
    T d2P0 = T(0), d2P1 = T(0);

    for (int k = 2; k <= n; ++k)
    {
        const T c = static_cast<T>(k);
        const T a = (T(2) * c - T(1)) / c;
        const T b = (c - T(1)) / c;
        const T P2   = a * x * P1 - b * P0;
        const T dP2  = a * (P1 + x * dP1) - b * dP0;
        const T d2P2 = a * (T(2) * dP1 + x * d2P1) - b * d2P0;
        P0 = P1;   P1 = P2;
        dP0 = dP1; dP1 = dP2;
        d2P0 = d2P1; d2P1 = d2P2;
    }
    return {P1, dP1, d2P1};
}

/// Return p-1 interior Gauss–Lobatto points in (0, 1), excluding endpoints.
///
/// For Pk geometry order p the GL quadrature on [-1,1] has p+1 nodes:
///   {-1, roots of P'_p(x) on (-1,1), +1}.
/// This function returns the p-1 roots of P'_p(x), mapped to [0, 1].
template <std::floating_point T>
inline std::vector<T> gauss_lobatto_interior_points_1d(int geom_order)
{
    const int n_interior = geom_order - 1;
    if (n_interior <= 0)
        return {};

    std::vector<T> pts(static_cast<std::size_t>(n_interior));
    const int p = geom_order;

    for (int i = 0; i < n_interior; ++i)
    {
        // Chebyshev-based initial guess
        T xi = -std::cos(static_cast<T>(i + 1) * std::numbers::pi_v<T>
                         / static_cast<T>(p));

        // Newton iteration: solve P'_p(xi) = 0 using P''_p(xi)
        for (int iter = 0; iter < 100; ++iter)
        {
            const auto [P, dP, d2P] = legendre_eval<T>(p, xi);
            if (std::abs(d2P) < T(1e-30))
                break;
            const T delta = dP / d2P;
            xi -= delta;
            if (std::abs(delta) < T(1e-15))
                break;
        }
        // Map [-1, 1] → [0, 1]
        pts[static_cast<std::size_t>(i)] = (xi + T(1)) / T(2);
    }

    std::sort(pts.begin(), pts.end());
    return pts;
}

} // namespace detail

// ============================================================================
// Curved interface entity nodes
// ============================================================================

/// Curve all interface entities for one level set.
///
/// Must be called after build_interface_entities(mesh, level_set_id).
/// For each 2D interface entity (segment):
///   1. Generate geom_order - 1 interior Gauss-Lobatto nodes on the straight segment
///   2. Compute the segment normal as the projection direction
///   3. Project each interior node to phi = 0 via find_ray_root
///   4. Store result in curved_iface_ref_nodes cache
///
/// Endpoints (root vertices) are NOT stored — they are shared with
/// adjacent entities and already exist as mesh vertices.
///
/// eval_phi:  T(const T* x_ref, int tdim) -> T
/// eval_grad: void(const T* x_ref, int tdim, T* grad_out)
template <std::floating_point T, typename EvalPhi, typename EvalGrad>
inline void curve_interface_entities(
    LocalMesh<T>&        mesh,
    EvalPhi&&            eval_phi,
    EvalGrad&&           eval_grad,
    int                  level_set_id,
    int                  geom_order,
    T                    tol = static_cast<T>(1e-12))
{
    const int n_entities = mesh.n_iface_entities();
    const int tdim = mesh.tdim;
    const int n_interior = geom_order - 1;

    if (n_entities == 0 || n_interior <= 0)
    {
        mesh.curved_iface_offsets.assign(
            static_cast<std::size_t>(n_entities + 1), 0);
        mesh.curved_iface_ref_nodes.clear();
        mesh.curved_iface_converged.clear();
        mesh.curved_iface_status.clear();
        mesh.curved_fallback_count = 0;
        mesh.curved_ambiguous_count = 0;
        mesh.curved_failed_count = 0;
        mesh.curved_geometry_order = geom_order;
        return;
    }

    // ---- allocate CSR offsets ----
    mesh.curved_iface_offsets.resize(
        static_cast<std::size_t>(n_entities + 1), 0);
    mesh.curved_iface_offsets[0] = 0;

    // Count interior nodes per entity
    for (int ie = 0; ie < n_entities; ++ie)
    {
        if (mesh.iface_level_set_id[static_cast<std::size_t>(ie)]
            == static_cast<int32_t>(level_set_id))
        {
            mesh.curved_iface_offsets[static_cast<std::size_t>(ie + 1)]
                = static_cast<int32_t>(n_interior);
        }
        else
        {
            mesh.curved_iface_offsets[static_cast<std::size_t>(ie + 1)] = 0;
        }
    }

    // Exclusive prefix sum
    for (int i = 0; i < n_entities; ++i)
    {
        mesh.curved_iface_offsets[static_cast<std::size_t>(i + 1)] +=
            mesh.curved_iface_offsets[static_cast<std::size_t>(i)];
    }

    const int total_nodes =
        mesh.curved_iface_offsets[static_cast<std::size_t>(n_entities)];
    mesh.curved_iface_ref_nodes.resize(
        static_cast<std::size_t>(total_nodes * tdim), T(0));
    mesh.curved_iface_converged.resize(
        static_cast<std::size_t>(total_nodes), 0);
    mesh.curved_iface_status.assign(
        static_cast<std::size_t>(total_nodes),
        static_cast<uint8_t>(CurveStatus::failed));
    mesh.curved_fallback_count = 0;
    mesh.curved_ambiguous_count = 0;
    mesh.curved_failed_count = 0;

    if (tdim == 2)
    {
        // 2D: each interface entity is a segment (2 vertices)
        const auto gl_pts = detail::gauss_lobatto_interior_points_1d<T>(geom_order);

        std::vector<T> x_s(static_cast<std::size_t>(tdim));
        std::vector<T> d_vec(static_cast<std::size_t>(tdim));
        std::vector<T> x_s_phys(static_cast<std::size_t>(mesh.gdim), T(0));
        std::vector<T> x_proj_phys(static_cast<std::size_t>(mesh.gdim), T(0));

        // Use the bounding-box cell type for domain checks inside find_ray_root.
        // This prevents probes near simplex diagonal edges from being pruned
        // when the normal direction points slightly outside the parent simplex.
        const cell::type bbox_type =
            (mesh.parent_cell_type == cell::type::triangle
             || mesh.parent_cell_type == cell::type::quadrilateral)
                ? cell::type::quadrilateral
                : cell::type::hexahedron;
        const auto J_cols = cell::jacobian_col_indices(mesh.parent_cell_type);
        const int n_p1 = cell::get_num_vertices(mesh.parent_cell_type);
        const bool has_parent_geom =
            static_cast<int>(mesh.parent_cell_coords_p1.size()) == n_p1 * mesh.gdim;

        auto ref_to_phys = [&](const T* x_ref, T* x_phys)
        {
            // Fallback: if parent geometry is unavailable, stay in reference
            // coordinates for diagnostics/ordering checks.
            if (!has_parent_geom)
            {
                for (int d = 0; d < mesh.gdim; ++d)
                    x_phys[d] = (d < tdim) ? x_ref[d] : T(0);
                return;
            }

            const T* x0 = mesh.parent_cell_coords_p1.data();
            for (int d = 0; d < mesh.gdim; ++d)
                x_phys[d] = x0[d];

            for (int k = 0; k < tdim; ++k)
            {
                const int vk = J_cols[static_cast<std::size_t>(k)];
                if (vk < 0 || vk >= n_p1)
                    continue;
                const T* xk = &mesh.parent_cell_coords_p1[static_cast<std::size_t>(vk * mesh.gdim)];
                const T xi = x_ref[static_cast<std::size_t>(k)];
                for (int d = 0; d < mesh.gdim; ++d)
                    x_phys[d] += xi * (xk[d] - x0[d]);
            }
        };

        for (int ie = 0; ie < n_entities; ++ie)
        {
            if (mesh.iface_level_set_id[static_cast<std::size_t>(ie)]
                != static_cast<int32_t>(level_set_id))
                continue;

            const int v0_off = mesh.iface_offsets[static_cast<std::size_t>(ie)];
            const int v0_id = mesh.iface_vertices[static_cast<std::size_t>(v0_off)];
            const int v1_id = mesh.iface_vertices[static_cast<std::size_t>(v0_off + 1)];

            const T* r0 = &mesh.vertex_ref_x[static_cast<std::size_t>(v0_id * tdim)];
            const T* r1 = &mesh.vertex_ref_x[static_cast<std::size_t>(v1_id * tdim)];

            // Segment tangent and length
            T h_local = T(0);
            for (int k = 0; k < tdim; ++k)
            {
                const T dk = r1[k] - r0[k];
                h_local += dk * dk;
            }
            h_local = std::sqrt(h_local);

            const T* p0_phys = &mesh.vertex_x[static_cast<std::size_t>(v0_id * mesh.gdim)];
            const T* p1_phys = &mesh.vertex_x[static_cast<std::size_t>(v1_id * mesh.gdim)];
            T h_local_phys = T(0);
            for (int d = 0; d < mesh.gdim; ++d)
            {
                const T dk = p1_phys[d] - p0_phys[d];
                h_local_phys += dk * dk;
            }
            h_local_phys = std::sqrt(h_local_phys);

            // Segment normal in 2D: perpendicular to tangent.
            // This remains a fallback when the level-set gradient is unavailable.
            T nx, ny;
            if (h_local > tol)
            {
                const T tx = (r1[0] - r0[0]) / h_local;
                const T ty = (r1[1] - r0[1]) / h_local;
                nx = -ty;
                ny = tx;
            }
            else
            {
                // Degenerate segment: fall back to gradient direction
                nx = T(0);
                ny = T(0);
            }

            const int node_start =
                mesh.curved_iface_offsets[static_cast<std::size_t>(ie)];
            std::vector<T> seg_phys(static_cast<std::size_t>(mesh.gdim), T(0));
            T seg_h2_phys = T(0);
            for (int d = 0; d < mesh.gdim; ++d)
            {
                seg_phys[static_cast<std::size_t>(d)] = p1_phys[d] - p0_phys[d];
                seg_h2_phys += seg_phys[static_cast<std::size_t>(d)]
                             * seg_phys[static_cast<std::size_t>(d)];
            }
            const T s_tol = T(1e-6);
            T prev_s = -std::numeric_limits<T>::max();

            cell::edge_root::RaySearchOptions<T> ray_opts;
            ray_opts.xtol = tol;
            ray_opts.ftol = tol;
            ray_opts.domain_tol = std::max<T>(tol, T(1e-10));
            ray_opts.max_probe_distance = T(2) * std::max(h_local, tol);

            for (int ni = 0; ni < n_interior; ++ni)
            {
                const T t_gl = gl_pts[static_cast<std::size_t>(ni)];
                const int node_idx = node_start + ni;
                T* out = &mesh.curved_iface_ref_nodes[
                    static_cast<std::size_t>(node_idx * tdim)];

                // Straight position on the segment at GL parameter
                x_s[0] = r0[0] + t_gl * (r1[0] - r0[0]);
                x_s[1] = r0[1] + t_gl * (r1[1] - r0[1]);

                // Prefer the level-set gradient direction at the straight point.
                // On coarse straight segments, the segment normal can intersect a
                // distant branch of the interface, while the gradient is the
                // correct local surface normal.
                std::vector<T> grad(static_cast<std::size_t>(tdim));
                eval_grad(x_s.data(), tdim, grad.data());
                T gnorm = T(0);
                for (int k = 0; k < tdim; ++k)
                    gnorm += grad[static_cast<std::size_t>(k)]
                           * grad[static_cast<std::size_t>(k)];
                gnorm = std::sqrt(gnorm);

                if (gnorm > tol)
                {
                    for (int k = 0; k < tdim; ++k)
                        d_vec[static_cast<std::size_t>(k)]
                            = grad[static_cast<std::size_t>(k)] / gnorm;
                }
                else if (h_local > tol)
                {
                    d_vec[0] = nx;
                    d_vec[1] = ny;
                }
                else
                {
                    // Fully degenerate: keep straight
                    for (int k = 0; k < tdim; ++k)
                        out[k] = x_s[static_cast<std::size_t>(k)];
                    mesh.curved_iface_converged[
                        static_cast<std::size_t>(node_idx)] = 0;
                    mesh.curved_iface_status[
                        static_cast<std::size_t>(node_idx)] = static_cast<uint8_t>(CurveStatus::failed);
                    ++mesh.curved_fallback_count;
                    ++mesh.curved_failed_count;
                    prev_s = t_gl;
                    continue;
                }

                // Solve phi(x_s + t * d) = 0 by searching both +d and -d and
                // selecting the nearest valid root.
                auto rr = cell::edge_root::find_ray_root_nearest<T>(
                    eval_phi, eval_grad,
                    std::span<const T>(x_s.data(),
                                       static_cast<std::size_t>(tdim)),
                    std::span<const T>(d_vec.data(),
                                       static_cast<std::size_t>(tdim)),
                    bbox_type,
                    tdim, h_local, ray_opts);
                rr.used_gradient_dir = (gnorm > tol);
                rr.used_fallback_dir = (gnorm <= tol);

                // Resolve two-sided ambiguity by nearest physical displacement.
                T dnorm = T(0);
                for (int k = 0; k < tdim; ++k)
                    dnorm += d_vec[static_cast<std::size_t>(k)]
                           * d_vec[static_cast<std::size_t>(k)];
                dnorm = std::sqrt(dnorm);
                if (dnorm > tol)
                {
                    for (int k = 0; k < tdim; ++k)
                        d_vec[static_cast<std::size_t>(k)] /= dnorm;
                }

                ref_to_phys(x_s.data(), x_s_phys.data());
                T chosen_dist_phys = std::numeric_limits<T>::max();
                bool have_chosen_candidate = false;
                std::array<T, 3> chosen_ref = {T(0), T(0), T(0)};

                auto choose_candidate = [&](T t_val)
                {
                    std::array<T, 3> cand_ref = {T(0), T(0), T(0)};
                    for (int k = 0; k < tdim; ++k)
                        cand_ref[static_cast<std::size_t>(k)]
                            = x_s[static_cast<std::size_t>(k)]
                            + t_val * d_vec[static_cast<std::size_t>(k)];
                    ref_to_phys(cand_ref.data(), x_proj_phys.data());

                    T dist2 = T(0);
                    for (int d = 0; d < mesh.gdim; ++d)
                    {
                        const T dd = x_proj_phys[static_cast<std::size_t>(d)]
                                   - x_s_phys[static_cast<std::size_t>(d)];
                        dist2 += dd * dd;
                    }
                    const T dist = std::sqrt(dist2);
                    if (!have_chosen_candidate || dist < chosen_dist_phys)
                    {
                        chosen_dist_phys = dist;
                        chosen_ref = cand_ref;
                        rr.t = t_val;
                        rr.dist = dist;
                        have_chosen_candidate = true;
                    }
                };

                if (rr.pos_valid)
                    choose_candidate(rr.pos_t);
                if (rr.neg_valid)
                    choose_candidate(rr.neg_t);

                if (have_chosen_candidate)
                {
                    for (int k = 0; k < tdim; ++k)
                        rr.x_ref[static_cast<std::size_t>(k)] = chosen_ref[static_cast<std::size_t>(k)];
                }

                bool accepted = rr.valid;
                const T max_dist = T(1.5) * std::max(h_local_phys, tol);
                if (accepted && (!have_chosen_candidate || rr.dist > max_dist))
                    accepted = false;

                if (rr.ambiguous_two_sided)
                    ++mesh.curved_ambiguous_count;

                if (!accepted)
                {
                    // Fallback projection: local Newton toward the nearest
                    // zero-level from x_s, then reuse standard acceptance checks.
                    std::array<T, 3> x_new = {x_s[0], x_s[1], T(0)};
                    bool newton_ok = false;
                    for (int it = 0; it < 8; ++it)
                    {
                        const T f = eval_phi(x_new.data(), tdim);
                        if (std::abs(f) <= std::max<T>(tol, T(1e-12)))
                        {
                            newton_ok = true;
                            break;
                        }
                        std::array<T, 3> g = {T(0), T(0), T(0)};
                        eval_grad(x_new.data(), tdim, g.data());
                        T g2 = T(0);
                        for (int k = 0; k < tdim; ++k)
                            g2 += g[static_cast<std::size_t>(k)] * g[static_cast<std::size_t>(k)];
                        if (g2 <= T(1e-20))
                            break;

                        T step = f / g2;
                        std::array<T, 3> x_trial = x_new;
                        for (int k = 0; k < tdim; ++k)
                            x_trial[static_cast<std::size_t>(k)]
                                = x_new[static_cast<std::size_t>(k)]
                                - step * g[static_cast<std::size_t>(k)];

                        for (int bt = 0; bt < 6; ++bt)
                        {
                            if (cell::edge_root::is_inside_reference_domain<T>(
                                    std::span<const T>(x_trial.data(), static_cast<std::size_t>(tdim)),
                                    bbox_type,
                                    ray_opts.domain_tol))
                                break;
                            step *= T(0.5);
                            for (int k = 0; k < tdim; ++k)
                                x_trial[static_cast<std::size_t>(k)]
                                    = x_new[static_cast<std::size_t>(k)]
                                    - step * g[static_cast<std::size_t>(k)];
                        }
                        x_new = x_trial;
                    }

                    if (newton_ok)
                    {
                        rr.valid = true;
                        for (int k = 0; k < tdim; ++k)
                            rr.x_ref[static_cast<std::size_t>(k)] = x_new[static_cast<std::size_t>(k)];
                        ref_to_phys(rr.x_ref.data(), x_proj_phys.data());
                        T dist2 = T(0);
                        for (int d = 0; d < mesh.gdim; ++d)
                        {
                            const T dd = x_proj_phys[static_cast<std::size_t>(d)]
                                       - x_s_phys[static_cast<std::size_t>(d)];
                            dist2 += dd * dd;
                        }
                        rr.dist = std::sqrt(dist2);
                        rr.converged = true;
                        accepted = (rr.dist <= max_dist);
                    }
                }

                if (accepted)
                {
                    // Reconstruct projected position from accepted root.
                    for (int k = 0; k < tdim; ++k)
                        out[k] = rr.x_ref[static_cast<std::size_t>(k)];

                    T s_curr = t_gl;
                    ref_to_phys(out, x_proj_phys.data());
                    if (seg_h2_phys > tol)
                    {
                        T s_num = T(0);
                        for (int d = 0; d < mesh.gdim; ++d)
                        {
                            const T dd = x_proj_phys[static_cast<std::size_t>(d)] - p0_phys[d];
                            s_num += dd * seg_phys[static_cast<std::size_t>(d)];
                        }
                        s_curr = s_num / seg_h2_phys;
                    }
                    const bool monotone_ok = (s_curr + s_tol >= prev_s)
                        && (s_curr >= -s_tol)
                        && (s_curr <= T(1) + s_tol);
                    if (!monotone_ok)
                    {
                        for (int k = 0; k < tdim; ++k)
                            out[k] = x_s[static_cast<std::size_t>(k)];
                        mesh.curved_iface_converged[
                            static_cast<std::size_t>(node_idx)] = 0;
                        mesh.curved_iface_status[
                            static_cast<std::size_t>(node_idx)] = static_cast<uint8_t>(CurveStatus::fallback_straight);
                        ++mesh.curved_fallback_count;
                        prev_s = t_gl;
                        continue;
                    }

                    mesh.curved_iface_converged[
                        static_cast<std::size_t>(node_idx)] = 1;
                    mesh.curved_iface_status[
                        static_cast<std::size_t>(node_idx)] =
                        static_cast<uint8_t>(rr.ambiguous_two_sided
                                                 ? CurveStatus::ambiguous
                                                 : CurveStatus::accepted);
                    prev_s = s_curr;
                }
                else
                {
                    // Fallback to straight position
                    for (int k = 0; k < tdim; ++k)
                        out[k] = x_s[static_cast<std::size_t>(k)];
                    mesh.curved_iface_converged[
                        static_cast<std::size_t>(node_idx)] = 0;
                    mesh.curved_iface_status[
                        static_cast<std::size_t>(node_idx)] =
                        static_cast<uint8_t>(rr.converged ? CurveStatus::fallback_straight
                                                          : CurveStatus::failed);
                    ++mesh.curved_fallback_count;
                    if (!rr.converged)
                        ++mesh.curved_failed_count;
                    prev_s = t_gl;
                }
            }
        }
    }
    else if (tdim == 3)
    {
        // 3D face curving: deferred — initialize empty cache
        // For now, curved_iface_ref_nodes stays all zeros and
        // converged flags stay 0 for any 3D entities.
    }

    mesh.curved_geometry_order = geom_order;
}

// ============================================================================
// curve_zero_entities_impl — project interior GL nodes for all zero entities
// ============================================================================

/// Project interior GL nodes onto the zero level set for all qualifying zero entities.
///
/// Iterates zero_entity_* directly (no iface_* intermediary).
/// Filters: zero_entity_dim[ze] == tdim - 1 AND (zero_entity_zero_mask[ze] >> level_set_id) & 1.
/// Writes curved node positions directly into curved_zero_ref_nodes / curved_zero_offsets /
/// curved_zero_converged / curved_zero_status indexed by zero entity id.
///
/// This is the canonical implementation used by curve_zero_entities_with_backend.
/// Eliminates the sync_curved_zero_compatibility_view count-mismatch failure mode.
///
/// eval_phi:  T(const T* x_ref, int tdim) -> T
/// eval_grad: void(const T* x_ref, int tdim, T* grad_out)
template <std::floating_point T, typename EvalPhi, typename EvalGrad>
inline void curve_zero_entities_impl(
    LocalMesh<T>&        mesh,
    EvalPhi&&            eval_phi,
    EvalGrad&&           eval_grad,
    int                  level_set_id,
    int                  geom_order,
    T                    tol = static_cast<T>(1e-12))
{
    const int n_zero   = mesh.n_zero_entities();
    const int tdim     = mesh.tdim;
    const int n_interior = geom_order - 1;
    const int codim1_dim = std::max(0, tdim - 1);

    // ---- (re)allocate curved_zero_offsets for every zero entity ----
    mesh.curved_zero_offsets.assign(
        static_cast<std::size_t>(n_zero + 1), 0);
    mesh.curved_zero_ref_nodes.clear();
    mesh.curved_zero_converged.clear();
    mesh.curved_zero_status.clear();
    mesh.curved_fallback_count  = 0;
    mesh.curved_ambiguous_count = 0;
    mesh.curved_failed_count    = 0;

    if (n_zero == 0 || n_interior <= 0)
    {
        mesh.curved_geometry_order = geom_order;
        ++mesh.curved_cache_version;
        mesh.curved_cache_zero_version = mesh.zero_entity_version;
        return;
    }

    // Count interior nodes per qualifying entity
    for (int ze = 0; ze < n_zero; ++ze)
    {
        if (mesh.zero_entity_dim[static_cast<std::size_t>(ze)] == codim1_dim
            && ((mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] >> level_set_id) & uint64_t(1)))
        {
            mesh.curved_zero_offsets[static_cast<std::size_t>(ze + 1)]
                = static_cast<int32_t>(n_interior);
        }
    }

    // Exclusive prefix sum
    for (int i = 0; i < n_zero; ++i)
    {
        mesh.curved_zero_offsets[static_cast<std::size_t>(i + 1)] +=
            mesh.curved_zero_offsets[static_cast<std::size_t>(i)];
    }

    const int total_nodes =
        mesh.curved_zero_offsets[static_cast<std::size_t>(n_zero)];
    mesh.curved_zero_ref_nodes.resize(
        static_cast<std::size_t>(total_nodes * tdim), T(0));
    mesh.curved_zero_converged.resize(
        static_cast<std::size_t>(total_nodes), 0);
    mesh.curved_zero_status.assign(
        static_cast<std::size_t>(total_nodes),
        static_cast<uint8_t>(CurveStatus::failed));

    if (tdim == 2)
    {
        // 2D: each qualifying zero entity is a segment (2 vertices)
        const auto gl_pts = detail::gauss_lobatto_interior_points_1d<T>(geom_order);

        std::vector<T> x_s(static_cast<std::size_t>(tdim));
        std::vector<T> d_vec(static_cast<std::size_t>(tdim));
        std::vector<T> x_s_phys(static_cast<std::size_t>(mesh.gdim), T(0));
        std::vector<T> x_proj_phys(static_cast<std::size_t>(mesh.gdim), T(0));

        const cell::type bbox_type =
            (mesh.parent_cell_type == cell::type::triangle
             || mesh.parent_cell_type == cell::type::quadrilateral)
                ? cell::type::quadrilateral
                : cell::type::hexahedron;
        const auto J_cols = cell::jacobian_col_indices(mesh.parent_cell_type);
        const int n_p1 = cell::get_num_vertices(mesh.parent_cell_type);
        const bool has_parent_geom =
            static_cast<int>(mesh.parent_cell_coords_p1.size()) == n_p1 * mesh.gdim;

        auto ref_to_phys = [&](const T* x_ref, T* x_phys)
        {
            if (!has_parent_geom)
            {
                for (int d = 0; d < mesh.gdim; ++d)
                    x_phys[d] = (d < tdim) ? x_ref[d] : T(0);
                return;
            }
            const T* x0 = mesh.parent_cell_coords_p1.data();
            for (int d = 0; d < mesh.gdim; ++d)
                x_phys[d] = x0[d];
            for (int k = 0; k < tdim; ++k)
            {
                const int vk = J_cols[static_cast<std::size_t>(k)];
                if (vk < 0 || vk >= n_p1)
                    continue;
                const T* xk = &mesh.parent_cell_coords_p1[
                    static_cast<std::size_t>(vk * mesh.gdim)];
                const T xi = x_ref[static_cast<std::size_t>(k)];
                for (int d = 0; d < mesh.gdim; ++d)
                    x_phys[d] += xi * (xk[d] - x0[d]);
            }
        };

        for (int ze = 0; ze < n_zero; ++ze)
        {
            if (mesh.zero_entity_dim[static_cast<std::size_t>(ze)] != codim1_dim)
                continue;
            if (!((mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] >> level_set_id)
                  & uint64_t(1)))
                continue;

            const int v0_id = mesh.zero_entity_endpoint_v0[static_cast<std::size_t>(ze)];
            const int v1_id = mesh.zero_entity_endpoint_v1[static_cast<std::size_t>(ze)];
            if (v0_id < 0 || v1_id < 0)
                continue;

            const T* r0 = &mesh.vertex_ref_x[static_cast<std::size_t>(v0_id * tdim)];
            const T* r1 = &mesh.vertex_ref_x[static_cast<std::size_t>(v1_id * tdim)];

            // Segment tangent length in reference space
            T h_local = T(0);
            for (int k = 0; k < tdim; ++k)
            {
                const T dk = r1[k] - r0[k];
                h_local += dk * dk;
            }
            h_local = std::sqrt(h_local);

            const T* p0_phys = &mesh.vertex_x[static_cast<std::size_t>(v0_id * mesh.gdim)];
            const T* p1_phys = &mesh.vertex_x[static_cast<std::size_t>(v1_id * mesh.gdim)];
            T h_local_phys = T(0);
            for (int d = 0; d < mesh.gdim; ++d)
            {
                const T dk = p1_phys[d] - p0_phys[d];
                h_local_phys += dk * dk;
            }
            h_local_phys = std::sqrt(h_local_phys);

            // Segment normal in reference 2D (fallback direction)
            T nx = T(0), ny = T(0);
            if (h_local > tol)
            {
                const T tx = (r1[0] - r0[0]) / h_local;
                const T ty = (r1[1] - r0[1]) / h_local;
                nx = -ty;
                ny =  tx;
            }

            const int node_start =
                mesh.curved_zero_offsets[static_cast<std::size_t>(ze)];
            std::vector<T> seg_phys(static_cast<std::size_t>(mesh.gdim), T(0));
            T seg_h2_phys = T(0);
            for (int d = 0; d < mesh.gdim; ++d)
            {
                seg_phys[static_cast<std::size_t>(d)] = p1_phys[d] - p0_phys[d];
                seg_h2_phys += seg_phys[static_cast<std::size_t>(d)]
                             * seg_phys[static_cast<std::size_t>(d)];
            }
            const T s_tol = T(1e-6);
            T prev_s = -std::numeric_limits<T>::max();

            cell::edge_root::RaySearchOptions<T> ray_opts;
            ray_opts.xtol = tol;
            ray_opts.ftol = tol;
            ray_opts.domain_tol = std::max<T>(tol, T(1e-10));
            ray_opts.max_probe_distance = T(2) * std::max(h_local, tol);

            for (int ni = 0; ni < n_interior; ++ni)
            {
                const T t_gl = gl_pts[static_cast<std::size_t>(ni)];
                const int node_idx = node_start + ni;
                T* out = &mesh.curved_zero_ref_nodes[
                    static_cast<std::size_t>(node_idx * tdim)];

                // Straight position on the segment at GL parameter
                x_s[0] = r0[0] + t_gl * (r1[0] - r0[0]);
                x_s[1] = r0[1] + t_gl * (r1[1] - r0[1]);

                // Prefer the level-set gradient direction at the straight point.
                std::vector<T> grad(static_cast<std::size_t>(tdim));
                eval_grad(x_s.data(), tdim, grad.data());
                T gnorm = T(0);
                for (int k = 0; k < tdim; ++k)
                    gnorm += grad[static_cast<std::size_t>(k)]
                           * grad[static_cast<std::size_t>(k)];
                gnorm = std::sqrt(gnorm);

                if (gnorm > tol)
                {
                    for (int k = 0; k < tdim; ++k)
                        d_vec[static_cast<std::size_t>(k)]
                            = grad[static_cast<std::size_t>(k)] / gnorm;
                }
                else if (h_local > tol)
                {
                    d_vec[0] = nx;
                    d_vec[1] = ny;
                }
                else
                {
                    // Fully degenerate: keep straight
                    for (int k = 0; k < tdim; ++k)
                        out[k] = x_s[static_cast<std::size_t>(k)];
                    mesh.curved_zero_converged[
                        static_cast<std::size_t>(node_idx)] = 0;
                    mesh.curved_zero_status[
                        static_cast<std::size_t>(node_idx)] =
                        static_cast<uint8_t>(CurveStatus::failed);
                    ++mesh.curved_fallback_count;
                    ++mesh.curved_failed_count;
                    prev_s = t_gl;
                    continue;
                }

                // Solve phi(x_s + t * d) = 0 searching both ±d directions
                auto rr = cell::edge_root::find_ray_root_nearest<T>(
                    eval_phi, eval_grad,
                    std::span<const T>(x_s.data(),
                                       static_cast<std::size_t>(tdim)),
                    std::span<const T>(d_vec.data(),
                                       static_cast<std::size_t>(tdim)),
                    bbox_type,
                    tdim, h_local, ray_opts);
                rr.used_gradient_dir = (gnorm > tol);
                rr.used_fallback_dir = (gnorm <= tol);

                // Resolve two-sided ambiguity by nearest physical displacement
                T dnorm = T(0);
                for (int k = 0; k < tdim; ++k)
                    dnorm += d_vec[static_cast<std::size_t>(k)]
                           * d_vec[static_cast<std::size_t>(k)];
                dnorm = std::sqrt(dnorm);
                if (dnorm > tol)
                {
                    for (int k = 0; k < tdim; ++k)
                        d_vec[static_cast<std::size_t>(k)] /= dnorm;
                }

                ref_to_phys(x_s.data(), x_s_phys.data());
                T chosen_dist_phys = std::numeric_limits<T>::max();
                bool have_chosen_candidate = false;
                std::array<T, 3> chosen_ref = {T(0), T(0), T(0)};

                auto choose_candidate = [&](T t_val)
                {
                    std::array<T, 3> cand_ref = {T(0), T(0), T(0)};
                    for (int k = 0; k < tdim; ++k)
                        cand_ref[static_cast<std::size_t>(k)]
                            = x_s[static_cast<std::size_t>(k)]
                            + t_val * d_vec[static_cast<std::size_t>(k)];
                    ref_to_phys(cand_ref.data(), x_proj_phys.data());

                    T dist2 = T(0);
                    for (int d = 0; d < mesh.gdim; ++d)
                    {
                        const T dd = x_proj_phys[static_cast<std::size_t>(d)]
                                   - x_s_phys[static_cast<std::size_t>(d)];
                        dist2 += dd * dd;
                    }
                    const T dist = std::sqrt(dist2);
                    if (!have_chosen_candidate || dist < chosen_dist_phys)
                    {
                        chosen_dist_phys = dist;
                        chosen_ref = cand_ref;
                        rr.t = t_val;
                        rr.dist = dist;
                        have_chosen_candidate = true;
                    }
                };

                if (rr.pos_valid)
                    choose_candidate(rr.pos_t);
                if (rr.neg_valid)
                    choose_candidate(rr.neg_t);

                if (have_chosen_candidate)
                {
                    for (int k = 0; k < tdim; ++k)
                        rr.x_ref[static_cast<std::size_t>(k)] =
                            chosen_ref[static_cast<std::size_t>(k)];
                }

                bool accepted = rr.valid;
                const T max_dist = T(1.5) * std::max(h_local_phys, tol);
                if (accepted && (!have_chosen_candidate || rr.dist > max_dist))
                    accepted = false;

                if (rr.ambiguous_two_sided)
                    ++mesh.curved_ambiguous_count;

                if (!accepted)
                {
                    // Fallback Newton projection toward the nearest zero level
                    std::array<T, 3> x_new = {x_s[0], x_s[1], T(0)};
                    bool newton_ok = false;
                    for (int it = 0; it < 8; ++it)
                    {
                        const T f = eval_phi(x_new.data(), tdim);
                        if (std::abs(f) <= std::max<T>(tol, T(1e-12)))
                        {
                            newton_ok = true;
                            break;
                        }
                        std::array<T, 3> g = {T(0), T(0), T(0)};
                        eval_grad(x_new.data(), tdim, g.data());
                        T g2 = T(0);
                        for (int k = 0; k < tdim; ++k)
                            g2 += g[static_cast<std::size_t>(k)]
                                * g[static_cast<std::size_t>(k)];
                        if (g2 <= T(1e-20))
                            break;

                        T step = f / g2;
                        std::array<T, 3> x_trial = x_new;
                        for (int k = 0; k < tdim; ++k)
                            x_trial[static_cast<std::size_t>(k)]
                                = x_new[static_cast<std::size_t>(k)]
                                - step * g[static_cast<std::size_t>(k)];

                        for (int bt = 0; bt < 6; ++bt)
                        {
                            if (cell::edge_root::is_inside_reference_domain<T>(
                                    std::span<const T>(x_trial.data(),
                                                       static_cast<std::size_t>(tdim)),
                                    bbox_type,
                                    ray_opts.domain_tol))
                                break;
                            step *= T(0.5);
                            for (int k = 0; k < tdim; ++k)
                                x_trial[static_cast<std::size_t>(k)]
                                    = x_new[static_cast<std::size_t>(k)]
                                    - step * g[static_cast<std::size_t>(k)];
                        }
                        x_new = x_trial;
                    }

                    if (newton_ok)
                    {
                        rr.valid = true;
                        for (int k = 0; k < tdim; ++k)
                            rr.x_ref[static_cast<std::size_t>(k)] =
                                x_new[static_cast<std::size_t>(k)];
                        ref_to_phys(rr.x_ref.data(), x_proj_phys.data());
                        T dist2 = T(0);
                        for (int d = 0; d < mesh.gdim; ++d)
                        {
                            const T dd = x_proj_phys[static_cast<std::size_t>(d)]
                                       - x_s_phys[static_cast<std::size_t>(d)];
                            dist2 += dd * dd;
                        }
                        rr.dist = std::sqrt(dist2);
                        rr.converged = true;
                        accepted = (rr.dist <= max_dist);
                    }
                }

                if (accepted)
                {
                    for (int k = 0; k < tdim; ++k)
                        out[k] = rr.x_ref[static_cast<std::size_t>(k)];

                    T s_curr = t_gl;
                    ref_to_phys(out, x_proj_phys.data());
                    if (seg_h2_phys > tol)
                    {
                        T s_num = T(0);
                        for (int d = 0; d < mesh.gdim; ++d)
                        {
                            const T dd = x_proj_phys[static_cast<std::size_t>(d)]
                                       - p0_phys[d];
                            s_num += dd * seg_phys[static_cast<std::size_t>(d)];
                        }
                        s_curr = s_num / seg_h2_phys;
                    }
                    const bool monotone_ok = (s_curr + s_tol >= prev_s)
                        && (s_curr >= -s_tol)
                        && (s_curr <= T(1) + s_tol);
                    if (!monotone_ok)
                    {
                        for (int k = 0; k < tdim; ++k)
                            out[k] = x_s[static_cast<std::size_t>(k)];
                        mesh.curved_zero_converged[
                            static_cast<std::size_t>(node_idx)] = 0;
                        mesh.curved_zero_status[
                            static_cast<std::size_t>(node_idx)] =
                            static_cast<uint8_t>(CurveStatus::fallback_straight);
                        ++mesh.curved_fallback_count;
                        prev_s = t_gl;
                        continue;
                    }

                    mesh.curved_zero_converged[
                        static_cast<std::size_t>(node_idx)] = 1;
                    mesh.curved_zero_status[
                        static_cast<std::size_t>(node_idx)] =
                        static_cast<uint8_t>(rr.ambiguous_two_sided
                                                 ? CurveStatus::ambiguous
                                                 : CurveStatus::accepted);
                    prev_s = s_curr;
                }
                else
                {
                    // Fallback to straight position
                    for (int k = 0; k < tdim; ++k)
                        out[k] = x_s[static_cast<std::size_t>(k)];
                    mesh.curved_zero_converged[
                        static_cast<std::size_t>(node_idx)] = 0;
                    mesh.curved_zero_status[
                        static_cast<std::size_t>(node_idx)] =
                        static_cast<uint8_t>(rr.converged
                                                 ? CurveStatus::fallback_straight
                                                 : CurveStatus::failed);
                    ++mesh.curved_fallback_count;
                    if (!rr.converged)
                        ++mesh.curved_failed_count;
                    prev_s = t_gl;
                }
            }
        }
    }
    else if (tdim == 3)
    {
        // 3D face curving: project interior face nodes onto the zero level set.
        //
        // For a triangular face of order p, the number of interior face nodes
        // is (p-1)(p-2)/2. For a quad face: (p-1)^2.
        //
        // Each interior node starts at its straight (bilinear/barycentric)
        // position within the face and is projected onto phi=0 via Newton
        // iteration along the gradient direction.

        // For the node count in the cache, we use:
        // - n_interior for 1D (edge) entities = geom_order - 1  (already handled by offset calculation above)
        // - n_interior for 2D (face) entities = depends on face shape

        // Re-do the offset calculation for 3D where dim-2 entities are faces
        // and may have varying interior node counts.
        // The initial offset calculation above used n_interior = geom_order - 1
        // which is correct for edge entities. For face entities we need to recompute.

        // Recompute offsets for qualifying face entities
        for (int ze = 0; ze < n_zero; ++ze)
        {
            const auto ze_dim = mesh.zero_entity_dim[static_cast<std::size_t>(ze)];
            if (ze_dim != codim1_dim) // codim1_dim = tdim - 1 = 2
                continue;
            if (!((mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] >> level_set_id)
                  & uint64_t(1)))
                continue;

            // Determine face shape from vertex count
            const int f0 = mesh.zero_entity_offsets[static_cast<std::size_t>(ze)];
            const int f1 = mesh.zero_entity_offsets[static_cast<std::size_t>(ze + 1)];
            const int n_face_verts = f1 - f0;

            int n_face_interior = 0;
            if (n_face_verts == 3) // triangle
                n_face_interior = (geom_order - 1) * (geom_order - 2) / 2;
            else if (n_face_verts == 4) // quad
                n_face_interior = (geom_order - 1) * (geom_order - 1);

            // Override the count set by the initial loop (which assumed n_interior)
            // We need to adjust: the initial loop set n_interior = geom_order - 1
            // but for faces we need n_face_interior.
            // Reconstruct the offset entry for this entity.
            mesh.curved_zero_offsets[static_cast<std::size_t>(ze + 1)]
                = mesh.curved_zero_offsets[static_cast<std::size_t>(ze)]
                + static_cast<int32_t>(n_face_interior);
        }

        // Fix up offsets for non-qualifying entities after the face entities
        // (the prefix sum needs to be re-done since we changed face entries)
        // Actually, let's redo the full prefix sum from scratch for correctness.
        mesh.curved_zero_offsets.assign(
            static_cast<std::size_t>(n_zero + 1), 0);
        for (int ze = 0; ze < n_zero; ++ze)
        {
            const auto ze_dim = mesh.zero_entity_dim[static_cast<std::size_t>(ze)];
            if (!((mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] >> level_set_id)
                  & uint64_t(1)))
            {
                mesh.curved_zero_offsets[static_cast<std::size_t>(ze + 1)] = 0;
                continue;
            }

            if (ze_dim == 1) // edge entity
            {
                mesh.curved_zero_offsets[static_cast<std::size_t>(ze + 1)]
                    = static_cast<int32_t>(n_interior);
            }
            else if (ze_dim == 2) // face entity
            {
                const int f0 = mesh.zero_entity_offsets[static_cast<std::size_t>(ze)];
                const int f1 = mesh.zero_entity_offsets[static_cast<std::size_t>(ze + 1)];
                const int n_face_verts = f1 - f0;

                int n_face_interior = 0;
                if (n_face_verts == 3)
                    n_face_interior = (geom_order - 1) * (geom_order - 2) / 2;
                else if (n_face_verts == 4)
                    n_face_interior = (geom_order - 1) * (geom_order - 1);

                mesh.curved_zero_offsets[static_cast<std::size_t>(ze + 1)]
                    = static_cast<int32_t>(n_face_interior);
            }
            else
            {
                mesh.curved_zero_offsets[static_cast<std::size_t>(ze + 1)] = 0;
            }
        }

        // Exclusive prefix sum
        for (int i = 0; i < n_zero; ++i)
        {
            mesh.curved_zero_offsets[static_cast<std::size_t>(i + 1)] +=
                mesh.curved_zero_offsets[static_cast<std::size_t>(i)];
        }

        const int total_nodes_3d =
            mesh.curved_zero_offsets[static_cast<std::size_t>(n_zero)];
        mesh.curved_zero_ref_nodes.resize(
            static_cast<std::size_t>(total_nodes_3d * tdim), T(0));
        mesh.curved_zero_converged.resize(
            static_cast<std::size_t>(total_nodes_3d), 0);
        mesh.curved_zero_status.assign(
            static_cast<std::size_t>(total_nodes_3d),
            static_cast<uint8_t>(CurveStatus::failed));

        // Bounding box type for domain checks
        const cell::type bbox_type = cell::type::hexahedron;

        // Project face interior nodes
        for (int ze = 0; ze < n_zero; ++ze)
        {
            if (mesh.zero_entity_dim[static_cast<std::size_t>(ze)] != 2)
                continue;
            if (!((mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] >> level_set_id)
                  & uint64_t(1)))
                continue;

            const int node_start =
                mesh.curved_zero_offsets[static_cast<std::size_t>(ze)];
            const int node_end =
                mesh.curved_zero_offsets[static_cast<std::size_t>(ze + 1)];
            const int n_face_nodes = node_end - node_start;
            if (n_face_nodes <= 0)
                continue;

            // Get face vertices in reference space
            const int f0 = mesh.zero_entity_offsets[static_cast<std::size_t>(ze)];
            const int f1 = mesh.zero_entity_offsets[static_cast<std::size_t>(ze + 1)];
            const int n_fv = f1 - f0;

            // Gather face vertex reference coordinates
            std::vector<T> face_ref(static_cast<std::size_t>(n_fv * tdim));
            for (int i = 0; i < n_fv; ++i)
            {
                const int vid = mesh.zero_entity_vertices[static_cast<std::size_t>(f0 + i)];
                for (int k = 0; k < tdim; ++k)
                    face_ref[static_cast<std::size_t>(i * tdim + k)]
                        = mesh.vertex_ref_x[static_cast<std::size_t>(vid * tdim + k)];
            }

            // Generate interior node positions on a regular grid in face param space
            // For triangles: use barycentric coordinates (L1, L2, L3) with
            //   L_i = i_k / p for i_k = 1..p-1, sum < p
            // For quads: tensor product of GL interior points

            int node_idx = 0;
            if (n_fv == 3) // triangular face
            {
                // Interior nodes with barycentric coordinates (i/p, j/p, k/p)
                // where i,j,k >= 1 and i+j+k = p
                const T dp = static_cast<T>(geom_order);
                for (int j = 1; j < geom_order; ++j)
                {
                    for (int i = 1; i < geom_order - j; ++i)
                    {
                        const T L0 = static_cast<T>(i) / dp;
                        const T L1 = static_cast<T>(j) / dp;
                        // L2 = 1 - L0 - L1
                        // Reference position = L0*v0 + L1*v1 + (1-L0-L1)*v2
                        // Actually: barycentric -> face coords
                        // x_ref = (1-L0-L1)*fv0 + L0*fv1 + L1*fv2
                        T* out = &mesh.curved_zero_ref_nodes[
                            static_cast<std::size_t>((node_start + node_idx) * tdim)];
                        for (int k = 0; k < tdim; ++k)
                        {
                            out[k] = (T(1) - L0 - L1) * face_ref[static_cast<std::size_t>(0 * tdim + k)]
                                   + L0 * face_ref[static_cast<std::size_t>(1 * tdim + k)]
                                   + L1 * face_ref[static_cast<std::size_t>(2 * tdim + k)];
                        }
                        ++node_idx;
                    }
                }
            }
            else if (n_fv == 4) // quad face
            {
                // Tensor product of equispaced interior points
                for (int j = 1; j < geom_order; ++j)
                {
                    const T v_param = static_cast<T>(j) / static_cast<T>(geom_order);
                    for (int i = 1; i < geom_order; ++i)
                    {
                        const T u_param = static_cast<T>(i) / static_cast<T>(geom_order);
                        T* out = &mesh.curved_zero_ref_nodes[
                            static_cast<std::size_t>((node_start + node_idx) * tdim)];
                        // Bilinear interpolation: (1-u)(1-v)*v0 + u(1-v)*v1 + uv*v2 + (1-u)v*v3
                        for (int k = 0; k < tdim; ++k)
                        {
                            out[k] = (T(1) - u_param) * (T(1) - v_param) * face_ref[static_cast<std::size_t>(0 * tdim + k)]
                                   + u_param * (T(1) - v_param) * face_ref[static_cast<std::size_t>(1 * tdim + k)]
                                   + u_param * v_param * face_ref[static_cast<std::size_t>(2 * tdim + k)]
                                   + (T(1) - u_param) * v_param * face_ref[static_cast<std::size_t>(3 * tdim + k)];
                        }
                        ++node_idx;
                    }
                }
            }

            // Now project each interior node onto phi=0 via Newton iteration
            for (int ni = 0; ni < n_face_nodes; ++ni)
            {
                const int abs_idx = node_start + ni;
                T* out = &mesh.curved_zero_ref_nodes[
                    static_cast<std::size_t>(abs_idx * tdim)];

                // Newton projection: x_{n+1} = x_n - (phi(x_n) / |grad|^2) * grad
                std::array<T, 3> x_cur = {T(0), T(0), T(0)};
                for (int k = 0; k < tdim; ++k)
                    x_cur[static_cast<std::size_t>(k)] = out[k];

                bool converged = false;
                for (int it = 0; it < 20; ++it)
                {
                    const T f = eval_phi(x_cur.data(), tdim);
                    if (std::abs(f) <= std::max<T>(tol, static_cast<T>(1e-12)))
                    {
                        converged = true;
                        break;
                    }

                    std::array<T, 3> g = {T(0), T(0), T(0)};
                    eval_grad(x_cur.data(), tdim, g.data());
                    T g2 = T(0);
                    for (int k = 0; k < tdim; ++k)
                        g2 += g[static_cast<std::size_t>(k)]
                            * g[static_cast<std::size_t>(k)];
                    if (g2 <= static_cast<T>(1e-20))
                        break;

                    T step = f / g2;

                    // Line search with backtracking for domain safety
                    std::array<T, 3> x_trial = x_cur;
                    for (int k = 0; k < tdim; ++k)
                        x_trial[static_cast<std::size_t>(k)] -= step * g[static_cast<std::size_t>(k)];

                    for (int bt = 0; bt < 6; ++bt)
                    {
                        if (cell::edge_root::is_inside_reference_domain<T>(
                                std::span<const T>(x_trial.data(),
                                                   static_cast<std::size_t>(tdim)),
                                bbox_type,
                                std::max<T>(tol, static_cast<T>(1e-10))))
                            break;
                        step *= T(0.5);
                        for (int k = 0; k < tdim; ++k)
                            x_trial[static_cast<std::size_t>(k)]
                                = x_cur[static_cast<std::size_t>(k)]
                                - step * g[static_cast<std::size_t>(k)];
                    }
                    x_cur = x_trial;
                }

                if (converged)
                {
                    for (int k = 0; k < tdim; ++k)
                        out[k] = x_cur[static_cast<std::size_t>(k)];
                    mesh.curved_zero_converged[static_cast<std::size_t>(abs_idx)] = 1;
                    mesh.curved_zero_status[static_cast<std::size_t>(abs_idx)]
                        = static_cast<uint8_t>(CurveStatus::accepted);
                }
                else
                {
                    // Keep straight position (already in out)
                    mesh.curved_zero_converged[static_cast<std::size_t>(abs_idx)] = 0;
                    mesh.curved_zero_status[static_cast<std::size_t>(abs_idx)]
                        = static_cast<uint8_t>(CurveStatus::failed);
                    ++mesh.curved_fallback_count;
                    ++mesh.curved_failed_count;
                }
            }
        }
    }

    mesh.curved_geometry_order = geom_order;
    ++mesh.curved_cache_version;
    mesh.curved_cache_zero_version = mesh.zero_entity_version;
}

} // namespace cutcells
