// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "cell_types.h"
#include "mesh_view.h"

namespace cutcells
{

// ---------------------------------------------------------------
// Tag enumerations for certification
// ---------------------------------------------------------------

/// Per-edge root classification for one level set.
enum class EdgeRootTag : std::uint8_t
{
    not_classified = 0,
    no_root        = 1,
    one_root       = 2,
    multiple_roots = 3,
    zero           = 4
};

/// Per-leaf-cell certification state for one level set.
enum class CellCertTag : std::uint8_t
{
    not_classified = 0,
    positive       = 1,
    negative       = 2,
    cut            = 3,
    zero           = 4,
    ambiguous      = 5,
    ready_to_cut   = 6
};

/// Per-face certification state for one level set.
/// A face is a 2D entity; its sign/cut state is analogous to a cell.
enum class FaceCertTag : std::uint8_t
{
    not_classified = 0,
    positive       = 1,
    negative       = 2,
    cut            = 3,
    zero           = 4,
    ambiguous      = 5
};

/// Why a leaf cell was created from an earlier leaf cell.
enum class CellRefinementReason : std::uint8_t
{
    none             = 0,
    green_edge       = 1,
    red_cell         = 2,
    graph_green_edge = 3,
    graph_red_cell   = 4,
    cut_level_set    = 5
};

/// Compact Sparse Row (CSR) adjacency map between entities.
///
/// Stores a variable-length list of adjacent entity indices for each source entity.
/// `offsets[i]` and `offsets[i+1]` delimit the slice of `indices` for entity `i`.
///
/// Used for:
///   - `entity_to_vertex[d]`:  maps entity i of dimension d to its vertex indices
///   - `connectivity[d0][d1]`: maps entities of dimension d0 to entities of dimension d1
struct EntityAdjacency
{
    /// offsets into `indices`; size = n_entities + 1.  offsets[0] == 0 always.
    std::vector<std::int32_t> offsets;

    /// Flat list of adjacent entity ids.
    std::vector<std::int32_t> indices;

    /// Number of source entities.
    std::int32_t size() const noexcept
    {
        return offsets.empty() ? 0 : static_cast<std::int32_t>(offsets.size()) - 1;
    }

    /// Return the adjacency list for source entity `i`.
    std::span<const std::int32_t> operator[](std::int32_t i) const
    {
        return {indices.data() + offsets[i],
                static_cast<std::size_t>(offsets[i + 1] - offsets[i])};
    }

    /// Return a mutable adjacency list for source entity `i`.
    std::span<std::int32_t> operator[](std::int32_t i)
    {
        return {indices.data() + offsets[i],
                static_cast<std::size_t>(offsets[i + 1] - offsets[i])};
    }
};

/// AdaptCell is the main local container for one background cell.
/// It is a dynamic mesh in reference coordinates: vertices carry the rich
/// semantic data (coordinates, level-set signs, provenance), and entities of
/// dimension 1, 2, 3 are defined purely through their vertex references.
/// Only leaf entities are stored — no refinement tree.
template <std::floating_point T>
struct AdaptCell
{
    // ---------------------------------------------------------------
    // Metadata
    // ---------------------------------------------------------------

    int gdim = 0;         ///< geometric dimension of the embedding space
    int tdim = 0;         ///< topological dimension (2 for triangles/quads, 3 for tets/hexes)

    cell::type parent_cell_type = cell::type::point;  ///< type of the original background cell
    int parent_cell_id = -1;                          ///< index of the background cell in the mesh

    ///< bit i set → level set i is active on this cell (cuts through it)
    std::uint64_t active_level_set_mask = 0;

    /// Per-leaf-cell active level-set mask.
    /// Parallel to entity_types[tdim]/entity_to_vertex[tdim]:
    /// bit i set ↔ level set i intersects that leaf cell.
    std::vector<std::uint64_t> cell_active_level_set_mask;

    /// Persistent leaf-cell provenance. These arrays are parallel to the
    /// top-dimensional entity pool and are preserved through repeated
    /// adaptcell refinements.
    std::vector<std::int32_t> cell_source_cell_id;
    std::vector<std::int32_t> cell_refinement_generation;
    std::vector<CellRefinementReason> cell_refinement_reason;
    std::vector<std::int32_t> cell_host_parent_cell_id;

    // ---------------------------------------------------------------
    // Vertices (dimension 0)
    //
    // Vertices are the primary data carriers. They are NOT stored in
    // the entity pools below. All arrays are parallel, indexed by
    // vertex id v = 0, 1, ..., n_vertices()-1.
    // ---------------------------------------------------------------

    /// Reference coordinates of each vertex, flattened: [v0_xi0, v0_xi1, ..., v1_xi0, ...].
    /// These are coordinates in the reference cell of the parent background cell.
    std::vector<T> vertex_coords;

    /// Tri-state sign storage per vertex and per level set.
    /// For level set i (0-indexed):
    ///   bit i set in zero_mask[v]     → φ_i(v) = 0
    ///   bit i set in negative_mask[v] → φ_i(v) < 0
    ///   neither                       → φ_i(v) > 0
    std::vector<std::uint64_t> zero_mask_per_vertex;
    std::vector<std::uint64_t> negative_mask_per_vertex;

    /// Provenance: where in the ORIGINAL parent background cell did this vertex come from?
    /// This is the key to global stitching without geometric deduplication.
    ///
    /// parent_entity_dim: 0 = parent vertex, 1 = parent edge, 2 = parent face, 3 = parent cell interior
    /// parent_entity_id:  GLOBAL mesh entity id (e.g. global edge 47, global vertex 12).
    ///                    Using global ids makes stitching across background cells trivial:
    ///                    two AdaptCells sharing an edge both reference the same global edge id,
    ///                    so matching vertices requires only (parent_entity_id, parent_param) comparison.
    /// parent_local_param: parametric coordinate(s) within the parent entity.
    ///                     For a vertex born on a parent edge: one parameter t ∈ [0,1].
    ///                     For a vertex born on a parent face: two parameters (u,v).
    ///                     For a vertex coinciding with a parent vertex: empty.
    ///                     For a vertex in the cell interior: tdim parameters.
    ///                     Stored flat; use parent_entity_dim to know how many params per vertex.
    std::vector<std::int8_t>  vertex_parent_dim;
    std::vector<std::int32_t> vertex_parent_id;
    std::vector<T>            vertex_parent_param;         ///< flat, variable-length per vertex
    std::vector<std::int32_t> vertex_parent_param_offset; ///< size = n_vertices + 1
    std::vector<std::int32_t> vertex_source_edge_id;      ///< leaf-edge provenance; -1 if not created from a leaf edge

    // ---------------------------------------------------------------
    // Entity pools (dimensions 1 .. tdim)
    //
    // Each entity of dimension d has:
    //   - a cell type (e.g. interval for d=1, triangle/quad for d=2, tet/hex for d=3)
    //   - an ordered list of vertex indices (via entity_to_vertex[d])
    //
    // Entities are numbered 0, 1, ..., n_entities(d)-1 within each dimension.
    // Only leaf (final, active) entities are stored.
    // ---------------------------------------------------------------

    static constexpr int max_dim = 4;

    /// Cell type of each entity. entity_types[d][i] is the type of entity i at dimension d.
    /// For d=1 this is always cell::type::interval, but stored for uniformity.
    std::array<std::vector<cell::type>, max_dim> entity_types;

    /// Canonical topology: entity → vertex.
    /// entity_to_vertex[d] stores the vertex indices for each entity of dimension d.
    /// This is CSR-style: offsets + flat indices.
    std::array<EntityAdjacency, max_dim> entity_to_vertex;

    /// Structured host provenance for entities. For zero faces this stores
    /// the uncut adaptcell leaf cell that produced the extracted surface
    /// element. For zero edges on a zero face boundary, host_face_id is filled
    /// by the face context when available.
    std::array<std::vector<std::int32_t>, max_dim> entity_host_cell_id;
    std::array<std::vector<cell::type>, max_dim> entity_host_cell_type;
    std::array<std::vector<std::int32_t>, max_dim> entity_host_face_id;
    std::array<std::vector<std::int32_t>, max_dim> entity_source_level_set;
    std::array<EntityAdjacency, max_dim> entity_host_cell_vertices;

    // ---------------------------------------------------------------
    // Connectivity cache (d0 → d1)
    //
    // Generic adjacency between entities of different dimensions.
    // Built on demand, not stored by default.
    // The canonical source is always entity_to_vertex; everything
    // else is derived from it.
    // ---------------------------------------------------------------

    std::array<std::array<EntityAdjacency, max_dim>, max_dim> connectivity;
    std::array<std::array<std::uint8_t, max_dim>, max_dim> has_connectivity = {};

    // ---------------------------------------------------------------
    // Zero-entity inventory (semantic layer for interfaces and intersections)
    //
    // A zero entity is an entity where one or more level sets vanish.
    // Each entry references an entity in the topology pools above.
    // This is the layer used for selections like "phi = 0" and
    // "phi1 = 0 and phi2 = 0".
    //
    // All arrays are parallel, indexed by zero-entity id z = 0, 1, ...
    // ---------------------------------------------------------------

    std::vector<std::uint8_t>  zero_entity_dim;         ///< dimension of the entity (1, 2, ...)
    std::vector<std::int32_t>  zero_entity_id;            ///< index into entity_types[dim] / entity_to_vertex[dim]
    std::vector<std::uint64_t> zero_entity_zero_mask;    ///< which level sets vanish: bit i set → φ_i = 0 on this entity
    std::vector<std::uint8_t>  zero_entity_is_owned;     ///< ownership flag for global assembly (avoids double-counting)

    /// Provenance of the zero entity back to the parent background cell.
    /// A zero entity born on parent face 2 gets parent_dim=2, parent_id=2.
    /// This enables global stitching of interface meshes.
    std::vector<std::int8_t>  zero_entity_parent_dim;    ///< -1 if not on a shared parent entity
    std::vector<std::int32_t> zero_entity_parent_id;     ///< -1 if not on a shared parent entity

    /// Provenance back to the uncut adaptcell leaf cell used to generate the
    /// zero entity. These arrays mirror entity_host_* for the zero-entity
    /// inventory and survive inventory rebuilds as long as the underlying
    /// entity provenance is available.
    std::vector<std::int32_t> zero_entity_host_cell_id;
    std::vector<cell::type> zero_entity_host_cell_type;
    std::vector<std::int32_t> zero_entity_host_face_id;
    std::vector<std::int32_t> zero_entity_source_level_set;
    EntityAdjacency zero_entity_host_cell_vertices;

    std::uint32_t zero_entity_version = 0;  ///< bumped when zero-entity inventory changes

    // ---------------------------------------------------------------
    // Per-level-set edge tags
    //
    // Flat storage: edge_root_tag[ls_id * n_edges(1) + edge_id]
    // ---------------------------------------------------------------

    std::vector<EdgeRootTag> edge_root_tag;
    int edge_root_tag_num_level_sets = 0;

    EdgeRootTag get_edge_root_tag(int level_set_id, int edge_id) const
    {
        return edge_root_tag[static_cast<std::size_t>(
            level_set_id * n_entities(1) + edge_id)];
    }

    void set_edge_root_tag(int level_set_id, int edge_id, EdgeRootTag tag)
    {
        edge_root_tag[static_cast<std::size_t>(
            level_set_id * n_entities(1) + edge_id)] = tag;
    }

    void resize_edge_root_tags(int num_level_sets)
    {
        const int old_num_level_sets = edge_root_tag_num_level_sets;
        const int old_n_edges = (old_num_level_sets > 0)
                                    ? static_cast<int>(edge_root_tag.size())
                                        / old_num_level_sets
                                    : 0;
        const int new_n_edges = n_entities(1);

        std::vector<EdgeRootTag> new_tags(
            static_cast<std::size_t>(num_level_sets * new_n_edges),
            EdgeRootTag::not_classified);

        const int copy_ls = std::min(old_num_level_sets, num_level_sets);
        const int copy_edges = std::min(old_n_edges, new_n_edges);
        for (int ls = 0; ls < copy_ls; ++ls)
        {
            for (int e = 0; e < copy_edges; ++e)
            {
                new_tags[static_cast<std::size_t>(ls * new_n_edges + e)] =
                    edge_root_tag[static_cast<std::size_t>(ls * old_n_edges + e)];
            }
        }

        edge_root_tag_num_level_sets = num_level_sets;
        edge_root_tag = std::move(new_tags);
    }

    // ---------------------------------------------------------------
    // Per-level-set cell certification tags
    //
    // Flat storage: cell_cert_tag[ls_id * n_entities(tdim) + cell_id]
    // ---------------------------------------------------------------

    std::vector<CellCertTag> cell_cert_tag;
    int cell_cert_tag_num_level_sets = 0;

    CellCertTag get_cell_cert_tag(int level_set_id, int cell_id) const
    {
        return cell_cert_tag[static_cast<std::size_t>(
            level_set_id * n_entities(tdim) + cell_id)];
    }

    void set_cell_cert_tag(int level_set_id, int cell_id, CellCertTag tag)
    {
        cell_cert_tag[static_cast<std::size_t>(
            level_set_id * n_entities(tdim) + cell_id)] = tag;
    }

    void resize_cell_cert_tags(int num_level_sets)
    {
        const int old_num_level_sets = cell_cert_tag_num_level_sets;
        const int old_n_cells = (old_num_level_sets > 0)
                                    ? static_cast<int>(cell_cert_tag.size())
                                        / old_num_level_sets
                                    : 0;
        const int new_n_cells = n_entities(tdim);

        std::vector<CellCertTag> new_tags(
            static_cast<std::size_t>(num_level_sets * new_n_cells),
            CellCertTag::not_classified);

        const int copy_ls = std::min(old_num_level_sets, num_level_sets);
        const int copy_cells = std::min(old_n_cells, new_n_cells);
        for (int ls = 0; ls < copy_ls; ++ls)
        {
            for (int c = 0; c < copy_cells; ++c)
            {
                new_tags[static_cast<std::size_t>(ls * new_n_cells + c)] =
                    cell_cert_tag[static_cast<std::size_t>(ls * old_n_cells + c)];
            }
        }

        cell_cert_tag_num_level_sets = num_level_sets;
        cell_cert_tag = std::move(new_tags);
    }

    // ---------------------------------------------------------------
    // Green-split metadata per edge per level set
    //
    // Flat storage: edge_green_split_param[ls_id * n_edges + edge_id]
    // edge_green_split_has_value: nonzero if the param is valid
    // ---------------------------------------------------------------

    std::vector<T> edge_green_split_param;
    std::vector<std::uint8_t> edge_green_split_has_value;

    void resize_green_split_data(int num_level_sets)
    {
        const int new_n_edges = n_entities(1);
        const int old_num_level_sets = (!edge_green_split_param.empty() && new_n_edges > 0)
                                           ? static_cast<int>(edge_green_split_param.size())
                                               / new_n_edges
                                           : 0;
        const int old_n_edges = new_n_edges;
        const auto n = static_cast<std::size_t>(num_level_sets * new_n_edges);

        std::vector<T> new_params(n, T(0));
        std::vector<std::uint8_t> new_has_value(n, std::uint8_t(0));

        const int copy_ls = std::min(old_num_level_sets, num_level_sets);
        const int copy_edges = std::min(old_n_edges, new_n_edges);
        for (int ls = 0; ls < copy_ls; ++ls)
        {
            for (int e = 0; e < copy_edges; ++e)
            {
                new_params[static_cast<std::size_t>(ls * new_n_edges + e)] =
                    edge_green_split_param[static_cast<std::size_t>(ls * old_n_edges + e)];
                new_has_value[static_cast<std::size_t>(ls * new_n_edges + e)] =
                    edge_green_split_has_value[static_cast<std::size_t>(ls * old_n_edges + e)];
            }
        }

        edge_green_split_param = std::move(new_params);
        edge_green_split_has_value = std::move(new_has_value);
    }

    // ---------------------------------------------------------------
    // One-root metadata per edge per level set
    //
    // Stores the localized root parameter on leaf edges tagged one_root
    // together with the AdaptCell vertex id created for that root.
    // ---------------------------------------------------------------

    std::vector<T> edge_one_root_param;
    std::vector<std::int32_t> edge_one_root_vertex_id;
    std::vector<std::uint8_t> edge_one_root_has_value;

    void resize_one_root_data(int num_level_sets)
    {
        const int new_n_edges = n_entities(1);
        const int old_num_level_sets = (!edge_one_root_param.empty() && new_n_edges > 0)
                                           ? static_cast<int>(edge_one_root_param.size())
                                               / new_n_edges
                                           : 0;
        const int old_n_edges = new_n_edges;
        const auto n = static_cast<std::size_t>(num_level_sets * new_n_edges);

        std::vector<T> new_params(n, T(0));
        std::vector<std::int32_t> new_vertex_ids(n, std::int32_t(-1));
        std::vector<std::uint8_t> new_has_value(n, std::uint8_t(0));

        const int copy_ls = std::min(old_num_level_sets, num_level_sets);
        const int copy_edges = std::min(old_n_edges, new_n_edges);
        for (int ls = 0; ls < copy_ls; ++ls)
        {
            for (int e = 0; e < copy_edges; ++e)
            {
                new_params[static_cast<std::size_t>(ls * new_n_edges + e)] =
                    edge_one_root_param[static_cast<std::size_t>(ls * old_n_edges + e)];
                new_vertex_ids[static_cast<std::size_t>(ls * new_n_edges + e)] =
                    edge_one_root_vertex_id[static_cast<std::size_t>(ls * old_n_edges + e)];
                new_has_value[static_cast<std::size_t>(ls * new_n_edges + e)] =
                    edge_one_root_has_value[static_cast<std::size_t>(ls * old_n_edges + e)];
            }
        }

        edge_one_root_param = std::move(new_params);
        edge_one_root_vertex_id = std::move(new_vertex_ids);
        edge_one_root_has_value = std::move(new_has_value);
    }

    // ---------------------------------------------------------------
    // Convenience accessors
    // ---------------------------------------------------------------

    /// Number of vertices in this cell
    int n_vertices() const
    {
        return tdim > 0 ? static_cast<int>(vertex_coords.size()) / tdim : 0;
    }

    /// Number of entities of dimension d
    int n_entities(int d) const
    {
        return (d >= 1 && d < max_dim)
                   ? static_cast<int>(entity_types[d].size())
                   : 0;
    }

    /// Check if connectivity map between dimensions d0 and d1 exists
    bool has_connectivity_map(int d0, int d1) const
    {
        return (d0 >= 1 && d0 < max_dim && d1 >= 1 && d1 < max_dim)
                   ? (has_connectivity[d0][d1] != 0)
                   : false;
    }

    /// Number of zero entities
    int n_zero_entities() const
    {
        return static_cast<int>(zero_entity_id.size());
    }

    // ---------------------------------------------------------------
    // Per-level-set face certification tags (3D only)
    //
    // Flat storage: face_cert_tag[ls_id * n_entities(2) + face_id]
    // ---------------------------------------------------------------

    std::vector<FaceCertTag> face_cert_tag;
    int face_cert_tag_num_level_sets = 0;

    FaceCertTag get_face_cert_tag(int level_set_id, int face_id) const
    {
        return face_cert_tag[static_cast<std::size_t>(
            level_set_id * n_entities(2) + face_id)];
    }

    void set_face_cert_tag(int level_set_id, int face_id, FaceCertTag tag)
    {
        face_cert_tag[static_cast<std::size_t>(
            level_set_id * n_entities(2) + face_id)] = tag;
    }

    void resize_face_cert_tags(int num_level_sets)
    {
        const int old_num_level_sets = face_cert_tag_num_level_sets;
        const int old_n_faces = (old_num_level_sets > 0)
                                    ? static_cast<int>(face_cert_tag.size())
                                        / old_num_level_sets
                                    : 0;
        const int new_n_faces = n_entities(2);

        std::vector<FaceCertTag> new_tags(
            static_cast<std::size_t>(num_level_sets * new_n_faces),
            FaceCertTag::not_classified);

        const int copy_ls = std::min(old_num_level_sets, num_level_sets);
        const int copy_faces = std::min(old_n_faces, new_n_faces);
        for (int ls = 0; ls < copy_ls; ++ls)
        {
            for (int f = 0; f < copy_faces; ++f)
            {
                new_tags[static_cast<std::size_t>(ls * new_n_faces + f)] =
                    face_cert_tag[static_cast<std::size_t>(ls * old_n_faces + f)];
            }
        }

        face_cert_tag_num_level_sets = num_level_sets;
        face_cert_tag = std::move(new_tags);
    }
};

/// Create an AdaptCell initialised from a single background cell in a MeshView.
///
/// The returned cell owns its own storage. Every vertex is tagged with
/// provenance dimension 0 (= parent vertex) and the global mesh node id.
/// Reference coordinates follow the canonical ordering defined in mapping.h.
/// Level-set sign masks are zeroed (all-positive) and must be filled
/// separately with fill_vertex_signs().
///
/// @param mesh     Background mesh. MeshView::cell_types must be present
///                 (cutcells integer codes).
/// @param cell_id  Index of the cell to extract.
/// @return         Initialised AdaptCell<T> in reference space.
template <std::floating_point T, std::integral I = int>
AdaptCell<T> make_adapt_cell(const MeshView<T, I>& mesh, I cell_id);

/// Fill the vertex-level-set sign masks of an AdaptCell for one level set.
///
/// @param ac                AdaptCell whose masks are updated in place.
/// @param vertex_ls_values  Level-set values at the vertices of the cell,
///                          in the same vertex ordering as ac.vertex_coords.
///                          Must have at least ac.n_vertices() entries.
/// @param ls_index          Which level set (bit position, 0–63).
/// @param tol               Absolute tolerance for the zero test.
template <std::floating_point T>
void fill_vertex_signs(AdaptCell<T>& ac,
                       std::span<const T> vertex_ls_values,
                       int ls_index,
                       T tol = T(1e-14));

/// Build (or rebuild) entity_to_vertex[1] from all top-dimensional leaf cells.
///
/// Clears any existing 1D entity pool, then re-derives edges by iterating
/// every cell in entity_to_vertex[tdim] and adding vertex-pairs from the
/// cell type's edge table (cell_topology.h ordering). Duplicate edges are
/// suppressed; each unique (min_v, max_v) pair appears once.
///
/// @param ac  AdaptCell to update in place.
template <std::floating_point T>
void build_edges(AdaptCell<T>& ac);

/// Build (or rebuild) entity_to_vertex[2] from all top-dimensional leaf cells.
///
/// Only meaningful for 3D cells (tdim == 3). Clears any existing 2D entity
/// pool, then re-derives faces by iterating every cell in entity_to_vertex[3]
/// and adding vertex-tuples from the cell type's face table
/// (cell_topology.h ordering). Duplicate faces are suppressed; each unique
/// sorted vertex set appears once.
///
/// @param ac  AdaptCell to update in place.
template <std::floating_point T>
void build_faces(AdaptCell<T>& ac);

/// Recompute per-leaf and per-cell active level-set masks from vertex masks.
///
/// For each top-dimensional leaf cell, bit i is set if the cell has a
/// sign change for level set i or contains at least one zero vertex for i.
///
/// @param ac               AdaptCell to update in place.
/// @param num_level_sets   Number of level sets to scan (bits 0..num_level_sets-1).
template <std::floating_point T>
void recompute_active_level_set_masks(AdaptCell<T>& ac, int num_level_sets);

/// Rebuild the zero-entity inventory from current topology and vertex masks.
///
/// Zero entities are registered for:
///   - dim 0: vertices with nonzero zero-mask
///   - dim 1..tdim-1: entities whose all vertices share at least one zero bit
///     (mask = intersection of vertex zero masks).
///
/// @param ac  AdaptCell to update in place.
template <std::floating_point T>
void rebuild_zero_entity_inventory(AdaptCell<T>& ac);

} // namespace cutcells
