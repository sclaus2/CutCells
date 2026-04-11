// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <span>
#include <vector>

#include "adapt_cell.h"

namespace cutcells
{

// =====================================================================
// Green refinement — split multiple-roots edges
// =====================================================================

/// Perform green refinement on cells incident to an edge tagged multiple_roots.
///
/// Each call refines the first multiple_roots edge with a stored green_split_t.
/// Repeated calls, or certify_and_refine(...), therefore resolve multiple such
/// edges recursively.
///
/// For the selected multiple_roots edge:
///   1. Insert a new vertex on the edge at parameter t.
///   2. Split all cells incident to that edge by connecting the new vertex.
///   3. Replace the old edge with two child edges.
///   4. New edges and new cells get tag not_classified.
///   5. Unaffected edges retain their tags.
///
/// @param adapt_cell     The AdaptCell (modified in place).
/// @param level_set_id   Which level set.
/// @return true if any refinement was performed.
template <std::floating_point T>
bool refine_green_on_multiple_root_edges(AdaptCell<T>& adapt_cell,
                                         int level_set_id);

// =====================================================================
// Red refinement — subdivide ambiguous cells
// =====================================================================

/// Perform red (full) refinement on cells tagged ambiguous.
///
/// Each ambiguous cell is subdivided uniformly: edge midpoints become
/// new vertices and the cell is replaced by 2^tdim children (or the
/// appropriate simplex subdivision).
///
/// @param adapt_cell     The AdaptCell (modified in place).
/// @param level_set_id   Which level set.
/// @return true if any refinement was performed.
template <std::floating_point T>
bool refine_red_on_ambiguous_cells(AdaptCell<T>& adapt_cell,
                                   int level_set_id);

/// Replace the current leaf-cell pool and rebuild the leaf-edge pool while
/// preserving certification state on surviving leaf entities. New leaf edges
/// and cells are left not_classified.
template <std::floating_point T>
void apply_topology_update_preserve_certification(
    AdaptCell<T>& adapt_cell,
    std::vector<cell::type>&& new_types,
    EntityAdjacency&& new_cells,
    std::span<const int> old_cell_ids_for_new_cells);

// =====================================================================
// Invalidation helpers
// =====================================================================

/// Reset edge tags to not_classified for a set of edge ids.
template <std::floating_point T>
void invalidate_edge_tags_for_new_edges(AdaptCell<T>& adapt_cell,
                                        std::span<const int> new_edge_ids);

/// Reset cell tags to not_classified for a set of cell ids.
template <std::floating_point T>
void invalidate_cell_tags_for_new_cells(AdaptCell<T>& adapt_cell,
                                        std::span<const int> new_cell_ids);

/// Reset face tags to not_classified for a set of face ids (3D cells).
template <std::floating_point T>
void invalidate_face_tags_for_new_faces(AdaptCell<T>& adapt_cell,
                                        std::span<const int> new_face_ids);

} // namespace cutcells
