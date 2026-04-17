// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "adapt_cell.h"
#include "cell_flags.h"
#include "level_set.h"
#include "level_set_cell.h"
#include "mesh_view.h"
#include "selection_expr.h"

namespace cutcells
{

// =====================================================================
// BackgroundMeshData — mesh-wide metadata produced during cutting
// =====================================================================

/// Mesh-wide metadata produced during cutting.
///
/// Holds a non-owning reference to the background mesh plus:
///   - level-set name registry
///   - per-cell, per-level-set domain classification
///   - cell_to_cut_index lookup
template <std::floating_point T, std::integral I = int>
struct BackgroundMeshData
{
    /// Non-owning reference to the background mesh.
    const MeshView<T, I>* mesh = nullptr;

    /// Level-set name registry.
    /// Index in this vector = bit position in AdaptCell bitmasks.
    std::vector<std::string> level_set_names;

    /// Per-cell, per-level-set domain classification.
    /// Flat storage: cell_domains[ls_index * num_cells + cell_id].
    std::vector<cell::domain> cell_domains;
    int num_cells      = 0;  ///< = mesh->num_cells()
    int num_level_sets = 0;

    /// Lookup: background cell id → index into HOCutCells arrays.
    /// -1 for uncut cells.  Size = num_cells.
    std::vector<int> cell_to_cut_index;

    /// Access the domain of a specific (level_set, cell) pair.
    cell::domain domain(int ls_index, I cell_id) const
    {
        return cell_domains[static_cast<std::size_t>(
            ls_index * num_cells + static_cast<int>(cell_id))];
    }
};

// =====================================================================
// HOCutCells — pure cut-cell storage (intersected cells only)
// =====================================================================

/// Pure storage of intersected cells.
/// No mesh reference, no uncut cells, no classification.
template <std::floating_point T, std::integral I = int>
struct HOCutCells
{
    int gdim = 0;
    int tdim = 0;

    /// Background-mesh cell ids for the cut cells.
    /// parent_cell_ids[k] indexes the background mesh for the k-th entry.
    std::vector<I> parent_cell_ids;

    /// AdaptCell (reference-space topology + vertex signs) for each cut cell.
    std::vector<AdaptCell<T>> adapt_cells;

    /// Per cut-cell level-set data, flat CSR layout:
    /// level_set_cells[ls_offsets[k] .. ls_offsets[k+1]] for k-th cut cell.
    std::vector<LevelSetCell<T, I>> level_set_cells;
    std::vector<int> ls_offsets;  ///< size = num_cut_cells + 1

    /// Bitmask of active (intersecting) level sets for each cut cell.
    /// Bit li set ↔ level set li has domain::intersected for the k-th cut cell.
    /// Size = num_cut_cells.
    std::vector<std::uint64_t> active_level_set_mask;

    /// Number of intersected cells.
    int num_cut_cells() const
    {
        return static_cast<int>(parent_cell_ids.size());
    }
};

// =====================================================================
// HOMeshPart — pure view/extraction layer over HOCutCells + BackgroundMeshData
// =====================================================================

/// View over HOCutCells + BackgroundMeshData, filtered by a SelectionExpr.
///
/// Stores no data of its own beyond the selection results.
template <std::floating_point T, std::integral I = int>
struct HOMeshPart
{
    /// Non-owning references.
    const HOCutCells<T, I>* cut_cells = nullptr;
    const BackgroundMeshData<T, I>* bg = nullptr;

    /// The selection expression (compiled).
    SelectionExpr expr;

    /// Inferred entity dimension of the selection.
    int dim = -1;

    /// Indices into cut_cells->adapt_cells for cells with matching sub-entities.
    std::vector<std::int32_t> cut_cell_ids;

    /// Background cell ids of uncut cells matching the selection.
    std::vector<I> uncut_cell_ids;

    /// true: iterate only cut cells (quadrature).
    /// false: iterate both cut + uncut (visualization).
    bool cut_only = false;
};

// =====================================================================
// Factory: cut() — produces both HOCutCells + BackgroundMeshData
// =====================================================================

/// Build HOCutCells and BackgroundMeshData from a mesh and a single level set.
///
/// Iterates all background cells, classifies each by a cheap Bernstein-cell
/// sign test when available (falling back to vertex signs otherwise), and
/// for intersected cells creates LevelSetCell + AdaptCell pairs.
///
/// @param mesh  Background mesh. Must have cell types.
/// @param ls    Level-set function.
/// @return pair of (HOCutCells, BackgroundMeshData).
template <std::floating_point T, std::integral I = int>
std::pair<HOCutCells<T, I>, BackgroundMeshData<T, I>>
cut(const MeshView<T, I>& mesh,
    const LevelSetFunction<T, I>& ls);

/// Build HOCutCells and BackgroundMeshData from a mesh and multiple level sets.
///
/// @param mesh        Background mesh.
/// @param level_sets  Vector of level-set functions.
/// @return pair of (HOCutCells, BackgroundMeshData).
template <std::floating_point T, std::integral I = int>
std::pair<HOCutCells<T, I>, BackgroundMeshData<T, I>>
cut(const MeshView<T, I>& mesh,
    const std::vector<LevelSetFunction<T, I>>& level_sets);

// =====================================================================
// select_part() — builds HOMeshPart
// =====================================================================

/// Build an HOMeshPart by filtering HOCutCells + BackgroundMeshData
/// with a selection expression string.
///
/// @param cut_cells  Intersected cell storage.
/// @param bg         Mesh-wide metadata.
/// @param expr_str   Selection expression, e.g. "phi < 0", "phi1 = 0 and phi2 < 0".
/// @return populated HOMeshPart.
template <std::floating_point T, std::integral I = int>
HOMeshPart<T, I> select_part(const HOCutCells<T, I>& cut_cells,
                              const BackgroundMeshData<T, I>& bg,
                              std::string_view expr_str);

} // namespace cutcells
