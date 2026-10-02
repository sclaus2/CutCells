// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "../cell_flags.h"
#include "../cell_topology.h"
#include "../cell_types.h"
#include "../level_set.h"
#include "../mesh_view.h"

/// The front end: every cell classified by every level set, mesh parts
/// selected by expressions, and their quadrature and visualisation by backend.
namespace cutcells::part
{

/// The facets of a cell: its faces in 3D, its edges in 2D.
inline int num_facets(cell::type type)
{
    return cell::get_tdim(type) == 3 ? cell::num_faces(type) : cell::num_edges(type);
}

/// The vertices of facet @p f (Basix numbering of the cell).
inline std::span<const int> facet_vertices(cell::type type, int f)
{
    if (cell::get_tdim(type) == 3)
        return cell::face_vertices(type, f);
    return std::span<const int>(cell::edges(type)[static_cast<std::size_t>(f)]);
}

inline cell::type facet_type(cell::type type, int f)
{
    return cell::get_tdim(type) == 3 ? cell::face_type(type, f) : cell::type::interval;
}

/// Options of cut().
struct ClassifyOptions
{
    /// Bisections of a cell's box before its sign is given up as unknown; the
    /// cell then counts as cut. Bernstein bounds rarely need more than a few,
    /// Taylor models of distance functions more.
    int max_depth = 12;
};

/// Every cell classified by every level set as inside (phi < 0), outside
/// (phi > 0) or cut, and the faces lying in a level set's zero set.
template <std::floating_point T, std::integral I = int>
struct CutResult
{
    const MeshView<T, I>* mesh = nullptr;                  ///< not owned
    std::vector<const LevelSetFunction<T, I>*> level_sets; ///< not owned
    std::vector<std::string> level_set_names;
    int num_cells = 0;
    /// domains[l * num_cells + c]: inside, outside or intersected
    std::vector<cell::domain> domains;
    /// Cells that some level set cuts, ascending.
    std::vector<I> cut_cells;
    /// Facets (faces in 3D, edges in 2D) in the zero set of a level set, each
    /// once: level set, owning cell, and the facet in that cell (Basix
    /// numbering). The owner is the cell on the negative side, or the lower
    /// cell index if neither or both sides are negative.
    std::vector<int> zero_face_level_sets;
    std::vector<I> zero_face_cells;
    std::vector<std::int8_t> zero_face_local;

    int n_level_sets() const { return static_cast<int>(level_sets.size()); }
    int n_zero_faces() const { return static_cast<int>(zero_face_cells.size()); }
    cell::domain domain(int level_set, I cell_id) const
    {
        return domains[static_cast<std::size_t>(level_set) * static_cast<std::size_t>(num_cells)
                       + static_cast<std::size_t>(cell_id)];
    }
};

/// @brief Bit l set if level set l cuts @p cell_id.
template <std::floating_point T, std::integral I>
std::uint64_t cut_mask(const CutResult<T, I>& result, I cell_id)
{
    std::uint64_t mask = 0;
    for (int l = 0; l < result.n_level_sets(); ++l)
        if (result.domain(l, cell_id) == cell::domain::intersected)
            mask |= std::uint64_t(1) << l;
    return mask;
}

/// @brief Classify every cell of @p mesh by every level set, and find the faces
/// lying in their zero sets.
///
/// Tetrahedra and hexahedra are classified by bounds of the level set itself:
/// Bernstein coefficients of a Pk level set, or Taylor models of an analytic
/// one, over the cell and its halves (cell_sign in classify.h). A cell whose
/// sign the bounds cannot prove counts as cut, which the backends handle.
/// Other cells use the signs of their Bernstein coefficients.
///
/// @param level_sets  Pk level sets (dof values) or analytic ones; they and
///                    @p mesh must outlive the result
/// @throws std::invalid_argument for level sets with neither, more than 64
///         level sets, or analytic level sets on cells other than
///         tetrahedra and hexahedra
template <std::floating_point T, std::integral I>
CutResult<T, I> cut(const MeshView<T, I>& mesh, std::span<const LevelSetFunction<T, I>> level_sets,
                    const ClassifyOptions& options = {});

} // namespace cutcells::part
