// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <concepts>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "../selection_expr.h"
#include "cut_result.h"

namespace cutcells::part
{

/// A part of the mesh selected by an expression such as "phi1 < 0 and
/// phi2 = 0": the cells wholly in it, the cut cells holding a piece of it,
/// which a backend integrates, and the facets, edges and vertices wholly in
/// it, whose rules are made once, on the entity (entity_rules).
template <std::floating_point T, std::integral I = int>
struct MeshPart
{
    const CutResult<T, I>* result = nullptr; ///< not owned
    SelectionExpr expr;
    int dim = -1;                ///< tdim for volume parts, tdim - 1 for interfaces
    std::vector<I> uncut_cells;  ///< cells wholly in a volume part, ascending
    std::vector<I> cut_cells;    ///< cells holding a piece the backend integrates, ascending
    /// Curve parts in cells of dimension 3: the cut cells that take a curve lying in
    /// one of their facets, where one of a term's two level sets vanishes on the
    /// facet and the other cuts the cell. The first cell sharing the facet takes it,
    /// so it counts once, with that level set moved off the facet into the cell
    /// (quadrature_rules). Ascending, and parallel: the level sets moved, by the side
    /// of them the cell lies on.
    std::vector<I> nudged_cells;
    std::vector<std::uint64_t> nudged_negative, nudged_positive;
    /// The entities wholly in a part of lower dimension, by dimension: indices
    /// into result->entities[d], ascending. A piece of the part is a whole
    /// entity where the zero sets it lies in all vanish there (and change sign
    /// across it), and the term's other clauses hold on all of it.
    std::array<std::vector<std::int32_t>, 3> entities;
    /// For each of them, the zero level sets of the terms that select it.
    std::array<std::vector<std::uint64_t>, 3> entity_zero;
    /// The entities of the part's dimension it holds a part of, by dimension:
    /// a term's zero clauses all count there, but a level set of its other
    /// clauses cuts every sharing cell, so the piece is bounded inside the
    /// entity (TermEntity::part). Indices into result->entities[d], ascending.
    std::array<std::vector<std::int32_t>, 3> partial_entities;
    /// For each of them, the zero level sets of the terms that select a part of it.
    std::array<std::vector<std::uint64_t>, 3> partial_entity_zero;
    /// The result's zero faces whose facet and level set the part holds (the
    /// facets among entities, with one owner each; CutResult::zero_face_cells).
    std::vector<int> zero_faces;

    int n_uncut_cells() const { return static_cast<int>(uncut_cells.size()); }
    int n_cut_cells() const { return static_cast<int>(cut_cells.size()); }
};

/// How a selection term meets a cell.
enum class TermCell : std::uint8_t
{
    none,  ///< the term selects nothing in the cell
    whole, ///< the whole cell
    piece  ///< a piece bounded by the level sets that cut the cell
};

/// How a selection term meets a facet, an edge or a vertex of the mesh.
enum class TermEntity : std::uint8_t
{
    none,  ///< the term selects nothing on it
    whole, ///< the whole entity
    part   ///< a part of it, bounded by level sets that cut every sharing cell
};

/// @brief How @p term meets entity @p e of @p table: its zero clauses must all
/// count there (they vanish there) and none of its strict clauses be on a
/// level set vanishing there; a strict clause holds on all of it where a
/// sharing cell holds it whole, on none where one holds the opposite, and on a
/// part where every sharing cell is cut by its level set.
/// @param cutting  output: for TermEntity::part, the level sets of the strict
///                 clauses that cut every sharing cell
template <std::floating_point T, std::integral I>
TermEntity term_on_entity(const SelectionTerm& term, const CutResult<T, I>& result, const EntityTable<I>& table,
                          std::int32_t e, std::uint64_t& cutting);

/// @brief How @p term meets cell @p cell_id, from the cell's domains.
template <std::floating_point T, std::integral I>
TermCell term_on_cell(const SelectionTerm& term, const CutResult<T, I>& result, I cell_id)
{
    bool piece = false;
    const std::uint64_t all = term.negative_required | term.positive_required | term.zero_required;
    for (int l = 0; l < result.n_level_sets(); ++l)
    {
        const std::uint64_t bit = std::uint64_t(1) << l;
        if (!(all & bit))
            continue;
        const cell::domain d = result.domain(l, cell_id);
        if ((term.negative_required & bit) && (term.positive_required & bit))
            return TermCell::none;
        if (term.zero_required & bit)
        {
            if (d != cell::domain::intersected)
                return TermCell::none;
            piece = true;
        }
        else if (term.negative_required & bit)
        {
            if (d == cell::domain::outside)
                return TermCell::none;
            piece |= d == cell::domain::intersected;
        }
        else
        {
            if (d == cell::domain::inside)
                return TermCell::none;
            piece |= d == cell::domain::intersected;
        }
    }
    return piece ? TermCell::piece : TermCell::whole;
}

/// @brief The part of @p result that @p expr selects.
///
/// A cell is wholly in a volume part if some term holds on all of it; it holds
/// a piece if no term holds on all of it and some term holds on part of it.
/// A part of lower dimension also holds the entities of its dimension on which
/// a term's zero clauses all count (they vanish there and change sign across
/// it, or it lies on the mesh boundary), none of its strict clauses is on a
/// level set vanishing there, and its other clauses hold on a sharing cell
/// whole (entities). Where every sharing cell is cut by the level set of such
/// a clause, the part holds the piece of the entity where it holds
/// (partial_entities).
///
/// Level sets made of level sets expand into terms on theirs
/// (compile_selection_expr with the result's insides). The dimension is that
/// of @p expr as written: "body = 0" selects an interface, also where it
/// expands to no term.
/// @throws std::runtime_error on syntax errors or unknown level-set names
template <std::floating_point T, std::integral I>
MeshPart<T, I> select(const CutResult<T, I>& result, std::string_view expr);

/// @brief The faces of @p part on the level sets named @p names: the part's
/// terms whose zero sets include one of them (no point is filtered).
///
/// A name is that of a level set of the result ("body.hole"); that of a
/// level set made of level sets, for all its faces ("body"); or a short name
/// ("hole"), the level set "X.hole" of a level set X made of level sets,
/// among the part's faces if one is.
/// @throws std::invalid_argument for parts of full dimension, no names,
///         unknown or ambiguous names, and names that are no face of the part
template <std::floating_point T, std::integral I>
MeshPart<T, I> on(const MeshPart<T, I>& part, const std::vector<std::string>& names);

} // namespace cutcells::part
