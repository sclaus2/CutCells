// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <string_view>
#include <vector>

#include "../selection_expr.h"
#include "cut_result.h"

namespace cutcells::part
{

/// A part of the mesh selected by an expression such as "phi1 < 0 and
/// phi2 = 0": the cells wholly in it, the cut cells holding a piece of it,
/// which a backend integrates, and the zero faces in it.
template <std::floating_point T, std::integral I = int>
struct MeshPart
{
    const CutResult<T, I>* result = nullptr; ///< not owned
    SelectionExpr expr;
    int dim = -1;                ///< tdim for volume parts, tdim - 1 for interfaces
    std::vector<I> uncut_cells;  ///< cells wholly in a volume part, ascending
    std::vector<I> cut_cells;    ///< cells holding a piece the backend integrates, ascending
    std::vector<int> zero_faces; ///< zero faces of the result in an interface part

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
/// a piece if no term holds on all of it and some term holds on part of it. An
/// interface part gets the cut cells of its level set and the zero faces of
/// that level set whose owning cells satisfy the term's other clauses.
/// @throws std::runtime_error on syntax errors or unknown level-set names
template <std::floating_point T, std::integral I>
MeshPart<T, I> select(const CutResult<T, I>& result, std::string_view expr);

} // namespace cutcells::part
