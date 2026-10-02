// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "mesh_part.h"

#include "../cell_types.h"

namespace cutcells::part
{

template <std::floating_point T, std::integral I>
MeshPart<T, I> select(const CutResult<T, I>& result, std::string_view expr)
{
    MeshPart<T, I> part;
    part.result = &result;
    part.expr = parse_selection_expr(expr);
    compile_selection_expr(part.expr, result.level_set_names);
    const int tdim = result.num_cells > 0 ? cell::get_tdim(result.mesh->cell_type(I(0))) : 3;
    part.dim = infer_selection_dim(part.expr, tdim);

    for (I c = 0; c < static_cast<I>(result.num_cells); ++c)
    {
        bool whole = false, piece = false;
        for (const SelectionTerm& term : part.expr.terms)
        {
            const TermCell tc = term_on_cell(term, result, c);
            whole |= tc == TermCell::whole;
            piece |= tc == TermCell::piece;
        }
        if (whole && part.dim == tdim)
            part.uncut_cells.push_back(c);
        else if (piece)
            part.cut_cells.push_back(c);
    }

    if (part.dim == tdim - 1)
    {
        for (int z = 0; z < result.n_zero_faces(); ++z)
        {
            const std::uint64_t bit = std::uint64_t(1) << result.zero_face_level_sets[static_cast<std::size_t>(z)];
            const I owner = result.zero_face_cells[static_cast<std::size_t>(z)];
            for (const SelectionTerm& term : part.expr.terms)
            {
                if (!(term.zero_required & bit))
                    continue;
                // the term without its zero clause must hold on the whole owner
                SelectionTerm rest = term;
                rest.zero_required &= ~bit;
                if (term_on_cell(rest, result, owner) == TermCell::whole)
                {
                    part.zero_faces.push_back(z);
                    break;
                }
            }
        }
    }
    return part;
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template MeshPart<float, int> select<float, int>(const CutResult<float, int>&, std::string_view);
template MeshPart<double, int> select<double, int>(const CutResult<double, int>&, std::string_view);

} // namespace cutcells::part
