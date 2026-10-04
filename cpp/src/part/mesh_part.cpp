// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "mesh_part.h"

#include <algorithm>
#include <bit>
#include <stdexcept>

#include "../cell_types.h"

namespace cutcells::part
{

namespace
{

template <std::floating_point T, std::integral I>
int result_tdim(const CutResult<T, I>& result)
{
    return result.num_cells > 0 ? cell::get_tdim(result.mesh->cell_type(I(0))) : 3;
}

/// The part of @p result that the compiled @p expr selects, of dimension @p dim.
template <std::floating_point T, std::integral I>
MeshPart<T, I> select_compiled(const CutResult<T, I>& result, SelectionExpr expr, int dim)
{
    MeshPart<T, I> part;
    part.result = &result;
    part.expr = std::move(expr);
    part.dim = dim;
    const int tdim = result_tdim(result);

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

/// The level sets of @p result that @p name names exactly: a level set, or
/// all of those a level set made of level sets is made of; 0 if none.
template <std::floating_point T, std::integral I>
std::uint64_t named_level_sets(const CutResult<T, I>& result, const std::string& name)
{
    const std::vector<std::string>& names = result.level_set_names;
    if (const auto it = std::find(names.begin(), names.end(), name); it != names.end())
        return std::uint64_t(1) << (it - names.begin());
    if (result.insides.count(name) == 0)
        return 0;
    std::uint64_t mask = 0;
    const std::string prefix = name + ".";
    for (std::size_t l = 0; l < names.size(); ++l)
        if (names[l].starts_with(prefix))
            mask |= std::uint64_t(1) << l;
    return mask;
}

std::string joined(const std::vector<std::string>& names)
{
    std::string s;
    for (const std::string& n : names)
        s += (s.empty() ? "" : ", ") + n;
    return s;
}

} // namespace

template <std::floating_point T, std::integral I>
MeshPart<T, I> select(const CutResult<T, I>& result, std::string_view expr)
{
    const SelectionExpr written = parse_selection_expr(expr);
    SelectionExpr compiled = written;
    compile_selection_expr(compiled, result.level_set_names, result.insides);
    // the dimension of the selection as written: "body = 0" expands to faces
    return select_compiled(result, std::move(compiled), infer_selection_dim(written, result_tdim(result)));
}

template <std::floating_point T, std::integral I>
MeshPart<T, I> on(const MeshPart<T, I>& part, const std::vector<std::string>& names)
{
    const CutResult<T, I>& r = *part.result;
    if (part.dim >= result_tdim(r))
        throw std::invalid_argument("on: a volume part has no faces; select a surface first, such as 'body = 0'");
    if (names.empty())
        throw std::invalid_argument("on: name one or more faces");
    std::uint64_t faces = 0;
    for (const SelectionTerm& term : part.expr.terms)
        faces |= term.zero_required;
    std::vector<std::string> face_names;
    for (std::uint64_t bits = faces; bits != 0; bits &= bits - 1)
        face_names.push_back(r.level_set_names[static_cast<std::size_t>(std::countr_zero(bits))]);

    std::uint64_t keep = 0;
    for (const std::string& name : names)
    {
        std::uint64_t mask = named_level_sets(r, name);
        if (mask == 0)
        {
            // a short name: "X.name" for X made of level sets, among the faces if one is
            std::vector<std::string> found, among_faces;
            auto consider = [&](const std::string& q)
            {
                if (q.size() > name.size() + 1 && q.ends_with(name) && q[q.size() - name.size() - 1] == '.'
                    && r.insides.count(q.substr(0, q.size() - name.size() - 1)) > 0)
                {
                    found.push_back(q);
                    if (named_level_sets(r, q) & faces)
                        among_faces.push_back(q);
                }
            };
            for (const std::string& q : r.level_set_names)
                consider(q);
            for (const auto& entry : r.insides)
                consider(entry.first);
            const std::vector<std::string>& candidates = among_faces.empty() ? found : among_faces;
            if (candidates.empty())
                throw std::invalid_argument("on: no level set is named '" + name + "'");
            if (candidates.size() > 1)
                throw std::invalid_argument("on: the short name '" + name + "' could be any of " + joined(candidates)
                                            + "; give the full name");
            mask = named_level_sets(r, candidates.front());
        }
        if ((mask & faces) == 0)
            throw std::invalid_argument("on: '" + name + "' is none of the part's faces (" + joined(face_names) + ")");
        keep |= mask;
    }

    SelectionExpr expr;
    for (const SelectionTerm& term : part.expr.terms)
        if (term.zero_required & keep)
            expr.terms.push_back(term);
    // the legacy single-term view
    expr.clauses = expr.terms.front().clauses;
    expr.zero_required = expr.terms.front().zero_required;
    expr.negative_required = expr.terms.front().negative_required;
    expr.positive_required = expr.terms.front().positive_required;
    return select_compiled(r, std::move(expr), part.dim);
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template MeshPart<float, int> select<float, int>(const CutResult<float, int>&, std::string_view);
template MeshPart<double, int> select<double, int>(const CutResult<double, int>&, std::string_view);
template MeshPart<float, int> on<float, int>(const MeshPart<float, int>&, const std::vector<std::string>&);
template MeshPart<double, int> on<double, int>(const MeshPart<double, int>&, const std::vector<std::string>&);

} // namespace cutcells::part
