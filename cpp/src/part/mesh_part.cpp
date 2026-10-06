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

/// Whether @p term has a piece on entity @p e of @p table that is the whole
/// entity: its zero clauses all count there, none of its strict clauses is on
/// a level set vanishing there, and a sharing cell holds each of its other
/// clauses whole (no sharing cell holding the opposite).
template <std::floating_point T, std::integral I>
bool term_on_entity(const SelectionTerm& term, const CutResult<T, I>& result, const EntityTable<I>& table,
                    std::int32_t e)
{
    const std::uint64_t strict = term.negative_required | term.positive_required;
    if (term.zero_required == 0 || (term.zero_required & ~table.zero[static_cast<std::size_t>(e)]) != 0
        || (strict & table.vanishing[static_cast<std::size_t>(e)]) != 0
        || (term.negative_required & term.positive_required) != 0)
        return false;
    for (std::uint64_t bits = strict; bits != 0; bits &= bits - 1)
    {
        const int l = std::countr_zero(bits);
        const bool negative = (term.negative_required >> l) & 1;
        const cell::domain side = negative ? cell::domain::inside : cell::domain::outside;
        const cell::domain opposite = negative ? cell::domain::outside : cell::domain::inside;
        bool holds = false;
        for (const I c : table.entity_cells(e))
        {
            const cell::domain d = result.domain(l, c);
            if (d == opposite)
                return false;
            holds |= d == side;
        }
        if (!holds)
            return false; // every sharing cell is cut by l: a part of the entity
    }
    return true;
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

    // the entities of the part's dimension that it holds whole
    if (part.dim >= 0 && part.dim < tdim)
    {
        const std::size_t d = static_cast<std::size_t>(part.dim);
        const EntityTable<I>& table = result.entities[d];
        for (std::int32_t e = 0; e < table.n_entities(); ++e)
        {
            if (table.zero[static_cast<std::size_t>(e)] == 0)
                continue;
            std::uint64_t zero = 0;
            for (const SelectionTerm& term : part.expr.terms)
                if (term_on_entity(term, result, table, e))
                    zero |= term.zero_required;
            if (zero == 0)
                continue;
            part.entities[d].push_back(e);
            part.entity_zero[d].push_back(zero);
        }
    }
    // the zero faces among them
    if (part.dim == tdim - 1)
    {
        const EntityTable<I>& facets = result.entities[static_cast<std::size_t>(part.dim)];
        const std::vector<std::int32_t>& held = part.entities[static_cast<std::size_t>(part.dim)];
        std::vector<I> nodes, scratch, local;
        for (int z = 0; z < result.n_zero_faces(); ++z)
        {
            const I owner = result.zero_face_cells[static_cast<std::size_t>(z)];
            const cell::type type = result.mesh->cell_type(owner);
            cell_vertex_nodes(*result.mesh, owner, nodes, scratch);
            local.clear();
            for (const int v : facet_vertices(type, result.zero_face_local[static_cast<std::size_t>(z)]))
                local.push_back(nodes[static_cast<std::size_t>(v)]);
            const std::int32_t e = find_entity(facets, entity_key(std::span<const I>(local)));
            const auto it = std::lower_bound(held.begin(), held.end(), e);
            if (it == held.end() || *it != e)
                continue;
            const std::uint64_t zero
                = part.entity_zero[static_cast<std::size_t>(part.dim)][static_cast<std::size_t>(it - held.begin())];
            if ((zero >> result.zero_face_level_sets[static_cast<std::size_t>(z)]) & 1)
                part.zero_faces.push_back(z);
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
