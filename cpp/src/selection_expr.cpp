// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "selection_expr.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace cutcells
{
namespace
{

/// Trim leading and trailing whitespace from a string_view.
std::string_view trim(std::string_view s)
{
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front())))
        s.remove_prefix(1);
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back())))
        s.remove_suffix(1);
    return s;
}

/// Parse a single clause: "name < 0", "name > 0", or "name = 0".
Clause parse_clause(std::string_view text)
{
    text = trim(text);
    if (text.empty())
        throw std::runtime_error("parse_selection_expr: empty clause");

    // Find the relation operator: '<', '>', or '='
    auto pos_lt = text.find('<');
    auto pos_gt = text.find('>');
    auto pos_eq = text.find('=');

    std::size_t op_pos = std::string_view::npos;
    Relation rel{};

    if (pos_lt != std::string_view::npos
        && (pos_lt < pos_gt || pos_gt == std::string_view::npos)
        && (pos_lt < pos_eq || pos_eq == std::string_view::npos))
    {
        op_pos = pos_lt;
        rel = Relation::LessThan;
    }
    else if (pos_gt != std::string_view::npos
             && (pos_gt < pos_eq || pos_eq == std::string_view::npos))
    {
        op_pos = pos_gt;
        rel = Relation::GreaterThan;
    }
    else if (pos_eq != std::string_view::npos)
    {
        op_pos = pos_eq;
        rel = Relation::EqualTo;
    }
    else
    {
        throw std::runtime_error(
            "parse_selection_expr: no operator (<, >, =) in clause: '"
            + std::string(text) + "'");
    }

    std::string_view name_part = trim(text.substr(0, op_pos));
    std::string_view rhs_part  = trim(text.substr(op_pos + 1));

    if (name_part.empty())
        throw std::runtime_error(
            "parse_selection_expr: missing level-set name in clause: '"
            + std::string(text) + "'");

    if (rhs_part != "0")
        throw std::runtime_error(
            "parse_selection_expr: right-hand side must be '0', got '"
            + std::string(rhs_part) + "' in clause: '" + std::string(text) + "'");

    Clause c;
    c.name = std::string(name_part);
    c.relation = rel;
    return c;
}

/// Case-insensitive check for " and " delimiter between clauses.
/// Returns the position of the 'a' in " and ", or npos.
std::size_t find_and_delimiter(std::string_view text, std::size_t start = 0)
{
    // Search for " and " (space-delimited, case-insensitive)
    for (std::size_t i = start; i + 4 < text.size(); ++i)
    {
        if (std::isspace(static_cast<unsigned char>(text[i]))
            && (text[i + 1] == 'a' || text[i + 1] == 'A')
            && (text[i + 2] == 'n' || text[i + 2] == 'N')
            && (text[i + 3] == 'd' || text[i + 3] == 'D')
            && std::isspace(static_cast<unsigned char>(text[i + 4])))
        {
            return i;
        }
    }
    return std::string_view::npos;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// parse_selection_expr
// ---------------------------------------------------------------------------

SelectionExpr parse_selection_expr(std::string_view text)
{
    text = trim(text);
    if (text.empty())
        throw std::runtime_error("parse_selection_expr: empty expression");

    SelectionExpr expr;

    std::size_t pos = 0;
    while (pos < text.size())
    {
        std::size_t and_pos = find_and_delimiter(text, pos);
        std::string_view clause_text;
        if (and_pos == std::string_view::npos)
        {
            clause_text = text.substr(pos);
            pos = text.size();
        }
        else
        {
            clause_text = text.substr(pos, and_pos - pos);
            pos = and_pos + 5; // skip " and "
        }

        expr.clauses.push_back(parse_clause(clause_text));
    }

    if (expr.clauses.empty())
        throw std::runtime_error("parse_selection_expr: no clauses found");

    return expr;
}

// ---------------------------------------------------------------------------
// compile_selection_expr
// ---------------------------------------------------------------------------

void compile_selection_expr(SelectionExpr& expr,
                            const std::vector<std::string>& level_set_names)
{
    expr.zero_required     = 0;
    expr.negative_required = 0;
    expr.positive_required = 0;

    for (auto& clause : expr.clauses)
    {
        auto it = std::find(level_set_names.begin(), level_set_names.end(),
                            clause.name);
        if (it == level_set_names.end())
        {
            throw std::runtime_error(
                "compile_selection_expr: unknown level-set name '"
                + clause.name + "'");
        }

        int idx = static_cast<int>(std::distance(level_set_names.begin(), it));
        if (idx >= 64)
            throw std::runtime_error(
                "compile_selection_expr: level-set index >= 64 not supported");

        clause.level_set_index = idx;
        const std::uint64_t bit = std::uint64_t(1) << idx;

        switch (clause.relation)
        {
            case Relation::LessThan:    expr.negative_required |= bit; break;
            case Relation::GreaterThan: expr.positive_required |= bit; break;
            case Relation::EqualTo:     expr.zero_required     |= bit; break;
        }
    }
}

// ---------------------------------------------------------------------------
// infer_selection_dim
// ---------------------------------------------------------------------------

int infer_selection_dim(const SelectionExpr& expr, int tdim)
{
    int n_eq = 0;
    for (const auto& clause : expr.clauses)
        if (clause.relation == Relation::EqualTo)
            ++n_eq;
    return tdim - n_eq;
}

} // namespace cutcells
