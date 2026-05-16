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

/// Case-insensitive check for a space-delimited keyword between clauses/terms.
/// Returns the position of the leading whitespace before the keyword, or npos.
std::size_t find_keyword_delimiter(std::string_view text,
                                   std::string_view keyword,
                                   std::size_t start = 0)
{
    for (std::size_t i = start; i + keyword.size() + 1 < text.size(); ++i)
    {
        if (!std::isspace(static_cast<unsigned char>(text[i])))
            continue;

        bool match = true;
        for (std::size_t k = 0; k < keyword.size(); ++k)
        {
            const char a = static_cast<char>(
                std::tolower(static_cast<unsigned char>(text[i + 1 + k])));
            const char b = static_cast<char>(
                std::tolower(static_cast<unsigned char>(keyword[k])));
            if (a != b)
            {
                match = false;
                break;
            }
        }

        if (match
            && std::isspace(static_cast<unsigned char>(
                text[i + 1 + keyword.size()])))
        {
            return i;
        }
    }
    return std::string_view::npos;
}

SelectionTerm parse_term(std::string_view text)
{
    text = trim(text);
    if (text.empty())
        throw std::runtime_error("parse_selection_expr: empty term");

    SelectionTerm term;
    std::size_t pos = 0;
    while (pos < text.size())
    {
        const std::size_t and_pos = find_keyword_delimiter(text, "and", pos);
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

        term.clauses.push_back(parse_clause(clause_text));
    }

    if (term.clauses.empty())
        throw std::runtime_error("parse_selection_expr: no clauses found in term");
    return term;
}

void mirror_first_term(SelectionExpr& expr)
{
    if (expr.terms.empty())
    {
        expr.clauses.clear();
        expr.zero_required = 0;
        expr.negative_required = 0;
        expr.positive_required = 0;
        return;
    }

    expr.clauses = expr.terms.front().clauses;
    expr.zero_required = expr.terms.front().zero_required;
    expr.negative_required = expr.terms.front().negative_required;
    expr.positive_required = expr.terms.front().positive_required;
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
        const std::size_t or_pos = find_keyword_delimiter(text, "or", pos);
        std::string_view term_text;
        if (or_pos == std::string_view::npos)
        {
            term_text = text.substr(pos);
            pos = text.size();
        }
        else
        {
            term_text = text.substr(pos, or_pos - pos);
            pos = or_pos + 4; // skip " or "
        }

        expr.terms.push_back(parse_term(term_text));
    }

    if (expr.terms.empty())
        throw std::runtime_error("parse_selection_expr: no terms found");

    mirror_first_term(expr);

    return expr;
}

// ---------------------------------------------------------------------------
// compile_selection_expr
// ---------------------------------------------------------------------------

void compile_selection_expr(SelectionExpr& expr,
                            const std::vector<std::string>& level_set_names)
{
    for (auto& term : expr.terms)
    {
        term.zero_required     = 0;
        term.negative_required = 0;
        term.positive_required = 0;

        for (auto& clause : term.clauses)
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
                case Relation::LessThan:
                    term.negative_required |= bit;
                    break;
                case Relation::GreaterThan:
                    term.positive_required |= bit;
                    break;
                case Relation::EqualTo:
                    term.zero_required |= bit;
                    break;
            }
        }
    }

    mirror_first_term(expr);
}

// ---------------------------------------------------------------------------
// infer_selection_dim
// ---------------------------------------------------------------------------

int infer_selection_dim(const SelectionExpr& expr, int tdim)
{
    if (expr.terms.empty())
        throw std::runtime_error("infer_selection_dim: empty selection expression");

    auto infer_term_dim = [tdim](const SelectionTerm& term)
    {
        int n_eq = 0;
        for (const auto& clause : term.clauses)
            if (clause.relation == Relation::EqualTo)
                ++n_eq;
        return std::max(0, tdim - n_eq);
    };

    const int dim = infer_term_dim(expr.terms.front());
    for (std::size_t i = 1; i < expr.terms.size(); ++i)
    {
        if (infer_term_dim(expr.terms[i]) != dim)
        {
            throw std::runtime_error(
                "infer_selection_dim: 'or' terms must select the same entity dimension");
        }
    }
    return dim;
}

} // namespace cutcells
