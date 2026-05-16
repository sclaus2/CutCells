// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace cutcells
{

/// Relation in a selection clause: φ < 0, φ > 0, or φ = 0.
enum class Relation
{
    LessThan,
    GreaterThan,
    EqualTo
};

/// One clause of a selection expression, e.g. "phi1 < 0".
struct Clause
{
    std::string name;              ///< level-set name (e.g. "phi1")
    int         level_set_index = -1; ///< resolved after compile()
    Relation    relation;          ///< the comparison operator
};

/// One conjunction term in a selection expression.
///
/// A term is an AND of clauses:
///   "phi1 < 0 and phi2 = 0" -> two clauses in one term.
///
/// After compilation against a name registry, bitmasks are populated:
///   - bit i set in zero_required     → φ_i = 0
///   - bit i set in negative_required → φ_i < 0
///   - bit i set in positive_required → φ_i > 0
struct SelectionTerm
{
    std::vector<Clause> clauses;          ///< implicitly AND-connected

    // Compiled bitmasks (set by compile_selection_expr())
    std::uint64_t zero_required     = 0;
    std::uint64_t negative_required = 0;
    std::uint64_t positive_required = 0;
};

/// A compiled selection expression.
///
/// The expression is a union (OR) of terms:
///   "phi1 < 0 or phi2 < 0" -> two terms.
///   "phi1 < 0 and phi2 > 0" -> one term with two clauses.
struct SelectionExpr
{
    std::vector<SelectionTerm> terms;

    // Legacy single-term view. These fields mirror terms.front() after
    // parsing/compilation so existing single-term code can be migrated
    // incrementally. New selection code should use terms.
    std::vector<Clause> clauses;
    std::uint64_t zero_required     = 0;
    std::uint64_t negative_required = 0;
    std::uint64_t positive_required = 0;
};

/// Parse a selection expression string.
///
/// Grammar:  term ('or' term)*
/// term:     clause ('and' clause)*
/// clause:   name ('<'|'>'|'=') '0'
///
/// @throws std::runtime_error on syntax error.
SelectionExpr parse_selection_expr(std::string_view text);

/// Compile a parsed expression against a level-set name registry.
///
/// Resolves each clause's name to its index in @p level_set_names
/// and populates the bitmask fields.
///
/// @throws std::runtime_error if a name is not found.
void compile_selection_expr(SelectionExpr& expr,
                            const std::vector<std::string>& level_set_names);

/// Infer the entity dimension selected by a compiled expression.
///
/// Each OR term must select the same dimension. Within a term, each EqualTo
/// clause reduces the dimension by one; a term with no EqualTo clauses selects
/// volume (tdim).
///
/// @param tdim  topological dimension of the background cells.
/// @return inferred entity dimension.
int infer_selection_dim(const SelectionExpr& expr, int tdim);

} // namespace cutcells
