// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace cutcells
{

/// Relation in a selection clause: φ < 0, φ <= 0, φ > 0, φ >= 0, or φ = 0
/// (also written φ == 0).
enum class Relation
{
    LessThan,
    LessEqual,
    GreaterThan,
    GreaterEqual,
    EqualTo
};

/// One clause of a selection expression, e.g. "phi1 < 0".
struct Clause
{
    std::string name;              ///< level-set name (e.g. "phi1")
    int         level_set_index = -1; ///< resolved after compile()
    Relation    relation = Relation::LessThan; ///< the comparison operator
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
///   - bit i set in normal_flipped    → on φ_i = 0, the outward normal of the
///     selected side is -∇φ_i/|∇φ_i|: the side lies where φ_i > 0. Clear for
///     plain clauses "φ_i = 0" (+∇φ_i/|∇φ_i|); set on the faces of a level
///     set made of level sets whose inside holds on the positive side of φ_i.
struct SelectionTerm
{
    std::vector<Clause> clauses;          ///< implicitly AND-connected

    // Compiled bitmasks (set by compile_selection_expr())
    std::uint64_t zero_required     = 0;
    std::uint64_t negative_required = 0;
    std::uint64_t positive_required = 0;
    std::uint64_t normal_flipped    = 0;
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
/// Grammar:  expression: term ('or' term)*
/// term:     factor ('and' factor)*
/// factor:   clause | '(' expression ')'
/// clause:   name ('<'|'<='|'>'|'>='|'='|'==') '0'
///
/// 'and' binds more tightly than 'or'; both are case-insensitive. A name is a
/// word without whitespace, parentheses, '<', '>' or '=' (is_selection_name).
/// Parentheses distribute: "(a < 0 or b < 0) and c > 0" gives the terms
/// "a < 0 and c > 0" and "b < 0 and c > 0", clauses in the order written.
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

/// @brief Compile a parsed expression against the level sets of a cut and the
/// level sets made of level sets it replaced by theirs.
///
/// A level set X made of level sets is replaced by its own, named "X.a",
/// "X.b", ... in @p level_set_names. @p insides gives, for each such X by its
/// name, its inside F: the selection of its own level sets, named without
/// "X.", where X is negative. Its level sets may be made of level sets in
/// turn, with names and insides of the same form ("X.a" in @p insides).
///
/// A name resolves to the level set or the level set made of level sets of
/// that name, or else, as a short name, to the one named "X.name" for some X
/// made of level sets, if there is exactly one.
///
/// Terms that name no level set made of level sets compile as the overload
/// above. In the others, "X < 0" (or <= 0) expands to F, "X > 0" (or >= 0)
/// to its opposite (and and or, < and > swapped), and "X = 0" to the faces
/// of X: for each level set a that F depends on, "a = 0" where F with a
/// negative differs from F with a positive, that is "a = 0 and F(a < 0) and
/// not F(a > 0)", and, with a's bit in normal_flipped, "a = 0 and F(a > 0)
/// and not F(a < 0)". A plane that F uses with both signs, such as the face
/// two boxes of a union share, gives no face where both sides hold. The
/// expanded terms drop contradictions and terms that others contain, and
/// get clauses from their bitmasks: the zero clause first, then the level
/// sets in order.
///
/// @throws std::runtime_error for unknown or ambiguous names, and for
///         expansions beyond 4096 terms
/// @throws std::invalid_argument for insides that do not parse, name other
///         level sets or ask for "= 0"
void compile_selection_expr(SelectionExpr& expr, const std::vector<std::string>& level_set_names,
                            const std::map<std::string, std::string>& insides);

/// Infer the entity dimension selected by a compiled expression.
///
/// Each OR term must select the same dimension. Within a term, each EqualTo
/// clause reduces the dimension by one; a term with no EqualTo clauses selects
/// volume (tdim).
///
/// @param tdim  topological dimension of the background cells.
/// @return inferred entity dimension.
int infer_selection_dim(const SelectionExpr& expr, int tdim);

/// @brief Whether @p name can name a level set in a selection: a word, not
/// empty, without whitespace, parentheses, '<', '>' or '=', and other than
/// 'and' and 'or'.
bool is_selection_name(std::string_view name);

/// @brief Check the inside of a level set made of level sets: a selection
/// whose clauses name level sets in @p names and ask for "< 0" or "> 0"
/// ("<= 0" and ">= 0" read the same).
/// @param owner  the level set's name, for messages
/// @throws std::invalid_argument otherwise, also for syntax errors
void check_inside(std::string_view inside, const std::vector<std::string>& names, std::string_view owner);

} // namespace cutcells
