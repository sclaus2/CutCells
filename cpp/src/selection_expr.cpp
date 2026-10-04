// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "selection_expr.h"

#include <algorithm>
#include <bit>
#include <cctype>
#include <stdexcept>

namespace cutcells
{
namespace
{

/// Expansions of level sets made of level sets stop beyond this many terms.
constexpr std::size_t max_terms = 4096;

bool is_space(char c)
{
    return std::isspace(static_cast<unsigned char>(c)) != 0;
}

/// Whitespace and the symbols of the grammar end a word.
bool ends_word(char c)
{
    return is_space(c) || c == '(' || c == ')' || c == '<' || c == '>' || c == '=';
}

/// Case-insensitive comparison of two words.
bool same_word(std::string_view a, std::string_view b)
{
    if (a.size() != b.size())
        return false;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        if (std::tolower(static_cast<unsigned char>(a[i])) != std::tolower(static_cast<unsigned char>(b[i])))
            return false;
    }
    return true;
}

bool is_keyword(std::string_view word)
{
    return same_word(word, "and") || same_word(word, "or");
}

bool is_less(Relation r)
{
    return r == Relation::LessThan || r == Relation::LessEqual;
}

std::string joined(const std::vector<std::string>& names)
{
    std::string s;
    for (const std::string& n : names)
        s += (s.empty() ? "" : ", ") + n;
    return s;
}

// ============================================================================
// Parsing
// ============================================================================

struct Token
{
    enum class Kind : std::uint8_t
    {
        word,
        open,
        close,
        relation,
        end
    };
    Kind kind = Kind::end;
    std::string_view text;
    Relation relation = Relation::LessThan;
};

std::vector<Token> tokenize(std::string_view text)
{
    std::vector<Token> tokens;
    std::size_t i = 0;
    while (i < text.size())
    {
        const char c = text[i];
        if (is_space(c))
        {
            ++i;
            continue;
        }
        if (c == '(' || c == ')')
        {
            tokens.push_back({c == '(' ? Token::Kind::open : Token::Kind::close, text.substr(i, 1)});
            ++i;
            continue;
        }
        if (c == '<' || c == '>' || c == '=')
        {
            // '<=', '>=' and '==' as well; '=' and '==' are the same relation
            const std::size_t n = i + 1 < text.size() && text[i + 1] == '=' ? 2 : 1;
            Relation r = Relation::EqualTo;
            if (c == '<')
                r = n == 2 ? Relation::LessEqual : Relation::LessThan;
            else if (c == '>')
                r = n == 2 ? Relation::GreaterEqual : Relation::GreaterThan;
            tokens.push_back({Token::Kind::relation, text.substr(i, n), r});
            i += n;
            continue;
        }
        std::size_t j = i;
        while (j < text.size() && !ends_word(text[j]))
            ++j;
        tokens.push_back({Token::Kind::word, text.substr(i, j - i)});
        i = j;
    }
    tokens.push_back({Token::Kind::end, text.substr(text.size())});
    return tokens;
}

/// A parsed selection: a clause, or the and (all) or or (any) of selections.
struct Node
{
    enum class Kind : std::uint8_t
    {
        clause,
        all,
        any
    };
    Kind kind = Kind::clause;
    Clause clause;
    std::vector<Node> children;
};

/// The tokens of a selection and the next one to read.
struct Tokens
{
    std::vector<Token> tokens;
    std::size_t next = 0;
    std::string_view text; ///< for messages
};

[[noreturn]] void syntax_error(const Tokens& t, const std::string& what)
{
    throw std::runtime_error("parse_selection_expr: " + what + " in '" + std::string(t.text) + "'");
}

bool at_keyword(const Tokens& t, std::string_view keyword)
{
    const Token& tok = t.tokens[t.next];
    return tok.kind == Token::Kind::word && same_word(tok.text, keyword);
}

Node parse_any(Tokens& t);

/// factor: clause | '(' expression ')'
Node parse_factor(Tokens& t)
{
    const Token& first = t.tokens[t.next];
    if (first.kind == Token::Kind::open)
    {
        ++t.next;
        Node inner = parse_any(t);
        if (t.tokens[t.next].kind != Token::Kind::close)
            syntax_error(t, "missing ')'");
        ++t.next;
        return inner;
    }
    if (first.kind != Token::Kind::word || is_keyword(first.text))
        syntax_error(t, "missing level-set name in a clause");
    Node node;
    node.clause.name = std::string(first.text);
    ++t.next;
    const Token& op = t.tokens[t.next];
    if (op.kind != Token::Kind::relation)
        syntax_error(t, "no operator (<, <=, >, >=, =, ==) after '" + node.clause.name + "'");
    node.clause.relation = op.relation;
    ++t.next;
    const Token& rhs = t.tokens[t.next];
    if (rhs.kind != Token::Kind::word || rhs.text != "0")
        syntax_error(t, "right-hand side must be '0', got '" + std::string(rhs.text) + "'");
    ++t.next;
    return node;
}

/// term: factor ('and' factor)*
Node parse_all(Tokens& t)
{
    Node first = parse_factor(t);
    if (!at_keyword(t, "and"))
        return first;
    Node all;
    all.kind = Node::Kind::all;
    all.children.push_back(std::move(first));
    while (at_keyword(t, "and"))
    {
        ++t.next;
        all.children.push_back(parse_factor(t));
    }
    return all;
}

/// expression: term ('or' term)*
Node parse_any(Tokens& t)
{
    Node first = parse_all(t);
    if (!at_keyword(t, "or"))
        return first;
    Node any;
    any.kind = Node::Kind::any;
    any.children.push_back(std::move(first));
    while (at_keyword(t, "or"))
    {
        ++t.next;
        any.children.push_back(parse_all(t));
    }
    return any;
}

Node parse_tree(std::string_view text)
{
    Tokens t{tokenize(text), 0, text};
    if (t.tokens.size() == 1)
        throw std::runtime_error("parse_selection_expr: empty expression");
    Node root = parse_any(t);
    const Token& rest = t.tokens[t.next];
    if (rest.kind == Token::Kind::close)
        syntax_error(t, "unexpected ')'");
    if (rest.kind != Token::Kind::end)
        syntax_error(t, "expected 'and' or 'or' before '" + std::string(rest.text) + "'");
    return root;
}

/// The or-terms of and-clauses of a parse tree, parentheses distributed, the
/// clauses in the order written.
std::vector<SelectionTerm> distribute(const Node& node)
{
    if (node.kind == Node::Kind::clause)
    {
        SelectionTerm term;
        term.clauses.push_back(node.clause);
        return {term};
    }
    std::vector<SelectionTerm> out;
    if (node.kind == Node::Kind::any)
    {
        for (const Node& child : node.children)
        {
            std::vector<SelectionTerm> terms = distribute(child);
            out.insert(out.end(), terms.begin(), terms.end());
        }
        return out;
    }
    out.resize(1);
    for (const Node& child : node.children)
    {
        const std::vector<SelectionTerm> right = distribute(child);
        if (out.size() * right.size() > max_terms)
            throw std::runtime_error("parse_selection_expr: the expression has more than "
                                     + std::to_string(max_terms) + " terms once its parentheses distribute");
        std::vector<SelectionTerm> next;
        next.reserve(out.size() * right.size());
        for (const SelectionTerm& a : out)
            for (const SelectionTerm& b : right)
            {
                SelectionTerm term = a;
                term.clauses.insert(term.clauses.end(), b.clauses.begin(), b.clauses.end());
                next.push_back(std::move(term));
            }
        out = std::move(next);
    }
    return out;
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

// ============================================================================
// Terms as bitmasks: and, or and simplification
// ============================================================================

bool contradictory(const SelectionTerm& t)
{
    return (t.negative_required & t.positive_required) != 0
           || (t.zero_required & (t.negative_required | t.positive_required)) != 0;
}

/// Whether term @p a holds wherever term @p b does: they lie in the same zero
/// sets, with the same normals, and a's signs are among b's.
bool contains(const SelectionTerm& a, const SelectionTerm& b)
{
    return a.zero_required == b.zero_required && a.normal_flipped == b.normal_flipped
           && (a.negative_required & ~b.negative_required) == 0 && (a.positive_required & ~b.positive_required) == 0;
}

/// Adds @p term to the or of @p terms, unless one of them contains it, and
/// drops those it contains.
void add_term(std::vector<SelectionTerm>& terms, const SelectionTerm& term)
{
    for (const SelectionTerm& t : terms)
        if (contains(t, term))
            return;
    std::erase_if(terms, [&term](const SelectionTerm& t) { return contains(term, t); });
    terms.push_back(term);
    if (terms.size() > max_terms)
        throw std::runtime_error("compile_selection_expr: the selection expands to more than "
                                 + std::to_string(max_terms) + " terms");
}

/// The and of two ors of terms, without contradictions.
std::vector<SelectionTerm> product(const std::vector<SelectionTerm>& a, const std::vector<SelectionTerm>& b)
{
    std::vector<SelectionTerm> out;
    for (const SelectionTerm& s : a)
        for (const SelectionTerm& t : b)
        {
            SelectionTerm u;
            u.zero_required = s.zero_required | t.zero_required;
            u.negative_required = s.negative_required | t.negative_required;
            u.positive_required = s.positive_required | t.positive_required;
            u.normal_flipped = s.normal_flipped | t.normal_flipped;
            if (!contradictory(u))
                add_term(out, u);
        }
    return out;
}

/// Clauses from a term's bitmasks: its zero sets first, then the other level
/// sets in order.
void write_clauses(SelectionTerm& term, const std::vector<std::string>& names)
{
    term.clauses.clear();
    for (std::uint64_t bits = term.zero_required; bits != 0; bits &= bits - 1)
    {
        const int l = std::countr_zero(bits);
        term.clauses.push_back({names[static_cast<std::size_t>(l)], l, Relation::EqualTo});
    }
    for (std::uint64_t bits = term.negative_required | term.positive_required; bits != 0; bits &= bits - 1)
    {
        const int l = std::countr_zero(bits);
        const Relation r = (term.negative_required >> l) & 1 ? Relation::LessThan : Relation::GreaterThan;
        term.clauses.push_back({names[static_cast<std::size_t>(l)], l, r});
    }
}

// ============================================================================
// Names and the expansion of level sets made of level sets
// ============================================================================

/// What a selection may name: the level sets of a cut, and the insides of
/// the level sets made of level sets it replaced by theirs, parsed once. A
/// name in the inside of X stands for "X.name", so no level set is made of
/// itself.
struct Registry
{
    const std::vector<std::string>& level_sets;
    const std::map<std::string, std::string>& insides;
    std::map<std::string, Node> parsed;
};

int index_of(const std::vector<std::string>& names, std::string_view name)
{
    const auto it = std::find(names.begin(), names.end(), name);
    return it == names.end() ? -1 : static_cast<int>(it - names.begin());
}

/// A name as written in a selection, resolved: the level set or level set
/// made of level sets of that name, or else the one named "X.name", X made of
/// level sets, if there is exactly one.
std::string resolve(const Registry& r, const std::string& name)
{
    if (index_of(r.level_sets, name) >= 0 || r.insides.count(name) > 0)
        return name;
    std::vector<std::string> found;
    auto consider = [&](const std::string& q)
    {
        if (q.size() > name.size() + 1 && q.ends_with(name) && q[q.size() - name.size() - 1] == '.'
            && r.insides.count(q.substr(0, q.size() - name.size() - 1)) > 0)
            found.push_back(q);
    };
    for (const std::string& q : r.level_sets)
        consider(q);
    for (const auto& entry : r.insides)
        consider(entry.first);
    if (found.empty())
        throw std::runtime_error("compile_selection_expr: unknown level-set name '" + name + "'");
    if (found.size() > 1)
        throw std::runtime_error("compile_selection_expr: the short name '" + name + "' could be any of "
                                 + joined(found) + "; give the full name");
    return found.front();
}

std::uint64_t level_set_bit(int index)
{
    if (index >= 64)
        throw std::runtime_error("compile_selection_expr: level-set index >= 64 not supported");
    return std::uint64_t(1) << index;
}

std::vector<SelectionTerm> expand_inside(Registry& r, const std::string& name, bool negated, int fixed, int sign);

/// The terms where @p node, a selection in the inside of @p owner, holds, or
/// with @p negated where it does not. Its names are those of owner's level
/// sets, without "owner.". Level set @p fixed, unless -1, has the sign
/// @p sign there.
std::vector<SelectionTerm> expand(Registry& r, const Node& node, const std::string& owner, bool negated, int fixed,
                                  int sign)
{
    if (node.kind == Node::Kind::clause)
    {
        const Clause& c = node.clause;
        if (c.relation == Relation::EqualTo)
            throw std::invalid_argument("compile_selection_expr: the inside of '" + owner + "' asks for '" + c.name
                                        + " = 0'; an inside is a region, of clauses < 0 and > 0");
        // where the clause holds (or, negated, fails), its level set is negative
        const bool negative = is_less(c.relation) != negated;
        const std::string name = owner + "." + c.name;
        const int l = index_of(r.level_sets, name);
        if (l >= 0)
        {
            if (l == fixed)
                return negative == (sign < 0) ? std::vector<SelectionTerm>(1) : std::vector<SelectionTerm>();
            SelectionTerm t;
            (negative ? t.negative_required : t.positive_required) = level_set_bit(l);
            return {t};
        }
        if (r.insides.count(name) > 0)
            return expand_inside(r, name, !negative, fixed, sign);
        throw std::invalid_argument("compile_selection_expr: the inside of '" + owner + "' names '" + c.name
                                    + "', none of its level sets");
    }
    // and of the children, or their or; negation swaps the two
    const bool all = (node.kind == Node::Kind::all) != negated;
    std::vector<SelectionTerm> out = all ? std::vector<SelectionTerm>(1) : std::vector<SelectionTerm>();
    for (const Node& child : node.children)
    {
        const std::vector<SelectionTerm> terms = expand(r, child, owner, negated, fixed, sign);
        if (all)
            out = product(out, terms);
        else
            for (const SelectionTerm& t : terms)
                add_term(out, t);
    }
    return out;
}

/// The terms where the level set made of level sets @p name is negative
/// (where its inside holds), or with @p negated positive.
std::vector<SelectionTerm> expand_inside(Registry& r, const std::string& name, bool negated, int fixed, int sign)
{
    auto it = r.parsed.find(name);
    if (it == r.parsed.end())
    {
        try
        {
            it = r.parsed.emplace(name, parse_tree(r.insides.at(name))).first;
        }
        catch (const std::runtime_error& e)
        {
            throw std::invalid_argument("compile_selection_expr: the inside of '" + name + "': " + e.what());
        }
    }
    return expand(r, it->second, name, negated, fixed, sign);
}

/// The faces of the level set made of level sets @p name, "name = 0": for each
/// level set a that its inside F depends on, a = 0 where F with a negative
/// differs from F with a positive, normal_flipped where F holds with a
/// positive.
std::vector<SelectionTerm> faces(Registry& r, const std::string& name)
{
    std::uint64_t support = 0;
    for (const SelectionTerm& t : expand_inside(r, name, false, -1, 0))
        support |= t.negative_required | t.positive_required;
    std::vector<SelectionTerm> out;
    for (std::uint64_t bits = support; bits != 0; bits &= bits - 1)
    {
        const int a = std::countr_zero(bits);
        const std::uint64_t bit = std::uint64_t(1) << a;
        const std::vector<SelectionTerm> below = expand_inside(r, name, false, a, -1);
        const std::vector<SelectionTerm> above = expand_inside(r, name, false, a, 1);
        for (SelectionTerm t : product(below, expand_inside(r, name, true, a, 1)))
        {
            t.zero_required |= bit;
            add_term(out, t);
        }
        for (SelectionTerm t : product(above, expand_inside(r, name, true, a, -1)))
        {
            t.zero_required |= bit;
            t.normal_flipped |= bit;
            add_term(out, t);
        }
    }
    return out;
}

/// Bitmasks of a term on level sets alone, from its resolved clauses.
void compile_term(SelectionTerm& term)
{
    term.zero_required = 0;
    term.negative_required = 0;
    term.positive_required = 0;
    term.normal_flipped = 0;
    for (const Clause& clause : term.clauses)
    {
        const std::uint64_t bit = level_set_bit(clause.level_set_index);
        switch (clause.relation)
        {
        case Relation::LessThan:
        case Relation::LessEqual:
            term.negative_required |= bit;
            break;
        case Relation::GreaterThan:
        case Relation::GreaterEqual:
            term.positive_required |= bit;
            break;
        case Relation::EqualTo:
            term.zero_required |= bit;
            break;
        }
    }
}

/// Checks the clauses of an inside against the names of its level sets.
void check_clauses(const Node& node, const std::vector<std::string>& names, std::string_view owner)
{
    if (node.kind != Node::Kind::clause)
    {
        for (const Node& child : node.children)
            check_clauses(child, names, owner);
        return;
    }
    const Clause& c = node.clause;
    if (index_of(names, c.name) < 0)
        throw std::invalid_argument("the inside of '" + std::string(owner) + "' names '" + c.name
                                    + "', none of its level sets (" + joined(names) + ")");
    if (c.relation == Relation::EqualTo)
        throw std::invalid_argument("the inside of '" + std::string(owner) + "' asks for '" + c.name
                                    + " = 0'; an inside is a region, of clauses < 0 and > 0");
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// parse_selection_expr
// ---------------------------------------------------------------------------

SelectionExpr parse_selection_expr(std::string_view text)
{
    SelectionExpr expr;
    expr.terms = distribute(parse_tree(text));
    mirror_first_term(expr);
    return expr;
}

// ---------------------------------------------------------------------------
// compile_selection_expr
// ---------------------------------------------------------------------------

void compile_selection_expr(SelectionExpr& expr,
                            const std::vector<std::string>& level_set_names)
{
    compile_selection_expr(expr, level_set_names, {});
}

void compile_selection_expr(SelectionExpr& expr, const std::vector<std::string>& level_set_names,
                            const std::map<std::string, std::string>& insides)
{
    Registry r{level_set_names, insides, {}};
    std::vector<SelectionTerm> terms;
    for (const SelectionTerm& written : expr.terms)
    {
        SelectionTerm term = written;
        bool expands = false;
        for (Clause& clause : term.clauses)
        {
            clause.name = resolve(r, clause.name);
            clause.level_set_index = index_of(level_set_names, clause.name);
            expands |= clause.level_set_index < 0;
        }
        if (!expands)
        {
            compile_term(term);
            terms.push_back(std::move(term));
            continue;
        }
        // a term on level sets made of level sets: the and of its clauses' expansions
        std::vector<SelectionTerm> expanded(1);
        for (const Clause& clause : term.clauses)
        {
            std::vector<SelectionTerm> right;
            if (clause.level_set_index >= 0)
            {
                SelectionTerm t;
                t.clauses.push_back(clause);
                compile_term(t);
                t.clauses.clear();
                right.push_back(t);
            }
            else if (clause.relation == Relation::EqualTo)
                right = faces(r, clause.name);
            else
                right = expand_inside(r, clause.name, !is_less(clause.relation), -1, 0);
            expanded = product(expanded, right);
        }
        for (SelectionTerm& t : expanded)
        {
            write_clauses(t, level_set_names);
            terms.push_back(std::move(t));
        }
        if (terms.size() > max_terms)
            throw std::runtime_error("compile_selection_expr: the selection expands to more than "
                                     + std::to_string(max_terms) + " terms");
    }
    expr.terms = std::move(terms);
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

// ---------------------------------------------------------------------------
// Names and insides
// ---------------------------------------------------------------------------

bool is_selection_name(std::string_view name)
{
    if (name.empty() || is_keyword(name))
        return false;
    return std::none_of(name.begin(), name.end(), ends_word);
}

void check_inside(std::string_view inside, const std::vector<std::string>& names, std::string_view owner)
{
    Node tree;
    try
    {
        tree = parse_tree(inside);
    }
    catch (const std::runtime_error& e)
    {
        throw std::invalid_argument("the inside of '" + std::string(owner) + "': " + e.what());
    }
    check_clauses(tree, names, owner);
}

} // namespace cutcells
