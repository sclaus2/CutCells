// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Level sets made of level sets (LevelSetFunction::level_sets and inside):
//  - "X < 0", "X > 0" and "X = 0" expand into terms that hold exactly where
//    the inside holds, fails, and changes across one of its level sets (with
//    normal_flipped where it holds on the positive side), at every assignment
//    of signs: insides with nested parentheses, level sets made of level sets
//    made of level sets, a plane used with both signs; and the terms of the
//    napkin ring and of a nested inside, written out;
//  - selections without them compile as before; parentheses distribute in
//    the order written; short names resolve when unambiguous; syntax errors,
//    bad insides and expansions beyond the cap raise;
//  - part::cut replaces them by their own level sets, named after them, and
//    keeps their insides; bad names and insides raise;
//  - on [-1, 1]^3 with hexahedra and tetrahedra: the napkin ring made of a
//    ball and a cylinder has the totals of its hand-written selections, and
//    on("cyl") gives 2 pi r h, on("ball") 2 pi R h; two boxes side by side
//    whose shared plane is one level set used with both signs have no surface
//    on it, exact volumes and surfaces, and a smaller second box leaves the
//    ring around it on the plane.
// Exits non-zero on failure.

#include <cutcells/level_set.h>
#include <cutcells/part/cut_result.h>
#include <cutcells/part/mesh_part.h>
#include <cutcells/part/output.h>
#include <cutcells/quadrays/analytic.h>
#include <cutcells/selection_expr.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "support/box_mesh.h"

using namespace cutcells;
using namespace cutcells::part;
using namespace cutcells::part::support;

namespace
{

int failures = 0;

void check(bool ok, const std::string& what)
{
    if (!ok)
    {
        std::printf("FAILED: %s\n", what.c_str());
        ++failures;
    }
}

using Names = std::vector<std::string>;
using Insides = std::map<std::string, std::string>;

std::vector<SelectionTerm> compiled(const std::string& text, const Names& names, const Insides& insides)
{
    SelectionExpr expr = parse_selection_expr(text);
    compile_selection_expr(expr, names, insides);
    return expr.terms;
}

/// A term as text, its clauses as written, a flipped zero set marked.
std::string text_of(const SelectionTerm& t)
{
    std::string s;
    for (const Clause& c : t.clauses)
    {
        const char* op = c.relation == Relation::LessThan       ? " < 0"
                         : c.relation == Relation::LessEqual    ? " <= 0"
                         : c.relation == Relation::GreaterThan  ? " > 0"
                         : c.relation == Relation::GreaterEqual ? " >= 0"
                                                                : " = 0";
        s += (s.empty() ? "" : " and ") + c.name + op;
        if (c.relation == Relation::EqualTo && ((t.normal_flipped >> c.level_set_index) & 1))
            s += " (flipped)";
    }
    return s;
}

std::string joined(const std::vector<std::string>& texts)
{
    std::string s;
    for (const std::string& t : texts)
        s += (s.empty() ? "" : " | ") + t;
    return s;
}

/// The terms of @p text, compiled, against @p expected, in any order.
void expect_terms(const std::string& text, const Names& names, const Insides& insides,
                  std::vector<std::string> expected)
{
    std::vector<std::string> got;
    for (const SelectionTerm& t : compiled(text, names, insides))
        got.push_back(text_of(t));
    std::sort(got.begin(), got.end());
    std::sort(expected.begin(), expected.end());
    check(got == expected, "'" + text + "' gives " + joined(got) + ", expected " + joined(expected));
}

bool term_holds(const SelectionTerm& t, const std::vector<int>& signs)
{
    for (std::size_t l = 0; l < signs.size(); ++l)
    {
        const std::uint64_t bit = std::uint64_t(1) << l;
        if ((t.zero_required & bit) && signs[l] != 0)
            return false;
        if ((t.negative_required & bit) && signs[l] != -1)
            return false;
        if ((t.positive_required & bit) && signs[l] != 1)
            return false;
    }
    return true;
}

/// "x < 0", "x > 0" and "x = 0" against @p holds, the inside of x, at every
/// assignment of signs to the level sets @p names: below where it holds,
/// above where it fails, on a = 0 where it changes as a does, flipped where
/// it holds with a positive.
void check_expansion(const std::string& x, const Names& names, const Insides& insides,
                     const std::function<bool(const std::vector<int>&)>& holds)
{
    const int n = static_cast<int>(names.size());
    const std::vector<SelectionTerm> below = compiled(x + " < 0", names, insides);
    const std::vector<SelectionTerm> above = compiled(x + " > 0", names, insides);
    const std::vector<SelectionTerm> faces = compiled(x + " = 0", names, insides);
    auto any = [](const std::vector<SelectionTerm>& terms, const std::vector<int>& s)
    { return std::any_of(terms.begin(), terms.end(), [&s](const SelectionTerm& t) { return term_holds(t, s); }); };
    bool ok = true;
    std::vector<int> s(static_cast<std::size_t>(n));
    for (int m = 0; m < (1 << n); ++m)
    {
        for (int l = 0; l < n; ++l)
            s[l] = (m >> l) & 1 ? 1 : -1;
        const bool in = holds(s);
        ok &= any(below, s) == in && any(above, s) == !in && !any(faces, s);
    }
    for (int a = 0; a < n; ++a)
        for (int m = 0; m < (1 << n); ++m)
        {
            if ((m >> a) & 1)
                continue;
            for (int l = 0; l < n; ++l)
                s[l] = l == a ? 0 : ((m >> l) & 1 ? 1 : -1);
            auto with = [&](int sign)
            {
                std::vector<int> t = s;
                t[a] = sign;
                return holds(t);
            };
            bool found = false;
            for (const SelectionTerm& t : faces)
                if (term_holds(t, s))
                {
                    found = true;
                    ok &= ((t.normal_flipped >> a) & 1) == static_cast<std::uint64_t>(with(1));
                }
            ok &= found == (with(-1) != with(1));
        }
    for (const SelectionTerm& t : faces)
        ok &= std::popcount(t.zero_required) == 1 && (t.normal_flipped & ~t.zero_required) == 0;
    std::printf("expansion of %s: %zu terms below, %zu above, %zu faces\n", x.c_str(), below.size(), above.size(),
                faces.size());
    check(ok, "the expansion of " + x + " against its inside at every assignment of signs");
}

template <typename E>
bool raises(const std::function<void()>& f)
{
    try
    {
        f();
    }
    catch (const E&)
    {
        return true;
    }
    catch (...)
    {
        return false;
    }
    return false;
}

void test_expansions()
{
    // the napkin ring, written out
    const Names ring = {"ring.ball", "ring.cyl", "fluid"};
    const Insides ring_inside = {{"ring", "ball < 0 and cyl > 0"}};
    expect_terms("ring < 0", ring, ring_inside, {"ring.ball < 0 and ring.cyl > 0"});
    expect_terms("ring <= 0", ring, ring_inside, {"ring.ball < 0 and ring.cyl > 0"});
    expect_terms("ring > 0", ring, ring_inside, {"ring.ball > 0", "ring.cyl < 0"});
    expect_terms("ring = 0", ring, ring_inside,
                 {"ring.ball = 0 and ring.cyl > 0", "ring.cyl = 0 (flipped) and ring.ball < 0"});
    expect_terms("ring == 0 and fluid < 0", ring, ring_inside,
                 {"ring.ball = 0 and ring.cyl > 0 and fluid < 0",
                  "ring.cyl = 0 (flipped) and ring.ball < 0 and fluid < 0"});
    expect_terms("(ring > 0 and fluid < 0) or fluid = 0", ring, ring_inside,
                 {"ring.ball > 0 and fluid < 0", "ring.cyl < 0 and fluid < 0", "fluid = 0"});
    // short names resolve, also in terms without level sets made of level sets
    expect_terms("cyl = 0 and ball < 0", ring, ring_inside, {"ring.cyl = 0 and ring.ball < 0"});
    check_expansion("ring", {"ring.ball", "ring.cyl"}, ring_inside,
                    [](const std::vector<int>& s) { return s[0] < 0 && s[1] > 0; });

    // a nested inside, written out
    const Names abc = {"body.a", "body.b", "body.c"};
    const Insides body = {{"body", "(a < 0 or b < 0) and c > 0"}};
    expect_terms("body < 0", abc, body, {"body.a < 0 and body.c > 0", "body.b < 0 and body.c > 0"});
    expect_terms("body > 0", abc, body, {"body.a > 0 and body.b > 0", "body.c < 0"});
    expect_terms("body = 0", abc, body,
                 {"body.a = 0 and body.b > 0 and body.c > 0", "body.b = 0 and body.a > 0 and body.c > 0",
                  "body.c = 0 (flipped) and body.a < 0", "body.c = 0 (flipped) and body.b < 0"});

    // deeper parentheses, against the inside itself
    const Names abcd = {"x.a", "x.b", "x.c", "x.d"};
    check_expansion("x", abcd, {{"x", "((a < 0 and b < 0) or c < 0) and (d > 0 or a > 0)"}},
                    [](const std::vector<int>& s)
                    { return ((s[0] < 0 && s[1] < 0) || s[2] < 0) && (s[3] > 0 || s[0] > 0); });
    check_expansion("x", abcd, {{"x", "a < 0 and (b < 0 or (c > 0 and (d < 0 or a > 0)))"}},
                    [](const std::vector<int>& s)
                    { return s[0] < 0 && (s[1] < 0 || (s[2] > 0 && (s[3] < 0 || s[0] > 0))); });
    check_expansion("x", abcd, {{"x", "(a < 0 or b > 0) and (b < 0 or c > 0) and (c < 0 or d > 0) or d < 0 and a > 0"}},
                    [](const std::vector<int>& s)
                    {
                        return ((s[0] < 0 || s[1] > 0) && (s[1] < 0 || s[2] > 0) && (s[2] < 0 || s[3] > 0))
                               || (s[3] < 0 && s[0] > 0);
                    });

    // level sets made of level sets made of level sets, and their short names
    const Names nested = {"outer.lens.a", "outer.lens.b", "outer.c", "outer.hole"};
    const Insides outer = {{"outer", "(lens < 0 or c < 0) and hole > 0"}, {"outer.lens", "a < 0 and b < 0"}};
    check_expansion("outer", nested, outer, [](const std::vector<int>& s)
                    { return ((s[0] < 0 && s[1] < 0) || s[2] < 0) && s[3] > 0; });
    check_expansion("outer.lens", nested, outer, [](const std::vector<int>& s) { return s[0] < 0 && s[1] < 0; });
    expect_terms("lens < 0", nested, outer, {"outer.lens.a < 0 and outer.lens.b < 0"});
    expect_terms("lens.a = 0", nested, outer, {"outer.lens.a = 0"});

    // two boxes side by side: their shared plane p is used with both signs
    const Names boxes = {"u.x1", "u.p", "u.x2", "u.y"};
    const Insides both = {{"u", "(x1 > 0 and p < 0 and y < 0) or (p > 0 and x2 < 0 and y < 0)"}};
    check_expansion("u", boxes, both, [](const std::vector<int>& s)
                    { return (s[0] > 0 && s[1] < 0 && s[3] < 0) || (s[1] > 0 && s[2] < 0 && s[3] < 0); });
    expect_terms("u = 0", boxes, both,
                 {"u.x1 = 0 (flipped) and u.p < 0 and u.y < 0", "u.p = 0 and u.x1 > 0 and u.x2 > 0 and u.y < 0",
                  "u.p = 0 (flipped) and u.x1 < 0 and u.x2 < 0 and u.y < 0", "u.x2 = 0 and u.p > 0 and u.y < 0",
                  "u.y = 0 and u.x1 > 0 and u.p < 0", "u.y = 0 and u.p > 0 and u.x2 < 0"});
    // a complement alone: its face has the opposite normal
    expect_terms("h = 0", {"h.a"}, {{"h", "a > 0"}}, {"h.a = 0 (flipped)"});
}

void test_plain_selections()
{
    const Names abc = {"a", "b", "c"};
    SelectionExpr e = parse_selection_expr("a < 0 and b > 0 or c = 0");
    compile_selection_expr(e, abc);
    bool same = e.terms.size() == 2 && e.terms[0].clauses.size() == 2 && e.terms[1].clauses.size() == 1
                && e.terms[0].negative_required == 1 && e.terms[0].positive_required == 2
                && e.terms[1].zero_required == 4 && e.terms[0].normal_flipped == 0 && e.terms[1].normal_flipped == 0
                && e.clauses.size() == 2 && e.negative_required == 1 && e.positive_required == 2;
    check(same, "a plain selection compiles as before");
    // parentheses distribute in the order written; relations stay as written
    expect_terms("(a < 0 or b <= 0) and c = 0", abc, {}, {"a < 0 and c = 0", "b <= 0 and c = 0"});
    e = parse_selection_expr("(a < 0 or b <= 0) and c = 0");
    check(e.terms.size() == 2 && e.terms[0].clauses[0].name == "a" && e.terms[0].clauses[1].name == "c"
              && e.terms[1].clauses[0].relation == Relation::LessEqual,
          "parentheses distribute in the order written");
    // contradictions stay in plain selections, as before
    expect_terms("a < 0 and a > 0", abc, {}, {"a < 0 and a > 0"});
    expect_terms("(a<0)and(b>0) OR c>=0", abc, {}, {"a < 0 and b > 0", "c >= 0"});
    expect_terms("((a < 0))", abc, {}, {"a < 0"});
    for (const char* bad : {"", "  ", "()", "(a < 0", "a < 0)", "a < 0 b < 0", "a <", "< 0", "a < 1", "a === 0",
                            "a < 0 and", "and < 0", "a < 0 or or b < 0", "a b < 0", "a < 00"})
        check(raises<std::runtime_error>([bad] { parse_selection_expr(bad); }),
              std::string("'") + bad + "' does not parse");
    check(raises<std::runtime_error>([&] { compiled("z < 0", abc, {}); }), "an unknown name raises");
    for (const char* name : {"a", "block.x-", "leg.2", "phi1", "x+"})
        check(is_selection_name(name), std::string("'") + name + "' is a selection name");
    for (const char* name : {"", "a b", "f(x)", "a<", "x=1", "and", "OR"})
        check(!is_selection_name(name), std::string("'") + name + "' is no selection name");
}

void test_bad_insides()
{
    const Names holes = {"s.hole", "t.hole"};
    const Insides two = {{"s", "hole < 0"}, {"t", "hole > 0"}};
    check(raises<std::runtime_error>([&] { compiled("hole < 0", holes, two); }), "an ambiguous short name raises");
    expect_terms("s.hole < 0 and t = 0", holes, two, {"t.hole = 0 (flipped) and s.hole < 0"});
    check(raises<std::invalid_argument>([] { compiled("x < 0", {"x.a"}, {{"x", "a = 0"}}); }),
          "an inside asking for = 0 raises");
    check(raises<std::invalid_argument>([] { compiled("x < 0", {"x.a"}, {{"x", "b < 0"}}); }),
          "an inside naming another level set raises");
    check(raises<std::invalid_argument>([] { compiled("x > 0", {"x.a"}, {{"x", "a < 0 and"}}); }),
          "an inside that does not parse raises");
    check(raises<std::invalid_argument>([] { check_inside("a < 0 and c > 0", {"a", "b"}, "x"); }),
          "check_inside: unknown names");
    check(raises<std::invalid_argument>([] { check_inside("a == 0", {"a"}, "x"); }), "check_inside: = 0");
    check(raises<std::invalid_argument>([] { check_inside("(a < 0", {"a"}, "x"); }), "check_inside: syntax");
    check(!raises<std::exception>([] { check_inside("(a <= 0 or b >= 0) and a < 0", {"a", "b"}, "x"); }),
          "check_inside: a good inside");
    // the union of five boxes: its outside has 6^5 terms, beyond the cap
    Names planes;
    std::string inside;
    for (int b = 0; b < 5; ++b)
    {
        std::string box;
        for (int f = 0; f < 6; ++f)
        {
            const std::string name = "f" + std::to_string(b) + std::to_string(f);
            planes.push_back("u." + name);
            box += (box.empty() ? "" : " and ") + name + " < 0";
        }
        inside += (inside.empty() ? "(" : " or (") + box + ")";
    }
    check(compiled("u < 0", planes, {{"u", inside}}).size() == 5, "the union of five boxes");
    check(raises<std::runtime_error>([&] { compiled("u > 0", planes, {{"u", inside}}); }),
          "an expansion beyond the cap raises");
}

/// The plane x_axis = offset.
struct AxisPlane
{
    int axis = 0;
    double offset = 0;

    template <typename V>
    V operator()(const std::array<V, 3>& x) const
    {
        return x[axis] - offset;
    }
};

/// |x - c|^2 - r^2, or, with axis >= 0, the cylinder about that axis.
struct Quadric
{
    V3 c = {0, 0, 0};
    double r = 0;
    int axis = -1;

    template <typename V>
    V operator()(const std::array<V, 3>& x) const
    {
        V q = V(-r * r);
        for (int i = 0; i < 3; ++i)
            if (i != axis)
                q = q + (x[i] - c[i]) * (x[i] - c[i]);
        return q;
    }
};

/// An analytic level set that keeps its functor, which it points to, alive.
template <typename F>
LevelSetFunction<double, int> analytic_ls(const F& functor, const std::string& name)
{
    struct Holder
    {
        F functor;
        quadrays::AnalyticLevelSet phi;
    };
    auto holder = std::make_shared<Holder>();
    holder->functor = functor;
    holder->phi = quadrays::analytic_level_set(holder->functor);
    const quadrays::AnalyticLevelSet* phi = &holder->phi;
    return create_level_set_function<double, int>(std::shared_ptr<const quadrays::AnalyticLevelSet>(holder, phi), 3,
                                                  name);
}

LevelSetFunction<double, int> made_of(std::vector<LevelSetFunction<double, int>> level_sets,
                                      const std::string& inside, const std::string& name)
{
    LevelSetFunction<double, int> ls;
    ls.name = name;
    ls.level_sets = std::move(level_sets);
    ls.inside = inside;
    return ls;
}

double total(const MeshPart<double, int>& part, int order, bool full)
{
    const auto rules = quadrature_rules(part, order, full);
    double s = 0;
    for (const double w : rules._weights)
        s += w;
    return s;
}

void test_cut()
{
    BoxMesh mesh;
    make_box_mesh("hex", 2, {0, 0, 0}, mesh);
    const LevelSetFunction<double, int> a = analytic_ls(AxisPlane{0, 0.1}, "a");
    const LevelSetFunction<double, int> b = analytic_ls(AxisPlane{1, 0.1}, "b");
    const std::vector<LevelSetFunction<double, int>> level_sets
        = {made_of({made_of({a, b}, "a < 0 and b < 0", "corner"), analytic_ls(AxisPlane{2, 0.1}, "c")},
                   "corner < 0 or c < 0", "body"),
           analytic_ls(AxisPlane{2, -0.3}, "fluid")};
    const CutResult<double, int> r = cut<double, int>(mesh.view, level_sets);
    check(r.level_set_names == Names{"body.corner.a", "body.corner.b", "body.c", "fluid"} && r.n_level_sets() == 4
              && r.insides == Insides{{"body", "corner < 0 or c < 0"}, {"body.corner", "a < 0 and b < 0"}},
          "cut replaces level sets made of level sets by their own, named after them");
    check(select(r, "body = 0").dim == 2 && select(r, "corner < 0 and fluid > 0").dim == 3,
          "selections of level sets made of level sets");

    auto cut_raises = [&mesh](std::vector<LevelSetFunction<double, int>> ls)
    { return raises<std::invalid_argument>([&] { cut<double, int>(mesh.view, ls); }); };
    check(cut_raises({a, a}), "two level sets of one name raise");
    check(cut_raises({made_of({a, a}, "a < 0", "x")}), "two level sets of one name within one raise");
    check(cut_raises({made_of({a, b}, "a < 0", "my body")}), "a name that selections cannot use raises");
    check(cut_raises({made_of({a, analytic_ls(AxisPlane{1, 0.1}, "b c")}, "a < 0", "x")}),
          "a name within that selections cannot use raises");
    check(cut_raises({made_of({a, b}, "a < 0", "x"), analytic_ls(AxisPlane{2, 0.1}, "x.c")}),
          "a name that reads as one of another's raises");
    check(cut_raises({made_of({a, b}, "", "x")}), "level sets without an inside raise");
    LevelSetFunction<double, int> lonely = a;
    lonely.inside = "a < 0";
    check(cut_raises({lonely}), "an inside without level sets raises");
    check(cut_raises({made_of({a, b}, "a < 0 and z > 0", "x")}), "an inside naming another level set raises");
    std::vector<LevelSetFunction<double, int>> many;
    for (int i = 0; i < 65; ++i)
        many.push_back(analytic_ls(AxisPlane{0, 0.01 * i}, "p" + std::to_string(i)));
    check(cut_raises({made_of(many, "p0 < 0", "x")}), "more than 64 level sets within raise");
    many.pop_back();
    check(!raises<std::exception>([&] { cut<double, int>(mesh.view, std::vector{made_of(many, "p0 < 0", "x")}); }),
          "64 level sets within");
}

/// The napkin ring, a ball of radius R outside a cylinder of radius r, made of
/// the two, against its hand-written selections and exact values (order 8:
/// 5 Gauss points per segment).
void test_napkin_ring()
{
    const V3 c = {0.0123, -0.0371, 0.0217};
    const double R = 0.8, r = 0.45, h = 2 * std::sqrt(R * R - r * r);
    const double volume = M_PI * h * h * h / 6, zone = 2 * M_PI * R * h, inner = 2 * M_PI * r * h;
    for (const char* kind : {"hex", "tet"})
    {
        BoxMesh mesh;
        make_box_mesh(kind, 8, {0, 0, 0}, mesh);
        const LevelSetFunction<double, int> ball = analytic_ls(Quadric{c, R, -1}, "ball");
        const LevelSetFunction<double, int> cyl = analytic_ls(Quadric{c, r, 2}, "cyl");
        const std::vector<LevelSetFunction<double, int>> plain_ls = {ball, cyl};
        const std::vector<LevelSetFunction<double, int>> ring_ls
            = {made_of({ball, cyl}, "ball < 0 and cyl > 0", "ring")};
        const CutResult<double, int> plain = cut<double, int>(mesh.view, plain_ls);
        const CutResult<double, int> made = cut<double, int>(mesh.view, ring_ls);
        const double v = total(select(made, "ring < 0"), 8, true);
        const double v0 = total(select(plain, "ball < 0 and cyl > 0"), 8, true);
        const double out = total(select(made, "ring > 0"), 8, true);
        const MeshPart<double, int> faces = select(made, "ring == 0");
        const double s = total(faces, 8, false);
        const double s0 = total(select(plain, "ball = 0 and cyl > 0"), 8, false)
                          + total(select(plain, "cyl = 0 and ball < 0"), 8, false);
        const double on_cyl = total(on(faces, {"cyl"}), 8, false), on_ball = total(on(faces, {"ball"}), 8, false);
        const double qualified = total(on(faces, {"ring.cyl"}), 8, false);
        const double both = total(on(faces, {"ball", "ring.cyl"}), 8, false);
        const double all = total(on(faces, {"ring"}), 8, false);
        std::printf("napkin ring %s: volume %.1e (hand-written %.1e), outside %.1e, faces %.1e (%.1e), "
                    "on cyl %.1e, on ball %.1e\n",
                    kind, v / volume - 1, v / v0 - 1, out / (8 - volume) - 1, s / (zone + inner) - 1, s / s0 - 1,
                    on_cyl / inner - 1, on_ball / zone - 1);
        check(std::abs(v / v0 - 1) < 1e-12 && std::abs(s / s0 - 1) < 1e-12 && std::abs(v / volume - 1) < 1e-6
                  && std::abs(out / (8 - volume) - 1) < 1e-7 && std::abs(on_cyl / inner - 1) < 1e-6
                  && std::abs(on_ball / zone - 1) < 1e-6 && qualified == on_cyl
                  && std::abs(both / s - 1) < 1e-12 && all == s,
              std::string("the napkin ring made of a ball and a cylinder, ") + kind);
        check(raises<std::invalid_argument>([&] { on(select(made, "ring < 0"), {"cyl"}); })
                  && raises<std::invalid_argument>([&] { on(faces, {"hole"}); })
                  && raises<std::invalid_argument>([&] { on(on(faces, {"ball"}), {"cyl"}); })
                  && raises<std::invalid_argument>([&] { on(faces, {}); }),
              "on raises for volumes, unknown names and names that are no faces");
    }
}

/// Two boxes side by side, [-0.7, 0.1] x Y x Z and [0.1, 0.7] x Y x Z with
/// Y = Z = [-0.6, 0.45], made of planes off the grid (h = 0.25); again with
/// the boxes meeting at x = 0, in mesh faces; and with a smaller second box,
/// Y = Z = [-0.35, 0.3].
void test_boxes()
{
    auto plane = [](int axis, double offset, const std::string& name)
    { return analytic_ls(AxisPlane{axis, offset}, name); };
    const std::string first = "xlo > 0 and p < 0 and ylo > 0 and yhi < 0 and zlo > 0 and zhi < 0";
    const double w = 1.05; // the width in y and z
    for (const char* kind : {"hex", "tet"})
    {
        BoxMesh mesh;
        make_box_mesh(kind, 8, {0, 0, 0}, mesh);
        for (const double shared : {0.1, 0.0})
        {
            const std::vector<LevelSetFunction<double, int>> ls
                = {made_of({plane(0, -0.7, "xlo"), plane(0, shared, "p"), plane(0, 0.7, "xhi"), plane(1, -0.6, "ylo"),
                            plane(1, 0.45, "yhi"), plane(2, -0.6, "zlo"), plane(2, 0.45, "zhi")},
                           "(" + first + ") or (p > 0 and xhi < 0 and ylo > 0 and yhi < 0 and zlo > 0 and zhi < 0)",
                           "union")};
            const CutResult<double, int> r = cut<double, int>(mesh.view, ls);
            const MeshPart<double, int> faces = select(r, "union = 0");
            const double v = total(select(r, "union < 0"), 2, true), s = total(faces, 2, false);
            const auto shared_rules = quadrature_rules(on(faces, {"p"}), 2, false);
            const double ends = total(on(faces, {"xlo", "xhi"}), 2, false);
            const double volume = 1.4 * w * w, surface = 2 * (2 * 1.4 * w + w * w);
            std::printf("two boxes %s, shared plane x = %.1f: volume %.1e, surface %.1e, %zu points on the plane, "
                        "ends %.1e\n",
                        kind, shared, v / volume - 1, s / surface - 1, shared_rules._weights.size(),
                        ends / (2 * w * w) - 1);
            check(std::abs(v / volume - 1) < 1e-12 && std::abs(s / surface - 1) < 1e-12 && shared_rules._weights.empty()
                      && std::abs(ends / (2 * w * w) - 1) < 1e-12,
                  std::string("two boxes sharing a plane have no surface on it, ") + kind);
        }
        // a smaller second box: the plane keeps the ring around it
        const std::vector<LevelSetFunction<double, int>> ls
            = {made_of({plane(0, -0.7, "xlo"), plane(0, 0.1, "p"), plane(0, 0.7, "xhi"), plane(1, -0.6, "ylo"),
                        plane(1, 0.45, "yhi"), plane(2, -0.6, "zlo"), plane(2, 0.45, "zhi"), plane(1, -0.35, "ylo2"),
                        plane(1, 0.3, "yhi2"), plane(2, -0.35, "zlo2"), plane(2, 0.3, "zhi2")},
                       "(" + first + ") or (p > 0 and xhi < 0 and ylo2 > 0 and yhi2 < 0 and zlo2 > 0 and zhi2 < 0)",
                       "union")};
        const CutResult<double, int> r = cut<double, int>(mesh.view, ls);
        const MeshPart<double, int> faces = select(r, "union = 0");
        const double v = total(select(r, "union < 0"), 2, true), s = total(faces, 2, false);
        const double ring = total(on(faces, {"p"}), 2, false);
        const double a2 = 0.65 * 0.65; // the second box's face on the plane
        const double volume = 0.8 * w * w + 0.6 * a2;
        const double surface = 2 * (2 * 0.8 * w + w * w) + 2 * (2 * 0.6 * 0.65 + a2) - 2 * a2;
        std::printf("two boxes %s, the second smaller: volume %.1e, surface %.1e, ring on the plane %.1e\n", kind,
                    v / volume - 1, s / surface - 1, ring / (w * w - a2) - 1);
        check(std::abs(v / volume - 1) < 1e-12 && std::abs(s / surface - 1) < 1e-12
                  && std::abs(ring / (w * w - a2) - 1) < 1e-12,
              std::string("a smaller second box leaves a ring on the shared plane, ") + kind);
    }
}

} // namespace

int main()
{
    test_expansions();
    test_plain_selections();
    test_bad_insides();
    test_cut();
    test_napkin_ring();
    test_boxes();
    return failures == 0 ? 0 : 1;
}
