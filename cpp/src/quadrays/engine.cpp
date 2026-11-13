// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "engine.h"
#include "curves.h"
#include "gauss.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <deque>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cutcells::quadrays
{

Part part_of(const SelectionTerm& term, int level_set)
{
    if (level_set < 0 || level_set >= 64)
        throw std::invalid_argument("quadrays: level-set index out of range");
    const std::uint64_t bit = std::uint64_t(1) << level_set;
    if ((term.zero_required | term.negative_required | term.positive_required) & ~bit)
    {
        throw std::invalid_argument(
            "quadrays: the selection term constrains another level set; "
            "one level set per term is supported");
    }
    const bool zero = term.zero_required & bit, negative = term.negative_required & bit,
               positive = term.positive_required & bit;
    if (int(zero) + int(negative) + int(positive) > 1)
        throw std::invalid_argument("quadrays: the selection term asks for contradictory signs");
    if (zero)
        return Part::interface;
    if (negative)
        return Part::negative;
    if (positive)
        return Part::positive;
    return Part::whole;
}

SelectionTerm term_of(Part part, int level_set)
{
    if (level_set < 0 || level_set >= 64)
        throw std::invalid_argument("quadrays: level-set index out of range");
    const std::uint64_t bit = std::uint64_t(1) << level_set;
    SelectionTerm term;
    switch (part)
    {
    case Part::negative:
        term.negative_required = bit;
        break;
    case Part::positive:
        term.positive_required = bit;
        break;
    case Part::interface:
        term.zero_required = bit;
        break;
    default:
        break;
    }
    return term;
}

int point_limit(const Options& opt, int q, int tdim)
{
    if (opt.max_points != 0)
        return opt.max_points < 0 ? std::numeric_limits<int>::max() : opt.max_points;
    return static_cast<int>(std::min(1000.0 * std::pow(static_cast<double>(q), tdim), 300000.0));
}

namespace
{

template <typename T, int D>
using VecD = std::array<T, D>;

// Tolerances of the prototype in double precision, scaled for float.
template <std::floating_point T>
constexpr T tiny = scaled_tolerance<T>(1e-14);
template <std::floating_point T>
constexpr T vertex_det_tol = scaled_tolerance<T>(1e-13);
template <std::floating_point T>
constexpr T feasible_tol = scaled_tolerance<T>(1e-11);
template <std::floating_point T>
constexpr T extent_tol = scaled_tolerance<T>(1e-13);
template <std::floating_point T>
constexpr T compare_tol = scaled_tolerance<T>(1e-12);
template <std::floating_point T>
constexpr T segment_tol = scaled_tolerance<T>(1e-15);
template <std::floating_point T>
constexpr T zero_function_tol = scaled_tolerance<T>(1e-12);
/// Accuracy checks: differences below this times the cell's measure are rounding.
template <std::floating_point T>
constexpr T accuracy_floor = scaled_tolerance<T>(1e-13);
/// level_margin: the margin of a level with D free coordinates is at most margin_cap / sqrt(D).
template <std::floating_point T>
constexpr T margin_cap = T(0.9);
/// Accuracy checks of boxes split at folds and crossings (Options::fold_splits):
/// the largest relative difference between q and q + 1 points they keep.
template <std::floating_point T>
constexpr T split_tolerance = scaled_tolerance<T>(1e-12);
/// Bisections of a level-2 box for the moment of the accuracy check alone
/// (Options::moment_tolerance).
constexpr int moment_depth = 6;
/// Those bisections in a cell: at a singular point (a cone's apex) more points do not
/// help, and neither do more boxes by much.
constexpr int moment_bisections = 32;
/// Points a level-2 box may add per direction for that moment before it is bisected.
constexpr int moment_raises = 4;
/// Level-2 bisections allowed below a 3D box without a margin along its height
/// direction (Options::silhouette_splits); beyond, the box itself is bisected.
constexpr int unmargined_bisections = 8;
/// A fold anchors the rule of a segment that ends within this many segment
/// lengths of it (anchored_nodes, fold_anchors).
template <std::floating_point T>
constexpr T fold_reach = T(4);
/// Points of a curve part on a 2D cell: a root at an end of a line counts
/// within this Newton distance |f| / |f'| (box units), and points closer than
/// twice it are one.
template <std::floating_point T>
constexpr T point_tol = std::max(static_cast<T>(1e-10), T(100) * std::numeric_limits<T>::epsilon());

template <std::floating_point T>
constexpr T infinity = std::numeric_limits<T>::infinity();

/// Entries of a level map below tiny count as zero, as for the degrees.
template <std::floating_point T>
T snap(T a)
{
    return std::abs(a) > tiny<T> ? a : T(0);
}

/// a and b non-zero with opposite signs: a root lies between them.
template <std::floating_point T>
bool opposite(T a, T b)
{
    return a != T(0) && b != T(0) && (a > T(0)) != (b > T(0));
}

// ============================================================================
// Problems at one level
// ============================================================================

/// psi(y) = phi_ls(A y + b) (curved), or a . y + c (linear), on level coordinates y.
template <typename T, int D>
struct Func
{
    std::uint32_t origin = 1; ///< diagnostics: chain of 4-bit codes, see origin_name
    bool linear = false;
    int ls = 0; ///< the level set of a curved function
    std::array<VecD<T, D>, 3> A{};
    Vec3<T> b{};
    VecD<T, D> a{};
    T c = 0;
};

/// A level set restricted to another's zero set along the height lines of the
/// level above: s(y) = on(y, r(y)), with r(y) the root of t -> under(y, t). It
/// vanishes where the roots of under and on cross. The inner functions are
/// curved functions of the level above, whose coordinate k is t; under is
/// monotone in t on the box's height range [t_lo, t_hi]. The root is followed
/// half that range beyond it, so that s stays smooth where the root leaves the
/// box (it would have a kink there if clamped at the range, which no box
/// along that curve could certify); beyond, it is clamped. Crossings outside
/// [t_lo, t_hi] are not in the box: they only add breakpoints.
template <typename T, int D>
struct Surface
{
    std::uint32_t origin = 11;
    Func<T, D + 1> under, on;
    int k = 0;
    T t_lo = 0, t_hi = 0;
    bool silhouette = false; ///< under is d_t on: it vanishes where the zero set of on turns along t
    /// under has two roots per line, on either side of its extreme point: this
    /// surface follows the one below it (-1) or above it (1). 0: under is monotone.
    int branch = 0;
    int slope_ls = -1;       ///< branch != 0: the level set d_t under
    int curvature = 0;       ///< branch != 0: the sign of d_t d_t under
};

/// A curved function of the level above, of this level's coordinates and t,
/// its coordinate k: its roots on the height lines above, t in [t_lo, t_hi],
/// are the height function of the integrand along them. Where its zero set
/// turns back along t (f = d_t f = 0; a surface's silhouette seen along t),
/// beyond the box, that height function has a square-root branch point
/// (Options::anchor_folds).
template <typename T, int D>
struct Height
{
    Func<T, D + 1> f;
    int k = 0;
    T t_lo = 0, t_hi = 0;
    bool two = false; ///< two roots per height line, else monotone along t on the box
};

/// c . y <= d
template <typename T, int D>
struct Half
{
    VecD<T, D> c{};
    T d = 0;
};

template <typename T, int D>
struct Problem
{
    VecD<T, D> lo{}, hi{};
    std::vector<Func<T, D>> funcs;
    std::vector<Surface<T, D>> surfaces;
    std::vector<Half<T, D>> clips;
    /// Below a certified 3D box: the curved functions whose roots its height
    /// lines take, at whose folds the height lines here are anchored
    /// (Options::anchor_folds).
    std::vector<Height<T, D>> heights;
    /// Top level: the map from the level's coordinates to box coordinates (the
    /// identity, or a rotation of it).
    Func<T, D> frame;
};

/// y_k = alpha + beta . y' on the base coordinates y'
template <typename T, int D>
struct Bound
{
    T alpha = 0;
    VecD<T, D - 1> beta{};
    int kind = 7; ///< diagnostics: 2/3 box face below/above, 4/5 clip below/above, 6 linear root
};

/// Origin codes: 1 phi; restricted to 2 lower box face, 3 upper box face, 4 lower
/// clip, 5 upper clip, 6 linear root, 7 unchanged; new linear functions: 8 switch
/// of lower bounds, 9 switch of upper bounds, 10 difference of linear roots; 11
/// a surface function (one level set on another's zero set).
std::string origin_name(std::uint32_t origin)
{
    static const char* names[] = {"?",        "phi",       "box-lo",    "box-hi",   "clip-lo", "clip-up",
                                  "lin-root", "same",      "switch-lo", "switch-up", "root-diff", "surface",
                                  "corner"};
    std::vector<std::string> parts;
    for (; origin != 0; origin >>= 4)
        parts.push_back(names[std::min<std::uint32_t>(origin & 15u, 12u)]);
    std::string s;
    for (auto it = parts.rbegin(); it != parts.rend(); ++it)
        s += (s.empty() ? "" : "|") + *it;
    return s;
}

/// The sign conditions of a selection term on the cell's level sets.
struct Requirement
{
    std::uint64_t negative = 0, positive = 0, zero = 0;
};

// ============================================================================
// Context: the cell, the options and scratch buffers per level
// ============================================================================

/// A fold the first pass of an accuracy check found near a height line of a
/// level-2 box (Options::anchor_folds): the index of its height function
/// (Problem::heights) and the point (tau, t) in the plane through the line.
template <typename T>
struct ColumnFold
{
    int height = 0;
    VecD<T, 2> y{};
};

/// A level set's form with one box coordinate fixed, kept between lines along one axis:
/// the lines of a slice of a base share that coordinate (line_form).
template <std::floating_point T>
struct PlaneForm
{
    int axis = -1;  ///< the axis of the lines
    int fixed = -1; ///< the coordinate fixed in coeffs
    T value = 0;    ///< and its value
    std::array<T, 3> last{}; ///< the origin of the last line
    std::vector<T> coeffs;   ///< over the other two coordinates, the lower one fastest
};

template <std::floating_point T>
struct Context
{
    std::vector<Source<T>> phis; ///< the level sets on the cell's unit box, then the derived ones
    /// Options::silhouette_splits: the derivatives of the cell's level sets along
    /// the axes of its box, level sets that no term names, whose forms these are.
    std::deque<BoxBernstein<T>> derived_forms;
    std::deque<DirectionalDerivative> derived_directions; ///< of analytic level sets
    std::deque<AnalyticLevelSet> derived_analytic;
    std::vector<int> derived; ///< index in phis of d_a phi_l at 3 l + a, or -1
    std::vector<BoxBernstein<T>> linear_forms; ///< of analytic level sets linear on the cell
    std::vector<T> scales;       ///< size of each (reference_magnitude)
    mutable std::vector<PlaneForm<T>> plane_forms; ///< by level set (line_form)
    bool any_analytic = false;   ///< some level set is analytic
    bool silhouettes = false;    ///< Options::silhouette_splits applies to this run
    int abandoned = 0;           ///< top boxes tried without a margin whose base did not settle, so far
    bool measure_only = false;   ///< the second run of an accuracy check: only the measure counts
    bool moment_only = false;    ///< the next bisection is for the moment alone
    T moment_floor = 0;          ///< > 0: below such a bisection, the tolerance on the moment it failed
    T measure_floor = 0;         ///< and the tolerance on the measure it passed
    int moment_degree = 0;       ///< of the Legendre polynomial in Options::moment_tolerance
    // what those tolerances ask of a rule with q points: all they say from q = 6,
    // less below (tolerance^(q / 6)), which a rule with fewer points cannot reach
    T moment_relative = 0;       ///< Options::moment_tolerance
    Vec3<T> ref_origin{};        ///< box coordinates to the cell's reference coordinates
    Mat3<T> ref_jacobian{};
    T split_relative = 0;        ///< split_tolerance
    std::vector<Requirement> terms;
    std::uint64_t used = 0; ///< level sets some term constrains
    int surface = -1;       ///< interface parts: the level set whose zero set is integrated
    int partner = -1;       ///< curve parts: the second level set, whose zero set meets surface's
    Options opt;
    T margin = 0;
    Stats* stats = nullptr;
    int vis_order = 0; ///< > 0: leaf nodes (vis_order + 1 per segment) instead of Gauss points
    T detj = 0;        ///< |det| of the box-to-physical map
    Mat3<T> inv = {};  ///< its inverse (interface weights)
    Mat3<T> jac = {};  ///< the map itself (curve weights)
    std::array<int, 3> box_ids{};
    int bisections = 0; ///< bisections so far in this cell
    int splits_blocked = 0; ///< > 0: inside a box whose fold-split search did not settle
    CellPoints<T>* emitted = nullptr; ///< the cell's points, as ctx.out but also without accuracy checks
    int first_point = 0;    ///< the first of this run in it
    int points_allowed = 0; ///< Options::max_points of this run
    bool over_limit = false; ///< it exceeded them: the recursion unwinds, and the cell is integrated again
    BoxTrace<T>* trace = nullptr; ///< not null: the boxes are recorded
    int trace_box = -1;     ///< the box being integrated

    // Gauss-Legendre rule on [0, 1], and the one with a point more for accuracy checks
    int q = 0;
    int cell_q = 0;   ///< the q of this run, which a box short of points raises for itself
    int q_raised = 0; ///< points more per direction in use (moment_raises)
    std::vector<T> gauss_x, gauss_w, check_x, check_w;
    // Gauss-Legendre rules by their number of points, for raise_rules
    std::vector<std::vector<T>> rule_x, rule_w;
    // the rules of segments anchored at a fold and at folds beyond both ends
    // (Options::anchor_folds, anchored_nodes), and those of the checks
    std::vector<T> fold_x, fold_w, fold_check_x, fold_check_w, cubic_x, cubic_w, cubic_check_x, cubic_check_w;
    // the same with twice the points, for boxes the accuracy check finds short of them
    std::vector<T> more_fold_x, more_fold_w, more_fold_check_x, more_fold_check_w, more_cubic_x, more_cubic_w,
        more_cubic_check_x, more_cubic_check_w;
    bool more_points = false;          ///< those are in use
    mutable bool anchored_any = false; ///< a segment took the anchored rule since this was last reset
    CellPoints<T>* out = nullptr; ///< the points emitted so far, which accuracy checks roll back
    std::vector<char> flat;       ///< per level set: affine on the cell
    bool curved = false;          ///< the current top box has a curved function that may vanish
    T measure_scale = 0;          ///< the cell's measure in the part's dimension (accuracy checks)
    T cell_size = 0;              ///< |det J|^(1 / tdim) (accuracy checks)
    T kappa = 0;                  ///< curvature at the centre of the current top box (accuracy checks)
    T top_extent = 0;             ///< its extent along its height direction (accuracy checks)
    bool top_unmargined = false;  ///< the current top box has no margin along it (Options::silhouette_splits)
    const std::vector<Half<T, 3>>* top_clips = nullptr; ///< its clip planes, in box coordinates; null in a rotated frame
    int level2_bisections = 0;    ///< bisections below it so far
    int moment_splits = 0;        ///< bisections for the moment alone in this cell (moment_bisections)
    bool unmargined_failed = false; ///< its base does not settle: its points are dropped, and it is bisected
    int tdim = 0;

    // scratch, indexed by level (number of free coordinates)
    std::array<BoxBernstein<T>, 4> form, deriv, second, sub, line;
    std::array<std::vector<T>, 4> form_work, line_work, root_work, nodes, slope_line;
    // scratch of surface functions, indexed by the level of their inner functions
    std::array<BoxBernstein<T>, 4> inner_form, inner_deriv, inner_line;
    std::array<std::vector<T>, 4> inner_work;
    // the level sets' Bernstein forms on the current top-level height line
    std::vector<BoxBernstein<T>> top_lines;
    std::vector<char> top_line_set;
    std::vector<T> top_values;

    // The folds of an accuracy-checked level-2 box (Options::anchor_folds): its
    // own (fold_anchors, bound_folds) and those of its height lines
    // (height_folds), found on the first pass, which the second takes on.
    int fold_pass = 0; ///< 0: search; 1: search and record; 2: take the records on
    std::vector<T> box_folds;
    std::vector<T> column_y;       ///< base coordinate of each recorded height line
    std::vector<int> column_start; ///< its folds: column_fold[column_start[i], column_start[i + 1])
    std::vector<ColumnFold<T>> column_fold;
    std::vector<int> column_order; ///< the recorded lines by base coordinate
    // scratch of the fold search
    Problem<T, 2> fold_plane;
    std::vector<VecD<T, 2>> fold_starts, fold_found;
    std::array<std::vector<T>, 4> line_folds; ///< by level: the folds above a height line
};

/// The margin required on boxes with D free coordinates. Some axis has
/// |n_k| >= 1/sqrt(D) for every unit normal n, but no axis more where n points
/// between them, so larger margins bisect until the limits and leave boxes
/// uncertified: they are capped below 1/sqrt(D), which the bounds' overestimate
/// of |grad psi| leaves room for.
template <std::floating_point T, int D>
T level_margin(const Context<T>& ctx)
{
    return std::min(ctx.margin, margin_cap<T> / std::sqrt(T(D)));
}

template <std::floating_point T>
bool is_analytic(const Context<T>& ctx, int ls)
{
    return ctx.phis[static_cast<std::size_t>(ls)].is_analytic();
}

// ============================================================================
// Evaluation
// ============================================================================

template <typename T, int D>
Vec3<T> to_u(const Func<T, D>& f, const VecD<T, D>& y)
{
    Vec3<T> u = f.b;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < D; ++j)
            u[i] += f.A[i][j] * y[j];
    return u;
}

/// Box coordinates of a top-level point (the identity map, u_j = 0 for j >= D).
template <typename T, int D>
Vec3<T> pad(const VecD<T, D>& y)
{
    Vec3<T> u = {0, 0, 0};
    for (int j = 0; j < D; ++j)
        u[j] = y[j];
    return u;
}

template <std::floating_point T, int D>
T value(const Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& y)
{
    if (f.linear)
    {
        T s = f.c;
        for (int j = 0; j < D; ++j)
            s += f.a[j] * y[j];
        return s;
    }
    const Vec3<T> u = to_u<T, D>(f, y);
    return evaluate(ctx.phis[static_cast<std::size_t>(f.ls)], std::span<const T>(u));
}

template <std::floating_point T, int D>
T derivative_along(const Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& y, int k)
{
    if (f.linear)
        return f.a[k];
    const Vec3<T> u = to_u<T, D>(f, y);
    Vec3<T> g = {0, 0, 0};
    gradient(ctx.phis[static_cast<std::size_t>(f.ls)], std::span<const T>(u), std::span<T>(g));
    T s = T(0);
    for (int i = 0; i < 3; ++i)
        s += g[i] * f.A[i][k];
    return s;
}

/// The gradient of a function in its level's coordinates.
template <std::floating_point T, int D>
void func_gradient(const Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& y, VecD<T, D>& g)
{
    if (f.linear)
    {
        g = f.a;
        return;
    }
    const Vec3<T> u = to_u<T, D>(f, y);
    Vec3<T> gu = {0, 0, 0};
    gradient(ctx.phis[static_cast<std::size_t>(f.ls)], std::span<const T>(u), std::span<T>(gu));
    for (int d = 0; d < D; ++d)
    {
        g[d] = T(0);
        for (int i = 0; i < 3; ++i)
            g[d] += gu[i] * f.A[i][d];
    }
}

template <typename T, int D>
VecD<T, D> insert(const VecD<T, D - 1>& yb, int k, T t)
{
    VecD<T, D> y{};
    for (int j = 0, jb = 0; j < D; ++j)
        y[j] = j == k ? t : yb[jb++];
    return y;
}

template <typename T, int D>
T bound_value(const Bound<T, D>& b, const VecD<T, D - 1>& yb)
{
    T s = b.alpha;
    for (int j = 0; j < D - 1; ++j)
        s += b.beta[j] * yb[j];
    return s;
}

/// f restricted to y_k = alpha + beta . y'
template <typename T, int D>
Func<T, D - 1> restrict_to(const Func<T, D>& f, int k, const Bound<T, D>& b)
{
    Func<T, D - 1> r;
    r.origin = (f.origin << 4) | static_cast<std::uint32_t>(b.kind);
    r.linear = f.linear;
    r.ls = f.ls;
    for (int jb = 0; jb < D - 1; ++jb)
    {
        const int j = jb < k ? jb : jb + 1;
        if (f.linear)
            r.a[jb] = f.a[j] + f.a[k] * b.beta[jb];
        else
            for (int i = 0; i < 3; ++i)
                r.A[i][jb] = f.A[i][j] + f.A[i][k] * b.beta[jb];
    }
    if (f.linear)
        r.c = f.c + f.a[k] * b.alpha;
    else
        for (int i = 0; i < 3; ++i)
            r.b[i] = f.b[i] + f.A[i][k] * b.alpha;
    return r;
}

/// Linear base function b1(y') - b2(y')
template <typename T, int D>
Func<T, D - 1> difference(const Bound<T, D>& b1, const Bound<T, D>& b2)
{
    Func<T, D - 1> r;
    r.origin = b1.kind == 6 ? 10u : (b1.kind == 2 || b1.kind == 4 ? 8u : 9u);
    r.linear = true;
    r.c = b1.alpha - b2.alpha;
    for (int j = 0; j < D - 1; ++j)
        r.a[j] = b1.beta[j] - b2.beta[j];
    return r;
}

/// A function of the level above, whose coordinate kt is t, restricted to
/// y_k = alpha + beta . y' of this level: along the same coordinate among its
/// inner coordinates, which leaves t at @p rk.
template <typename T, int D>
Func<T, D> restrict_inner(const Func<T, D + 1>& f, int kt, int k, const Bound<T, D>& b, int& rk)
{
    // y_k among the inner coordinates insert(y, kt, t), and t after removing it
    const int K = k < kt ? k : k + 1;
    rk = k < kt ? kt - 1 : kt;
    Bound<T, D + 1> inner;
    inner.alpha = b.alpha;
    inner.kind = b.kind;
    for (int jb = 0; jb < D - 1; ++jb)
        inner.beta[jb < rk ? jb : jb + 1] = b.beta[jb];
    inner.beta[rk] = T(0);
    return restrict_to<T, D + 1>(f, K, inner);
}

/// A surface function restricted to y_k = alpha + beta . y': its inner
/// functions restricted along the same coordinate.
template <typename T, int D>
Surface<T, D - 1> restrict_surface(const Surface<T, D>& s, int k, const Bound<T, D>& b)
{
    Surface<T, D - 1> r;
    r.origin = (s.origin << 4) | static_cast<std::uint32_t>(b.kind);
    r.t_lo = s.t_lo;
    r.t_hi = s.t_hi;
    r.silhouette = s.silhouette;
    r.branch = s.branch;
    r.slope_ls = s.slope_ls;
    r.curvature = s.curvature;
    r.under = restrict_inner<T, D>(s.under, s.k, k, b, r.k);
    r.on = restrict_inner<T, D>(s.on, s.k, k, b, r.k);
    return r;
}

/// A height function restricted to y_k = alpha + beta . y', likewise.
template <typename T, int D>
Height<T, D - 1> restrict_height(const Height<T, D>& h, int k, const Bound<T, D>& b)
{
    Height<T, D - 1> r;
    r.t_lo = h.t_lo;
    r.t_hi = h.t_hi;
    r.two = h.two;
    r.f = restrict_inner<T, D>(h.f, h.k, k, b, r.k);
    return r;
}

/// A function of the level above, whose coordinate kt is t, on the plane
/// through the line y = insert(yb, k, tau) of this level and t: a function of
/// (tau, t).
template <typename T, int D>
Func<T, 2> line_plane(const Func<T, D + 1>& f, const VecD<T, D - 1>& yb, int k, int kt)
{
    const int K = k < kt ? k : k + 1; // tau among the inner coordinates
    // the inner coordinates at tau = 0 and t = 0
    const VecD<T, D + 1> y0 = insert<T, D + 1>(insert<T, D>(yb, k, T(0)), kt, T(0));
    Func<T, 2> g;
    g.origin = f.origin;
    g.linear = f.linear;
    g.ls = f.ls;
    if (f.linear)
    {
        g.a = {f.a[K], f.a[kt]};
        g.c = f.c;
        for (int j = 0; j < D + 1; ++j)
            g.c += f.a[j] * y0[j];
    }
    else
    {
        g.b = to_u<T, D + 1>(f, y0);
        for (int i = 0; i < 3; ++i)
            g.A[i] = {f.A[i][K], f.A[i][kt]};
    }
    return g;
}

/// A surface function along the line y = insert(yb, k, tau), as a surface
/// function of tau: inner coordinates (tau, t).
template <typename T, int D>
Surface<T, 1> surface_on_line(const Surface<T, D>& s, const VecD<T, D - 1>& yb, int k)
{
    Surface<T, 1> r;
    r.origin = s.origin;
    r.t_lo = s.t_lo;
    r.t_hi = s.t_hi;
    r.silhouette = s.silhouette;
    r.branch = s.branch;
    r.slope_ls = s.slope_ls;
    r.curvature = s.curvature;
    r.k = 1;
    r.under = line_plane<T, D>(s.under, yb, k, s.k);
    r.on = line_plane<T, D>(s.on, yb, k, s.k);
    return r;
}

// ============================================================================
// Bernstein forms of curved functions
// ============================================================================

/// Bernstein form of a curved function on the box [lo, hi] of its level: exact
/// restriction of its level set to the affine image of the box.
template <std::floating_point T, int D>
void bernstein_form(const Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi,
                    BoxBernstein<T>& out, std::vector<T>& work)
{
    const BoxBernstein<T>& phi = *ctx.phis[static_cast<std::size_t>(f.ls)].bernstein;
    std::array<T, 3> origin;
    std::array<T, 3 * D> matrix;
    for (int i = 0; i < 3; ++i)
    {
        T o = f.b[i];
        for (int j = 0; j < D; ++j)
        {
            const T a = snap(f.A[i][j]);
            o += a * lo[j];
            matrix[i * D + j] = a * (hi[j] - lo[j]);
        }
        origin[i] = o;
    }
    // a form of a 2D cell reads the first two rows
    restrict_affine(phi, std::span<const T>(origin.data(), static_cast<std::size_t>(phi.dim)),
                    std::span<const T>(matrix.data(), static_cast<std::size_t>(phi.dim * D)), D, out, work);
}

/// Bernstein basis of degree n at s.
template <std::floating_point T>
void bernstein_basis(int n, T s, std::array<T, max_box_degree + 1>& b)
{
    const T t = T(1) - s;
    b[0] = T(1);
    for (int r = 1; r <= n; ++r)
    {
        T carry = T(0);
        for (int j = 0; j < r; ++j)
        {
            const T v = b[static_cast<std::size_t>(j)];
            b[static_cast<std::size_t>(j)] = carry + t * v;
            carry = s * v;
        }
        b[static_cast<std::size_t>(r)] = carry;
    }
}

/// The Bernstein form on [0, 1] of a 3D form along a line of direction @p step along
/// @p axis from @p origin. Lines come in slices sharing one coordinate: the contraction
/// along it is kept in @p plane, and only the other is done for each line.
template <std::floating_point T>
void axis_line_form(const BoxBernstein<T>& phi, PlaneForm<T>& plane, const std::array<T, 3>& origin, int axis, T step,
                    BoxBernstein<T>& out)
{
    const int i = axis == 0 ? 1 : 0, j = axis == 2 ? 1 : 2; // the fixed coordinates, i < j
    const std::array<int, 3> extent = {phi.degree[0] + 1, phi.degree[1] + 1, phi.degree[2] + 1};
    const std::array<int, 3> stride = {1, extent[0], extent[0] * extent[1]};
    std::array<T, max_box_degree + 1> basis;
    const auto index = [](int c) { return static_cast<std::size_t>(c); };
    if (!(plane.axis == axis && ((plane.fixed == i && plane.value == origin[index(i)])
                                 || (plane.fixed == j && plane.value == origin[index(j)]))))
    {
        // fix the coordinate that stayed from the last line
        const int fixed = plane.axis == axis && plane.last[index(i)] != origin[index(i)]
                                  && plane.last[index(j)] == origin[index(j)]
                              ? j
                              : i;
        const int low = fixed == 0 ? 1 : 0, high = fixed == 2 ? 1 : 2; // the kept coordinates
        bernstein_basis(phi.degree[index(fixed)], origin[index(fixed)], basis);
        plane.coeffs.resize(index(extent[index(low)]) * index(extent[index(high)]));
        std::size_t at = 0;
        for (int b = 0; b < extent[index(high)]; ++b)
            for (int a = 0; a < extent[index(low)]; ++a)
            {
                const T* source = phi.coeffs.data() + a * stride[index(low)] + b * stride[index(high)];
                T sum = T(0);
                for (int m = 0; m < extent[index(fixed)]; ++m)
                    sum += basis[index(m)] * source[m * stride[index(fixed)]];
                plane.coeffs[at++] = sum;
            }
        plane.axis = axis;
        plane.fixed = fixed;
        plane.value = origin[index(fixed)];
    }
    plane.last = origin;
    // the other fixed coordinate, in the kept plane
    const int other = plane.fixed == i ? j : i;
    const int low = plane.fixed == 0 ? 1 : 0;
    const int n = phi.degree[index(axis)];
    const int plane_low = extent[index(low)];
    const int along = axis == low ? 1 : plane_low, across = other == low ? 1 : plane_low;
    bernstein_basis(phi.degree[index(other)], origin[index(other)], basis);
    std::array<T, max_box_degree + 1> line;
    for (int t = 0; t <= n; ++t)
    {
        const T* source = plane.coeffs.data() + t * along;
        T sum = T(0);
        for (int m = 0; m < extent[index(other)]; ++m)
            sum += basis[index(m)] * source[m * across];
        line[index(t)] = sum;
    }
    out.dim = 1;
    out.degree = {n, 0, 0};
    out.coeffs.resize(index(n + 1));
    reparametrise_1d(std::span<T>(line.data(), index(n + 1)), origin[index(axis)], origin[index(axis)] + step,
                     std::span<T>(out.coeffs));
}

/// Bernstein form on [L, U] of t -> psi(y0 + (t - L) e_k), y0 on the line at t = L.
template <std::floating_point T, int D>
void line_form(const Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& y0, int k, T L, T U,
               BoxBernstein<T>& out, std::vector<T>& work)
{
    const BoxBernstein<T>& phi = *ctx.phis[static_cast<std::size_t>(f.ls)].bernstein;
    std::array<T, 3> origin, column;
    int axis = -1, moving = 0;
    for (int i = 0; i < 3; ++i)
    {
        T o = f.b[i];
        for (int j = 0; j < D; ++j)
            o += snap(f.A[i][j]) * y0[j];
        origin[i] = o;
        column[i] = snap(f.A[i][k]) * (U - L);
        if (column[i] != T(0))
        {
            axis = i;
            ++moving;
        }
    }
    if (phi.dim == 3 && moving == 1)
    {
        if (ctx.plane_forms.size() <= static_cast<std::size_t>(f.ls))
            ctx.plane_forms.resize(static_cast<std::size_t>(f.ls) + 1);
        axis_line_form(phi, ctx.plane_forms[static_cast<std::size_t>(f.ls)], origin, axis,
                       column[static_cast<std::size_t>(axis)], out);
        return;
    }
    restrict_affine(phi, std::span<const T>(origin.data(), static_cast<std::size_t>(phi.dim)),
                    std::span<const T>(column.data(), static_cast<std::size_t>(phi.dim)), 1, out, work);
}

/// The line y = y0 + t e_k of a level as u = u0 + t dir in box coordinates.
template <typename T, int D>
void box_line(const Func<T, D>& f, const VecD<T, D>& y0, int k, Vec3<T>& u0, Vec3<T>& dir)
{
    u0 = to_u<T, D>(f, y0);
    for (int i = 0; i < 3; ++i)
        dir[i] = f.A[i][k];
}

/// The affine image of [lo, hi] for an analytic level set: origin and matrix
/// (3 x D, row-major) of the map s -> A (lo + diag(hi - lo) s) + b.
template <typename T, int D>
void affine_image(const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi, std::array<T, 3>& origin,
                  std::array<T, 3 * D>& matrix)
{
    for (int i = 0; i < 3; ++i)
    {
        T o = f.b[i];
        for (int j = 0; j < D; ++j)
        {
            o += f.A[i][j] * lo[j];
            matrix[i * D + j] = f.A[i][j] * (hi[j] - lo[j]);
        }
        origin[i] = o;
    }
}

/// Bounds of a curved function on [lo, hi] for an analytic level set, from
/// Taylor models of phi on the affine image of the box (the map to physical
/// space is affine, so the models represent it exactly). Sets the sign of the
/// function if it is certain, the direction margins as margins() does from
/// Bernstein coefficients, and the largest |psi|. Without derivative bounds
/// no direction is certified; without any bound nothing is certain.
template <std::floating_point T, int D>
void analytic_bounds(const Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi,
                     int& sign, VecD<T, D>& ratio, T& magnitude, VecD<T, D>& upper)
{
    upper.fill(infinity<T>);
    std::array<T, 3> origin;
    std::array<T, 3 * D> matrix;
    affine_image<T, D>(f, lo, hi, origin, matrix);
    VecD<T, D> lengths;
    for (int j = 0; j < D; ++j)
        lengths[j] = hi[j] - lo[j];
    AffineBounds<T> b;
    if (!affine_bounds(ctx.phis[static_cast<std::size_t>(f.ls)], std::span<const T>(origin),
                       std::span<const T>(matrix), D, b))
    {
        sign = 0;
        ratio.fill(T(0));
        magnitude = infinity<T>;
        return;
    }
    sign = b.sign;
    magnitude = b.magnitude;
    if (!b.has_derivatives)
    {
        ratio.fill(T(0));
        return;
    }
    // derivatives with respect to y_j = lo_j + lengths_j s_j
    VecD<T, D> lower;
    for (int j = 0; j < D; ++j)
    {
        lower[j] = b.lower[j] / lengths[j];
        upper[j] = b.upper[j] / lengths[j];
    }
    const T norm = scaled_norm(std::span<const T>(upper));
    for (int j = 0; j < D; ++j)
        ratio[j] = norm > T(0) ? lower[j] / norm : T(0);
}

/// analytic_bounds with margins from Taylor models over the M^D sub-boxes of
/// [lo, hi] that may meet the clipped region. The sign and the size come from
/// the whole box: a zero set on the boundary between sub-boxes touches each of
/// them with one sign, and would be lost. The margins count the sub-boxes on
/// which psi may vanish, as local_margins does for Bernstein forms: if d_k psi
/// has one strict sign on all of them, psi has at most one root on every
/// height line in the region. Without such a sub-box (the zero set on their
/// boundaries) or a bound, the whole box's margins stand.
template <std::floating_point T, int D>
void analytic_local_bounds(const Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi,
                           const std::vector<Half<T, D>>& clips, int M, int& sign, VecD<T, D>& ratio,
                           T& magnitude, VecD<T, D>& upper_out)
{
    analytic_bounds<T, D>(ctx, f, lo, hi, sign, ratio, magnitude, upper_out);
    if (sign != 0 || !std::isfinite(magnitude))
        return;
    int cells = 1;
    for (int j = 0; j < D; ++j)
        cells *= M;
    bool relevant = false;
    std::array<int, D> dsign{};
    VecD<T, D> lower, upper;
    lower.fill(infinity<T>);
    upper.fill(T(0));
    const Source<T>& phi = ctx.phis[static_cast<std::size_t>(f.ls)];
    for (int n = 0; n < cells; ++n)
    {
        VecD<T, D> ylo, yhi;
        for (int j = 0, m = n; j < D; ++j, m /= M)
        {
            ylo[j] = lo[j] + (hi[j] - lo[j]) * T(m % M) / T(M);
            yhi[j] = lo[j] + (hi[j] - lo[j]) * T(m % M + 1) / T(M);
        }
        bool meets = true;
        for (const Half<T, D>& h : clips)
        {
            T vmin = T(0);
            for (int j = 0; j < D; ++j)
                vmin += h.c[j] * (h.c[j] >= T(0) ? ylo[j] : yhi[j]);
            meets &= vmin <= h.d + compare_tol<T>;
        }
        if (!meets)
            continue;
        std::array<T, 3> origin;
        std::array<T, 3 * D> matrix;
        affine_image<T, D>(f, ylo, yhi, origin, matrix);
        AffineBounds<T> b;
        if (!affine_bounds(phi, std::span<const T>(origin), std::span<const T>(matrix), D, b) || !b.has_derivatives)
            return; // the whole box's margins stand
        if (b.sign != 0)
            continue; // no root inside (at most on its boundary)
        relevant = true;
        for (int j = 0; j < D; ++j)
        {
            const T length = yhi[j] - ylo[j];
            const int s = b.dmin[j] > T(0) ? 1 : (b.dmax[j] < T(0) ? -1 : 0);
            // one strict sign shared by all counted sub-boxes
            if (s == 0 || (dsign[j] != 0 && dsign[j] != s))
                dsign[j] = 2;
            else if (dsign[j] == 0)
                dsign[j] = s;
            lower[j] = std::min(lower[j], s != 0 ? b.lower[j] / length : T(0));
            upper[j] = std::max(upper[j], b.upper[j] / length);
        }
    }
    if (!relevant)
        return; // the zero set on sub-box boundaries: the whole box's margins
    const T norm = scaled_norm(std::span<const T>(upper));
    VecD<T, D> local;
    for (int j = 0; j < D; ++j)
        local[j] = (norm > T(0) && (dsign[j] == 1 || dsign[j] == -1)) ? lower[j] / norm : T(0);
    ratio = local;
    upper_out = upper;
}

/// Margins from M^D sub-cells. Only sub-cells that may meet the clipped region and
/// on which psi may vanish count. If d_k psi has one strict sign on all of them,
/// psi has at most one root on every height line in the region, and ratio[k] is the
/// smallest local min|d_k psi| / max|grad psi|. Returns false if no sub-cell
/// counts: the function does not vanish in the cell and can be dropped.
template <std::floating_point T, int D>
bool local_margins(Context<T>& ctx, const BoxBernstein<T>& p, const VecD<T, D>& lo, const VecD<T, D>& hi,
                   const std::vector<Half<T, D>>& clips, int M, VecD<T, D>& ratio, VecD<T, D>& upper)
{
    upper.fill(T(0));
    VecD<T, D> local, scale;
    std::array<int, D> sign{};
    local.fill(infinity<T>);
    for (int j = 0; j < D; ++j)
        scale[j] = M / (hi[j] - lo[j]);
    bool relevant = false;
    int cells = 1;
    for (int j = 0; j < D; ++j)
        cells *= M;
    BoxBernstein<T>& sub = ctx.sub[D];
    BoxBernstein<T>& dk = ctx.deriv[D];
    for (int n = 0; n < cells; ++n)
    {
        VecD<T, D> a, b, ylo, yhi;
        for (int j = 0, m = n; j < D; ++j, m /= M)
        {
            a[j] = T(m % M) / M;
            b[j] = T(m % M + 1) / M;
            ylo[j] = lo[j] + (hi[j] - lo[j]) * a[j];
            yhi[j] = lo[j] + (hi[j] - lo[j]) * b[j];
        }
        bool meets = true;
        for (const Half<T, D>& h : clips)
        {
            T vmin = T(0);
            for (int j = 0; j < D; ++j)
                vmin += h.c[j] * (h.c[j] >= T(0) ? ylo[j] : yhi[j]);
            meets &= vmin <= h.d + compare_tol<T>;
        }
        if (!meets)
            continue;
        subdivide(p, std::span<const T>(a), std::span<const T>(b), sub, ctx.form_work[D]);
        if (!may_vanish(std::span<const T>(sub.coeffs)))
            continue;
        relevant = true;
        VecD<T, D> dmin{}, dmax{};
        std::array<int, D> s{};
        for (int k = 0; k < D; ++k)
        {
            if (sub.degree[k] == 0)
                continue;
            derivative(sub, k, dk);
            bool pos = true, neg = true;
            T amin = infinity<T>, amax = T(0);
            for (const T v : dk.coeffs)
            {
                pos &= v > T(0);
                neg &= v < T(0);
                amin = std::min(amin, std::abs(v));
                amax = std::max(amax, std::abs(v));
            }
            s[k] = pos ? 1 : (neg ? -1 : 0);
            dmin[k] = (pos || neg) ? amin * scale[k] : T(0);
            dmax[k] = amax * scale[k];
        }
        const T norm = scaled_norm(std::span<const T>(dmax));
        for (int k = 0; k < D; ++k)
            upper[k] = std::max(upper[k], dmax[k]);
        for (int k = 0; k < D; ++k)
        {
            // one strict sign shared by all counted sub-cells
            if (s[k] == 0 || (sign[k] != 0 && sign[k] != s[k]))
                sign[k] = 2;
            else if (sign[k] == 0)
                sign[k] = s[k];
            local[k] = std::min(local[k], norm > T(0) ? dmin[k] / norm : T(0));
        }
    }
    for (int k = 0; k < D; ++k)
        ratio[k] = (relevant && (sign[k] == 1 || sign[k] == -1)) ? local[k] : T(0);
    return relevant;
}

// ============================================================================
// Bounds of curved functions and surface functions as intervals
// ============================================================================

/// Bounds of a curved function over a box: the sign of its value if certain,
/// and intervals of its derivatives with respect to the level's coordinates.
template <typename T, int D>
struct Ranges
{
    bool valid = false; ///< false: no bound holds
    int sign = 0;
    bool has_derivatives = false;
    VecD<T, D> dmin{}, dmax{};
};

template <std::floating_point T, int D>
Ranges<T, D> curved_ranges(Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi,
                           BoxBernstein<T>& form, BoxBernstein<T>& deriv, std::vector<T>& work)
{
    Ranges<T, D> r;
    if (is_analytic(ctx, f.ls))
    {
        std::array<T, 3> origin;
        std::array<T, 3 * D> matrix;
        affine_image<T, D>(f, lo, hi, origin, matrix);
        AffineBounds<T> b;
        if (!affine_bounds(ctx.phis[static_cast<std::size_t>(f.ls)], std::span<const T>(origin),
                           std::span<const T>(matrix), D, b))
            return r;
        r.valid = true;
        r.sign = b.sign;
        r.has_derivatives = b.has_derivatives;
        for (int j = 0; j < D && r.has_derivatives; ++j)
        {
            r.dmin[j] = b.dmin[j] / (hi[j] - lo[j]);
            r.dmax[j] = b.dmax[j] / (hi[j] - lo[j]);
        }
        return r;
    }
    bernstein_form<T, D>(ctx, f, lo, hi, form, work);
    r.valid = true;
    r.sign = coefficient_sign(std::span<const T>(form.coeffs));
    r.has_derivatives = true;
    for (int j = 0; j < D; ++j)
    {
        if (form.degree[j] == 0)
            continue; // constant along y_j
        derivative(form, j, deriv);
        T vmin = infinity<T>, vmax = -infinity<T>;
        for (const T v : deriv.coeffs)
        {
            vmin = std::min(vmin, v);
            vmax = std::max(vmax, v);
        }
        r.dmin[j] = vmin / (hi[j] - lo[j]);
        r.dmax[j] = vmax / (hi[j] - lo[j]);
    }
    return r;
}

template <std::floating_point T>
std::array<T, 2> interval_mul(T a, T b, T c, T d)
{
    const T p[4] = {a * c, a * d, b * c, b * d};
    return {*std::min_element(p, p + 4), *std::max_element(p, p + 4)};
}

/// The range in which a surface function follows the root of under.
template <typename T, int D>
std::array<T, 2> surface_range(const Surface<T, D>& s)
{
    const T e = T(0.5) * (s.t_hi - s.t_lo);
    return {s.t_lo - e, s.t_hi + e};
}

/// under at the inner point Y.
template <std::floating_point T, int D>
T under_value(const Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D + 1>& Y)
{
    return value<T, D + 1>(ctx, s.under, Y);
}

/// The root of t -> under(y, t) in [a, b], where its values @p fa and @p fb
/// have opposite signs.
template <std::floating_point T, int D>
T under_root(Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D>& y, T a, T b, T fa, T fb)
{
    if (is_analytic(ctx, s.under.ls))
    {
        Vec3<T> u0, dir;
        box_line<T, D + 1>(s.under, insert<T, D + 1>(y, s.k, T(0)), s.k, u0, dir);
        return line_root(ctx.phis[static_cast<std::size_t>(s.under.ls)], std::span<const T>(u0),
                         std::span<const T>(dir), a, b, fa, fb);
    }
    BoxBernstein<T>& line = ctx.inner_line[D + 1];
    line_form<T, D + 1>(ctx, s.under, insert<T, D + 1>(y, s.k, a), s.k, a, b, line, ctx.inner_work[D + 1]);
    return ctx.vis_order == 0 ? newton_root(std::span<const T>(line.coeffs), a, b, a, b, fa, fb)
                              : bracketed_root(std::span<const T>(line.coeffs), a, b, a, b, fa, fb);
}

/// The root of a surface function's under on its branch (Surface::branch), within
/// surface_range: between the extreme point of the line and the end of the range on
/// that side. @p clamped: there is none, and this is the end of that interval where
/// |under| is smaller, the extreme point where the line has passed the silhouette.
template <std::floating_point T, int D>
T branch_root(Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D>& y, bool& clamped)
{
    const std::array<T, 2> range = surface_range(s);
    Func<T, D + 1> slope = s.under;
    slope.ls = s.slope_ls;
    const auto slope_at = [&](T t) { return value<T, D + 1>(ctx, slope, insert<T, D + 1>(y, s.k, t)); };
    const T dl = slope_at(range[0]), dh = slope_at(range[1]);
    T c;
    if (opposite(dl, dh))
        c = illinois_root(slope_at, range[0], range[1], dl, dh);
    else if (dl == T(0) || dh == T(0))
        c = dl == T(0) ? range[0] : range[1];
    else
        c = (dl > T(0)) == (s.curvature > 0) ? range[0] : range[1]; // the extreme point is beyond that end
    const T a = s.branch < 0 ? range[0] : c, b = s.branch < 0 ? c : range[1];
    clamped = true;
    if (!(b > a))
        return c;
    const T fa = under_value<T, D>(ctx, s, insert<T, D + 1>(y, s.k, a)),
            fb = under_value<T, D>(ctx, s, insert<T, D + 1>(y, s.k, b));
    if (fa == T(0) || fb == T(0))
    {
        clamped = false;
        return fa == T(0) ? a : b;
    }
    if (opposite(fa, fb))
    {
        clamped = false;
        return under_root<T, D>(ctx, s, y, a, b, fa, fb);
    }
    return std::abs(fa) <= std::abs(fb) ? a : b;
}

/// The root of t -> under(y, t): in [t_lo, t_hi], where under is monotone,
/// else beyond the end where |under| is smaller, within surface_range, else
/// that end of the range. On a branch, branch_root.
template <std::floating_point T, int D>
T surface_root(Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D>& y)
{
    if (s.branch != 0)
    {
        bool clamped = false;
        return branch_root<T, D>(ctx, s, y, clamped);
    }
    const T flo = under_value<T, D>(ctx, s, insert<T, D + 1>(y, s.k, s.t_lo)),
            fhi = under_value<T, D>(ctx, s, insert<T, D + 1>(y, s.k, s.t_hi));
    if (flo == T(0))
        return s.t_lo;
    if (fhi == T(0))
        return s.t_hi;
    if (opposite(flo, fhi))
        return under_root<T, D>(ctx, s, y, s.t_lo, s.t_hi, flo, fhi);
    const std::array<T, 2> range = surface_range(s);
    if (std::abs(flo) <= std::abs(fhi))
    {
        const T fa = under_value<T, D>(ctx, s, insert<T, D + 1>(y, s.k, range[0]));
        return opposite(fa, flo) ? under_root<T, D>(ctx, s, y, range[0], s.t_lo, fa, flo) : range[0];
    }
    const T fb = under_value<T, D>(ctx, s, insert<T, D + 1>(y, s.k, range[1]));
    return opposite(fhi, fb) ? under_root<T, D>(ctx, s, y, s.t_hi, range[1], fhi, fb) : range[1];
}

/// The value of on at the root of under: the surface function at y.
template <std::floating_point T, int D>
T surface_value(Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D>& y)
{
    return value<T, D + 1>(ctx, s.on, insert<T, D + 1>(y, s.k, surface_root<T, D>(ctx, s, y)));
}

/// The gradient of a surface function at y: d_j on - d_t on d_j under / d_t
/// under at the root, or d_j on where the root is clamped.
template <std::floating_point T, int D>
void surface_gradient(Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D>& y, VecD<T, D>& g)
{
    bool clamped = false;
    const T r = s.branch != 0 ? branch_root<T, D>(ctx, s, y, clamped) : surface_root<T, D>(ctx, s, y);
    const VecD<T, D + 1> yr = insert<T, D + 1>(y, s.k, r);
    const std::array<T, 2> range = surface_range(s);
    const T ut = derivative_along<T, D + 1>(ctx, s.under, yr, s.k);
    const T ot = derivative_along<T, D + 1>(ctx, s.on, yr, s.k);
    const bool interior = !clamped && r > range[0] && r < range[1] && ut != T(0);
    for (int j = 0; j < D; ++j)
    {
        const int J = j < s.k ? j : j + 1;
        g[j] = derivative_along<T, D + 1>(ctx, s.on, yr, J);
        if (interior)
            g[j] -= ot * derivative_along<T, D + 1>(ctx, s.under, yr, J) / ut;
    }
}

/// What bounds say about a surface function on a box.
enum class SurfaceCase
{
    none,     ///< it has no root in the box (its on function has one sign), or
              ///< under has none, so its zeros are those of on on a box face
    bounded,  ///< it may vanish; dmin, dmax bound its derivatives
    unbounded ///< it may vanish; no derivative bound
};

/// Bounds of a surface function on [lo, hi] from those of its inner
/// functions. The root r of under lies, over the box, within r(centre) plus
/// its slope bound |d_j r| <= max |d_j under| / min |d_t under| (bounds on
/// [lo, hi] x [t_lo, t_hi], or on a slab about r(centre) that holds the root
/// where d_t under changes sign on that range) times the half widths; the
/// inner functions are then bounded on that slab, which bisection of the base
/// thins. Where r is
/// interior, d_j s = d_j on - d_t on d_j under / d_t under; where it is
/// clamped, d_j s = d_j on. Intervals of both enclose every one-sided
/// derivative. @p centre_value receives s at the centre.
template <std::floating_point T, int D>
SurfaceCase surface_bounds(Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D>& lo, const VecD<T, D>& hi,
                           VecD<T, D>& dmin, VecD<T, D>& dmax, T& centre_value)
{
    VecD<T, D> centre;
    for (int j = 0; j < D; ++j)
        centre[j] = T(0.5) * (lo[j] + hi[j]);
    const T rc = surface_root<T, D>(ctx, s, centre);
    centre_value = value<T, D + 1>(ctx, s.on, insert<T, D + 1>(centre, s.k, rc));

    // the slab of the root over the box
    const std::array<T, 2> range = surface_range(s);
    VecD<T, D + 1> ilo = insert<T, D + 1>(lo, s.k, range[0]), ihi = insert<T, D + 1>(hi, s.k, range[1]);
    const auto under_ranges = [&](const VecD<T, D + 1>& a, const VecD<T, D + 1>& b)
    {
        return curved_ranges<T, D + 1>(ctx, s.under, a, b, ctx.inner_form[D + 1], ctx.inner_deriv[D + 1],
                                       ctx.inner_work[D + 1]);
    };
    // A branch's root leaves any slab near the silhouette, where its slope is
    // unbounded. Its crossings need both functions to vanish in a slice of the
    // box's height range, which holds there too.
    if (s.branch != 0)
    {
        constexpr int pieces = 8;
        bool crossing = false;
        VecD<T, D + 1> slo = ilo, shi = ihi;
        for (int i = 0; i < pieces && !crossing; ++i)
        {
            slo[s.k] = s.t_lo + (s.t_hi - s.t_lo) * T(i) / T(pieces);
            shi[s.k] = s.t_lo + (s.t_hi - s.t_lo) * T(i + 1) / T(pieces);
            const Ranges<T, D + 1> u = under_ranges(slo, shi);
            if (u.valid && u.sign != 0)
                continue;
            const Ranges<T, D + 1> o = curved_ranges<T, D + 1>(ctx, s.on, slo, shi, ctx.inner_form[D + 1],
                                                               ctx.inner_deriv[D + 1], ctx.inner_work[D + 1]);
            crossing = !o.valid || o.sign == 0;
        }
        if (!crossing)
            return SurfaceCase::none;
    }
    // how far the root strays from rc over the box, from bounds of under on a
    // slab where d_t under has one sign: |d_j r| <= max |d_j under| / min |d_t under|
    const auto slab_reach = [&](const Ranges<T, D + 1>& r, T& reach)
    {
        const T wl = r.dmin[s.k], wh = r.dmax[s.k];
        if (!(wl > T(0) || wh < T(0)))
            return false;
        const T slope_min = std::min(std::abs(wl), std::abs(wh));
        reach = T(0);
        for (int j = 0; j < D; ++j)
        {
            const int J = j < s.k ? j : j + 1;
            reach += std::max(std::abs(r.dmin[J]), std::abs(r.dmax[J])) / slope_min * T(0.5) * (hi[j] - lo[j]);
        }
        return true;
    };
    const Ranges<T, D + 1> wide = under_ranges(ilo, ihi);
    if (wide.valid && wide.sign != 0)
        return SurfaceCase::none; // no root: zeros of on beyond the box
    if (!wide.valid)
        return SurfaceCase::unbounded;
    T reach = T(0);
    if (!wide.has_derivatives || !slab_reach(wide, reach))
    {
        // d_t under changes sign on the range (the box above was certified on
        // its parts where under may vanish only, Taylor sub-boxes), or its
        // derivatives have no bounds there (a distance function, whose Taylor
        // models fail near its centre): a slab about rc on which d_t under has
        // one sign holds the root over the box if the reach its bounds give is
        // at most its half width. Grow it from a thin one.
        if (!(rc > range[0] && rc < range[1]))
        {
            // no root at the centre: none if under keeps one sign on thinner
            // slabs of the range, whose bounds are tighter
            const int pieces = 8;
            int sign = 0;
            for (int i = 0; i < pieces; ++i)
            {
                ilo[s.k] = range[0] + (range[1] - range[0]) * T(i) / T(pieces);
                ihi[s.k] = range[0] + (range[1] - range[0]) * T(i + 1) / T(pieces);
                const Ranges<T, D + 1> piece = under_ranges(ilo, ihi);
                if (!piece.valid || piece.sign == 0 || (sign != 0 && piece.sign != sign))
                    return SurfaceCase::unbounded;
                sign = piece.sign;
            }
            return SurfaceCase::none;
        }
        T half = std::max(T(1e-12), tiny<T>) * (s.t_hi - s.t_lo);
        bool held = false;
        for (int it = 0; it < 6 && !held; ++it)
        {
            ilo[s.k] = std::max(range[0], rc - half);
            ihi[s.k] = std::min(range[1], rc + half);
            const Ranges<T, D + 1> local = under_ranges(ilo, ihi);
            T r = T(0);
            if (!local.valid || !local.has_derivatives || !slab_reach(local, r))
                return SurfaceCase::unbounded;
            held = r <= half;
            if (held)
                reach = r;
            else
                half = T(2) * r;
        }
        if (!held)
            return SurfaceCase::unbounded;
    }
    // a root outside the box's range everywhere: its crossings are not in the box
    if (!(rc - reach <= s.t_hi && rc + reach >= s.t_lo))
        return SurfaceCase::none;
    const T rlo = std::max(range[0], rc - reach), rhi = std::min(range[1], rc + reach);
    // the root stays inside the range it is followed in: no clamping
    const bool interior = rc - reach > range[0] && rc + reach < range[1];
    ilo[s.k] = rlo;
    ihi[s.k] = std::max(rhi, rlo + std::max(T(1e-12), tiny<T>) * (s.t_hi - s.t_lo));

    const Ranges<T, D + 1> on = curved_ranges<T, D + 1>(ctx, s.on, ilo, ihi, ctx.inner_form[D + 1],
                                                        ctx.inner_deriv[D + 1], ctx.inner_work[D + 1]);
    if (on.valid && on.sign != 0)
        return SurfaceCase::none;
    const Ranges<T, D + 1> under = under_ranges(ilo, ihi);
    if (!on.valid || !under.valid || !on.has_derivatives || !under.has_derivatives)
        return SurfaceCase::unbounded;
    const T ul = under.dmin[s.k], uh = under.dmax[s.k];
    if (!(ul > T(0) || uh < T(0)))
        return SurfaceCase::unbounded;
    for (int j = 0; j < D; ++j)
    {
        const int J = j < s.k ? j : j + 1;
        // d_j under / d_t under, times d_t on
        const std::array<T, 2> q = interval_mul(under.dmin[J], under.dmax[J], T(1) / uh, T(1) / ul);
        const std::array<T, 2> p = interval_mul(on.dmin[s.k], on.dmax[s.k], q[0], q[1]);
        dmin[j] = on.dmin[J] - p[1];
        dmax[j] = on.dmax[J] - p[0];
        if (!interior)
        {
            dmin[j] = std::min(dmin[j], on.dmin[J]);
            dmax[j] = std::max(dmax[j], on.dmax[J]);
        }
    }
    return SurfaceCase::bounded;
}

/// Margins of a surface function on [lo, hi] as margins() gives them from
/// derivative bounds; false if it has no root there. With derivative bounds,
/// its value at the centre and the mean-value bound may show one sign.
template <std::floating_point T, int D>
bool surface_margins(Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D>& lo, const VecD<T, D>& hi,
                     VecD<T, D>& ratio)
{
    ratio.fill(T(0));
    VecD<T, D> dmin{}, dmax{};
    T centre_value = T(0);
    const SurfaceCase sc = surface_bounds<T, D>(ctx, s, lo, hi, dmin, dmax, centre_value);
    if (sc == SurfaceCase::none)
        return false;
    if (sc == SurfaceCase::unbounded)
        return true;
    VecD<T, D> upper;
    T reach = T(0);
    for (int j = 0; j < D; ++j)
    {
        upper[j] = std::max(std::abs(dmin[j]), std::abs(dmax[j]));
        reach += upper[j] * T(0.5) * (hi[j] - lo[j]);
    }
    if (std::abs(centre_value) > reach)
        return false;
    const T norm = scaled_norm(std::span<const T>(upper));
    for (int j = 0; j < D; ++j)
    {
        const bool strict = dmin[j] > T(0) || dmax[j] < T(0);
        ratio[j] = (strict && norm > T(0)) ? std::min(std::abs(dmin[j]), std::abs(dmax[j])) / norm : T(0);
    }
    return true;
}

/// True if a surface function's under level set has one sign on [lo, hi]
/// times the height range [t_lo, t_hi] of the box above: its root, and a curve
/// on its zero set, do not cross that box over [lo, hi]. (The surface function
/// follows the root beyond the box, where its bounds may never settle.) Where
/// the whole range shows no sign, thinner slabs of it, whose bounds are
/// tighter (Taylor models of an analytic level set), may.
template <std::floating_point T, int D>
bool under_signed(Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D>& lo, const VecD<T, D>& hi)
{
    VecD<T, D + 1> ilo = insert<T, D + 1>(lo, s.k, s.t_lo), ihi = insert<T, D + 1>(hi, s.k, s.t_hi);
    const auto sign = [&]()
    {
        const Ranges<T, D + 1> r = curved_ranges<T, D + 1>(ctx, s.under, ilo, ihi, ctx.inner_form[D + 1],
                                                           ctx.inner_deriv[D + 1], ctx.inner_work[D + 1]);
        return r.valid ? r.sign : 0;
    };
    const int whole = sign();
    if (whole != 0 || !is_analytic(ctx, s.under.ls))
        return whole != 0;
    const int pieces = 8;
    int first = 0;
    for (int i = 0; i < pieces; ++i)
    {
        ilo[s.k] = s.t_lo + (s.t_hi - s.t_lo) * T(i) / T(pieces);
        ihi[s.k] = s.t_lo + (s.t_hi - s.t_lo) * T(i + 1) / T(pieces);
        const int piece = sign();
        if (piece == 0 || (first != 0 && piece != first))
            return false;
        first = piece;
    }
    return true;
}

/// True if a surface function is 0 up to rounding, relative to the size of its
/// on level set, at the centre and the corners of [lo, hi]: the two zero sets
/// coincide there, and its roots would be rounding noise.
template <std::floating_point T, int D>
bool surface_vanishes_identically(Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D>& lo,
                                  const VecD<T, D>& hi)
{
    const T tol = zero_function_tol<T> * ctx.scales[static_cast<std::size_t>(s.on.ls)];
    VecD<T, D> y;
    for (int j = 0; j < D; ++j)
        y[j] = T(0.5) * (lo[j] + hi[j]);
    if (!(std::abs(surface_value<T, D>(ctx, s, y)) <= tol))
        return false;
    for (int corner = 0; corner < (1 << D); ++corner)
    {
        for (int j = 0; j < D; ++j)
            y[j] = ((corner >> j) & 1) ? hi[j] : lo[j];
        if (!(std::abs(surface_value<T, D>(ctx, s, y)) <= tol))
            return false;
    }
    return true;
}

/// Whether a point of a zero curve in the base of a 3D box may bound the part: some term
/// of the selection holds at it, or nearly (within @p margin times the size of each level
/// set), and on an interface part it lies on the interface. Elsewhere the curve is
/// outside the part: what happens to it there leaves the integrand smooth.
template <std::floating_point T>
bool may_bound_part(Context<T>& ctx, const Vec3<T>& u, T margin)
{
    const auto value_at = [&](int l) { return evaluate(ctx.phis[static_cast<std::size_t>(l)], std::span<const T>(u)); };
    const auto reach = [&](int l) { return margin * ctx.scales[static_cast<std::size_t>(l)]; };
    if (ctx.tdim == 3 && ctx.top_clips != nullptr)
        for (const Half<T, 3>& h : *ctx.top_clips)
            if (h.c[0] * u[0] + h.c[1] * u[1] + h.c[2] * u[2] > h.d + margin)
                return false;
    if (ctx.surface >= 0 && std::abs(value_at(ctx.surface)) > reach(ctx.surface))
        return false;
    for (const Requirement& r : ctx.terms)
    {
        bool ok = true;
        for (std::uint64_t bits = r.negative | r.positive; bits != 0 && ok; bits &= bits - 1)
        {
            const int l = std::countr_zero(bits);
            const T v = value_at(l);
            ok = (r.negative & (std::uint64_t(1) << l)) ? v <= reach(l) : v >= -reach(l);
        }
        if (ok)
            return true;
    }
    return false;
}

/// The 3D point of a surface function's curve above y.
template <std::floating_point T>
Vec3<T> surface_point(Context<T>& ctx, const Surface<T, 2>& s, const VecD<T, 2>& y)
{
    return to_u<T, 3>(s.on, insert<T, 3>(y, s.k, surface_root<T, 2>(ctx, s, y)));
}

/// Whether a zero of a surface function at y is a point of the curve it stands for:
/// under has a root on the line there, in the height range of the box or just beyond
/// it, not the continuation by a clamped root.
template <std::floating_point T, int D>
bool surface_point_in_box(Context<T>& ctx, const Surface<T, D>& s, const VecD<T, D>& y)
{
    bool clamped = false;
    const T r = s.branch != 0 ? branch_root<T, D>(ctx, s, y, clamped) : surface_root<T, D>(ctx, s, y);
    const T height = s.t_hi - s.t_lo;
    if (clamped || r < s.t_lo - T(0.05) * height || r > s.t_hi + T(0.05) * height)
        return false;
    const VecD<T, D + 1> Y = insert<T, D + 1>(y, s.k, r);
    return !(std::abs(under_value<T, D>(ctx, s, Y))
             > T(1e-6) * std::abs(derivative_along<T, D + 1>(ctx, s.under, Y, s.k)) * height);
}

/// Roots of a surface function of one coordinate on (a, b), appended to
/// @p roots: none where bounds show one sign, one where they show a monotone
/// function whose ends differ in sign, bisection otherwise (down to 2^-24 of
/// the interval, at most @p budget intervals).
template <std::floating_point T>
void surface_line_roots(Context<T>& ctx, const Surface<T, 1>& s, T a, T b, T sa, T sb, int depth, int& budget,
                        std::vector<T>& roots)
{
    --budget;
    VecD<T, 1> dmin{}, dmax{};
    T centre_value = T(0);
    const SurfaceCase sc = surface_bounds<T, 1>(ctx, s, VecD<T, 1>{a}, VecD<T, 1>{b}, dmin, dmax, centre_value);
    if (sc == SurfaceCase::none)
        return;
    if (sc == SurfaceCase::bounded
        && std::abs(centre_value) > std::max(std::abs(dmin[0]), std::abs(dmax[0])) * T(0.5) * (b - a))
        return; // one sign on (a, b)
    const bool monotone = sc == SurfaceCase::bounded && (dmin[0] > T(0) || dmax[0] < T(0));
    if (monotone || depth >= 24 || budget <= 0)
    {
        if (opposite(sa, sb))
            roots.push_back(illinois_root([&](T t) { return surface_value<T, 1>(ctx, s, VecD<T, 1>{t}); }, a, b, sa,
                                          sb));
        return;
    }
    const T m = T(0.5) * (a + b);
    const T sm = surface_value<T, 1>(ctx, s, VecD<T, 1>{m});
    if (sm == T(0))
        roots.push_back(m);
    surface_line_roots<T>(ctx, s, a, m, sa, sm, depth + 1, budget, roots);
    surface_line_roots<T>(ctx, s, m, b, sm, sb, depth + 1, budget, roots);
}

/// A surface function is sampled along the line: bounds of its root prove little, and
/// cost far more than the roots do (restricting forms to a box for every line of a
/// base, most of the time of a cell). The boxes below it are checked, which has what
/// this misses. Curve parts and leaf cells, which have no check, keep the bounds.
template <std::floating_point T>
void surface_line_roots(Context<T>& ctx, const Surface<T, 1>& s, T a, T b, std::vector<T>& roots)
{
    const T sa = surface_value<T, 1>(ctx, s, VecD<T, 1>{a}), sb = surface_value<T, 1>(ctx, s, VecD<T, 1>{b});
    if (s.branch != 0 || (ctx.partner < 0 && ctx.vis_order == 0))
    {
        constexpr int pieces = 16;
        const auto at = [&](T t) { return surface_value<T, 1>(ctx, s, VecD<T, 1>{t}); };
        const auto root_between = [&](T t0, T t1, T s0, T s1)
        {
            const T r = illinois_root(at, t0, t1, s0, s1);
            if (surface_point_in_box<T, 1>(ctx, s, VecD<T, 1>{r}))
                roots.push_back(r);
        };
        std::array<T, pieces + 1> ts, vs;
        for (int i = 0; i <= pieces; ++i)
        {
            ts[static_cast<std::size_t>(i)] = i == pieces ? b : a + (b - a) * T(i) / T(pieces);
            vs[static_cast<std::size_t>(i)] = i == 0 ? sa : (i == pieces ? sb : at(ts[static_cast<std::size_t>(i)]));
        }
        for (std::size_t i = 0; i < pieces; ++i)
        {
            if (opposite(vs[i], vs[i + 1]))
                root_between(ts[i], ts[i + 1], vs[i], vs[i + 1]);
            else if (vs[i + 1] == T(0) && i + 1 < pieces)
                roots.push_back(ts[i + 1]);
        }
        // Two roots between samples: |s| has a small minimum there, which takes the
        // opposite sign. A few parabolas through three samples find it.
        for (std::size_t i = 1; i < pieces; ++i)
        {
            T t0 = ts[i - 1], t1 = ts[i], t2 = ts[i + 1], s0 = vs[i - 1], s1 = vs[i], s2 = vs[i + 1];
            if (s1 == T(0) || opposite(s0, s1) || opposite(s1, s2) || std::abs(s1) > std::abs(s0)
                || std::abs(s1) > std::abs(s2))
                continue;
            for (int it = 0; it < 4; ++it)
            {
                const T below = (t1 - t0) * (s1 - s2) - (t1 - t2) * (s1 - s0);
                if (!(std::abs(below) > T(0)))
                    break;
                const T t = t1 - T(0.5) * ((t1 - t0) * (t1 - t0) * (s1 - s2) - (t1 - t2) * (t1 - t2) * (s1 - s0)) / below;
                if (!(t > t0 && t < t2) || t == t1)
                    break;
                const T v = at(t);
                if (opposite(v, s1))
                {
                    root_between(t0, t, s0, v);
                    root_between(t, t2, v, s2);
                    break;
                }
                if (std::abs(v) <= std::abs(s1))
                {
                    (t < t1 ? t2 : t0) = t1;
                    (t < t1 ? s2 : s0) = s1;
                    t1 = t;
                    s1 = v;
                }
                else
                {
                    (t < t1 ? t0 : t2) = t;
                    (t < t1 ? s0 : s2) = v;
                }
            }
        }
        return;
    }
    int budget = 256;
    surface_line_roots<T>(ctx, s, a, b, sa, sb, 0, budget, roots);
}

// ============================================================================
// Clipped boxes at one level
// ============================================================================

template <typename T, int D>
T determinant(const std::array<VecD<T, D>, D>& m)
{
    if constexpr (D == 1)
        return m[0][0];
    else if constexpr (D == 2)
        return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    else
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
               - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
               + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

/// Vertices of {y in [lo, hi] : clips} by enumeration of plane D-tuples.
template <std::floating_point T, int D>
void polytope_vertices(const VecD<T, D>& lo, const VecD<T, D>& hi, const std::vector<Half<T, D>>& clips,
                       std::vector<VecD<T, D>>& vertices)
{
    vertices.clear();
    std::vector<Half<T, D>> planes;
    planes.reserve(2 * D + clips.size());
    for (int j = 0; j < D; ++j)
    {
        Half<T, D> l, u;
        l.c[j] = T(-1);
        l.d = -lo[j];
        u.c[j] = T(1);
        u.d = hi[j];
        planes.push_back(l);
        planes.push_back(u);
    }
    planes.insert(planes.end(), clips.begin(), clips.end());
    const int np = static_cast<int>(planes.size());

    auto vertex = [&](const std::array<int, D>& idx)
    {
        std::array<VecD<T, D>, D> m;
        for (int r = 0; r < D; ++r)
            m[r] = planes[idx[r]].c;
        const T det = determinant<T, D>(m);
        if (std::abs(det) < vertex_det_tol<T>)
            return;
        VecD<T, D> x;
        for (int col = 0; col < D; ++col)
        {
            auto mc = m;
            for (int r = 0; r < D; ++r)
                mc[r][col] = planes[idx[r]].d;
            x[col] = determinant<T, D>(mc) / det;
        }
        for (const Half<T, D>& h : planes)
        {
            T s = T(0);
            for (int j = 0; j < D; ++j)
                s += h.c[j] * x[j];
            if (s > h.d + feasible_tol<T>)
                return;
        }
        vertices.push_back(x);
    };
    if constexpr (D == 1)
    {
        for (int a = 0; a < np; ++a)
            vertex({a});
    }
    else if constexpr (D == 2)
    {
        for (int a = 0; a < np; ++a)
            for (int b = a + 1; b < np; ++b)
                vertex({a, b});
    }
    else
    {
        for (int a = 0; a < np; ++a)
            for (int b = a + 1; b < np; ++b)
                for (int c = b + 1; c < np; ++c)
                    vertex({a, b, c});
    }
}

/// Shrink [lo, hi] to the bounding box of {y in [lo, hi] : clips}; false if empty.
template <std::floating_point T, int D>
bool tighten(VecD<T, D>& lo, VecD<T, D>& hi, const std::vector<Half<T, D>>& clips)
{
    if (clips.empty())
        return true;
    std::vector<VecD<T, D>> vertices;
    polytope_vertices<T, D>(lo, hi, clips, vertices);
    if (vertices.empty())
        return false;
    VecD<T, D> nlo, nhi;
    nlo.fill(infinity<T>);
    nhi.fill(-infinity<T>);
    for (const VecD<T, D>& x : vertices)
        for (int j = 0; j < D; ++j)
        {
            nlo[j] = std::min(nlo[j], x[j]);
            nhi[j] = std::max(nhi[j], x[j]);
        }
    for (int j = 0; j < D; ++j)
    {
        lo[j] = std::max(lo[j], nlo[j]);
        hi[j] = std::min(hi[j], nhi[j]);
        if (!(hi[j] - lo[j] > extent_tol<T>))
            return false;
    }
    return true;
}

/// True if a clip plane cuts the box [lo, hi].
template <std::floating_point T, int D>
bool box_clipped(const VecD<T, D>& lo, const VecD<T, D>& hi, const std::vector<Half<T, D>>& clips)
{
    for (const Half<T, D>& h : clips)
    {
        T vmax = T(0);
        for (int j = 0; j < D; ++j)
            vmax += h.c[j] * (h.c[j] >= T(0) ? hi[j] : lo[j]);
        if (vmax > h.d + compare_tol<T>)
            return true;
    }
    return false;
}

template <std::floating_point T, int D>
bool linear_may_vanish(const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi)
{
    T vmin = f.c, vmax = f.c;
    for (int j = 0; j < D; ++j)
    {
        vmin += std::min(f.a[j] * lo[j], f.a[j] * hi[j]);
        vmax += std::max(f.a[j] * lo[j], f.a[j] * hi[j]);
    }
    return vmin <= T(0) && vmax >= T(0);
}

// ============================================================================
// Analysis of one box: which functions may vanish, and in which direction
// ============================================================================

template <typename T, int D>
struct Analysis
{
    std::vector<Func<T, D>> funcs;
    std::vector<Surface<T, D>> surfaces;
    std::vector<VecD<T, D>> ratios;         ///< margins of each curved function in curved
    std::vector<VecD<T, D>> uppers;         ///< upper bounds of its |d_j psi|
    std::vector<VecD<T, D>> surface_ratios; ///< margins of each surface function
    std::vector<int> curved;                ///< index in funcs of each entry of ratios
    /// Certified without a margin along k: the lines take all roots, and the base
    /// has the silhouettes of the two-root functions (Options::silhouette_splits).
    bool unmargined = false;
    std::vector<char> two_roots;            ///< per entry of funcs: two roots per line along k
    std::vector<char> independent;          ///< per entry of funcs: constant along k
    /// per entry of funcs: top level of an interface part, another level set
    /// than the interface's; only evaluated at the interface's roots, so it
    /// needs no margin along k
    std::vector<char> passive;
    int k = 0;
    T best = -1;
    bool certified = false;
    /// Top level: each level set's sign on the box (0: may vanish or unknown),
    /// the terms left once those signs are known, and whether they decide the
    /// box (+1 all of it, -1 none of it, 0 the level sets decide).
    std::vector<int> signs;
    std::vector<Requirement> terms;
    int decided = 0;
    /// Top level: the level set whose roots the height lines take, the
    /// interface's; in a curve part either of its two level sets, the second
    /// where only it has a margin.
    int surface = -1;
};

/// A curve part's surface function: one of its two level sets on the other's
/// zero set, which vanishes where the curve crosses the height lines.
template <std::floating_point T, int D>
bool partner_surface(const Context<T>& ctx, const Surface<T, D>& s)
{
    return ctx.partner >= 0
           && ((s.under.ls == ctx.surface && s.on.ls == ctx.partner)
               || (s.under.ls == ctx.partner && s.on.ls == ctx.surface));
}

/// The terms that the level sets' signs on a box leave: terms contradicted by
/// a known sign go, conditions met by one are dropped.
inline int residual_terms(const std::vector<Requirement>& terms, const std::vector<int>& signs, int surface,
                          bool surface_vanishes, std::vector<Requirement>& out)
{
    out.clear();
    if (surface >= 0 && !surface_vanishes)
        return -1; // the zero set does not cross the box
    for (const Requirement& t : terms)
    {
        Requirement r;
        bool holds = true;
        for (std::uint64_t bits = t.negative | t.positive; bits != 0 && holds; bits &= bits - 1)
        {
            const int l = std::countr_zero(bits);
            const std::uint64_t bit = std::uint64_t(1) << l;
            const int s = signs[static_cast<std::size_t>(l)];
            const int wanted = (t.negative & bit) ? -1 : 1;
            if (s == 0)
                (wanted < 0 ? r.negative : r.positive) |= bit;
            else
                holds = s == wanted;
        }
        if (!holds)
            continue;
        r.zero = t.zero;
        if ((r.negative | r.positive) == 0)
        {
            out.assign(1, r);
            return 1;
        }
        out.push_back(r);
    }
    return out.empty() ? -1 : 0;
}

/// Second-derivative margin of a curved function along k on [lo, hi]:
/// min |d_k d_k psi| / max |grad d_k psi|, from Bernstein coefficients or
/// Hessian bounds; 0 if d_k d_k psi may vanish or nothing bounds it.
template <std::floating_point T, int D>
T second_margin(Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi, int k,
                int& curvature)
{
    curvature = 0;
    VecD<T, D> lower{}, upper{};
    if (is_analytic(ctx, f.ls))
    {
        std::array<T, 3> origin;
        std::array<T, 3 * D> matrix;
        affine_image<T, D>(f, lo, hi, origin, matrix);
        std::array<T, 9> hlo, hhi;
        if (!affine_hessian(ctx.phis[static_cast<std::size_t>(f.ls)], std::span<const T>(origin),
                            std::span<const T>(matrix), D, hlo, hhi))
            return T(0);
        for (int j = 0; j < D; ++j)
        {
            const T scale = T(1) / ((hi[j] - lo[j]) * (hi[k] - lo[k]));
            const T a = hlo[k * D + j] * scale, b = hhi[k * D + j] * scale;
            upper[j] = std::max(std::abs(a), std::abs(b));
            if (j == k)
            {
                curvature = a > T(0) ? 1 : (b < T(0) ? -1 : 0);
                lower[j] = curvature != 0 ? std::min(std::abs(a), std::abs(b)) : T(0);
            }
        }
    }
    else
    {
        bernstein_form<T, D>(ctx, f, lo, hi, ctx.form[D], ctx.form_work[D]);
        if (ctx.form[D].degree[k] < 2)
            return T(0);
        derivative(ctx.form[D], k, ctx.deriv[D]);
        for (int j = 0; j < D; ++j)
        {
            if (ctx.deriv[D].degree[j] == 0)
                continue;
            derivative(ctx.deriv[D], j, ctx.second[D]);
            bool pos = true, neg = true;
            T amin = infinity<T>, amax = T(0);
            for (const T v : ctx.second[D].coeffs)
            {
                pos &= v > T(0);
                neg &= v < T(0);
                amin = std::min(amin, std::abs(v));
                amax = std::max(amax, std::abs(v));
            }
            const T scale = T(1) / ((hi[j] - lo[j]) * (hi[k] - lo[k]));
            upper[j] = amax * scale;
            if (j == k)
            {
                curvature = pos ? 1 : (neg ? -1 : 0);
                lower[j] = curvature != 0 ? amin * scale : T(0);
            }
        }
    }
    const T norm = scaled_norm(std::span<const T>(upper));
    return norm > T(0) && curvature != 0 ? lower[k] / norm : T(0);
}

/// The sign of psi on the plane y_k = beta . y' + alpha over the base box
/// [blo, bhi] (Bernstein coefficients, strictly, or Taylor models); 0 if the
/// bounds do not show one.
template <std::floating_point T, int D>
int plane_sign(Context<T>& ctx, const Func<T, D>& f, int k, const Bound<T, D>& plane, const VecD<T, D - 1>& blo,
               const VecD<T, D - 1>& bhi)
{
    const Func<T, D - 1> slice = restrict_to<T, D>(f, k, plane);
    if (is_analytic(ctx, f.ls))
    {
        std::array<T, 3> origin;
        std::array<T, 3 * (D - 1)> matrix;
        affine_image<T, D - 1>(slice, blo, bhi, origin, matrix);
        AffineBounds<T> b;
        if (!affine_bounds(ctx.phis[static_cast<std::size_t>(f.ls)], std::span<const T>(origin),
                           std::span<const T>(matrix), D - 1, b))
            return 0;
        return b.sign;
    }
    bernstein_form<T, D - 1>(ctx, slice, blo, bhi, ctx.form[D - 1], ctx.form_work[D - 1]);
    bool pos = true, neg = true;
    for (const T v : ctx.form[D - 1].coeffs)
    {
        pos &= v > T(0);
        neg &= v < T(0);
    }
    return pos ? 1 : (neg ? -1 : 0);
}

/// Two roots per line along k for a function psi whose derivative along k is
/// monotone (second-derivative margin; curvature its sign), judged from the
/// roots themselves rather than by bisection. On a grid of height lines over
/// the base, the extreme point of psi (the root of d_k psi, or the box end
/// nearer to it) lies between the two sheets; the centre line must cross
/// both sheets inside the box (else psi has one root there, a case for the
/// margin: -1). Over each cell of the grid, a
/// plane through the extreme points of its corners does too. If psi has the
/// sign -curvature on every such plane (Bernstein coefficients or Taylor
/// models), the extreme value of psi on every height line has that sign, so
/// its two roots never merge in the box. Returns the smallest
/// |d_k psi| / |grad psi| at the roots on the lines through the base's
/// corners and centre (1 if none; 0 where another axis is much better there,
/// near a silhouette along k), or -1 if a plane does not show the sign.
template <std::floating_point T, int D>
T two_root_score(Context<T>& ctx, const Func<T, D>& f, const VecD<T, D>& lo, const VecD<T, D>& hi, int k,
                 int curvature)
{
    // |d_k psi| / |grad psi| at a root, or 0 where another axis is much better
    auto slope_ratio = [&](const VecD<T, D>& y)
    {
        VecD<T, D> g;
        for (int j = 0; j < D; ++j)
            g[j] = derivative_along<T, D>(ctx, f, y, j);
        const T norm = scaled_norm(std::span<const T>(g));
        if (!(norm > T(0)))
            return T(0);
        T best = T(0);
        for (int j = 0; j < D; ++j)
            best = std::max(best, std::abs(g[j]));
        return std::abs(g[k]) >= T(0.75) * best ? std::abs(g[k]) / norm : T(0);
    };
    auto line_point = [&](const VecD<T, D - 1>& yb, T t)
    {
        VecD<T, D> y = insert<T, D>(yb, k, t);
        return y;
    };
    // the extreme point on the line through yb, and the roots on either side
    // (their slopes lower score)
    T score = T(1);
    int found = 0; // roots found on the last line examined
    auto extreme_on = [&](const VecD<T, D - 1>& yb, bool roots) -> T
    {
        found = 0;
        auto slope = [&](T t) { return derivative_along<T, D>(ctx, f, line_point(yb, t), k); };
        auto val = [&](T t) { return value<T, D>(ctx, f, line_point(yb, t)); };
        const T dl = slope(lo[k]), dh = slope(hi[k]);
        // without an extreme point in the box, the end nearer to it (a NaN
        // marks the line for the caller)
        if (!opposite(dl, dh))
            return std::numeric_limits<T>::quiet_NaN();
        const T c = illinois_root(slope, lo[k], hi[k], dl, dh);
        if (roots)
        {
            const std::array<T, 3> ends = {lo[k], c, hi[k]};
            for (int piece = 0; piece < 2; ++piece)
            {
                const T a = ends[static_cast<std::size_t>(piece)], b = ends[static_cast<std::size_t>(piece + 1)];
                if (!(b > a))
                    continue;
                const T fa = val(a), fb = val(b);
                if (opposite(fa, fb))
                {
                    score = std::min(score, slope_ratio(line_point(yb, illinois_root(val, a, b, fa, fb))));
                    ++found;
                }
            }
        }
        return c;
    };
    // the grid of lines over the base: M cells per base direction
    constexpr int M = 4;
    constexpr int nodes = D == 3 ? (M + 1) * (M + 1) : M + 1;
    VecD<T, D - 1> blo, bhi;
    for (int jb = 0; jb < D - 1; ++jb)
    {
        blo[static_cast<std::size_t>(jb)] = lo[jb < k ? jb : jb + 1];
        bhi[static_cast<std::size_t>(jb)] = hi[jb < k ? jb : jb + 1];
    }
    std::array<VecD<T, D - 1>, nodes> base{};
    std::array<T, nodes> extreme{};
    for (int n = 0; n < nodes; ++n)
    {
        int m = n;
        for (int jb = 0; jb < D - 1; ++jb, m /= (M + 1))
            base[static_cast<std::size_t>(n)][static_cast<std::size_t>(jb)]
                = blo[static_cast<std::size_t>(jb)]
                  + (bhi[static_cast<std::size_t>(jb)] - blo[static_cast<std::size_t>(jb)]) * T(m % (M + 1)) / T(M);
        // roots on the lines through the base's corners and centre
        bool corner = true;
        m = n;
        for (int jb = 0; jb < D - 1; ++jb, m /= (M + 1))
            corner &= m % (M + 1) == 0 || m % (M + 1) == M;
        bool centre = true;
        m = n;
        for (int jb = 0; jb < D - 1; ++jb, m /= (M + 1))
            centre &= m % (M + 1) == M / 2;
        T c = extreme_on(base[static_cast<std::size_t>(n)], corner || centre);
        // two sheets: both crossed by the centre line inside the box; else one
        // root there, a case for the margin
        if (centre && (!std::isfinite(c) || found < 2))
            return T(-1);
        if (!std::isfinite(c))
        {
            const VecD<T, D - 1>& yb = base[static_cast<std::size_t>(n)];
            const T dl = derivative_along<T, D>(ctx, f, line_point(yb, lo[k]), k),
                    dh = derivative_along<T, D>(ctx, f, line_point(yb, hi[k]), k);
            c = std::abs(dl) <= std::abs(dh) ? lo[k] : hi[k];
        }
        extreme[static_cast<std::size_t>(n)] = c;
    }
    // a plane per grid cell through the extreme points of three of its corners
    constexpr int cells = D == 3 ? M * M : M;
    for (int cell = 0; cell < cells; ++cell)
    {
        const int i = cell % M, j = D == 3 ? cell / M : 0;
        const int n00 = D == 3 ? i + (M + 1) * j : i;
        VecD<T, D - 1> clo = base[static_cast<std::size_t>(n00)], chi;
        Bound<T, D> plane;
        plane.kind = 7;
        if constexpr (D == 3)
        {
            const int n10 = n00 + 1, n01 = n00 + (M + 1);
            chi = {base[static_cast<std::size_t>(n10)][0], base[static_cast<std::size_t>(n01)][1]};
            plane.beta[0] = (extreme[static_cast<std::size_t>(n10)] - extreme[static_cast<std::size_t>(n00)])
                            / (chi[0] - clo[0]);
            plane.beta[1] = (extreme[static_cast<std::size_t>(n01)] - extreme[static_cast<std::size_t>(n00)])
                            / (chi[1] - clo[1]);
        }
        else
        {
            const int n10 = n00 + 1;
            chi = {base[static_cast<std::size_t>(n10)][0]};
            plane.beta[0] = (extreme[static_cast<std::size_t>(n10)] - extreme[static_cast<std::size_t>(n00)])
                            / (chi[0] - clo[0]);
        }
        plane.alpha = extreme[static_cast<std::size_t>(n00)];
        for (int jb = 0; jb < D - 1; ++jb)
            plane.alpha -= plane.beta[static_cast<std::size_t>(jb)] * clo[static_cast<std::size_t>(jb)];
        if (plane_sign<T, D>(ctx, f, k, plane, clo, chi) != -curvature)
            return T(-1);
    }
    return score;
}

/// The height direction of a 3D box to integrate without bisecting it
/// (Options::silhouette_splits): an axis along which each of its curved functions
/// is monotone (margin_of positive: constant, or passive, counts) or has a monotone
/// derivative, so that its lines take at most two roots of each. Among them the
/// largest score, on a 7^3 grid: near each function's zero set the mean |d_k f| /
/// |grad f| as in Saye's score, near a pair's meeting curve the mean length of its
/// tangent T = grad f x grad g seen along k, relative to |T|, for lines that cross the
/// curve rather than run along it; less a little per function that needs a
/// silhouette. @p two receives, per entry of
/// Analysis::funcs, whether it does. False if the box does not qualify.
template <std::floating_point T, typename Margin>
bool silhouette_direction(Context<T>& ctx, const Problem<T, 3>& p, const Analysis<T, 3>& an, Margin&& margin_of,
                          int& k, std::vector<char>& two)
{
    // A failed attempt is thrown away, so a cell gets few. Volumes gain most from what
    // succeeds lower down; on interfaces a bisected box costs little anyway.
    if (!ctx.silhouettes || an.curved.empty() || ctx.abandoned >= (ctx.surface >= 0 ? 1 : 4))
        return false;
    for (const int c : an.curved)
    {
        const Func<T, 3>& f = an.funcs[static_cast<std::size_t>(c)];
        // the derived level sets are along the axes of the cell's box
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                if (f.A[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] != (i == j ? T(1) : T(0)))
                    return false;
    }
    const std::size_t n = an.curved.size();
    // which axes qualify: this needs bounds only
    std::array<std::vector<char>, 3> candidates;
    std::array<int, 3> needed{};
    std::array<int, 3> axes{};
    int n_axes = 0;
    for (int kk = 0; kk < 3; ++kk)
    {
        std::vector<char>& candidate = candidates[static_cast<std::size_t>(kk)];
        candidate.assign(an.funcs.size(), 0);
        bool qualifies = true;
        for (std::size_t c = 0; c < n && qualifies; ++c)
        {
            if (margin_of(c, kk) > T(0))
                continue;
            const std::size_t i = static_cast<std::size_t>(an.curved[c]);
            const Func<T, 3>& f = an.funcs[i];
            int curvature = 0;
            if (ctx.flat[static_cast<std::size_t>(f.ls)] || ctx.derived[static_cast<std::size_t>(3 * f.ls + kk)] < 0
                || !(second_margin<T, 3>(ctx, f, p.lo, p.hi, kk, curvature) > T(0)))
                qualifies = false;
            candidate[i] = 1;
            ++needed[static_cast<std::size_t>(kk)];
        }
        if (qualifies)
            axes[static_cast<std::size_t>(n_axes++)] = kk;
    }
    if (n_axes == 0)
        return false;
    if (n_axes == 1)
    {
        k = axes[0];
        two = candidates[static_cast<std::size_t>(k)];
        return true;
    }
    // Several: the largest score on a grid
    static constexpr int G = 7, total = G * G * G;
    const auto grid_point = [&p](int m)
    {
        VecD<T, 3> y;
        for (std::size_t j = 0; j < 3; ++j, m /= G)
            y[j] = p.lo[j] + (p.hi[j] - p.lo[j]) * (T(m % G) + T(0.5)) / T(G);
        return y;
    };
    std::vector<T> values(n * total);
    std::vector<T> largest(n, T(0));
    for (std::size_t c = 0; c < n; ++c)
        for (int m = 0; m < total; ++m)
        {
            const T v = value<T, 3>(ctx, an.funcs[static_cast<std::size_t>(an.curved[c])], grid_point(m));
            values[c * total + static_cast<std::size_t>(m)] = v;
            largest[c] = std::max(largest[c], std::abs(v));
        }
    const auto near_zero = [&](std::size_t c, int m)
    { return std::abs(values[c * total + static_cast<std::size_t>(m)]) <= T(0.2) * largest[c]; };
    std::vector<VecD<T, 3>> mean(n, VecD<T, 3>{});
    std::vector<int> count(n, 0);
    std::vector<VecD<T, 3>> grads(n);
    std::vector<int> close;
    close.reserve(n);
    VecD<T, 3> pair_sum{};
    int pair_count = 0;
    for (int m = 0; m < total; ++m)
    {
        const VecD<T, 3> y = grid_point(m);
        bool in_region = true;
        for (const Half<T, 3>& h : p.clips)
            in_region = in_region && h.c[0] * y[0] + h.c[1] * y[1] + h.c[2] * y[2] <= h.d;
        if (!in_region)
            continue;
        close.clear();
        for (std::size_t c = 0; c < n; ++c)
            if (near_zero(c, m))
            {
                func_gradient<T, 3>(ctx, an.funcs[static_cast<std::size_t>(an.curved[c])], y, grads[c]);
                const T norm = scaled_norm(std::span<const T>(grads[c]));
                if (!(norm > T(0)))
                    continue;
                for (std::size_t j = 0; j < 3; ++j)
                {
                    grads[c][j] /= norm;
                    mean[c][j] += std::abs(grads[c][j]);
                }
                ++count[c];
                close.push_back(static_cast<int>(c));
            }
        for (std::size_t a = 0; a < close.size(); ++a)
            for (std::size_t b = a + 1; b < close.size(); ++b)
            {
                const VecD<T, 3>& f = grads[static_cast<std::size_t>(close[a])];
                const VecD<T, 3>& g = grads[static_cast<std::size_t>(close[b])];
                const VecD<T, 3> t = {f[1] * g[2] - f[2] * g[1], f[2] * g[0] - f[0] * g[2], f[0] * g[1] - f[1] * g[0]};
                const T norm = std::sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
                if (!(norm > T(0)))
                    continue;
                // the length of the curve seen along each axis
                for (std::size_t j = 0; j < 3; ++j)
                    pair_sum[j] += std::sqrt(std::max(T(0), T(1) - t[j] * t[j] / (norm * norm)));
                ++pair_count;
            }
    }
    T best = -infinity<T>;
    bool found = false;
    for (int a = 0; a < n_axes; ++a)
    {
        const int kk = axes[static_cast<std::size_t>(a)];
        const std::size_t K = static_cast<std::size_t>(kk);
        T score = T(1);
        for (std::size_t c = 0; c < n; ++c)
            if (count[c] > 0)
                score = std::min(score, mean[c][K] / T(count[c]));
        if (pair_count > 0)
            score *= pair_sum[K] / T(pair_count);
        score -= T(0.1) * T(needed[K]);
        if (score > best)
        {
            best = score;
            k = kk;
            two = candidates[K];
            found = true;
        }
    }
    return found;
}

template <std::floating_point T, int D>
Analysis<T, D> analyse(Context<T>& ctx, const Problem<T, D>& p, bool top, int depth)
{
    Analysis<T, D> an;
    an.surface = ctx.surface;
    BoxBernstein<T>& form = ctx.form[D];
    VecD<T, D> lengths;
    for (int j = 0; j < D; ++j)
        lengths[j] = p.hi[j] - p.lo[j];
    if (top)
        an.signs.assign(ctx.phis.size(), 0);
    const bool clipped = ctx.opt.taylor_subdivisions > 1 && box_clipped<T, D>(p.lo, p.hi, p.clips);
    bool surface_vanishes = false, partner_vanishes = false;
    for (const Func<T, D>& f : p.funcs)
    {
        if (f.linear)
        {
            if (linear_may_vanish<T, D>(f, p.lo, p.hi))
                an.funcs.push_back(f);
            continue;
        }
        const T scale = ctx.scales[static_cast<std::size_t>(f.ls)];
        if (is_analytic(ctx, f.ls))
        {
            int sign = 0;
            VecD<T, D> ratio{}, upper{};
            T magnitude = T(0);
            if (clipped)
                analytic_local_bounds<T, D>(ctx, f, p.lo, p.hi, p.clips, ctx.opt.taylor_subdivisions, sign, ratio,
                                            magnitude, upper);
            else
                analytic_bounds<T, D>(ctx, f, p.lo, p.hi, sign, ratio, magnitude, upper);
            // phi = 0 on a face, say: no root inside the region
            if (magnitude <= zero_function_tol<T> * scale)
                continue;
            if (sign != 0)
            {
                if (top)
                    an.signs[static_cast<std::size_t>(f.ls)] = sign;
                continue;
            }
            an.curved.push_back(static_cast<int>(an.funcs.size()));
            an.funcs.push_back(f);
            an.ratios.push_back(ratio);
            an.uppers.push_back(upper);
            surface_vanishes |= f.ls == ctx.surface;
            partner_vanishes |= f.ls == ctx.partner;
            continue;
        }
        bernstein_form<T, D>(ctx, f, p.lo, p.hi, form, ctx.form_work[D]);
        const std::span<const T> coeffs(form.coeffs);
        // phi vanishes identically here (phi = 0 on a face, say): only rounding
        // noise is left, which would block certification; it has no root inside
        if (max_abs(coeffs) <= zero_function_tol<T> * scale)
            continue;
        VecD<T, D> ratio{}, upper{};
        if (ctx.opt.mask_subdivisions > 1)
        {
            if (!local_margins<T, D>(ctx, form, p.lo, p.hi, p.clips, ctx.opt.mask_subdivisions, ratio, upper))
                continue; // no zero of this function inside the cell
        }
        else
        {
            if (!may_vanish(coeffs))
            {
                // one sign on the box, zeros at most where it touches zero
                const T largest = *std::max_element(coeffs.begin(), coeffs.end(), [](T a, T b)
                                                    { return std::abs(a) < std::abs(b); });
                if (top)
                    an.signs[static_cast<std::size_t>(f.ls)] = largest > T(0) ? 1 : -1;
                continue;
            }
            margins(form, std::span<const T>(lengths), std::span<T>(ratio), ctx.deriv[D], std::span<T>(upper));
        }
        an.curved.push_back(static_cast<int>(an.funcs.size()));
        an.funcs.push_back(f);
        an.ratios.push_back(ratio);
        an.uppers.push_back(upper);
        surface_vanishes |= f.ls == ctx.surface;
        partner_vanishes |= f.ls == ctx.partner;
    }

    // Top level: the level sets' signs may decide the terms on the box, and
    // level sets the remaining terms do not name need no breakpoints. A curve
    // crosses the box only where both its level sets may vanish.
    if (top)
    {
        an.decided = residual_terms(ctx.terms, an.signs, ctx.surface,
                                    surface_vanishes && (ctx.partner < 0 || partner_vanishes), an.terms);
        std::uint64_t needed = 0;
        if (an.decided == 0)
            for (const Requirement& r : an.terms)
                needed |= r.negative | r.positive;
        if (ctx.surface >= 0)
            needed |= std::uint64_t(1) << ctx.surface;
        if (ctx.partner >= 0)
            needed |= std::uint64_t(1) << ctx.partner;
        if (an.decided == -1)
            needed = 0;
        std::vector<Func<T, D>> funcs;
        std::vector<VecD<T, D>> ratios, uppers;
        std::vector<int> curved;
        for (std::size_t i = 0, c = 0; i < an.funcs.size(); ++i)
        {
            const bool is_curved = c < an.curved.size() && an.curved[c] == static_cast<int>(i);
            if (!is_curved || ((needed >> an.funcs[i].ls) & 1))
            {
                if (is_curved)
                {
                    curved.push_back(static_cast<int>(funcs.size()));
                    ratios.push_back(an.ratios[c]);
                    uppers.push_back(an.uppers[c]);
                }
                funcs.push_back(an.funcs[i]);
            }
            c += is_curved;
        }
        an.funcs = std::move(funcs);
        an.ratios = std::move(ratios);
        an.uppers = std::move(uppers);
        an.curved = std::move(curved);
    }

    // A function constant along y_k has no root on its height lines and needs no
    // margin there; its restrictions to the bounds carry it to the base. Not the
    // level set of an interface, whose roots the top level integrates.
    const auto is_passive = [&](std::size_t c)
    { return top && an.surface >= 0 && an.funcs[static_cast<std::size_t>(an.curved[c])].ls != an.surface; };
    const auto margin_of = [&](std::size_t c, int kk)
    {
        if (is_passive(c))
            return T(1);
        const VecD<T, D>& upper = an.uppers[c];
        const T norm = scaled_norm(std::span<const T>(upper));
        const bool constant = norm > T(0) && norm < infinity<T> && upper[kk] <= tiny<T> * norm
                              && !(top && an.funcs[static_cast<std::size_t>(an.curved[c])].ls == an.surface);
        return constant ? T(1) : an.ratios[c][kk];
    };

    // surface functions: dropped where they have no root, where their two zero
    // sets coincide (no kink on the height lines, and the function is rounding
    // noise, whose sign changes would each split them), and a curve part's
    // where the curve does not cross the box above
    if constexpr (D <= 2)
    {
        for (const Surface<T, D>& s : p.surfaces)
        {
            if (surface_vanishes_identically<T, D>(ctx, s, p.lo, p.hi)
                || (partner_surface<T, D>(ctx, s) && under_signed<T, D>(ctx, s, p.lo, p.hi)))
                continue;
            VecD<T, D> ratio{};
            if (!surface_margins<T, D>(ctx, s, p.lo, p.hi, ratio))
                continue;
            an.surfaces.push_back(s);
            an.surface_ratios.push_back(ratio);
        }
    }

    // height direction: the best certified margin over all curved functions
    const auto choose_direction = [&]()
    {
        an.best = T(-1);
        an.k = 0;
        for (int kk = 0; kk < D; ++kk)
        {
            T best = T(1);
            for (std::size_t c = 0; c < an.ratios.size(); ++c)
                best = std::min(best, margin_of(c, kk));
            for (const VecD<T, D>& r : an.surface_ratios)
                best = std::min(best, r[kk]);
            if (best > an.best || (best == an.best && p.hi[kk] - p.lo[kk] > p.hi[an.k] - p.lo[an.k]))
            {
                an.best = best;
                an.k = kk;
            }
        }
        an.certified
            = (an.ratios.empty() && an.surface_ratios.empty()) || (an.best > T(0) && an.best >= level_margin<T, D>(ctx));
    };
    choose_direction();
    // A curve part takes the roots of either of its level sets: the second
    // where only it has a margin (the first may have two roots per line in the
    // box, which no surface function follows).
    if (!an.certified && top && ctx.partner >= 0)
    {
        an.surface = ctx.partner;
        choose_direction();
        if (!an.certified)
        {
            an.surface = ctx.surface;
            choose_direction();
        }
    }
    an.two_roots.assign(an.funcs.size(), 0);
    an.independent.assign(an.funcs.size(), 0);
    an.passive.assign(an.funcs.size(), 0);
    const auto mark_independent = [&]()
    {
        for (std::size_t c = 0; c < an.ratios.size(); ++c)
        {
            const std::size_t i = static_cast<std::size_t>(an.curved[c]);
            an.passive[i] = is_passive(c);
            an.independent[i] = !an.passive[i] && margin_of(c, an.k) == T(1) && an.ratios[c][an.k] < T(1);
        }
    };
    // curve parts: no two roots per line, which no surface function follows
    if (an.certified || depth < ctx.opt.two_roots_depth || an.ratios.empty() || (top && ctx.partner >= 0))
    {
        mark_independent();
        return an;
    }

    // Two roots per line, where no direction gives one: bisection rarely
    // separates two sheets of one level set, so the roots themselves decide. A
    // direction qualifies where every curved function has a margin, or a
    // monotone derivative along it (second-derivative margin) with roots that
    // never merge in the box and are not steep (two_root_score). Each line
    // then splits at the extreme point.
    T best2 = -1;
    int k2 = -1;
    std::vector<int> chosen;
    for (int kk = 0; kk < D; ++kk)
    {
        T score = T(1);
        for (const VecD<T, D>& r : an.surface_ratios)
            score = std::min(score, r[kk]);
        if (!(score >= level_margin<T, D>(ctx)))
            continue;
        std::vector<int> two;
        for (std::size_t c = 0; c < an.ratios.size() && score >= level_margin<T, D>(ctx); ++c)
        {
            if (margin_of(c, kk) >= level_margin<T, D>(ctx))
            {
                score = std::min(score, margin_of(c, kk));
                continue;
            }
            const Func<T, D>& f = an.funcs[static_cast<std::size_t>(an.curved[c])];
            int cv = 0;
            if (!(second_margin<T, D>(ctx, f, p.lo, p.hi, kk, cv) >= level_margin<T, D>(ctx)))
            {
                score = T(-1);
                break;
            }
            score = std::min(score, two_root_score<T, D>(ctx, f, p.lo, p.hi, kk, cv));
            two.push_back(static_cast<int>(c));
        }
        if (score >= level_margin<T, D>(ctx) && score > best2)
        {
            best2 = score;
            k2 = kk;
            chosen = two;
        }
    }
    if (k2 < 0)
    {
        if constexpr (D == 3)
        {
            int k3 = 0;
            std::vector<char> two;
            if (top && silhouette_direction<T>(ctx, p, an, margin_of, k3, two))
            {
                an.two_roots = two;
                an.unmargined = true;
                an.k = k3;
                an.best = T(0);
                an.certified = true;
                mark_independent();
                ++ctx.stats->silhouettes;
            }
        }
        return an;
    }
    for (const int c : chosen)
        an.two_roots[static_cast<std::size_t>(an.curved[static_cast<std::size_t>(c)])] = 1;
    an.k = k2;
    an.best = best2;
    an.certified = true;
    mark_independent();
    ++ctx.stats->two_roots;
    return an;
}

/// Margins of @p f near its zero set, sampled: on a 9^D grid of the clipped
/// region, at the points where |f| is within a fifth of its largest value, the
/// smallest |d_j f| / |grad f| along each axis j (infinity if no point is
/// that close). False if f takes one sign at every point (its zero set is
/// outside the cell, or between the points).
template <std::floating_point T, int D>
bool sampled_margins(const Context<T>& ctx, const Problem<T, D>& p, const Func<T, D>& f, VecD<T, D>& rmin)
{
    const int G = 9;
    std::vector<VecD<T, D>> points;
    int total = 1;
    for (int j = 0; j < D; ++j)
        total *= G;
    for (int n = 0; n < total; ++n)
    {
        int m = n;
        VecD<T, D> y;
        for (int j = 0; j < D; ++j)
        {
            y[j] = p.lo[j] + (p.hi[j] - p.lo[j]) * (T(m % G) + T(0.5)) / G;
            m /= G;
        }
        bool inside = true;
        for (const Half<T, D>& h : p.clips)
        {
            T sum = T(0);
            for (int j = 0; j < D; ++j)
                sum += h.c[j] * y[j];
            inside &= sum <= h.d;
        }
        if (inside)
            points.push_back(y);
    }
    std::vector<T> values(points.size());
    T vmax = T(0);
    bool pos = false, neg = false;
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        values[i] = value<T, D>(ctx, f, points[i]);
        vmax = std::max(vmax, std::abs(values[i]));
        pos |= values[i] > T(0);
        neg |= values[i] < T(0);
    }
    rmin.fill(infinity<T>);
    if (!(pos && neg))
        return false;
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        if (std::abs(values[i]) > T(0.2) * vmax)
            continue;
        VecD<T, D> g;
        for (int j = 0; j < D; ++j)
            g[j] = derivative_along<T, D>(ctx, f, points[i], j);
        const T norm = scaled_norm(std::span<const T>(g));
        for (int j = 0; j < D && norm > T(0); ++j)
            rmin[j] = std::min(rmin[j], std::abs(g[j]) / norm);
    }
    return true;
}

/// The height direction of a box integrated uncertified: no axis has a
/// certified margin, but its lines must still cross the zero sets, as
/// transversally as the box allows. Where the samples of a curved function
/// change sign, the axis with the largest sampled margin near the zero sets
/// (sampled_margins). Otherwise @p k, the longest side, unless the bounds
/// show that a function hardly varies along it (its |d_k psi| below a tenth
/// of that along the best axis): then the longest of the other axes. At a
/// pyramid's apex, where the form (1 - z)^n phi has a degenerate zero, every
/// bound is 0 and the longest side may be an axis along which the level set
/// does not vary at all: its lines would find no roots.
template <std::floating_point T, int D>
int uncertified_direction(const Context<T>& ctx, const Problem<T, D>& p, const Analysis<T, D>& an, int k)
{
    const auto passive = [&](std::size_t i) { return i < an.passive.size() && an.passive[i]; };
    VecD<T, D> score;
    score.fill(infinity<T>);
    bool seen = false;
    for (std::size_t i = 0; i < p.funcs.size(); ++i)
    {
        if (p.funcs[i].linear || passive(i))
            continue;
        VecD<T, D> rmin;
        if (!sampled_margins<T, D>(ctx, p, p.funcs[i], rmin))
            continue;
        seen = true;
        for (int j = 0; j < D; ++j)
            score[j] = std::min(score[j], rmin[j]);
    }
    if (seen)
    {
        for (int j = 0; j < D; ++j)
            if (score[j] > score[k])
                k = j;
        return k;
    }
    // the zero sets fall between the samples: how much each function may vary
    // along each axis, relative to its gradient, from the bounds
    VecD<T, D> variation;
    variation.fill(T(1));
    for (std::size_t c = 0; c < an.uppers.size(); ++c)
    {
        if (passive(static_cast<std::size_t>(an.curved[c])))
            continue;
        const VecD<T, D>& upper = an.uppers[c];
        const T norm = scaled_norm(std::span<const T>(upper));
        for (int j = 0; j < D; ++j)
            variation[j] = std::min(variation[j], norm > T(0) && norm < infinity<T> ? upper[j] / norm : T(0));
    }
    T best = T(0);
    for (int j = 0; j < D; ++j)
        best = std::max(best, variation[j]);
    if (variation[k] >= T(0.1) * best)
        return k;
    int longest = k;
    for (int j = 0; j < D; ++j)
        if (variation[j] >= T(0.1) * best && (longest == k || p.hi[j] - p.lo[j] > p.hi[longest] - p.lo[longest]))
            longest = j;
    return longest;
}

/// Why did certification fail? Sample the clipped region on a 9^D grid and look
/// at the function with the smallest margin for the chosen axis near its zero set.
template <std::floating_point T, int D>
void diagnose_failure(Context<T>& ctx, const Problem<T, D>& p, const Analysis<T, D>& an, int depth)
{
    std::map<std::string, std::int64_t>& causes = ctx.stats->causes;
    const std::string level = "level " + std::to_string(D);
    // the origin of a function, with its level set's index
    const auto name = [](const Func<T, D>& f)
    { return origin_name(f.origin) + (f.linear ? std::string() : "[" + std::to_string(f.ls) + "]"); };
    const std::vector<VecD<T, D>>& ratios = an.ratios;
    if (ratios.empty())
    {
        ++causes[level + " | surface function | depth " + std::to_string(depth)];
        return;
    }
    const int k = an.k;
    std::size_t blocker = 0;
    for (std::size_t i = 1; i < ratios.size(); ++i)
        if (ratios[i][k] < ratios[blocker][k])
            blocker = i;
    const Func<T, D>& f = p.funcs[an.curved[blocker]];
    VecD<T, D> rmin;
    std::string category;
    if (!sampled_margins<T, D>(ctx, p, f, rmin))
        category = "zero set outside the cell";
    else
    {
        T other = T(0);
        for (int j = 0; j < D; ++j)
            other = std::max(other, rmin[j]);
        if (rmin[k] >= level_margin<T, D>(ctx))
            category = "margin holds near the zero set in the cell";
        else if (other >= level_margin<T, D>(ctx))
            category = "another axis has the margin there";
        else
            category = "margin fails near the zero set in the cell";
    }
    ++causes[level + " | " + name(f) + " | "
             + (an.best > T(0) ? "bound below margin" : "no certain sign") + " | " + category];
    ++causes[level + " | depth " + std::to_string(depth)];

    // From the certified bounds alone: does one function fail on every axis, or
    // does each function have an axis but no axis suits all (a conflict)?
    std::string axes;
    for (std::size_t i = 0; i < ratios.size() && axes.empty(); ++i)
    {
        T top = T(0);
        for (int j = 0; j < D; ++j)
            top = std::max(top, ratios[i][j]);
        if (top < level_margin<T, D>(ctx))
            axes = "one function fails on every axis: " + name(p.funcs[an.curved[i]]);
    }
    if (axes.empty())
    {
        std::vector<std::string> names;
        for (int j = 0; j < D; ++j)
        {
            std::size_t b = 0;
            for (std::size_t i = 1; i < ratios.size(); ++i)
                if (ratios[i][j] < ratios[b][j])
                    b = i;
            names.push_back(name(p.funcs[an.curved[b]]));
        }
        std::sort(names.begin(), names.end());
        names.erase(std::unique(names.begin(), names.end()), names.end());
        axes = "conflict:";
        for (const std::string& n : names)
            axes += " " + n;
    }
    // is the box cut by the boundary of the clipped region?
    bool cut = false;
    for (const Half<T, D>& h : p.clips)
    {
        T vmin = T(0), vmax_c = T(0);
        for (int j = 0; j < D; ++j)
        {
            vmin += h.c[j] * (h.c[j] >= T(0) ? p.lo[j] : p.hi[j]);
            vmax_c += h.c[j] * (h.c[j] >= T(0) ? p.hi[j] : p.lo[j]);
        }
        cut |= vmin < h.d - compare_tol<T> && vmax_c > h.d + compare_tol<T>;
    }
    ++causes["axes | " + level + (depth >= ctx.opt.max_depth ? " | at depth limit" : " | before limit") + " | "
             + axes + (cut ? " | box cut by a face" : " | box inside")];
}

/// A level-2 problem in the diagonal frame y = T z, T = [[1/2, -1/2], [1/2, 1/2]]:
/// z_0 runs along (1, 1) and z_1 along (-1, 1). The box becomes four clips.
template <typename T>
Problem<T, 2> rotate_diagonal(const Problem<T, 2>& p)
{
    auto transpose_times = [](const VecD<T, 2>& c)
    { return VecD<T, 2>{T(0.5) * (c[0] + c[1]), T(0.5) * (c[1] - c[0])}; };
    Problem<T, 2> r;
    for (const Half<T, 2>& h : p.clips)
        r.clips.push_back({transpose_times(h.c), h.d});
    for (int j = 0; j < 2; ++j)
    {
        VecD<T, 2> e{};
        e[j] = T(1);
        r.clips.push_back({transpose_times(e), p.hi[j]});
        e[j] = T(-1);
        r.clips.push_back({transpose_times(e), -p.lo[j]});
    }
    // z_0 = y_0 + y_1, z_1 = y_1 - y_0
    r.lo = {p.lo[0] + p.lo[1], p.lo[1] - p.hi[0]};
    r.hi = {p.hi[0] + p.hi[1], p.hi[1] - p.lo[0]};
    for (const Func<T, 2>& f : p.funcs)
    {
        Func<T, 2> g = f;
        if (f.linear)
            g.a = transpose_times(f.a);
        else
            for (int i = 0; i < 3; ++i)
                g.A[i] = transpose_times(f.A[i]);
        r.funcs.push_back(g);
    }
    // functions of the level above (surface functions, heights): the two
    // columns of y among their inner coordinates, t being their coordinate k
    const auto rotate_inner = [&transpose_times](Func<T, 3>& f, int k)
    {
        const int j0 = k == 0 ? 1 : 0, j1 = k == 2 ? 1 : 2;
        for (int i = 0; i < 3; ++i)
        {
            const VecD<T, 2> c = transpose_times(VecD<T, 2>{f.A[i][j0], f.A[i][j1]});
            f.A[i][j0] = c[0];
            f.A[i][j1] = c[1];
        }
    };
    for (const Surface<T, 2>& s : p.surfaces)
    {
        Surface<T, 2> g = s;
        rotate_inner(g.under, s.k);
        rotate_inner(g.on, s.k);
        r.surfaces.push_back(g);
    }
    for (const Height<T, 2>& h : p.heights)
    {
        Height<T, 2> g = h;
        rotate_inner(g.f, h.k);
        r.heights.push_back(g);
    }
    return r;
}

/// A problem in a rotated frame y = R z, R orthogonal with the new axes as its
/// columns. The box becomes 2 D clips; the functions, the clips and the frame
/// map are composed with R.
template <typename T, int D>
Problem<T, D> rotate_frame(const Problem<T, D>& p, const std::array<VecD<T, D>, D>& R)
{
    // a row vector c of y, as a row vector of z: c R
    auto times = [&R](const VecD<T, D>& c)
    {
        VecD<T, D> r{};
        for (int m = 0; m < D; ++m)
            for (int j = 0; j < D; ++j)
                r[m] += c[j] * R[j][m];
        return r;
    };
    Problem<T, D> r;
    r.frame = p.frame;
    for (int i = 0; i < 3; ++i)
        r.frame.A[i] = times(p.frame.A[i]);
    for (const Half<T, D>& h : p.clips)
        r.clips.push_back({times(h.c), h.d});
    for (int j = 0; j < D; ++j)
    {
        VecD<T, D> e{};
        e[j] = T(1);
        r.clips.push_back({times(e), p.hi[j]});
        e[j] = T(-1);
        r.clips.push_back({times(e), -p.lo[j]});
    }
    // the bounding box of the box's corners, z = R^T y
    r.lo.fill(infinity<T>);
    r.hi.fill(-infinity<T>);
    for (int corner = 0; corner < (1 << D); ++corner)
        for (int m = 0; m < D; ++m)
        {
            T z = T(0);
            for (int j = 0; j < D; ++j)
                z += R[j][m] * (((corner >> j) & 1) ? p.hi[j] : p.lo[j]);
            r.lo[m] = std::min(r.lo[m], z);
            r.hi[m] = std::max(r.hi[m], z);
        }
    for (const Func<T, D>& f : p.funcs)
    {
        Func<T, D> g = f;
        if (f.linear)
            g.a = times(f.a);
        else
            for (int i = 0; i < 3; ++i)
                g.A[i] = times(f.A[i]);
        r.funcs.push_back(g);
    }
    if constexpr (D <= 2)
    {
        // functions of the level above (surface functions, heights): the
        // columns of y among their inner coordinates, t being their coordinate k
        const auto rotate_inner = [&times](Func<T, D + 1>& f, int k)
        {
            for (int i = 0; i < 3; ++i)
            {
                VecD<T, D> c;
                for (int j = 0; j < D; ++j)
                    c[j] = f.A[i][j < k ? j : j + 1];
                c = times(c);
                for (int j = 0; j < D; ++j)
                    f.A[i][j < k ? j : j + 1] = c[j];
            }
        };
        for (const Surface<T, D>& s : p.surfaces)
        {
            Surface<T, D> g = s;
            rotate_inner(g.under, s.k);
            rotate_inner(g.on, s.k);
            r.surfaces.push_back(g);
        }
        for (const Height<T, D>& h : p.heights)
        {
            Height<T, D> g = h;
            rotate_inner(g.f, h.k);
            r.heights.push_back(g);
        }
    }
    return r;
}

/// Whether level-2 boxes try frames before bisecting: only where they are not
/// split at their folds and crossings, which needs no frame.
template <std::floating_point T>
bool level2_frames(const Context<T>& ctx)
{
    return ctx.opt.diagonal_frames && !(ctx.opt.fold_splits && ctx.opt.anchor_folds > 0 && ctx.vis_order == 0);
}

/// Whether the zero sets of two level sets may meet in a box: on a grid of
/// its sub-boxes, one shows a sign change of both at its corners. Where they do
/// not meet, bisection separates them, and a rotated frame only costs points.
template <std::floating_point T>
bool zero_sets_meet(Context<T>& ctx, const Problem<T, 3>& p, const std::vector<int>& curved)
{
    constexpr int n = 8;
    constexpr int nodes = n + 1;
    std::vector<std::uint64_t> changes(n * n * n, 0);
    std::vector<T> values(nodes * nodes * nodes);
    for (const int c : curved)
    {
        const Func<T, 3>& f = p.funcs[static_cast<std::size_t>(c)];
        const std::uint64_t bit = std::uint64_t(1) << f.ls;
        for (int i = 0; i < nodes; ++i)
            for (int j = 0; j < nodes; ++j)
                for (int l = 0; l < nodes; ++l)
                {
                    const VecD<T, 3> y = {p.lo[0] + (p.hi[0] - p.lo[0]) * T(i) / T(n),
                                          p.lo[1] + (p.hi[1] - p.lo[1]) * T(j) / T(n),
                                          p.lo[2] + (p.hi[2] - p.lo[2]) * T(l) / T(n)};
                    values[(i * nodes + j) * nodes + l] = value<T, 3>(ctx, f, y);
                }
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int l = 0; l < n; ++l)
                {
                    bool negative = false, positive = false;
                    for (int corner = 0; corner < 8; ++corner)
                    {
                        const T v = values[((i + (corner & 1)) * nodes + j + ((corner >> 1) & 1)) * nodes + l
                                           + ((corner >> 2) & 1)];
                        negative = negative || v < T(0);
                        positive = positive || v > T(0);
                    }
                    if (negative && positive)
                        changes[(i * n + j) * n + l] |= bit;
                }
    }
    return std::any_of(changes.begin(), changes.end(), [](std::uint64_t m) { return std::popcount(m) >= 2; });
}

/// A rotated frame for a box where no axis suits all curved functions, as
/// where the zero sets of two level sets meet: the height direction among the
/// unit normals of the functions (curved ones and surface functions) at the
/// box centre and the bisectors of pairs of them that makes the smallest
/// |n_i . d| largest. False if none reaches the margin there. The curved
/// functions are those of @p p that @p an lists.
template <std::floating_point T, int D>
bool choose_frame(Context<T>& ctx, const Problem<T, D>& p, const Analysis<T, D>& an,
                  std::array<VecD<T, D>, D>& R)
{
    VecD<T, D> centre;
    for (int j = 0; j < D; ++j)
        centre[j] = T(0.5) * (p.lo[j] + p.hi[j]);
    std::vector<VecD<T, D>> normals;
    auto add = [&normals](VecD<T, D> g)
    {
        const T norm = scaled_norm(std::span<const T>(g));
        if (!(norm > T(0)) || !std::isfinite(norm))
            return;
        for (int j = 0; j < D; ++j)
            g[j] /= norm;
        normals.push_back(g);
    };
    for (const int c : an.curved)
    {
        const Func<T, D>& f = p.funcs[static_cast<std::size_t>(c)];
        VecD<T, D> g;
        for (int j = 0; j < D; ++j)
            g[j] = derivative_along<T, D>(ctx, f, centre, j);
        add(g);
    }
    if constexpr (D <= 2)
    {
        for (const Surface<T, D>& s : p.surfaces)
        {
            VecD<T, D> g;
            surface_gradient<T, D>(ctx, s, centre, g);
            add(g);
        }
    }
    if (normals.size() < 2)
        return false;
    std::vector<VecD<T, D>> candidates = normals;
    for (std::size_t a = 0; a < normals.size(); ++a)
        for (std::size_t b = a + 1; b < normals.size(); ++b)
            for (const T sign : {T(1), T(-1)})
            {
                VecD<T, D> d;
                for (int j = 0; j < D; ++j)
                    d[j] = normals[a][j] + sign * normals[b][j];
                const T norm = scaled_norm(std::span<const T>(d));
                if (!(norm > T(0.1)))
                    continue;
                for (int j = 0; j < D; ++j)
                    d[j] /= norm;
                candidates.push_back(d);
            }
    T best = -1;
    VecD<T, D> dir{};
    for (const VecD<T, D>& d : candidates)
    {
        T score = T(1);
        for (const VecD<T, D>& n : normals)
        {
            T dot = T(0);
            for (int j = 0; j < D; ++j)
                dot += n[j] * d[j];
            score = std::min(score, std::abs(dot));
        }
        if (score > best)
        {
            best = score;
            dir = d;
        }
    }
    if (!(best >= level_margin<T, D>(ctx)))
        return false;
    // an orthonormal frame with dir as its first axis (columns of R)
    for (int j = 0; j < D; ++j)
        R[j][0] = dir[j];
    if constexpr (D == 2)
    {
        R[0][1] = -dir[1];
        R[1][1] = dir[0];
    }
    else if constexpr (D == 3)
    {
        int e = 0;
        for (int j = 1; j < 3; ++j)
            if (std::abs(dir[j]) < std::abs(dir[e]))
                e = j;
        VecD<T, 3> a{};
        a[e] = T(1);
        T dot = dir[e];
        T norm = T(0);
        for (int j = 0; j < 3; ++j)
        {
            a[j] -= dot * dir[j];
            norm += a[j] * a[j];
        }
        norm = std::sqrt(norm);
        for (int j = 0; j < 3; ++j)
            R[j][1] = a[j] / norm;
        const VecD<T, 3> b = {dir[1] * R[2][1] - dir[2] * R[1][1], dir[2] * R[0][1] - dir[0] * R[2][1],
                              dir[0] * R[1][1] - dir[1] * R[0][1]};
        for (int j = 0; j < 3; ++j)
            R[j][2] = b[j];
    }
    else
        return false;
    return true;
}

// ============================================================================
// Integration by dimension reduction
// ============================================================================

/// Nodes on a segment [a, b]: Gauss-Legendre points and weights for quadrature,
/// or vis_order + 1 equispaced points pulled 1e-5 inside the ends for leaf cells
/// (closer to a breakpoint, a root may fall on either side of a bound by rounding).
template <std::floating_point T, typename F>
void segment_nodes(const Context<T>& ctx, T a, T b, F&& f)
{
    if (ctx.vis_order > 0)
    {
        const T eps = T(1e-5);
        for (int j = 0; j <= ctx.vis_order; ++j)
            f(j, a + (b - a) * (eps + (T(1) - 2 * eps) * T(j) / T(ctx.vis_order)), T(0));
        return;
    }
    for (int j = 0; j < ctx.q; ++j)
        f(j, a + (b - a) * ctx.gauss_x[j], (b - a) * ctx.gauss_w[j]);
}

/// The parameter of the Bernstein ellipse through a point at distance d (in
/// lengths of a segment) beyond an end of the segment: m Gauss-Legendre points
/// on it err by about rho^(-2m) at a branch point there.
template <std::floating_point T>
T ellipse_parameter(T d)
{
    const T u = 1 + 2 * d;
    return u + std::sqrt(u * u - 1);
}

/// Nodes on a base segment [a, b] with the folds of the level above
/// (fold_anchors): where one lies at or beyond an end, within fold_reach
/// segment lengths, the anchored rule, Gauss-Legendre points in s with
/// x = x_f -+ s^2, which makes the square-root branch point there analytic.
/// With folds beyond both ends, the rule anchored at one of them if its
/// n q points in s converge at the other as fast as q points do at
/// fold_reach lengths (rho^n >= ellipse_parameter(fold_reach), rho the other
/// fold's in s); otherwise x = c + r (3t - t^3) / 2, whose derivative
/// vanishes at both, x = c -+ r, with ceil(3 n q / 2) points in t, which
/// integrate polynomials of the same degree n q - 1 in x exactly. Otherwise
/// segment_nodes.
template <std::floating_point T, typename F>
void anchored_nodes(const Context<T>& ctx, T a, T b, const std::vector<T>& folds, F&& f)
{
    const T len = b - a, tol = T(1e-12) * len;
    T left = -infinity<T>, right = infinity<T>;
    for (const T x : folds)
    {
        if (x >= b - tol && x - b <= fold_reach<T> * len)
            right = std::min(right, x);
        if (x <= a + tol && a - x <= fold_reach<T> * len)
            left = std::max(left, x);
    }
    const bool has_left = left > -infinity<T>, has_right = right < infinity<T>;
    if (!has_left && !has_right)
    {
        segment_nodes(ctx, a, b, f);
        return;
    }
    ctx.anchored_any = true;
    // x = x_f - s^2 (a fold above b) or x_f + s^2 (below a) on [a, b]
    const auto anchored = [&](T fold)
    {
        const bool above = fold >= b - tol;
        const T s0 = std::sqrt(std::max(above ? fold - b : a - fold, T(0)));
        const T s1 = std::sqrt(std::max(above ? fold - a : b - fold, T(0)));
        for (std::size_t i = 0; i < ctx.fold_x.size(); ++i)
        {
            const T s = s0 + (s1 - s0) * ctx.fold_x[i];
            f(static_cast<int>(i), above ? fold - s * s : fold + s * s, (s1 - s0) * ctx.fold_w[i] * 2 * s);
        }
    };
    if (!has_left || !has_right)
    {
        anchored(has_right ? right : left);
        return;
    }
    // the other fold beyond the far end of the segment in s, in its lengths,
    // anchored at the right fold (s^2 = right - x) or the left one
    const T width = right - left, sw = std::sqrt(width);
    const T sa = std::sqrt(std::max(right - a, T(0))), sb = std::sqrt(std::max(right - b, T(0)));
    const T ta = std::sqrt(std::max(a - left, T(0))), tb = std::sqrt(std::max(b - left, T(0)));
    const T other_left = std::max(a - left, T(0)) * (sa + sb) / (len * (sw + sa));
    const T other_right = std::max(right - b, T(0)) * (ta + tb) / (len * (sw + tb));
    if (std::pow(ellipse_parameter(std::max(other_left, other_right)), T(ctx.fold_x.size() / ctx.gauss_x.size()))
        >= ellipse_parameter(fold_reach<T>))
    {
        anchored(other_left >= other_right ? right : left);
        return;
    }
    // x = c + r phi(t), phi(t) = (3t - t^3) / 2 = sin(3 theta) for t = 2 sin(theta).
    // Inverted from the distance to the nearer fold: asin((x - c) / r) would turn
    // rounding at a fold into sqrt(eps), where the map's derivative vanishes, and
    // a surface's weight, 1 / sqrt(x - x_f), is large enough to lose that sliver.
    const T c = T(0.5) * (left + right), r = T(0.5) * width;
    const auto t_of = [&](T x)
    {
        const T from_left = std::max(x - left, T(0)) / width, from_right = std::max(right - x, T(0)) / width;
        const T half_pi = T(1.57079632679489661923);
        const T theta = from_left <= from_right ? 2 * std::asin(std::sqrt(from_left)) - half_pi
                                                : half_pi - 2 * std::asin(std::sqrt(from_right));
        return 2 * std::sin(theta / 3);
    };
    const T t0 = t_of(a), t1 = t_of(b);
    for (std::size_t i = 0; i < ctx.cubic_x.size(); ++i)
    {
        const T t = t0 + (t1 - t0) * ctx.cubic_x[i];
        f(static_cast<int>(i), c + r * T(0.5) * (3 * t - t * t * t),
          (t1 - t0) * ctx.cubic_w[i] * r * T(1.5) * (1 - t * t));
    }
}

/// Roots in (L, U) of a curved function on the height line y = insert(yb, k, t),
/// appended to @p roots: on a certified box the one root the ends bracket, or
/// with two roots per line the roots on either side of the extreme point;
/// otherwise all roots, isolated. Bernstein forms leave the line's form in
/// @p form.
template <std::floating_point T, int D>
void curved_roots(Context<T>& ctx, const Func<T, D>& f, const VecD<T, D - 1>& yb, int k, T L, T U, bool certified,
                  bool two, BoxBernstein<T>& form, std::vector<T>& roots)
{
    const Source<T>& phi = ctx.phis[static_cast<std::size_t>(f.ls)];
    const VecD<T, D> yL = insert<T, D>(yb, k, L);
    if (phi.is_analytic())
    {
        Vec3<T> u0, dir;
        box_line<T, D>(f, insert<T, D>(yb, k, T(0)), k, u0, dir);
        if (!certified)
        {
            line_roots(phi, std::span<const T>(u0), std::span<const T>(dir), L, U,
                       zero_function_tol<T> * ctx.scales[static_cast<std::size_t>(f.ls)], roots);
            return;
        }
        const T gl = value<T, D>(ctx, f, yL), gu = value<T, D>(ctx, f, insert<T, D>(yb, k, U));
        T c = U, gc = gu;
        if (two)
        {
            // the extreme point: the root of d_k psi, which is monotone
            auto slope = [&](T t) { return derivative_along<T, D>(ctx, f, insert<T, D>(yb, k, t), k); };
            const T dl = slope(L), du = slope(U);
            if (opposite(dl, du))
            {
                c = illinois_root(slope, L, U, dl, du);
                gc = value<T, D>(ctx, f, insert<T, D>(yb, k, c));
                if (opposite(gl, gc))
                    roots.push_back(line_root(phi, std::span<const T>(u0), std::span<const T>(dir), L, c, gl, gc));
                if (opposite(gc, gu))
                    roots.push_back(line_root(phi, std::span<const T>(u0), std::span<const T>(dir), c, U, gc, gu));
                return;
            }
        }
        if (opposite(gl, gu))
            roots.push_back(line_root(phi, std::span<const T>(u0), std::span<const T>(dir), L, U, gl, gu));
        return;
    }
    line_form<T, D>(ctx, f, yL, k, L, U, form, ctx.line_work[D]);
    const std::span<const T> c(form.coeffs);
    if (!certified)
    {
        isolate_roots(c, L, U, roots, ctx.root_work[D]);
        return;
    }
    const T gl = c.front(), gu = c.back();
    const auto solve = [&ctx](std::span<const T> poly, T a, T b, T lo, T hi, T glo, T ghi)
    {
        return ctx.vis_order == 0 ? newton_root(poly, a, b, lo, hi, glo, ghi)
                                  : bracketed_root(poly, a, b, lo, hi, glo, ghi);
    };
    if (two && c.size() > 2)
    {
        // the extreme point: the root of the line's derivative, which is monotone
        std::vector<T>& dc = ctx.slope_line[D];
        const int n = static_cast<int>(c.size()) - 1;
        dc.resize(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i)
            dc[static_cast<std::size_t>(i)] = T(n) * (c[static_cast<std::size_t>(i + 1)] - c[static_cast<std::size_t>(i)]);
        const T dl = dc.front(), du = dc.back();
        if (opposite(dl, du))
        {
            const T x = solve(std::span<const T>(dc), L, U, L, U, dl, du);
            const T gx = evaluate_1d(c, (x - L) / (U - L));
            if (opposite(gl, gx))
                roots.push_back(solve(c, L, U, L, x, gl, gx));
            if (opposite(gx, gu))
                roots.push_back(solve(c, L, U, x, U, gx, gu));
            return;
        }
    }
    if (opposite(gl, gu))
        roots.push_back(solve(c, L, U, L, U, gl, gu));
}

/// The gradient and Hessian of f at y (forward differences of the gradient,
/// step eps), divided by the gradient's largest component, so that level sets
/// scaled by 1e+-200 give the same folds: that component, or 0 if singular.
template <std::floating_point T, typename Gradient>
T fold_local(Gradient&& gradient_at, const VecD<T, 2>& y, T eps, VecD<T, 2>& g, std::array<VecD<T, 2>, 2>& h)
{
    gradient_at(y, g);
    const T scale = std::max(std::abs(g[0]), std::abs(g[1]));
    if (!(scale > T(0)) || !(scale < infinity<T>))
        return T(0);
    for (int d = 0; d < 2; ++d)
    {
        VecD<T, 2> yp = y, gp;
        yp[d] += eps;
        gradient_at(yp, gp);
        for (int e = 0; e < 2; ++e)
            h[d][e] = (gp[e] - g[e]) / (eps * scale);
    }
    h[0][1] = h[1][0] = T(0.5) * (h[0][1] + h[1][0]);
    g[0] /= scale;
    g[1] /= scale;
    return scale;
}

/// Newton's method on (f, d_k f) in (y_j, y_k), j = 1 - k, from y: true if its
/// steps fall below 1e-12 @p size within 12 iterations, y then the fold.
template <std::floating_point T, typename Value, typename Gradient>
bool fold_newton(Value&& value_at, Gradient&& gradient_at, int k, T size, VecD<T, 2>& y)
{
    const int j = 1 - k;
    VecD<T, 2> g;
    std::array<VecD<T, 2>, 2> h;
    for (int it = 0; it < 12; ++it)
    {
        const T scale = fold_local<T>(gradient_at, y, T(1e-5) * size, g, h);
        if (scale == T(0))
            return false;
        const T r1 = value_at(y) / scale, r2 = g[k];
        const T det = g[j] * h[k][k] - g[k] * h[k][j];
        if (!(std::abs(det) > T(0)))
            return false;
        const T dj = (r1 * h[k][k] - g[k] * r2) / det, dk = (g[j] * r2 - h[k][j] * r1) / det;
        y[j] -= dj;
        y[k] -= dk;
        if (!(std::abs(y[j]) + std::abs(y[k]) < infinity<T>))
            return false;
        if (std::abs(dj) + std::abs(dk) <= T(1e-12) * size)
            return true;
    }
    return false;
}

/// Sort folds and merge those closer than 1e-10 @p size.
template <std::floating_point T>
void unique_folds(std::vector<T>& folds, T size)
{
    std::sort(folds.begin(), folds.end());
    folds.erase(std::unique(folds.begin(), folds.end(), [&](T a, T b) { return b - a <= T(1e-10) * size; }),
                folds.end());
}

/// From where a zero curve leaves a level-2 box (@p starts), the folds along k it
/// reaches: Newton's method on (f, d_k f) from the points of the osculating circle
/// there whose tangent runs along k. Their base coordinates within fold_reach times
/// the box's extent of it are appended; @p found, if given, receives all it reaches.
template <std::floating_point T, typename Value, typename Gradient>
void osculating_folds(Value&& value_at, Gradient&& gradient_at, const std::vector<VecD<T, 2>>& starts,
                      const VecD<T, 2>& lo, const VecD<T, 2>& hi, int k, std::vector<T>& folds,
                      std::vector<VecD<T, 2>>* found)
{
    const int j = 1 - k;
    const T extent = hi[j] - lo[j], size = std::max(extent, hi[k] - lo[k]);
    const T reach = fold_reach<T> * extent;
    for (const VecD<T, 2>& start : starts)
    {
        VecD<T, 2> g;
        std::array<VecD<T, 2>, 2> h;
        if (fold_local<T>(gradient_at, start, T(1e-5) * size, g, h) == T(0))
            continue;
        const T gn = std::sqrt(g[0] * g[0] + g[1] * g[1]);
        // signed curvature div(grad f / |grad f|): the centre is at -n / kappa
        const T kappa = (h[j][j] * g[k] * g[k] - 2 * h[j][k] * g[j] * g[k] + h[k][k] * g[j] * g[j]) / (gn * gn * gn);
        if (!(std::abs(kappa) * size > T(1e-6)))
            continue; // straight there: no fold near
        const T radius = T(1) / std::abs(kappa);
        VecD<T, 2> centre;
        for (int d = 0; d < 2; ++d)
            centre[d] = start[d] - g[d] / (gn * kappa);
        for (const T side : {T(-1), T(1)})
        {
            VecD<T, 2> y = centre;
            y[j] += side * radius;
            if (y[j] < lo[j] - reach - extent || y[j] > hi[j] + reach + extent)
                continue;
            if (!fold_newton<T>(value_at, gradient_at, k, size, y))
                continue;
            if (found != nullptr)
                found->push_back(y);
            if (y[j] >= lo[j] - reach && y[j] <= hi[j] + reach)
                folds.push_back(y[j]);
        }
    }
}

/// The folds near a certified level-2 box with height direction k: points where
/// the zero set of a curved function turns back along k (f = d_k f = 0). Its
/// height function has a square-root branch point at the fold's coordinate on
/// the base line, which slows Gauss-Legendre points on the base segments that
/// end near it. Each zero set is followed from where it leaves the box: Newton's
/// method on (f, d_k f) starts at the points of the osculating circle there
/// whose tangent runs along k. The base coordinates of the folds within
/// fold_reach times the box's extent of it along the base line are appended
/// (Options::anchor_folds). With @p monotone, the functions are monotone along
/// k on the box, so that the ends of a side along k bracket its one root; @p
/// found, if given, receives every fold Newton's method reaches, near or not.
template <std::floating_point T>
void fold_anchors(Context<T>& ctx, const Problem<T, 2>& p, int k, std::vector<T>& folds, bool monotone = false,
                  std::vector<VecD<T, 2>>* found = nullptr)
{
    const int j = 1 - k;
    const T extent = p.hi[j] - p.lo[j], size = std::max(extent, p.hi[k] - p.lo[k]);
    if (!(size > T(0)))
        return;
    const T reach = fold_reach<T> * extent;
    std::vector<T>& exits = ctx.nodes[2];
    std::vector<VecD<T, 2>>& starts = ctx.fold_starts;
    for (const Func<T, 2>& f : p.funcs)
    {
        if (f.linear || ctx.flat[static_cast<std::size_t>(f.ls)])
            continue;
        // where the zero set leaves the box: its roots on the four sides
        starts.clear();
        for (int side = 0; side < 2; ++side)
            for (const int along : {j, k})
            {
                const int across = 1 - along;
                const T c = side == 0 ? p.lo[across] : p.hi[across];
                exits.clear();
                curved_roots<T, 2>(ctx, f, VecD<T, 1>{c}, along, p.lo[along], p.hi[along], monotone && along == k,
                                   false, ctx.line[2], exits);
                for (const T t : exits)
                    starts.push_back(insert<T, 2>(VecD<T, 1>{c}, along, t));
            }
        osculating_folds<T>([&](const VecD<T, 2>& y) { return value<T, 2>(ctx, f, y); },
                            [&](const VecD<T, 2>& y, VecD<T, 2>& g) { func_gradient<T, 2>(ctx, f, y, g); }, starts,
                            p.lo, p.hi, k, folds, found);
    }
    unique_folds(folds, size);
}

/// fold_anchors for a surface function in a level-2 box, a silhouette or a meeting
/// curve: its folds just beyond the box, which the box's split (fold_split_points)
/// does not reach. There its height lines have a square-root branch point, and so
/// has the integrand of the base line.
template <std::floating_point T>
void silhouette_anchors(Context<T>& ctx, const Surface<T, 2>& s, const Problem<T, 2>& p, int k, std::vector<T>& folds)
{
    const int j = 1 - k;
    std::vector<T>& exits = ctx.nodes[2];
    std::vector<VecD<T, 2>>& starts = ctx.fold_starts;
    starts.clear();
    for (int side = 0; side < 2; ++side)
        for (const int along : {j, k})
        {
            const int across = 1 - along;
            const VecD<T, 1> c = {side == 0 ? p.lo[across] : p.hi[across]};
            exits.clear();
            surface_line_roots<T>(ctx, surface_on_line<T, 2>(s, c, along), p.lo[along], p.hi[along], exits);
            for (const T t : exits)
                starts.push_back(insert<T, 2>(c, along, t));
        }
    osculating_folds<T>([&](const VecD<T, 2>& y) { return surface_value<T, 2>(ctx, s, y); },
                        [&](const VecD<T, 2>& y, VecD<T, 2>& g) { surface_gradient<T, 2>(ctx, s, y, g); }, starts, p.lo,
                        p.hi, k, folds, nullptr);
}

/// The folds of the level above a level-2 box on its height line y =
/// insert(yb, k, tau), tau in [L, U]: where the zero set of a height function
/// (Problem::heights), in the plane through the line and the height direction
/// above, turns back along that direction, beyond the box, where |d_t f| no
/// longer has the margin (the silhouette of the level above, seen along it).
/// The integrand along the line has a square-root branch point there; the
/// coordinates tau of those near the line (fold_anchors on the plane) are
/// appended. The plane lies in the certified box above, where a function with
/// one root per height line is monotone along t. On the first pass of an
/// accuracy check the folds are recorded (Context::fold_pass).
template <std::floating_point T>
void height_folds(Context<T>& ctx, const std::vector<Height<T, 2>>& heights, const VecD<T, 1>& yb, int k, T L, T U,
                  std::vector<T>& folds)
{
    Problem<T, 2>& plane = ctx.fold_plane;
    plane.funcs.resize(1);
    const bool record = ctx.fold_pass == 1;
    for (std::size_t i = 0; i < heights.size(); ++i)
    {
        const Height<T, 2>& h = heights[i];
        plane.lo = {L, h.t_lo};
        plane.hi = {U, h.t_hi};
        plane.funcs[0] = line_plane<T, 2>(h.f, yb, k, h.k);
        ctx.fold_found.clear();
        fold_anchors<T>(ctx, plane, 1, folds, !h.two, record ? &ctx.fold_found : nullptr);
        for (const VecD<T, 2>& y : ctx.fold_found)
            ctx.column_fold.push_back({static_cast<int>(i), y});
    }
    if (record)
    {
        ctx.column_y.push_back(yb[0]);
        ctx.column_start.push_back(static_cast<int>(ctx.column_fold.size()));
    }
}

/// height_folds on the second pass of an accuracy check: Newton's method from
/// the folds the first pass found on the nearest height lines on either side,
/// whose rule had other points, instead of the search. A fold no line of the
/// first pass came near is missed, which leaves the second rule less accurate
/// and so the check stricter, never laxer.
template <std::floating_point T>
void continued_height_folds(Context<T>& ctx, const std::vector<Height<T, 2>>& heights, const VecD<T, 1>& yb, int k,
                            T L, T U, std::vector<T>& folds)
{
    const std::vector<int>& order = ctx.column_order;
    const auto above = std::lower_bound(order.begin(), order.end(), yb[0],
                                        [&](int c, T y) { return ctx.column_y[static_cast<std::size_t>(c)] < y; });
    const T extent = U - L, reach = fold_reach<T> * extent;
    T size = extent;
    for (const Height<T, 2>& h : heights)
        size = std::max(size, h.t_hi - h.t_lo);
    for (auto it = above == order.begin() ? above : above - 1; it != order.end() && it <= above; ++it)
    {
        const std::size_t c = static_cast<std::size_t>(*it);
        for (int m = ctx.column_start[c]; m < ctx.column_start[c + 1]; ++m)
        {
            const ColumnFold<T>& fold = ctx.column_fold[static_cast<std::size_t>(m)];
            const Height<T, 2>& h = heights[static_cast<std::size_t>(fold.height)];
            VecD<T, 2> y = fold.y;
            const Func<T, 2> f = line_plane<T, 2>(h.f, yb, k, h.k);
            if (fold_newton<T>([&](const VecD<T, 2>& z) { return value<T, 2>(ctx, f, z); },
                               [&](const VecD<T, 2>& z, VecD<T, 2>& g) { func_gradient<T, 2>(ctx, f, z, g); }, 1,
                               std::max(extent, h.t_hi - h.t_lo), y)
                && y[0] >= L - reach && y[0] <= U + reach)
                folds.push_back(y[0]);
        }
    }
    unique_folds(folds, size);
}

/// The folds of the level above a level-2 box where its silhouette (see
/// height_folds) crosses a bound of the height lines, y_k = alpha + beta x:
/// there the branch point of the lines' integrand reaches their end, and the
/// integrand of the base line has one too. The base coordinates of those near
/// the base line [lo, hi] (fold_anchors on the plane through the bound and
/// the height direction above) are appended.
template <std::floating_point T>
void bound_folds(Context<T>& ctx, const std::vector<Height<T, 2>>& heights, int k,
                 const std::vector<Bound<T, 2>>& bounds, T lo, T hi, std::vector<T>& folds)
{
    Problem<T, 2>& plane = ctx.fold_plane;
    plane.funcs.resize(1);
    for (const Height<T, 2>& h : heights)
        for (const Bound<T, 2>& b : bounds)
        {
            const Height<T, 1> r = restrict_height<T, 2>(h, k, b);
            plane.lo[1 - r.k] = lo;
            plane.hi[1 - r.k] = hi;
            plane.lo[r.k] = r.t_lo;
            plane.hi[r.k] = r.t_hi;
            plane.funcs[0] = r.f;
            fold_anchors<T>(ctx, plane, r.k, folds);
        }
}

/// Does a top-level point satisfy one of the terms left on its box? Values of
/// the level sets come from their forms on the current line where they exist.
template <std::floating_point T>
bool holds(Context<T>& ctx, const std::vector<Requirement>& terms, const Vec3<T>& u, T s)
{
    std::uint64_t evaluated = 0;
    for (const Requirement& r : terms)
    {
        bool ok = true;
        for (std::uint64_t bits = r.negative | r.positive; bits != 0 && ok; bits &= bits - 1)
        {
            const int l = std::countr_zero(bits);
            const std::uint64_t bit = std::uint64_t(1) << l;
            T& v = ctx.top_values[static_cast<std::size_t>(l)];
            if (!(evaluated & bit))
            {
                evaluated |= bit;
                if (ctx.top_line_set[static_cast<std::size_t>(l)])
                    v = evaluate_1d(std::span<const T>(ctx.top_lines[static_cast<std::size_t>(l)].coeffs), s);
                else
                    v = evaluate(ctx.phis[static_cast<std::size_t>(l)], std::span<const T>(u));
            }
            ok = (r.negative & bit) ? v < T(0) : v > T(0);
        }
        if (ok)
            return true;
    }
    return false;
}

/// The points where a curve part's two level sets cross, on a 2D cell: the
/// roots of their surface function on the base line [L, U] of a box, each with
/// weight 1. An end of the line counts within the Newton distance point_tol
/// (a point on the face two boxes share is found by both, and merged once the
/// cell is done); where the two zero sets coincide there are none.
template <std::floating_point T, typename Emit>
void crossing_points(Context<T>& ctx, const Problem<T, 1>& p, T L, T U, const Emit& emit)
{
    std::vector<T>& nodes = ctx.nodes[1];
    NodeTag tag;
    tag.box[0] = ctx.box_ids[0]++;
    int index = 0;
    for (const Surface<T, 1>& s : p.surfaces)
    {
        if (!partner_surface<T, 1>(ctx, s) || surface_vanishes_identically<T, 1>(ctx, s, {L}, {U}))
            continue;
        nodes.clear();
        surface_line_roots<T>(ctx, s, L, U, nodes);
        for (const T end : {L, U})
        {
            VecD<T, 1> g;
            surface_gradient<T, 1>(ctx, s, {end}, g);
            if (std::abs(surface_value<T, 1>(ctx, s, {end})) <= point_tol<T> * std::abs(g[0]))
                nodes.push_back(end);
        }
        std::sort(nodes.begin(), nodes.end());
        for (const T t : nodes)
        {
            tag.segment[0] = index++;
            tag.node[0] = 0;
            emit(VecD<T, 1>{t}, T(1), tag, VecD<T, 1>{T(0)});
        }
    }
}

/// One free coordinate: Gauss-Legendre on the segments between bounds and
/// roots, anchored at the given folds of the level above (anchored_nodes); the
/// derivative of the points along their coordinate is 1. Below the top of a 2D
/// cell (Points), a curve part's crossing points instead.
template <std::floating_point T, bool Points, typename Emit>
void integrate_line(Context<T>& ctx, const Problem<T, 1>& p, const Emit& emit, const std::vector<T>* folds = nullptr)
{
    T L = p.lo[0], U = p.hi[0];
    for (const Half<T, 1>& h : p.clips)
    {
        if (h.c[0] > tiny<T>)
            U = std::min(U, h.d / h.c[0]);
        else if (h.c[0] < -tiny<T>)
            L = std::max(L, h.d / h.c[0]);
        else if (h.d < T(0))
            return;
    }
    if (!(U > L))
        return;
    if constexpr (Points)
    {
        if (ctx.partner >= 0)
        {
            crossing_points<T>(ctx, p, L, U, emit);
            return;
        }
    }
    std::vector<T>& nodes = ctx.nodes[1];
    nodes.clear();
    nodes.push_back(L);
    nodes.push_back(U);
    for (const Func<T, 1>& f : p.funcs)
    {
        if (f.linear)
        {
            if (std::abs(f.a[0]) > tiny<T>)
            {
                const T t = -f.c / f.a[0];
                if (t > L && t < U)
                    nodes.push_back(t);
            }
            continue;
        }
        if (is_analytic(ctx, f.ls))
        {
            Vec3<T> u0, dir;
            box_line<T, 1>(f, VecD<T, 1>{T(0)}, 0, u0, dir);
            line_roots(ctx.phis[static_cast<std::size_t>(f.ls)], std::span<const T>(u0), std::span<const T>(dir), L,
                       U, zero_function_tol<T> * ctx.scales[static_cast<std::size_t>(f.ls)], nodes);
            continue;
        }
        line_form<T, 1>(ctx, f, VecD<T, 1>{L}, 0, L, U, ctx.line[1], ctx.line_work[1]);
        isolate_roots(std::span<const T>(ctx.line[1].coeffs), L, U, nodes, ctx.root_work[1]);
    }
    // Where a silhouette crosses a bound of the height lines above, the integral
    // along them has a square root in this coordinate: bound_folds does not see it
    // where the zero set in the plane through the bound never leaves the box.
    std::vector<T>& anchors = ctx.line_folds[1];
    anchors.clear();
    if (folds != nullptr)
        anchors = *folds;
    for (const Surface<T, 1>& s : p.surfaces)
        if (!surface_vanishes_identically<T, 1>(ctx, s, {L}, {U}))
        {
            const std::size_t found = nodes.size();
            surface_line_roots<T>(ctx, s, L, U, nodes);
            if (s.silhouette)
                anchors.insert(anchors.end(), nodes.begin() + static_cast<std::ptrdiff_t>(found), nodes.end());
        }
    std::sort(nodes.begin(), nodes.end());
    NodeTag tag;
    tag.box[0] = ctx.box_ids[0]++;
    int segment = 0;
    for (std::size_t s = 0; s + 1 < nodes.size(); ++s)
    {
        const T a = nodes[s], b = nodes[s + 1];
        if (b - a <= segment_tol<T>)
            continue;
        tag.segment[0] = segment++;
        const auto node = [&](int j, T t, T w)
        {
            tag.node[0] = j;
            emit(VecD<T, 1>{t}, w, tag, VecD<T, 1>{T(1)});
        };
        if (!anchors.empty() && ctx.vis_order == 0)
            anchored_nodes(ctx, a, b, anchors, node);
        else
            segment_nodes(ctx, a, b, node);
    }
}

/// Where the roots of two surface functions on the height lines of a level-2
/// box cross: three level sets meet there (a corner), and the integrand of the
/// line below has a kink. Both are monotone along k on the certified box. A
/// root exists where its function differs in sign at the two ends of the line,
/// so between the points where either function vanishes on an end, both roots
/// exist throughout or nowhere; the difference of the roots is sampled on each
/// such interval and its sign changes are refined. The base coordinates of the
/// crossings are appended.
template <std::floating_point T>
void surface_crossings(Context<T>& ctx, const Surface<T, 2>& a, const Surface<T, 2>& b, int k, T lo, T hi, T blo,
                       T bhi, std::vector<T>& out)
{
    auto root = [&](const Surface<T, 2>& s, T y) -> T
    {
        const T fl = surface_value<T, 2>(ctx, s, insert<T, 2>(VecD<T, 1>{y}, k, lo)),
                fh = surface_value<T, 2>(ctx, s, insert<T, 2>(VecD<T, 1>{y}, k, hi));
        if (!opposite(fl, fh))
            return std::numeric_limits<T>::quiet_NaN();
        return illinois_root([&](T t) { return surface_value<T, 2>(ctx, s, insert<T, 2>(VecD<T, 1>{y}, k, t)); }, lo, hi,
                             fl, fh);
    };
    auto gap = [&](T y) { return root(a, y) - root(b, y); };

    // where a root enters or leaves the line through one of its ends
    std::vector<T> ends = {blo, bhi};
    for (const Surface<T, 2>* s : {&a, &b})
        for (int side = 0; side < 2; ++side)
        {
            Bound<T, 2> end;
            end.alpha = side == 0 ? lo : hi;
            end.kind = side == 0 ? 2 : 3;
            surface_line_roots<T>(ctx, restrict_surface<T, 2>(*s, k, end), blo, bhi, ends);
        }
    std::sort(ends.begin(), ends.end());

    constexpr int samples = 9;
    const T tol = segment_tol<T> * (bhi - blo);
    for (std::size_t e = 0; e + 1 < ends.size(); ++e)
    {
        const T e0 = ends[e], e1 = ends[e + 1];
        if (!(e1 - e0 > tol))
            continue;
        // the first and last samples just inside, where a root may sit on an end
        auto sample = [&](int i)
        {
            const T t = std::clamp(T(i) / T(samples - 1), T(1e-6), T(1) - T(1e-6));
            return e0 + (e1 - e0) * t;
        };
        T y0 = sample(0), g0 = gap(y0);
        for (int i = 1; i < samples; ++i)
        {
            const T y1 = sample(i), g1 = gap(y1);
            if (std::isfinite(g0) && std::isfinite(g1) && opposite(g0, g1))
            {
                // the roots may still vanish inside: bisect while both exist
                T l = y0, h = y1, gl = g0;
                for (int it = 0; it < 60 && h - l > tol; ++it)
                {
                    const T m = T(0.5) * (l + h), gm = gap(m);
                    if (!std::isfinite(gm))
                        break;
                    if (opposite(gl, gm))
                        h = m;
                    else
                    {
                        l = m;
                        gl = gm;
                    }
                }
                out.push_back(T(0.5) * (l + h));
            }
            y0 = y1;
            g0 = g1;
        }
    }
}

template <std::floating_point T, int D, bool Rotated, int Top, typename Emit>
void integrate_box(Context<T>& ctx, Problem<T, D> p, const Emit& emit, int depth, bool interface);

/// True if the interface's level set vanishes on the whole plane y_axis = mid
/// of the box: halves split there would share its zero set as a face, on which
/// neither has a root.
template <std::floating_point T, int D>
bool interface_on_plane(Context<T>& ctx, const Problem<T, D>& p, int axis, T mid)
{
    if constexpr (D < 2)
        return false;
    else
    {
        thread_local BoxBernstein<T> form;
        thread_local std::vector<T> work;
        for (const Func<T, D>& f : p.funcs)
        {
            if (f.linear || f.ls != ctx.surface)
                continue;
            Bound<T, D> plane;
            plane.alpha = mid;
            const Func<T, D - 1> g = restrict_to<T, D>(f, axis, plane);
            VecD<T, D - 1> lo, hi;
            for (int jb = 0; jb < D - 1; ++jb)
            {
                const int j = jb < axis ? jb : jb + 1;
                lo[jb] = p.lo[j];
                hi[jb] = p.hi[j];
            }
            const T tol = zero_function_tol<T> * ctx.scales[static_cast<std::size_t>(f.ls)];
            if (is_analytic(ctx, f.ls))
            {
                int sign = 0;
                VecD<T, D - 1> ratio{}, upper{};
                T magnitude = T(0);
                analytic_bounds<T, D - 1>(ctx, g, lo, hi, sign, ratio, magnitude, upper);
                return magnitude <= tol;
            }
            bernstein_form<T, D - 1>(ctx, g, lo, hi, form, work);
            return max_abs(std::span<const T>(form.coeffs)) <= tol;
        }
        return false;
    }
}

/// Physical curvature of the zero sets of the curved level sets the part
/// names at box coordinates @p u: the largest |d n / dx| along the cell's
/// axes, n the unit normal, from central differences; the largest over the
/// level sets.
template <std::floating_point T>
T curvature(const Context<T>& ctx, const Vec3<T>& u)
{
    const auto normal = [&ctx](const Source<T>& phi, const Vec3<T>& v, Vec3<T>& n)
    {
        Vec3<T> g{};
        gradient(phi, std::span<const T>(v), std::span<T>(g));
        // level sets scaled by 1e+-200 must give the same normals: no squares
        // before dividing by the largest component
        T largest = max_abs(std::span<const T>(g));
        if (!(largest > T(0)) || !(largest < infinity<T>))
            return false;
        n = {};
        for (int i = 0; i < 3; ++i)
            for (int m = 0; m < 3; ++m)
                n[i] += ctx.inv[m][i] * (g[m] / largest);
        const T norm = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        if (!(norm > T(0)))
            return false;
        for (T& c : n)
            c /= norm;
        return true;
    };
    const T eps = T(1e-4);
    T kappa = 0;
    for (std::size_t l = 0; l < ctx.derived.size() / 3; ++l)
    {
        if (!((ctx.used >> l) & 1) || ctx.flat[l])
            continue;
        for (int j = 0; j < ctx.tdim; ++j)
        {
            Vec3<T> up = u, down = u, np, nm;
            up[j] += eps;
            down[j] -= eps;
            if (!normal(ctx.phis[l], up, np) || !normal(ctx.phis[l], down, nm))
                return infinity<T>; // a singular point: no check
            T dn = 0, dx = 0;
            for (int i = 0; i < 3; ++i)
            {
                dn += (np[i] - nm[i]) * (np[i] - nm[i]);
                dx += ctx.jac[i][j] * ctx.jac[i][j];
            }
            kappa = std::max(kappa, std::sqrt(dn / dx) / (2 * eps));
        }
    }
    return kappa;
}

/// Enters box @p p in the trace (Context::trace) as the one being integrated.
/// Returns its index, or -1 if nothing is recorded.
template <std::floating_point T, int D, int Top>
int trace_enter(Context<T>& ctx, const Problem<T, D>& p, int k, int depth, bool rotated, BoxFate fate)
{
    if (ctx.trace == nullptr)
        return -1;
    BoxTrace<T>& t = *ctx.trace;
    const int id = t.n_boxes();
    t.parent.push_back(ctx.trace_box);
    t.level.push_back(D);
    t.depth.push_back(depth);
    t.height.push_back(k);
    t.fate.push_back(fate);
    t.rotated.push_back(rotated ? 1 : 0);
    for (int j = 0; j < 3; ++j)
    {
        t.lo.push_back(j < D ? p.lo[j] : T(0));
        t.hi.push_back(j < D ? p.hi[j] : T(0));
    }
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            t.frame.push_back(D == Top && j < D ? p.frame.A[i][j] : T(0));
    for (int i = 0; i < 3; ++i)
        t.frame.push_back(D == Top ? p.frame.b[i] : T(0));
    t.points.push_back(ctx.emitted->n_points());
    ctx.trace_box = id;
    return id;
}

/// Leaves box @p id of the trace (trace_enter): its points are those emitted since.
template <std::floating_point T>
void trace_leave(Context<T>& ctx, int id)
{
    if (id < 0)
        return;
    BoxTrace<T>& t = *ctx.trace;
    t.points[static_cast<std::size_t>(id)] = ctx.emitted->n_points() - t.points[static_cast<std::size_t>(id)];
    ctx.trace_box = t.parent[static_cast<std::size_t>(id)];
}

/// Sets what became of the box being integrated.
template <std::floating_point T>
void trace_fate(Context<T>& ctx, BoxFate fate)
{
    if (ctx.trace != nullptr)
        ctx.trace->fate[static_cast<std::size_t>(ctx.trace_box)] = fate;
}

/// Bisect a box along its longest side and integrate the halves.
template <std::floating_point T, int D, bool Rotated, int Top, typename Emit>
void bisect(Context<T>& ctx, Problem<T, D> p, const Emit& emit, int depth, bool interface)
{
    if constexpr (D == 2 && Top == 3)
    {
        const bool for_moment = ctx.moment_only;
        ctx.moment_only = false;
        if (ctx.top_unmargined && !for_moment && ++ctx.level2_bisections > unmargined_bisections)
        {
            ctx.unmargined_failed = true;
            return;
        }
    }
    ++ctx.bisections;
    ++ctx.stats->bisections;
    int axis = 0;
    for (int j = 1; j < D; ++j)
        if (p.hi[j] - p.lo[j] > p.hi[axis] - p.lo[axis])
            axis = j;
    // A curve in the mid-plane of a box (its partner's zero set, say a
    // plane through a cell's centre, or two balls symmetric about it), or
    // whose base lies on the mid-line of a base, would lie on the face the
    // halves share, where neither finds it: curve parts split off-centre.
    T mid = ctx.partner >= 0 ? p.lo[axis] + T(0.53) * (p.hi[axis] - p.lo[axis]) : T(0.5) * (p.lo[axis] + p.hi[axis]);
    if constexpr (D == Top)
        if (interface && interface_on_plane<T, D>(ctx, p, axis, mid))
            mid = p.lo[axis] + T(0.625) * (p.hi[axis] - p.lo[axis]);
    Problem<T, D> first = p, second = std::move(p);
    first.hi[axis] = mid;
    second.lo[axis] = mid;
    integrate_box<T, D, Rotated, Top>(ctx, std::move(first), emit, depth + 1, interface);
    integrate_box<T, D, Rotated, Top>(ctx, std::move(second), emit, depth + 1, interface);
}

/// Where a meeting curve of two level sets touches the silhouette of one of them in
/// a level-2 box: f = g = d_t h = 0, h one of the two. Seen along t both curves pass
/// through the point with one tangent, so that Newton's method on them fails; in 3D
/// the three functions cross transversally. Subdivision with bounds, then Newton's
/// method; the base coordinates of the points inside the box are appended.
template <std::floating_point T>
void contact_points(Context<T>& ctx, const Func<T, 3>& f, const Func<T, 3>& g, const Func<T, 3>& slope, int k, T t_lo,
                    T t_hi, const VecD<T, 2>& lo, const VecD<T, 2>& hi, std::vector<VecD<T, 2>>& out)
{
    constexpr int max_depth = 7, max_cells = 2000;
    struct Cell
    {
        VecD<T, 3> lo, hi;
        int depth;
    };
    thread_local BoxBernstein<T> form;
    thread_local std::vector<T> work;
    const std::array<const Func<T, 3>*, 3> functions = {&f, &g, &slope};
    const auto may_vanish_on = [&](const Func<T, 3>& h, const VecD<T, 3>& a, const VecD<T, 3>& b)
    {
        if (is_analytic(ctx, h.ls))
        {
            const Ranges<T, 3> r = curved_ranges<T, 3>(ctx, h, a, b, ctx.inner_form[3], ctx.inner_deriv[3], ctx.inner_work[3]);
            return !r.valid || r.sign == 0;
        }
        bernstein_form<T, 3>(ctx, h, a, b, form, work);
        return may_vanish(std::span<const T>(form.coeffs));
    };
    std::vector<Cell> stack;
    stack.push_back({insert<T, 3>(lo, k, t_lo), insert<T, 3>(hi, k, t_hi), 0});
    int cells = 0;
    const T size = std::max({hi[0] - lo[0], hi[1] - lo[1], t_hi - t_lo});
    while (!stack.empty() && ++cells <= max_cells)
    {
        const Cell c = stack.back();
        stack.pop_back();
        bool keep = true;
        for (const Func<T, 3>* h : functions)
            keep = keep && may_vanish_on(*h, c.lo, c.hi);
        if (!keep)
            continue;
        if (c.depth < max_depth)
        {
            for (int corner = 0; corner < 8; ++corner)
            {
                Cell child = c;
                child.depth = c.depth + 1;
                for (std::size_t j = 0; j < 3; ++j)
                {
                    const T m = T(0.5) * (c.lo[j] + c.hi[j]);
                    ((corner >> j) & 1 ? child.lo[j] : child.hi[j]) = m;
                }
                stack.push_back(child);
            }
            continue;
        }
        // Newton's method on the three functions
        VecD<T, 3> y;
        for (std::size_t j = 0; j < 3; ++j)
            y[j] = T(0.5) * (c.lo[j] + c.hi[j]);
        bool converged = false;
        for (int it = 0; it < 12 && !converged; ++it)
        {
            std::array<VecD<T, 3>, 3> m;
            VecD<T, 3> r;
            bool regular = true;
            for (std::size_t i = 0; i < 3; ++i)
            {
                func_gradient<T, 3>(ctx, *functions[i], y, m[i]);
                const T scale = max_abs(std::span<const T>(m[i]));
                regular = regular && scale > T(0) && scale < infinity<T>;
                if (!regular)
                    break;
                r[i] = value<T, 3>(ctx, *functions[i], y) / scale;
                for (T& v : m[i])
                    v /= scale;
            }
            const T det = regular ? determinant<T, 3>(m) : T(0);
            if (!(std::abs(det) > T(0)))
                break;
            T step = T(0);
            VecD<T, 3> next = y;
            for (std::size_t col = 0; col < 3; ++col)
            {
                std::array<VecD<T, 3>, 3> mc = m;
                for (std::size_t i = 0; i < 3; ++i)
                    mc[i][col] = r[i];
                const T d = determinant<T, 3>(mc) / det;
                next[col] -= d;
                step += std::abs(d);
            }
            y = next;
            converged = step <= T(1e-12) * size;
        }
        if (!converged)
            continue;
        const T tol = T(1e-9) * size;
        const VecD<T, 2> base = {y[k == 0 ? 1 : 0], y[k == 2 ? 1 : 2]};
        if (base[0] < lo[0] - tol || base[0] > hi[0] + tol || base[1] < lo[1] - tol || base[1] > hi[1] + tol
            || y[static_cast<std::size_t>(k)] < t_lo - tol || y[static_cast<std::size_t>(k)] > t_hi + tol)
            continue;
        bool seen = false;
        for (const VecD<T, 2>& z : out)
            seen = seen || std::abs(z[0] - base[0]) + std::abs(z[1] - base[1]) <= tol;
        if (!seen)
            out.push_back(base);
    }
}

/// The height direction of a level-2 box below a 3D box: the axis
/// whose lines cross all zero curves best, functions and surface functions, by the
/// smallest over the curves of the mean |d_j f| / |grad f| near them on a grid. Lines
/// along a nearly straight curve would leave nothing to split at, and its fold
/// search nothing to prune. @p k if no curve shows on the grid.
template <std::floating_point T>
int crossing_direction(Context<T>& ctx, const Problem<T, 2>& p, int k)
{
    static constexpr int G = 9;
    VecD<T, 2> score;
    score.fill(infinity<T>);
    bool seen = false;
    const auto accumulate = [&](auto&& value_at, auto&& gradient_at)
    {
        std::array<T, G * G> values;
        T largest = T(0);
        bool pos = false, neg = false;
        const auto point = [&p](int i, int j)
        {
            return VecD<T, 2>{p.lo[0] + (p.hi[0] - p.lo[0]) * (T(i) + T(0.5)) / T(G),
                              p.lo[1] + (p.hi[1] - p.lo[1]) * (T(j) + T(0.5)) / T(G)};
        };
        const auto in_region = [&p](const VecD<T, 2>& y)
        {
            for (const Half<T, 2>& h : p.clips)
                if (h.c[0] * y[0] + h.c[1] * y[1] > h.d)
                    return false;
            return true;
        };
        for (int i = 0; i < G; ++i)
            for (int j = 0; j < G; ++j)
            {
                const VecD<T, 2> y = point(i, j);
                const bool in = in_region(y);
                const T v = in ? value_at(y) : infinity<T>;
                values[static_cast<std::size_t>(i * G + j)] = v;
                if (!in)
                    continue;
                largest = std::max(largest, std::abs(v));
                pos = pos || v > T(0);
                neg = neg || v < T(0);
            }
        if (!(pos && neg))
            return;
        VecD<T, 2> mean{};
        int count = 0;
        for (int i = 0; i < G; ++i)
            for (int j = 0; j < G; ++j)
            {
                if (std::abs(values[static_cast<std::size_t>(i * G + j)]) > T(0.2) * largest)
                    continue;
                VecD<T, 2> g;
                gradient_at(point(i, j), g);
                const T norm = scaled_norm(std::span<const T>(g));
                if (!(norm > T(0)))
                    continue;
                mean[0] += std::abs(g[0]) / norm;
                mean[1] += std::abs(g[1]) / norm;
                ++count;
            }
        if (count == 0)
            return;
        seen = true;
        score[0] = std::min(score[0], mean[0] / T(count));
        score[1] = std::min(score[1], mean[1] / T(count));
    };
    for (const Func<T, 2>& f : p.funcs)
        if (!f.linear)
            accumulate([&](const VecD<T, 2>& y) { return value<T, 2>(ctx, f, y); },
                       [&](const VecD<T, 2>& y, VecD<T, 2>& g) { func_gradient<T, 2>(ctx, f, y, g); });
    for (const Surface<T, 2>& s : p.surfaces)
        accumulate([&](const VecD<T, 2>& y) { return surface_value<T, 2>(ctx, s, y); },
                   [&](const VecD<T, 2>& y, VecD<T, 2>& g) { surface_gradient<T, 2>(ctx, s, y, g); });
    if (!seen)
        return k;
    return score[1] > score[0] ? 1 : 0;
}

/// How near the part a point of a curve must lie to split the base at it, relative to
/// the box's size (may_bound_part).
template <std::floating_point T>
inline constexpr T part_margin = T(0.02);

/// Fold splits (Options::fold_splits) of a level-2 box with height direction
/// k: the points where the zero curve of a curved function or surface function
/// turns along k (a fold, f = d_k f = 0) and where two zero curves cross
/// (f = g = 0). Sub-boxes where bounds show that one of the pair has one sign
/// are dropped; at the depth limit Newton's method runs from the sub-box's
/// centre. On the base segments between them, the roots of the functions on
/// the height lines neither merge nor cross, so that the integrand is smooth.
/// Folds of passive functions (top level of an interface part: only evaluated
/// at the interface's roots) are not needed. @p folds and @p crossings receive
/// the points; false if the search does not settle (its box budget, too many
/// points, or a sub-box at the depth limit where Newton's method does not
/// converge: a singular point, or nearly tangent curves), and the box is
/// bisected instead. Below a 3D box the curves of surface functions are walked on
/// a grid instead of bounded, and only points that may bound the part split it.
template <std::floating_point T>
bool fold_split_points(Context<T>& ctx, const Problem<T, 2>& p, int k, bool top, int surface,
                       std::vector<VecD<T, 2>>& folds, std::vector<VecD<T, 2>>& crossings)
{
    constexpr int max_depth = 9, max_points = 24;
    const bool Top3 = !top;
    // Below a 3D box the curves of surface functions are walked on a grid (see below):
    // the search need not be complete, the accuracy check has what it misses. A curve
    // part has no such check for what it misses (its one curve is the integration
    // domain), and the segments of leaf cells must line up across all curves.
    const bool lenient = Top3 && ctx.partner < 0 && ctx.vis_order == 0;
    constexpr int max_boxes = 400;
    folds.clear();
    crossings.clear();
    struct Member
    {
        int func = -1, surf = -1, ls = 0;
        bool passive = false;
        bool straight = false; ///< a level set affine on the cell: crossings, but no folds
        bool branch = false;   ///< of a surface on a branch (Surface::branch)
    };
    std::vector<Member> members;
    for (int i = 0; i < static_cast<int>(p.funcs.size()); ++i)
    {
        const Func<T, 2>& f = p.funcs[static_cast<std::size_t>(i)];
        if (f.linear)
            continue;
        const bool flat = ctx.flat[static_cast<std::size_t>(f.ls)];
        // below a box without a margin, a straight zero curve's crossings are
        // where a meeting curve leaves through a face
        if (flat && !(Top3 && ctx.top_unmargined))
            continue;
        members.push_back({i, -1, f.ls, top && surface >= 0 && f.ls != surface, flat, false});
    }
    for (int i = 0; i < static_cast<int>(p.surfaces.size()); ++i)
        members.push_back({-1, i, p.surfaces[static_cast<std::size_t>(i)].on.ls, false, false,
                           p.surfaces[static_cast<std::size_t>(i)].branch != 0});

    const auto value_of = [&](const Member& m, const VecD<T, 2>& y)
    {
        return m.func >= 0 ? value<T, 2>(ctx, p.funcs[static_cast<std::size_t>(m.func)], y)
                           : surface_value<T, 2>(ctx, p.surfaces[static_cast<std::size_t>(m.surf)], y);
    };
    const auto gradient_of = [&](const Member& m, const VecD<T, 2>& y, VecD<T, 2>& g)
    {
        if (m.func >= 0)
            func_gradient<T, 2>(ctx, p.funcs[static_cast<std::size_t>(m.func)], y, g);
        else
            surface_gradient<T, 2>(ctx, p.surfaces[static_cast<std::size_t>(m.surf)], y, g);
    };
    // bounds on a sub-box: may it vanish, and the sign of its derivative along k (0: unknown, 2: none along k)
    thread_local BoxBernstein<T> form, deriv;
    thread_local std::vector<T> work;
    // and whether its gradient may vanish there (both derivatives' bounds hold 0)
    const auto bounds_of = [&](const Member& m, const VecD<T, 2>& lo, const VecD<T, 2>& hi, int& dk,
                               bool* singular = nullptr)
    {
        dk = 0;
        if (singular != nullptr)
            *singular = true;
        if (m.func >= 0)
        {
            const Ranges<T, 2> r
                = curved_ranges<T, 2>(ctx, p.funcs[static_cast<std::size_t>(m.func)], lo, hi, form, deriv, work);
            if (!r.valid)
                return true;
            if (r.has_derivatives)
            {
                dk = r.dmin[k] > T(0) ? 1 : (r.dmax[k] < T(0) ? -1 : 0);
                // a straight curve along k: no fold, its point in the base is its trace's
                const int o = 1 - k;
                if (std::max(std::abs(r.dmin[k]), std::abs(r.dmax[k]))
                    <= T(1e-12) * std::max(std::abs(r.dmin[o]), std::abs(r.dmax[o])))
                    dk = 2;
                if (singular != nullptr)
                    *singular = !(r.dmin[0] > T(0) || r.dmax[0] < T(0) || r.dmin[1] > T(0) || r.dmax[1] < T(0));
            }
            return r.sign == 0;
        }
        VecD<T, 2> dmin{}, dmax{};
        T centre = T(0);
        const SurfaceCase c
            = surface_bounds<T, 2>(ctx, p.surfaces[static_cast<std::size_t>(m.surf)], lo, hi, dmin, dmax, centre);
        if (c == SurfaceCase::none)
            return false;
        if (c == SurfaceCase::unbounded)
            return true;
        dk = dmin[k] > T(0) ? 1 : (dmax[k] < T(0) ? -1 : 0);
        if (singular != nullptr)
            *singular = !(dmin[0] > T(0) || dmax[0] < T(0) || dmin[1] > T(0) || dmax[1] < T(0));
        T reach = T(0);
        for (int j = 0; j < 2; ++j)
            reach += std::max(std::abs(dmin[j]), std::abs(dmax[j])) * T(0.5) * (hi[j] - lo[j]);
        return !(std::abs(centre) > reach);
    };
    const T size = std::max(p.hi[0] - p.lo[0], p.hi[1] - p.lo[1]);
    // Newton's method on (F0, F1) from y, F and its Jacobian given by eval
    const auto newton = [&](auto&& eval, VecD<T, 2>& y)
    {
        for (int it = 0; it < 20; ++it)
        {
            VecD<T, 2> f;
            std::array<VecD<T, 2>, 2> j;
            if (!eval(y, f, j))
                return false;
            // each row by its largest entry: level sets scaled by 1e+-200 take the same steps
            for (int r = 0; r < 2; ++r)
            {
                const T largest = std::max(std::abs(j[r][0]), std::abs(j[r][1]));
                if (!(largest > T(0)) || !(largest < infinity<T>))
                    return false;
                f[r] /= largest;
                j[r][0] /= largest;
                j[r][1] /= largest;
            }
            const T det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
            if (!(std::abs(det) > T(0)))
                return false;
            const T d0 = (f[0] * j[1][1] - f[1] * j[0][1]) / det, d1 = (j[0][0] * f[1] - j[1][0] * f[0]) / det;
            y[0] -= d0;
            y[1] -= d1;
            if (!(std::abs(y[0] - p.lo[0]) + std::abs(y[1] - p.lo[1]) < T(4) * size))
                return false;
            if (std::abs(d0) + std::abs(d1) <= T(1e-12) * size)
                return true;
        }
        return false;
    };
    const auto inside = [&](const VecD<T, 2>& y, const VecD<T, 2>& lo, const VecD<T, 2>& hi)
    {
        const T tol = T(1e-9) * size;
        if (lenient)
            for (const Half<T, 2>& h : p.clips)
                if (h.c[0] * y[0] + h.c[1] * y[1] > h.d + tol)
                    return false;
        return y[0] >= lo[0] - tol && y[0] <= hi[0] + tol && y[1] >= lo[1] - tol && y[1] <= hi[1] + tol;
    };
    const auto add = [&](std::vector<VecD<T, 2>>& list, const VecD<T, 2>& y)
    {
        for (const VecD<T, 2>& z : list)
            if (std::abs(z[0] - y[0]) + std::abs(z[1] - y[1]) <= T(1e-9) * size)
                return;
        list.push_back(y);
    };
    // subdivide [lo, hi] while keep(lo, hi) holds; at the depth limit, solve from the centre
    struct Cell
    {
        VecD<T, 2> lo, hi;
        int depth;
    };
    std::vector<Cell> stack;
    int boxes = 0;
    bool unresolved = false;
    const auto search = [&](auto&& keep, auto&& solve)
    {
        stack.assign(1, Cell{p.lo, p.hi, 0});
        boxes = 0;
        while (!stack.empty())
        {
            const Cell c = stack.back();
            stack.pop_back();
            if (++boxes > max_boxes || unresolved)
                return false;
            if (!keep(c.lo, c.hi))
                continue;
            if (c.depth == max_depth)
            {
                solve(c.lo, c.hi);
                continue;
            }
            const T m0 = T(0.5) * (c.lo[0] + c.hi[0]), m1 = T(0.5) * (c.lo[1] + c.hi[1]);
            stack.push_back({{c.lo[0], c.lo[1]}, {m0, m1}, c.depth + 1});
            stack.push_back({{m0, c.lo[1]}, {c.hi[0], m1}, c.depth + 1});
            stack.push_back({{c.lo[0], m1}, {m0, c.hi[1]}, c.depth + 1});
            stack.push_back({{m0, m1}, {c.hi[0], c.hi[1]}, c.depth + 1});
        }
        return true;
    };

    // (f, d_k f) with the Hessian row from differences of the gradient
    const auto fold_system = [&](const Member& m)
    {
        return [&, m](const VecD<T, 2>& y, VecD<T, 2>& f, std::array<VecD<T, 2>, 2>& j)
        {
            const T eps = T(1e-6) * size;
            VecD<T, 2> g;
            gradient_of(m, y, g);
            f = {value_of(m, y), g[k]};
            j[0] = g;
            for (int d = 0; d < 2; ++d)
            {
                VecD<T, 2> yp = y, gp;
                yp[d] += eps;
                gradient_of(m, yp, gp);
                j[1][d] = (gp[k] - g[k]) / eps;
            }
            return std::isfinite(f[0]) && std::isfinite(f[1]);
        };
    };
    const auto on_its_curve = [&](const Member& m, const VecD<T, 2>& y)
    {
        return m.surf < 0 || surface_point_in_box<T, 2>(ctx, p.surfaces[static_cast<std::size_t>(m.surf)], y);
    };
    // the 3D point of a member's curve above y
    const auto point_of = [&](const Member& m, const VecD<T, 2>& y)
    {
        return m.func >= 0 ? to_u<T, 2>(p.funcs[static_cast<std::size_t>(m.func)], y)
                           : surface_point<T>(ctx, p.surfaces[static_cast<std::size_t>(m.surf)], y);
    };
    // Below a 3D box a point splits the base only if it may bound the part, and two
    // curves cross only where they meet in 3D: where they cross at different heights
    // along the lines of the box, seen from above, each adds its own smooth term to the
    // integrand and neither sees the other.
    const auto splits_at = [&](const Member& m, const VecD<T, 2>& y)
    { return !lenient || may_bound_part<T>(ctx, point_of(m, y), part_margin<T> * size); };
    const auto meet_at = [&](const Member& a, const Member& b, const VecD<T, 2>& y)
    {
        if (!lenient)
            return true;
        const T margin = part_margin<T> * size;
        const auto vanishes = [&](int l, const Vec3<T>& u)
        {
            return std::abs(evaluate(ctx.phis[static_cast<std::size_t>(l)], std::span<const T>(u)))
                   <= margin * ctx.scales[static_cast<std::size_t>(l)];
        };
        // a function's zero set may reach the other's point from above its trace
        const Vec3<T> ua = point_of(a, y), ub = point_of(b, y);
        return (vanishes(a.ls, ub) && may_bound_part<T>(ctx, ub, margin))
               || (vanishes(b.ls, ua) && may_bound_part<T>(ctx, ua, margin));
    };

    // Below a 3D box, bounds of a surface function's derivative over a sub-box prove
    // little (its root moves over the sub-box, and its slope is unbounded at a
    // silhouette), and every sub-box along a curve is kept. So walk the curves on a
    // grid of the base instead: a curve passes through the cells whose corners differ in
    // sign, a fold lies where the slope along k differs in sign on the way through a
    // cell, and two curves may cross only in a cell they share. Newton's method from
    // there. What the grid does not resolve is left to the accuracy check.
    constexpr int grid = 12;
    constexpr int n_nodes = (grid + 1) * (grid + 1);
    std::vector<T> samples;
    const auto node = [&p](int i, int j)
    {
        return VecD<T, 2>{p.lo[0] + (p.hi[0] - p.lo[0]) * T(i) / T(grid),
                          p.lo[1] + (p.hi[1] - p.lo[1]) * T(j) / T(grid)};
    };
    const auto sample_members = [&]()
    {
        if (!samples.empty())
            return;
        samples.resize(members.size() * n_nodes);
        for (std::size_t a = 0; a < members.size(); ++a)
            for (int i = 0; i <= grid; ++i)
                for (int j = 0; j <= grid; ++j)
                    samples[a * n_nodes + static_cast<std::size_t>(i * (grid + 1) + j)] = value_of(members[a], node(i, j));
    };
    // where a member's curve leaves a cell, from the values at its corners
    const auto curve_in_cell = [&](std::size_t member, int i, int j, std::array<VecD<T, 2>, 4>& at)
    {
        const std::array<std::array<int, 2>, 4> corners = {{{i, j}, {i + 1, j}, {i + 1, j + 1}, {i, j + 1}}};
        const auto sample = [&](const std::array<int, 2>& c)
        { return samples[member * n_nodes + static_cast<std::size_t>(c[0] * (grid + 1) + c[1])]; };
        int n = 0;
        for (int e = 0; e < 4; ++e)
        {
            const std::array<int, 2>& a = corners[static_cast<std::size_t>(e)];
            const std::array<int, 2>& b = corners[static_cast<std::size_t>((e + 1) % 4)];
            const T va = sample(a), vb = sample(b);
            if (!(va == T(0) || opposite(va, vb)))
                continue;
            const T w = va == T(0) ? T(0) : va / (va - vb);
            const VecD<T, 2> ya = node(a[0], a[1]), yb = node(b[0], b[1]);
            at[static_cast<std::size_t>(n++)] = {ya[0] + w * (yb[0] - ya[0]), ya[1] + w * (yb[1] - ya[1])};
        }
        return n;
    };
    const auto walk_folds = [&](std::size_t a)
    {
        sample_members();
        const Member& m = members[a];
        std::array<VecD<T, 2>, 4> at;
        for (int i = 0; i < grid; ++i)
            for (int j = 0; j < grid; ++j)
            {
                const int n = curve_in_cell(a, i, j, at);
                if (n < 2)
                    continue;
                int up = -1, down = -1;
                for (int e = 0; e < n; ++e)
                {
                    VecD<T, 2> g;
                    gradient_of(m, at[static_cast<std::size_t>(e)], g);
                    const T norm = scaled_norm(std::span<const T>(g));
                    if (g[k] > T(1e-10) * norm)
                        up = e;
                    else if (g[k] < -T(1e-10) * norm)
                        down = e;
                }
                if (up < 0 || down < 0)
                    continue;
                const VecD<T, 2>&yu = at[static_cast<std::size_t>(up)], &yd = at[static_cast<std::size_t>(down)];
                VecD<T, 2> y = {T(0.5) * (yu[0] + yd[0]), T(0.5) * (yu[1] + yd[1])};
                if (newton(fold_system(m), y) && inside(y, p.lo, p.hi) && on_its_curve(m, y) && splits_at(m, y))
                    add(folds, y);
            }
    };
    // A closed curve smaller than a cell of the grid shows in no sign change, and the
    // check cannot see it either: no point of either rule lands on it, say where two
    // balls barely overlap. It shows as a stationary point where the function takes the
    // opposite sign of the samples around it. Its folds along k lie to both sides.
    const auto walk_loops = [&](std::size_t a)
    {
        sample_members();
        const Member& m = members[a];
        const T eps = T(1e-6) * size;
        const auto hessian = [&](const VecD<T, 2>& y, VecD<T, 2>& g, std::array<VecD<T, 2>, 2>& h)
        {
            gradient_of(m, y, g);
            for (int d = 0; d < 2; ++d)
            {
                VecD<T, 2> yp = y, gp;
                yp[d] += eps;
                gradient_of(m, yp, gp);
                h[0][d] = (gp[0] - g[0]) / eps;
                h[1][d] = (gp[1] - g[1]) / eps;
            }
            return std::isfinite(g[0]) && std::isfinite(g[1]);
        };
        const auto at_node = [&](int i, int j)
        { return samples[a * n_nodes + static_cast<std::size_t>(i * (grid + 1) + j)]; };
        const T reach = T(2) * std::max(p.hi[0] - p.lo[0], p.hi[1] - p.lo[1]) / T(grid);
        for (int i = 0; i <= grid; ++i)
            for (int j = 0; j <= grid; ++j)
            {
                const T v = at_node(i, j);
                if (v == T(0))
                    continue;
                bool smallest = true;
                for (int di = -1; di <= 1 && smallest; ++di)
                    for (int dj = -1; dj <= 1 && smallest; ++dj)
                    {
                        const int ni = i + di, nj = j + dj;
                        if ((di == 0 && dj == 0) || ni < 0 || ni > grid || nj < 0 || nj > grid)
                            continue;
                        const T w = at_node(ni, nj);
                        smallest = !opposite(v, w) && std::abs(w) >= std::abs(v);
                    }
                if (!smallest)
                    continue;
                const VecD<T, 2> start = node(i, j);
                VecD<T, 2> y = start;
                if (!newton(hessian, y) || std::abs(y[0] - start[0]) + std::abs(y[1] - start[1]) > reach)
                    continue;
                const T centre = value_of(m, y);
                if (!opposite(v, centre))
                    continue;
                VecD<T, 2> g;
                std::array<VecD<T, 2>, 2> h;
                hessian(y, g, h);
                const int o = 1 - k;
                if (!(std::abs(h[k][k]) > T(0)))
                    continue;
                const T curvature = h[o][o] - h[o][k] * h[k][o] / h[k][k];
                if (!(std::abs(curvature) > T(0)))
                    continue;
                const T radius = std::sqrt(T(2) * std::abs(centre / curvature));
                for (const T side : {T(-1), T(1)})
                {
                    VecD<T, 2> z = y;
                    z[o] += side * radius;
                    z[k] -= side * radius * h[k][o] / h[k][k];
                    if (newton(fold_system(m), z) && inside(z, p.lo, p.hi) && on_its_curve(m, z) && splits_at(m, z))
                        add(folds, z);
                }
            }
    };
    const auto walk_crossings = [&](std::size_t a, std::size_t b, auto&& system)
    {
        sample_members();
        std::array<VecD<T, 2>, 4> at_a, at_b;
        for (int i = 0; i < grid; ++i)
            for (int j = 0; j < grid; ++j)
            {
                const int na = curve_in_cell(a, i, j, at_a);
                if (na == 0)
                    continue;
                const int nb = curve_in_cell(b, i, j, at_b);
                if (nb == 0)
                    continue;
                VecD<T, 2> y{};
                for (int e = 0; e < na; ++e)
                    for (int d = 0; d < 2; ++d)
                        y[d] += at_a[static_cast<std::size_t>(e)][d] / T(2 * na);
                for (int e = 0; e < nb; ++e)
                    for (int d = 0; d < 2; ++d)
                        y[d] += at_b[static_cast<std::size_t>(e)][d] / T(2 * nb);
                if (newton(system, y) && inside(y, p.lo, p.hi) && on_its_curve(members[a], y)
                    && on_its_curve(members[b], y) && meet_at(members[a], members[b], y))
                    add(crossings, y);
            }
    };

    // folds: f = d_k f = 0
    for (std::size_t a = 0; a < members.size(); ++a)
    {
        const Member& m = members[a];
        if (m.passive || m.straight)
            continue;
        if (lenient && m.surf >= 0)
        {
            walk_folds(a);
            walk_loops(a);
            continue;
        }
        const auto keep = [&](const VecD<T, 2>& lo, const VecD<T, 2>& hi)
        {
            int dk = 0;
            return bounds_of(m, lo, hi, dk) && dk == 0;
        };
        const auto solve = [&](const VecD<T, 2>& lo, const VecD<T, 2>& hi)
        {
            // a sub-box at the depth limit where the gradient may vanish may
            // hold a singular point (two branches of the zero set crossing)
            int dk = 0;
            bool singular = true;
            bounds_of(m, lo, hi, dk, &singular);
            VecD<T, 2> y = {T(0.5) * (lo[0] + hi[0]), T(0.5) * (lo[1] + hi[1])};
            if (singular || !newton(fold_system(m), y))
                unresolved = true; // bisect instead
            else if (inside(y, p.lo, p.hi) && splits_at(m, y))
                add(folds, y);
        };
        if (!search(keep, solve))
            return false;
    }
    // crossings: f = g = 0, for two functions of different level sets or a surface function
    for (std::size_t a = 0; a < members.size(); ++a)
        for (std::size_t b = a + 1; b < members.size(); ++b)
        {
            const Member &ma = members[a], &mb = members[b];
            if (ma.func >= 0 && mb.func >= 0 && ma.ls == mb.ls)
                continue;
            if (ma.surf >= 0 && mb.surf >= 0)
            {
                const Surface<T, 2>& sa = p.surfaces[static_cast<std::size_t>(ma.surf)];
                const Surface<T, 2>& sb = p.surfaces[static_cast<std::size_t>(mb.surf)];
                // the two branches of a pair join at the contacts with a silhouette
                if (sa.branch != 0 && sb.branch != 0 && sa.under.ls == sb.under.ls && sa.on.ls == sb.on.ls)
                    continue;
                // a meeting curve and the silhouette of one of its level sets touch
                const Surface<T, 2>* curve = sa.silhouette ? (sb.silhouette ? nullptr : &sb) : (sb.silhouette ? &sa : nullptr);
                const Surface<T, 2>* outline = curve == &sa ? &sb : &sa;
                if (curve != nullptr && (outline->on.ls == curve->under.ls || outline->on.ls == curve->on.ls))
                {
                    const std::size_t found = folds.size();
                    contact_points<T>(ctx, curve->under, curve->on, outline->under, curve->k, curve->t_lo, curve->t_hi,
                                      p.lo, p.hi, folds);
                    folds.erase(std::remove_if(folds.begin() + static_cast<std::ptrdiff_t>(found), folds.end(),
                                               [&](const VecD<T, 2>& y)
                                               { return !splits_at(sa.silhouette ? mb : ma, y); }),
                                folds.end());
                    continue;
                }
            }
            // A silhouette touches the trace of its level set on a bound of the
            // lines where it crosses that bound (where the extreme point of a
            // line reaches its end): tangent curves, on which Newton's method
            // fails. There f and d_t f on the bound cross transversally.
            const Member& mf = ma.func >= 0 ? ma : mb;
            const Member& ms = ma.func >= 0 ? mb : ma;
            if (mf.func >= 0 && ms.surf >= 0 && p.surfaces[static_cast<std::size_t>(ms.surf)].silhouette)
            {
                const Surface<T, 2>& s = p.surfaces[static_cast<std::size_t>(ms.surf)];
                const Func<T, 2>& f = p.funcs[static_cast<std::size_t>(mf.func)];
                if (f.ls == s.on.ls)
                {
                    Func<T, 2> slope = f;
                    slope.ls = s.under.ls;
                    const auto keep = [&](const VecD<T, 2>& lo, const VecD<T, 2>& hi)
                    {
                        const Ranges<T, 2> r0 = curved_ranges<T, 2>(ctx, f, lo, hi, form, deriv, work);
                        const Ranges<T, 2> r1 = curved_ranges<T, 2>(ctx, slope, lo, hi, form, deriv, work);
                        return (!r0.valid || r0.sign == 0) && (!r1.valid || r1.sign == 0);
                    };
                    const auto solve = [&](const VecD<T, 2>& lo, const VecD<T, 2>& hi)
                    {
                        const auto eval = [&](const VecD<T, 2>& y, VecD<T, 2>& v, std::array<VecD<T, 2>, 2>& j)
                        {
                            v = {value<T, 2>(ctx, f, y), value<T, 2>(ctx, slope, y)};
                            func_gradient<T, 2>(ctx, f, y, j[0]);
                            func_gradient<T, 2>(ctx, slope, y, j[1]);
                            return std::isfinite(v[0]) && std::isfinite(v[1]);
                        };
                        VecD<T, 2> y = {T(0.5) * (lo[0] + hi[0]), T(0.5) * (lo[1] + hi[1])};
                        if (!newton(eval, y))
                            unresolved = unresolved || !lenient;
                        else if (inside(y, p.lo, p.hi) && splits_at(mf, y))
                            add(crossings, y);
                    };
                    if (!search(keep, solve) && !lenient)
                        return false;
                    continue;
                }
            }
            const auto crossing_system = [&](const VecD<T, 2>& y, VecD<T, 2>& f, std::array<VecD<T, 2>, 2>& j)
            {
                f = {value_of(ma, y), value_of(mb, y)};
                gradient_of(ma, y, j[0]);
                gradient_of(mb, y, j[1]);
                return std::isfinite(f[0]) && std::isfinite(f[1]);
            };
            if (lenient && (ma.surf >= 0 || mb.surf >= 0))
            {
                walk_crossings(a, b, crossing_system);
                continue;
            }
            const auto keep = [&](const VecD<T, 2>& lo, const VecD<T, 2>& hi)
            {
                int dk = 0;
                return bounds_of(ma, lo, hi, dk) && bounds_of(mb, lo, hi, dk);
            };
            const auto solve = [&](const VecD<T, 2>& lo, const VecD<T, 2>& hi)
            {
                VecD<T, 2> y = {T(0.5) * (lo[0] + hi[0]), T(0.5) * (lo[1] + hi[1])};
                if (!newton(crossing_system, y))
                    unresolved = unresolved || !lenient; // nearly tangent curves, say: bisect instead
                else if (inside(y, p.lo, p.hi) && meet_at(ma, mb, y))
                    add(crossings, y);
            };
            if (!search(keep, solve))
            {
                if (!lenient)
                    return false;
                walk_crossings(a, b, crossing_system);
            }
        }
    return !unresolved && static_cast<int>(folds.size() + crossings.size()) <= max_points;
}

/// Swap the anchored rules with those of twice the points.
template <std::floating_point T>
void swap_more_points(Context<T>& ctx)
{
    std::swap(ctx.fold_x, ctx.more_fold_x);
    std::swap(ctx.fold_w, ctx.more_fold_w);
    std::swap(ctx.fold_check_x, ctx.more_fold_check_x);
    std::swap(ctx.fold_check_w, ctx.more_fold_check_w);
    std::swap(ctx.cubic_x, ctx.more_cubic_x);
    std::swap(ctx.cubic_w, ctx.more_cubic_w);
    std::swap(ctx.cubic_check_x, ctx.more_cubic_check_x);
    std::swap(ctx.cubic_check_w, ctx.more_cubic_check_w);
    ctx.more_points = !ctx.more_points;
}

/// Swap the Gauss-Legendre rule with the one of the accuracy checks.
template <std::floating_point T>
void swap_rules(Context<T>& ctx)
{
    std::swap(ctx.gauss_x, ctx.check_x);
    std::swap(ctx.gauss_w, ctx.check_w);
    std::swap(ctx.fold_x, ctx.fold_check_x);
    std::swap(ctx.fold_w, ctx.fold_check_w);
    std::swap(ctx.cubic_x, ctx.cubic_check_x);
    std::swap(ctx.cubic_w, ctx.cubic_check_w);
    ctx.q = static_cast<int>(ctx.gauss_x.size());
}

/// All the rules of a box, with @p level points more per direction than the run's.
template <std::floating_point T>
void raise_rules(Context<T>& ctx, int level)
{
    const auto rule = [&ctx](int n, std::vector<T>& x, std::vector<T>& w)
    {
        if (static_cast<int>(x.size()) == n)
            return;
        const std::size_t i = static_cast<std::size_t>(n);
        if (ctx.rule_x.size() <= i)
        {
            ctx.rule_x.resize(i + 1);
            ctx.rule_w.resize(i + 1);
        }
        if (ctx.rule_x[i].empty())
            gauss_legendre<T>(n, ctx.rule_x[i], ctx.rule_w[i]);
        x = ctx.rule_x[i];
        w = ctx.rule_w[i];
    };
    const bool more = ctx.more_points;
    if (more)
        swap_more_points(ctx);
    const int q = ctx.cell_q + level;
    const int af = ctx.opt.anchor_folds;
    rule(q, ctx.gauss_x, ctx.gauss_w);
    rule(q + 1, ctx.check_x, ctx.check_w);
    if (af > 0)
    {
        rule(af * q, ctx.fold_x, ctx.fold_w);
        rule(af * (q + 1), ctx.fold_check_x, ctx.fold_check_w);
        rule((3 * af * q + 1) / 2, ctx.cubic_x, ctx.cubic_w);
        rule((3 * af * (q + 1) + 1) / 2, ctx.cubic_check_x, ctx.cubic_check_w);
    }
    if (af == 1)
    {
        rule(2 * q, ctx.more_fold_x, ctx.more_fold_w);
        rule(2 * (q + 1), ctx.more_fold_check_x, ctx.more_fold_check_w);
        rule(3 * q, ctx.more_cubic_x, ctx.more_cubic_w);
        rule(3 * (q + 1), ctx.more_cubic_check_x, ctx.more_cubic_check_w);
    }
    if (more)
        swap_more_points(ctx);
    ctx.q = q;
    ctx.q_raised = level;
}

/// The dimension reduction of a tightened and analysed box: bisect if it is not
/// certified, otherwise integrate along height lines over the base level. With
/// @p check, a certified box with two free coordinates whose piece is curved
/// is integrated with q and q + 1 points, and bisected while they disagree
/// (Options::tolerance).
template <std::floating_point T, int D, bool Rotated, int Top, typename Emit>
void reduce(Context<T>& ctx, Problem<T, D> p, const Analysis<T, D>& an, const Emit& emit, int depth,
            bool interface, bool check = true)
{
    int k = an.k;
    const bool certified = an.certified;
    if constexpr (D == Top)
    {
        // the height functions below are curved where a curved function may vanish
        ctx.curved = !p.surfaces.empty();
        const Func<T, D>* curved = nullptr;
        for (const Func<T, D>& f : p.funcs)
            if (!f.linear && !ctx.flat[static_cast<std::size_t>(f.ls)])
                curved = &f;
        ctx.curved = ctx.curved || curved != nullptr;
        // what the accuracy checks below take from the box, the same for all parts
        ctx.top_extent = p.hi[k] - p.lo[k];
        ctx.top_unmargined = an.unmargined;
        if constexpr (D == 3)
            ctx.top_clips = Rotated ? nullptr : &p.clips;
        ctx.kappa = 0;
        if (curved != nullptr && ctx.out != nullptr && ctx.opt.tolerance == 0)
        {
            VecD<T, D> centre;
            for (int j = 0; j < D; ++j)
                centre[j] = T(0.5) * (p.lo[j] + p.hi[j]);
            ctx.kappa = curvature(ctx, to_u<T, D>(*curved, centre));
        }
    }
    if (!certified && ctx.opt.diagnose)
        diagnose_failure<T, D>(ctx, p, an, depth);
    // whether the accuracy check below runs on this box (a 2D cell's points
    // where two curves cross have weight 1 each)
    const bool checked = D == 2 && check && ctx.curved && ctx.out != nullptr && !(Top == 2 && ctx.partner >= 0);
    // Options::fold_splits: a level-2 box that no axis certifies splits its
    // base at the folds and crossings along the axis with the best sampled
    // margins instead of bisecting
    bool split = false, block = false;
    std::vector<VecD<T, 2>> split_folds, split_crossings;
    if constexpr (D == 2)
    {
        if (!certified && ctx.opt.fold_splits && ctx.opt.anchor_folds > 0 && ctx.vis_order == 0
            && ctx.splits_blocked == 0)
        {
            const int ks = Top == 3 && ctx.partner < 0 && ctx.vis_order == 0 ? crossing_direction<T>(ctx, p, k)
                                                                             : uncertified_direction<T, D>(ctx, p, an, k);
            split = fold_split_points<T>(ctx, p, ks, D == Top, an.surface, split_folds, split_crossings);
            if (split)
            {
                k = ks;
                ctx.stats->fold_splits += !checked;
                trace_fate(ctx, BoxFate::split);
            }
            else
            {
                block = true;
                // a singular point, or nearly tangent curves: in the box above
                if (Top == 3 && ctx.top_unmargined)
                {
                    ctx.unmargined_failed = true;
                    return;
                }
            }
        }
    }
    if (!certified && !split)
    {
        if (depth < ctx.opt.max_depth && ctx.bisections < ctx.opt.max_bisections)
        {
            // where the search does not settle (a singular point, nearly
            // tangent curves), the boxes within are bisected as before: the
            // lines of a split box would pass arbitrarily close to it
            trace_fate(ctx, BoxFate::bisected);
            ctx.splits_blocked += block;
            bisect<T, D, Rotated, Top>(ctx, std::move(p), emit, depth, interface);
            ctx.splits_blocked -= block;
            return;
        }
        ++ctx.stats->uncertified;
        trace_fate(ctx, BoxFate::uncertified);
        // a curve part's top box makes no surface functions uncertified: its
        // piece of the curve, if any, is lost
        if (D == Top && ctx.partner >= 0)
            ++ctx.stats->curve_lost;
        k = uncertified_direction<T, D>(ctx, p, an, k);
    }
    if constexpr (D == 2)
    {
        if (checked && (certified || split))
        {
            CellPoints<T>& out = *ctx.out;
            const auto measure = [&out](int from, int to)
            {
                T sum = 0;
                for (int i = from; i < to; ++i)
                    sum += out.weights[static_cast<std::size_t>(i)];
                return sum;
            };
            const auto truncate = [&out](int n)
            {
                out.points.resize(3 * static_cast<std::size_t>(n));
                out.weights.resize(static_cast<std::size_t>(n));
            };
            const int start = out.n_points();
            ctx.anchored_any = false;
            // the second run takes on the folds the first finds
            const bool anchored = ctx.opt.anchor_folds > 0;
            if (anchored)
            {
                ctx.fold_pass = 1;
                ctx.box_folds.clear();
                ctx.column_y.clear();
                ctx.column_start.assign(1, 0);
                ctx.column_fold.clear();
            }
            reduce<T, D, Rotated, Top>(ctx, p, an, emit, depth, interface, false);
            const int mid = out.n_points();
            if (anchored)
            {
                ctx.column_order.resize(ctx.column_y.size());
                for (std::size_t i = 0; i < ctx.column_order.size(); ++i)
                    ctx.column_order[i] = static_cast<int>(i);
                std::sort(ctx.column_order.begin(), ctx.column_order.end(),
                          [&](int a, int b)
                          { return ctx.column_y[static_cast<std::size_t>(a)] < ctx.column_y[static_cast<std::size_t>(b)]; });
                ctx.fold_pass = 2;
            }
            // the second run counts nothing
            const Stats stats = *ctx.stats;
            const int bisections = ctx.bisections;
            swap_rules(ctx);
            const bool moments = ctx.opt.moment_tolerance > 0 && ctx.partner < 0;
            ctx.measure_only = !interface && !ctx.any_analytic && !moments;
            reduce<T, D, Rotated, Top>(ctx, p, an, emit, depth, interface, false);
            ctx.measure_only = false;
            swap_rules(ctx);
            ctx.fold_pass = 0;
            *ctx.stats = stats;
            ctx.bisections = bisections;
            const T a = measure(start, mid), b = measure(mid, out.n_points());
            const auto moment = [&out, &ctx](int from, int to)
            {
                T sum = 0;
                for (int i = from; i < to; ++i)
                {
                    T shape = 0;
                    for (int j = 0; j < ctx.tdim; ++j)
                    {
                        T xi = ctx.ref_origin[static_cast<std::size_t>(j)];
                        for (std::size_t c = 0; c < 3; ++c)
                            xi += ctx.ref_jacobian[static_cast<std::size_t>(j)][c] * out.points[3 * static_cast<std::size_t>(i) + c];
                        const T t = 2 * xi - 1;
                        T low = T(1), high = t;
                        for (int m = 1; m < ctx.moment_degree; ++m)
                        {
                            const T next = (T(2 * m + 1) * t * high - T(m) * low) / T(m + 1);
                            low = high;
                            high = next;
                        }
                        shape += high * high;
                    }
                    sum += out.weights[static_cast<std::size_t>(i)] * shape;
                }
                return sum;
            };
            const T moment_a = moments ? moment(start, mid) : T(0), moment_b = moments ? moment(mid, out.n_points()) : T(0);
            truncate(mid);
            if (mid == start)
                return;
            T rel = static_cast<T>(ctx.opt.tolerance);
            if (rel == T(0))
            {
                // the rule's own rate: q Gauss-Legendre points on a segment of
                // length h err by about (h / (4 r))^(2q) where the integrand is
                // analytic within r of it; r = 1 / (2 kappa) leaves room for the
                // branch points of the height functions, which come closer
                // than the radius of curvature
                rel = std::pow(ctx.cell_size * ctx.kappa / T(2), T(2 * ctx.cell_q));
            }
            // A box split at its folds and crossings is not bisected for
            // certification, so it is larger than a certified one and its q
            // points need not reach that accuracy: it must agree with q + 1 to
            // split_tolerance, or be bisected. Likewise below a 3D box without a
            // margin along its height direction.
            if ((split || (Top == 3 && ctx.top_unmargined)) && ctx.opt.tolerance == 0)
                rel = std::min(rel, ctx.split_relative);
            // Relative to the piece's measure on interfaces and curves; a volume
            // part's to the whole box's, as the difference: both sides of a
            // level set then bisect alike and still add up to the box exactly.
            T scale = std::abs(b);
            if (!interface)
            {
                scale = ctx.measure_scale * (Top == 3 ? ctx.top_extent : T(1));
                for (int j = 0; j < D; ++j)
                    scale *= p.hi[j] - p.lo[j];
            }
            // Below a bisection for the moment, a box is held to the tolerances of the box
            // that failed: not to shares by its size, which the part with a singular point
            // would never pass, nor to a certified box's rate, which would undo the measure.
            const T tol = ctx.measure_floor > T(0) ? ctx.measure_floor : rel * scale + accuracy_floor<T> * ctx.measure_scale;
            const T moment_tol = ctx.moment_floor > T(0) ? ctx.moment_floor : ctx.moment_relative * scale + accuracy_floor<T> * ctx.measure_scale;
            const bool measure_holds = !(std::abs(a - b) > tol);
            const bool moment_holds = !(std::abs(moment_a - moment_b) > moment_tol);
            if (measure_holds && moment_holds)
                return;
            // Anchored segments with one factor are short of points before the box is too
            // large: twice as many in that direction first, which is often all it takes.
            if (ctx.opt.anchor_folds == 1 && !ctx.more_points && ctx.anchored_any)
            {
                truncate(start);
                swap_more_points(ctx);
                reduce<T, D, Rotated, Top>(ctx, p, an, emit, depth, interface);
                swap_more_points(ctx);
                return;
            }
            // q points are exact to degree 2 q - 1 on a flat piece, and curvature takes the
            // moment past that: a point more per direction gains digits where a bisection
            // would gain little for four times the points. The measure failing is something
            // else, which only bisection resolves.
            const int level = ctx.q_raised;
            if (measure_holds && level < moment_raises)
            {
                truncate(start);
                raise_rules(ctx, level + 1);
                reduce<T, D, Rotated, Top>(ctx, p, an, emit, depth, interface);
                raise_rules(ctx, level);
                return;
            }
            if (depth >= ctx.opt.max_depth || ctx.bisections >= ctx.opt.max_bisections)
                return;
            if (measure_holds)
            {
                // The moment alone bisects a few times, as far as bisection helps a rule
                // that is right about the measure: not without end at a singular point (a
                // cone's apex), and it uses up no attempt below a box without a margin.
                if (depth >= moment_depth || ctx.moment_splits >= moment_bisections)
                    return;
                ++ctx.moment_splits;
                ctx.moment_only = true;
                const T saved_moment = ctx.moment_floor, saved_measure = ctx.measure_floor;
                ctx.moment_floor = moment_tol;
                ctx.measure_floor = std::min(tol, ctx.split_relative * scale + accuracy_floor<T> * ctx.measure_scale);
                truncate(start);
                trace_fate(ctx, BoxFate::refined);
                if (level > 0)
                    raise_rules(ctx, 0);
                bisect<T, D, Rotated, Top>(ctx, std::move(p), emit, depth, interface);
                if (level > 0)
                    raise_rules(ctx, level);
                ctx.moment_floor = saved_moment;
                ctx.measure_floor = saved_measure;
                return;
            }
            // below a bisection for the moment, the same few bisections answer for the measure
            if (ctx.moment_floor > T(0))
            {
                if (ctx.moment_splits >= moment_bisections)
                    return;
                ++ctx.moment_splits;
            }
            truncate(start);
            trace_fate(ctx, BoxFate::refined);
            if (level > 0)
                raise_rules(ctx, 0);
            bisect<T, D, Rotated, Top>(ctx, std::move(p), emit, depth, interface);
            if (level > 0)
                raise_rules(ctx, level);
            return;
        }
    }

    if (ctx.trace != nullptr)
        ctx.trace->height[static_cast<std::size_t>(ctx.trace_box)] = k;

    // bounds of the height lines: box faces and clip planes
    std::vector<Bound<T, D>> lowers, uppers;
    Problem<T, D - 1> base_box;
    for (int jb = 0; jb < D - 1; ++jb)
    {
        const int j = jb < k ? jb : jb + 1;
        base_box.lo[jb] = p.lo[j];
        base_box.hi[jb] = p.hi[j];
    }
    {
        Bound<T, D> l, u;
        l.alpha = p.lo[k];
        u.alpha = p.hi[k];
        l.kind = 2;
        u.kind = 3;
        lowers.push_back(l);
        uppers.push_back(u);
    }
    for (const Half<T, D>& h : p.clips)
    {
        const T ck = h.c[k];
        if (std::abs(ck) <= tiny<T>)
        {
            Half<T, D - 1> hb;
            for (int jb = 0; jb < D - 1; ++jb)
                hb.c[jb] = h.c[jb < k ? jb : jb + 1];
            hb.d = h.d;
            base_box.clips.push_back(hb);
            continue;
        }
        Bound<T, D> b;
        b.alpha = h.d / ck;
        b.kind = ck > T(0) ? 5 : 4;
        for (int jb = 0; jb < D - 1; ++jb)
            b.beta[jb] = -h.c[jb < k ? jb : jb + 1] / ck;
        (ck > T(0) ? uppers : lowers).push_back(b);
    }
    // Fourier-Motzkin: the height line is non-empty where every lower <= every upper
    for (std::size_t i = 0; i < lowers.size(); ++i)
        for (std::size_t j = 0; j < uppers.size(); ++j)
        {
            if (i == 0 && j == 0)
                continue;
            Half<T, D - 1> hb;
            for (int jb = 0; jb < D - 1; ++jb)
                hb.c[jb] = lowers[i].beta[jb] - uppers[j].beta[jb];
            hb.d = uppers[j].alpha - lowers[i].alpha;
            base_box.clips.push_back(hb);
        }

    // Drop bounds that are never active on the base region: a lower bound below
    // another lower bound at every vertex of the base polytope (bounds are affine),
    // or an upper bound above another upper bound. Their restrictions and switch
    // functions would only cause needless bisections.
    if (ctx.opt.prune_bounds && lowers.size() + uppers.size() > 2)
    {
        std::vector<VecD<T, D - 1>> corners;
        polytope_vertices<T, D - 1>(base_box.lo, base_box.hi, base_box.clips, corners);
        const T tol = compare_tol<T>;
        auto prune = [&](std::vector<Bound<T, D>>& list, bool lower)
        {
            std::vector<Bound<T, D>> kept;
            for (std::size_t i = 0; i < list.size(); ++i)
            {
                bool dominated = false;
                for (std::size_t j = 0; j < list.size() && !dominated; ++j)
                {
                    if (j == i)
                        continue;
                    bool all = !corners.empty();
                    for (const VecD<T, D - 1>& y : corners)
                    {
                        const T di = bound_value<T, D>(list[i], y), dj = bound_value<T, D>(list[j], y);
                        // j at least as tight everywhere; ties keep the earlier bound
                        const bool tighter = lower ? (dj > di + tol || (std::abs(dj - di) <= tol && j < i))
                                                   : (dj < di - tol || (std::abs(dj - di) <= tol && j < i));
                        all &= tighter || std::abs(dj - di) <= tol;
                        if (!all)
                            break;
                    }
                    // all corners tied: keep the earlier one
                    if (all)
                    {
                        bool strictly = false;
                        for (const VecD<T, D - 1>& y : corners)
                        {
                            const T di = bound_value<T, D>(list[i], y), dj = bound_value<T, D>(list[j], y);
                            strictly |= lower ? dj > di + tol : dj < di - tol;
                        }
                        dominated = strictly || j < i;
                    }
                }
                if (!dominated)
                    kept.push_back(list[i]);
            }
            list = std::move(kept);
        };
        prune(lowers, true);
        prune(uppers, false);
    }

    // Duplicates (a clip in a box face, say) would be counted twice below.
    auto dedupe = [](std::vector<Bound<T, D>>& list)
    {
        std::vector<Bound<T, D>> kept;
        for (const Bound<T, D>& b : list)
        {
            bool same = false;
            for (const Bound<T, D>& c : kept)
            {
                bool equal = std::abs(b.alpha - c.alpha) <= compare_tol<T>;
                for (int jb = 0; jb < D - 1; ++jb)
                    equal &= std::abs(b.beta[jb] - c.beta[jb]) <= compare_tol<T>;
                same |= equal;
            }
            if (!same)
                kept.push_back(b);
        }
        list = std::move(kept);
    };
    dedupe(lowers);
    dedupe(uppers);

    std::vector<Bound<T, D>> linear_roots;
    for (const Func<T, D>& f : p.funcs)
        if (f.linear && std::abs(f.a[k]) > tiny<T>)
        {
            Bound<T, D> r;
            r.kind = 6;
            r.alpha = -f.c / f.a[k];
            for (int jb = 0; jb < D - 1; ++jb)
                r.beta[jb] = -f.a[jb < k ? jb : jb + 1] / f.a[k];
            linear_roots.push_back(r);
        }

    // corners: the roots of two surface functions crossing (three level sets)
    std::vector<T> corners;
    if constexpr (D == 2)
    {
        for (const std::vector<VecD<T, 2>>* points : {&split_folds, &split_crossings})
            for (const VecD<T, 2>& y : *points)
                corners.push_back(y[1 - k]);
        if (certified)
            for (std::size_t i = 0; i < p.surfaces.size(); ++i)
                for (std::size_t j = i + 1; j < p.surfaces.size(); ++j)
                    surface_crossings<T>(ctx, p.surfaces[i], p.surfaces[j], k, p.lo[k], p.hi[k], base_box.lo[0],
                                         base_box.hi[0], corners);
    }

    // the folds of the curved functions near a certified level-2 box, and of
    // the level above where they cross its bounds, at which the segments of its
    // base line are anchored (Options::anchor_folds)
    std::vector<T> folds;
    if constexpr (D == 2)
        if (ctx.opt.anchor_folds > 0 && (certified || split) && ctx.vis_order == 0)
        {
            if (ctx.fold_pass == 2)
                folds = ctx.box_folds;
            else
            {
                for (const VecD<T, 2>& y : split_folds)
                    folds.push_back(y[1 - k]);
                fold_anchors<T>(ctx, p, k, folds);
                if (Top == 3 && ctx.top_unmargined)
                    for (const Surface<T, 2>& s : p.surfaces)
                        silhouette_anchors<T>(ctx, s, p, k, folds);
                bound_folds<T>(ctx, p.heights, k, lowers, base_box.lo[0], base_box.hi[0], folds);
                bound_folds<T>(ctx, p.heights, k, uppers, base_box.lo[0], base_box.hi[0], folds);
                if (ctx.fold_pass == 1)
                    ctx.box_folds = folds;
            }
        }

    // The base of the height lines between the given lower and upper bounds, on
    // the base polytope cut by the given clips: its functions are where the
    // integrand along the lines changes form.
    const auto integrate_base = [&](const std::vector<Bound<T, D>>& lows, const std::vector<Bound<T, D>>& ups,
                                    const std::vector<Half<T, D - 1>>& clips)
    {
        Problem<T, D - 1> base;
        base.lo = base_box.lo;
        base.hi = base_box.hi;
        base.clips = base_box.clips;
        base.clips.insert(base.clips.end(), clips.begin(), clips.end());
        std::vector<Bound<T, D>> bounds = lows;
        bounds.insert(bounds.end(), ups.begin(), ups.end());
        // A passive function has no roots that matter here, but its restrictions
        // to the bounds pair with the interface's: where the interface's zero set
        // leaves through a face and crosses the other zero set there, the
        // integrand below has a kink, which only their surface function marks.
        for (std::size_t i = 0; i < p.funcs.size(); ++i)
        {
            const Func<T, D>& f = p.funcs[i];
            if (f.linear && std::abs(f.a[k]) <= tiny<T>)
            {
                base.funcs.push_back(restrict_to<T, D>(f, k, Bound<T, D>{})); // independent of y_k
                continue;
            }
            for (const Bound<T, D>& b : bounds) // a root reaches a bound
                base.funcs.push_back(restrict_to<T, D>(f, k, b));
            if (!f.linear)
                for (const Bound<T, D>& r : linear_roots) // a curved root meets a linear one
                    base.funcs.push_back(restrict_to<T, D>(f, k, r));
        }
        for (std::size_t i = 0; i < linear_roots.size(); ++i)
            for (std::size_t j = i + 1; j < linear_roots.size(); ++j)
                base.funcs.push_back(difference<T, D>(linear_roots[i], linear_roots[j]));
        for (std::size_t i = 0; i < lows.size(); ++i) // the active lower bound changes
            for (std::size_t j = i + 1; j < lows.size(); ++j)
                base.funcs.push_back(difference<T, D>(lows[i], lows[j]));
        for (std::size_t i = 0; i < ups.size(); ++i)
            for (std::size_t j = i + 1; j < ups.size(); ++j)
                base.funcs.push_back(difference<T, D>(ups[i], ups[j]));
        if constexpr (D <= 2)
        {
            for (const Surface<T, D>& s : p.surfaces)
            {
                for (const Bound<T, D>& b : bounds)
                    base.surfaces.push_back(restrict_surface<T, D>(s, k, b));
                for (const Bound<T, D>& r : linear_roots)
                    base.surfaces.push_back(restrict_surface<T, D>(s, k, r));
            }
        }
        if constexpr (D == 2)
        {
            for (const T y : corners)
            {
                Func<T, 1> f;
                f.origin = 12;
                f.linear = true;
                f.a[0] = T(1);
                f.c = -y;
                base.funcs.push_back(f);
            }
        }
        // the height functions of the integrand along the lines, at whose folds
        // the height lines of the base are anchored (Options::anchor_folds): the
        // roots of the curved functions the lines take, not of passive ones
        if constexpr (D == 3)
        {
            if (ctx.opt.anchor_folds > 0 && certified && ctx.vis_order == 0 && ctx.partner < 0)
                for (std::size_t i = 0; i < p.funcs.size(); ++i)
                {
                    const Func<T, D>& f = p.funcs[i];
                    if (f.linear || ctx.flat[static_cast<std::size_t>(f.ls)] || an.independent[i] || an.passive[i]
                        || (an.surface >= 0 && f.ls != an.surface))
                        continue;
                    base.heights.push_back({f, k, p.lo[k], p.hi[k], an.two_roots[i] != 0});
                }
            // the silhouette: where the two roots of a line merge, unless d_k f
            // keeps its sign on the box, up to 0 on a face (an equator in a face)
            if (an.unmargined)
                for (std::size_t i = 0; i < p.funcs.size(); ++i)
                {
                    if (!an.two_roots[i])
                        continue;
                    Surface<T, 2> s;
                    s.silhouette = true;
                    s.on = p.funcs[i];
                    s.under = s.on;
                    s.under.ls = ctx.derived[static_cast<std::size_t>(3 * s.on.ls + k)];
                    s.k = k;
                    s.t_lo = p.lo[k];
                    s.t_hi = p.hi[k];
                    const Ranges<T, 3> r = curved_ranges<T, 3>(ctx, s.under, p.lo, p.hi, ctx.inner_form[3],
                                                               ctx.inner_deriv[3], ctx.inner_work[3]);
                    if (!r.valid || r.sign == 0)
                        base.surfaces.push_back(s);
                }
        }
        // Below a box without a margin, a pair's surface runs along one of its
        // functions whose roots the lines take, with one root per line if it can;
        // with two, one surface per branch: every crossing lies on a branch of one.
        if constexpr (D == 3)
        {
            if (an.unmargined)
                for (std::size_t i = 0; i < p.funcs.size(); ++i)
                    for (std::size_t j = i + 1; j < p.funcs.size(); ++j)
                    {
                        std::size_t u = i, v = j;
                        if (p.funcs[u].linear || p.funcs[v].linear || p.funcs[u].ls == p.funcs[v].ls
                            || an.independent[u] || an.independent[v] || (an.passive[u] && an.passive[v]))
                            continue;
                        if (an.passive[u] || (an.two_roots[u] && !an.two_roots[v] && !an.passive[v]))
                            std::swap(u, v);
                        Surface<T, 2> s;
                        s.under = p.funcs[u];
                        s.on = p.funcs[v];
                        s.k = k;
                        s.t_lo = p.lo[k];
                        s.t_hi = p.hi[k];
                        ++ctx.stats->surfaces;
                        if (!an.two_roots[u])
                        {
                            base.surfaces.push_back(s);
                            continue;
                        }
                        s.slope_ls = ctx.derived[static_cast<std::size_t>(3 * s.under.ls + k)];
                        second_margin<T, 3>(ctx, s.under, p.lo, p.hi, k, s.curvature);
                        for (const int branch : {-1, 1})
                        {
                            s.branch = branch;
                            base.surfaces.push_back(s);
                        }
                    }
        }
        // Where two level sets cut a certified box, their roots on the height lines
        // cross on a ridge: the base splits there, where one level set vanishes on
        // the other's zero set. The surface runs along one with a single root per line.
        if (certified && !an.unmargined)
        {
            for (std::size_t i = 0; i < p.funcs.size(); ++i)
                for (std::size_t j = 0; j < p.funcs.size(); ++j)
                {
                    const Func<T, D>& f = p.funcs[i];
                    const Func<T, D>& g = p.funcs[j];
                    // a function constant along y_k: no roots there to cross; a
                    // passive one is only crossed
                    if (f.linear || g.linear || f.ls == g.ls || an.two_roots[i] || an.independent[i]
                        || an.independent[j] || an.passive[i])
                        continue;
                    // one surface per pair: along the first with one root per line
                    if (j < i && !an.two_roots[j] && !an.passive[j])
                        continue;
                    Surface<T, D - 1> s;
                    s.under = f;
                    s.on = g;
                    s.k = k;
                    s.t_lo = p.lo[k];
                    s.t_hi = p.hi[k];
                    base.surfaces.push_back(s);
                    ++ctx.stats->surfaces;
                }
        }

        // the rule along the height lines: the integrand of the base level; dyb
        // is the derivative of yb along a curve part's parameter
        const int box_id = ctx.box_ids[D - 1]++;
        const auto line = [&, k, certified, box_id](const VecD<T, D - 1>& yb, T w, const NodeTag& below,
                                                    const VecD<T, D - 1>& dyb)
        {
            T L = -infinity<T>, U = infinity<T>;
            for (const Bound<T, D>& b : lows)
                L = std::max(L, bound_value<T, D>(b, yb));
            for (const Bound<T, D>& b : ups)
                U = std::min(U, bound_value<T, D>(b, yb));
            if (!(U > L))
            {
                // a crossing point at a corner of the region, where the line has
                // shrunk to it (a triangle's acute corner), still counts
                bool corner = false;
                if constexpr (D == Top && Top == 2)
                    corner = ctx.partner >= 0 && U >= L - point_tol<T>;
                if (!corner)
                    return;
                U = L;
            }
            // the folds of the level above on the line, at which its segments
            // are anchored (height_folds, which takes nodes as scratch)
            std::vector<T>& folds_above = ctx.line_folds[D];
            folds_above.clear();
            if constexpr (D == 2)
                if ((certified || split) && !p.heights.empty() && U > L)
                {
                    if (ctx.fold_pass == 2)
                        continued_height_folds<T>(ctx, p.heights, yb, k, L, U, folds_above);
                    else
                        height_folds<T>(ctx, p.heights, yb, k, L, U, folds_above);
                }
            std::vector<T>& nodes = ctx.nodes[D];
            nodes.clear();
            if constexpr (D == 2 && Top == 3)
            {
                // A curve part, on the base of the top box: the points where the
                // curve crosses the line, the roots of its surface function,
                // instead of segments. Along the curve the surface function stays
                // 0, so its root moves by dt/ds = -(grad_b s . dyb) / d_t s.
                if (ctx.partner >= 0)
                {
                    NodeTag tag = below;
                    tag.box[D - 1] = box_id;
                    tag.node[D - 1] = 0;
                    int index = 0;
                    for (const Surface<T, D>& s : p.surfaces)
                    {
                        if (!partner_surface<T, D>(ctx, s))
                            continue;
                        nodes.clear();
                        if (certified)
                        {
                            const T sl = surface_value<T, D>(ctx, s, insert<T, D>(yb, k, L)),
                                    su = surface_value<T, D>(ctx, s, insert<T, D>(yb, k, U));
                            if (opposite(sl, su))
                                nodes.push_back(illinois_root(
                                    [&](T t) { return surface_value<T, D>(ctx, s, insert<T, D>(yb, k, t)); }, L, U,
                                    sl, su));
                        }
                        else
                            surface_line_roots<T>(ctx, surface_on_line<T, D>(s, yb, k), L, U, nodes);
                        std::sort(nodes.begin(), nodes.end());
                        for (const T t : nodes)
                        {
                            const VecD<T, D> y = insert<T, D>(yb, k, t);
                            // leaf cells: only crossings in the box above, as below
                            if (ctx.vis_order > 0)
                            {
                                const T r = surface_root<T, D>(ctx, s, y);
                                if (r < s.t_lo || r > s.t_hi)
                                    continue;
                            }
                            VecD<T, D> g;
                            surface_gradient<T, D>(ctx, s, y, g);
                            if (!(std::abs(g[k]) > T(0)))
                                continue;
                            T slope = T(0);
                            for (int jb = 0; jb < D - 1; ++jb)
                                slope += g[jb < k ? jb : jb + 1] * dyb[jb];
                            tag.segment[D - 1] = index++;
                            emit(y, w, tag, insert<T, D>(dyb, k, -slope / g[k]));
                        }
                    }
                    return;
                }
            }
            if constexpr (D == Top)
                std::fill(ctx.top_line_set.begin(), ctx.top_line_set.end(), char(0));
            for (std::size_t i = 0; i < p.funcs.size(); ++i)
            {
                const Func<T, D>& f = p.funcs[i];
                if (f.linear)
                {
                    if (std::abs(f.a[k]) > tiny<T>)
                    {
                        T t = f.c;
                        for (int jb = 0; jb < D - 1; ++jb)
                            t += f.a[jb < k ? jb : jb + 1] * yb[jb];
                        t = -t / f.a[k];
                        if (t > L && t < U)
                            nodes.push_back(t);
                    }
                    continue;
                }
                // An interface part: the roots of the interface's level set (in a
                // curve part, of the box's); the others are evaluated at them.
                // Leaf cells do not split at the others' restrictions below
                // either, which only pair with the interface's in surface
                // functions: their roots may cross a surface function's inside a
                // base segment, and the segments of a leaf would no longer line up.
                if (an.surface >= 0 && f.ls != an.surface && (D == Top || ctx.vis_order > 0))
                    continue;
                BoxBernstein<T>* form = &ctx.line[D];
                if (U > L)
                {
                    if constexpr (D == Top)
                    {
                        if (!is_analytic(ctx, f.ls))
                        {
                            form = &ctx.top_lines[static_cast<std::size_t>(f.ls)];
                            ctx.top_line_set[static_cast<std::size_t>(f.ls)] = 1;
                        }
                    }
                    const std::size_t found = nodes.size();
                    curved_roots<T, D>(ctx, f, yb, k, L, U, certified, an.two_roots[i] != 0, *form, nodes);
                    // the roots of a curve outside the part do not split the line
                    if constexpr (D == 2 && Top == 3)
                        if (ctx.partner < 0 && ctx.vis_order == 0)
                            nodes.erase(std::remove_if(nodes.begin() + static_cast<std::ptrdiff_t>(found), nodes.end(),
                                                       [&](T t)
                                                       {
                                                           return !may_bound_part<T>(ctx, to_u<T, D>(f, insert<T, D>(yb, k, t)),
                                                                                     part_margin<T> * (U - L));
                                                       }),
                                        nodes.end());
                }
                if constexpr (D == Top && Top == 2)
                {
                    // a crossing point at an end of the line counts there too
                    if (ctx.partner >= 0)
                        for (const T end : {L, U})
                        {
                            const VecD<T, D> y = insert<T, D>(yb, k, end);
                            if (std::abs(value<T, D>(ctx, f, y))
                                <= point_tol<T> * std::abs(derivative_along<T, D>(ctx, f, y, k)))
                                nodes.push_back(end);
                            if (!(U > L))
                                break;
                        }
                }
            }
            if constexpr (D <= 2)
            {
                for (const Surface<T, D>& s : p.surfaces)
                {
                    const std::size_t found = nodes.size();
                    if (certified)
                    {
                        const T sl = surface_value<T, D>(ctx, s, insert<T, D>(yb, k, L)),
                                su = surface_value<T, D>(ctx, s, insert<T, D>(yb, k, U));
                        if (opposite(sl, su))
                            nodes.push_back(illinois_root(
                                [&](T t) { return surface_value<T, D>(ctx, s, insert<T, D>(yb, k, t)); }, L, U, sl, su));
                    }
                    else
                        surface_line_roots<T>(ctx, surface_on_line<T, D>(s, yb, k), L, U, nodes);
                    if constexpr (D == 2 && Top == 3)
                        if (ctx.partner < 0 && ctx.vis_order == 0)
                            nodes.erase(std::remove_if(nodes.begin() + static_cast<std::ptrdiff_t>(found), nodes.end(),
                                                       [&](T t)
                                                       {
                                                           return !may_bound_part<T>(ctx, surface_point<T>(ctx, s, insert<T, D>(yb, k, t)),
                                                                                     part_margin<T> * (U - L));
                                                       }),
                                        nodes.end());
                    // A silhouette's roots are the folds of the level above on this
                    // line, also where its zero set in the plane through the line
                    // never leaves the box, which height_folds would not see.
                    if (s.silhouette)
                        folds_above.insert(folds_above.end(), nodes.begin() + static_cast<std::ptrdiff_t>(found),
                                           nodes.end());
                    // Leaf cells: only crossings in the box above. Beyond it the
                    // surface function follows the root of under, then clamps it;
                    // such roots may vanish inside a base segment, and the
                    // segments of a leaf would no longer line up.
                    if (ctx.vis_order > 0)
                        nodes.erase(std::remove_if(nodes.begin() + static_cast<std::ptrdiff_t>(found), nodes.end(),
                                                   [&](T t)
                                                   {
                                                       const T r = surface_root<T, D>(ctx, s, insert<T, D>(yb, k, t));
                                                       return r < s.t_lo || r > s.t_hi;
                                                   }),
                                    nodes.end());
                }
            }

            NodeTag tag = below;
            tag.box[D - 1] = box_id;
            if constexpr (D == Top)
            {
                if (interface)
                {
                    // the roots of the level set on each line, weighted by the
                    // physical surface measure |det J| |J^-T grad phi| / |d_k phi|;
                    // on a curve by the physical speed |J dU/ds|, the root moving
                    // by dr/ds = -(grad_b phi . dyb) / d_k phi along it; a 2D curve
                    // part's crossing points weigh 1
                    std::sort(nodes.begin(), nodes.end());
                    const Source<T>& phi = ctx.phis[static_cast<std::size_t>(an.surface)];
                    for (std::size_t r = 0; r < nodes.size(); ++r)
                    {
                        tag.segment[D - 1] = static_cast<int>(r);
                        const VecD<T, D> y = insert<T, D>(yb, k, nodes[r]);
                        const Vec3<T> u = to_u<T, D>(p.frame, y);
                        if (ctx.vis_order == 0 && an.decided == 0 && !holds<T>(ctx, an.terms, u, T(0)))
                            continue;
                        Vec3<T> g = {0, 0, 0};
                        gradient(phi, std::span<const T>(u), std::span<T>(g));
                        // the cell's own directions only: on a 2D cell, the
                        // derivative of an analytic level set normal to it (off
                        // the plane z = 0, or off a facet's plane) is no part of
                        // the curve's gradient
                        for (int i = Top; i < 3; ++i)
                            g[i] = T(0);
                        // grad phi in the level's (possibly rotated) coordinates
                        VecD<T, D> dphi{};
                        for (int j = 0; j < D; ++j)
                            for (int i = 0; i < 3; ++i)
                                dphi[j] += p.frame.A[i][j] * g[i];
                        const T dk = std::abs(dphi[k]);
                        if (!(dk > T(0)))
                            continue;
                        if (ctx.partner < 0)
                        {
                            const T gnorm = scaled_norm(std::span<const T>(g));
                            Vec3<T> gx = {0, 0, 0};
                            for (int i = 0; i < 3; ++i)
                                for (int m = 0; m < 3; ++m)
                                    gx[i] += ctx.inv[m][i] * g[m] / gnorm;
                            emit(y, w * gnorm / dk * ctx.detj * scaled_norm(std::span<const T>(gx)), tag,
                                 VecD<T, D>{});
                        }
                        else if constexpr (Top == 2)
                            emit(y, w, tag, VecD<T, D>{});
                        else
                        {
                            T slope = T(0);
                            for (int jb = 0; jb < D - 1; ++jb)
                                slope += dphi[jb < k ? jb : jb + 1] * dyb[jb];
                            const VecD<T, D> dy = insert<T, D>(dyb, k, -slope / dphi[k]);
                            Vec3<T> dx = {0, 0, 0};
                            for (int i = 0; i < 3; ++i)
                            {
                                T du = T(0);
                                for (int j = 0; j < D; ++j)
                                    du += p.frame.A[i][j] * dy[j];
                                for (int m = 0; m < 3; ++m)
                                    dx[m] += ctx.jac[m][i] * du;
                            }
                            emit(y, w * scaled_norm(std::span<const T>(dx)), tag, dy);
                        }
                    }
                    return;
                }
            }
            nodes.push_back(L);
            nodes.push_back(U);
            std::sort(nodes.begin(), nodes.end());
            for (std::size_t s = 0; s + 1 < nodes.size(); ++s)
            {
                const T a = nodes[s], b = nodes[s + 1];
                // leaf cells keep degenerate segments so segment numbers stay aligned
                if (ctx.vis_order == 0 && b - a <= segment_tol<T>)
                    continue;
                tag.segment[D - 1] = static_cast<int>(s);
                // the roots of Bernstein forms on this line are all found: no level set
                // changes sign within a segment, its middle decides for all its points
                // (but for a sliver of rounding, whose points are tested one by one)
                bool whole = false;
                if constexpr (D == Top)
                    if (ctx.vis_order == 0 && !ctx.any_analytic && b - a > T(1e-12) * (U - L))
                    {
                        const T middle = T(0.5) * (a + b);
                        const VecD<T, D> y = insert<T, D>(yb, k, middle);
                        if (an.decided == 0 && !holds<T>(ctx, an.terms, to_u<T, D>(p.frame, y), (middle - L) / (U - L)))
                            continue;
                        whole = true;
                        if (ctx.measure_only)
                        {
                            emit(y, w * (b - a) * ctx.detj, tag, VecD<T, D>{});
                            continue;
                        }
                    }
                const auto node = [&](int j, T t, T wt)
                {
                    tag.node[D - 1] = j;
                    const VecD<T, D> y = insert<T, D>(yb, k, t);
                    if constexpr (D == Top)
                    {
                        if (ctx.vis_order == 0)
                        {
                            // top level: keep the points of the selected part
                            if (!whole && an.decided == 0
                                && !holds<T>(ctx, an.terms, to_u<T, D>(p.frame, y), (t - L) / (U - L)))
                                return;
                            emit(y, w * wt * ctx.detj, tag, VecD<T, D>{});
                            return;
                        }
                    }
                    emit(y, w * wt, tag, VecD<T, D>{});
                };
                if (folds_above.empty())
                    segment_nodes(ctx, a, b, node);
                else
                    anchored_nodes(ctx, a, b, folds_above, node);
            }
        };
        if constexpr (D == 2)
            integrate_line<T, Top == 2>(ctx, base, line, folds.empty() ? nullptr : &folds);
        else
            integrate_box<T, D - 1, false, Top>(ctx, std::move(base), line, 0, false);
    };

    // With several active bounds, the base splits into the regions where one
    // lower and one upper bound are active: each needs only the restrictions to
    // those two, which are fewer functions to certify together. Where several
    // level sets meet, or in a rotated frame (whose box faces are clips), that
    // keeps the levels below tractable; with one level set the joint base
    // bisects more and is a little more accurate.
    if constexpr (D == 2 && Top == 2)
    {
        if (split && ctx.partner >= 0)
        {
            // a 2D cell's curve part on a split box: the crossings of its two
            // level sets, in the region and the part, with weight 1 (merged
            // with those of neighbouring boxes once the cell is done)
            NodeTag tag;
            tag.box[0] = ctx.box_ids[0]++;
            int index = 0;
            for (const VecD<T, 2>& y : split_crossings)
            {
                std::array<bool, 2> on = {false, false};
                for (const Func<T, 2>& f : p.funcs)
                {
                    if (f.linear || (f.ls != ctx.surface && f.ls != ctx.partner))
                        continue;
                    VecD<T, 2> g;
                    func_gradient<T, 2>(ctx, f, y, g);
                    if (std::abs(value<T, 2>(ctx, f, y)) <= point_tol<T> * (std::abs(g[0]) + std::abs(g[1])))
                        on[f.ls == ctx.partner ? 1 : 0] = true;
                }
                bool in = on[0] && on[1];
                for (const Half<T, 2>& h : p.clips)
                    in = in && h.c[0] * y[0] + h.c[1] * y[1] <= h.d + point_tol<T>;
                if (!in || (an.decided == 0 && !holds<T>(ctx, an.terms, to_u<T, 2>(p.frame, y), T(0))))
                    continue;
                tag.segment[0] = index++;
                emit(y, T(1), tag, VecD<T, 2>{});
            }
            return;
        }
    }
    bool several = !p.surfaces.empty();
    for (std::size_t i = 1; i < p.funcs.size() && !several; ++i)
        several = !p.funcs[i].linear && !p.funcs[0].linear && p.funcs[i].ls != p.funcs[0].ls;
    for (std::size_t i = 0; i < p.funcs.size() && !several; ++i)
        for (std::size_t j = i + 1; j < p.funcs.size() && !several; ++j)
            several = !p.funcs[i].linear && !p.funcs[j].linear && p.funcs[i].ls != p.funcs[j].ls;
    if (!ctx.opt.split_bounds || !(several || Rotated) || lowers.size() * uppers.size() == 1)
    {
        integrate_base(lowers, uppers, {});
        return;
    }
    std::vector<Half<T, D - 1>> clips;
    for (std::size_t i = 0; i < lowers.size(); ++i)
        for (std::size_t j = 0; j < uppers.size(); ++j)
        {
            clips.clear();
            // lower_m <= lower_i and upper_j <= upper_m
            for (std::size_t m = 0; m < lowers.size(); ++m)
                if (m != i)
                {
                    Half<T, D - 1> h;
                    for (int jb = 0; jb < D - 1; ++jb)
                        h.c[jb] = lowers[m].beta[jb] - lowers[i].beta[jb];
                    h.d = lowers[i].alpha - lowers[m].alpha;
                    clips.push_back(h);
                }
            for (std::size_t m = 0; m < uppers.size(); ++m)
                if (m != j)
                {
                    Half<T, D - 1> h;
                    for (int jb = 0; jb < D - 1; ++jb)
                        h.c[jb] = uppers[j].beta[jb] - uppers[m].beta[jb];
                    h.d = uppers[m].alpha - uppers[j].alpha;
                    clips.push_back(h);
                }
            integrate_base({lowers[i]}, {uppers[j]}, clips);
        }
}

template <std::floating_point T, int D, bool Rotated, int Top, typename Emit>
void integrate_box(Context<T>& ctx, Problem<T, D> p, const Emit& emit, int depth, bool interface)
{
    if (ctx.emitted->n_points() - ctx.first_point > ctx.points_allowed)
        ctx.over_limit = true;
    if (ctx.over_limit || ctx.unmargined_failed || !tighten<T, D>(p.lo, p.hi, p.clips))
        return;

    // keep the functions that may vanish on the box, with their direction margins
    Analysis<T, D> an = analyse<T, D>(ctx, p, D == Top, depth);
    if constexpr (D == Top)
    {
        if (an.decided < 0)
            return; // the selection holds nowhere in the box
        if (an.decided > 0 && !interface)
        {
            // it holds on the whole box: no breakpoints needed
            an.funcs.clear();
            an.ratios.clear();
            an.uppers.clear();
            an.curved.clear();
            an.two_roots.clear();
            an.independent.clear();
            an.passive.clear();
            an.certified = true;
        }
    }
    else
    {
        if (interface && an.funcs.empty())
            return;
        // a curve part, on the base of the top box: nothing where its surface
        // function has no root, nor need the other functions certify there
        if constexpr (D == 2 && Top == 3)
            if (ctx.partner >= 0
                && std::none_of(an.surfaces.begin(), an.surfaces.end(),
                                [&ctx](const Surface<T, 2>& s) { return partner_surface<T, 2>(ctx, s); }))
                return;
    }
    p.funcs = std::move(an.funcs);
    p.surfaces = std::move(an.surfaces);

    // Options::silhouette_splits: as long as the base settles
    if constexpr (D == 3 && Top == 3)
    {
        if (an.unmargined)
        {
            const int box = trace_enter<T, D, Top>(ctx, p, an.k, depth, Rotated, BoxFate::certified);
            const int start = ctx.emitted->n_points(), bisections = ctx.bisections;
            ctx.level2_bisections = 0;
            ctx.unmargined_failed = false;
            reduce<T, D, Rotated, Top>(ctx, p, an, emit, depth, interface);
            const bool failed = ctx.unmargined_failed;
            ctx.unmargined_failed = false;
            ctx.top_unmargined = false;
            if (!failed)
            {
                trace_leave(ctx, box);
                return;
            }
            trace_fate(ctx, BoxFate::abandoned);
            ++ctx.abandoned;
            ctx.bisections = bisections;
            ctx.emitted->points.resize(3 * static_cast<std::size_t>(start));
            ctx.emitted->weights.resize(static_cast<std::size_t>(start));
            trace_leave(ctx, box);
            an.unmargined = false;
            an.certified = false;
            std::fill(an.two_roots.begin(), an.two_roots.end(), char(0));
        }
    }

    // Level 2: two functions whose zero curves meet (on a face of a tet, say) may
    // each need a different axis, and no bisection separates them. A diagonal
    // frame often suits both. Not at the top level of 2D cells, whose points
    // are the box's own coordinates.
    if constexpr (D == 2 && !Rotated && Top != 2)
    {
        if (!an.certified && level2_frames(ctx))
        {
            Problem<T, 2> r = rotate_diagonal<T>(p);
            if (tighten<T, 2>(r.lo, r.hi, r.clips))
            {
                Analysis<T, 2> ar = analyse<T, 2>(ctx, r, false, depth);
                if (ar.certified)
                {
                    ++ctx.stats->rotations;
                    r.funcs = std::move(ar.funcs);
                    r.surfaces = std::move(ar.surfaces);
                    // the frame's factor 1/2 is an area's; a curve's points get
                    // theirs from the derivative along it, which the map carries
                    const bool curve = ctx.partner >= 0;
                    const auto back = [&emit, curve](const VecD<T, 2>& z, T w, const NodeTag& tag, const VecD<T, 2>& dz)
                    {
                        emit(VecD<T, 2>{T(0.5) * (z[0] - z[1]), T(0.5) * (z[0] + z[1])}, curve ? w : T(0.5) * w, tag,
                             VecD<T, 2>{T(0.5) * (dz[0] - dz[1]), T(0.5) * (dz[0] + dz[1])});
                    };
                    const int region = trace_enter<T, 2, Top>(ctx, p, an.k, depth, Rotated, BoxFate::rotated);
                    const int box = trace_enter<T, 2, Top>(ctx, r, ar.k, depth, true, BoxFate::certified);
                    reduce<T, 2, true, Top>(ctx, std::move(r), ar, back, depth, interface);
                    trace_leave(ctx, box);
                    trace_leave(ctx, region);
                    return;
                }
            }
        }
    }
    // Where the zero sets of two level sets meet, their normals may differ so
    // much that no axis suits both: a height direction between them does.
    if constexpr (!Rotated)
    {
        // two level sets: curved functions of both, or surface functions
        bool two = !p.surfaces.empty();
        for (std::size_t c = 1; c < an.curved.size() && !two; ++c)
            two = p.funcs[static_cast<std::size_t>(an.curved[c])].ls
                  != p.funcs[static_cast<std::size_t>(an.curved[0])].ls;
        std::array<VecD<T, D>, D> R;
        // in 3D only once bisection has not helped and where the zero sets meet:
        // a rotated box has its faces as clip planes, which make the levels
        // below costly
        bool allowed = D < 3 ? level2_frames(ctx) : depth >= ctx.opt.rotation_depth;
        if constexpr (D == 3)
            allowed = allowed && !an.certified && two && zero_sets_meet<T>(ctx, p, an.curved);
        if (!an.certified && two && allowed && choose_frame<T, D>(ctx, p, an, R))
        {
            Problem<T, D> r = rotate_frame<T, D>(p, R);
            if (tighten<T, D>(r.lo, r.hi, r.clips))
            {
                Analysis<T, D> ar = analyse<T, D>(ctx, r, D == Top, depth);
                if constexpr (D == Top)
                {
                    // the rotated box fits the region differently: its bounds may decide it
                    if (ar.decided < 0)
                        return;
                    if (ar.decided > 0 && !interface)
                    {
                        ar.funcs.clear();
                        ar.ratios.clear();
                        ar.uppers.clear();
                        ar.curved.clear();
                        ar.two_roots.clear();
                        ar.independent.clear();
                        ar.passive.clear();
                        ar.certified = true;
                    }
                }
                if (ar.certified)
                {
                    ++ctx.stats->rotations;
                    r.funcs = std::move(ar.funcs);
                    r.surfaces = std::move(ar.surfaces);
                    const auto back = [&emit, R](const VecD<T, D>& z, T w, const NodeTag& tag, const VecD<T, D>& dz)
                    {
                        VecD<T, D> y{}, dy{};
                        for (int j = 0; j < D; ++j)
                            for (int m = 0; m < D; ++m)
                            {
                                y[j] += R[j][m] * z[m];
                                dy[j] += R[j][m] * dz[m];
                            }
                        emit(y, w, tag, dy);
                    };
                    const int region = trace_enter<T, D, Top>(ctx, p, an.k, depth, Rotated, BoxFate::rotated);
                    const int box = trace_enter<T, D, Top>(ctx, r, ar.k, depth, true, BoxFate::certified);
                    reduce<T, D, true, Top>(ctx, std::move(r), ar, back, depth, interface);
                    trace_leave(ctx, box);
                    trace_leave(ctx, region);
                    return;
                }
            }
        }
    }
    const int box = trace_enter<T, D, Top>(ctx, p, an.k, depth, Rotated, BoxFate::certified);
    reduce<T, D, Rotated, Top>(ctx, std::move(p), an, emit, depth, interface);
    trace_leave(ctx, box);
}

/// Whether a Bernstein form is affine: its coefficients are then its values at
/// the Greville points, c_ijk = c_000 + i/p_0 d_0 + j/p_1 d_1 + k/p_2 d_2.
template <std::floating_point T>
bool affine(const BoxBernstein<T>& f)
{
    const std::array<int, 3> p = f.degree;
    const std::array<int, 3> stride = {1, p[0] + 1, (p[0] + 1) * (p[1] + 1)};
    const T c0 = f.coeffs[0];
    std::array<T, 3> d{};
    for (int j = 0; j < 3; ++j)
        d[j] = p[j] > 0 ? f.coeffs[static_cast<std::size_t>(p[j] * stride[j])] - c0 : T(0);
    const T tol = compare_tol<T> * max_abs(std::span<const T>(f.coeffs));
    for (int k = 0; k <= p[2]; ++k)
        for (int j = 0; j <= p[1]; ++j)
            for (int i = 0; i <= p[0]; ++i)
            {
                T value = c0;
                if (p[0] > 0)
                    value += d[0] * T(i) / T(p[0]);
                if (p[1] > 0)
                    value += d[1] * T(j) / T(p[1]);
                if (p[2] > 0)
                    value += d[2] * T(k) / T(p[2]);
                if (std::abs(f.coeffs[static_cast<std::size_t>(i * stride[0] + j * stride[1] + k * stride[2])] - value)
                    > tol)
                    return false;
            }
    return true;
}

/// The top level: the unit box with the cell's clips and one function per
/// level set the terms name.
template <std::floating_point T, int Top>
void run_top(Context<T>& ctx, const ClippedBox<T>& cell, CellPoints<T>& out, bool interface)
{
    Problem<T, Top> top;
    top.lo.fill(T(0));
    top.hi.fill(T(1));
    for (int i = 0; i < Top; ++i)
        top.frame.A[i][i] = T(1);
    for (int l = 0; l < static_cast<int>(ctx.phis.size()); ++l)
    {
        if (!((ctx.used >> l) & 1))
            continue;
        Func<T, Top> f;
        f.ls = l;
        for (int i = 0; i < Top; ++i)
            f.A[i][i] = T(1);
        top.funcs.push_back(f);
    }
    for (const HalfSpace<T>& h : cell.clips)
    {
        Half<T, Top> c;
        for (int j = 0; j < Top; ++j)
            c.c[j] = h.c[j];
        c.d = h.d;
        top.clips.push_back(c);
    }
    const int vis_order = ctx.vis_order;
    const auto emit = [&out, vis_order](const VecD<T, Top>& y, T w, const NodeTag& tag, const VecD<T, Top>&)
    {
        const Vec3<T> u = pad<T, Top>(y);
        out.points.insert(out.points.end(), u.begin(), u.end());
        out.weights.push_back(w);
        if (vis_order > 0)
            out.tags.push_back(tag);
    };
    integrate_box<T, Top, false, Top>(ctx, std::move(top), emit, 0, interface);
}

/// The top level of an interval: the roots of its level sets in (0, 1) split
/// it into segments on which each has one sign. A volume part takes the
/// segments where a term holds at their midpoints; an interface part the
/// roots of its level set where a term holds, each a point of weight 1, as
/// the crossing points of 2D cells. A root at an end lies in a facet (a
/// vertex) of the cell, which the front end gives to one cell as a zero face.
template <std::floating_point T>
void run_interval(Context<T>& ctx, CellPoints<T>& out, bool interface)
{
    const Vec3<T> origin = {0, 0, 0}, direction = {1, 0, 0};
    std::vector<T>& nodes = ctx.nodes[1];
    nodes.clear();
    // the roots of level set l in (0, 1), appended to nodes
    const auto roots = [&](int l)
    {
        const Source<T>& phi = ctx.phis[static_cast<std::size_t>(l)];
        if (phi.is_analytic())
            line_roots(phi, std::span<const T>(origin), std::span<const T>(direction), T(0), T(1),
                       zero_function_tol<T> * ctx.scales[static_cast<std::size_t>(l)], nodes);
        else
            isolate_roots(std::span<const T>(phi.bernstein->coeffs), T(0), T(1), nodes, ctx.root_work[1]);
    };
    const auto emit = [&](T t, T w, const NodeTag& tag)
    {
        out.points.insert(out.points.end(), {t, T(0), T(0)});
        out.weights.push_back(w);
        if (ctx.vis_order > 0)
            out.tags.push_back(tag);
    };
    NodeTag tag;
    tag.box[0] = 0;
    if (interface)
    {
        roots(ctx.surface);
        std::sort(nodes.begin(), nodes.end());
        int index = 0;
        for (const T t : nodes)
        {
            if (!holds<T>(ctx, ctx.terms, Vec3<T>{t, 0, 0}, T(0)))
                continue;
            tag.segment[0] = index++;
            emit(t, ctx.vis_order > 0 ? T(0) : T(1), tag);
        }
        return;
    }
    nodes.push_back(T(0));
    nodes.push_back(T(1));
    for (int l = 0; l < static_cast<int>(ctx.phis.size()); ++l)
        if ((ctx.used >> l) & 1)
            roots(l);
    std::sort(nodes.begin(), nodes.end());
    int segment = 0;
    for (std::size_t s = 0; s + 1 < nodes.size(); ++s)
    {
        const T a = nodes[s], b = nodes[s + 1];
        if (b - a <= segment_tol<T> || !holds<T>(ctx, ctx.terms, Vec3<T>{T(0.5) * (a + b), 0, 0}, T(0)))
            continue;
        tag.segment[0] = segment++;
        segment_nodes(ctx, a, b,
                      [&](int j, T t, T w)
                      {
                          tag.node[0] = j;
                          emit(t, w * ctx.detj, tag);
                      });
    }
}

/// Points of @p out from @p start on within @p tol (box coordinates) of an
/// earlier one are dropped: a crossing point on the face two boxes share is
/// found by both.
template <std::floating_point T>
void merge_points(CellPoints<T>& out, int start, T tol)
{
    const bool tagged = out.tags.size() == out.weights.size();
    int kept = start;
    for (int i = start; i < out.n_points(); ++i)
    {
        bool duplicate = false;
        for (int j = start; j < kept && !duplicate; ++j)
        {
            T d = T(0);
            for (int c = 0; c < 3; ++c)
                d = std::max(d, std::abs(out.points[3 * i + c] - out.points[3 * j + c]));
            duplicate = d <= tol;
        }
        if (duplicate)
            continue;
        for (int c = 0; c < 3; ++c)
            out.points[3 * kept + c] = out.points[3 * i + c];
        out.weights[kept] = out.weights[i];
        if (tagged)
            out.tags[kept] = out.tags[i];
        ++kept;
    }
    out.points.resize(3 * static_cast<std::size_t>(kept));
    out.weights.resize(static_cast<std::size_t>(kept));
    if (tagged)
        out.tags.resize(static_cast<std::size_t>(kept));
}

/// The engine on one cell.
template <std::floating_point T>
void run(const ClippedBox<T>& cell, std::span<const Source<T>> phis, std::span<const SelectionTerm> terms, int q,
         int vis_order, const Options& opt, CellPoints<T>& out, Stats& stats, BoxTrace<T>* trace = nullptr)
{
    if (phis.empty() || phis.size() > 64)
        throw std::invalid_argument("quadrays: give 1 to 64 level sets");
    if (cell.tdim < 1 || cell.tdim > 3)
        throw std::invalid_argument("quadrays: cells of dimension 1, 2 or 3");
    for (const Source<T>& phi : phis)
    {
        if ((phi.bernstein == nullptr) == (phi.analytic == nullptr))
            throw std::invalid_argument("quadrays: a source is a Bernstein form or an analytic level set");
        if (phi.bernstein != nullptr && phi.bernstein->dim != cell.tdim)
            throw std::invalid_argument("quadrays: the level set must be a form in as many variables as the cell has");
        if (phi.analytic != nullptr && phi.tdim != cell.tdim)
            throw std::invalid_argument("quadrays: the analytic level set belongs to a cell of another dimension");
    }
    if (vis_order == 0 && q < 1)
        throw std::invalid_argument("quadrays: at least one Gauss point per segment is needed");
    if (opt.anchor_folds < 0 || opt.anchor_folds > 4)
        throw std::invalid_argument("quadrays: anchor_folds is 0 (off) to 4");

    thread_local Context<T> ctx;
    const std::size_t n = phis.size();
    const std::uint64_t all = n == 64 ? ~std::uint64_t(0) : (std::uint64_t(1) << n) - 1;
    ctx.terms.clear();
    ctx.used = 0;
    std::uint64_t zero = 0;
    bool first = true;
    for (const SelectionTerm& t : terms)
    {
        const Requirement r = {t.negative_required, t.positive_required, t.zero_required};
        if ((r.negative | r.positive | r.zero) & ~all)
            throw std::invalid_argument("quadrays: a term names a level set the cell does not have");
        if (std::popcount(r.zero) > 2)
            throw std::invalid_argument("quadrays: the points where three level sets vanish are not integrated");
        if (cell.tdim == 1 && std::popcount(r.zero) > 1)
            throw std::invalid_argument("quadrays: on intervals, the points where two level sets vanish are not "
                                        "integrated");
        if (!first && r.zero != zero)
            throw std::invalid_argument("quadrays: the terms ask for different zero sets");
        zero = r.zero;
        first = false;
        // a term that asks for two signs of one level set selects nothing
        if ((r.negative & r.positive) | (r.zero & (r.negative | r.positive)))
            continue;
        ctx.terms.push_back(r);
        ctx.used |= r.negative | r.positive | r.zero;
    }
    if (ctx.terms.empty())
        return;
    // a curve part (two zero sets; points on 2D cells) takes the roots of the
    // first on its height lines and finds the second on its zero set
    ctx.surface = zero != 0 ? std::countr_zero(zero) : -1;
    ctx.partner = std::popcount(zero) == 2 ? 63 - std::countl_zero(zero) : -1;
    ctx.phis.assign(phis.begin(), phis.end());
    for (PlaneForm<T>& plane : ctx.plane_forms)
        plane.axis = -1;
    ctx.measure_only = false;
    // An analytic level set that is linear on the cell (a plane) is read as its
    // form of degree 1, as a level set of degree 1 is: its zero set then needs
    // no surface function where it meets another's, whose bounds would come
    // from Taylor models.
    ctx.linear_forms.resize(n);
    for (std::size_t l = 0; l < n; ++l)
        if (((ctx.used >> l) & 1) && phis[l].is_analytic() && linear_form(phis[l], ctx.linear_forms[l]))
            ctx.phis[l] = bernstein_source(ctx.linear_forms[l]);
    if (vis_order == 0 && cell.tdim == 3 && ctx.partner >= 0 && opt.follow_curves
        && curve_points<T>(cell, std::span<const Source<T>>(ctx.phis.data(), n), terms, ctx.surface, ctx.partner, q, opt, out, stats))
        return;
    ctx.scales.resize(n);
    for (std::size_t l = 0; l < n; ++l)
        ctx.scales[l] = ((ctx.used >> l) & 1) ? reference_magnitude(ctx.phis[l]) : T(0);
    ctx.top_lines.resize(n);
    ctx.top_line_set.assign(n, 0);
    ctx.top_values.assign(n, T(0));
    ctx.opt = opt;
    // A curve part's check has no moment to catch the rule of one factor, and has
    // few anchored segments to save points on: it keeps two.
    if (ctx.partner >= 0 && ctx.opt.anchor_folds == 1)
        ctx.opt.anchor_folds = 2;
    const int anchor_folds = ctx.opt.anchor_folds;
    ctx.margin = static_cast<T>(opt.margin);
    ctx.stats = &stats;
    ctx.vis_order = vis_order;
    const bool interface = ctx.surface >= 0;
    ctx.detj = std::abs(jacobian_determinant(cell));
    ctx.inv = interface ? inverse_jacobian(cell) : Mat3<T>{};
    ctx.jac = cell.jacobian;
    ctx.box_ids = {0, 0, 0};
    ctx.bisections = 0;
    ctx.splits_blocked = 0;
    if (vis_order == 0 && ctx.q != q)
    {
        gauss_legendre<T>(q, ctx.gauss_x, ctx.gauss_w);
        ctx.q = q;
    }
    if (vis_order == 0 && anchor_folds > 0)
    {
        if (static_cast<int>(ctx.fold_x.size()) != anchor_folds * q)
            gauss_legendre<T>(anchor_folds * q, ctx.fold_x, ctx.fold_w);
        if (static_cast<int>(ctx.fold_check_x.size()) != anchor_folds * (q + 1))
            gauss_legendre<T>(anchor_folds * (q + 1), ctx.fold_check_x, ctx.fold_check_w);
        if (static_cast<int>(ctx.cubic_x.size()) != (3 * anchor_folds * q + 1) / 2)
            gauss_legendre<T>((3 * anchor_folds * q + 1) / 2, ctx.cubic_x, ctx.cubic_w);
        if (static_cast<int>(ctx.cubic_check_x.size()) != (3 * anchor_folds * (q + 1) + 1) / 2)
            gauss_legendre<T>((3 * anchor_folds * (q + 1) + 1) / 2, ctx.cubic_check_x, ctx.cubic_check_w);
        if (anchor_folds == 1)
        {
            if (static_cast<int>(ctx.more_fold_x.size()) != 2 * q)
                gauss_legendre<T>(2 * q, ctx.more_fold_x, ctx.more_fold_w);
            if (static_cast<int>(ctx.more_fold_check_x.size()) != 2 * (q + 1))
                gauss_legendre<T>(2 * (q + 1), ctx.more_fold_check_x, ctx.more_fold_check_w);
            if (static_cast<int>(ctx.more_cubic_x.size()) != 3 * q)
                gauss_legendre<T>(3 * q, ctx.more_cubic_x, ctx.more_cubic_w);
            if (static_cast<int>(ctx.more_cubic_check_x.size()) != 3 * (q + 1))
                gauss_legendre<T>(3 * (q + 1), ctx.more_cubic_check_x, ctx.more_cubic_check_w);
        }
    }
    ctx.more_points = false;
    ctx.flat.assign(n, 0);
    for (std::size_t l = 0; l < n; ++l)
        ctx.flat[l] = ((ctx.used >> l) & 1) && !ctx.phis[l].is_analytic() && affine(*ctx.phis[l].bernstein);
    ctx.out = nullptr;
    if (vis_order == 0 && opt.tolerance >= 0 && cell.tdim >= 2)
    {
        if (static_cast<int>(ctx.check_x.size()) != q + 1)
            gauss_legendre<T>(q + 1, ctx.check_x, ctx.check_w);
        ctx.out = &out;
        // h^m for a part of dimension m
        const int m = cell.tdim - std::popcount(zero);
        ctx.moment_degree = std::max(1, (2 * q - m) / 2);
        ctx.moment_relative = static_cast<T>(std::pow(opt.moment_tolerance, std::min(q, 6) / 6.0));
        ctx.ref_origin = cell.ref_origin;
        ctx.ref_jacobian = cell.ref_jacobian;
        ctx.cell_q = q;
        ctx.q_raised = 0;
        ctx.split_relative = static_cast<T>(std::pow(static_cast<double>(split_tolerance<T>), std::min(q, 6) / 6.0));
        ctx.measure_scale = std::pow(ctx.detj, T(m) / T(cell.tdim));
        ctx.cell_size = std::pow(ctx.detj, T(1) / T(cell.tdim));
        ctx.inv = inverse_jacobian(cell);
    }
    ctx.tdim = cell.tdim;
    ctx.derived.assign(3 * n, -1);
    ctx.derived_forms.clear();
    ctx.derived_directions.clear();
    ctx.derived_analytic.clear();
    ctx.silhouettes = opt.silhouette_splits && opt.fold_splits && opt.anchor_folds > 0 && vis_order == 0 && cell.tdim == 3;
    ctx.abandoned = 0;
    if (ctx.silhouettes)
        for (std::size_t l = 0; l < n; ++l)
        {
            if (!((ctx.used >> l) & 1) || ctx.flat[l])
                continue;
            if (ctx.phis[l].is_analytic())
            {
                if (ctx.phis[l].analytic->hessian_bounds == nullptr)
                    continue;
                for (int a = 0; a < 3; ++a)
                {
                    DirectionalDerivative& d = ctx.derived_directions.emplace_back();
                    d.phi = ctx.phis[l].analytic;
                    for (int i = 0; i < 3; ++i)
                        d.direction[i] = static_cast<double>(ctx.phis[l].jacobian[static_cast<std::size_t>(i)]
                                                                                  [static_cast<std::size_t>(a)]);
                    ctx.derived_analytic.push_back(directional_derivative(d));
                    Source<T> source = ctx.phis[l];
                    source.analytic = &ctx.derived_analytic.back();
                    ctx.derived[3 * l + static_cast<std::size_t>(a)] = static_cast<int>(ctx.phis.size());
                    ctx.phis.push_back(source);
                    ctx.scales.push_back(reference_magnitude(ctx.phis.back()));
                    ctx.flat.push_back(0);
                }
                continue;
            }
            for (int a = 0; a < 3; ++a)
            {
                if (ctx.phis[l].bernstein->degree[static_cast<std::size_t>(a)] < 2)
                    continue;
                derivative(*ctx.phis[l].bernstein, a, ctx.derived_forms.emplace_back());
                ctx.derived[3 * l + static_cast<std::size_t>(a)] = static_cast<int>(ctx.phis.size());
                ctx.phis.push_back(bernstein_source(ctx.derived_forms.back()));
                ctx.scales.push_back(reference_magnitude(ctx.phis.back()));
                ctx.flat.push_back(0);
            }
        }
    // after an analytic level set linear on the cell became its form of degree 1
    ctx.any_analytic = false;
    for (const Source<T>& phi : ctx.phis)
        ctx.any_analytic = ctx.any_analytic || phi.is_analytic();
    const int start = out.n_points();
    ctx.emitted = &out;
    ctx.first_point = start;
    ctx.points_allowed = vis_order == 0 ? point_limit(opt, q, cell.tdim) : std::numeric_limits<int>::max();
    ctx.over_limit = false;
    ctx.top_unmargined = false;
    ctx.top_clips = nullptr;
    ctx.unmargined_failed = false;
    ctx.moment_only = false;
    ctx.moment_splits = 0;
    ctx.moment_floor = T(0);
    ctx.measure_floor = T(0);
    ctx.trace = trace;
    ctx.trace_box = -1;
    const auto integrate_cell = [&]()
    {
        if (cell.tdim == 3)
            run_top<T, 3>(ctx, cell, out, interface);
        else if (cell.tdim == 2)
            run_top<T, 2>(ctx, cell, out, interface);
        else
            run_interval<T>(ctx, out, interface);
    };
    integrate_cell();
    if (ctx.over_limit)
    {
        ++stats.points_limited;
        // the trace keeps the run that went wrong
        ctx.trace = nullptr;
        for (const int depth : {6, 4, 2, 1, 0})
        {
            if (depth >= opt.max_depth)
                continue;
            out.points.resize(3 * static_cast<std::size_t>(start));
            out.weights.resize(static_cast<std::size_t>(start));
            ctx.opt.max_depth = depth;
            ctx.over_limit = false;
            ctx.bisections = 0;
            ctx.splits_blocked = 0;
            ctx.box_ids = {0, 0, 0};
            if (depth == 0)
                ctx.points_allowed = std::numeric_limits<int>::max();
            integrate_cell();
            if (!ctx.over_limit)
                break;
        }
    }
    if (ctx.bisections >= opt.max_bisections)
        ++stats.bisections_limited;
    if (ctx.partner >= 0 && cell.tdim == 2)
        merge_points(out, start, T(2) * point_tol<T>);
}

} // namespace

template <std::floating_point T>
void integrate(const ClippedBox<T>& cell, std::span<const Source<T>> phis, std::span<const SelectionTerm> terms,
               int q, const Options& opt, CellPoints<T>& out, Stats& stats, BoxTrace<T>* trace)
{
    run(cell, phis, terms, q, 0, opt, out, stats, trace);
}

template <std::floating_point T>
void integrate(const ClippedBox<T>& cell, const Source<T>& phi, Part part, int q, const Options& opt,
               CellPoints<T>& out, Stats& stats)
{
    const SelectionTerm term = term_of(part);
    run(cell, std::span<const Source<T>>(&phi, 1), std::span<const SelectionTerm>(&term, 1), q, 0, opt, out, stats);
}

template <std::floating_point T>
void integrate(const ClippedBox<T>& cell, const BoxBernstein<T>& phi, Part part, int q,
               const Options& opt, CellPoints<T>& out, Stats& stats)
{
    integrate(cell, bernstein_source(phi), part, q, opt, out, stats);
}

template <std::floating_point T>
void leaf_nodes(const ClippedBox<T>& cell, std::span<const Source<T>> phis, std::span<const SelectionTerm> terms,
                int degree, const Options& opt, CellPoints<T>& out, Stats& stats)
{
    if (degree < 1)
        throw std::invalid_argument("quadrays: leaf cells need degree 1 or more");
    run(cell, phis, terms, 0, degree, opt, out, stats);
}

template <std::floating_point T>
void leaf_nodes(const ClippedBox<T>& cell, const Source<T>& phi, Part part, int degree, const Options& opt,
                CellPoints<T>& out, Stats& stats)
{
    const SelectionTerm term = term_of(part);
    leaf_nodes(cell, std::span<const Source<T>>(&phi, 1), std::span<const SelectionTerm>(&term, 1), degree, opt, out,
               stats);
}

template <std::floating_point T>
void leaf_nodes(const ClippedBox<T>& cell, const BoxBernstein<T>& phi, Part part, int degree,
                const Options& opt, CellPoints<T>& out, Stats& stats)
{
    leaf_nodes(cell, bernstein_source(phi), part, degree, opt, out, stats);
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template void integrate<float>(const ClippedBox<float>&, std::span<const Source<float>>,
                               std::span<const SelectionTerm>, int, const Options&, CellPoints<float>&, Stats&,
                               BoxTrace<float>*);
template void integrate<double>(const ClippedBox<double>&, std::span<const Source<double>>,
                               std::span<const SelectionTerm>, int, const Options&, CellPoints<double>&, Stats&,
                               BoxTrace<double>*);
template void integrate<float>(const ClippedBox<float>&, const Source<float>&, Part, int, const Options&,
                               CellPoints<float>&, Stats&);
template void integrate<double>(const ClippedBox<double>&, const Source<double>&, Part, int, const Options&,
                                CellPoints<double>&, Stats&);
template void integrate<float>(const ClippedBox<float>&, const BoxBernstein<float>&, Part, int,
                               const Options&, CellPoints<float>&, Stats&);
template void integrate<double>(const ClippedBox<double>&, const BoxBernstein<double>&, Part, int,
                                const Options&, CellPoints<double>&, Stats&);
template void leaf_nodes<float>(const ClippedBox<float>&, std::span<const Source<float>>,
                                std::span<const SelectionTerm>, int, const Options&, CellPoints<float>&, Stats&);
template void leaf_nodes<double>(const ClippedBox<double>&, std::span<const Source<double>>,
                                 std::span<const SelectionTerm>, int, const Options&, CellPoints<double>&, Stats&);
template void leaf_nodes<float>(const ClippedBox<float>&, const Source<float>&, Part, int, const Options&,
                                CellPoints<float>&, Stats&);
template void leaf_nodes<double>(const ClippedBox<double>&, const Source<double>&, Part, int, const Options&,
                                 CellPoints<double>&, Stats&);
template void leaf_nodes<float>(const ClippedBox<float>&, const BoxBernstein<float>&, Part, int,
                                const Options&, CellPoints<float>&, Stats&);
template void leaf_nodes<double>(const ClippedBox<double>&, const BoxBernstein<double>&, Part, int,
                                 const Options&, CellPoints<double>&, Stats&);

} // namespace cutcells::quadrays
