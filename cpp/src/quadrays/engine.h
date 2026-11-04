// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <concepts>
#include <cstdint>
#include <map>
#include <span>
#include <string>
#include <vector>

#include "../selection_expr.h"
#include "box_bernstein.h"
#include "clipped_box.h"
#include "source.h"

namespace cutcells::quadrays
{

/// Options of the quadrays engine.
///
/// The engine reduces the dimension one height direction at a time, as in
/// Saye's 2015 algorithm, with the clip planes of the cell handled directly.
/// At each level the functions are the level sets restricted to a flat,
/// phi(A y + b), linear functions, and where two level sets cut a box, one
/// level set restricted to the other's zero set along the height lines (a
/// surface function, which vanishes where their roots cross). A height
/// direction k is accepted only if every curved function satisfies
/// |d_k psi| >= margin |grad psi| on the level's box, from Bernstein bounds
/// or, for analytic level sets, Taylor models; otherwise a box with two free
/// coordinates is split at the folds and crossings of its zero curves
/// (fold_splits), and other boxes are bisected.
struct Options
{
    /// Required |d_k psi| / |grad psi| on a box. Larger margins bisect more.
    /// A box needs a direction with this margin for every function on it, so
    /// margins above 1/sqrt(D) on boxes with D free coordinates cannot be met
    /// where the normal points between the axes: the engine caps the margin of
    /// each level at 0.9/sqrt(D).
    double margin = 0.25;
    /// Relative accuracy of curved pieces. Where a box with two free coordinates
    /// (a 2D cell's box, or the base of a 3D cell's box) carries a curved zero
    /// set, its piece is integrated with q and q + 1 points per segment; while
    /// the two differ by more than this times its measure, the box is bisected,
    /// and the q points are kept. The margin bounds the slope of the height
    /// functions, not how close their branch points (where a zero set turns
    /// along the height direction) come to a box: where a box face cuts a
    /// sphere near its pole, say, the zero set on the face is a circle of any
    /// radius, whose branch points can lie just beyond the base segments, and
    /// the points then converge slowly in q and not under mesh refinement.
    /// anchor_folds anchors the segments near them; the check bisects where
    /// the points still converge slowly.
    /// 0: the rule's own rate, (h kappa / 2)^(2 q) with h = |det J|^(1/tdim)
    /// the cell's size and kappa the curvature of the zero sets in the box, so
    /// that the error falls as h^(2q) and pieces that converge as expected
    /// are not bisected. < 0: no check. Flat pieces are not checked.
    double tolerance = 0;
    /// On a box with two free coordinates, where the zero set of a curved
    /// function turns back along the height direction (a fold, where
    /// d_k psi = 0) at or within four segment lengths beyond an end of a base
    /// segment, its height function has a square-root branch point there, and
    /// Gauss-Legendre points in x converge slowly. n > 0: such segments take
    /// n q Gauss-Legendre points in s, x = x_f -+ s^2 anchored at the fold x_f,
    /// which makes the integrand analytic. With folds beyond both ends, one of
    /// them anchors the segment if its points converge at the other as fast as
    /// q points do at four lengths; otherwise both, by x = c + r (3t - t^3) / 2,
    /// whose derivative vanishes at both folds, with ceil(3 n q / 2) points.
    /// They integrate polynomials of degree n q - 1 in x exactly (q points:
    /// 2 q - 1), so n = 2 keeps the degree, and n = 1 halves the rate at which
    /// the error falls under refinement. Below a 3D box, the zero set of its
    /// level set turns back along its height direction where its normal
    /// turns perpendicular to it, beyond the box (its silhouette): the height
    /// lines of the base boxes are anchored likewise at the silhouette's
    /// coordinate on them, and their base segments where it crosses the
    /// lines' bounds. 0: off.
    int anchor_folds = 2;
    /// Bisections allowed per level along a branch. Bernstein bounds certify
    /// after a few; Taylor models of an analytic level set, a distance function
    /// in particular, need boxes small against the radius of curvature, which
    /// for a feature half a cell wide takes about 10. (The prototype allowed 8.)
    int max_depth = 12;
    /// Bisections allowed per cell, all levels and branches together. Beyond it
    /// boxes are integrated uncertified (roots isolated in full), so the cost
    /// stays bounded where certification cannot succeed (two sheets of the
    /// level set closer than the depth limit resolves, singular points). (The
    /// prototype allowed 256.)
    int max_bisections = 1024;
    /// Points allowed in a cell, per selection group. A cell that exceeds them
    /// is integrated again with at most 6, 4, 2, 1 or 0 bisections per level,
    /// the first that stays below, which costs accuracy; the front end warns.
    /// 0: 1000 q^tdim for q points per segment, at most 300000. < 0: no limit.
    int max_points = 0;
    /// Drop bounds of the height lines that are never active on the base region.
    bool prune_bounds = true;
    /// Split the base where the active bounds of the height lines change, so
    /// that each region carries only the restrictions to its two bounds.
    bool split_bounds = true;
    /// Level 2, where boxes are not split at their folds and crossings
    /// (fold_splits off, anchor_folds 0, or leaf nodes): if no axis certifies, try
    /// the diagonal frame and a frame rotated between the normals of the zero
    /// curves before bisecting. On a tet, phi on a box face and phi on the slanted
    /// face can need different axes where their zero curves meet; no bisection
    /// separates them. With fold splits these frames gain nothing and cost up to
    /// 40% of the points (a ball, 5 cells per direction), so they are not tried.
    bool diagonal_frames = true;
    /// Where the zero sets of two level sets meet in 3D, their normals may differ
    /// so much that no axis suits both, and bisection does not help. Boxes from
    /// this depth of bisection on, in which both zero sets change sign in one
    /// sub-box, are tried in a frame rotated between the normals; the box's faces
    /// become oblique clip planes, so that the rule covers the box exactly, and
    /// each such box costs about nine times the points of an axis-aligned one.
    /// The default is max_depth: a last resort before a box goes uncertified.
    /// Above max_depth: never.
    int rotation_depth = 12;
    /// M > 1: margins from M^D sub-cells that meet the clipped region and on
    /// which the function may vanish (enough for at most one root per height
    /// line). Fewer bisections, but less accurate at the same margin.
    /// 1: Bernstein bounds on the whole box. Analytic level sets ignore it.
    int mask_subdivisions = 1;
    /// Analytic level sets: Taylor models over the M^D sub-boxes that meet the
    /// clipped region, on boxes that a clip plane cuts (tetrahedra, prisms,
    /// pyramids, triangles); remainders shrink with the square of the size.
    /// 1: one model over the whole box.
    int taylor_subdivisions = 2;
    /// From this depth of bisection on, where no direction gives one root per
    /// height line, accept a direction in which the derivative along it is
    /// monotone and the two roots it allows never merge in the box (two sheets
    /// of one level set, a thin shell). Needs second derivatives: Bernstein
    /// forms, or analytic level sets with hessian_bounds. Earlier, bisection
    /// separates the roots and keeps the height functions flat; above
    /// max_depth, never.
    int two_roots_depth = 0;
    /// Record why each bisection happened in Stats::causes.
    bool diagnose = false;
    /// The front end (part::quadrature_rules), for studying cells that take
    /// many points. Not empty: each cell that reaches a limit or exceeds
    /// debug_points is written to this directory as a JSON file (its level sets
    /// on a grid, its boxes and points; write_cell_report), at most 16 per call,
    /// which `python -m cutcells.quadrays_debug` plots. Empty: nothing is
    /// recorded, at the cost of one test per box.
    std::string debug_dir;
    /// Points from which a cell is written to debug_dir.
    int debug_points = 10000;
    /// A box with two free coordinates that no axis certifies is not
    /// bisected: along the axis with the best sampled margins, the points where
    /// a zero curve turns along it (a fold, f = d_k f = 0) and where two zero
    /// curves cross (f = g = 0) are found by subdivision with bounds and
    /// Newton's method, the base line splits at them, the height lines take
    /// all roots, and segments anchor at the folds (as in Saye 2022, which
    /// removes these points from the base). Needs anchor_folds > 0, whose rule
    /// takes the square-root branch of the height functions at a fold at the
    /// end of a segment. Such boxes are larger than certified ones: with
    /// tolerance 0, their check asks q and q + 1 points to agree to 1e-12 and
    /// bisects otherwise. Where the search does not settle (a singular point,
    /// nearly tangent curves), the box and the boxes within are bisected.
    bool fold_splits = true;
    /// A 3D box that no axis certifies is not bisected where, along an axis, each
    /// of its level sets is monotone or has a monotone derivative: its lines take
    /// all roots, at most two of each per line, and the base splits at the
    /// silhouettes, where a zero set turns along the axis (f = d_k f = 0; Beck and
    /// Kummer 2023, Section 5.4, read d_k f as a second level set the same way), and
    /// at the curves where two zero sets meet. Segments anchor there (anchor_folds)
    /// and level-2 boxes split at these curves' folds and crossings (fold_splits);
    /// needs both. If the base does not settle, the box is bisected after all. For
    /// now: Bernstein level sets and no clip planes; other boxes bisect.
    bool silhouette_splits = true;
    /// The front end (part::quadrature_rules): Gauss-Legendre points per
    /// segment of each height line. 0: enough for its order, ceil((order + m)
    /// / 2) for a part of dimension m, so that flat pieces integrate
    /// polynomials of degree order exactly. The engine's own functions take
    /// the points as an argument and ignore it.
    int points_per_segment = 0;
};

/// Counters of one or more engine runs.
struct Stats
{
    int bisections = 0;        ///< boxes bisected, all levels
    int uncertified = 0;       ///< boxes integrated without a certified direction
    int rotations = 0;         ///< level-2 boxes integrated in the diagonal frame
    int incomplete_leaves = 0; ///< leaves of the part dropped because their nodes did not line up
    int two_roots = 0;         ///< boxes certified with two roots per height line
    int surfaces = 0;          ///< surface functions made where two level sets cut a box
    int fold_splits = 0;       ///< boxes split at folds and crossings (Options::fold_splits)
    int silhouettes = 0;       ///< 3D boxes integrated without certificate (Options::silhouette_splits)
    int points_limited = 0;    ///< runs that exceeded Options::max_points and were redone with fewer bisections
    int bisections_limited = 0; ///< runs that used all of Options::max_bisections
    /// top boxes of a curve part integrated uncertified (depth or bisection
    /// limit): no surface function finds the curve there, so its piece, if
    /// any, is missing
    int curve_lost = 0;
    /// With Options::diagnose: bisections counted by level, the function that
    /// blocked certification, and what sampling the clipped region says.
    std::map<std::string, std::int64_t> causes;
};

/// @brief The points Options::max_points allows in a cell of dimension @p tdim
/// with @p q points per segment; the largest int if none.
int point_limit(const Options& opt, int q, int tdim);

/// What became of a box in the engine's recursion (BoxTrace).
enum class BoxFate : std::uint8_t
{
    certified = 0,   ///< a height direction holds on all of it
    split = 1,       ///< split at folds and crossings (Options::fold_splits)
    bisected = 2,    ///< no direction holds on all of it: bisected
    refined = 3,     ///< certified or split, but its accuracy check failed: bisected
    rotated = 4,     ///< no axis holds: its one child is the box in a rotated frame
    uncertified = 5, ///< no direction holds and a limit is reached: roots isolated in full
    abandoned = 6    ///< tried without a margin, its base did not settle: its points were dropped
};

/// The boxes of the engine's recursion in a cell, in the order it visits them,
/// for studying cells that take many points. Boxes without a point of the part
/// that the bounds decide at once are not entered.
template <std::floating_point T>
struct BoxTrace
{
    std::vector<int> parent;       ///< the box it came from; -1: the cell
    std::vector<int> level;        ///< free coordinates
    std::vector<int> depth;        ///< bisections above it on its level
    std::vector<int> height;       ///< its height direction
    std::vector<BoxFate> fate;
    std::vector<std::uint8_t> rotated; ///< in a rotated frame: its faces are not the cell's
    std::vector<T> lo, hi;         ///< 3 per box, in the frame's coordinates (0 for unused ones)
    /// 12 per box, top boxes: the frame u = A y + b, A row-wise, then b; else 0.
    std::vector<T> frame;
    std::vector<int> points;       ///< points emitted in it and below it

    int n_boxes() const { return static_cast<int>(parent.size()); }
};

/// What a selection term asks of one level set.
enum class Part : std::uint8_t
{
    negative = 0,  ///< phi < 0
    positive = 1,  ///< phi > 0
    interface = 2, ///< phi = 0
    whole = 3      ///< no condition
};

/// @brief The part of level set @p level_set that a compiled selection term selects.
/// @throws std::invalid_argument if the term constrains another level set or
///         asks for contradictory signs.
Part part_of(const SelectionTerm& term, int level_set = 0);

/// @brief The selection term of a part of level set @p level_set.
SelectionTerm term_of(Part part, int level_set = 0);

/// Where an emitted point sits in the decomposition. Per level (index = number
/// of free coordinates - 1): the certified box, the segment along its height
/// line and the node within the segment. Leaf cells are assembled from it;
/// on 2D cells the third entries stay unused.
struct NodeTag
{
    std::array<int, 3> box = {-1, -1, -1};
    std::array<int, 3> segment = {0, 0, 0};
    std::array<int, 3> node = {0, 0, 0};
};

/// Points produced by the engine for one cell.
template <std::floating_point T>
struct CellPoints
{
    std::vector<T> points;     ///< box coordinates, 3 per point (u_j = 0 for j >= tdim)
    std::vector<T> weights;    ///< physical weights (0 for leaf nodes)
    std::vector<NodeTag> tags; ///< leaf nodes only

    int n_points() const { return static_cast<int>(weights.size()); }
};

/// @brief Quadrature points of the part of a cell that @p terms select,
/// appended to @p out.
///
/// The part is where some term holds; bit l of a term's masks refers to
/// phis[l]. All terms ask for the same zero sets: of none (a volume part), of
/// one level set (an interface part; weights are surface measures), or of two
/// (a curve part). A curve's points lie on the base lines of the boxes where
/// both level sets vanish, Gauss-Legendre in the line's coordinate, weighted by
/// the physical speed along the curve; on 2D cells the curve is the points
/// where the two zero sets cross, each with weight 1. On an interval, the
/// roots of the level sets split the cell into segments, and an interface is
/// the points where its level set vanishes, each with weight 1.
///
/// Weights are measures in the cell itself: on a facet in 3D, areas and
/// lengths on the facet. An analytic level set's derivative normal to the
/// cell does not count.
///
/// @param cell   the cell as a clipped box (intervals, 2D and 3D cells, in
///               any space of their dimension or higher)
/// @param phis   the cell's level sets: Bernstein forms on its box or analytic
///               level sets (bernstein_source, analytic_source), at most 64
/// @param terms  the selection, an or of terms
/// @param q      Gauss-Legendre points per segment of each height line
/// @param opt    engine options
/// @param out    points in box coordinates and physical weights
/// @param stats  counters, accumulated
/// @param trace  not null: the boxes of the recursion are appended (of the run
///               that exceeded Options::max_points, if one did)
/// @throws std::invalid_argument if terms ask for different zero sets, for
///         three zero sets at once, or for two on an interval
template <std::floating_point T>
void integrate(const ClippedBox<T>& cell, std::span<const Source<T>> phis, std::span<const SelectionTerm> terms,
               int q, const Options& opt, CellPoints<T>& out, Stats& stats, BoxTrace<T>* trace = nullptr);

/// @brief Quadrature points of one part of a cell with one level set,
/// appended to @p out.
///
/// @param phi    the cell's level set: a Bernstein form on its box or an
///               analytic level set (bernstein_source, analytic_source)
/// @param part   the selected part
template <std::floating_point T>
void integrate(const ClippedBox<T>& cell, const Source<T>& phi, Part part, int q, const Options& opt,
               CellPoints<T>& out, Stats& stats);

/// @brief integrate for a Bernstein form on the cell's box (cell_bernstein_on_box).
template <std::floating_point T>
void integrate(const ClippedBox<T>& cell, const BoxBernstein<T>& phi, Part part, int q,
               const Options& opt, CellPoints<T>& out, Stats& stats);

/// @brief Nodes of the leaf cells of the same decomposition, appended to @p out.
///
/// The decomposition is the certified one: the bisections of the accuracy
/// checks (Options::tolerance), which need Gauss points, are not made.
/// Every segment gets degree + 1 equispaced nodes, pulled 1e-5 of its length
/// inside its ends, with their tags (a curve part's: those of its base lines;
/// on 2D cells one node per crossing point, on intervals one per point of an
/// interface); weights are 0. Leaves are not filtered by
/// the terms where a box needs the level sets to decide: a leaf lies on one
/// side of each level set as a whole. The height lines below the top do not
/// split where a surface function vanishes beyond the box above, nor, in an
/// interface part, at the roots of the other level sets' restrictions: neither
/// is an edge of a leaf, and such roots may cross others inside a base
/// segment, so that the segments of a leaf would not line up.
template <std::floating_point T>
void leaf_nodes(const ClippedBox<T>& cell, std::span<const Source<T>> phis, std::span<const SelectionTerm> terms,
                int degree, const Options& opt, CellPoints<T>& out, Stats& stats);

/// @brief leaf_nodes for one level set.
template <std::floating_point T>
void leaf_nodes(const ClippedBox<T>& cell, const Source<T>& phi, Part part, int degree, const Options& opt,
                CellPoints<T>& out, Stats& stats);

/// @brief leaf_nodes for a Bernstein form on the cell's box.
template <std::floating_point T>
void leaf_nodes(const ClippedBox<T>& cell, const BoxBernstein<T>& phi, Part part, int degree,
                const Options& opt, CellPoints<T>& out, Stats& stats);

} // namespace cutcells::quadrays
