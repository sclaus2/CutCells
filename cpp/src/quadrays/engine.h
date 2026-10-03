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
/// or, for analytic level sets, Taylor models; otherwise the box is bisected.
struct Options
{
    /// Required |d_k psi| / |grad psi| on a box. It also controls the accuracy:
    /// larger margins bisect more and integrate more accurately.
    double margin = 0.25;
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
    /// Drop bounds of the height lines that are never active on the base region.
    bool prune_bounds = true;
    /// Split the base where the active bounds of the height lines change, so
    /// that each region carries only the restrictions to its two bounds.
    bool split_bounds = true;
    /// Level 2: if no axis certifies, try the diagonal frame before bisecting.
    /// On a tet, phi on a box face and phi on the slanted face can need
    /// different axes where their zero curves meet; no bisection separates them.
    bool diagonal_frames = true;
    /// Where the zero sets of two level sets meet, their normals may differ so
    /// much that no axis suits both, and bisection does not help. Such boxes
    /// are tried in a frame rotated between the normals: at level 2 always,
    /// in 3D from this depth of bisection on (the levels below a rotated box
    /// get its faces as oblique clip planes, which costs). Above max_depth:
    /// never in 3D.
    int rotation_depth = 6;
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
    /// With Options::diagnose: bisections counted by level, the function that
    /// blocked certification, and what sampling the clipped region says.
    std::map<std::string, std::int64_t> causes;
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
    std::vector<T> points;     ///< box coordinates, 3 per point (u2 = 0 on 2D cells)
    std::vector<T> weights;    ///< physical weights (0 for leaf nodes)
    std::vector<NodeTag> tags; ///< leaf nodes only

    int n_points() const { return static_cast<int>(weights.size()); }
};

/// @brief Quadrature points of the part of a cell that @p terms select,
/// appended to @p out.
///
/// The part is where some term holds; bit l of a term's masks refers to
/// phis[l]. All terms ask for the zero set of the same level set (an interface
/// part; weights are surface measures) or of none (a volume part). Curves
/// where two level sets vanish are not integrated.
///
/// @param cell   the cell as a clipped box
/// @param phis   the cell's level sets: Bernstein forms on its box or analytic
///               level sets (bernstein_source, analytic_source), at most 64
/// @param terms  the selection, an or of terms
/// @param q      Gauss-Legendre points per segment of each height line
/// @param opt    engine options
/// @param out    points in box coordinates and physical weights
/// @param stats  counters, accumulated
/// @throws std::invalid_argument if terms ask for different zero sets, or for
///         two zero sets at once
template <std::floating_point T>
void integrate(const ClippedBox<T>& cell, std::span<const Source<T>> phis, std::span<const SelectionTerm> terms,
               int q, const Options& opt, CellPoints<T>& out, Stats& stats);

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
/// Every segment gets degree + 1 equispaced nodes, pulled 1e-5 of its length
/// inside its ends, with their tags; weights are 0. Leaves are not filtered by
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
