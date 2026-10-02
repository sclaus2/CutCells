// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <concepts>
#include <cstdint>
#include <map>
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
/// At each level the functions are the level set restricted to a flat,
/// phi(A y + b), or linear functions. A height direction k is accepted only if
/// every curved function satisfies |d_k psi| >= margin |grad psi| on the
/// level's box, from Bernstein bounds or, for analytic level sets, Taylor
/// models; otherwise the box is bisected.
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
    /// Level 2: if no axis certifies, try the diagonal frame before bisecting.
    /// On a tet, phi on a box face and phi on the slanted face can need
    /// different axes where their zero curves meet; no bisection separates them.
    bool diagonal_frames = true;
    /// M > 1: margins from M^D sub-cells that meet the clipped region and on
    /// which the function may vanish (enough for at most one root per height
    /// line). Fewer bisections, but less accurate at the same margin.
    /// 1: Bernstein bounds on the whole box. Analytic level sets ignore it.
    int mask_subdivisions = 1;
    /// Record why each bisection happened in Stats::causes.
    bool diagnose = false;
};

/// Counters of one or more engine runs.
struct Stats
{
    int bisections = 0;        ///< boxes bisected, all levels
    int uncertified = 0;       ///< boxes integrated without a certified direction
    int rotations = 0;         ///< level-2 boxes integrated in the diagonal frame
    int incomplete_leaves = 0; ///< leaves dropped because their nodes did not line up
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

/// Where an emitted point sits in the decomposition. Per level (index = number
/// of free coordinates - 1): the certified box, the segment along its height
/// line and the node within the segment. Leaf cells are assembled from it.
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
    std::vector<T> points;     ///< box coordinates, 3 per point
    std::vector<T> weights;    ///< physical weights (0 for leaf nodes)
    std::vector<NodeTag> tags; ///< leaf nodes only

    int n_points() const { return static_cast<int>(weights.size()); }
};

/// @brief Quadrature points of one part of a cell, appended to @p out.
///
/// @param cell   the cell as a clipped box
/// @param phi    the cell's level set: a Bernstein form on its box or an
///               analytic level set (bernstein_source, analytic_source)
/// @param part   the selected part
/// @param q      Gauss-Legendre points per segment of each height line
/// @param opt    engine options
/// @param out    points in box coordinates and physical weights
/// @param stats  counters, accumulated
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
/// inside its ends, with their tags; weights are 0. Volume parts are not
/// filtered by sign: a leaf lies on one side of the level set as a whole.
template <std::floating_point T>
void leaf_nodes(const ClippedBox<T>& cell, const Source<T>& phi, Part part, int degree, const Options& opt,
                CellPoints<T>& out, Stats& stats);

/// @brief leaf_nodes for a Bernstein form on the cell's box.
template <std::floating_point T>
void leaf_nodes(const ClippedBox<T>& cell, const BoxBernstein<T>& phi, Part part, int degree,
                const Options& opt, CellPoints<T>& out, Stats& stats);

} // namespace cutcells::quadrays
