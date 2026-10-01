// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <map>
#include <string>

#include "clipped_box.h"
#include "generators.h"
#include "leaf_mesh.h"
#include "selection_expr.h"

namespace cutcells::proto
{

/// Certify-and-bisect quadrature on a clipped box (one level set).
///
/// Dimension reduction as in Saye's 2015 algorithm, with clip planes handled
/// directly. At each level (3, 2, 1 free coordinates) the functions are the
/// level set restricted to a flat, phi(A y + b), or linear functions. A height
/// direction k is accepted only if every curved function satisfies
/// |d_k psi| >= margin * |grad psi| on the level's box, from Bernstein bounds;
/// otherwise the box is bisected. Clip planes and box faces bound the height
/// lines; the base level gets their Fourier-Motzkin projections as clips, the
/// restrictions of each function to each bounding plane (an affine
/// substitution, no resultants) and the linear functions where the active bound
/// changes. Each level's box is first shrunk to the bounding box of its clipped
/// region, so decisions only see the cell.
struct CertifyOptions
{
    double margin = 0.25; ///< required |d_k psi| / |grad psi| on the box
    int max_depth = 8;    ///< bisections allowed per level along a branch
    bool prune_bounds = true; ///< drop bounds that are never active on the base region
    bool diagnose = false;    ///< record why each bisection happened (CertifyStats::causes)
    /// Level 2: if no axis certifies, try the diagonal frame before bisecting. On a
    /// tet, phi on a box face and phi on the slanted face can need different axes
    /// where their zero curves meet on the shared edge; no bisection separates them.
    bool diagonal_frames = true;
    /// M > 1: margins from M^D sub-cells, counting only those that meet the clipped
    /// region and on which the function may vanish (enough for at most one root per
    /// height line). Fewer bisections, but less accurate at the same margin: the box
    /// bounds also act as the accuracy control. 1 = Bernstein bounds on the whole box.
    int mask_subdivisions = 1;
};

struct CertifyStats
{
    int bisections = 0;        ///< boxes bisected, all levels
    int uncertified = 0;       ///< boxes integrated without a certified direction (depth limit)
    int incomplete_leaves = 0; ///< leaves dropped because their nodes did not line up
    int rotations = 0;         ///< level-2 boxes integrated in the diagonal frame
    /// With CertifyOptions::diagnose: bisections counted by level, the origin of the
    /// function that blocked certification, and what sampling the clipped region says.
    std::map<std::string, long> causes;
};

void certified_bisection(const ClippedBox& cell, const LevelSet& ls, const SelectionTerm& term, int q,
                         const CertifyOptions& opt, Rule& rule, CertifyStats& stats);

/// Leaf cells of the same decomposition, for visualisation: every piece that
/// certified_bisection integrates, as a Lagrange cell of the given degree
/// (hexahedra for volume parts, quadrilaterals for the interface) whose nodes are
/// the images of an equispaced grid under the piece's parametrisation. Appended to
/// @p mesh with @p parent as background cell.
void certified_leaves(const ClippedBox& cell, const LevelSet& ls, const SelectionTerm& term, int degree,
                      const CertifyOptions& opt, std::int32_t parent, LeafMesh& mesh, CertifyStats& stats);

} // namespace cutcells::proto
