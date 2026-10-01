// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include "clipped_box.h"
#include "generators.h"
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
};

struct CertifyStats
{
    int bisections = 0;  ///< boxes bisected, all levels
    int uncertified = 0; ///< boxes integrated without a certified direction (depth limit)
};

void certified_bisection(const ClippedBox& cell, const LevelSet& ls, const SelectionTerm& term, int q,
                         const CertifyOptions& opt, Rule& rule, CertifyStats& stats);

} // namespace cutcells::proto
