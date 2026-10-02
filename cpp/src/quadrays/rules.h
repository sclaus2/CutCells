// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

/// quadrays: height-function quadrature on cut cells.
///
/// Every cell is a box clipped by half-spaces; a tetrahedron is the unit box
/// with u0 + u1 + u2 <= 1.  The engine reduces the dimension one height
/// direction at a time, accepts a direction only where bounds of the level
/// set certify it, and bisects the box otherwise.
///
/// This header is the module's entry point: quadrature rules for one cell and
/// one selection term, with points in the parent cell's reference
/// coordinates and physical weights.
namespace cutcells::quadrays
{
} // namespace cutcells::quadrays
