// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <span>

#include "../quadrature.h"
#include "analytic.h"
#include "box_bernstein.h"
#include "clipped_box.h"
#include "engine.h"
#include "source.h"

/// quadrays: height-function quadrature on cut cells.
///
/// Every cell is a box clipped by half-spaces; a tetrahedron is the unit box
/// with u0 + u1 + u2 <= 1. The engine reduces the dimension one height
/// direction at a time, accepts a direction only where bounds of the level set
/// certify it, and bisects the box otherwise (engine.h). This header is the
/// module's entry point: quadrature rules for one cell and one selection term.
namespace cutcells::quadrays
{

/// @brief Quadrature rule of one part of one cell, appended to @p rules.
///
/// Points are in the parent cell's reference coordinates, weights are
/// physical, as in quadrature::QuadratureRules. The rule is appended only if it
/// has points.
///
/// @param cell         the cell as a clipped box (make_clipped_box)
/// @param phi          the cell's level set (bernstein_source, analytic_source)
/// @param part         the selected part (part_of)
/// @param q            Gauss-Legendre points per segment of each height line
/// @param opt          engine options
/// @param parent_cell  background cell of the rule
/// @param rules        output, appended to
/// @param stats        counters, accumulated
template <std::floating_point T>
void append_rules(const ClippedBox<T>& cell, const Source<T>& phi, Part part, int q, const Options& opt,
                  std::int32_t parent_cell, quadrature::QuadratureRules<T>& rules, Stats& stats);

/// @brief append_rules for a Bernstein form on the cell's box (cell_bernstein_on_box).
template <std::floating_point T>
void append_rules(const ClippedBox<T>& cell, const BoxBernstein<T>& phi, Part part, int q,
                  const Options& opt, std::int32_t parent_cell,
                  quadrature::QuadratureRules<T>& rules, Stats& stats);

/// @brief Quadrature rule of one part of one cell given by its type, vertices
/// and the Bernstein coefficients of its level set, appended to @p rules.
///
/// @param cell_type      tetrahedron or hexahedron
/// @param vertex_coords  vertices in Basix order, flat, 3 per vertex
/// @param degree         degree of the level set
/// @param coeffs         Bernstein coefficients in CutCells' order (bernstein.h)
/// @param term           compiled selection term on one level set
/// @param level_set      index of that level set in the term's bitmasks
template <std::floating_point T>
void append_cell_rules(cell::type cell_type, std::span<const T> vertex_coords, int degree,
                       std::span<const T> coeffs, const SelectionTerm& term, int level_set, int q,
                       const Options& opt, std::int32_t parent_cell,
                       quadrature::QuadratureRules<T>& rules, Stats& stats);

/// @brief Quadrature rule of one part of one cell given by its type and
/// vertices, for an analytic level set in physical coordinates, appended to
/// @p rules.
template <std::floating_point T>
void append_cell_rules(cell::type cell_type, std::span<const T> vertex_coords, const AnalyticLevelSet& phi,
                       const SelectionTerm& term, int level_set, int q, const Options& opt,
                       std::int32_t parent_cell, quadrature::QuadratureRules<T>& rules, Stats& stats);

} // namespace cutcells::quadrays
