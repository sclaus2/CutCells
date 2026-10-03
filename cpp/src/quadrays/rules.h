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
/// with u0 + u1 + u2 <= 1, a triangle the unit square with u0 + u1 <= 1. The
/// engine reduces the dimension one height direction at a time, accepts a
/// direction only where bounds of the level sets certify it, and bisects the
/// box otherwise (engine.h). This header is the module's entry point:
/// quadrature rules for one cell and a selection.
namespace cutcells::quadrays
{

/// @brief Quadrature rule of the part of one cell that @p terms select,
/// appended to @p rules.
///
/// Points are in the parent cell's reference coordinates (tdim of them),
/// weights are physical, as in quadrature::QuadratureRules. The rule is
/// appended only if it has points.
///
/// @param cell         the cell as a clipped box (make_clipped_box)
/// @param phis         the cell's level sets (bernstein_source, analytic_source)
/// @param terms        the selection; bit l of a term refers to phis[l] (integrate)
/// @param q            Gauss-Legendre points per segment of each height line
/// @param opt          engine options
/// @param parent_cell  background cell of the rule
/// @param rules        output, appended to
/// @param stats        counters, accumulated
template <std::floating_point T>
void append_rules(const ClippedBox<T>& cell, std::span<const Source<T>> phis, std::span<const SelectionTerm> terms,
                  int q, const Options& opt, std::int32_t parent_cell, quadrature::QuadratureRules<T>& rules,
                  Stats& stats);

/// @brief append_rules for one level set and one part of it (part_of).
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
/// @param cell_type      triangle, quadrilateral, tetrahedron, hexahedron, prism
///                       or pyramid
/// @param vertex_coords  vertices in Basix order, flat, tdim per vertex
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
/// vertices (tdim coordinates each), for an analytic level set in physical
/// coordinates (u2 = 0 on 2D cells), appended to @p rules. Any cell type of
/// make_clipped_box.
template <std::floating_point T>
void append_cell_rules(cell::type cell_type, std::span<const T> vertex_coords, const AnalyticLevelSet& phi,
                       const SelectionTerm& term, int level_set, int q, const Options& opt,
                       std::int32_t parent_cell, quadrature::QuadratureRules<T>& rules, Stats& stats);

} // namespace cutcells::quadrays
