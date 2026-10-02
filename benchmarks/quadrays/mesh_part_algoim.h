// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>

#include <cutcells/part/mesh_part.h>
#include <cutcells/quadrature.h>

/// algoim's quadrature on the parts of the front end, for comparisons. It was
/// the library's backend='algoim' on HOMeshPart until phase 5; the library no
/// longer includes algoim or links LAPACK.
namespace cutcells::benchmarks
{

/// @brief Rules on a part that one clause on one level set selects: on its cut
/// quadrilaterals and hexahedra algoim's multi-polynomial engine on the
/// Bernstein form (general = false) or its general engine with the Bernstein
/// form as a functor (general = true); on cut intervals the roots, weight 1
/// each; on whole cells and faces in the zero set the lookup tables' rules.
///
/// Points are in the cells' reference coordinates, weights are physical; the
/// cut cells' rules come first, cells ascending.
///
/// @param order  Gauss points per direction (algoim's q), 1 to 10
template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> algoim_rules(const part::MeshPart<T, I>& part, int order, bool include_uncut_cells,
                                            bool general);

} // namespace cutcells::benchmarks
