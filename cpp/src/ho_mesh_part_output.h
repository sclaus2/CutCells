// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "cut_mesh.h"
#include "ho_cut_mesh.h"
#include "quadrature.h"
#include "quadrays/engine.h"
#include "quadrays/leaves.h"

namespace cutcells::output
{

enum class QuadratureBackend
{
    Straight,
    AlgoimBernstein,
    AlgoimGeneral,
    Quadrays
};

struct SelectedZeroEntityInfo
{
    std::int32_t cut_cell_id = -1;
    std::int32_t parent_cell_id = -1;
    std::int32_t local_zero_entity_id = -1;
    std::int32_t dimension = -1;
};

template <std::floating_point T, std::integral I = int>
std::vector<SelectedZeroEntityInfo> selected_zero_entity_infos(
    const HOMeshPart<T, I>& part);

template <std::floating_point T, std::integral I = int>
mesh::CutMesh<T> visualization_mesh(const HOMeshPart<T, I>& part,
                                    bool include_uncut_cells);

template <std::floating_point T, std::integral I = int>
quadrature::QuadratureRules<T> quadrature_rules(const HOMeshPart<T, I>& part,
                                                int order,
                                                bool include_uncut_cells);

QuadratureBackend quadrature_backend_from_string(std::string_view backend);

template <std::floating_point T, std::integral I = int>
quadrature::QuadratureRules<T> quadrature_rules(const HOMeshPart<T, I>& part,
                                                int order,
                                                bool include_uncut_cells,
                                                QuadratureBackend backend);

/// @brief Quadrature rules of a part from the quadrays engine.
///
/// Cut cells get the engine's rules, with @p order Gauss-Legendre points per
/// segment of each height line (as `order` counts points for backend
/// "algoim"); uncut cells of volume parts get the straight rules of degree
/// @p order. The selection must be one term on one level set, whose cells
/// carry Bernstein coefficients; parent cells are tetrahedra or hexahedra.
template <std::floating_point T, std::integral I = int>
quadrature::QuadratureRules<T> quadrays_quadrature_rules(
    const HOMeshPart<T, I>& part,
    int order,
    bool include_uncut_cells,
    const quadrays::Options& options = {});

/// @brief Leaf cells of a part: every piece the quadrays engine integrates in
/// a cut cell, as a Lagrange cell of the given degree, and, for volume parts
/// with @p include_uncut_cells, the uncut cells as linear cells.
template <std::floating_point T, std::integral I = int>
quadrays::LeafMesh<T> quadrays_leaves(const HOMeshPart<T, I>& part,
                                      int degree,
                                      bool include_uncut_cells,
                                      const quadrays::Options& options = {});

template <std::floating_point T, std::integral I = int>
std::vector<std::pair<std::string, quadrature::QuadratureRules<T>>>
paired_quadrature_rules(
    const std::vector<std::pair<std::string, HOMeshPart<T, I>>>& parts,
    int order,
    bool include_uncut_cells,
    QuadratureBackend backend);

/// Compute the selected volume fraction for each parent cell represented by a
/// volume mesh part.
///
/// The mesh part already encodes the selection expression. For each returned
/// parent id, the matching value is the reference-space measure of the selected
/// part divided by the reference-space measure of the full parent cell. Multiple
/// selected pieces with the same parent are accumulated into one fraction.
///
/// @throws std::runtime_error if @p part is not a volume part.
template <std::floating_point T, std::integral I = int>
std::pair<std::vector<I>, std::vector<T>>
volume_fractions(const HOMeshPart<T, I>& part);

} // namespace cutcells::output
