// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "../cell_types.h"
#include "box_bernstein.h"
#include "clipped_box.h"
#include "engine.h"
#include "source.h"

namespace cutcells::quadrays
{

/// VTK cell types of leaf meshes.
inline constexpr std::uint8_t vtk_tetra = 10;
inline constexpr std::uint8_t vtk_hexahedron = 12;
inline constexpr std::uint8_t vtk_lagrange_quadrilateral = 70;
inline constexpr std::uint8_t vtk_lagrange_hexahedron = 72;

/// Cells for visualisation in CSR layout: the leaves of the engine's
/// decomposition as Lagrange cells in VTK's node order, and whole cells as
/// linear VTK cells.
template <std::floating_point T>
struct LeafMesh
{
    std::vector<T> points;                   ///< physical coordinates, 3 per node
    std::vector<std::int32_t> connectivity;  ///< node indices of all cells
    std::vector<std::int32_t> offsets = {0}; ///< size n_cells + 1, offsets[0] = 0
    std::vector<std::uint8_t> vtk_types;     ///< VTK cell type per cell
    std::vector<std::int32_t> parent;        ///< background cell per cell
    std::vector<std::int32_t> degree;        ///< polynomial degree per cell (1 for linear cells)

    int n_points() const { return static_cast<int>(points.size()) / 3; }
    int n_cells() const { return static_cast<int>(offsets.size()) - 1; }
};

/// @brief Leaf cells of one part of one cell, appended to @p mesh.
///
/// Every piece the engine integrates becomes a Lagrange cell of the given
/// degree whose nodes are the images of an equispaced grid under the piece's
/// parametrisation: hexahedra for volume parts, oriented to a positive
/// Jacobian, and quadrilaterals for the interface, with their normal along
/// grad phi. A picture of the leaves shows exactly what is integrated.
///
/// @param phi          the cell's level set (bernstein_source, analytic_source)
/// @param parent_cell  background cell of the leaves
/// @param stats        counts leaves dropped because their nodes did not line up
template <std::floating_point T>
void append_leaves(const ClippedBox<T>& cell, const Source<T>& phi, Part part, int degree, const Options& opt,
                   std::int32_t parent_cell, LeafMesh<T>& mesh, Stats& stats);

/// @brief append_leaves for a Bernstein form on the cell's box (cell_bernstein_on_box).
template <std::floating_point T>
void append_leaves(const ClippedBox<T>& cell, const BoxBernstein<T>& phi, Part part, int degree,
                   const Options& opt, std::int32_t parent_cell, LeafMesh<T>& mesh, Stats& stats);

/// @brief A whole cell as a linear VTK cell (tetrahedron or hexahedron),
/// appended to @p mesh.
///
/// @param vertex_coords  vertices in Basix order, flat, 3 per vertex
template <std::floating_point T>
void append_linear_cell(cell::type cell_type, std::span<const T> vertex_coords,
                        std::int32_t parent_cell, LeafMesh<T>& mesh);

/// @brief Write @p mesh to a .vtu file, with the cell data "parent_id".
template <std::floating_point T>
void write_leaves(const std::string& filename, const LeafMesh<T>& mesh);

} // namespace cutcells::quadrays
