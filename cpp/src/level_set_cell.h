// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <span>
#include <vector>

#include "adapt_cell.h"
#include "cell_types.h"
#include "level_set.h"
#include "mesh_view.h"

namespace cutcells
{

template <std::floating_point T, std::integral I = int>
struct LevelSetCell
{
    const LevelSetFunction<T, I>* global_level_set = nullptr;  ///< non-owning view of the global level set function

    int level_set_id = 0 ; ///< identifier for the level set function on this cell (e.g., for multi-level-set scenarios)

    // --- Cell geometry ---
    cell::type cell_type = cell::type::point;  ///< cell type of the background cell
    int tdim = 0;                              ///< topological dimension
    I cell_id = static_cast<I>(-1);            ///< background cell index

    // --- Evaluation interface (uniform for both backends) ---
    // Evaluate the level set at a point xi in reference coordinates.
    T value(std::span<const T> xi) const;

    // Evaluate the gradient at xi in reference coordinates.
    void grad(std::span<const T> xi, std::span<T> g) const;

    // --- Polynomial backend storage ---
    // Bernstein coefficients on this cell.
    // For a polynomial level set, the global nodal values are extracted for
    // this cell, then converted from the nodal (Lagrange) basis to the
    // Bernstein basis on the reference cell.
    std::vector<T> bernstein_coeffs;
    int bernstein_order = 0;

    // The original nodal values on this cell (before conversion), kept for
    // reference and debugging.
    std::vector<T> nodal_values;
    int nodal_order = 0;
};

/// Create a LevelSetCell and AdaptCell for a single background cell.
///
/// @param global_ls      The global level set function.
/// @param mesh           The background mesh.
/// @param cell_id        Index of the background cell to process.
///
/// @return A pair of (LevelSetCell, AdaptCell) initialized for the given cell.
template <std::floating_point T, std::integral I = int>
LevelSetCell<T, I>
make_cell_level_set(const LevelSetFunction<T, I>& global_ls,
                    I cell_id);

} // namespace cutcells