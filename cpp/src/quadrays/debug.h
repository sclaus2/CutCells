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

#include "clipped_box.h"
#include "engine.h"
#include "source.h"

namespace cutcells::quadrays
{

/// Who a cell is, for write_cell_report: what the engine does not know of it.
template <std::floating_point T>
struct CellLabel
{
    std::int64_t cell = -1;              ///< its index in the mesh
    std::string cell_type;
    int gdim = 0;
    std::vector<T> vertices;             ///< gdim per vertex
    std::vector<std::string> level_sets; ///< the names of the phis
    std::string selection;               ///< the part, as text
    int dim = 0;                         ///< the part's dimension
    int order = 0;                       ///< the order asked for
    int q = 0;                           ///< points per segment
    double seconds = 0;                  ///< how long the cell took
};

/// @brief Writes what the engine did in a cell to a JSON file, for studying
/// cells that take many points (Options::debug_dir;
/// `python -m cutcells.quadrays_debug` plots it): the cell, its level sets on a
/// uniform grid of the unit box, the boxes of the recursion and the points.
///
/// @param path    the file
/// @param label   who the cell is
/// @param box     the cell as a clipped box
/// @param phis    its level sets on the box
/// @param opt     the options of the run
/// @param trace   the boxes of the recursion
/// @param points  the cell's points, every n-th if there are more than 300000
/// @param stats   the cell's counters
/// @throws std::runtime_error if the file cannot be written
template <std::floating_point T>
void write_cell_report(const std::string& path, const CellLabel<T>& label, const ClippedBox<T>& box,
                       std::span<const Source<T>> phis, const Options& opt, const BoxTrace<T>& trace,
                       const CellPoints<T>& points, const Stats& stats);

} // namespace cutcells::quadrays
