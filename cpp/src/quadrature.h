// Copyright (c) 2024 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include "cut_cell.h"
#include "cell_types.h"

#include <vector>
#include <cstdint>
#include <concepts>
#include <span>
#include <string>

namespace cutcells::quadrature
{

/// Flat batch of quadrature rules for a collection of cut cells.
///
/// Conventions:
///   - _points are in the **parent reference element** (dimension _tdim each).
///   - _weights are **physical integration weights** (|det J_phys| scaled).
///   - Rule i covers parent cell _parent_map[i].
///   - Points for rule i live at _points[_offset[i]*_tdim .. _offset[i+1]*_tdim).
///   - Weights for rule i live at _weights[_offset[i] .. _offset[i+1]).
///
/// Preconditions before calling make_quadrature / append_quadrature:
///   - cut_cell._vertex_coords     must hold parent **reference** coordinates.
///   - cut_cell._vertex_coords_phys must hold **physical** coordinates.
///   - Both are filled after calling compute_physical_cut_vertices() or
///     complete_from_physical() from mapping.h.
template <std::floating_point T>
struct QuadratureRules
{
    /// Topological dimension of the parent reference element (= dim of each point)
    int _tdim = 0;

    /// Flat quadrature points in parent reference space.
    /// Size: total_num_points * _tdim
    std::vector<T> _points;

    /// Physical integration weights.
    /// Size: total_num_points
    std::vector<T> _weights;

    /// CSR offsets: rule i has points in [_offset[i], _offset[i+1]).
    /// Size: num_rules + 1   (_offset[0] = 0 always)
    std::vector<int32_t> _offset;

    /// Parent cell index for each rule.
    /// Size: num_rules
    std::vector<int32_t> _parent_map;
};

/// Append the quadrature contribution of one cut cell.
///
/// Iterates over the subcells of cut_cell, fetches the canonical rule for
/// each subcell type at the given polynomial order, maps canonical points to
/// the parent reference element (using _vertex_coords), computes physical
/// weights (using _vertex_coords_phys), and appends everything to rules.
///
/// @param cut_cell  enriched cut cell (both coord frames must be set)
/// @param order     polynomial order for the canonical quadrature rule
/// @param rules     output accumulator (appended to, not reset)
template <std::floating_point T>
void append_quadrature(const cutcells::cell::CutCell<T>& cut_cell,
                       int order,
                       QuadratureRules<T>& rules);

/// Build QuadratureRules for a vector of cut cells.
/// rules._offset and rules._parent_map are built alongside the point/weight data.
/// @param cut_cells  enriched cut-cell vector
/// @param order      polynomial order
/// @param rules      output (reset before use)
template <std::floating_point T>
void make_quadrature(const std::vector<cutcells::cell::CutCell<T>>& cut_cells,
                     int order,
                     QuadratureRules<T>& rules);

/// Convenience overload that returns a freshly constructed QuadratureRules.
template <std::floating_point T>
QuadratureRules<T> make_quadrature(
    const std::vector<cutcells::cell::CutCell<T>>& cut_cells, int order);

/// @brief Generate flat quadrature rules for all mesh cells matching the
///        requested level-set domain.
///
/// Iterates over all cells in the VTK-format mesh.  For each cell the level-set
/// domain is classified:
///   - *Full cells*  (entirely inside / outside): a standard reference-space
///     rule is used, scaled by |det J_phys|.  When @p triangulate is true and
///     the cell is not a simplex, the cell is first decomposed into simplices.
///   - *Cut cells*   (intersected): the cell is cut, the reference/physical
///     coordinates are completed via mapping.h, and append_quadrature is called.
///
/// All input arrays are flat and follow VTK conventions:
///   @p points        npts × 3  (gdim = 3 always for VTK meshes)
///   @p connectivity  concatenated vertex-index lists for all cells
///   @p offset        CSR offsets: cell i uses connectivity[offset[i]..offset[i+1])
///   @p vtk_type      VTK element-type code per cell
///
/// @param ls_vals       level-set values at all mesh vertices (size = npts)
/// @param points        flat mesh vertex coordinates (size = npts * 3)
/// @param connectivity  flat vertex connectivity list
/// @param offset        CSR offsets (size = ncells + 1)
/// @param vtk_type      VTK cell types (size = ncells)
/// @param cut_type_str  one of "phi<0", "phi>0", "phi=0"
/// @param triangulate   if true, non-simplex cells are split into simplices
/// @param order         polynomial quadrature order
/// @returns             flat QuadratureRules in parent reference space
template <std::floating_point T>
QuadratureRules<T> runtime_quadrature(
    std::span<const T>   ls_vals,
    std::span<const T>   points,
    std::span<const int> connectivity,
    std::span<const int> offset,
    std::span<const int> vtk_type,
    const std::string&   cut_type_str,
    bool triangulate,
    int  order);

/// @brief Push reference-space quadrature points to physical space.
///
/// For each rule i in @p rules, the reference-space points stored in
/// rules._points are mapped to physical space using the affine pushforward
/// of parent cell @c rules._parent_map[i].
///
/// @param rules        quadrature rules (output of runtime_quadrature or
///                     make_quadrature)
/// @param points        flat mesh vertex coordinates (size = npts * 3)
/// @param connectivity  flat vertex connectivity list
/// @param offset        CSR offsets (size = ncells + 1)
/// @param vtk_type      VTK cell types (size = ncells)
/// @returns             flat physical point coordinates
///                      (total_num_points * 3)
template <std::floating_point T>
std::vector<T> physical_points(
    const QuadratureRules<T>&  rules,
    std::span<const T>   points,
    std::span<const int> connectivity,
    std::span<const int> offset,
    std::span<const int> vtk_type);

} // namespace cutcells::quadrature
