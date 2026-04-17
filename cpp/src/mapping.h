// Copyright (c) 2024 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include "cut_cell.h"
#include "cell_types.h"

#include <span>
#include <vector>
#include <concepts>
#include <array>
#include <stdexcept>

namespace cutcells::cell
{

/// @brief Indices into parent_vertex_coords for the Jacobian columns.
///
/// Column k of J  =  parent_vertex[ col[k] ] - parent_vertex[0].
/// Up to tdim columns; unused entries are -1.
///
/// Choices follow VTK first-order vertex ordering so that J is diagonal
/// (and thus the identity) for the unit reference cell.
///
///  interval:       v0=0   v1=1                     → col[0]=1
///  triangle:       v0=(0,0)  v1=(1,0)  v2=(0,1)    → col={1,2}
///  tetrahedron:    v0=(0,0,0) v1=(1,0,0) v2=(0,1,0) v3=(0,0,1) → col={1,2,3}
///  quadrilateral:  v0=(0,0) v1=(1,0) v2=(1,1) v3=(0,1)         → col={1,3}
///  hexahedron:     v0=(0,0,0) v1=(1,0,0) v3=(0,1,0) v4=(0,0,1) → col={1,3,4}
///  prism:          v0=(0,0,0) v1=(1,0,0) v2=(0,1,0) v3=(0,0,1) → col={1,2,3}
///  pyramid:        v0=(0,0,0) v1=(1,0,0) v3=(0,1,0) v4=apex    → col={1,3,4}
inline std::array<int, 3> jacobian_col_indices(type cell_type)
{
    switch (cell_type)
    {
      case type::interval:      return {1, -1, -1};
      case type::triangle:      return {1,  2, -1};
      case type::tetrahedron:   return {1,  2,  3};
      case type::quadrilateral: return {1,  3, -1};
      case type::hexahedron:    return {1,  3,  4};
      case type::prism:         return {1,  2,  3};
      case type::pyramid:       return {1,  3,  4};
      default:
        throw std::invalid_argument("mapping: unsupported parent cell type");
    }
}

/// @brief Compute the affine volume scaling factor |det J| (or √det(JᵀJ) for
///        embedded cells) for a single cell.
///
/// For volume cells  (tdim == gdim)  returns |det J|, the absolute value of
/// the Jacobian determinant of the affine map from the reference element to
/// physical space.
///
/// For surface/embedded cells  (tdim < gdim)  returns the Gramian root
/// √(det(JᵀJ)) which gives the correct area / length scaling factor.
///
/// @param cell_type     VTK cell type
/// @param phys_verts    flat physical vertex coordinates in VTK ordering
///                      (nv * gdim values)
/// @param gdim          geometric dimension of the embedding space (1, 2, or 3)
/// @returns volume scaling factor ≥ 0
template <std::floating_point T>
T affine_volume_factor(type cell_type, const T* phys_verts, int gdim);

/// @brief Push all cut vertices from parent reference space into physical space.
///
/// Uses the affine map  x_phys = x0 + J * X_ref  where J is built from the
/// first-order parent cell vertex differences (VTK ordering) and x0 is the first
/// parent physical vertex.
///
/// Precondition  : _vertex_coords holds reference coordinates in the parent
///                 reference cell; _parent_vertex_coords and _parent_cell_type
///                 are set.
/// Postcondition : _vertex_coords_phys is allocated and filled.
///
/// Note: supported for gdim == topological dimension of parent (volume cells).
template <std::floating_point T>
void compute_physical_cut_vertices(CutCell<T>& cut_cell);

/// @brief Complete a cut cell that was computed in physical space.
///
/// The cutting routine wrote physical coordinates into _vertex_coords.
/// This function:
///   1. moves _vertex_coords  →  _vertex_coords_phys
///   2. pulls _vertex_coords_phys back through  J * X = x - x0  to overwrite
///      _vertex_coords with the reference coordinates.
///
/// After this call _vertex_coords is always in reference space and
/// _vertex_coords_phys is always in physical space.
///
/// Precondition  : _vertex_coords holds raw physical coordinates from the cutter;
///                 _parent_vertex_coords and _parent_cell_type are set.
/// Postcondition : _vertex_coords   = reference coordinates (parent ref. space)
///                 _vertex_coords_phys = physical coordinates
///
/// Note: supported for gdim == topological dimension of parent (volume cells).
template <std::floating_point T>
void complete_from_physical(CutCell<T>& cut_cell);

/// @brief Affine push-forward for a batch of n reference points.
///
/// Maps X_ref (flat, n*gdim) from parent reference space to physical space
/// x_phys (flat, n*gdim) via  x = x0 + J * X  where J is derived from
/// parent_vertex_coords and parent_type in VTK vertex ordering.
///
/// @param parent_type         cell type of the parent element
/// @param parent_vertex_coords flat physical vertex coords in VTK ordering
/// @param gdim                geometric / topological dimension (must equal
///                            the topological dimension of parent_type)
/// @param X_ref               flat input: n * gdim reference coordinates
/// @param x_phys              flat output: n * gdim physical coordinates
template <std::floating_point T>
void push_forward_affine(type parent_type,
                         const std::vector<T>& parent_vertex_coords,
                         int gdim,
                         std::span<const T> X_ref,
                         std::span<T> x_phys);

/// @brief Affine push-forward for parent-reference points when tdim <= gdim.
///
/// Maps points in the parent reference space (dimension = get_tdim(parent_type))
/// to physical space (dimension = gdim). This supports embedded cells, e.g.
/// triangles in 3D or intervals in 2D/3D.
///
/// @param parent_type          cell type of the parent element
/// @param parent_vertex_coords flat physical vertex coords in VTK ordering
/// @param gdim                 physical embedding dimension
/// @param X_ref                flat input: n * tdim reference coordinates
/// @returns                    flat output: n * gdim physical coordinates
template <std::floating_point T>
std::vector<T> push_forward_affine_map(type parent_type,
                                       const std::vector<T>& parent_vertex_coords,
                                       int gdim,
                                       std::span<const T> X_ref);

/// @brief Affine pullback for a batch of n physical points.
///
/// Maps x_phys (flat, n*gdim) from physical space to parent reference space
/// X_ref (flat, n*gdim) by solving  J * X = x - x0.
///
/// @param parent_type         cell type of the parent element
/// @param parent_vertex_coords flat physical vertex coords in VTK ordering
/// @param gdim                geometric / topological dimension
/// @param x_phys              flat input: n * gdim physical coordinates
/// @param X_ref               flat output: n * gdim reference coordinates
template <std::floating_point T>
void pull_back_affine(type parent_type,
                      const std::vector<T>& parent_vertex_coords,
                      int gdim,
                      std::span<const T> x_phys,
                      std::span<T> X_ref);

} // namespace cutcells::cell
