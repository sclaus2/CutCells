// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include "local_mesh.h"
#include "quadrature.h"

#include <concepts>
#include <vector>

namespace cutcells::quadrature
{

/// @brief Append curved Pk quadrature for one interface entity.
///
/// Integrates along the Pk-curved interface segment (in 2D) using the
/// curved node cache stored in mesh.curved_zero_ref_nodes.
/// Quadrature points are stored in parent background-cell reference space
/// (dimension tdim). Physical weights account for the Pk surface Jacobian.
///
/// If the curved cache is not built (curved_geometry_order < 2), falls back
/// to straight-segment quadrature using only the two endpoint vertices.
///
/// @param mesh            decomposed LocalMesh with curved interface cache filled
/// @param zero_entity_id  index of the zero entity to integrate over
/// @param order           polynomial order of the quadrature rule
/// @param rules           output accumulator (appended to, not reset)
template <std::floating_point T>
void append_interface_quadrature_curved(
    const LocalMesh<T>& mesh,
    int zero_entity_id,
    int order,
    QuadratureRules<T>& rules);

template <std::floating_point T>
void append_interface_quadrature_curved(
    const LocalMesh<T>& mesh,
    int zero_entity_id,
    int order,
    QuadratureRules<T>& rules,
    std::vector<T>& physical_points);

/// @brief Append curved volume quadrature for one inside sub-cell.
///
/// If the sub-cell has an interface face (curved_geometry_order >= 2), uses
/// a blended Pk mapping: F(xi_1, xi_n) = (1-xi_n)*gamma(xi_1) + xi_n*v_opp
/// where gamma is the Pk curve from the cache and v_opp is the opposite vertex.
/// Integration domain is the unit square [0,1]^2 (collapsing at xi_n=1).
///
/// Accepts both consistent orientations (positive or negative det(J) in 2D),
/// and uses |det(J)| in the weight. If det(J) changes sign across quadrature
/// points (fold-over), falls back to the affine path and increments
/// fallback_count.
///
/// If no curved cache exists or the cell does not touch the interface,
/// uses the standard affine path.
///
/// Currently implemented for 2D triangles only.
/// 3D cells and non-simplex types fall back to affine automatically.
///
/// @param mesh          LocalMesh with curved interface cache
/// @param cell_id       inside sub-cell index
/// @param level_set_id  which level set's interface to use
/// @param order         polynomial order of the quadrature rule
/// @param rules         output accumulator (appended to, not reset)
/// @param fallback_count output counter incremented when affine fallback is used
template <std::floating_point T>
void append_volume_quadrature_curved(
    const LocalMesh<T>& mesh,
    int cell_id,
    int level_set_id,
    int order,
    QuadratureRules<T>& rules,
    int& fallback_count);

template <std::floating_point T>
void append_volume_quadrature_curved(
    const LocalMesh<T>& mesh,
    int cell_id,
    int level_set_id,
    int order,
    QuadratureRules<T>& rules,
    int& fallback_count,
    std::vector<T>& physical_points);

/// @brief Build curved quadrature rules for all inside cells and interface
///        entities of one level set in a LocalMesh.
///
/// Iterates over:
///   - all codim-1 zero entities owned by level_set_id
///     → calls append_interface_quadrature_curved for each
///   - all inside sub-cells
///     → calls append_volume_quadrature_curved for each
///
/// Both output rules are reset before filling.
///
/// @param mesh              decomposed LocalMesh with curved cache
/// @param level_set_id      which level set to process
/// @param order             quadrature polynomial order
/// @param volume_rules      output: volume quadrature (appended from scratch)
/// @param interface_rules   output: interface quadrature (appended from scratch)
template <std::floating_point T>
void make_quadrature_curved(
    const LocalMesh<T>& mesh,
    int level_set_id,
    int order,
    QuadratureRules<T>& volume_rules,
    QuadratureRules<T>& interface_rules);

template <std::floating_point T>
void make_quadrature_curved(
    const LocalMesh<T>& mesh,
    int level_set_id,
    int order,
    QuadratureRules<T>& volume_rules,
    std::vector<T>& volume_physical_points,
    QuadratureRules<T>& interface_rules,
    std::vector<T>& interface_physical_points);

} // namespace cutcells::quadrature
