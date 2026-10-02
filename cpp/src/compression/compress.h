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

#include "../quadrature.h"

/// Compression of quadrature rules: a positive rule replaced by a subset of
/// its points, with new positive weights, that integrates a polynomial space
/// exactly as the original does (Caratheodory-Tchakaloff). A space of
/// dimension M needs at most M points, whatever the original rule had.
namespace cutcells::compression
{

/// Polynomial space whose moments a compressed rule keeps, in the reference
/// coordinates of the rule's points.
enum class MomentSpace : std::uint8_t
{
    tensor = 0, ///< degree <= p in each coordinate (Q_p): quadrilaterals, hexahedra
    total = 1   ///< total degree <= p (P_p): triangles, tetrahedra
};

/// @brief "tensor" or "total".
std::string moment_space_to_str(MomentSpace space);

/// @brief The MomentSpace named "tensor" or "total".
/// @throws std::invalid_argument for any other name
MomentSpace string_to_moment_space(const std::string& name);

/// @brief Dimension of the space: (p + 1)^dim for tensor, binom(p + dim, dim) for total.
int n_moments(MomentSpace space, int dim, int degree);

/// Counters of compress_rules, accumulated over calls.
struct CompressionStats
{
    int n_rules = 0;                 ///< rules read
    int n_compressed = 0;            ///< rules reduced
    int n_skipped = 0;               ///< rules with a negative weight, copied unchanged
    std::int64_t points_before = 0;  ///< points read
    std::int64_t points_after = 0;   ///< points written
    double max_residual = 0;         ///< largest moment error relative to the rule's total |weight|
};

/// @brief A positive rule with the moments of @p weights at @p points over
/// the space, on at most n_moments(space, dim, degree) of the points.
///
/// Points are recombined in groups (Tchernychova and Lyons) until a direct
/// Caratheodory reduction is cheap; the basis is Legendre on the points'
/// bounding box, and directions the points do not resolve (a piece flat in
/// one coordinate, say) are dropped, so such pieces keep fewer points.
///
/// @param points       dim coordinates per point
/// @param weights      one weight per point
/// @param out_points   selected points (overwritten)
/// @param out_weights  their weights (overwritten)
/// @param residual     largest moment error relative to sum |weights|
/// @return false, with the rule copied, if a weight is negative or the rule
///         has no more points than the space has moments
template <std::floating_point T>
bool compress_rule(std::span<const T> points, std::span<const T> weights, int dim, int degree,
                   MomentSpace space, std::vector<T>& out_points, std::vector<T>& out_weights,
                   double& residual);

/// @brief compress_rule for every rule of @p rules, into @p out (same parent
/// cells, same order); rules run in parallel with OpenMP when it is enabled.
///
/// For Q_k elements on affine hexahedra, stiffness and mass integrands lie in
/// Q_{2k}: degree = 2k keeps their assembly exact relative to the original
/// rule. On affine tetrahedra with P_k elements use space total, degree 2k.
template <std::floating_point T>
void compress_rules(const quadrature::QuadratureRules<T>& rules, int degree, MomentSpace space,
                    quadrature::QuadratureRules<T>& out, CompressionStats& stats);

} // namespace cutcells::compression
