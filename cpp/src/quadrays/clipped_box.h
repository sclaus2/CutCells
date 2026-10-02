// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <concepts>
#include <cstdint>
#include <span>
#include <vector>

#include "../cell_types.h"

namespace cutcells::quadrays
{

template <std::floating_point T>
using Vec3 = std::array<T, 3>;

/// Row-major 3 x 3 matrix, m[i][k].
template <std::floating_point T>
using Mat3 = std::array<std::array<T, 3>, 3>;

/// Half-space c . u <= d in the box coordinates u of a ClippedBox.
template <std::floating_point T>
struct HalfSpace
{
    Vec3<T> c = {0, 0, 0};
    T d = 0;
};

/// A cell, or a piece of one, as the unit box clipped by half-spaces.
///
///   x(u)  = origin + jacobian u           maps the unit box to physical space,
///   xi(u) = ref_origin + ref_jacobian u   maps it to the parent cell's reference
///                                         coordinates, where points are reported.
///
/// A hexahedron is its own box; a tetrahedron is the box of its affine map
/// clipped by u0 + u1 + u2 <= 1. Maps are affine, as everywhere in CutCells.
template <std::floating_point T>
struct ClippedBox
{
    Vec3<T> origin = {0, 0, 0};
    Mat3<T> jacobian = {};
    Vec3<T> ref_origin = {0, 0, 0};
    Mat3<T> ref_jacobian = {};
    std::vector<HalfSpace<T>> clips;
};

/// Clipped region of a ClippedBox (the unit box intersected with its clips) in
/// box coordinates. Faces list their edges in CSR layout.
template <std::floating_point T>
struct Polytope
{
    std::vector<Vec3<T>> vertices;
    std::vector<std::array<int, 2>> edges;
    std::vector<Vec3<T>> face_normals; ///< outward unit normals
    std::vector<int> face_edges;        ///< edges of all faces
    std::vector<int> face_offsets = {0}; ///< size n_faces + 1

    int n_vertices() const { return static_cast<int>(vertices.size()); }
    int n_faces() const { return static_cast<int>(face_offsets.size()) - 1; }
};

/// @brief The ClippedBox of a cell from its type and vertices.
///
/// @param cell_type      tetrahedron or hexahedron
/// @param vertex_coords  vertices in Basix order, flat, gdim per vertex
/// @param gdim           geometric dimension (3)
/// @param box            output; box coordinates equal the cell's Basix
///                       reference coordinates
template <std::floating_point T>
void make_clipped_box(cell::type cell_type, std::span<const T> vertex_coords, int gdim,
                      ClippedBox<T>& box);

/// @brief Determinant of the box-to-physical map.
template <std::floating_point T>
T jacobian_determinant(const ClippedBox<T>& box);

/// @brief Inverse of the box-to-physical Jacobian.
template <std::floating_point T>
Mat3<T> inverse_jacobian(const ClippedBox<T>& box);

/// @brief Physical point of the box coordinates @p u.
template <std::floating_point T>
Vec3<T> physical_point(const ClippedBox<T>& box, const Vec3<T>& u);

/// @brief Parent reference coordinates of the box coordinates @p u.
template <std::floating_point T>
Vec3<T> reference_point(const ClippedBox<T>& box, const Vec3<T>& u);

/// @brief Child covering [lo, hi] of the box coordinates of @p box; clips carry over.
template <std::floating_point T>
ClippedBox<T> sub_box(const ClippedBox<T>& box, const Vec3<T>& lo, const Vec3<T>& hi);

/// @brief Box axis with the longest physical edge.
template <std::floating_point T>
int longest_axis(const ClippedBox<T>& box);

/// @brief True if [lo, hi] (box coordinates) may meet the clipped region.
///
/// Tests each half-space separately: conservative, it may keep a sub-box that
/// misses the region.
template <std::floating_point T>
bool may_meet_clips(const ClippedBox<T>& box, const Vec3<T>& lo, const Vec3<T>& hi,
                    T tol = static_cast<T>(1e-12));

/// @brief True if @p u satisfies all clips up to @p tol.
template <std::floating_point T>
bool inside_clips(const ClippedBox<T>& box, const Vec3<T>& u, T tol = T(0));

/// @brief Exact clipped region by vertex enumeration; empty if nothing is left.
template <std::floating_point T>
void clipped_polytope(const ClippedBox<T>& box, Polytope<T>& poly);

} // namespace cutcells::quadrays
