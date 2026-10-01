// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <vector>

namespace cutcells::proto
{

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<Vec3, 3>; ///< row-major, m[i][k]

/// Half-space c . u <= d in the unit-box coordinates u of a ClippedBox.
struct HalfSpace
{
    Vec3 c = {0, 0, 0};
    double d = 0;
};

/// A cell, or a piece of one, as the unit box clipped by half-spaces.
///
///   x(u)  = origin + jacobian u          maps the unit box to physical space,
///   xi(u) = ref_origin + ref_jacobian u  maps it to the parent cell's reference
///                                        coordinates, where points are reported.
struct ClippedBox
{
    Vec3 origin = {0, 0, 0};
    Mat3 jacobian = {};
    Vec3 ref_origin = {0, 0, 0};
    Mat3 ref_jacobian = {};
    std::vector<HalfSpace> clips;
};

/// Clipped region of a ClippedBox (unit box intersected with its clips) in box
/// coordinates.
struct Polytope
{
    std::vector<Vec3> vertices;
    std::vector<std::array<int, 2>> edges;
    std::vector<Vec3> face_normals;           ///< outward unit normals
    std::vector<std::vector<int>> face_edges; ///< edges lying on each face

    int n_vertices() const { return static_cast<int>(vertices.size()); }
};

/// Grid hexahedron [lo, lo + h]^3; reference coordinates equal box coordinates.
ClippedBox hex_cell(const Vec3& lo, double h);

/// Tetrahedron X[0..3] in the reference box of its affine map: box corner at X[0],
/// box edges towards X[1], X[2], X[3], clipped by u0 + u1 + u2 <= 1. Box
/// coordinates equal the Basix reference coordinates of the tetrahedron.
ClippedBox tet_cell(const std::array<Vec3, 4>& X);

/// Determinant of the box-to-physical map.
double jacobian_determinant(const ClippedBox& box);

/// Inverse of the box-to-physical Jacobian.
Mat3 inverse_jacobian(const ClippedBox& box);

Vec3 physical_point(const ClippedBox& box, const Vec3& u);
Vec3 reference_point(const ClippedBox& box, const Vec3& u);

/// Child covering [lo, hi] of the box coordinates of @p box; clips carry over.
ClippedBox sub_box(const ClippedBox& box, const Vec3& lo, const Vec3& hi);

/// Box axis with the longest physical edge.
int longest_axis(const ClippedBox& box);

/// True if [lo, hi] (box coordinates) may meet the clipped region. Tests each
/// half-space separately: conservative, it may keep a sub-box that misses it.
bool may_meet_clips(const ClippedBox& box, const Vec3& lo, const Vec3& hi, double tol = 1e-12);

/// True if u satisfies all clips.
bool inside_clips(const ClippedBox& box, const Vec3& u);

/// Exact clipped region by vertex enumeration; empty if nothing is left.
Polytope clipped_polytope(const ClippedBox& box);

} // namespace cutcells::proto
