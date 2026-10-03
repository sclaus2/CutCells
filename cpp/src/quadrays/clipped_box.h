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
/// A hexahedron is its own box. In the reference frame (BoxFrame), a
/// tetrahedron is the box of its affine map clipped by u0 + u1 + u2 <= 1, a
/// prism the box clipped by u0 + u1 <= 1, a pyramid the box clipped by
/// u0 + u2 <= 1 and u1 + u2 <= 1; in the orthogonal frame, a simplex is its
/// bounding box in an orthonormal frame, clipped by its facets. Maps are
/// affine, as everywhere in CutCells.
///
/// Cells of dimension 2 (triangles, quadrilaterals) live in the plane u2 = 0:
/// their box is [0, 1]^2 x {0}, the third column of the maps is e2 (so their
/// determinants are the cells' areas), and points keep u2 = 0.
template <std::floating_point T>
struct ClippedBox
{
    int tdim = 3; ///< 2 or 3
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

/// How make_clipped_box chooses the box of a cell.
enum class BoxFrame : std::uint8_t
{
    /// The unit box of the cell's affine map from vertex 0 along its edges:
    /// box coordinates are the cell's reference coordinates.
    reference = 0,
    /// Triangles and tetrahedra: over the physical axes and the orthonormal
    /// frames two of the cell's edges span, the one in which the cell's
    /// bounding box is smallest, the physical axes on ties; the cell is that
    /// box clipped by its facets. Prisms: their bottom triangle's smallest
    /// rectangle in an orthonormal frame of its plane, times the lateral edge,
    /// clipped by the triangle's sides. Height directions along orthogonal
    /// axes: one of them always suits a smooth zero set, while the sheared
    /// reference frame of a Kuhn tetrahedron may offer none. Pyramids keep
    /// the reference frame: (1 - z)^n phi, their engine's function, varies
    /// wildly in the box outside the pyramid, of which an orthogonal box has
    /// more. Hexahedra, quadrilaterals and degenerate cells: reference.
    orthogonal = 1
};

/// @brief The ClippedBox of a cell from its type and vertices.
///
/// The maps are affine. Quadrilaterals, hexahedra, prisms and pyramids are
/// exact only if they are affine images of their reference cells
/// (parallelograms, parallelepipeds, ...).
///
/// @param cell_type      triangle or quadrilateral (gdim 2); tetrahedron,
///                       hexahedron, prism or pyramid (gdim 3)
/// @param vertex_coords  vertices in Basix order, flat, gdim per vertex
/// @param gdim           geometric dimension, the cell's dimension
/// @param box            output
/// @param frame          the frame of the box (BoxFrame)
template <std::floating_point T>
void make_clipped_box(cell::type cell_type, std::span<const T> vertex_coords, int gdim,
                      ClippedBox<T>& box, BoxFrame frame = BoxFrame::orthogonal);

/// @brief True if quadrays takes cells of this type (make_clipped_box).
inline bool supported_cell(cell::type cell_type)
{
    switch (cell_type)
    {
    case cell::type::triangle:
    case cell::type::quadrilateral:
    case cell::type::tetrahedron:
    case cell::type::hexahedron:
    case cell::type::prism:
    case cell::type::pyramid:
        return true;
    default:
        return false;
    }
}

/// @brief Determinant of the box-to-physical map.
template <std::floating_point T>
T jacobian_determinant(const ClippedBox<T>& box);

/// @brief Determinant of the map from the parent's reference coordinates to
/// physical space: reference measures times its absolute value are physical.
template <std::floating_point T>
T reference_jacobian_determinant(const ClippedBox<T>& box);

/// @brief Inverse of the box-to-physical Jacobian.
template <std::floating_point T>
Mat3<T> inverse_jacobian(const ClippedBox<T>& box);

/// @brief Physical point of the box coordinates @p u.
template <std::floating_point T>
Vec3<T> physical_point(const ClippedBox<T>& box, const Vec3<T>& u);

/// @brief Parent reference coordinates of the box coordinates @p u.
template <std::floating_point T>
Vec3<T> reference_point(const ClippedBox<T>& box, const Vec3<T>& u);

/// @brief Box coordinates of the parent reference coordinates @p xi (the
/// inverse of reference_point).
template <std::floating_point T>
Vec3<T> box_point(const ClippedBox<T>& box, const Vec3<T>& xi);

/// @brief True if the box coordinates of @p box are its parent's reference
/// coordinates (BoxFrame::reference, or a hexahedron or quadrilateral).
template <std::floating_point T>
bool reference_frame(const ClippedBox<T>& box);

/// @brief Child covering [lo, hi] of the box coordinates of @p box; clips carry
/// over. Boxes of 2D cells keep lo[2] = 0 and hi[2] = 1.
template <std::floating_point T>
ClippedBox<T> sub_box(const ClippedBox<T>& box, const Vec3<T>& lo, const Vec3<T>& hi);

/// @brief Box axis (of the first tdim) with the longest physical edge.
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
