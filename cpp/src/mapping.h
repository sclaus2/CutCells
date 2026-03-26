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
#include <algorithm>
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
/// Maps X_ref (flat, n*tdim) from parent reference space to physical space
/// x_phys (flat, n*gdim) via the parent-cell P1 mapping.
///
/// @param parent_type         cell type of the parent element
/// @param parent_vertex_coords flat physical vertex coords in VTK ordering
/// @param gdim                geometric dimension of the embedding space
/// @param X_ref               flat input: n * tdim reference coordinates
/// @param x_phys              flat output: n * gdim physical coordinates
template <std::floating_point T>
void push_forward_affine(type parent_type,
                         const std::vector<T>& parent_vertex_coords,
                         int gdim,
                         std::span<const T> X_ref,
                         std::span<T> x_phys);

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

/// @brief Push-forward one reference point to physical space using P1 shape
/// functions on the parent reference cell.
///
/// Works for simplex and non-simplex cells (quadrilateral/hexahedron/prism/
/// pyramid) using the corresponding first-order reference basis.
template <std::floating_point T>
inline void ref_to_phys_affine(type parent_type,
                               std::span<const T> parent_vertex_coords,
                               int gdim,
                               int tdim,
                               const T* x_ref,
                               T* x_phys);

/// @brief Compute Jacobian J = d x_phys / d x_ref at one reference point.
///
/// J is stored row-major with shape (gdim, tdim):
///   J[d * tdim + k] = d x_d / d xi_k.
template <std::floating_point T>
inline void compute_jacobian(type parent_type,
                             std::span<const T> parent_vertex_coords,
                             int gdim,
                             int tdim,
                             const T* x_ref,
                             T* J);

namespace detail
{
template <std::floating_point T>
inline void eval_p1_shape(type ct, std::span<const T> x_ref, std::span<T> N)
{
    const int n_vertices = get_num_vertices(ct);
    if (static_cast<int>(N.size()) != n_vertices)
        throw std::invalid_argument("eval_p1_shape: invalid N size");

    std::fill(N.begin(), N.end(), T(0));

    switch (ct)
    {
    case type::interval:
    {
        const T x = x_ref[0];
        N[0] = T(1) - x;
        N[1] = x;
        return;
    }
    case type::triangle:
    {
        const T x = x_ref[0];
        const T y = x_ref[1];
        N[0] = T(1) - x - y;
        N[1] = x;
        N[2] = y;
        return;
    }
    case type::tetrahedron:
    {
        const T x = x_ref[0];
        const T y = x_ref[1];
        const T z = x_ref[2];
        N[0] = T(1) - x - y - z;
        N[1] = x;
        N[2] = y;
        N[3] = z;
        return;
    }
    case type::quadrilateral:
    {
        const T x = x_ref[0];
        const T y = x_ref[1];
        N[0] = (T(1) - x) * (T(1) - y);
        N[1] = x * (T(1) - y);
        N[2] = x * y;
        N[3] = (T(1) - x) * y;
        return;
    }
    case type::hexahedron:
    {
        const T x = x_ref[0];
        const T y = x_ref[1];
        const T z = x_ref[2];
        const T mx = T(1) - x;
        const T my = T(1) - y;
        const T mz = T(1) - z;
        N[0] = mx * my * mz;
        N[1] = x * my * mz;
        N[2] = x * y * mz;
        N[3] = mx * y * mz;
        N[4] = mx * my * z;
        N[5] = x * my * z;
        N[6] = x * y * z;
        N[7] = mx * y * z;
        return;
    }
    case type::prism:
    {
        const T x = x_ref[0];
        const T y = x_ref[1];
        const T z = x_ref[2];
        const T a = T(1) - x - y;
        N[0] = a * (T(1) - z);
        N[1] = x * (T(1) - z);
        N[2] = y * (T(1) - z);
        N[3] = a * z;
        N[4] = x * z;
        N[5] = y * z;
        return;
    }
    case type::pyramid:
    {
        const T x = x_ref[0];
        const T y = x_ref[1];
        const T z = x_ref[2];
        const T r = T(1) - z;
        if (r <= T(1e-12))
        {
            N[0] = T(0);
            N[1] = T(0);
            N[2] = T(0);
            N[3] = T(0);
            N[4] = T(1);
            return;
        }
        const T inv_r = T(1) / r;
        const T xy_r = x * y * inv_r;
        N[0] = r - x - y + xy_r;
        N[1] = x - xy_r;
        N[2] = xy_r;
        N[3] = y - xy_r;
        N[4] = z;
        return;
    }
    default:
        throw std::invalid_argument("eval_p1_shape: unsupported parent cell type");
    }
}

template <std::floating_point T>
inline void eval_p1_shape_grad(type ct, std::span<const T> x_ref, std::span<T> dN)
{
    const int tdim = get_tdim(ct);
    const int n_vertices = get_num_vertices(ct);
    if (static_cast<int>(dN.size()) != n_vertices * tdim)
        throw std::invalid_argument("eval_p1_shape_grad: invalid dN size");

    std::fill(dN.begin(), dN.end(), T(0));

    auto set_grad = [&](int i, int k, T val)
    {
        dN[static_cast<std::size_t>(i * tdim + k)] = val;
    };

    switch (ct)
    {
    case type::interval:
        set_grad(0, 0, T(-1));
        set_grad(1, 0, T(1));
        return;
    case type::triangle:
        set_grad(0, 0, T(-1)); set_grad(0, 1, T(-1));
        set_grad(1, 0, T(1));  set_grad(1, 1, T(0));
        set_grad(2, 0, T(0));  set_grad(2, 1, T(1));
        return;
    case type::tetrahedron:
        set_grad(0, 0, T(-1)); set_grad(0, 1, T(-1)); set_grad(0, 2, T(-1));
        set_grad(1, 0, T(1));  set_grad(1, 1, T(0));  set_grad(1, 2, T(0));
        set_grad(2, 0, T(0));  set_grad(2, 1, T(1));  set_grad(2, 2, T(0));
        set_grad(3, 0, T(0));  set_grad(3, 1, T(0));  set_grad(3, 2, T(1));
        return;
    case type::quadrilateral:
    {
        const T x = x_ref[0];
        const T y = x_ref[1];
        set_grad(0, 0, -(T(1) - y)); set_grad(0, 1, -(T(1) - x));
        set_grad(1, 0,  (T(1) - y)); set_grad(1, 1, -x);
        set_grad(2, 0,  y);          set_grad(2, 1, x);
        set_grad(3, 0, -y);          set_grad(3, 1, T(1) - x);
        return;
    }
    case type::hexahedron:
    {
        const T x = x_ref[0];
        const T y = x_ref[1];
        const T z = x_ref[2];
        const T mx = T(1) - x;
        const T my = T(1) - y;
        const T mz = T(1) - z;

        set_grad(0, 0, -my * mz); set_grad(0, 1, -mx * mz); set_grad(0, 2, -mx * my);
        set_grad(1, 0,  my * mz); set_grad(1, 1, -x  * mz); set_grad(1, 2, -x  * my);
        set_grad(2, 0,  y  * mz); set_grad(2, 1,  x  * mz); set_grad(2, 2, -x  * y);
        set_grad(3, 0, -y  * mz); set_grad(3, 1,  mx * mz); set_grad(3, 2, -mx * y);
        set_grad(4, 0, -my * z ); set_grad(4, 1, -mx * z ); set_grad(4, 2,  mx * my);
        set_grad(5, 0,  my * z ); set_grad(5, 1, -x  * z ); set_grad(5, 2,  x  * my);
        set_grad(6, 0,  y  * z ); set_grad(6, 1,  x  * z ); set_grad(6, 2,  x  * y);
        set_grad(7, 0, -y  * z ); set_grad(7, 1,  mx * z ); set_grad(7, 2,  mx * y);
        return;
    }
    case type::prism:
    {
        const T x = x_ref[0];
        const T y = x_ref[1];
        const T z = x_ref[2];
        set_grad(0, 0, -(T(1) - z)); set_grad(0, 1, -(T(1) - z)); set_grad(0, 2, -(T(1) - x - y));
        set_grad(1, 0,  (T(1) - z)); set_grad(1, 1, T(0));        set_grad(1, 2, -x);
        set_grad(2, 0,  T(0));       set_grad(2, 1, (T(1) - z));  set_grad(2, 2, -y);
        set_grad(3, 0, -z);          set_grad(3, 1, -z);          set_grad(3, 2, (T(1) - x - y));
        set_grad(4, 0,  z);          set_grad(4, 1, T(0));        set_grad(4, 2, x);
        set_grad(5, 0,  T(0));       set_grad(5, 1, z);           set_grad(5, 2, y);
        return;
    }
    case type::pyramid:
    {
        const T x = x_ref[0];
        const T y = x_ref[1];
        const T z = x_ref[2];
        const T r = T(1) - z;
        if (r <= T(1e-12))
        {
            const auto cols = jacobian_col_indices(ct);
            const int ncols = get_tdim(ct);
            for (int k = 0; k < ncols; ++k)
            {
                for (int i = 0; i < n_vertices; ++i)
                    set_grad(i, k, T(0));
                set_grad(0, k, T(-1));
                const int vk = cols[static_cast<std::size_t>(k)];
                if (vk >= 0 && vk < n_vertices)
                    set_grad(vk, k, T(1));
            }
            return;
        }
        const T inv_r = T(1) / r;
        const T inv_r2 = inv_r * inv_r;
        set_grad(0, 0, -T(1) + y * inv_r);
        set_grad(0, 1, -T(1) + x * inv_r);
        set_grad(0, 2, -T(1) + x * y * inv_r2);

        set_grad(1, 0, T(1) - y * inv_r);
        set_grad(1, 1, -x * inv_r);
        set_grad(1, 2, -x * y * inv_r2);

        set_grad(2, 0, y * inv_r);
        set_grad(2, 1, x * inv_r);
        set_grad(2, 2, x * y * inv_r2);

        set_grad(3, 0, -y * inv_r);
        set_grad(3, 1, T(1) - x * inv_r);
        set_grad(3, 2, -x * y * inv_r2);

        set_grad(4, 0, T(0));
        set_grad(4, 1, T(0));
        set_grad(4, 2, T(1));
        return;
    }
    default:
        throw std::invalid_argument("eval_p1_shape_grad: unsupported parent cell type");
    }
}
} // namespace detail

template <std::floating_point T>
inline void ref_to_phys_affine(type parent_type,
                               std::span<const T> parent_vertex_coords,
                               int gdim,
                               int tdim,
                               const T* x_ref,
                               T* x_phys)
{
    if (gdim <= 0 || gdim > 3)
        throw std::invalid_argument("ref_to_phys_affine: gdim must be in [1, 3]");
    if (tdim != get_tdim(parent_type))
        throw std::invalid_argument("ref_to_phys_affine: tdim does not match parent type");

    const int n_vertices = get_num_vertices(parent_type);
    if (static_cast<int>(parent_vertex_coords.size()) != n_vertices * gdim)
        throw std::invalid_argument("ref_to_phys_affine: invalid parent_vertex_coords size");

    std::vector<T> N(static_cast<std::size_t>(n_vertices), T(0));
    detail::eval_p1_shape<T>(
        parent_type,
        std::span<const T>(x_ref, static_cast<std::size_t>(tdim)),
        std::span<T>(N.data(), N.size()));

    for (int d = 0; d < gdim; ++d)
    {
        T val = T(0);
        for (int i = 0; i < n_vertices; ++i)
        {
            val += N[static_cast<std::size_t>(i)]
                * parent_vertex_coords[static_cast<std::size_t>(i * gdim + d)];
        }
        x_phys[d] = val;
    }
}

template <std::floating_point T>
inline void compute_jacobian(type parent_type,
                             std::span<const T> parent_vertex_coords,
                             int gdim,
                             int tdim,
                             const T* x_ref,
                             T* J)
{
    if (gdim <= 0 || gdim > 3)
        throw std::invalid_argument("compute_jacobian: gdim must be in [1, 3]");
    if (tdim != get_tdim(parent_type))
        throw std::invalid_argument("compute_jacobian: tdim does not match parent type");

    const int n_vertices = get_num_vertices(parent_type);
    if (static_cast<int>(parent_vertex_coords.size()) != n_vertices * gdim)
        throw std::invalid_argument("compute_jacobian: invalid parent_vertex_coords size");

    std::vector<T> dN(static_cast<std::size_t>(n_vertices * tdim), T(0));
    detail::eval_p1_shape_grad<T>(
        parent_type,
        std::span<const T>(x_ref, static_cast<std::size_t>(tdim)),
        std::span<T>(dN.data(), dN.size()));

    for (int d = 0; d < gdim; ++d)
    {
        for (int k = 0; k < tdim; ++k)
        {
            T val = T(0);
            for (int i = 0; i < n_vertices; ++i)
            {
                val += parent_vertex_coords[static_cast<std::size_t>(i * gdim + d)]
                     * dN[static_cast<std::size_t>(i * tdim + k)];
            }
            J[static_cast<std::size_t>(d * tdim + k)] = val;
        }
    }
}

} // namespace cutcells::cell
