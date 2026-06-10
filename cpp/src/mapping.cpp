// Copyright (c) 2024 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include "mapping.h"

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <array>

namespace cutcells::cell
{

namespace
{
  // -------------------------------------------------------------------------
  // Per-cell-type helpers
  // -------------------------------------------------------------------------

  /// Number of first-order vertices per parent cell type.
  inline int num_parent_vertices(type cell_type)
  {
    switch (cell_type)
    {
      case type::interval:      return 2;
      case type::triangle:      return 3;
      case type::tetrahedron:   return 4;
      case type::quadrilateral: return 4;
      case type::hexahedron:    return 8;
      case type::prism:         return 6;
      case type::pyramid:       return 5;
      default:
        throw std::invalid_argument("mapping: unsupported parent cell type");
    }
  }

  // jacobian_col_indices is now a public inline in mapping.h

  // -------------------------------------------------------------------------
  // Small dense matrix-vector multiply  y = J * x
  // J is stored column-major: col k occupies J[k*gdim .. (k+1)*gdim - 1]
  // -------------------------------------------------------------------------

  template <std::floating_point T>
  inline void mat_vec_1x1(const T* J, const T* x, T* y)
  {
    y[0] = J[0] * x[0];
  }

  template <std::floating_point T>
  inline void mat_vec_2x2(const T* J, const T* x, T* y)
  {
    // J = | J[0]  J[2] |
    //     | J[1]  J[3] |
    y[0] = J[0]*x[0] + J[2]*x[1];
    y[1] = J[1]*x[0] + J[3]*x[1];
  }

  template <std::floating_point T>
  inline void mat_vec_3x3(const T* J, const T* x, T* y)
  {
    // J = | J[0]  J[3]  J[6] |
    //     | J[1]  J[4]  J[7] |
    //     | J[2]  J[5]  J[8] |
    y[0] = J[0]*x[0] + J[3]*x[1] + J[6]*x[2];
    y[1] = J[1]*x[0] + J[4]*x[1] + J[7]*x[2];
    y[2] = J[2]*x[0] + J[5]*x[1] + J[8]*x[2];
  }

  // -------------------------------------------------------------------------
  // Solve J * y = x for y  (J column-major, same layout as above)
  // -------------------------------------------------------------------------

  template <std::floating_point T>
  inline void solve_1x1(const T* J, const T* x, T* y)
  {
    y[0] = x[0] / J[0];
  }

  template <std::floating_point T>
  inline void solve_2x2(const T* J, const T* x, T* y)
  {
    // J = | a  b |    det = a*d - b*c
    //     | c  d |
    // column-major: a=J[0] c=J[1] b=J[2] d=J[3]
    const T a = J[0], c = J[1], b = J[2], d = J[3];
    const T inv_det = T(1) / (a*d - b*c);
    y[0] = inv_det * ( d*x[0] - b*x[1]);
    y[1] = inv_det * (-c*x[0] + a*x[1]);
  }

  template <std::floating_point T>
  inline void solve_3x3(const T* J, const T* x, T* y)
  {
    // Column-major layout:
    //   J = | J[0]  J[3]  J[6] |
    //       | J[1]  J[4]  J[7] |
    //       | J[2]  J[5]  J[8] |
    //
    // Cofactors of column 0 (used for det and first row of J^{-1}):
    const T C00 =  J[4]*J[8] - J[7]*J[5];   // cof(0,0)
    const T C10 = -(J[3]*J[8] - J[6]*J[5]); // cof(1,0)
    const T C20 =  J[3]*J[7] - J[6]*J[4];   // cof(2,0)

    const T det     = J[0]*C00 + J[1]*C10 + J[2]*C20;
    const T inv_det = T(1) / det;

    // Remaining cofactors:
    const T C01 = -(J[1]*J[8] - J[7]*J[2]); // cof(0,1)
    const T C11 =  J[0]*J[8] - J[6]*J[2];   // cof(1,1)
    const T C21 = -(J[0]*J[7] - J[6]*J[1]); // cof(2,1)
    const T C02 =  J[1]*J[5] - J[4]*J[2];   // cof(0,2)
    const T C12 = -(J[0]*J[5] - J[3]*J[2]); // cof(1,2)
    const T C22 =  J[0]*J[4] - J[3]*J[1];   // cof(2,2)

    // J^{-1}[row i, col j] = cof(j, i) / det
    // y = J^{-1} * x
    y[0] = inv_det * (C00*x[0] + C10*x[1] + C20*x[2]);
    y[1] = inv_det * (C01*x[0] + C11*x[1] + C21*x[2]);
    y[2] = inv_det * (C02*x[0] + C12*x[1] + C22*x[2]);
  }

  // -------------------------------------------------------------------------
  // Build the column-major Jacobian matrix (gdim x gdim) for the affine map.
  // J_out must point to storage for at least gdim*gdim values.
  // -------------------------------------------------------------------------
  template <std::floating_point T>
  void build_jacobian(type parent_type,
                      const std::vector<T>& pv,
                      int gdim,
                      T* J_out)
  {
    const auto cols = jacobian_col_indices(parent_type);
    const T* v0 = pv.data();
    for (int col = 0; col < gdim; ++col)
    {
      const T* vi = pv.data() + cols[col] * gdim;
      for (int row = 0; row < gdim; ++row)
        J_out[col * gdim + row] = vi[row] - v0[row];
    }
  }

} // anonymous namespace

// =============================================================================
// Public API implementations
// =============================================================================

template <std::floating_point T>
void push_forward_affine(type parent_type,
                         const std::vector<T>& parent_vertex_coords,
                         int gdim,
                         std::span<const T> X_ref,
                         std::span<T> x_phys)
{
  assert(X_ref.size() == x_phys.size());
  const int n = static_cast<int>(X_ref.size()) / gdim;
  const T* x0 = parent_vertex_coords.data();

  T J[9] = {};
  build_jacobian(parent_type, parent_vertex_coords, gdim, J);

  if (gdim == 1)
  {
    for (int i = 0; i < n; ++i)
    {
      T y[1];
      mat_vec_1x1(J, X_ref.data() + i, y);
      x_phys[i] = x0[0] + y[0];
    }
  }
  else if (gdim == 2)
  {
    for (int i = 0; i < n; ++i)
    {
      T y[2];
      mat_vec_2x2(J, X_ref.data() + 2*i, y);
      x_phys[2*i    ] = x0[0] + y[0];
      x_phys[2*i + 1] = x0[1] + y[1];
    }
  }
  else // gdim == 3
  {
    for (int i = 0; i < n; ++i)
    {
      T y[3];
      mat_vec_3x3(J, X_ref.data() + 3*i, y);
      x_phys[3*i    ] = x0[0] + y[0];
      x_phys[3*i + 1] = x0[1] + y[1];
      x_phys[3*i + 2] = x0[2] + y[2];
    }
  }
}

template <std::floating_point T>
std::vector<T> push_forward_affine_map(type parent_type,
                                       const std::vector<T>& parent_vertex_coords,
                                       int gdim,
                                       std::span<const T> X_ref)
{
  const int tdim = get_tdim(parent_type);
  if (tdim <= 0 || tdim > gdim)
    throw std::invalid_argument("push_forward_affine_map: invalid tdim/gdim combination");
  if (X_ref.size() % static_cast<std::size_t>(tdim) != 0)
  {
    throw std::invalid_argument(
        "push_forward_affine_map: X_ref size must be a multiple of the parent tdim");
  }

  const int n = static_cast<int>(X_ref.size()) / tdim;
  std::vector<T> x_phys(static_cast<std::size_t>(n * gdim), T(0));
  const T* x0 = parent_vertex_coords.data();
  const auto cols = jacobian_col_indices(parent_type);

  for (int i = 0; i < n; ++i)
  {
    const T* X = X_ref.data() + i * tdim;
    T* x = x_phys.data() + i * gdim;

    for (int d = 0; d < gdim; ++d)
      x[d] = x0[d];

    for (int col = 0; col < tdim; ++col)
    {
      const T* vi = parent_vertex_coords.data() + cols[col] * gdim;
      for (int d = 0; d < gdim; ++d)
        x[d] += X[col] * (vi[d] - x0[d]);
    }
  }

  return x_phys;
}

template <std::floating_point T>
void pull_back_affine(type parent_type,
                      const std::vector<T>& parent_vertex_coords,
                      int gdim,
                      std::span<const T> x_phys,
                      std::span<T> X_ref)
{
  assert(x_phys.size() == X_ref.size());
  const int n = static_cast<int>(x_phys.size()) / gdim;
  const T* x0 = parent_vertex_coords.data();

  T J[9] = {};
  build_jacobian(parent_type, parent_vertex_coords, gdim, J);

  if (gdim == 1)
  {
    for (int i = 0; i < n; ++i)
    {
      T rhs = x_phys[i] - x0[0];
      solve_1x1(J, &rhs, X_ref.data() + i);
    }
  }
  else if (gdim == 2)
  {
    T rhs[2];
    for (int i = 0; i < n; ++i)
    {
      rhs[0] = x_phys[2*i    ] - x0[0];
      rhs[1] = x_phys[2*i + 1] - x0[1];
      solve_2x2(J, rhs, X_ref.data() + 2*i);
    }
  }
  else // gdim == 3
  {
    T rhs[3];
    for (int i = 0; i < n; ++i)
    {
      rhs[0] = x_phys[3*i    ] - x0[0];
      rhs[1] = x_phys[3*i + 1] - x0[1];
      rhs[2] = x_phys[3*i + 2] - x0[2];
      solve_3x3(J, rhs, X_ref.data() + 3*i);
    }
  }
}

template <std::floating_point T>
void compute_physical_cut_vertices(CutCell<T>& cut_cell)
{
  const int gdim   = cut_cell._gdim;
  const int n_verts = static_cast<int>(cut_cell._vertex_coords.size()) / gdim;
  cut_cell._vertex_coords_phys.resize(n_verts * gdim);

  push_forward_affine(cut_cell._parent_cell_type,
                      cut_cell._parent_vertex_coords,
                      gdim,
                      std::span<const T>(cut_cell._vertex_coords),
                      std::span<T>(cut_cell._vertex_coords_phys));
}

template <std::floating_point T>
void complete_from_physical(CutCell<T>& cut_cell)
{
  const int gdim    = cut_cell._gdim;
  const int n_verts = static_cast<int>(cut_cell._vertex_coords.size()) / gdim;

  // Step 1: move the raw physical output from the cutter into _vertex_coords_phys.
  cut_cell._vertex_coords_phys = std::move(cut_cell._vertex_coords);

  // Step 2: pull back the physical coords to reference space and store in _vertex_coords.
  cut_cell._vertex_coords.resize(n_verts * gdim);
  pull_back_affine(cut_cell._parent_cell_type,
                   cut_cell._parent_vertex_coords,
                   gdim,
                   std::span<const T>(cut_cell._vertex_coords_phys),
                   std::span<T>(cut_cell._vertex_coords));
}

// =============================================================================
// affine_volume_factor
// =============================================================================

template <std::floating_point T>
T affine_volume_factor(type cell_type, const T* phys_verts, int gdim)
{
  const int tdim = get_tdim(cell_type);
  const auto cols = jacobian_col_indices(cell_type);
  const T* v0 = phys_verts;

  if (tdim == gdim)
  {
    // Square Jacobian — build column-major gdim×gdim matrix and take |det|.
    T J[9] = {};
    for (int c = 0; c < tdim; ++c)
    {
      const T* vi = phys_verts + cols[c] * gdim;
      for (int r = 0; r < gdim; ++r)
        J[c * gdim + r] = vi[r] - v0[r];
    }
    if (gdim == 1)
    {
      return std::abs(J[0]);
    }
    else if (gdim == 2)
    {
      // det of 2×2 column-major J: col0=(J[0],J[1]), col1=(J[2],J[3])
      const T det = J[0]*J[3] - J[2]*J[1];
      return std::abs(det);
    }
    else // gdim == 3
    {
      // det of 3×3 column-major J by cofactor expansion along first row
      const T det = J[0]*(J[4]*J[8] - J[7]*J[5])
                  - J[1]*(J[3]*J[8] - J[6]*J[5])
                  + J[2]*(J[3]*J[7] - J[6]*J[4]);
      return std::abs(det);
    }
  }
  else if (tdim == 1)
  {
    // Embedded edge: J is a single column vector of length gdim.
    const T* vi = phys_verts + cols[0] * gdim;
    T norm2 = T(0);
    for (int r = 0; r < gdim; ++r)
    {
      const T d = vi[r] - v0[r];
      norm2 += d * d;
    }
    return std::sqrt(norm2);
  }
  else if (tdim == 2 && gdim == 3)
  {
    // Embedded surface: J has two columns (3-vectors). Volume factor = |j0 × j1|.
    const T* vi0 = phys_verts + cols[0] * gdim;
    const T* vi1 = phys_verts + cols[1] * gdim;
    const T j0[3] = {vi0[0]-v0[0], vi0[1]-v0[1], vi0[2]-v0[2]};
    const T j1[3] = {vi1[0]-v0[0], vi1[1]-v0[1], vi1[2]-v0[2]};
    const T cx = j0[1]*j1[2] - j0[2]*j1[1];
    const T cy = j0[2]*j1[0] - j0[0]*j1[2];
    const T cz = j0[0]*j1[1] - j0[1]*j1[0];
    return std::sqrt(cx*cx + cy*cy + cz*cz);
  }
  else
  {
    throw std::invalid_argument("affine_volume_factor: unsupported tdim/gdim combination");
  }
}

// =============================================================================
// Explicit instantiations for double and float
// =============================================================================
template void push_forward_affine(type, const std::vector<double>&, int,
                                  std::span<const double>, std::span<double>);
template void push_forward_affine(type, const std::vector<float>&,  int,
                                  std::span<const float>,  std::span<float>);
template std::vector<double> push_forward_affine_map(
    type, const std::vector<double>&, int, std::span<const double>);
template std::vector<float> push_forward_affine_map(
    type, const std::vector<float>&, int, std::span<const float>);

template void pull_back_affine(type, const std::vector<double>&, int,
                               std::span<const double>, std::span<double>);
template void pull_back_affine(type, const std::vector<float>&,  int,
                               std::span<const float>,  std::span<float>);

template void compute_physical_cut_vertices(CutCell<double>&);
template void compute_physical_cut_vertices(CutCell<float>&);

template void complete_from_physical(CutCell<double>&);
template void complete_from_physical(CutCell<float>&);

template double affine_volume_factor(type, const double*, int);
template float  affine_volume_factor(type, const float*, int);

} // namespace cutcells::cell
