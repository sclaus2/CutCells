// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <functional>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "cell_types.h"
#include "mesh_view.h"

namespace cutcells
{

enum class LevelSetType
{
  Analytical,
  Polynomial
};

template <std::floating_point T, std::integral I = int>
struct LevelSetMeshData
{
  int gdim = 0;
  int tdim = 0;
  int degree = 1;

  // Flattened coordinates of the discrete level-set dofs:
  // [x0_0, x0_1, ..., x1_0, x1_1, ...]
  std::vector<T> dof_coordinates;

  // CSR-like cell -> global dof map. The local ordering on each cell follows
  // Basix ordering for the corresponding cell type and degree.
  std::vector<I> cell_dofs;
  std::vector<I> cell_offsets;

  // Optional flattened reference coordinates for one local cell dof row:
  // [xi0_0, xi0_1, ..., xi1_0, xi1_1, ...]. When empty, CutCells falls
  // back to its equispaced reference lattice for the requested cell type.
  std::vector<T> cell_reference_points;

  // Optional non-owning cell -> global dof map. When cell_offsets_view is empty,
  // cell_dofs_view is interpreted as fixed-width rows with cell_dof_stride
  // entries and cell_dofs_per_cell active entries.
  std::span<const I> cell_dofs_view;
  std::span<const I> cell_offsets_view;
  I cell_count = 0;
  I cell_dofs_per_cell = 0;
  I cell_dof_stride = 0;
  std::span<const StridedRowBlock<I>> cell_dof_blocks;

  // Optional callbacks for dof layouts that cannot expose contiguous rows.
  std::function<I(I)> cell_dof_row_size;
  std::function<I(I, I)> cell_dof_row_value;
  std::function<void(I, std::vector<I>&)> cell_dof_row_gather;

  // Optional keep-alive anchor for external dofmap memory.
  std::shared_ptr<const void> owner;

  // cell types, one per cell when available.
  std::vector<cell::type> cell_types;
  std::span<const CellTypeBlock<I>> cell_type_blocks;
  cell::type uniform_cell_type = cell::type::point;
  bool has_uniform_cell_type = false;

  // Provenance aligned with AdaptCell:
  // 0 = parent vertex, 1 = parent edge, 2 = parent face, 3 = parent cell interior.
  std::vector<int8_t> dof_parent_dim;
  std::vector<int32_t> dof_parent_id;
  std::vector<T> dof_parent_param;
  std::vector<int32_t> dof_parent_param_offset;

  RowAccess<I> cell_dof_row_access() const
  {
    RowAccess<I> rows;
    rows.values = !cell_dofs.empty()
                      ? std::span<const I>(cell_dofs.data(), cell_dofs.size())
                      : cell_dofs_view;
    rows.offsets = !cell_offsets.empty()
                       ? std::span<const I>(cell_offsets.data(),
                                            cell_offsets.size())
                       : cell_offsets_view;
    rows.row_count = cell_count;
    rows.row_width = cell_dofs_per_cell;
    rows.row_stride = cell_dof_stride;
    rows.blocks = cell_dof_blocks;
    rows.row_size_fn = cell_dof_row_size;
    rows.row_value_fn = cell_dof_row_value;
    rows.row_gather_fn = cell_dof_row_gather;
    return rows;
  }

  I num_dofs() const
  {
    if (gdim <= 0)
      return 0;
    return static_cast<I>(dof_coordinates.size() / static_cast<std::size_t>(gdim));
  }

  I num_cells() const
  {
    return cell_dof_row_access().num_rows();
  }

  bool has_cell_types() const
  {
    return !cell_types.empty() || !cell_type_blocks.empty()
        || has_uniform_cell_type;
  }

  cell::type cell_type(I cell_id) const
  {
    if (!cell_types.empty())
      return cell_types[static_cast<std::size_t>(cell_id)];
    for (const auto& block : cell_type_blocks)
    {
      if (block.contains(cell_id))
        return block.type;
    }
    if (has_uniform_cell_type)
      return uniform_cell_type;
    throw std::runtime_error("LevelSetMeshData::cell_type not available");
  }

  I cell_num_dofs(I cell_id) const
  {
    return cell_dof_row_access().row_size(cell_id);
  }

  std::span<const I> cell_dofs_span(I cell_id) const
  {
    return cell_dof_row_access().direct_row(cell_id);
  }

  std::span<const I> cell_dofs_span(I cell_id, std::vector<I>& scratch) const
  {
    return cell_dof_row_access().row(cell_id, scratch);
  }

  std::span<const I> cell_dofs_storage_span() const
  {
    if (!cell_dofs.empty())
      return std::span<const I>(cell_dofs.data(), cell_dofs.size());
    return cell_dofs_view;
  }

  std::span<const I> cell_offsets_span() const
  {
    if (!cell_offsets.empty())
      return std::span<const I>(cell_offsets.data(), cell_offsets.size());
    return cell_offsets_view;
  }

  const T* dof_coordinate(I dof_id) const
  {
    return dof_coordinates.data()
         + static_cast<std::size_t>(dof_id) * static_cast<std::size_t>(gdim);
  }

  std::span<const T> dof_parent_param_span(I dof_id) const
  {
    const std::size_t begin = static_cast<std::size_t>(
        dof_parent_param_offset[static_cast<std::size_t>(dof_id)]);
    const std::size_t end = static_cast<std::size_t>(
        dof_parent_param_offset[static_cast<std::size_t>(dof_id) + 1]);
    return std::span<const T>(dof_parent_param.data() + begin, end - begin);
  }
};

template <std::floating_point T, std::integral I = int>
struct LevelSetFunction
{

  std::string name = "phi";
  LevelSetType type = LevelSetType::Analytical;
  int gdim = 0;
  
  // cell_id == -1 means "unknown / not provided"
  // The value and gradient of the level set function in physical coordinates (x)
  // together with background cell id.
  std::function<T(const T*, I)> value_fn;
  std::function<void(const T*, I, T*)> grad_fn;

  // Legacy low-order nodal storage on mesh vertices.
  std::span<const T> nodal_values;

  // Optional higher-order mesh/dof layout and corresponding global dof values.
  LevelSetMeshData<T, I> mesh_data;
  bool has_mesh_data_storage = false;
  std::span<const T> dof_values;

  // Optional keep-alive anchor for Python / external memory.
  std::shared_ptr<const void> owner;

  bool has_value() const
  {
    return static_cast<bool>(value_fn);
  }

  bool has_gradient() const
  {
    return static_cast<bool>(grad_fn);
  }

  bool has_nodal_values() const
  {
    return !nodal_values.empty();
  }

  bool has_mesh_data() const
  {
    return has_mesh_data_storage;
  }

  bool has_dof_values() const
  {
    return !dof_values.empty();
  }

  T value(const T* x, I cell_id = static_cast<I>(-1)) const
  {
    if (!value_fn)
      throw std::runtime_error("LevelSetFunction::value not available");
    return value_fn(x, cell_id);
  }

  void grad(const T* x, I cell_id, T* g) const
  {
    if (!grad_fn)
      throw std::runtime_error("LevelSetFunction::grad not available");
    grad_fn(x, cell_id, g);
  }

  T value_at_node(I node_id) const
  {
    if (nodal_values.empty())
      throw std::runtime_error("LevelSetFunction::value_at_node not available");
    return nodal_values[static_cast<std::size_t>(node_id)];
  }
};

template <std::floating_point T, std::integral I = int>
LevelSetMeshData<T, I> create_level_set_mesh_data(
    const MeshView<T, I>& mesh, int degree, T merge_tol = T(-1));

template <std::floating_point T, std::integral I = int>
LevelSetMeshData<T, I> create_level_set_mesh_data(
    int gdim,
    int tdim,
    int degree,
    std::span<const T> dof_coordinates,
    std::span<const I> cell_dofs,
    std::span<const I> cell_offsets,
    std::span<const cell::type> cell_types = {},
    std::span<const T> cell_reference_points = {});

template <std::floating_point T, std::integral I = int>
LevelSetMeshData<T, I> create_level_set_mesh_data_view(
    int gdim,
    int tdim,
    int degree,
    std::vector<T>&& dof_coordinates,
    std::span<const I> cell_dofs,
    I num_cells,
    I dofs_per_cell,
    I dof_stride,
    std::span<const cell::type> cell_types = {});

template <std::floating_point T, std::integral I = int>
LevelSetMeshData<T, I> create_level_set_mesh_data_view(
    int gdim,
    int tdim,
    int degree,
    std::vector<T>&& dof_coordinates,
    std::span<const StridedRowBlock<I>> cell_dof_blocks,
    I num_cells,
    std::span<const cell::type> cell_types = {});

template <std::floating_point T, std::integral I = int>
LevelSetMeshData<T, I> create_level_set_mesh_data_view(
    int gdim,
    int tdim,
    int degree,
    std::vector<T>&& dof_coordinates,
    std::span<const StridedRowBlock<I>> cell_dof_blocks,
    I num_cells,
    std::span<const CellTypeBlock<I>> cell_type_blocks);

template <std::floating_point T, std::integral I = int>
LevelSetFunction<T, I> create_level_set_function(
    LevelSetMeshData<T, I> mesh_data,
    std::span<const T> dof_values,
    std::string name = "phi");

template <std::floating_point T, std::integral I = int>
LevelSetFunction<T, I> create_level_set_function_view(
    LevelSetMeshData<T, I> mesh_data,
    std::span<const T> dof_values,
    std::string name = "phi",
    std::shared_ptr<const void> owner = nullptr);

} // namespace cutcells
