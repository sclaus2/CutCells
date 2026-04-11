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

  // cell types, one per cell when available.
  std::vector<cell::type> cell_types;

  // Provenance aligned with AdaptCell:
  // 0 = parent vertex, 1 = parent edge, 2 = parent face, 3 = parent cell interior.
  std::vector<int8_t> dof_parent_dim;
  std::vector<int32_t> dof_parent_id;
  std::vector<T> dof_parent_param;
  std::vector<int32_t> dof_parent_param_offset;

  I num_dofs() const
  {
    if (gdim <= 0)
      return 0;
    return static_cast<I>(dof_coordinates.size() / static_cast<std::size_t>(gdim));
  }

  I num_cells() const
  {
    if (cell_offsets.empty())
      return 0;
    return static_cast<I>(cell_offsets.size() - 1);
  }

  I cell_num_dofs(I cell_id) const
  {
    return cell_offsets[static_cast<std::size_t>(cell_id) + 1]
         - cell_offsets[static_cast<std::size_t>(cell_id)];
  }

  std::span<const I> cell_dofs_span(I cell_id) const
  {
    const std::size_t begin = static_cast<std::size_t>(
        cell_offsets[static_cast<std::size_t>(cell_id)]);
    const std::size_t end = static_cast<std::size_t>(
        cell_offsets[static_cast<std::size_t>(cell_id) + 1]);
    return std::span<const I>(cell_dofs.data() + begin, end - begin);
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
  std::shared_ptr<LevelSetMeshData<T, I>> mesh_data;
  std::span<const T> dof_values;

  // Optional keep-alive anchor for Python / external memory.
  std::shared_ptr<void> owner;

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
    return static_cast<bool>(mesh_data);
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
    std::span<const cell::type> cell_types = {});

template <std::floating_point T, std::integral I = int>
LevelSetFunction<T, I> create_level_set_function(
    std::shared_ptr<LevelSetMeshData<T, I>> mesh_data,
    std::span<const T> dof_values,
    std::string name = "phi");

} // namespace cutcells
