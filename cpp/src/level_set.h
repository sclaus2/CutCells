// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <functional>
#include <memory>
#include <span>
#include <stdexcept>

namespace cutcells
{

template <std::floating_point T, std::integral I>
struct MeshView;

template <std::floating_point T, std::integral I>
struct LocalLevelSetFunction;

template <std::floating_point T, std::integral I = int>
struct LevelSetFunction
{
  enum class Kind : uint8_t
  {
    callable = 0,
    fem_nodal = 1
  };

  // cell_id == -1 means "unknown / not provided"
  std::function<T(const T*, I)> value_fn;
  std::function<void(const T*, I, T*)> grad_fn;

  // Optional nodal values
  std::span<const T> nodal_values;

  // Optional keep-alive anchor for Python / external memory
  std::shared_ptr<void> owner;

  // Optional global mesh for FEM-backed local extraction
  const MeshView<T, I>* mesh = nullptr;

  int gdim = 0;
  int degree = -1;
  Kind kind = Kind::callable;

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

  LocalLevelSetFunction<T, I> local(I cell_id) const;
};

template <std::floating_point T, std::integral I = int>
inline LevelSetFunction<T, I> make_level_set_function_from_callable(
    std::function<T(const T*, I)>            value_fn,
    std::function<void(const T*, I, T*)>     grad_fn = {},
    int                                      gdim = 0,
    int                                      degree = -1,
    std::shared_ptr<void>                    owner = nullptr)
{
  LevelSetFunction<T, I> level_set;
  level_set.value_fn = std::move(value_fn);
  level_set.grad_fn = std::move(grad_fn);
  level_set.owner = std::move(owner);
  level_set.gdim = gdim;
  level_set.degree = degree;
  level_set.kind = LevelSetFunction<T, I>::Kind::callable;
  return level_set;
}

template <std::floating_point T, std::integral I = int>
LevelSetFunction<T, I> make_level_set_function_from_fem(
    const MeshView<T, I>&                      mesh,
    std::span<const T>                         nodal_values,
    int                                        degree,
    std::shared_ptr<void>                      owner = nullptr);

} // namespace cutcells
