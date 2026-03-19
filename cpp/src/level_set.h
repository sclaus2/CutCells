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

namespace cutcells
{

template <std::floating_point T, std::integral I = int>
struct LevelSetFunction
{
  // cell_id == -1 means "unknown / not provided"
  std::function<T(const T*, I)> value_fn;
  std::function<void(const T*, I, T*)> grad_fn;

  // Optional nodal values
  std::span<const T> nodal_values;

  // Optional keep-alive anchor for Python / external memory
  std::shared_ptr<void> owner;

  int gdim = 0;
  int degree = -1;

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
};

} // namespace cutcells
