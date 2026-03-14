// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <cstddef>
#include <memory>
#include <span>
#include <stdexcept>

namespace cutcells
{

template <std::floating_point T, std::integral I = int>
struct MeshView
{
  int gdim = 0;
  int tdim = 0;

  // Flattened coordinates: [x0_0, x0_1, ..., x1_0, x1_1, ...]
  std::span<const T> coordinates;

  // CSR-like cell storage
  std::span<const I> connectivity;
  std::span<const I> offsets;

  // Cell types, one per cell
  std::span<const I> cell_types;

  // Optional keep-alive anchor for Python / external memory
  std::shared_ptr<void> owner;

  I num_nodes() const
  {
    if (gdim <= 0)
      return 0;
    return static_cast<I>(coordinates.size() / static_cast<std::size_t>(gdim));
  }

  I num_cells() const
  {
    if (offsets.empty())
      return 0;
    return static_cast<I>(offsets.size() - 1);
  }

  bool has_cell_types() const
  {
    return !cell_types.empty();
  }

  const T* node(I node_id) const
  {
    return coordinates.data() + static_cast<std::size_t>(node_id) * static_cast<std::size_t>(gdim);
  }

  I cell_num_nodes(I cell_id) const
  {
    return offsets[static_cast<std::size_t>(cell_id) + 1] - offsets[static_cast<std::size_t>(cell_id)];
  }

  std::span<const I> cell_nodes(I cell_id) const
  {
    const std::size_t begin = static_cast<std::size_t>(offsets[static_cast<std::size_t>(cell_id)]);
    const std::size_t end = static_cast<std::size_t>(offsets[static_cast<std::size_t>(cell_id) + 1]);
    return connectivity.subspan(begin, end - begin);
  }

  I cell_node(I cell_id, I local_id) const
  {
    const std::size_t pos = static_cast<std::size_t>(offsets[static_cast<std::size_t>(cell_id)])
                          + static_cast<std::size_t>(local_id);
    return connectivity[pos];
  }

  I cell_type(I cell_id) const
  {
    if (cell_types.empty())
      throw std::runtime_error("MeshView::cell_type not available");
    return cell_types[static_cast<std::size_t>(cell_id)];
  }
};

} // namespace cutcells