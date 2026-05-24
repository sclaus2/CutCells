// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <cstddef>
#include <functional>
#include <memory>
#include <span>
#include <stdexcept>
#include <vector>

#include "cell_types.h"
#include "row_view.h"

namespace cutcells
{

template <std::integral I>
struct CellTypeBlock
{
  I first_cell = 0;
  I cell_count = 0;
  cell::type type = cell::type::point;

  bool contains(I cell_id) const
  {
    return cell_id >= first_cell && cell_id < first_cell + cell_count;
  }
};

template <std::floating_point T, std::integral I = int>
struct MeshView
{
  int gdim = 0;
  int tdim = 0;

  // Flattened coordinates: [x0_0, x0_1, ..., x1_0, x1_1, ...]
  std::span<const T> coordinates;

  // Number of scalars between consecutive coordinate rows. Defaults to gdim.
  // This allows views of DOLFINx geometry storage, which is always padded to
  // xyz rows even when the geometric dimension is 1 or 2.
  std::size_t coordinate_stride = 0;

  // CSR-like cell storage
  std::span<const I> connectivity;
  std::span<const I> offsets;

  // Optional fixed-width cell storage metadata. When offsets is empty,
  // connectivity is interpreted as rows of cell_stride entries and the first
  // cell_width entries of each row are the parent vertices.
  I cell_count = 0;
  I cell_width = 0;
  I cell_stride = 0;
  std::span<const StridedRowBlock<I>> cell_blocks;

  // Optional callbacks for layouts that cannot expose a direct contiguous row.
  // The cutting pipeline passes scratch storage to avoid hidden allocations.
  std::function<I(I)> cell_row_size;
  std::function<I(I, I)> cell_row_value;
  std::function<void(I, std::vector<I>&)> cell_row_gather;

  // Cell types, one per cell in cutcells cell::type enum numbering (not VTK codes).
  std::span<const cell::type> cell_types;
  std::span<const CellTypeBlock<I>> cell_type_blocks;

  // Optional uniform cell type for meshes with one topology type.
  cell::type uniform_cell_type = cell::type::point;
  bool has_uniform_cell_type = false;

  // True when the vertex connectivity of this mesh uses VTK vertex ordering.
  // Set by the Python boundary when constructing from VTK data.
  // Used internally to apply the VTK→Basix vertex permutation.
  bool vtk_vertex_order = false;

  // Optional keep-alive anchor for Python / external memory
  std::shared_ptr<void> owner;

  RowAccess<I> cell_row_access() const
  {
    RowAccess<I> rows;
    rows.values = connectivity;
    rows.offsets = offsets;
    rows.row_count = cell_count;
    rows.row_width = cell_width;
    rows.row_stride = cell_stride;
    rows.blocks = cell_blocks;
    rows.row_size_fn = cell_row_size;
    rows.row_value_fn = cell_row_value;
    rows.row_gather_fn = cell_row_gather;
    if (rows.row_width <= 0 && has_uniform_cell_type)
      rows.row_width = static_cast<I>(cell::get_num_vertices(uniform_cell_type));
    return rows;
  }

  I num_nodes() const
  {
    if (gdim <= 0)
      return 0;
    const std::size_t stride = coordinate_stride == 0
                                   ? static_cast<std::size_t>(gdim)
                                   : coordinate_stride;
    if (stride == 0)
      return 0;
    return static_cast<I>(coordinates.size() / stride);
  }

  I num_cells() const
  {
    return cell_row_access().num_rows();
  }

  bool has_cell_types() const
  {
    return !cell_types.empty() || !cell_type_blocks.empty()
        || has_uniform_cell_type;
  }

  const T* node(I node_id) const
  {
    const std::size_t stride = coordinate_stride == 0
                                   ? static_cast<std::size_t>(gdim)
                                   : coordinate_stride;
    return coordinates.data() + static_cast<std::size_t>(node_id) * stride;
  }

  I cell_num_nodes(I cell_id) const
  {
    if (!cell_row_size && cell_width <= 0 && cell_blocks.empty()
        && offsets.empty() && has_cell_types())
      return static_cast<I>(cell::get_num_vertices(cell_type(cell_id)));
    return cell_row_access().row_size(cell_id);
  }

  std::span<const I> cell_nodes(I cell_id) const
  {
    if (!cell_row_size && cell_width <= 0 && cell_blocks.empty()
        && offsets.empty() && !cell_types.empty() && !connectivity.empty())
    {
      const I width = cell_num_nodes(cell_id);
      const I stride = cell_stride > 0 ? cell_stride : width;
      const std::size_t begin
          = static_cast<std::size_t>(cell_id) * static_cast<std::size_t>(stride);
      return std::span<const I>(connectivity.data() + begin,
                                static_cast<std::size_t>(width));
    }
    return cell_row_access().direct_row(cell_id);
  }

  std::span<const I> cell_nodes(I cell_id, std::vector<I>& scratch) const
  {
    if (!cell_row_size && cell_width <= 0 && cell_blocks.empty()
        && offsets.empty() && !cell_types.empty() && !connectivity.empty())
      return cell_nodes(cell_id);
    return cell_row_access().row(cell_id, scratch);
  }

  I cell_node(I cell_id, I local_id) const
  {
    return cell_row_access().value(cell_id, local_id);
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
    throw std::runtime_error("MeshView::cell_type not available");
  }
};

} // namespace cutcells
