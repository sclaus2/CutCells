// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#pragma once

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <functional>
#include <span>
#include <stdexcept>
#include <vector>

namespace cutcells
{

enum class RowLayout
{
  empty,
  csr,
  strided,
  blocked_strided,
  callback
};

template <std::integral I>
struct StridedRowBlock
{
  I first_row = 0;
  I row_count = 0;
  I row_width = 0;
  I row_stride = 0;
  std::span<const I> values;

  bool contains(I row) const
  {
    return row >= first_row && row < first_row + row_count;
  }

  I width() const
  {
    return row_width > 0 ? row_width : row_stride;
  }

  I stride() const
  {
    return row_stride > 0 ? row_stride : width();
  }

  std::span<const I> row(I row) const
  {
    const I local_row = row - first_row;
    const I w = width();
    const I s = stride();
    const std::size_t begin
        = static_cast<std::size_t>(local_row) * static_cast<std::size_t>(s);
    return std::span<const I>(values.data() + begin, static_cast<std::size_t>(w));
  }
};

template <std::integral I>
struct RowAccess
{
  std::span<const I> values;
  std::span<const I> offsets;

  I row_count = 0;
  I row_width = 0;
  I row_stride = 0;

  std::span<const StridedRowBlock<I>> blocks;

  std::function<I(I)> row_size_fn;
  std::function<I(I, I)> row_value_fn;
  std::function<void(I, std::vector<I>&)> row_gather_fn;

  RowLayout layout() const
  {
    if (!offsets.empty())
      return RowLayout::csr;
    if (!blocks.empty())
      return RowLayout::blocked_strided;
    if (!values.empty() && (row_count > 0 || row_width > 0 || row_stride > 0))
      return RowLayout::strided;
    if (row_gather_fn || row_value_fn || row_size_fn)
      return RowLayout::callback;
    return RowLayout::empty;
  }

  I num_rows() const
  {
    if (!offsets.empty())
      return static_cast<I>(offsets.size() - 1);
    if (row_count > 0)
      return row_count;
    if (!blocks.empty())
    {
      I n = 0;
      for (const auto& block : blocks)
        n = std::max(n, block.first_row + block.row_count);
      return n;
    }
    const I s = stride();
    if (!values.empty() && s > 0)
      return static_cast<I>(values.size() / static_cast<std::size_t>(s));
    return 0;
  }

  I row_size(I row_id) const
  {
    const RowLayout current_layout = layout();
    if (current_layout == RowLayout::callback)
    {
      if (row_size_fn)
        return row_size_fn(row_id);
      throw std::runtime_error("RowAccess callback layout requires a row_size callback");
    }
    if (!offsets.empty())
    {
      return offsets[static_cast<std::size_t>(row_id) + 1]
           - offsets[static_cast<std::size_t>(row_id)];
    }
    if (!blocks.empty())
      return block(row_id).width();
    if (row_width > 0)
      return row_width;
    return row_stride;
  }

  std::span<const I> direct_row(I row_id) const
  {
    if (!offsets.empty())
    {
      const std::size_t begin
          = static_cast<std::size_t>(offsets[static_cast<std::size_t>(row_id)]);
      const std::size_t end = static_cast<std::size_t>(
          offsets[static_cast<std::size_t>(row_id) + 1]);
      return values.subspan(begin, end - begin);
    }
    if (!blocks.empty())
      return block(row_id).row(row_id);
    if (!values.empty())
    {
      const I s = stride();
      const I w = row_size(row_id);
      const std::size_t begin
          = static_cast<std::size_t>(row_id) * static_cast<std::size_t>(s);
      return std::span<const I>(values.data() + begin,
                                static_cast<std::size_t>(w));
    }
    throw std::runtime_error("RowAccess::direct_row is not available");
  }

  std::span<const I> row(I row_id, std::vector<I>& scratch) const
  {
    switch (layout())
    {
    case RowLayout::csr:
    case RowLayout::strided:
    case RowLayout::blocked_strided:
      return direct_row(row_id);
    case RowLayout::callback:
      scratch.clear();
      if (row_gather_fn)
      {
        row_gather_fn(row_id, scratch);
      }
      else if (row_size_fn && row_value_fn)
      {
        const I n = row_size_fn(row_id);
        scratch.resize(static_cast<std::size_t>(n));
        for (I i = 0; i < n; ++i)
          scratch[static_cast<std::size_t>(i)] = row_value_fn(row_id, i);
      }
      else
      {
        throw std::runtime_error(
            "RowAccess callback layout requires gather or size+value callbacks");
      }
      return std::span<const I>(scratch.data(), scratch.size());
    case RowLayout::empty:
      break;
    }
    throw std::runtime_error("RowAccess::row is not available");
  }

  I value(I row_id, I local_id) const
  {
    switch (layout())
    {
    case RowLayout::csr:
    case RowLayout::strided:
    case RowLayout::blocked_strided:
      return direct_row(row_id)[static_cast<std::size_t>(local_id)];
    case RowLayout::callback:
      if (row_value_fn)
        return row_value_fn(row_id, local_id);
      break;
    case RowLayout::empty:
      break;
    }
    throw std::runtime_error("RowAccess::value is not available");
  }

private:
  I stride() const
  {
    return row_stride > 0 ? row_stride : row_width;
  }

  const StridedRowBlock<I>& block(I row_id) const
  {
    for (const auto& candidate : blocks)
    {
      if (candidate.contains(row_id))
        return candidate;
    }
    throw std::out_of_range("RowAccess row is not covered by any block");
  }
};

} // namespace cutcells
