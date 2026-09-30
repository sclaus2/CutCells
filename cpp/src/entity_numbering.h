// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <vector>

namespace cutcells::detail
{

/// Number visited local entities by first occurrence of their vertex key.
/// On return `entity_of[k]` is the id of the k-th visited entity (ids are
/// assigned in visiting order) and `is_first[k]` marks the visit that created
/// the entity. Returns the number of distinct entities.
template <typename Key>
int number_by_first_occurrence(const std::vector<Key>& keys,
                               std::vector<std::int32_t>& entity_of,
                               std::vector<std::uint8_t>& is_first,
                               std::vector<std::int32_t>& order)
{
    const std::size_t n = keys.size();
    order.resize(n);
    std::iota(order.begin(), order.end(), std::int32_t(0));
    std::sort(order.begin(), order.end(),
              [&keys](std::int32_t a, std::int32_t b)
              {
                  const auto& ka = keys[static_cast<std::size_t>(a)];
                  const auto& kb = keys[static_cast<std::size_t>(b)];
                  return ka < kb || (ka == kb && a < b);
              });

    // entity_of temporarily holds the first visit sharing the key.
    entity_of.resize(n);
    for (std::size_t i = 0; i < n;)
    {
        const std::int32_t first = order[i];
        std::size_t j = i;
        while (j < n
               && keys[static_cast<std::size_t>(order[j])]
                      == keys[static_cast<std::size_t>(first)])
        {
            entity_of[static_cast<std::size_t>(order[j])] = first;
            ++j;
        }
        i = j;
    }

    is_first.assign(n, std::uint8_t(0));
    int count = 0;
    for (std::size_t k = 0; k < n; ++k)
    {
        const std::int32_t first = entity_of[k];
        if (first == static_cast<std::int32_t>(k))
        {
            is_first[k] = 1;
            entity_of[k] = count++;
        }
        else
            entity_of[k] = entity_of[static_cast<std::size_t>(first)];
    }
    return count;
}

template <typename Key>
int number_by_first_occurrence(const std::vector<Key>& keys,
                               std::vector<std::int32_t>& entity_of,
                               std::vector<std::uint8_t>& is_first)
{
    std::vector<std::int32_t> order;
    return number_by_first_occurrence(keys, entity_of, is_first, order);
}

} // namespace cutcells::detail
