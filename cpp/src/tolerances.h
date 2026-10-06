// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <limits>

namespace cutcells
{

/// The one tolerance of touching zero, relative to a level set's magnitude (its
/// largest Bernstein coefficient on a cell, say): values within it count as 0.
/// Classification (quadrays::coefficient_sign), the lookup tables' values and
/// the decision that a level set vanishes on a facet, edge or vertex of the
/// mesh (part::cut) all read it.
template <std::floating_point T>
inline constexpr T touch_tolerance = T(64) * std::numeric_limits<T>::epsilon();

} // namespace cutcells
