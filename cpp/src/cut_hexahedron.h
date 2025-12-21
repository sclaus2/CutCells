// Copyright (c) 2025 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include "cut_cell.h"
#include "cell_types.h"

#include <concepts>
#include <span>
#include <string>
#include <vector>

namespace cutcells::cell::hexahedron
{
    template <std::floating_point T>
    void cut(const std::span<const T> vertex_coordinates, const int gdim,
             const std::span<const T> ls_values, const std::string& cut_type_str,
             CutCell<T>& cut_cell, bool triangulate = false);

    template <std::floating_point T>
    void cut(const std::span<const T> vertex_coordinates, const int gdim,
             const std::span<const T> ls_values, const std::vector<std::string>& cut_type_str,
             std::vector<CutCell<T>>& cut_cell, bool triangulate = false);
}
