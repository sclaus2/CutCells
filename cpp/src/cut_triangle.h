// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <span>
#include <vector>

#include "cell_types.h"
#include "cut_cell.h"

namespace cutcells::cell
{
    namespace triangle
    {
        int get_num_intersection_points(const int &flag);

        int get_num_sub_elements(const int &flag, bool triangulate);

        // get interface cut
        template <std::floating_point T>
        void cut(const std::span<const T> vertex_coordinates, const int gdim,
                 const std::span<const T> ls_values, const std::string& cut_type_str,
                 CutCell<T>& cut_cell, bool triangulate);

        template <std::floating_point T>
        void cut(const std::span<const T> vertex_coordinates, const int gdim,
             const std::span<const T> ls_values, const std::vector<std::string>& cut_type_str,
             std::vector<CutCell<T>>& cut_cell, bool triangulate);

        template <std::floating_point T>
        T volume(const std::span<const T> vertex_coordinates, const int gdim);
    }
}

