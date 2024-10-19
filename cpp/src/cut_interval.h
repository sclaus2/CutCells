// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <span>
#include <vector>

#include "cut_cell.h"

namespace cutcells::cell
{
    namespace interval{
        // cut interval by linear interpolation to obtain intersection point
        template <std::floating_point T>
        void compute_intersection_point(const T &level, std::span<const T> p0, std::span<const T> p1,
                 const T& v0, const T& v1, std::vector<T>& intersection_point, const int & offset=0);

        template <std::floating_point T>
        void cut(const std::span<const T> vertex_coordinates, const int gdim,
                const std::span<const T> ls_values, const std::string& cut_type_str,
                CutCell<T>& cut_cell);

        template <std::floating_point T>
        T volume(const std::span<const T> vertex_coordinates, const int gdim);
    }
}
