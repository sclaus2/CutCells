// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    LGPL-3.0-or-later
#pragma once

#include <span>
#include <vector>

#include "cut_cell.h"

namespace cutcells::cell
{
    namespace interval{
        // cut interval by linear interpolation to obtain intersection point
        void compute_intersection_point(const double &level, const std::span<const double> p0, const std::span<const double> p1,
                 const double& v0, const double& v1, std::vector<double>& intersection_point, const int & offset=0);

        void cut(const std::span<const double> vertex_coordinates, const int gdim, 
                const std::span<const double> ls_values, const std::string& cut_type_str,
                CutCell& cut_cell);
    }
}
