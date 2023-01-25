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
    namespace tetrahedron
    {
        int get_num_intersection_points(const int &flag);

        int get_num_sub_elements(const int &flag, bool triangulate);

        // get interface cut
        void cut(const std::span<const double> vertex_coordinates, const int gdim, 
                 const std::span<const double> ls_values, const std::string& cut_type_str,
                 CutCell& cut_cell, bool triangulate);
    }
}

