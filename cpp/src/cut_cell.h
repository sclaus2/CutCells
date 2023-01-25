// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <iostream>

#include "cell_types.h"
#include <vector>
#include <span>
#include <string>

namespace cutcells
{
    namespace cell
    {
        /// @brief Stores the sub-cells resulting from cutting a cell
        struct CutCell
        {
            /// Geometric Dimension of Cell 
            int _gdim; 

            /// Topological Dimension of Cell
            int _tdim;

            /// Coordinates of vertices of cut cell
            std::vector<double> _vertex_coords; 

            /// Vertex ids of cut cells
            std::vector<std::vector<int>> _connectivity;

            /// Cell type of cut cells
            std::vector<type> _types;
        
        };
        
        void str(const CutCell &cut_cell);  

        void cut(const type cell_type, const std::span<const double> vertex_coordinates, const int gdim, 
                 const std::span<const double> ls_values, const std::string& cut_type_str,
                CutCell& cut_cell, bool triangulate=false);

        void cut(const type cell_type, const std::span<const double> vertex_coordinates, const int gdim, 
             const std::span<const double> ls_values, const std::vector<std::string>& cut_type_str,
             std::vector<CutCell>& cut_cell, bool triangulate=false);
    }

}