// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    LGPL-3.0-or-later
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
        // Decompositon into interior, exterior and interface sub-elements of a cell 
        // given by the intersection of a background cell with a level set function
        struct CutCell
        {
            // Geometric Dimension of Cell 
            int _gdim; 

            // Topological Dimension of Cell
            int _tdim;

            // Coordinates of vertices of cut cell
            std::vector<double> _vertex_coords; 

           // Vertex ids of cut cells
            std::vector<std::vector<int>> _connectivity;
            // Cell type of cut cells
            std::vector<type> _types;
        
        };
        
        // One sub cell of a cut 
        // struct SubCell
        // {
        //     // Geometric Dimension of Cell 
        //     int _gdim; 
        //     // Coordinates of cell
        //     std::vector<double> _coords; 
        //     // Type of cell
        //     type _sub_cell_type;
        // };

        void str(const CutCell &cut_cell);  

        void cut(const type cell_type, const std::span<const double> vertex_coordinates, const int gdim, 
                 const std::span<const double> ls_values, const std::string& cut_type_str,
                CutCell& cut_cell, bool triangulate=false);

        void cut(const type cell_type, const std::span<const double> vertex_coordinates, const int gdim, 
             const std::span<const double> ls_values, const std::vector<std::string>& cut_type_str,
             std::vector<CutCell>& cut_cell, bool triangulate=false);

        // void extract_sub_cell(const CutCellDecomposition& cut_cell, const std::string cut_type_str, int sub_cell_index, SubCell& sub_cell);

        //void extract_cut_cell_part(const CutCell& cut_cell, const std::string cut_type_str, CutCellPart& cut_cell_part);
    }

}