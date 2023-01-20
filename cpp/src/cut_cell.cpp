// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include "cut_cell.h"
#include "cut_tetrahedron.h"
#include "cut_triangle.h"
#include "cut_interval.h"
#include "cell_flags.h"

#include <cassert>
#include <set>
#include <unordered_map>


namespace cutcells::cell{

    void str(const CutCell &cut_cell)
    {
        std::cout << "CutCell: " << std::endl;

        std::cout << "vertex coordinates=[";
        for(auto &i: cut_cell._vertex_coords)
        {
            std::cout << i << ", ";
        }
        std::cout << "]" << std::endl;

        std::cout << "connectivity=[";
        for(int i=0;i<cut_cell._connectivity.size();i++)
        {
            std::cout << i << ": ";
            for(int j=0;j<cut_cell._connectivity[i].size();j++)
            {
                    std::cout << cut_cell._connectivity[i][j]; 
                    if(j<cut_cell._connectivity[i].size()-1) 
                        std::cout << ", ";
            }
            if(i<cut_cell._connectivity.size()-1)
                std::cout << std::endl;
        }
        std::cout << "]" << std::endl;

        std::cout << "cell types=[";
        for(int i=0;i<cut_cell._types.size();i++)
        {
            std::cout << i << ": " << cell_type_to_str(cut_cell._types[i]); 
            if(i<cut_cell._types.size()-1)
                std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }

    /// Cut a cell of type cell_type with the level set values in its vertices and 
    /// return a list of elements and their cell types
    void cut(const type cell_type, const std::span<const double> vertex_coordinates, const int gdim, 
             const std::span<const double> ls_values, const std::string& cut_type_str,
             CutCell& cut_cell, bool triangulate)
    {
        switch(cell_type)
        {
            case type::interval: interval::cut(vertex_coordinates, gdim, ls_values, cut_type_str, cut_cell);
                                 break;
            case type::triangle: triangle::cut(vertex_coordinates, gdim, ls_values, cut_type_str, cut_cell, triangulate);
                                 break;
            case type::tetrahedron: tetrahedron::cut(vertex_coordinates, gdim, ls_values, cut_type_str, cut_cell, triangulate);
                                 break;
        }
    }


    void cut(const type cell_type, const std::span<const double> vertex_coordinates, const int gdim, 
             const std::span<const double> ls_values, const std::vector<std::string>& cut_type_str,
             std::vector<CutCell>& cut_cell, bool triangulate)
    {
        switch(cell_type)
        {
            case type::triangle: triangle::cut(vertex_coordinates, gdim, ls_values, cut_type_str, cut_cell, triangulate);
                                 break;
        }
    }

}