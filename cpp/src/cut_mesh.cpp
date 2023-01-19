// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include "cut_mesh.h"

namespace cutcells::mesh
{
    void str(CutMesh &cut_mesh)
    {
        std::cout << "CutMesh: " << std::endl;

        std::cout << "vertex coordinates=[";
        for(auto &i: cut_mesh._vertex_coordinates)
        {
            std::cout << i << ", ";
        }
        std::cout << "]" << std::endl;

        std::cout << "cut cells=[";
        for(int i=0;i<cut_mesh._cut_cells.size();i++)
        {
            std::cout << i << ": ";
            for(int j=0;j<cut_mesh._cut_cells[i].size();j++)
            {
                    std::cout << cut_mesh._cut_cells[i][j] << ", ";
            }
            std::cout << std::endl;
        }
        std::cout << "]" << std::endl;

        std::cout << "Parent map=[";
        for(int i=0;i<cut_mesh._parent_element_map.size();i++)
        {
            std::cout << i << ": " << cut_mesh._parent_element_map[i] << std::endl;
        }
        std::cout << "]" << std::endl;
    }

    // inverse map of parent_element_map
    std::unordered_map<int, std::vector<int>> create_parent_cut_cells_map(const std::span<int> parent_element_map)
    {
        std::unordered_map<int, std::vector<int>> parent_cut_cell_map;

        for(std::size_t i=0; i< parent_element_map.size();i++)
        {
            parent_cut_cell_map[parent_element_map[i]].push_back(i);
        }

        return parent_cut_cell_map;
    }
}