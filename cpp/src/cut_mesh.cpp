// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include "cut_mesh.h"

namespace cutcells::mesh
{
    void str(CutCells &cut_mesh)
    {
        std::cout << "CutCells: " << std::endl;
        int cnt = 0;

        for(auto &cell: cut_mesh._cut_cells)
        {
            std::cout << "Cut Cell ";
            std::cout << cnt << ": ";
            cnt++;
            std::cout << "vertex coordinates=[";
            for(auto &i: cell._vertex_coords)
            {
                std::cout << i << ", ";
            }
            std::cout << "]" << std::endl;
            std::cout << "connectivity=[";
            for(int i=0;i<cell._connectivity.size();i++)
            {
                std::cout << i << ": ";
                for(int j=0;j<cell._connectivity[i].size();j++)
                {
                    std::cout << cell._connectivity[i][j] << ", ";
                }
            }
            std::cout << "]" << std::endl;
        }

        std::cout << "Cell Types= ";
        for(auto &type : cut_mesh._types)
        {
            std::cout << cell_type_to_str(type) << ", ";
        }
        std::cout << std::endl;

        std::cout << "Parent map=[";
        for(int i=0;i<cut_mesh._parent_map.size();i++)
        {
            std::cout << i << ": " << cut_mesh._parent_map[i] << std::endl;
        }
        std::cout << "]" << std::endl;
    }

    std::unordered_map<int, std::vector<int>> create_parent_cut_cells_map(const std::span<int> parent_map)
    {
        std::unordered_map<int, std::vector<int>> parent_cut_cell_map;

        for(std::size_t i=0; i< parent_map.size();i++)
        {
            parent_cut_cell_map[parent_map[i]].push_back(i);
        }

        return parent_cut_cell_map;
    }
}