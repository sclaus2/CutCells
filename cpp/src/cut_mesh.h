// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <iostream>
#include <span>
#include <vector>
#include <unordered_map>

#include "cell_types.h"
#include "cut_cell.h"

namespace cutcells::mesh
{
    /// Collection of cut cells that have been cut with regards to a parent mesh/entities
    struct CutCells
    {
            /// vector of all cut cells 
            std::vector<cell::CutCell> _cut_cells;
            /// map of cut cell id to parent cell id 
            /// this vector contains all cells that are cut
            std::vector<std::int32_t> _parent_map;
            /// the types of elements contained in all cut_cells
            std::vector<cell::type> _types;

            /// Total number of vertices in CutMesh
            std::size_t _num_vertices;
    };

    /// Collection of cut cells that have been cut with regards to a parent mesh/entities
    struct CutInterface
    {
            /// vector of all cut cells 
            std::vector<cell::CutCell> _cut_cells;
            /// map of cut cell id to parent cell ids (2 for interface)
            /// this vector contains all cells that are cut
            std::vector<std::pair<std::int32_t, std::int32_t>> _parent_map;

            /// the types of elements contained in all cut_cells
            std::vector<cell::type> _types;

            /// Total number of vertices in CutMesh
            std::size_t _num_vertices;
    };

    /// @brief  Print information about cut_mesh to screen
    /// @param cut_mesh 
    void str(CutCells &cut_mesh);

    /// @brief Get inverse map of parent_map 
    /// @param parent_map: map from cutcell index to parent cell index. 
    /// @return inverse map from parent cell index to cut cell index 
    std::unordered_map<int, std::vector<int>> create_parent_cut_cells_map(std::span<int> parent_map);

    //Get number of cells in CutMesh
    int get_num_cells(const cutcells::mesh::CutCells& cut_mesh);
}