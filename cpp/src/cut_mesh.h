// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    LGPL-3.0-or-later
#pragma once

#include <iostream>
#include <span>
#include <vector>
#include <unordered_map>

#include "cell_types.h"

namespace cutcells::mesh
{
    // Collection of cut cells that have been cut with regards to a parent mesh/entities
    struct CutMesh
    {
            std::vector<double> _vertex_coordinates;
            std::vector<std::vector<int>> _cut_cells;
            std::vector<cell::type> _cell_types;
            std::vector<int> _parent_element_map;
    };

    void str(CutMesh &cut_mesh);

    // inverse map of parent_element_map
    std::unordered_map<int, std::vector<int>> create_parent_cut_cells_map(std::span<int> parent_element_map);
}