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
#include "cut_cell.h"

namespace cutcells::mesh
{
    // QuadratureRule on reference element for CutCell
    struct QuadratureRule
    {
        std::vector<double> _points; 
        std::vector<double> _weights;
    };

    // Collection of cut cells that have been cut with regards to a parent mesh/entities
    struct CutMesh
    {
            // vector of all cut cells 
            std::vector<cell::CutCell> _cut_cells;
            // map of cut cell id to parent cell id 
            // this vector contains all cells that are cut
            std::vector<int> _parent_map;
            // the types of elements contained in all cut_cells
            std::vector<cell::type> _types;
            // vector that are to be filled with quadrature rules in reference coordinates for each cut_cell
            std::vector<QuadratureRule> _quadrature_rules;

            //Total number of vertices 
            int _num_vertices;
    };

    void str(CutMesh &cut_mesh);

    // inverse map of parent_element_map
    std::unordered_map<int, std::vector<int>> create_parent_cut_cells_map(std::span<int> parent_element_map);
}