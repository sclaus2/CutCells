// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <string>
#include "cut_cell.h"

namespace cutcells::io
{
    void write_vtk(std::string filename, const std::span<const double> element_vertex_coords,  
                    const std::vector<std::vector<int>> elements,
                    const std::span<cell::type> element_types, 
                    const int gdim);

    void write_vtk(std::string filename, cell::CutCell<double>& cut_cell);
}