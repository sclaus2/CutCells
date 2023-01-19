// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    LGPL-3.0-or-later
#pragma once

#include "cell_types.h"

namespace cutcells::cell
{
    inline void triangulation(const type cell_type, int* vertices, std::vector<std::vector<int>>& tris)
    {
        switch(cell_type)
        {
            case type::quadrilateral: tris.resize(2, std::vector<int>(3));
                                      tris = {{vertices[0],vertices[1],vertices[2]}, {vertices[0],vertices[2],vertices[3]}};
                                      break;
            case type::prism:   // Tetrahedron 0 (original vertices): { 0, 2, 1, 3 }
                                // Tetrahedron 1 (original vertices): { 1, 3, 5, 4 }
                                // Tetrahedron 2 (original vertices): { 1, 2, 5, 3 }
                                tris.resize(3, std::vector<int>(3));
                                tris = {{vertices[0],vertices[2],vertices[1],vertices[3]}, 
                                        {vertices[1],vertices[3],vertices[5],vertices[4]},
                                        {vertices[1],vertices[2],vertices[5],vertices[3]}};
                                break;

            default: throw std::invalid_argument("triangulation not implemented for given cell type");
                    break;
        }
    };
}
