// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    LGPL-3.0-or-later
#pragma once

namespace cutcells
{
    namespace cell
    {
        /// Cell type
        enum class type
        {
            point = 0,
            interval = 1,
            triangle = 2,
            tetrahedron = 3,
            quadrilateral = 4,
            hexahedron = 5,
            prism = 6,
            pyramid = 7
        };

        enum class vtk_types
        {
            VTK_EMPTY_CELL = 0,
            VTK_VERTEX = 1,
            VTK_POLY_VERTEX = 2,
            VTK_LINE = 3,
            VTK_POLY_LINE = 4,
            VTK_TRIANGLE = 5,
            VTK_TRIANGLE_STRIP = 6,
            VTK_POLYGON = 7,
            VTK_PIXEL = 8,
            VTK_QUAD = 9,
            VTK_TETRA = 10,
            VTK_VOXEL = 11,
            VTK_HEXAHEDRON = 12,
            VTK_WEDGE = 13,
            VTK_PYRAMID = 14,
            VTK_PENTAGONAL_PRISM = 15,
            VTK_HEXAGONAL_PRISM = 16,
        };

        inline vtk_types map_cell_type_to_vtk(type cell_type)
        {
            vtk_types vtk_type; 

            switch(cell_type)
            {
                case type::point: vtk_type = vtk_types::VTK_VERTEX;
                                  break;
                case type::interval: vtk_type=vtk_types::VTK_LINE;
                                     break;
                case type::triangle: vtk_type=vtk_types::VTK_TRIANGLE;
                                     break;
                case type::quadrilateral: vtk_type=vtk_types::VTK_QUAD;
                                          break;
                case type::tetrahedron: vtk_type=vtk_types::VTK_TETRA;
                                     break;
                case type::prism: vtk_type=vtk_types::VTK_WEDGE;
                                          break;
            }

            return vtk_type;
        };

        inline int get_num_vertices(type cell_type)
        {
            int num_vertices = 0;
            switch(cell_type)
            {
                case type::interval:      num_vertices = 2;
                                          break;
                case type::triangle:      num_vertices = 3;
                                          break;
                case type::quadrilateral: num_vertices = 4;
                                          break;
                case type::tetrahedron:   num_vertices = 4;
                                          break;
                case type::prism:         num_vertices = 6;
                                          break;
            }

            return num_vertices;
        }
    }
}