// Copyright (c) 2022 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <string>
#include <stdexcept>

#include <vector>

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
                case type::point:           vtk_type = vtk_types::VTK_VERTEX;
                                            break;
                case type::interval:        vtk_type=vtk_types::VTK_LINE;
                                            break;
                case type::triangle:        vtk_type=vtk_types::VTK_TRIANGLE;
                                            break;
                case type::quadrilateral:   vtk_type=vtk_types::VTK_QUAD;
                                            break;
                case type::tetrahedron:     vtk_type=vtk_types::VTK_TETRA;
                                            break;
                case type::hexahedron:      vtk_type=vtk_types::VTK_HEXAHEDRON;
                                            break;
                case type::prism:           vtk_type=vtk_types::VTK_WEDGE;
                                            break;
                case type::pyramid:         vtk_type=vtk_types::VTK_PYRAMID;
                                            break;
                default: throw std::invalid_argument("cell type not recognised in map_cell_type_to_vtk of cell_types.h");
                                            break;
            }

            return vtk_type;
        };

        inline type map_vtk_type_to_cell_type(vtk_types vtk_type)
        {
            type cell_type;

            switch(vtk_type)
            {
                case vtk_types::VTK_VERTEX: cell_type = type::point;
                                            break;
                case vtk_types::VTK_LINE:   cell_type = type::interval;
                                            break;
                case vtk_types::VTK_TRIANGLE: cell_type = type::triangle;
                                            break;
                case vtk_types::VTK_QUAD:   cell_type =  type::quadrilateral;
                                            break;
                case vtk_types::VTK_TETRA:  cell_type = type::tetrahedron;
                                            break;
                case vtk_types::VTK_HEXAHEDRON: cell_type = type::hexahedron;
                                            break;
                case vtk_types::VTK_WEDGE:  cell_type = type::prism;
                                            break;
                case  vtk_types::VTK_PYRAMID: cell_type = type::pyramid;
                                              break;
                default: throw std::invalid_argument("vtk cell type not recognised in map_vtk_type_to_cell_type of cell_types.h");
                                            break;
            }

            return cell_type;
        };

        inline int get_num_vertices(type cell_type)
        {
            int num_vertices = 0;
            switch(cell_type)
            {
                case type::point:         num_vertices = 1;
                                          break;
                case type::interval:      num_vertices = 2;
                                          break;
                case type::triangle:      num_vertices = 3;
                                          break;
                case type::quadrilateral: num_vertices = 4;
                                          break;
                case type::tetrahedron:   num_vertices = 4;
                                          break;
                case type::hexahedron:    num_vertices = 8;
                                          break;
                case type::prism:         num_vertices = 6;
                                          break;
                case type::pyramid:       num_vertices = 5;
                                          break;
                default: throw std::invalid_argument("cell type not recognised in get_num_vertices of cell_types.h");
                         break;
            }

            return num_vertices;
        }

        inline int get_tdim(type cell_type)
        {
            int tdim = 0;
            switch(cell_type)
            {
                case type::point:         tdim = 0;
                                          break;
                case type::interval:      tdim = 1;
                                          break;
                case type::triangle:      tdim = 2;
                                          break;
                case type::quadrilateral: tdim = 2;
                                          break;
                case type::tetrahedron:   tdim = 3;
                                          break;
                case type::hexahedron:    tdim = 3;
                                          break;
                case type::prism:         tdim = 3;
                                          break;
                case type::pyramid:       tdim = 3;
                                          break;
                default: throw std::invalid_argument("cell type not recognised in get_num_vertices of cell_types.h");
                         break;
            }

            return tdim;
        }

        inline std::string cell_type_to_str(type cell_type)
        {
            std::string type_str;
            switch(cell_type)
            {
                case type::point:         type_str = "point";
                                          break;
                case type::interval:      type_str = "interval";
                                          break;
                case type::triangle:      type_str = "triangle";
                                          break;
                case type::quadrilateral: type_str = "quadrilateral";
                                          break;
                case type::tetrahedron:   type_str = "tetrahedron";
                                          break;
                case type::hexahedron:    type_str = "hexahedron";
                                          break;
                case type::prism:         type_str = "prism";
                                          break;
                case type::pyramid:       type_str = "pyramid";
                                          break;
                default: throw std::invalid_argument("cell type not recognised in cell_type_to_str of cell_types.h");
                         break;
            }
            return type_str;
        }

        //get possible cell types from type of cut for cell type cell_type
        inline std::vector<type> cut_cell_types(type cell_type, const std::string& ls_part)
        {
            std::vector<type> cell_types;

            switch(cell_type)
            {
                case type::point:         {cell_types.resize(1);
                                          cell_types[0]= type::point;
                                          break;}
                case type::interval:      {if(ls_part=="phi=0"){
                                            cell_types.resize(1);
                                            cell_types[0]= type::point;}
                                          else{
                                            cell_types.resize(1);
                                            cell_types[0]= type::interval;}
                                            break;}
                case type::triangle:      {if(ls_part=="phi=0")
                                          {
                                            cell_types.resize(1);
                                            cell_types[0]=type::interval;
                                          }
                                          else
                                          {
                                            cell_types.resize(2);
                                            cell_types[0]= type::triangle;
                                            cell_types[1]= type::quadrilateral;
                                          }
                                          break;}
                case type::quadrilateral: {throw std::invalid_argument("quadrilateral yet not implemented");
                                          break;}
                case type::tetrahedron:   {if(ls_part=="phi=0")
                                          {
                                            cell_types.resize(2);
                                            cell_types[0]= type::triangle;
                                            cell_types[1]= type::quadrilateral;
                                          }
                                          else
                                          {
                                            cell_types.resize(2);
                                            cell_types[0]= type::tetrahedron;
                                            cell_types[1]= type::prism;
                                          }
                                          break;}
                case type::hexahedron:    {throw std::invalid_argument("hexahedron yet not implemented");
                                          break;}
                case type::prism:         {throw std::invalid_argument("prism yet not implemented");
                                          break;}
                case type::pyramid:       {throw std::invalid_argument("pyramid yet not implemented");
                                          break;}
                default: {throw std::invalid_argument("cell type not recognised in cell_type_to_str of cell_types.h");
                         break;}
            }
            return cell_types;
        }
    }
}
