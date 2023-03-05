// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include "cut_triangle.h"
#include "cut_interval.h"
#include "cell_flags.h"
#include "triangulation.h"

#include <cassert>
#include <stdexcept>
#include <iostream>
#include <unordered_map>

namespace cutcells::cell
{
// Look up tables for intersection
namespace{
    // Choose numbering of edges for cutting in accordance with numbering of vtk
    int edges[3][2] = {{0,1}, {1,2}, {2,0}};

    // List of ids of intersected edges these will be used below to obtain 
    // the intersection points (one intersection point is assumed per edge) 
    int triangle_intersected_edges[8][2] = 
    {   {-1, -1},                    // 0
        {0, 2},                      // 1 
        {1, 0},                      // 2
        {1, 2},                      // 3 
        {2, 1},                      // 4
        {0, 1},                      // 5
        {2, 0},                      // 6
        {-1,-1}                      // 7
    };

    // List of vertex ids that form sub elements created by intersection. 
    // Intersection points are numbered with 0,1,2...
    // original vertices are numbered by 100, 101, ...
    // The triangle vertex and edge numbering used is 
    //       102
    //     /    \
    //    2       1
    //   /         \
    // 100 -- 0 -- 101
    // The list gives back either the triangle or quadrilateral produced by the intersection
    // This list can be used to reconstruct the coordinates of the sub-elements.
    int triangle_sub_element[8][4] = 
    {   {-1, -1, -1 ,-1 },                    // 0
        {0, 2, 100, -1},                      // 1 
        {1, 0, 101, -1},                      // 2
        {1, 2, 100, 101},                     // 3 
        {2, 1, 102, -1},                      // 4
        {0, 1, 102, 100},                     // 5
        {2, 0, 101, 102},                     // 6
        {100, 101, 102, -1}                   // 7
    };

    // List of cell types for each sub-element produced in each cut case
    type triangle_sub_element_cell_types[8] = 
    {   type::triangle,                       // 0
        type::triangle,                       // 1
        type::triangle,                       // 2
        type::quadrilateral,                  // 3
        type::triangle,                       // 4
        type::quadrilateral,                  // 5
        type::quadrilateral,                  // 6
        type::triangle                        // 7
    };
}

namespace triangle{

    int get_num_intersection_points(const int &flag)
    {
        // number of intersection points for intersected triangle is always 2
        return 2; 
    }

    //for triangles the number of sub-elements without triangulation is 1
    int get_num_sub_elements(const int &flag, bool triangulate)
    {   
        type cell_type = triangle_sub_element_cell_types[flag];
        int num_sub_elements = 1;

        if(cell_type == type::quadrilateral && triangulate == true)
        {
            num_sub_elements = 2;
        }
        return num_sub_elements;
    }

    void compute_intersection_points(const std::span<const double> vertex_coordinates, const int gdim, 
             const std::span<const double> ls_values, const int flag, std::vector<double>& intersection_points, 
             std::unordered_map<int,int>& vertex_case_map)
    {
        // vertex ids will be edges[edge_0][0] edges[edge_0][1]
        // get vertex coordinates and interpolate 
        int num_intersection_points = get_num_intersection_points(flag);
        
        intersection_points.resize(num_intersection_points*gdim);

        std::vector<double> v0(gdim);
        std::vector<double> v1(gdim);
        std::vector<double> intersection_point(gdim);

        for(int ip=0; ip<num_intersection_points; ip++)
        {
            // edge has two vertices
            // intersection points are listed first in triangle case table 
            // therefore the first two entries are the intersected edges
            int vertex_id_0 = edges[triangle_intersected_edges[flag][ip]][0];
            int vertex_id_1 = edges[triangle_intersected_edges[flag][ip]][1]; 

            for(int j=0;j<gdim;j++)
            {
                v0[j] = vertex_coordinates[vertex_id_0*gdim+j];
                v1[j] = vertex_coordinates[vertex_id_1*gdim+j];
            }

            double ls0 = ls_values[vertex_id_0];
            double ls1 = ls_values[vertex_id_1];

            interval::compute_intersection_point(0.0, v0, v1,ls0, ls1, intersection_point);

            for(int j=0;j<gdim;j++)
            {
                intersection_points[ip*gdim+j] = intersection_point[j];
            }
        }

        vertex_case_map[triangle_intersected_edges[flag][0]] = 0;
        vertex_case_map[triangle_intersected_edges[flag][1]] = 1;
    }

    void create_sub_cell_vertex_coords(const int& flag, const std::span<const double> vertex_coordinates, const int gdim, 
                                       const std::span<const double> intersection_points, 
                                       std::vector<double>& coords, std::unordered_map<int,int>& vertex_case_map)
    {
        int num_intersection_points = intersection_points.size()/gdim;
        type cell_type = triangle_sub_element_cell_types[flag];
        int total_num_points = get_num_vertices(cell_type);

        // Copy vertex coordinates into cut_cell
        coords.resize(total_num_points*gdim);

        for(int i=0;i<intersection_points.size();i++)
        {
            coords[i] = intersection_points[i];
        }

        int ncols = 4;
        int offset = num_intersection_points; 
        int cnt = num_intersection_points; 
        int k = 0;

        for(int i=0;i<ncols;i++)
        {
           if( triangle_sub_element[flag][i] < 100 || triangle_sub_element[flag][i]==-1)
           {
            //intersection point which is alreay added or is not a vertex
           }
           else 
           {    
                int vertex_id = triangle_sub_element[flag][i];
                vertex_case_map[triangle_sub_element[flag][i]] = cnt;
                cnt++;

                vertex_id -= 100;
                for(int j=0;j<gdim;j++)
                {
                    coords[(offset+k)*gdim+j] = vertex_coordinates[vertex_id*gdim+j];
                }
                k++;
           }
        }
    }

    void create_interface_cells(std::vector<std::vector<int>>& interface_cells, std::vector<type>& interface_cell_types)
    {
        interface_cells.resize(1);
        interface_cells[0].resize(2);
        interface_cells[0][0] = 0;
        interface_cells[0][1] = 1;

        interface_cell_types.resize(1);
        interface_cell_types[0] = type::interval;
    }

    void create_sub_cells(const int& flag, bool triangulate, std::vector<std::vector<int>>& sub_cells, 
                          std::vector<type>& sub_cell_types, std::unordered_map<int,int>& vertex_case_map)
    {
        // Allocate memory
        int num_sub_elements = get_num_sub_elements(flag, triangulate);
        sub_cells.resize(num_sub_elements);
        sub_cell_types.resize(num_sub_elements);

        for(int i=0;i<num_sub_elements;i++)
        {
            type sub_cell_type;

            if(triangulate)
            {
                sub_cell_type = type::triangle;
            }
            else
            {
                sub_cell_type = triangle_sub_element_cell_types[flag];
            }
            sub_cell_types[i] = sub_cell_type;
            int num_vertices = get_num_vertices(sub_cell_type);
            sub_cells[i].resize(num_vertices);
        }

        // Fill in vertex ids
        //int vertices[4] = triangle_sub_element[flag];
        type sub_cell_type = triangle_sub_element_cell_types[flag];

        //collect all vertices of sub-elements (elements are separated by -1 in list)
        if(sub_cell_type == type::quadrilateral && triangulate == true)
        {
            std::vector<std::vector<int>> triangles;
            triangulation(sub_cell_type, triangle_sub_element[flag], triangles);

            for(int i=0;i<num_sub_elements;i++)
                for(int j=0;j<3;j++)
                {
                    sub_cells[i][j] = vertex_case_map[triangles[i][j]];
                }
        }
        else
        {
            //number of entries in sub_element table
            int ncols = 4;
            int sub_element = 0;

            for(int i=0;i<num_sub_elements;i++)
            {
                for(int j=0;j<ncols;j++)
                {
                    if(triangle_sub_element[flag][j]!=-1)
                    {
                        sub_cells[i][j] = vertex_case_map[triangle_sub_element[flag][j]];
                    }
                    else
                    {
                        //vertex is -1, move on to next element
                    }
                } 
            }
        }
    }

    void create_cut_cell(const std::span<const double> vertex_coordinates, 
                         const int gdim,  const std::span<const double> ls_values,
                         const std::string& cut_type_str,
                         CutCell& cut_cell, bool triangulate, 
                         const std::span<const double> intersection_points, 
                         std::unordered_map<int,int>& vertex_case_map)
    {
        cut_cell._gdim = gdim;
        
        if(cut_type_str=="phi=0")
        {
            cut_cell._tdim = 1;
            int flag_interior = get_entity_flag(ls_values, false);
            // Copy vertex coordinates into cut_cell
            int num_intersection_points = get_num_intersection_points(flag_interior);
            cut_cell._vertex_coords.resize(num_intersection_points*gdim);
            for(int i=0;i<intersection_points.size();i++)
            {
                cut_cell._vertex_coords[i] = intersection_points[i];
            }
            // Determine interface cells for triangle this is an interval with two points (the intersection points)
            create_interface_cells(cut_cell._connectivity, cut_cell._types);
        }
        else if(cut_type_str=="phi<0")
        {
            cut_cell._tdim = 2;
            int flag_interior = get_entity_flag(ls_values, false);
            create_sub_cell_vertex_coords(flag_interior, vertex_coordinates, gdim, intersection_points, 
                        cut_cell._vertex_coords, vertex_case_map);
            //Determine interior sub-cells
            create_sub_cells(flag_interior, triangulate, cut_cell._connectivity, 
                        cut_cell._types, vertex_case_map);
        }
        else if(cut_type_str=="phi>0")
        {
            cut_cell._tdim = 2;
            //Determine exterior sub-cells
            int flag_exterior = get_entity_flag(ls_values, true);
            create_sub_cell_vertex_coords(flag_exterior, vertex_coordinates, gdim, intersection_points, 
                        cut_cell._vertex_coords, vertex_case_map);
            create_sub_cells(flag_exterior, triangulate, cut_cell._connectivity, 
                    cut_cell._types, vertex_case_map);
        }
        else
        {
            throw std::invalid_argument("cutting type unknown");
        }
    }

    void str(CutCell& cut_cell, std::unordered_map<int,int>& vertex_case_map)
    {
        std::cout << "vertex case map=[";
        for(auto& it: vertex_case_map)
        {
        std::cout << it.first << ": " << it.second << std::endl;
        }
        std::cout << "]" << std::endl;

        std::cout << "connectivity=[";
        for(auto &i: cut_cell._connectivity)
        {
            for(auto &j : i)
            {
                std::cout << j << ", ";
            }
        }
        std::cout << "]" << std::endl;

        std::cout << "types=[";
        for(auto &i: cut_cell._types)
        {
            std::cout << static_cast<int>(i) << ", ";
        }
        std::cout << "]" << std::endl;

        std::cout << "vertex coordinates=[";
        for(auto &i : cut_cell._vertex_coords)
        {
            std::cout << i << ",";
        }
        std::cout << "]" << std::endl;
    }

    // cut triangle
    void cut(const std::span<const double> vertex_coordinates, const int gdim, 
             const std::span<const double> ls_values, const std::string& cut_type_str,
             CutCell& cut_cell, bool triangulate)
    {
        int flag_interior = get_entity_flag(ls_values, false);

        //throw error if cell is not intersected, only intersected cells should land here 
        if(flag_interior ==0 || flag_interior == 7)
        {
            throw std::invalid_argument("triangle is not intersected and therefore cannot be cut");
        }

        // Compute intersection points these are required for any cut cell part (interface, interior, exterior)
        // get the number of intersection points
        std::vector<double> intersection_points; 
        // the vertex case map, 
        // first few entries map from intersected edge to intersection point number 
        // next entries map from orginal vertex id to number of vertex in vertex_coordinates of CutCell (renumbered to go from 0,...,N)
        // example: intersected edges 0 -> 0, 2 -> 1 
        //          then orginal vertex 101 -> 2 , 102 -> 3 etc.
        std::unordered_map<int,int> vertex_case_map;       
        compute_intersection_points(vertex_coordinates, gdim, ls_values, flag_interior, intersection_points, vertex_case_map);

        create_cut_cell(vertex_coordinates, gdim, ls_values, cut_type_str, cut_cell, 
                        triangulate, intersection_points, vertex_case_map);
    };

    // cut triangle version with vector of string in case multiple parts of the cut-cell are needed (very common)
    void cut(const std::span<const double> vertex_coordinates, const int gdim, 
             const std::span<const double> ls_values, const std::vector<std::string>& cut_type_str,
             std::vector<CutCell>& cut_cell, bool triangulate)
    {
        int flag_interior = get_entity_flag(ls_values, false);

        // throw error if cell is not intersected, only intersected cells should land here 
        if(flag_interior ==0 || flag_interior == 7)
        {
            throw std::invalid_argument("triangle is not registered as intersected and therefore cannot be cut (possible corner case not implemented yet)");
        }

        // Compute intersection points these are required for any cut cell part (interface, interior, exterior)
        // get the number of intersection points
        std::vector<double> intersection_points; 
        std::unordered_map<int,int> vertex_case_map;       
        compute_intersection_points(vertex_coordinates, gdim, ls_values, flag_interior, intersection_points, vertex_case_map);

        cut_cell.resize(cut_type_str.size());

        for(int i=0;i<cut_type_str.size();i++)
        {
            create_cut_cell(vertex_coordinates, gdim, ls_values, cut_type_str[i], cut_cell[i], 
                            triangulate, intersection_points, vertex_case_map);
        }
    };
}
}