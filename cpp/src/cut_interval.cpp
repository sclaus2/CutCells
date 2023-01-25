// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#include "cut_interval.h"
#include "cell_flags.h"
#include <unordered_map>

namespace cutcells::cell
{
    namespace{
        int interval_sub_element[4][2] = 
        {   {-1 -1},           // 0
            {100,0},           // 1
            {0, 101},          // 2
            {100,101}          // 3
        };
    }

    namespace interval{

        int get_num_intersection_points(const int &flag)
        {
            // number of intersection points for intersected triangle is always 2
            return 1; 
        }

        // cut interval by linear interpolation to obtain intersection point
        void compute_intersection_point(const double& level, const std::span<const double> p0, const std::span<const double> p1,
                 const double& v0, const double& v1, std::vector<double>& intersection_point, const int & offset)
        {
            // TODO: Catch precision errors and almost alignments.
            for(int i=0;i<p0.size();i++)
            {
                intersection_point[i+offset] = p0[i] +  (p1[i] - p0[i])*(level - v0)/(v1-v0);
            }
        };

        void create_interface_cells(std::vector<std::vector<int>>& interface_cells, std::vector<type>& interface_cell_types)
        {
            interface_cells.resize(1);
            interface_cells[0].resize(1);
            interface_cells[0][0] = 0;

            interface_cell_types.resize(1);
            interface_cell_types[0] = type::point;
        }

        void create_sub_cell_vertex_coords(const int& flag, const std::span<const double> vertex_coordinates, const int gdim, 
                                        const std::span<const double> intersection_points, 
                                        std::vector<double>& coords, std::unordered_map<int,int>& vertex_case_map)
        {
            int num_intersection_points = intersection_points.size()/gdim;
            type cell_type = cell::type::interval;
            int total_num_points = get_num_vertices(cell_type);

            // Copy vertex coordinates into cut_cell
            coords.resize(total_num_points*gdim);

            for(int i=0;i<intersection_points.size();i++)
            {
                coords[i] = intersection_points[i];
            }

            vertex_case_map[0] = 0;

            int ncols = 2;
            int offset = num_intersection_points; 
            int cnt = num_intersection_points; 
            int k = 0;

            for(int i=0;i<ncols;i++)
            {
            if( interval_sub_element[flag][i] < 100 || interval_sub_element[flag][i]==-1)
            {
                //intersection point which is alreay added or is not a vertex
            }
            else 
            {    
                    int vertex_id = interval_sub_element[flag][i];
                    vertex_case_map[interval_sub_element[flag][i]] = cnt;
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

        void create_sub_cells(const int& flag, std::vector<std::vector<int>>& sub_cells, 
                          std::vector<type>& sub_cell_types, std::unordered_map<int,int>& vertex_case_map)
        {
            // Allocate memory
            int num_sub_elements = 1;
            sub_cells.resize(num_sub_elements);
            sub_cell_types.resize(num_sub_elements);

            for(int i=0;i<num_sub_elements;i++)
            {
                type sub_cell_type = type::interval;
                sub_cell_types[i] = sub_cell_type;
                int num_vertices = get_num_vertices(sub_cell_type);
                sub_cells[i].resize(num_vertices);

                for(int j=0;j<num_vertices;j++)
                {
                    sub_cells[i][j] = vertex_case_map[interval_sub_element[flag][j]];
                }
            }
        }

        void create_cut_cell(const std::span<const double> vertex_coordinates, 
                        const int gdim,  const std::span<const double> ls_values,
                        const std::string& cut_type_str,
                        CutCell& cut_cell, 
                        const std::span<const double> intersection_points)
        {
            cut_cell._gdim = gdim;
            std::unordered_map<int,int> vertex_case_map;
            
            if(cut_type_str=="phi=0")
            {
                cut_cell._tdim = 0;
                cut_cell._vertex_coords.resize(intersection_points.size());
                for(int i=0;i<intersection_points.size();i++)
                {
                    cut_cell._vertex_coords[i] = intersection_points[i];
                }
                //Fill in cut cell with intersection point for interval this is just the intersection point
                create_interface_cells(cut_cell._connectivity, cut_cell._types);
            }
            else if(cut_type_str=="phi<0")
            {
                cut_cell._tdim = 1;
                int flag_interior = get_entity_flag(ls_values, false);
                create_sub_cell_vertex_coords(flag_interior, vertex_coordinates, gdim, intersection_points, 
                            cut_cell._vertex_coords, vertex_case_map);
                //Determine interior sub-cells
                create_sub_cells(flag_interior, cut_cell._connectivity, 
                            cut_cell._types, vertex_case_map);
            }
            else if(cut_type_str=="phi>0")
            {
                cut_cell._tdim = 1;
                //Determine exterior sub-cells
                int flag_exterior = get_entity_flag(ls_values, true);
                create_sub_cell_vertex_coords(flag_exterior, vertex_coordinates, gdim, intersection_points, 
                            cut_cell._vertex_coords, vertex_case_map);
                create_sub_cells(flag_exterior, cut_cell._connectivity, 
                        cut_cell._types, vertex_case_map);
            }
            else
            {
                throw std::invalid_argument("cutting type unknown");
            }
        }

        // cut interval
        void cut(const std::span<const double> vertex_coordinates, const int gdim, 
                const std::span<const double> ls_values, const std::string& cut_type_str,
                CutCell& cut_cell)
        {
            int flag_interior = get_entity_flag(ls_values, false);

            //throw error if cell is not intersected, only intersected cells should land here 
            if(flag_interior ==0 || flag_interior == 3)
            {
                throw std::invalid_argument("interval is not intersected and therefore cannot be cut");
            }

            // Compute intersection point these are required for any cut cell part (interface, interior, exterior)
            std::vector<double> intersection_point(gdim); 
            std::vector<double> p0(gdim); 
            std::vector<double> p1(gdim);
            double level = 0.0;

            for(auto i=0;i<gdim;i++)
                p0[i] = vertex_coordinates[i];
            
            for(auto i=0;i<gdim;i++)
                p1[i] = vertex_coordinates[gdim+i];

            compute_intersection_point(level, p0, p1,ls_values[0],ls_values[1], intersection_point);

            //Create the cut cell depending on which cut is requested
            create_cut_cell(vertex_coordinates, gdim, ls_values, cut_type_str, cut_cell, intersection_point);
        };
    }
}
