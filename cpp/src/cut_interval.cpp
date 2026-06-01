// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#include "cut_interval.h"
#include "cell_flags.h"
#include "span_math.h"
#include "utils.h"

namespace cutcells::cell
{
    using VertexCaseMap = std::array<int, cutcells::utils::MAX_TOKEN_LOOKUP>;

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

template <std::floating_point T>
        void create_sub_cell_vertex_coords(const int& flag, const std::span<const T> vertex_coordinates, const int gdim,
                                        const std::span<const T> intersection_points,
                                        std::vector<T>& coords, VertexCaseMap& vertex_case_map)
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

        template <std::floating_point T>
        void create_cut_cell(const std::span<const T> vertex_coordinates,
                        const int gdim,  const std::span<const T> ls_values,
                        const std::string& cut_type_str,
                        CutCell<T>& cut_cell,
                        const std::span<const T> intersection_points)
        {
            cut_cell._gdim = gdim;
            cutcells::cell::clear_cell_topology(cut_cell);
            VertexCaseMap vertex_case_map;
            vertex_case_map.fill(-1);

            if(cut_type_str=="phi=0")
            {
                cut_cell._tdim = 0;
                cut_cell._vertex_coords.resize(intersection_points.size());
                for(int i=0;i<intersection_points.size();i++)
                {
                    cut_cell._vertex_coords[i] = intersection_points[i];
                }
                // Interface is the single intersection point (index 0)
                const int v = 0;
                cutcells::cell::append_cell(cut_cell, type::point, &v, 1);
            }
            else if(cut_type_str=="phi<0")
            {
                cut_cell._tdim = 1;
                int flag_interior = get_entity_flag(ls_values, false);
                create_sub_cell_vertex_coords(flag_interior, vertex_coordinates, gdim, intersection_points,
                            cut_cell._vertex_coords, vertex_case_map);
                // Interval always produces exactly one sub-interval with 2 vertices
                std::array<int, 2> verts;
                verts[0] = vertex_case_map[interval_sub_element[flag_interior][0]];
                verts[1] = vertex_case_map[interval_sub_element[flag_interior][1]];
                cutcells::cell::append_cell(cut_cell, type::interval, verts, 2);
            }
            else if(cut_type_str=="phi>0")
            {
                cut_cell._tdim = 1;
                int flag_exterior = get_entity_flag(ls_values, true);
                create_sub_cell_vertex_coords(flag_exterior, vertex_coordinates, gdim, intersection_points,
                            cut_cell._vertex_coords, vertex_case_map);
                // Interval always produces exactly one sub-interval with 2 vertices
                std::array<int, 2> verts;
                verts[0] = vertex_case_map[interval_sub_element[flag_exterior][0]];
                verts[1] = vertex_case_map[interval_sub_element[flag_exterior][1]];
                cutcells::cell::append_cell(cut_cell, type::interval, verts, 2);
            }
            else
            {
                throw std::invalid_argument("cutting type unknown");
            }

            cutcells::utils::create_vertex_parent_entity_map<T>(vertex_case_map, cut_cell._vertex_parent_entity);
        }

        // cut interval
        template <std::floating_point T>
        void cut(const std::span<const T> vertex_coordinates, const int gdim,
                const std::span<const T> ls_values, const std::string& cut_type_str,
                CutCell<T>& cut_cell)
        {
            int flag_interior = get_entity_flag(ls_values, false);

            //throw error if cell is not intersected, only intersected cells should land here 
            if(flag_interior ==0 || flag_interior == 3)
            {
                throw std::invalid_argument("interval is not intersected and therefore cannot be cut");
            }

            // Compute intersection point these are required for any cut cell part (interface, interior, exterior)
            thread_local std::vector<T> intersection_point;
            intersection_point.resize(gdim);
            std::array<T, 3> p0 = {};
            std::array<T, 3> p1 = {};
            T level = 0.0;

            for(auto i=0;i<gdim;i++)
                p0[i] = vertex_coordinates[i];

            for(auto i=0;i<gdim;i++)
                p1[i] = vertex_coordinates[gdim+i];

            compute_intersection_point<T>(level, std::span<const T>(p0.data(), gdim),
                                          std::span<const T>(p1.data(), gdim),
                                          ls_values[0], ls_values[1], intersection_point);

            //Create the cut cell depending on which cut is requested
            create_cut_cell<T>(vertex_coordinates, gdim, ls_values, cut_type_str, cut_cell, intersection_point);
        };

      template <std::floating_point T>
      T volume(const std::span<const T> vertex_coordinates, const int gdim)
      {
          const auto p0 = vertex_coordinates.subspan(0, gdim);
          const auto p1 = vertex_coordinates.subspan(gdim, gdim);

          T length =  cutcells::math::distance(p0, p1);
          return length;
      }

//-----------------------------------------------------------------------------
    template void cut(const std::span<const double> vertex_coordinates, const int gdim,
            const std::span<const double> ls_values, const std::string& cut_type_str,
            CutCell<double>& cut_cell);
    template void cut(const std::span<const float> vertex_coordinates, const int gdim,
              const std::span<const float> ls_values, const std::string& cut_type_str,
              CutCell<float>& cut_cell);

    template double volume(const std::span<const double> vertex_coordinates, const int gdim);
    template float volume(const std::span<const float> vertex_coordinates, const int gdim);

//-----------------------------------------------------------------------------
    }
}
