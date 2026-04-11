// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include "cut_tetrahedron.h"
#include "cut_interval.h"
#include "cell_flags.h"
#include "triangulation.h"
#include "span_math.h"
#include "utils.h"

#include <cassert>
#include <stdexcept>
#include <iostream>
#include <array>
#include <cmath>

namespace cutcells::cell
{

// Look up tables for intersection
namespace{
    using VertexCaseMap = std::array<int, cutcells::utils::MAX_TOKEN_LOOKUP>;
    constexpr int reserve_vertex_coords = 16;
    constexpr int reserve_connectivity = 16;
    constexpr int reserve_types = 16;

    // Choose numbering of edges for cutting in accordance with numbering in vtk
    // Tetrahedron numbering:
    //
    //      2
    //     /|\
    //    / | \
    //   / /3\ \
    //   |/___\|
    //   0     1
    //
    constexpr int edges[6][2] = 
    {
        { 0, 1 }, // 0
        { 1, 2 }, // 1
        { 2, 0 }, // 2
        { 0, 3 }, // 3
        { 1, 3 }, // 4
        { 2, 3 }, // 5
    };

    // List of ids of intersected edges these will be used below to obtain 
    // the intersection points (one intersection point is assumed per edge) 
    int tetrahedron_intersected_edges[16][4] = 
    {   {-1, -1, -1, -1},                    // 0
        {3, 0, 2, -1},                       // 1 
        {1, 0, 4, -1},                       // 2
        {2, 3, 4, 1},                        // 3 
        {2, 1, 5, -1},                       // 4 
        {5, 3, 0, 1},                        // 5
        {2, 0, 4, 5},                        // 6 
        {5, 3, 4, -1},                       // 7 
        {4, 3, 5, -1},                       // 8
        {4, 0, 2, 5},                        // 9 
        {5, 3, 0, 1},                        // 10 
        {2, 5, 1, -1},                       // 11
        {4, 3, 2, 1},                        // 12
        {4, 0, 1, -1},                       // 13 
        {2, 0, 3, -1},                       // 14
        {-1, -1, -1, -1},                    // 15
    };

    // List of cell types for each interface produced in each cut case
    type interface_sub_element_cell_types[16] = 
    {   type::quadrilateral,                  // 0
        type::triangle,                       // 1
        type::triangle,                       // 2
        type::quadrilateral,                  // 3
        type::triangle,                       // 4
        type::quadrilateral,                  // 5
        type::quadrilateral,                  // 6
        type::triangle,                       // 7
        type::triangle,                       // 8
        type::quadrilateral,                  // 9
        type::quadrilateral,                  // 10
        type::triangle,                       // 11
        type::quadrilateral,                  // 12
        type::triangle,                       // 13
        type::triangle,                       // 14
        type::quadrilateral,                  // 15
    };

    // List of vertex ids that form sub elements created by intersection. 
    // Intersection points are numbered with 0,1,2...
    // original vertices are numbered by 100, 101, ...
    //
    //      102
    //     / | \
    //    /  |  \
            103
    //   / /   \ \
    //  / /_____\ \
    //   100   101
    //
    // The list gives back either the tetrahedron or prism produced by the intersection
    // This list can be used to reconstruct the coordinates of the sub-elements.
    int tetrahedron_sub_element[16][6] = 
    {
        {  -1, -1, -1, -1, -1, -1  },    // 0
        {  0, 3, 2, 100, -1, -1  },      // 1
        {  0, 1, 4, 101, -1, -1  },      // 2
        {  100, 2, 3, 101, 1, 4  },      // 3  (prism: swapped bottom/top)
        {  1, 2, 5, 102, -1, -1  },      // 4
        {  100, 3, 0, 102, 5, 1  },      // 5  (prism: swapped bottom/top)
        {  101, 0, 4, 102, 2, 5  },      // 6  (prism: swapped bottom/top)
        {  100, 101, 102, 3, 4, 5},      // 7  (prism: swapped bottom/top)
        {  3, 4, 5, 103, -1, -1  },      // 8
        {  100, 0, 2, 103, 4, 5  },      // 9  (prism: swapped bottom/top)
        {  101, 1, 0, 103, 5, 3  },      // 10 (prism: swapped bottom/top)
        {  2, 1, 5, 100, 101, 103},      // 11 (prism: swapped bottom/top)
        {  3, 103, 4, 2, 102, 1  },      // 12 (prism: swapped bottom/top)
        {  100, 102, 103, 4, 0, 1},      // 13 (prism: swapped bottom/top)
        {  101, 103, 102, 0, 3, 2},      // 14 (prism: swapped bottom/top)
        {  100, 101, 102, 103, -1, -1}   // 15
    };

    // List of cell types for each sub-element produced in each cut case
    type tetrahedron_sub_element_cell_types[16] = 
    {   type::tetrahedron,                    // 0
        type::tetrahedron,                    // 1
        type::tetrahedron,                    // 2
        type::prism,                          // 3
        type::tetrahedron,                    // 4
        type::prism,                          // 5
        type::prism,                          // 6
        type::prism,                          // 7
        type::tetrahedron,                    // 8
        type::prism,                          // 9
        type::prism,                          // 10
        type::prism,                          // 11
        type::prism,                          // 12
        type::prism,                          // 13
        type::prism,                          // 14
        type::tetrahedron,                    // 15
    };
}

namespace tetrahedron{

    int get_num_intersection_points(const int &flag)
    {
        // number of intersection points depends if intersection is triangle or quadrilateral
        return cell::get_num_vertices(interface_sub_element_cell_types[flag]);
    }

    // for triangles the number of sub-elements without triangulation is 1
    int get_num_sub_elements(const int &flag, bool triangulate)
    {
        type cell_type = tetrahedron_sub_element_cell_types[flag];
        int num_sub_elements = 1;

        if(cell_type == type::prism && triangulate == true)
        {
            num_sub_elements = 3;
        }
        return num_sub_elements;
    }

    int get_num_interface_elements(const int &flag, bool triangulate)
    {
        type cell_type = interface_sub_element_cell_types[flag];
        int num_sub_elements = 1;

        if(cell_type == type::quadrilateral && triangulate == true)
        {
            num_sub_elements = 2;
        }
        return num_sub_elements;
    }

    template <std::floating_point T>
    void compute_intersection_points(std::span<const T> vertex_coordinates, const int gdim,
             std::span<const T> ls_values, const int flag, std::vector<T>& intersection_points,
             VertexCaseMap& vertex_case_map)
    {
        // vertex ids will be edges[edge_0][0] edges[edge_0][1]
        // get vertex coordinates and interpolate
        int num_intersection_points = get_num_intersection_points(flag);

        intersection_points.resize(num_intersection_points*gdim);

        // gdim is at most 3; use stack arrays to avoid heap allocation per intersection point
        std::array<T, 3> v0 = {};
        std::array<T, 3> v1 = {};

        for(int ip=0; ip<num_intersection_points; ip++)
        {
            // edge has two vertices
            // intersection points are listed first in triangle case table
            // therefore the first two entries are the intersected edges
            int vertex_id_0 = edges[tetrahedron_intersected_edges[flag][ip]][0];
            int vertex_id_1 = edges[tetrahedron_intersected_edges[flag][ip]][1];

            for(int j=0;j<gdim;j++)
            {
                v0[j] = vertex_coordinates[vertex_id_0*gdim+j];
                v1[j] = vertex_coordinates[vertex_id_1*gdim+j];
            }

            T ls0 = ls_values[vertex_id_0];
            T ls1 = ls_values[vertex_id_1];

            interval::compute_intersection_point<T>(T(0), std::span<const T>(v0.data(), gdim),
                                                    std::span<const T>(v1.data(), gdim),
                                                    ls0, ls1, intersection_points, ip * gdim);
        }

        int cnt = 0;
        for(int ip=0; ip<num_intersection_points; ip++)
        {
            vertex_case_map[tetrahedron_intersected_edges[flag][ip]] = cnt;
            cnt++;
        }

    }

    template <std::floating_point T>
    void create_sub_cell_vertex_coords(const int& flag, const std::span<const T> vertex_coordinates, const int gdim,
                                       const std::span<const T> intersection_points,
                                       std::vector<T>& coords, VertexCaseMap& vertex_case_map)
    {
        int num_intersection_points = intersection_points.size()/gdim;
        type cell_type = tetrahedron_sub_element_cell_types[flag];
        int total_num_points = get_num_vertices(cell_type);

        // Copy vertex coordinates into cut_cell
        coords.resize(total_num_points*gdim);

        for(int i=0;i<intersection_points.size();i++)
        {
            coords[i] = intersection_points[i];
        }

        int ncols = 6;
        int offset = num_intersection_points;
        int cnt = num_intersection_points;
        int k = 0;

        for(int i=0;i<ncols;i++)
        {
           if( tetrahedron_sub_element[flag][i] < 100 || tetrahedron_sub_element[flag][i]==-1)
           {
            //intersection point which is alreay added or is not a vertex
           }
           else
           {
                int vertex_id = tetrahedron_sub_element[flag][i];
                vertex_case_map[vertex_id] = cnt;
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
                         CutCell<T>& cut_cell, bool triangulate,
                         const std::span<const T> intersection_points,
                         VertexCaseMap& vertex_case_map)
    {
        cut_cell._gdim = gdim;
        cutcells::cell::clear_cell_topology(cut_cell);
        const cut_type cut_kind = string_to_cut_type(cut_type_str);

        if(cut_kind == cut_type::phieq0)
        {
            cut_cell._tdim = 2;
            int flag_interior = get_entity_flag(ls_values, false);
            // Copy vertex coordinates into cut_cell
            int num_intersection_points = get_num_intersection_points(flag_interior);
            cut_cell._vertex_coords.resize(num_intersection_points*gdim);
            for(int i=0;i<intersection_points.size();i++)
            {
                cut_cell._vertex_coords[i] = intersection_points[i];
            }
            // Append interface cells directly
            type iface_type = interface_sub_element_cell_types[flag_interior];
            if(iface_type == type::quadrilateral && triangulate)
            {
                std::vector<std::vector<int>> triangles;
                triangulation(iface_type, tetrahedron_intersected_edges[flag_interior], triangles);
                for(int i = 0; i < static_cast<int>(triangles.size()); ++i)
                {
                    std::array<int, 3> t;
                    t[0] = vertex_case_map[triangles[i][0]];
                    t[1] = vertex_case_map[triangles[i][1]];
                    t[2] = vertex_case_map[triangles[i][2]];
                    cutcells::cell::append_cell(cut_cell, type::triangle, t, 3);
                }
            }
            else
            {
                int num_vertices = get_num_vertices(iface_type);
                std::array<int, 4> verts;
                for(int j = 0; j < num_vertices; ++j)
                    verts[j] = vertex_case_map[tetrahedron_intersected_edges[flag_interior][j]];
                cutcells::cell::append_cell(cut_cell, iface_type, verts, num_vertices);
            }
        }
        else if(cut_kind == cut_type::philt0)
        {
            cut_cell._tdim = 3;
            int flag_interior = get_entity_flag(ls_values, false);
            create_sub_cell_vertex_coords(flag_interior, vertex_coordinates, gdim, intersection_points, 
                        cut_cell._vertex_coords, vertex_case_map);
            // Append sub-cells directly
            type sub_cell_type = tetrahedron_sub_element_cell_types[flag_interior];
            if(sub_cell_type == type::prism && triangulate)
            {
                std::vector<std::vector<int>> tets;
                triangulation(sub_cell_type, tetrahedron_sub_element[flag_interior], tets);
                for(int i = 0; i < static_cast<int>(tets.size()); ++i)
                {
                    std::array<int, 4> t;
                    t[0] = vertex_case_map[tets[i][0]];
                    t[1] = vertex_case_map[tets[i][1]];
                    t[2] = vertex_case_map[tets[i][2]];
                    t[3] = vertex_case_map[tets[i][3]];
                    cutcells::cell::append_cell(cut_cell, type::tetrahedron, t, 4);
                }
            }
            else
            {
                int num_vertices = get_num_vertices(sub_cell_type);
                std::array<int, 6> verts;
                for(int j = 0; j < num_vertices; ++j)
                    verts[j] = vertex_case_map[tetrahedron_sub_element[flag_interior][j]];
                cutcells::cell::append_cell(cut_cell, sub_cell_type, verts, num_vertices);
            }
        }
        else if(cut_kind == cut_type::phigt0)
        {
            cut_cell._tdim = 3;
            int flag_exterior = get_entity_flag(ls_values, true);
            create_sub_cell_vertex_coords(flag_exterior, vertex_coordinates, gdim, intersection_points, 
                        cut_cell._vertex_coords, vertex_case_map);
            // Append sub-cells directly
            type sub_cell_type = tetrahedron_sub_element_cell_types[flag_exterior];
            if(sub_cell_type == type::prism && triangulate)
            {
                std::vector<std::vector<int>> tets;
                triangulation(sub_cell_type, tetrahedron_sub_element[flag_exterior], tets);
                for(int i = 0; i < static_cast<int>(tets.size()); ++i)
                {
                    std::array<int, 4> t;
                    t[0] = vertex_case_map[tets[i][0]];
                    t[1] = vertex_case_map[tets[i][1]];
                    t[2] = vertex_case_map[tets[i][2]];
                    t[3] = vertex_case_map[tets[i][3]];
                    cutcells::cell::append_cell(cut_cell, type::tetrahedron, t, 4);
                }
            }
            else
            {
                int num_vertices = get_num_vertices(sub_cell_type);
                std::array<int, 6> verts;
                for(int j = 0; j < num_vertices; ++j)
                    verts[j] = vertex_case_map[tetrahedron_sub_element[flag_exterior][j]];
                cutcells::cell::append_cell(cut_cell, sub_cell_type, verts, num_vertices);
            }
        }
        else
        {
            throw std::invalid_argument("cutting type unknown");
        }
    }

    template <std::floating_point T>
    void str(CutCell<T>& cut_cell, const VertexCaseMap& vertex_case_map)
    {
        std::cout << "vertex case map=[";
        for(std::size_t token = 0; token < vertex_case_map.size(); ++token)
        {
            if (vertex_case_map[token] >= 0)
                std::cout << token << ": " << vertex_case_map[token] << std::endl;
        }
        std::cout << "]" << std::endl;

        std::cout << "connectivity=[";
        for(auto &i: cut_cell._connectivity)
        {
            std::cout << i << ", ";
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

    // cut tetrahedron
    template <std::floating_point T>
    void cut(const std::span<const T> vertex_coordinates, const int gdim,
             const std::span<const T> ls_values, const std::string& cut_type_str,
             CutCell<T>& cut_cell, bool triangulate)
    {
        int flag_interior = get_entity_flag(ls_values, false);

        // throw error if cell is not intersected, only intersected cells should land here
        if(flag_interior ==0 || flag_interior == 15)
        {
            throw std::invalid_argument("tetrahedron is not intersected and therefore cannot be cut");
        }

        // Compute intersection points these are required for any cut cell part (interface, interior, exterior)
        // get the number of intersection points
        thread_local std::vector<T> intersection_points;
        intersection_points.clear();
        intersection_points.reserve(4 * gdim);
        // the vertex case map,
        // first few entries map from intersected edge to intersection point number
        // next entries map from orginal vertex id to number of vertex in vertex_coordinates of CutCell (renumbered to go from 0,...,N)
        // example: intersected edges 0 -> 0, 2 -> 1
        //          then orginal vertex 101 -> 2 , 102 -> 3 etc.
        VertexCaseMap vertex_case_map;
        vertex_case_map.fill(-1);
        compute_intersection_points<T>(vertex_coordinates, gdim, ls_values, flag_interior, intersection_points, vertex_case_map);

        cut_cell._vertex_coords.reserve(reserve_vertex_coords * gdim);
        cutcells::cell::reserve_cell_topology(cut_cell, reserve_connectivity, reserve_types);

        create_cut_cell<T>(vertex_coordinates, gdim, ls_values, cut_type_str, cut_cell,
                        triangulate, intersection_points, vertex_case_map);

        cutcells::utils::create_vertex_parent_entity_map<T>(vertex_case_map, cut_cell._vertex_parent_entity, 6, 4);
    }

    template <std::floating_point T>
    T volume(const std::span<const T> vertex_coordinates, const int gdim)
    {
            const auto a = cutcells::math::to_vec3<T>(vertex_coordinates.subspan(0, gdim));
            const auto b = cutcells::math::to_vec3<T>(vertex_coordinates.subspan(gdim, gdim));
            const auto c = cutcells::math::to_vec3<T>(vertex_coordinates.subspan(2*gdim, gdim));
            const auto d = cutcells::math::to_vec3<T>(vertex_coordinates.subspan(3*gdim, gdim));

      const auto ad = cutcells::math::subtract(a,d);
      const auto bd = cutcells::math::subtract(b,d);
      const auto cd = cutcells::math::subtract(c,d);

      const auto bdxcd = cutcells::math::cross<T>(bd,cd);

      T vol = fabs(cutcells::math::dot<T>(ad,bdxcd))/6.0;
      return vol;
    }

    //-----------------------------------------------------------------------------
    template void cut(const std::span<const double> vertex_coordinates, const int gdim,
            const std::span<const double> ls_values, const std::string& cut_type_str,
            CutCell<double>& cut_cell, bool triangulate);
    template void cut(const std::span<const float> vertex_coordinates, const int gdim,
              const std::span<const float> ls_values, const std::string& cut_type_str,
              CutCell<float>& cut_cell, bool triangulate);

    template void compute_intersection_points(std::span<const double> vertex_coordinates, const int gdim,
             std::span<const double> ls_values, const int flag, std::vector<double>& intersection_points,
             VertexCaseMap& vertex_case_map);
    template void compute_intersection_points(std::span<const float> vertex_coordinates, const int gdim,
             std::span<const float> ls_values, const int flag, std::vector<float>& intersection_points,
             VertexCaseMap& vertex_case_map);

    template double volume(const std::span<const double> vertex_coordinates, const int gdim);
    template float volume(const std::span<const float> vertex_coordinates, const int gdim);

    //-----------------------------------------------------------------------------
}
}
