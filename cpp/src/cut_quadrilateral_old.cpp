// Copyright (c) 2022-2025 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include "cut_quadrilateral.h"
#include "cell_flags.h"
#include "cell_topology.h"
#include "cut_interval.h"
#include "cut_triangle.h"
#include "generated/cut_quadrilateral_inside_tables.h"
#include "generated/cut_quadrilateral_outside_tables.h"
#include "generated/cut_quadrilateral_interface_tables.h"
#include "triangulation.h"
#include "utils.h"

#include <array>
#include <concepts>
#include <stdexcept>
#include <unordered_map>

namespace cutcells::cell::quadrilateral
{
    using namespace cutcells::cell::generated;

    namespace
    {
        template <std::floating_point T>
        void compute_intersection_points(const std::span<const T> vertex_coordinates, const int gdim,
                                          const std::span<const T> ls_values, const int flag,
                                          std::vector<T>& intersection_points,
                                          std::unordered_map<int, int>& vertex_case_map)
        {
            // Per-case edge flags come from the generated table; use them to locate linear edge/edge intersections.
            const auto& intersected = cut_quadrilateral_inside_intersected_edges[flag];

            int num_intersection_points = 0;
            for (int e = 0; e < 4; ++e)
            {
                if (intersected[e] == 1)
                    ++num_intersection_points;
            }

            intersection_points.resize(num_intersection_points * gdim);

            std::vector<T> v0(gdim);
            std::vector<T> v1(gdim);
            std::vector<T> ip(gdim);

            int ip_index = 0;
            for (int e = 0; e < 4; ++e)
            {
                if (intersected[e] == 0)
                    continue;

                const int v0_id = quadrilateral_edges[e][0];
                const int v1_id = quadrilateral_edges[e][1];

                for (int j = 0; j < gdim; ++j)
                {
                    v0[j] = vertex_coordinates[v0_id * gdim + j];
                    v1[j] = vertex_coordinates[v1_id * gdim + j];
                }

                const T ls0 = ls_values[v0_id];
                const T ls1 = ls_values[v1_id];

                interval::compute_intersection_point<T>(0.0, v0, v1, ls0, ls1, ip);

                for (int j = 0; j < gdim; ++j)
                {
                    intersection_points[ip_index * gdim + j] = ip[j];
                }

                vertex_case_map[e] = ip_index;
                ++ip_index;
            }
        }

        // sub_element stores a flat stream: [n, v0, ..., v{n-1}] repeated per emitted cell.
        // sub_element_offset tells us which cell stream entries belong to a case; this helper finds
        // the starting token index for a given emitted cell id.
        inline int cell_token_start(const int target_cell_idx, const int* sub_element, const int total_cells)
        {
            int token_idx = 0;
            for (int c = 0; c < target_cell_idx && c < total_cells; ++c)
            {
                const int nverts = sub_element[token_idx];
                token_idx += 1 + nverts;
            }
            return token_idx;
        }

        template <std::floating_point T>
        void append_vertex_if_needed(int token, const std::span<const T> vertex_coordinates, int gdim,
                                     std::unordered_map<int, int>& vertex_case_map,
                                     std::vector<T>& out_coords)
        {
            if (vertex_case_map.find(token) != vertex_case_map.end())
                return;

            // token represents original vertex 100+vid; store it once and reuse via vertex_case_map
            const int vid = token - 100;
            const int local_idx = static_cast<int>(out_coords.size() / gdim);
            vertex_case_map[token] = local_idx;
            out_coords.resize(out_coords.size() + gdim);
            for (int j = 0; j < gdim; ++j)
            {
                out_coords[local_idx * gdim + j] = vertex_coordinates[vid * gdim + j];
            }
        }

        template <std::floating_point T>
        void decode_case(const std::span<const T> vertex_coordinates, const int gdim, int flag,
                         const int* sub_element_offset, const type* sub_element_cell_types, const int* sub_element,
                         bool triangulate, CutCell<T>& cut_cell, std::unordered_map<int, int>& vertex_case_map)
        {
            const int cell_begin = sub_element_offset[flag];
            const int cell_end = sub_element_offset[flag + 1];
            // sub_element_offset has length NCases+1 (16+1 here); the last entry stores the total emitted cells.
            const int total_cells = sub_element_offset[16];

            // Ensure original vertices referenced in this case are appended
            int token_idx = cell_token_start(cell_begin, sub_element, total_cells);
            for (int cell_idx = cell_begin; cell_idx < cell_end; ++cell_idx)
            {
                const int nverts = sub_element[token_idx];
                for (int i = 1; i <= nverts; ++i)
                {
                    const int token = sub_element[token_idx + i];
                    if (token >= 100)
                    {
                        append_vertex_if_needed(token, vertex_coordinates, gdim, vertex_case_map, cut_cell._vertex_coords);
                    }
                }
                token_idx += 1 + nverts;
            }

            // Build connectivity
            token_idx = cell_token_start(cell_begin, sub_element, total_cells);
            for (int cell_idx = cell_begin; cell_idx < cell_end; ++cell_idx)
            {
                const int nverts = sub_element[token_idx];
                const type sub_type = sub_element_cell_types[cell_idx];

                std::vector<int> verts_local(nverts);
                for (int i = 0; i < nverts; ++i)
                {
                    const int token = sub_element[token_idx + 1 + i];
                    if (token >= 100)
                    {
                        verts_local[i] = vertex_case_map[token];
                    }
                    else
                    {
                        // edge intersection point
                        verts_local[i] = vertex_case_map[token];
                    }
                }

                if (triangulate && sub_type == type::quadrilateral)
                {
                    int verts_arr[4] = {verts_local[0], verts_local[1], verts_local[2], verts_local[3]};
                    std::vector<std::vector<int>> tris;
                    triangulation(sub_type, verts_arr, tris);
                    for (auto& tri : tris)
                    {
                        cut_cell._types.push_back(type::triangle);
                        cut_cell._connectivity.push_back(tri);
                    }
                }
                else
                {
                    cut_cell._types.push_back(sub_type);
                    cut_cell._connectivity.push_back(std::move(verts_local));
                }

                token_idx += 1 + nverts;
            }
        }
    } // namespace

    template <std::floating_point T>
    void cut(const std::span<const T> vertex_coordinates, const int gdim,
             const std::span<const T> ls_values, const std::string& cut_type_str,
             CutCell<T>& cut_cell, bool triangulate)
    {
        const int flag_interior = get_entity_flag(ls_values, false);

        if (flag_interior == 0 || flag_interior == 15)
        {
            throw std::invalid_argument("quadrilateral is not intersected and therefore cannot be cut");
        }

        // Compute intersections (shared for all parts)
        std::vector<T> intersection_points;
        std::unordered_map<int, int> vertex_case_map;
        compute_intersection_points(vertex_coordinates, gdim, ls_values, flag_interior, intersection_points, vertex_case_map);

        cut_cell._gdim = gdim;
        cut_cell._vertex_coords = intersection_points;
        cut_cell._types.clear();
        cut_cell._connectivity.clear();

        // Ambiguous opposite-corner cases (masks 0b0101, 0b1010) are already disambiguated in the generated tables
        // into two disjoint pieces (two triangles for volume, two segments for interface).

        if (cut_type_str == "phi=0")
        {
            cut_cell._tdim = 1;
            const int* offsets = cut_quadrilateral_interface_sub_element_offset;
            const type* cell_types = cut_quadrilateral_interface_sub_element_cell_types;
            const int* sub_el = cut_quadrilateral_interface_sub_element;
            decode_case(vertex_coordinates, gdim, flag_interior, offsets, cell_types, sub_el, triangulate, cut_cell, vertex_case_map);
        }
        else if (cut_type_str == "phi<0")
        {
            cut_cell._tdim = 2;
            const int* offsets = cut_quadrilateral_inside_sub_element_offset;
            const type* cell_types = cut_quadrilateral_inside_sub_element_cell_types;
            const int* sub_el = cut_quadrilateral_inside_sub_element;
            decode_case(vertex_coordinates, gdim, flag_interior, offsets, cell_types, sub_el, triangulate, cut_cell, vertex_case_map);
        }
        else if (cut_type_str == "phi>0")
        {
            cut_cell._tdim = 2;
            const int flag_exterior = get_entity_flag(ls_values, true);
            const int* offsets = cut_quadrilateral_outside_sub_element_offset;
            const type* cell_types = cut_quadrilateral_outside_sub_element_cell_types;
            const int* sub_el = cut_quadrilateral_outside_sub_element;
            decode_case(vertex_coordinates, gdim, flag_exterior, offsets, cell_types, sub_el, triangulate, cut_cell, vertex_case_map);
        }
        else
        {
            throw std::invalid_argument("cutting type unknown");
        }

        cutcells::utils::create_vertex_parent_entity_map<T>(vertex_case_map, cut_cell._vertex_parent_entity);
    }

    template <std::floating_point T>
    void cut(const std::span<const T> vertex_coordinates, const int gdim,
             const std::span<const T> ls_values, const std::vector<std::string>& cut_type_str,
             std::vector<CutCell<T>>& cut_cell, bool triangulate)
    {
        const int flag_interior = get_entity_flag(ls_values, false);
        if (flag_interior == 0 || flag_interior == 15)
        {
            throw std::invalid_argument("quadrilateral is not intersected and therefore cannot be cut");
        }

        cut_cell.resize(cut_type_str.size());
        for (std::size_t i = 0; i < cut_type_str.size(); ++i)
        {
            cut(vertex_coordinates, gdim, ls_values, cut_type_str[i], cut_cell[i], triangulate);
        }
    }

    template <std::floating_point T>
    T volume(const std::span<const T> vertex_coordinates, const int gdim)
    {
        // Split quad into two triangles (0-1-2) and (0-2-3)
        std::vector<T> t0;
        t0.reserve(3 * gdim);
        for (int j = 0; j < gdim; ++j)
        {
            t0.push_back(vertex_coordinates[0 * gdim + j]);
        }
        for (int j = 0; j < gdim; ++j)
        {
            t0.push_back(vertex_coordinates[1 * gdim + j]);
        }
        for (int j = 0; j < gdim; ++j)
        {
            t0.push_back(vertex_coordinates[2 * gdim + j]);
        }

        std::vector<T> t1;
        t1.reserve(3 * gdim);
        for (int j = 0; j < gdim; ++j)
        {
            t1.push_back(vertex_coordinates[0 * gdim + j]);
        }
        for (int j = 0; j < gdim; ++j)
        {
            t1.push_back(vertex_coordinates[2 * gdim + j]);
        }
        for (int j = 0; j < gdim; ++j)
        {
            t1.push_back(vertex_coordinates[3 * gdim + j]);
        }

        return triangle::volume<T>(t0, gdim) + triangle::volume<T>(t1, gdim);
    }

    template void cut(const std::span<const double> vertex_coordinates, const int gdim,
                      const std::span<const double> ls_values, const std::string& cut_type_str,
                      CutCell<double>& cut_cell, bool triangulate);
    template void cut(const std::span<const float> vertex_coordinates, const int gdim,
                      const std::span<const float> ls_values, const std::string& cut_type_str,
                      CutCell<float>& cut_cell, bool triangulate);

    template void cut(const std::span<const double> vertex_coordinates, const int gdim,
                      const std::span<const double> ls_values, const std::vector<std::string>& cut_type_str,
                      std::vector<CutCell<double>>& cut_cell, bool triangulate);
    template void cut(const std::span<const float> vertex_coordinates, const int gdim,
                      const std::span<const float> ls_values, const std::vector<std::string>& cut_type_str,
                      std::vector<CutCell<float>>& cut_cell, bool triangulate);

    template double volume(const std::span<const double> vertex_coordinates, const int gdim);
    template float volume(const std::span<const float> vertex_coordinates, const int gdim);
}
