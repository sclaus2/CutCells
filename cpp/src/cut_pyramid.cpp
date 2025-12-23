// Copyright (c) 2025 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include "cut_pyramid.h"

#include "cell_flags.h"
#include "cut_cell.h"
#include "cut_interval.h"
#include "generated/cut_pyramid_inside_tables.h"
#include "generated/cut_pyramid_interface_tables.h"
#include "generated/cut_pyramid_outside_tables.h"
#include "utils.h"

#include <concepts>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace cutcells::cell::pyramid
{
    namespace
    {
        // VTK_PYRAMID vertex ordering assumed:
        // base quad: 0,1,2,3 and apex: 4.
        // Edge ids must match the VTK TableBasedClip case stream.
        constexpr int edges[8][2] = {
            {0, 1}, {1, 2}, {2, 3}, {3, 0}, // base
            {0, 4}, {1, 4}, {2, 4}, {3, 4}  // sides
        };

        [[noreturn]] void throw_missing_token(const int flag, const int token, const char* where)
        {
            throw std::runtime_error(std::string(where) + ": missing token=" + std::to_string(token)
                                     + " for flag=" + std::to_string(flag));
        }

        int lookup_token_or_throw(const std::unordered_map<int, int>& vertex_case_map,
                                  const int flag,
                                  const int token,
                                  const char* where)
        {
            const auto it = vertex_case_map.find(token);
            if (it == vertex_case_map.end())
                throw_missing_token(flag, token, where);
            return it->second;
        }

        template <std::floating_point T>
        void compute_intersection_points(const std::span<const T> vertex_coordinates, const int gdim,
                                          const std::span<const T> ls_values, const int flag,
                                          std::vector<T>& intersection_points,
                                          std::unordered_map<int, int>& vertex_case_map)
        {
            intersection_points.clear();
            vertex_case_map.clear();

            std::vector<T> v0(gdim);
            std::vector<T> v1(gdim);
            std::vector<T> ip(gdim);

            int ip_index = 0;
            for (int e = 0; e < 8; ++e)
            {
                if (intersected_edges[flag][e] == 0)
                    continue;

                const int v0_id = edges[e][0];
                const int v1_id = edges[e][1];

                for (int j = 0; j < gdim; ++j)
                {
                    v0[j] = vertex_coordinates[v0_id * gdim + j];
                    v1[j] = vertex_coordinates[v1_id * gdim + j];
                }

                const T ls0 = ls_values[v0_id];
                const T ls1 = ls_values[v1_id];

                interval::compute_intersection_point<T>(0.0, v0, v1, ls0, ls1, ip);

                for (int j = 0; j < gdim; ++j)
                    intersection_points.push_back(ip[j]);

                vertex_case_map[e] = ip_index;
                ++ip_index;
            }
        }

        template <std::floating_point T>
        void ensure_vertex_token(const std::span<const T> vertex_coordinates, const int gdim,
                                 const int token,
                                 CutCell<T>& cut_cell,
                                 std::unordered_map<int, int>& vertex_case_map)
        {
            if (vertex_case_map.find(token) != vertex_case_map.end())
                return;

            const int vid = token - 100;
            const int local_idx = static_cast<int>(cut_cell._vertex_coords.size() / gdim);
            vertex_case_map[token] = local_idx;
            cut_cell._vertex_coords.resize(cut_cell._vertex_coords.size() + gdim);
            for (int j = 0; j < gdim; ++j)
                cut_cell._vertex_coords[local_idx * gdim + j] = vertex_coordinates[vid * gdim + j];
        }

        template <std::floating_point T>
        void ensure_special_token(const std::span<const T> vertex_coordinates, const int gdim,
                                  const int flag, const int token,
                                  const int* special_point_count,
                                  const int* special_point_offset,
                                  const int* special_point_data,
                                  CutCell<T>& cut_cell,
                                  std::unordered_map<int, int>& vertex_case_map)
        {
            if (vertex_case_map.find(token) != vertex_case_map.end())
                return;

            if (special_point_count == nullptr || special_point_offset == nullptr || special_point_data == nullptr)
                throw std::runtime_error("Special points referenced but not available for this table");

            const int sp_id = token - 200;
            const int count = special_point_count[flag];
            if (sp_id < 0 || sp_id >= count)
                throw std::runtime_error("Invalid special point id in table");

            int idx = special_point_offset[flag];
            for (int sid = 0; sid < sp_id; ++sid)
            {
                const int nrefs = special_point_data[idx++];
                idx += nrefs;
            }

            const int nrefs = special_point_data[idx++];
            if (nrefs <= 0)
                throw std::runtime_error("Malformed special point definition");

            std::vector<T> coord(gdim, T(0));
            for (int r = 0; r < nrefs; ++r)
            {
                const int ref = special_point_data[idx++];

                if (ref < 100)
                {
                    const int local = lookup_token_or_throw(vertex_case_map, flag, ref,
                                                           "ensure_special_token(ref<100)");
                    for (int j = 0; j < gdim; ++j)
                        coord[j] += cut_cell._vertex_coords[local * gdim + j];
                }
                else if (ref < 200)
                {
                    const int vid = ref - 100;
                    for (int j = 0; j < gdim; ++j)
                        coord[j] += vertex_coordinates[vid * gdim + j];
                }
                else
                {
                    ensure_special_token<T>(vertex_coordinates, gdim, flag, ref,
                                            special_point_count, special_point_offset, special_point_data,
                                            cut_cell, vertex_case_map);
                    const int local = lookup_token_or_throw(vertex_case_map, flag, ref,
                                                           "ensure_special_token(ref>=200)");
                    for (int j = 0; j < gdim; ++j)
                        coord[j] += cut_cell._vertex_coords[local * gdim + j];
                }
            }

            for (int j = 0; j < gdim; ++j)
                coord[j] /= static_cast<T>(nrefs);

            const int local_idx = static_cast<int>(cut_cell._vertex_coords.size() / gdim);
            vertex_case_map[token] = local_idx;
            cut_cell._vertex_coords.resize(cut_cell._vertex_coords.size() + gdim);
            for (int j = 0; j < gdim; ++j)
                cut_cell._vertex_coords[local_idx * gdim + j] = coord[j];
        }

        template <std::floating_point T, int MaxVerts>
        void decode_case(const std::span<const T> vertex_coordinates, const int gdim,
                         const int flag,
                         const int* case_offsets,
                         const type* cell_types,
                         const int (*subcell_verts)[MaxVerts],
                         const int* special_point_count,
                         const int* special_point_offset,
                         const int* special_point_data,
                         bool triangulate,
                         CutCell<T>& cut_cell,
                         std::unordered_map<int, int>& vertex_case_map)
        {
            const int cell_begin = case_offsets[flag];
            const int cell_end = case_offsets[flag + 1];

            for (int cell_idx = cell_begin; cell_idx < cell_end; ++cell_idx)
            {
                const type sub_type = cell_types[cell_idx];
                const int* verts = subcell_verts[cell_idx];

                std::vector<int> verts_local;
                verts_local.reserve(MaxVerts);
                for (int i = 0; i < MaxVerts; ++i)
                {
                    const int token = verts[i];
                    if (token == -1)
                        break;

                    if (token < 100)
                    {
                        // Edge intersection token: already present
                    }
                    else if (token < 200)
                    {
                        ensure_vertex_token<T>(vertex_coordinates, gdim, token, cut_cell, vertex_case_map);
                    }
                    else
                    {
                        ensure_special_token<T>(vertex_coordinates, gdim, flag, token,
                                                special_point_count, special_point_offset, special_point_data,
                                                cut_cell, vertex_case_map);
                    }

                    verts_local.push_back(lookup_token_or_throw(vertex_case_map, flag, token, "decode_case"));
                }

                if (triangulate && sub_type == type::quadrilateral && verts_local.size() == 4)
                {
                    cut_cell._types.push_back(type::triangle);
                    cut_cell._connectivity.push_back({verts_local[0], verts_local[1], verts_local[2]});
                    cut_cell._types.push_back(type::triangle);
                    cut_cell._connectivity.push_back({verts_local[0], verts_local[2], verts_local[3]});
                }
                else
                {
                    cut_cell._types.push_back(sub_type);
                    cut_cell._connectivity.push_back(std::move(verts_local));
                }
            }
        }
    }

    template <std::floating_point T>
    void cut(const std::span<const T> vertex_coordinates, const int gdim,
             const std::span<const T> ls_values, const std::string& cut_type_str,
             CutCell<T>& cut_cell, bool triangulate)
    {
        if (vertex_coordinates.size() != static_cast<std::size_t>(5 * gdim))
            throw std::invalid_argument("pyramid::cut expects 5 vertices");
        if (ls_values.size() != 5)
            throw std::invalid_argument("pyramid::cut expects 5 level set values");

        // CutCells convention: case mask bit i is set when phi_i < 0.
        const int flag_lt0 = get_entity_flag(ls_values, false);
        if (flag_lt0 == 0 || flag_lt0 == 31)
            throw std::invalid_argument("pyramid is not intersected and therefore cannot be cut");

        // Store parent geometry
        cut_cell._parent_cell_type = type::pyramid;
        cut_cell._parent_vertex_coords.resize(5 * gdim);
        cut_cell._parent_vertex_ids.resize(5);
        for (int i = 0; i < 5; ++i)
        {
            for (int j = 0; j < gdim; ++j)
                cut_cell._parent_vertex_coords[i * gdim + j] = vertex_coordinates[i * gdim + j];
            cut_cell._parent_vertex_ids[i] = i;
        }

        // Compute intersections (shared for all parts)
        std::vector<T> intersection_points;
        std::unordered_map<int, int> vertex_case_map;
        compute_intersection_points<T>(vertex_coordinates, gdim, ls_values, flag_lt0,
                                       intersection_points, vertex_case_map);

        cut_cell._gdim = gdim;
        cut_cell._vertex_coords = std::move(intersection_points);
        cut_cell._types.clear();
        cut_cell._connectivity.clear();

        if (cut_type_str == "phi=0")
        {
            cut_cell._tdim = 2;
            decode_case<T, 4>(vertex_coordinates, gdim, flag_lt0,
                              case_subcell_offset_interface, subcell_type_interface, subcell_verts_interface,
                              nullptr, nullptr, nullptr,
                              triangulate, cut_cell, vertex_case_map);
        }
        else if (cut_type_str == "phi<0")
        {
            cut_cell._tdim = 3;
            decode_case<T, 8>(vertex_coordinates, gdim, flag_lt0,
                              case_subcell_offset_inside, subcell_type_inside, subcell_verts_inside,
                              special_point_count_inside, special_point_offset_inside, special_point_data_inside,
                              triangulate, cut_cell, vertex_case_map);
        }
        else if (cut_type_str == "phi>0")
        {
            cut_cell._tdim = 3;
            decode_case<T, 8>(vertex_coordinates, gdim, flag_lt0,
                              case_subcell_offset_outside, subcell_type_outside, subcell_verts_outside,
                              special_point_count_outside, special_point_offset_outside, special_point_data_outside,
                              triangulate, cut_cell, vertex_case_map);
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
        cut_cell.resize(cut_type_str.size());
        for (std::size_t i = 0; i < cut_type_str.size(); ++i)
            cut<T>(vertex_coordinates, gdim, ls_values, cut_type_str[i], cut_cell[i], triangulate);
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
}
