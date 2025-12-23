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
    namespace
    {
        inline bool case_is_ambiguous(int flag)
        {
            // Opposite-corner patterns 0b0101 (5) and 0b1010 (10)
            return flag == 5 || flag == 10;
        }

        inline bool asymptotic_decider(double f0, double f1, double f2, double f3)
        {
            // Bilinear form for marching squares; true -> choose diagonal (0,2), false -> (1,3)
            // This form is invariant to scaling of phi.
            const double d = f0 * f2 - f1 * f3;
            // Tie-break: when d==0 (center exactly on interface), prefer the alternate variant
            // to keep the opposite-corners case disconnected.
            return d > 0.0;
        }

        template <std::floating_point T>
        void compute_intersection_points(const std::span<const T> vertex_coordinates, const int gdim,
                                          const std::span<const T> ls_values, const int flag,
                                          std::vector<T>& intersection_points,
                                          std::unordered_map<int, int>& vertex_case_map)
        {
            const auto& intersected = intersected_edges[flag];

            intersection_points.clear();
            vertex_case_map.clear();

            std::vector<T> v0(gdim);
            std::vector<T> v1(gdim);
            std::vector<T> ip(gdim);

            int ip_index = 0;
            for (int e = 0; e < 4; ++e)
            {
                if (intersected[e] == 0)
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
                {
                    intersection_points.push_back(ip[j]);
                }

                vertex_case_map[e] = ip_index;
                ++ip_index;
            }
        }
        template <std::floating_point T, int MaxVerts>
        void decode_range(const std::span<const T> vertex_coordinates, const int gdim,
                          const int cell_begin, const int cell_end,
                          const type* cell_types,
                          const int (*subcell_verts)[MaxVerts],
                          const std::span<const T> intersection_points,
                          bool triangulate, CutCell<T>& cut_cell,
                          const std::unordered_map<int, int>& edge_ip_map,
                          std::unordered_map<int, int>& token_to_local_out)
        {
            // Map from token (edge id or 100+vid) to local vertex index in cut_cell._vertex_coords
            std::unordered_map<int, int> token_to_local;

            // Prefill edge intersections so they get stable indices
            for (int e = 0; e < 4; ++e)
            {
                const auto it = edge_ip_map.find(e);
                if (it == edge_ip_map.end())
                    continue;

                const int token = it->first; // edge id
                const int ip_idx = it->second;
                const int local_idx = static_cast<int>(cut_cell._vertex_coords.size() / gdim);
                token_to_local[token] = local_idx;
                cut_cell._vertex_coords.resize(cut_cell._vertex_coords.size() + gdim);
                for (int j = 0; j < gdim; ++j)
                    cut_cell._vertex_coords[local_idx * gdim + j] = intersection_points[ip_idx * gdim + j];
            }

            auto get_local = [&](int token) -> int
            {
                auto it = token_to_local.find(token);
                if (it != token_to_local.end())
                    return it->second;

                const int local_idx = static_cast<int>(cut_cell._vertex_coords.size() / gdim);
                token_to_local[token] = local_idx;
                cut_cell._vertex_coords.resize(cut_cell._vertex_coords.size() + gdim);

                if (token < 100)
                {
                    // edge intersection point
                    const int ip_idx = edge_ip_map.at(token);
                    for (int j = 0; j < gdim; ++j)
                        cut_cell._vertex_coords[local_idx * gdim + j] = intersection_points[ip_idx * gdim + j];
                }
                else
                {
                    const int vid = token - 100;
                    for (int j = 0; j < gdim; ++j)
                        cut_cell._vertex_coords[local_idx * gdim + j] = vertex_coordinates[vid * gdim + j];
                }
                return local_idx;
            };

            for (int cell_idx = cell_begin; cell_idx < cell_end; ++cell_idx)
            {
                const type sub_type = cell_types[cell_idx];
                const int* verts = subcell_verts[cell_idx];

                std::vector<int> verts_local;
                for (int i = 0; i < MaxVerts; ++i)
                {
                    const int token = verts[i];
                    if (token == -1)
                        break;
                    verts_local.push_back(get_local(token));
                }

                if (triangulate && sub_type == type::quadrilateral && verts_local.size() == 4)
                {
                    // split into two triangles (0,1,2) and (0,2,3)
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

            token_to_local_out = std::move(token_to_local);
        }

        template <std::floating_point T, int MaxVerts>
        void decode_case(const std::span<const T> vertex_coordinates, const int gdim, int flag,
                         const int* case_offsets, const type* cell_types,
                         const int (*subcell_verts)[MaxVerts],
                         const std::span<const T> intersection_points,
                         bool triangulate, CutCell<T>& cut_cell,
                         const std::unordered_map<int, int>& edge_ip_map,
                         std::unordered_map<int, int>& token_to_local_out)
        {
            decode_range<T, MaxVerts>(vertex_coordinates, gdim,
                                      case_offsets[flag], case_offsets[flag + 1],
                                      cell_types, subcell_verts,
                                      intersection_points, triangulate, cut_cell, edge_ip_map, token_to_local_out);
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

        // Store parent geometry
        cut_cell._parent_cell_type = type::quadrilateral;
        cut_cell._parent_vertex_coords.resize(4 * gdim);
        cut_cell._parent_vertex_ids.resize(4);
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < gdim; ++j)
            {
                cut_cell._parent_vertex_coords[i * gdim + j] = vertex_coordinates[i * gdim + j];
            }
            cut_cell._parent_vertex_ids[i] = i;
        }

        // Compute intersections (shared for all parts)
        std::vector<T> intersection_points;
        std::unordered_map<int, int> edge_ip_map;
        compute_intersection_points(vertex_coordinates, gdim, ls_values, flag_interior, intersection_points, edge_ip_map);

        // Final token->local vertex index map used to populate CutCell::_vertex_parent_entity
        std::unordered_map<int, int> token_to_local;

        cut_cell._gdim = gdim;
        cut_cell._vertex_coords.clear();
        cut_cell._types.clear();
        cut_cell._connectivity.clear();

        const bool is_amb = case_is_ambiguous(flag_interior);
        const int variant = asymptotic_decider(ls_values[0], ls_values[1], ls_values[2], ls_values[3]) ? 0 : 1;

        if (cut_type_str == "phi=0")
        {
            cut_cell._tdim = 1;
            std::span<const T> ip_span(intersection_points.data(), intersection_points.size());
            if (is_amb)
            {
                const int amb_id = amb_case_id[flag_interior];
                const int begin = amb_range_interface[4 * amb_id + 2 * variant + 0];
                const int end = amb_range_interface[4 * amb_id + 2 * variant + 1];
                decode_range<T, 2>(vertex_coordinates, gdim, begin, end,
                                   subcell_type_interface, subcell_verts_interface,
                                   ip_span, triangulate, cut_cell, edge_ip_map, token_to_local);
            }
            else
            {
                decode_case<T, 2>(vertex_coordinates, gdim, flag_interior,
                                  case_subcell_offset_interface, subcell_type_interface,
                                  subcell_verts_interface, ip_span, triangulate, cut_cell, edge_ip_map, token_to_local);
            }
        }
        else if (cut_type_str == "phi<0")
        {
            cut_cell._tdim = 2;
            std::span<const T> ip_span(intersection_points.data(), intersection_points.size());
            if (is_amb)
            {
                const int amb_id = amb_case_id[flag_interior];
                const int begin = amb_range_inside[4 * amb_id + 2 * variant + 0];
                const int end = amb_range_inside[4 * amb_id + 2 * variant + 1];
                decode_range<T, 4>(vertex_coordinates, gdim, begin, end,
                                   subcell_type_inside, subcell_verts_inside,
                                   ip_span, triangulate, cut_cell, edge_ip_map, token_to_local);
            }
            else
            {
                decode_case<T, 4>(vertex_coordinates, gdim, flag_interior,
                                  case_subcell_offset_inside, subcell_type_inside,
                                  subcell_verts_inside, ip_span, triangulate, cut_cell, edge_ip_map, token_to_local);
            }
        }
        else if (cut_type_str == "phi>0")
        {
            cut_cell._tdim = 2;
            std::span<const T> ip_span(intersection_points.data(), intersection_points.size());
            // Note: the generated "outside" tables are keyed by the *interior* mask (phi<0),
            // i.e. they already represent the complement region.
            if (is_amb)
            {
                const int amb_id = amb_case_id[flag_interior];
                const int begin = amb_range_outside[4 * amb_id + 2 * variant + 0];
                const int end = amb_range_outside[4 * amb_id + 2 * variant + 1];
                decode_range<T, 4>(vertex_coordinates, gdim, begin, end,
                                   subcell_type_outside, subcell_verts_outside,
                                   ip_span, triangulate, cut_cell, edge_ip_map, token_to_local);
            }
            else
            {
                decode_case<T, 4>(vertex_coordinates, gdim, flag_interior,
                                  case_subcell_offset_outside, subcell_type_outside,
                                  subcell_verts_outside, ip_span, triangulate, cut_cell, edge_ip_map, token_to_local);
            }
        }
        else
        {
            throw std::invalid_argument("cutting type unknown");
        }

        cutcells::utils::create_vertex_parent_entity_map<T>(token_to_local, cut_cell._vertex_parent_entity);
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
