// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <iostream>
#include <cstdint>

#include "cell_types.h"
#include <vector>
#include <span>
#include <string>
#include <concepts>

namespace cutcells
{
    // Structure to represent mesh local to each cell
    namespace cell
    {
        /// @brief Stores the sub-cells resulting from cutting a cell
        // Note: can also be used to merge range of cutcells into one mesh representation
        template <std::floating_point T>
        struct CutCell
        {
            /// Geometric Dimension of Cell
            int _gdim;

            /// Topological Dimension of Cell
            int _tdim;

            /// Coordinates of vertices of cut cell
            std::vector<T> _vertex_coords;

            /// Vertex ids of cut cells
            /// Flattened cell-to-vertex connectivity in CSR layout
            std::vector<int> _connectivity;

            /// Offsets into _connectivity (size = num_cells + 1)
            std::vector<int> _offset;

            /// Cell type of cut cells
            std::vector<type> _types;

            /// parent facet or vertex of vertex coordinates
            /// these correspond to the local numbering of intersected facets or vertices
            /// vertices are indicated by a shift of 100 in the numbering  vertex 1 -> 101 etc.
            std::vector<int32_t> _vertex_parent_entity;

            /// Parent cell index
            int32_t _parent_cell_index;

            /// Parent cell type before cutting
            type _parent_cell_type;

            /// Parent vertex coordinates (V*gdim)
            std::vector<T> _parent_vertex_coords;

            /// Parent vertex IDs (context-global ids)
            std::vector<int> _parent_vertex_ids;
        };

        template <std::floating_point T>
        void str(const CutCell<T> &cut_cell);

        template <std::floating_point T>
        void sub_cell_vertices(const CutCell<T> &cut_cell, const int& id, std::vector<T>& vertex_coordinates);

        template <std::floating_point T>
        T volume(const CutCell<T> &cut_cell);

        template <std::floating_point T>
        void cut(const type cell_type, const std::span<const T> vertex_coordinates, const int gdim,
                 const std::span<const T> ls_values, const std::string& cut_type_str,
                CutCell<T>& cut_cell, bool triangulate=false);

        template <std::floating_point T>
        void cut(const type cell_type, const std::span<const T> vertex_coordinates, const int gdim,
             const std::span<const T> ls_values, const std::vector<std::string>& cut_type_str,
             std::vector<CutCell<T>>& cut_cell, bool triangulate=false);

        template <std::floating_point T>
        CutCell<T> higher_order_cut(const type cell_type, const std::span<const T> vertex_coordinates, const int gdim,
             const std::span<const T> ls_values, const std::string& cut_type_str,
             bool triangulate=false);

        template <std::floating_point T>
        CutCell<T> merge(std::vector<CutCell<T>> cut_cell_vec);

        template <std::floating_point T>
        CutCell<T> create_cut_cell(const type& cell_type, std::span<const T> vertex_coords, const int& gdim);

        template <std::floating_point T>
        void recursive_cut(cutcells::cell::CutCell<T> &cut_cell,
                  std::span<const T> ls_vals_all,
                  const std::string& cut_type_str,
                  bool triangulate);

                template <std::floating_point T>
                inline int num_cells(const CutCell<T>& cut_cell)
                {
                    return cut_cell._offset.empty() ? 0 : static_cast<int>(cut_cell._offset.size()) - 1;
                }

                template <std::floating_point T>
                inline int num_cell_vertices(const CutCell<T>& cut_cell, const int cell_id)
                {
                    return cut_cell._offset[cell_id + 1] - cut_cell._offset[cell_id];
                }

                template <std::floating_point T>
                inline std::span<const int> cell_vertices(const CutCell<T>& cut_cell, const int cell_id)
                {
                    const int begin = cut_cell._offset[cell_id];
                    const int end = cut_cell._offset[cell_id + 1];
                    return std::span<const int>(cut_cell._connectivity.data() + begin, end - begin);
                }

                template <std::floating_point T>
                inline void clear_cell_topology(CutCell<T>& cut_cell)
                {
                    cut_cell._connectivity.clear();
                    cut_cell._offset.clear();
                    cut_cell._offset.push_back(0);
                    cut_cell._types.clear();
                }

                template <std::floating_point T>
                inline void reserve_cell_topology(CutCell<T>& cut_cell,
                                                                                    const int connectivity_capacity,
                                                                                    const int cell_capacity)
                {
                    cut_cell._connectivity.reserve(connectivity_capacity);
                    cut_cell._offset.reserve(cell_capacity + 1);
                    cut_cell._types.reserve(cell_capacity);
                }

                template <std::floating_point T>
                inline void append_cell(CutCell<T>& cut_cell, const type cell_type, std::span<const int> vertices)
                {
                    cut_cell._types.push_back(cell_type);
                    cut_cell._connectivity.insert(cut_cell._connectivity.end(), vertices.begin(), vertices.end());
                    cut_cell._offset.push_back(static_cast<int>(cut_cell._connectivity.size()));
                }

                /// Append a subcell given as a raw pointer + count.
                /// Avoids forcing a temporary std::vector<int> at call sites.
                template <std::floating_point T>
                inline void append_cell(CutCell<T>& cut_cell, const type cell_type, const int* vertices, int n)
                {
                    cut_cell._types.push_back(cell_type);
                    cut_cell._connectivity.insert(cut_cell._connectivity.end(), vertices, vertices + n);
                    cut_cell._offset.push_back(static_cast<int>(cut_cell._connectivity.size()));
                }

                /// Append a subcell given as a fixed-size std::array using the first n entries.
                template <std::floating_point T, std::size_t N>
                inline void append_cell(CutCell<T>& cut_cell, const type cell_type, const std::array<int, N>& vertices, int n)
                {
                    cut_cell._types.push_back(cell_type);
                    cut_cell._connectivity.insert(cut_cell._connectivity.end(), vertices.begin(), vertices.begin() + n);
                    cut_cell._offset.push_back(static_cast<int>(cut_cell._connectivity.size()));
                }
    }

}