// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <iostream>

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
            /// @todo: maybe change this to connectivity and offset vectors
            std::vector<std::vector<int>> _connectivity;

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
    }

}