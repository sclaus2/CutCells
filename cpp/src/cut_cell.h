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

namespace cutcells
{
    namespace cell
    {
        /// @brief Stores the sub-cells resulting from cutting a cell
        struct CutCell
        {
            /// Geometric Dimension of Cell
            int _gdim;

            /// Topological Dimension of Cell
            int _tdim;

            /// Coordinates of vertices of cut cell
            std::vector<double> _vertex_coords;

            /// Vertex ids of cut cells
            std::vector<std::vector<int>> _connectivity;

            /// Cell type of cut cells
            std::vector<type> _types;

            /// Parent index for cell, pair of indices for interfaces
            std::vector<std::int32_t> _parent_cell_index;

        };

        void str(const CutCell &cut_cell);

        void sub_cell_vertices(const CutCell &cut_cell, const int& id, std::vector<double>& vertex_coordinates);
        double volume(const CutCell &cut_cell);

        void cut(const type cell_type, const std::span<const double> vertex_coordinates, const int gdim, 
                 const std::span<const double> ls_values, const std::string& cut_type_str,
                CutCell& cut_cell, bool triangulate=false);

        void cut(const type cell_type, const std::span<const double> vertex_coordinates, const int gdim, 
             const std::span<const double> ls_values, const std::vector<std::string>& cut_type_str,
             std::vector<CutCell>& cut_cell, bool triangulate=false);

        CutCell higher_order_cut(const type cell_type, const std::span<const double> vertex_coordinates, const int gdim, 
             const std::span<const double> ls_values, const std::string& cut_type_str,
             bool triangulate=false);

        CutCell merge(std::vector<CutCell> cut_cell_vec);

        CutCell create_cut_cell(const type& cell_type, std::span<const double> vertex_coords, const int& gdim);

         void cut_cut_cell(cutcells::cell::CutCell &cut_cell,
                    std::span<const double> ls_vals_all,
                    const int& parent_cell_index,
                    const std::string& cut_type_str,
                    bool triangulate);
    }

}