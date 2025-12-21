// Copyright (c) 2025 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include "cut_hexahedron.h"

#include "cell_types.h"
#include "cut_cell.h"

#include <array>
#include <concepts>
#include <stdexcept>
#include <vector>

namespace cutcells::cell::hexahedron
{
    namespace
    {
        // Minimal fallback: subdivide the VTK-ordered hexahedron into tetrahedra.
        // VTK vertex numbering (from cell_topology.h):
        // 0:(0,0,0) 1:(1,0,0) 2:(1,1,0) 3:(0,1,0) 4:(0,0,1) 5:(1,0,1) 6:(1,1,1) 7:(0,1,1)
        //
        // Decomposition using the body diagonal (0-6), yielding 6 tets:
        // (0,1,2,6), (0,2,3,6), (0,3,7,6), (0,7,4,6), (0,4,5,6), (0,5,1,6)
        inline constexpr std::array<std::array<int, 4>, 6> hex_to_tets = {{
            {0, 1, 2, 6},
            {0, 2, 3, 6},
            {0, 3, 7, 6},
            {0, 7, 4, 6},
            {0, 4, 5, 6},
            {0, 5, 1, 6},
        }};
    }

    template <std::floating_point T>
    void cut(const std::span<const T> vertex_coordinates, const int gdim,
             const std::span<const T> ls_values, const std::string& cut_type_str,
             CutCell<T>& cut_cell, bool triangulate)
    {
        if (vertex_coordinates.size() != static_cast<std::size_t>(8 * gdim))
            throw std::invalid_argument("hexahedron::cut expects 8 vertices");
        if (ls_values.size() != 8)
            throw std::invalid_argument("hexahedron::cut expects 8 level set values");

        CutCell<T> tet_mesh;
        tet_mesh._gdim = gdim;
        tet_mesh._tdim = 3;
        tet_mesh._vertex_coords.assign(vertex_coordinates.begin(), vertex_coordinates.end());

        tet_mesh._connectivity.resize(hex_to_tets.size());
        tet_mesh._types.resize(hex_to_tets.size(), type::tetrahedron);

        for (std::size_t i = 0; i < hex_to_tets.size(); ++i)
        {
            tet_mesh._connectivity[i] = {hex_to_tets[i][0], hex_to_tets[i][1], hex_to_tets[i][2], hex_to_tets[i][3]};
        }

        // Reuse existing per-tet cutting + merge pipeline.
        // This is a correctness-first fallback until table-driven hex clipping is implemented.
        cutcells::cell::recursive_cut<T>(tet_mesh, ls_values, cut_type_str, triangulate);
        cut_cell = std::move(tet_mesh);
    }

    template <std::floating_point T>
    void cut(const std::span<const T> vertex_coordinates, const int gdim,
             const std::span<const T> ls_values, const std::vector<std::string>& cut_type_str,
             std::vector<CutCell<T>>& cut_cell, bool triangulate)
    {
        cut_cell.resize(cut_type_str.size());
        for (std::size_t i = 0; i < cut_type_str.size(); ++i)
        {
            cut<T>(vertex_coordinates, gdim, ls_values, cut_type_str[i], cut_cell[i], triangulate);
        }
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
