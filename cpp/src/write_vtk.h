// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <string>
#include <concepts>
#include <cstdint>
#include <vector>
#include "cut_cell.h"
#include "cut_mesh.h"
#include "reference_cell.h"
#include "level_set.h"

namespace cutcells::io
{
    std::vector<int> basix_to_vtk_lagrange_permutation(cell::type cell_type,
                                                       int local_dofs,
                                                       int degree);

    void write_vtk(std::string filename, const std::span<const double> element_vertex_coords,  
                    const std::span<const int> connectivity,
                    const std::span<const int> offsets,
                    const std::span<cell::type> element_types, 
                    const int gdim);

    void write_vtk(std::string filename, cell::CutCell<double>& cut_cell);

    /// Write a CutMesh (basix-ordered internally) to a VTU file.
    /// Non-simplex cells are permuted to VTK vertex ordering before writing.
    template <std::floating_point T>
    void write_vtk(std::string filename, const mesh::CutMesh<T>& cut_mesh)
    {
        const int n_cells = cut_mesh._num_cells;
        const int gdim = cut_mesh._gdim;

        // Permute non-simplex cell connectivity from basix to VTK ordering.
        std::vector<int> vtk_connectivity;
        vtk_connectivity.reserve(cut_mesh._connectivity.size());

        for (int c = 0; c < n_cells; ++c)
        {
            const int start = (c == 0) ? 0 : cut_mesh._offset[static_cast<std::size_t>(c)];
            const int end   = cut_mesh._offset[static_cast<std::size_t>(c + 1)];
            const int nv    = end - start;
            const cell::type ctype = cut_mesh._types[static_cast<std::size_t>(c)];

            if (ctype == cell::type::point
                || ctype == cell::type::interval
                || ctype == cell::type::triangle
                || ctype == cell::type::tetrahedron)
            {
                for (int k = start; k < end; ++k)
                    vtk_connectivity.push_back(cut_mesh._connectivity[static_cast<std::size_t>(k)]);
            }
            else
            {
                const auto perm = cell::basix_to_vtk_vertex_permutation(ctype);
                if (static_cast<int>(perm.size()) != nv)
                    throw std::runtime_error("write_vtk(CutMesh): cell vertex count mismatch");
                for (int j = 0; j < nv; ++j)
                    vtk_connectivity.push_back(
                        cut_mesh._connectivity[static_cast<std::size_t>(
                            start + perm[static_cast<std::size_t>(j)])]);
            }
        }

        // write_vtk low-level function accepts only double coords; convert if needed.
        std::vector<double> coords_d(cut_mesh._vertex_coords.begin(),
                                     cut_mesh._vertex_coords.end());
        write_vtk(
            filename,
            std::span<const double>(coords_d.data(), coords_d.size()),
            std::span<const int>(vtk_connectivity.data(), vtk_connectivity.size()),
            std::span<const int>(cut_mesh._offset.data(), cut_mesh._offset.size()),
            std::span<cell::type>(const_cast<cell::type*>(cut_mesh._types.data()),
                                  cut_mesh._types.size()),
            gdim);
    }

    void write_level_set_vtu(std::string filename,
                             const cutcells::LevelSetFunction<double>& ls,
                             std::string field_name = "phi");

    void write_lagrange_vtk(std::string filename,
                            const std::span<const double> point_coords,
                            const std::span<const int> connectivity,
                            const std::span<const int> offsets,
                            const std::span<const int> vtk_types,
                            int gdim,
                            const std::span<const std::int32_t> parent_map = {},
                            const std::span<const std::int32_t> subdivision_depth = {});

    template <std::floating_point T>
    void write_lagrange_vtk(std::string filename,
                            const std::span<const T> point_coords,
                            const std::span<const int> connectivity,
                            const std::span<const int> offsets,
                            const std::span<const int> vtk_types,
                            int gdim,
                            const std::span<const std::int32_t> parent_map = {},
                            const std::span<const std::int32_t> subdivision_depth = {})
    {
        std::vector<double> coords_d(point_coords.begin(), point_coords.end());
        write_lagrange_vtk(
            filename,
            std::span<const double>(coords_d.data(), coords_d.size()),
            connectivity,
            offsets,
            vtk_types,
            gdim,
            parent_map,
            subdivision_depth);
    }
}
