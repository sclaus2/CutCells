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

    template <std::floating_point T>
    void write_vtk(std::string filename,
                   const std::span<const T> element_vertex_coords,
                   const std::span<const int> connectivity,
                   const std::span<const int> offsets,
                   const std::span<cell::type> element_types,
                   const int gdim)
    {
        std::vector<double> coords_d(element_vertex_coords.begin(),
                                     element_vertex_coords.end());
        write_vtk(
            filename,
            std::span<const double>(coords_d.data(), coords_d.size()),
            connectivity,
            offsets,
            element_types,
            gdim);
    }

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

    /// @brief Position of node (i, j, k), 0 <= i, j, k <= p, of a degree-p
    /// Lagrange hexahedron in VTK's node order (vertices, edges, faces,
    /// interior), as in VTK 9.1 and later.
    inline int vtk_lagrange_hexahedron_index(int i, int j, int k, int p)
    {
        const bool ib = i == 0 || i == p, jb = j == 0 || j == p, kb = k == 0 || k == p;
        const int nb = int(ib) + int(jb) + int(kb);
        if (nb == 3)
            return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
        int offset = 8;
        if (nb == 2)
        {
            if (!ib)
                return (i - 1) + (j ? 2 * (p - 1) : 0) + (k ? 4 * (p - 1) : 0) + offset;
            if (!jb)
                return (j - 1) + (i ? (p - 1) : 3 * (p - 1)) + (k ? 4 * (p - 1) : 0) + offset;
            offset += 8 * (p - 1);
            return (k - 1) + (p - 1) * (i ? (j ? 2 : 1) : (j ? 3 : 0)) + offset;
        }
        offset += 12 * (p - 1);
        const int f = (p - 1) * (p - 1);
        if (nb == 1)
        {
            if (ib)
                return (j - 1) + (p - 1) * (k - 1) + (i ? f : 0) + offset;
            offset += 2 * f;
            if (jb)
                return (i - 1) + (p - 1) * (k - 1) + (j ? f : 0) + offset;
            offset += 2 * f;
            return (i - 1) + (p - 1) * (j - 1) + (k ? f : 0) + offset;
        }
        offset += 6 * f;
        return offset + (i - 1) + (p - 1) * ((j - 1) + (p - 1) * (k - 1));
    }

    /// @brief Position of node (i, j), 0 <= i, j <= p, of a degree-p Lagrange
    /// quadrilateral in VTK's node order.
    inline int vtk_lagrange_quadrilateral_index(int i, int j, int p)
    {
        const bool ib = i == 0 || i == p, jb = j == 0 || j == p;
        if (ib && jb)
            return i ? (j ? 2 : 1) : (j ? 3 : 0);
        const int offset = 4;
        if (!ib && jb)
            return (i - 1) + (j ? 2 * (p - 1) : 0) + offset;
        if (ib && !jb)
            return (j - 1) + (i ? (p - 1) : 3 * (p - 1)) + offset;
        return offset + 4 * (p - 1) + (i - 1) + (p - 1) * (j - 1);
    }

    /// @brief Write cells in VTK's Lagrange (or linear) node order to a .vtu file.
    ///
    /// Files with Lagrange quadrilaterals, hexahedra or wedges are written as
    /// version 2.2, which tells VTK that their node order is the one of VTK 9.1
    /// and later. @p degrees, when given, is written as HigherOrderDegrees.
    void write_lagrange_vtk(std::string filename,
                            const std::span<const double> point_coords,
                            const std::span<const int> connectivity,
                            const std::span<const int> offsets,
                            const std::span<const int> vtk_types,
                            int gdim,
                            const std::span<const std::int32_t> parent_map = {},
                            const std::span<const std::int32_t> subdivision_depth = {},
                            const std::span<const std::int32_t> degrees = {});

    template <std::floating_point T>
    void write_lagrange_vtk(std::string filename,
                            const std::span<const T> point_coords,
                            const std::span<const int> connectivity,
                            const std::span<const int> offsets,
                            const std::span<const int> vtk_types,
                            int gdim,
                            const std::span<const std::int32_t> parent_map = {},
                            const std::span<const std::int32_t> subdivision_depth = {},
                            const std::span<const std::int32_t> degrees = {})
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
            subdivision_depth,
            degrees);
    }
}
