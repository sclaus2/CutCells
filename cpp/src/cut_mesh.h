// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <iostream>
#include <span>
#include <vector>
#include <unordered_map>

#include "cut_cell.h"

#include "cell_types.h"
#include "cell_flags.h"

/// Structures to represent mesh formed from several local cell meshes
namespace cutcells::mesh
{
    /// Collection of cut cells that have been cut with regards to a parent mesh/entities
    template <std::floating_point T>
    struct CutCells
    {
      /// vector of all cut cells
      std::vector<cell::CutCell<T>> _cut_cells;

      /// map of cut cell id to parent cell id
      /// this vector points to parent cell for each cut cell in vector _cut_cells
      /// for an interface their are two parent cells
      std::vector<std::int32_t> _parent_map;

      /// the types of elements contained in all cut_cells
      std::vector<cell::type> _types;

    };

    // Class to represent mesh formed by merging CutCells (removal of double vertices + flattening of connectivties)
    template <std::floating_point T>
    struct CutMesh
    {
      /// Geometric Dimension of Cell
      int _gdim;

      /// Topological Dimension of Cell
      int _tdim;

      /// Number of cells
      int _num_cells;

      /// Number of vertices
      int _num_vertices;

      /// Coordinates of vertices of cut cell
      std::vector<T> _vertex_coords;

      /// Vertex ids of cut cells
      /// @todo: maybe change this to connectivity and offset vectors
      std::vector<int> _connectivity;
      std::vector<int> _offset; //offsets in connectivity vector

      /// Cell types
      std::vector<cell::type> _types;

      /// Parent index for cell, pair of indices for interfaces
      std::vector<std::int32_t> _parent_map;
    };

    /// @brief  Print information about cut_mesh to screen
    /// @param cut_mesh
    template <std::floating_point T>
    void str(const CutCells<T> &cut_mesh);

    /// @brief Get inverse map of parent_map
    /// @param parent_map: map from cutcell index to parent cell index.
    /// @return inverse map from parent cell index to cut cell index
    std::unordered_map<int, std::vector<int>> create_parent_cut_cells_map(std::span<int> parent_map);

    //Get number of cells in CutMesh
    template <std::floating_point T>
    int get_num_cells(const cutcells::mesh::CutCells<T>& cut_mesh);

    template <std::floating_point T>
    cutcells::mesh::CutMesh<T> create_cut_mesh(CutCells<T>& cut_cells);

    template <std::floating_point T>
    std::vector<int> locate_cells(std::span<const T> ls_vals, std::span<const T> points,
                                    std::span<const int> connectivity, std::span<const int> offset,
                                    std::span<const int> vtk_type,
                                    cell::cut_type ctype);

    template <std::floating_point T>
    cutcells::mesh::CutMesh<T> cut_vtk_mesh(std::span<const T> ls_vals, std::span<const T> points,
                                            std::span<const int> connectivity, std::span<const int> offset,
                                            std::span<const int> vtk_type,
                                            const std::string& cut_type_str);
}