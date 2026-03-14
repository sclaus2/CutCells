// Copyright (c) 2022-2023 ONERA
// Authors: Susanne Claus
// This file is part of CutFEMx
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <array>

//These tables describe the subdivision of cells by inserting one node at each edge 
//of the cell
namespace cutcells::cell
{
  inline constexpr std::array<std::array<int, 2>, 2> interval_subdivision_table = {{
    {0, 2},
    {2, 1}
  }};

  inline constexpr std::array<std::array<int, 3>, 4> triangle_subdivision_table = {{
    {0, 5, 4},
    {5, 1, 3},
    {4, 3, 2},
    {5, 3, 4}
  }};

//Notation: 0,1,2,3,e01(9),e02(8),e03(7),e12(6),e13(5),e23(4)
//where exy denotes middle point on edge between x and y
//Node ordering is chosen according to basix (fenicsx)
  inline constexpr std::array<std::array<int, 4>, 8> tetrahedron_subdivision_table = {{
    {0, 9, 8, 7}, {1, 9, 5, 6}, {2, 6, 4, 8}, {3, 7, 4, 5},
    {9, 6, 8, 7}, {9, 5, 7, 8}, {6, 8, 7, 4}, {6, 7, 5, 4}
  }};

// Quadrilateral P2 node numbering in Basix ordering:
//   0,1,2,3 (vertices)
//   e01=4, e02=5, e13=6, e23=7 (edge midpoints; Basix edge order)
//   f0123=8 (cell center)
//
// Subdivision into 4 quads (each keeps hexa/quad vertex layout used in this codebase).
  inline constexpr std::array<std::array<int, 4>, 4> quadrilateral_subdivision_table = {{
    {0, 4, 5, 8},
    {4, 1, 8, 6},
    {5, 8, 2, 7},
    {8, 6, 7, 3}
  }};

// Hexahedron P2 node numbering in Basix ordering:
//   Vertices: 0..7
//   Edge midpoints (Basix edge order in cell_topology.h): 8..19
//   Face centers (Basix face order): 20..25
//   Cell center: 26
//
// 8-way octree-like subdivision into hexahedra.
  inline constexpr std::array<std::array<int, 8>, 8> hexahedron_subdivision_table = {{
    {0, 8, 9, 20, 10, 21, 22, 26},
    {8, 1, 20, 11, 21, 12, 26, 23},
    {9, 20, 2, 13, 22, 26, 14, 24},
    {20, 11, 13, 3, 26, 23, 24, 15},
    {10, 21, 22, 26, 4, 16, 17, 25},
    {21, 12, 26, 23, 16, 5, 25, 18},
    {22, 26, 14, 24, 17, 25, 6, 19},
    {26, 23, 24, 15, 25, 18, 19, 7}
  }};

// basix -> vtk map   0, 1, 2, 3, 9, 8, 5, 7, 6, 4
//  basix: e01: 9, e12: 6, e02: 8, e03: 7, e13: 5, e23: 4
//  vtk: e01 : 4 , e12: 5, e02: 6, e03: 7, e13: 8, e23: 9
//  vtk subdivision: {0,4,6,7},{1,4,8,5},{2,5,9,6},{3,7,9,8},{4,5,6,7},{4,5,7,8},{5,6,7,9},{5,7,8,9}
// basix subdivision: {0,9,8,7},{1,9,5,6},{2,6,4,8},{3,7,4,5},{9,6,8,7},{9,5,7,8},{6,8,7,4},{6,7,5,4}
}
