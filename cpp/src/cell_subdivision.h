// Copyright (c) 2022-2023 ONERA 
// Authors: Susanne Claus 
// This file is part of CutFEMx
//
// SPDX-License-Identifier:    MIT
#include <array>

//These tables describe the subdivision of cells by inserting one node at each edge 
//of the cell
namespace cutcells::cell
{
  std::array<std::array<int, 2>, 2>  interval_subdivision_table{
    {{0,2},
     {2,1}}
  };

  std::array<std::array<int, 3>, 4> triangle_subdivision_table{
    {{0,5,4},
     {5,1,3},
     {4,3,2},
     {5,3,4}}
  };

//Notation: 0,1,2,3,e01(9),e02(8),e03(7),e12(6),e13(5),e23(4)
//where exy denotes middle point on edge between x and y
//Node ordering is chosen according to basix (fenicsx)
  std::array<std::array<int, 4>, 8> tetrahedron_subdivision_table{
    {// four tetrahedra from each vertex
      {0,9,8,7},{1,9,5,6},{2,6,4,8},{3,7,4,5},{9,6,8,7},{9,5,7,8},{6,8,7,4},{6,7,5,4}
    }
  };

// basix -> vtk map   0, 1, 2, 3, 9, 8, 5, 7, 6, 4
//  basix: e01: 9, e12: 6, e02: 8, e03: 7, e13: 5, e23: 4
//  vtk: e01 : 4 , e12: 5, e02: 6, e03: 7, e13: 8, e23: 9
//  vtk subdivision: {0,4,6,7},{1,4,8,5},{2,5,9,6},{3,7,9,8},{4,5,6,7},{4,5,7,8},{5,6,7,9},{5,7,8,9}
// basix subdivision: {0,9,8,7},{1,9,5,6},{2,6,4,8},{3,7,4,5},{9,6,8,7},{9,5,7,8},{6,8,7,4},{6,7,5,4}
}