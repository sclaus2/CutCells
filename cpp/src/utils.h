// Copyright (c) 2022-2023 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <cmath>
#include <concepts>

namespace cutcells::utils
{
    //Check if two vertices are equal
    template <std::floating_point T>
    static bool equal(std::span<const T> coord1, const int &id1,
                      std::span<const T> coord2, const int &id2, const int& gdim)
    {
      T tol = 1e-15;
      T distance = 0;

      //Take the distance between two points with id1 and id2 and return
      for(std::size_t j=0;j<gdim;j++)
      {
          distance +=(coord1[id1*gdim+j] - coord2[id2*gdim+j])*(coord1[id1*gdim+j] - coord2[id2*gdim+j]);
      }

      distance = sqrt(distance);

      if(distance<tol)
      {
        return true;
      }
      else
      {
        return false;
      }
    }

    // check if vertex coordinate exists in geom vector
    // and if vertex exists return its id in geom otherwise return -1
    template <std::floating_point T>
    static int vertex_exists(std::span<const T> geom, std::span<const T> coords, const int &id2, const int &gdim)
    {
      int num_vertices_geom = geom.size()/gdim;

      for(int i=0;i<num_vertices_geom;i++)
      {
        if(equal(geom, i, coords, id2,gdim))
        {
          return i;
        }
      }
      return -1;
    }
}