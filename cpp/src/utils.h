// Copyright (c) 2022-2023 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <cmath>

namespace cutcells::utils
{
    //Check if two vertices are equal
    static bool equal(std::span<const double> coord1, const int &id1,
                      std::span<const double> coord2, const int &id2, const int& gdim)
    {
      double tol = 1e-12;
      double distance = 0;

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
    static int vertex_exists(std::span<const double> geom, std::span<const double> coords, const int &id2, const int &gdim)
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