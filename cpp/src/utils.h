// Copyright (c) 2022-2023 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <cmath>
#include <concepts>
#include <span>
#include <unordered_map>
#include <vector>

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

    // create map from vertex coordinates to local parent entity, i.e. intersected edges and vertices
    template <std::floating_point T>
    void create_vertex_parent_entity_map(const std::unordered_map<int,int>& vertex_case_map, std::vector<int32_t>& vertex_parent_entity)
    {
      vertex_parent_entity.resize(vertex_case_map.size());

      for (auto &[first,second]: vertex_case_map) //using structured binding
      {
        vertex_parent_entity[second] = first;  // values can also be manipulated as they are refrences
      }
    }
}