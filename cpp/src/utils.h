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
      // vertex_case_map maps (token -> local_vertex_index). Tokens encode the origin:
      // - edge intersections: 0..(n_edges-1)
      // - original vertices:  100+vid
      // - special points:     200+sid
      // Local indices are expected to refer to indices in CutCell::_vertex_coords.
      int max_local = -1;
      for (const auto& kv : vertex_case_map)
        max_local = std::max(max_local, kv.second);

      vertex_parent_entity.assign(static_cast<std::size_t>(max_local + 1), -1);
      for (const auto& kv : vertex_case_map)
        vertex_parent_entity[static_cast<std::size_t>(kv.second)] = kv.first;
    }
}