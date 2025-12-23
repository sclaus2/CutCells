// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include "cut_mesh.h"




#include "utils.h"
#include "cell_topology.h"
#include <map>
#include <algorithm>

namespace
{
  struct VertexKey
  {
    // 0: original vertex (global vertex id)
    // 1: edge intersection (global vertex id pair)
    // 2: special point (parent cell id + sid)
    uint8_t kind = 0;
    int32_t a = 0;
    int32_t b = 0;

    bool operator==(const VertexKey& other) const noexcept
    {
      return kind == other.kind && a == other.a && b == other.b;
    }
  };

  struct VertexKeyHash
  {
    std::size_t operator()(const VertexKey& k) const noexcept
    {
      // Basic mix; good enough for small integer keys.
      const std::size_t h0 = std::hash<int32_t>{}(k.a);
      const std::size_t h1 = std::hash<int32_t>{}(k.b);
      const std::size_t hk = std::hash<uint8_t>{}(k.kind);
      return h0 ^ (h1 + 0x9e3779b97f4a7c15ULL + (h0 << 6) + (h0 >> 2)) ^ (hk << 1);
    }
  };

  template <std::floating_point T>
  inline void append_vertex_coords(std::vector<T>& out, const std::vector<T>& in,
                   int local_vertex_id, int gdim)
  {
    const int base = local_vertex_id * gdim;
    for (int j = 0; j < gdim; ++j)
      out.push_back(in[base + j]);
  }
}

namespace cutcells::mesh
{
    template <std::floating_point T>
    void str(const CutCells<T> &cut_mesh)
    {
        std::cout << "CutCells: " << std::endl;
        int cnt = 0;

        for(auto &cell: cut_mesh._cut_cells)
        {
            std::cout << "Cut Cell ";
            std::cout << cnt << ": ";
            cnt++;
            std::cout << "vertex coordinates=[";
            for(auto &i: cell._vertex_coords)
            {
                std::cout << i << ", ";
            }
            std::cout << "]" << std::endl;
            std::cout << "connectivity=[";
            for(int i=0;i<cell._connectivity.size();i++)
            {
                std::cout << i << ": ";
                for(int j=0;j<cell._connectivity[i].size();j++)
                {
                    std::cout << cell._connectivity[i][j] << ", ";
                }
            }
            std::cout << "]" << std::endl;
        }

        std::cout << "Cell Types= ";
        for(auto &type : cut_mesh._types)
        {
            std::cout << cell_type_to_str(type) << ", ";
        }
        std::cout << std::endl;

        std::cout << "Parent map=[";
        for(int i=0;i<cut_mesh._parent_map.size();i++)
        {
            std::cout << i << ": " << cut_mesh._parent_map[i] << std::endl;
        }
        std::cout << "]" << std::endl;
    }

    std::unordered_map<int, std::vector<int>> create_parent_cut_cells_map(const std::span<int> parent_map)
    {
        std::unordered_map<int, std::vector<int>> parent_cut_cell_map;

        for(std::size_t i=0; i< parent_map.size();i++)
        {
            parent_cut_cell_map[parent_map[i]].push_back(i);
        }

        return parent_cut_cell_map;
    }

    template <std::floating_point T>
    int get_num_cells(const cutcells::mesh::CutCells<T>& cut_mesh)
    {
      int num_cells = 0;

      for(auto & cut_cell : cut_mesh._cut_cells)
      {
        num_cells += cut_cell._connectivity.size();
      }
      return num_cells;
    }

    //create cutmesh from cutcells by merging cut cells
    template <std::floating_point T>
    cutcells::mesh::CutMesh<T> create_cut_mesh(CutCells<T>& cut_cells)
    {
      CutMesh<T> cut_mesh;

      if (cut_cells._cut_cells.empty())
      {
        cut_mesh._gdim = 0;
        cut_mesh._tdim = 0;
        cut_mesh._num_cells = 0;
        cut_mesh._num_vertices = 0;
        return cut_mesh;
      }

      std::size_t gdim = cut_cells._cut_cells[0]._gdim;
      std::size_t tdim = cut_cells._cut_cells[0]._tdim;

      cut_mesh._gdim = gdim;
      cut_mesh._tdim = tdim;

      //Count the total number of cells in vector
      int num_cells =0;
      for(auto & cut_cell :  cut_cells._cut_cells)
      {
        num_cells += cut_cell._connectivity.size();
      }

      int num_connectivity=0;
      for(auto & cut_cell :  cut_cells._cut_cells)
        for(int i=0;i<cut_cell._connectivity.size();i++)
          num_connectivity += cut_cell._connectivity[i].size();

      cut_mesh._offset.resize(num_cells+1);
      cut_mesh._connectivity.resize(num_connectivity);

      // either two or one; allow missing parent map
      int num_parents = 0;
      if (!cut_cells._parent_map.empty())
        num_parents = static_cast<int>(cut_cells._parent_map.size() / cut_cells._cut_cells.size());
      if (num_parents <= 0)
        num_parents = 1;

      cut_mesh._num_cells = num_cells;
      cut_mesh._types.resize(num_cells);
      cut_mesh._parent_map.resize(num_cells * num_parents, -1);

      int sub_cell_offset = 0;
      int element_offset = 0;
      int cnt = 0;

      // Fast global dedup: uses CutCell::_vertex_parent_entity tokens and CutCell::_parent_vertex_ids
      // (which should be context-global vertex ids for the parent mesh).
      std::unordered_map<VertexKey, int, VertexKeyHash> global_vertex_ids;
      global_vertex_ids.reserve(static_cast<std::size_t>(num_connectivity));

      //all cutcells in vector above should have the same gdim and tdim
      for(auto & cut_cell : cut_cells._cut_cells)
      {
        if(cut_cell._vertex_coords.size()==0)
        {
          continue;
        }
        //check that current cut_cell has same dimensions and parent
        if((cut_cell._gdim!=cut_mesh._gdim)||(cut_cell._tdim!=cut_mesh._tdim))
        {
          std::cout << "gdim: (" << gdim << ", " << cut_cell._gdim << ")" << std::endl;
          std::cout << "tdim: (" << tdim << ", " << cut_cell._tdim << ")" << std::endl;
          throw std::runtime_error ("Error in merging cutcells into mesh as they have differing dimensions");
        }

        int num_cut_cell_vertices = cut_cell._vertex_coords.size()/gdim;

        int local_num_cells = cut_cell._connectivity.size();

        // Map from vertex id in current cutcell to merged cutmesh
        std::vector<int> local_merged_vertex_ids(num_cut_cell_vertices, -1);

        // Determine whether we can use fast token-based dedup for this cut cell.
        const bool has_tokens = (static_cast<int>(cut_cell._vertex_parent_entity.size()) == num_cut_cell_vertices);
        const int parent_vertices = cell::get_num_vertices(cut_cell._parent_cell_type);
        const bool has_parent_ids = (static_cast<int>(cut_cell._parent_vertex_ids.size()) == parent_vertices);
        const int parent_cell_id = cut_cells._parent_map.empty()
                                       ? -1
                                       : cut_cells._parent_map[cnt * num_parents + 0];
        const bool can_fast_dedup = has_tokens && has_parent_ids && (parent_cell_id >= 0);

        if (can_fast_dedup)
        {
          const auto parent_edges = cell::edges(cut_cell._parent_cell_type);

          for (int local_id = 0; local_id < num_cut_cell_vertices; ++local_id)
          {
            const int32_t token = cut_cell._vertex_parent_entity[local_id];
            VertexKey key;

            if (token >= 100 && token < 200)
            {
              const int ref_vid = static_cast<int>(token - 100);
              if (ref_vid < 0 || ref_vid >= parent_vertices)
              {
                // Fallback for unexpected token.
                key.kind = 2;
                key.a = parent_cell_id;
                key.b = token;
              }
              else
              {
                key.kind = 0;
                key.a = static_cast<int32_t>(cut_cell._parent_vertex_ids[ref_vid]);
                key.b = 0;
              }
            }
            else if (token >= 200)
            {
              key.kind = 2;
              key.a = parent_cell_id;
              key.b = static_cast<int32_t>(token - 200);
            }
            else
            {
              const int edge_id = static_cast<int>(token);
              if (edge_id < 0 || edge_id >= static_cast<int>(parent_edges.size()))
              {
                key.kind = 2;
                key.a = parent_cell_id;
                key.b = token;
              }
              else
              {
                const int v0 = parent_edges[edge_id][0];
                const int v1 = parent_edges[edge_id][1];
                const int32_t gv0 = static_cast<int32_t>(cut_cell._parent_vertex_ids[v0]);
                const int32_t gv1 = static_cast<int32_t>(cut_cell._parent_vertex_ids[v1]);
                key.kind = 1;
                key.a = std::min(gv0, gv1);
                key.b = std::max(gv0, gv1);
              }
            }

            auto it = global_vertex_ids.find(key);
            if (it == global_vertex_ids.end())
            {
              const int merged_vertex_id = static_cast<int>(cut_mesh._vertex_coords.size() / gdim);
              append_vertex_coords<T>(cut_mesh._vertex_coords, cut_cell._vertex_coords, local_id, static_cast<int>(gdim));
              global_vertex_ids.emplace(key, merged_vertex_id);
              local_merged_vertex_ids[local_id] = merged_vertex_id;
            }
            else
            {
              local_merged_vertex_ids[local_id] = it->second;
            }
          }
        }
        else
        {
          // Slow fallback: coordinate-based dedup (kept for compatibility).
          for (int local_id = 0; local_id < num_cut_cell_vertices; ++local_id)
          {
            const int id = cutcells::utils::vertex_exists<T>(cut_mesh._vertex_coords, cut_cell._vertex_coords,
                                                            local_id, static_cast<int>(gdim));
            if (id == -1)
            {
              const int merged_vertex_id = static_cast<int>(cut_mesh._vertex_coords.size() / gdim);
              append_vertex_coords<T>(cut_mesh._vertex_coords, cut_cell._vertex_coords, local_id, static_cast<int>(gdim));
              local_merged_vertex_ids[local_id] = merged_vertex_id;
            }
            else
            {
              local_merged_vertex_ids[local_id] = id;
            }
          }
        }

        for(int i=0;i<local_num_cells;i++)
        {
          int num_vertices = cut_cell._connectivity[i].size();
          for(int j=0;j<num_vertices;j++)
          {
              int64_t index = cut_cell._connectivity[i][j];
              cut_mesh._connectivity[element_offset+j] = local_merged_vertex_ids[static_cast<int>(index)];
          }

          cut_mesh._offset[sub_cell_offset+i] = element_offset;
          element_offset += num_vertices;
          // one type per cell
          cut_mesh._types[sub_cell_offset+i]=cut_cell._types[i];

          for(int j=0;j<num_parents;j++)
            cut_mesh._parent_map[sub_cell_offset+i*num_parents+j] = cut_cells._parent_map[cnt*num_parents+j];
        }

        sub_cell_offset+=local_num_cells;
        cnt++;
      }

      cut_mesh._offset[num_cells]=element_offset; //last offset

      cut_mesh._num_vertices = cut_mesh._vertex_coords.size()/cut_mesh._gdim;

      return cut_mesh;
    }

    //extract the cell ids that are in a given domain
    // either phi<0, phi = 0, or phi>0
    template <std::floating_point T>
    std::vector<int> locate_cells(std::span<const T> ls_vals, std::span<const T> points,
                                    std::span<const int> connectivity, std::span<const int> offset,
                                    std::span<const int> vtk_type,
                                    cell::cut_type ctype)
    {
      std::vector<int> cells;

      int num_cells = vtk_type.size();
      for(std::size_t i=0;i<num_cells;i++)
      {
        int cell_offset = offset[i];
        cell::vtk_types vtype = static_cast<cell::vtk_types>(vtk_type[i]);
        cell::type cell_type = cell::map_vtk_type_to_cell_type(vtype);

        int num_vertices = cell::get_num_vertices(cell_type);

        std::vector<T> level_set_values(num_vertices);
        for(std::size_t j=0;j<num_vertices;j++)
        {
          int vertex_id = connectivity[cell_offset+j];
          level_set_values[j] = ls_vals[vertex_id];
        }

        cell::domain domain_type = cell::classify_cell_domain<T>(std::span(level_set_values.data(),level_set_values.size()));

        if((domain_type == cell::domain::inside) && (ctype == cell::cut_type::philt0))
        {
          cells.push_back(i);
        }
        if((domain_type == cell::domain::outside) && (ctype == cell::cut_type::phigt0))
        {
          cells.push_back(i);
        }
        if((domain_type == cell::domain::intersected) && (ctype == cell::cut_type::phieq0))
        {
          cells.push_back(i);
        }
      }

      return cells;
    }

    // cut a vtk mesh
    // input: level set values at points
    // connectivity of cells with offsets and type
    template <std::floating_point T>
    cutcells::mesh::CutMesh<T> cut_vtk_mesh(std::span<const T> ls_vals, std::span<const T> points,
                                            std::span<const int> connectivity, std::span<const int> offset,
                                            std::span<const int> vtk_type,
                                            const std::string& cut_type_str)
    {
      return cut_vtk_mesh<T>(ls_vals, points,
                             connectivity, offset,
                             vtk_type,
                             cut_type_str,
                             true);
    }

    template <std::floating_point T>
    cutcells::mesh::CutMesh<T> cut_vtk_mesh(std::span<const T> ls_vals, std::span<const T> points,
                                            std::span<const int> connectivity, std::span<const int> offset,
                                            std::span<const int> vtk_type,
                                            const std::string& cut_type_str,
                                            bool triangulate)
    {
      cutcells::mesh::CutCells<T> cut_cells;

      auto intersected_cells = locate_cells<T>(ls_vals, points,
                                               connectivity, offset,
                                               vtk_type,
                                               cell::cut_type::phieq0);

      cut_cells._cut_cells.resize(intersected_cells.size());
      cut_cells._parent_map.resize(intersected_cells.size());

      for(std::size_t i=0;i<intersected_cells.size();i++)
      {
        int cell_index = intersected_cells[i];
        int cell_offset = offset[cell_index];
        cell::vtk_types cell_vtk_type = static_cast<cell::vtk_types>(vtk_type[cell_index]);
        cell::type cell_type = cell::map_vtk_type_to_cell_type(cell_vtk_type);
        int num_vertices = cell::get_num_vertices(cell_type);

        std::vector<T> vertex_coords(num_vertices*3);
        std::vector<T> level_set_values(num_vertices);

        for(std::size_t j=0;j<num_vertices;j++)
        {
          int vertex_id = connectivity[cell_offset+j];
          for(std::size_t k=0;k<3;k++)
            vertex_coords[j*3+k] = points[vertex_id*3+k];

          level_set_values[j] = ls_vals[vertex_id];
        }

        cell::CutCell<T> cut_cell;

        //cut the cell
        cell::cut<T>(cell_type, vertex_coords,3,level_set_values,cut_type_str,cut_cell,triangulate);

        // Populate context-global parent vertex IDs to enable fast merging in create_cut_mesh.
        // The cutters themselves default these IDs to 0..V-1.
        cut_cell._parent_cell_type = cell_type;
        cut_cell._parent_vertex_ids.resize(num_vertices);
        for (int j = 0; j < num_vertices; ++j)
        {
          const int vertex_id = connectivity[cell_offset + j];
          cut_cell._parent_vertex_ids[j] = vertex_id;
        }

        cut_cells._cut_cells[i] = cut_cell;
        cut_cells._parent_map[i] = cell_index;

        if (std::find(cut_cells._types.begin(), cut_cells._types.end(), cell_type) == cut_cells._types.end())
          cut_cells._types.push_back(cell_type);
      }

      return create_cut_mesh<T>(cut_cells);
    }

//-----------------------------------------------------------------------------
    template void str(const CutCells<double> &cut_mesh);
    template void str(const CutCells<float> &cut_mesh);

    template cutcells::mesh::CutMesh<double> create_cut_mesh(CutCells<double>& cut_cells);
    template cutcells::mesh::CutMesh<float> create_cut_mesh(CutCells<float>& cut_cells);

    template
    std::vector<int> locate_cells(std::span<const double> ls_vals, std::span<const double> points,
                                    std::span<const int> connectivity, std::span<const int> offset,
                                    std::span<const int> vtk_type,
                                    cell::cut_type ctype);

    template
    std::vector<int> locate_cells(std::span<const float> ls_vals, std::span<const float> points,
                                    std::span<const int> connectivity, std::span<const int> offset,
                                    std::span<const int> vtk_type,
                                    cell::cut_type ctype);

    template
    cutcells::mesh::CutMesh<double> cut_vtk_mesh(std::span<const double> ls_vals, std::span<const double> points,
                                            std::span<const int> connectivity, std::span<const int> offset,
                                            std::span<const int> vtk_type,
                                            const std::string& cut_type_str);
    template
    cutcells::mesh::CutMesh<float> cut_vtk_mesh(std::span<const float> ls_vals, std::span<const float> points,
                                            std::span<const int> connectivity, std::span<const int> offset,
                                            std::span<const int> vtk_type,
                                            const std::string& cut_type_str);

    template
    cutcells::mesh::CutMesh<double> cut_vtk_mesh(std::span<const double> ls_vals, std::span<const double> points,
                        std::span<const int> connectivity, std::span<const int> offset,
                        std::span<const int> vtk_type,
                        const std::string& cut_type_str,
                        bool triangulate);
    template
    cutcells::mesh::CutMesh<float> cut_vtk_mesh(std::span<const float> ls_vals, std::span<const float> points,
                        std::span<const int> connectivity, std::span<const int> offset,
                        std::span<const int> vtk_type,
                        const std::string& cut_type_str,
                        bool triangulate);

//-----------------------------------------------------------------------------
}
