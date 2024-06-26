// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include "cut_mesh.h"

#include "utils.h"
#include <map>

namespace cutcells::mesh
{
    void str(const CutCells &cut_mesh)
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

    int get_num_cells(const cutcells::mesh::CutCells& cut_mesh)
    {
      int num_cells = 0;

      for(auto & cut_cell : cut_mesh._cut_cells)
      {
        num_cells += cut_cell._connectivity.size();
      }
      return num_cells;
    }

    //create cutmesh from cutcells by merging all cut cells
    cutcells::mesh::CutMesh create_cut_mesh(std::vector<const cell::CutCell>& cut_cells)
    {
      CutMesh cut_mesh;

      std::size_t gdim = cut_cells[0]._gdim;
      std::size_t tdim = cut_cells[0]._tdim;

      cut_mesh._gdim = cut_cells[0]._gdim;
      cut_mesh._tdim = cut_cells[0]._tdim;

      //Count the total number of cells in vector
      int num_cells =0;
      for(auto & cut_cell :  cut_cells)
      {
        num_cells += cut_cell._connectivity.size();
      }

      cut_mesh._num_cells = num_cells;
      cut_mesh._connectivity.resize(num_cells);
      cut_mesh._parent_cell_index.resize(num_cells);
      cut_mesh._types.resize(num_cells);

      int merged_vertex_id = 0;
      int sub_cell_offset = 0;
      int vertex_counter = 0;
      int element_offset = 0;

      //all cutcells in vector above should have the same gdim and tdim
      for(auto & cut_cell : cut_cells)
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

        //Map from vertex id in current cutcell to merged cutmesh
        std::map<int,int> local_merged_vertex_ids;

        for(int local_id=0;local_id<num_cut_cell_vertices;local_id++)
        {
          //check if vertex already exists to avoid doubling of vertices
          int id = cutcells::utils::vertex_exists(cut_mesh._vertex_coords, cut_cell._vertex_coords, local_id, gdim);

          if(id==-1) //not found
          {
            //add vertex
            merged_vertex_id=vertex_counter;
            vertex_counter++;

            for(int j=0;j<gdim;j++)
            {
              cut_mesh._vertex_coords.push_back(cut_cell._vertex_coords[local_id*gdim+j]);
            }
          }
          else //found
          {
            //take already existing vertex for local mapping
            merged_vertex_id = id;
          }
          //offset is vertex_id
          local_merged_vertex_ids[local_id] = merged_vertex_id;
        }

        for(int i=0;i<local_num_cells;i++)
        {
          int num_vertices = cut_cell._connectivity[i].size();
          for(int j=0;j<num_vertices;j++)
          {
              int64_t index = cut_cell._connectivity[i][j];
              cut_mesh._connectivity[sub_cell_offset+i].push_back(local_merged_vertex_ids[index]);
          }

          // cut_mesh._offsets[sub_cell_offset+i] = element_offset;
          // element_offset += num_vertices;
          // one type per cell
          cut_mesh._types[sub_cell_offset+i]=cut_cell._types[i];
          //@todo: extend this for interface mesh here only first parent cell is kept!
          if(cut_cell._parent_cell_index.size()>0) //check if parent cell index is set
            cut_mesh._parent_cell_index[sub_cell_offset+i] = cut_cell._parent_cell_index[0];
        }

        sub_cell_offset+=local_num_cells;
      }

      cut_mesh._num_vertices = cut_mesh._vertex_coords.size()/cut_mesh._gdim;

      return cut_mesh;
    }
}