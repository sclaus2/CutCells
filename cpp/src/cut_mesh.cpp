// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include "cut_mesh.h"




#include "utils.h"
#include <map>
#include <algorithm>

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

      //either two or one
      int num_parents = cut_cells._parent_map.size()/cut_cells._cut_cells.size();

      cut_mesh._num_cells = num_cells;
      cut_mesh._types.resize(num_cells);
      cut_mesh._parent_map.resize(num_cells*num_parents);

      int merged_vertex_id = 0;
      int sub_cell_offset = 0;
      int vertex_counter = 0;
      int element_offset = 0;
      int cnt = 0;

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

        //Map from vertex id in current cutcell to merged cutmesh
        std::map<int,int> local_merged_vertex_ids;

        for(int local_id=0;local_id<num_cut_cell_vertices;local_id++)
        {
          //check if vertex already exists to avoid doubling of vertices
          int id = cutcells::utils::vertex_exists<T>(cut_mesh._vertex_coords, cut_cell._vertex_coords, local_id, gdim);

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
              cut_mesh._connectivity[element_offset+j] = local_merged_vertex_ids[index];
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
      cutcells::mesh::CutCells<T> cut_cells;

      bool triangulate = true;

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

//-----------------------------------------------------------------------------
}
