// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include "cut_cell.h"
#include "cut_tetrahedron.h"
#include "cut_triangle.h"
#include "cut_interval.h"
#include "cell_flags.h"
#include "cell_subdivision.h"
#include "utils.h"

#include <cassert>
#include <set>
#include <unordered_map>
#include <map>


namespace cutcells::cell{

    void str(const CutCell &cut_cell)
    {
        std::cout << "CutCell: " << std::endl;

        std::cout << "vertex coordinates=[";
        for(auto &i: cut_cell._vertex_coords)
        {
            std::cout << i << ", ";
        }
        std::cout << "]" << std::endl;

        std::cout << "connectivity=[";
        for(int i=0;i<cut_cell._connectivity.size();i++)
        {
            std::cout << i << ": ";
            for(int j=0;j<cut_cell._connectivity[i].size();j++)
            {
                    std::cout << cut_cell._connectivity[i][j]; 
                    if(j<cut_cell._connectivity[i].size()-1) 
                        std::cout << ", ";
            }
            if(i<cut_cell._connectivity.size()-1)
                std::cout << std::endl;
        }
        std::cout << "]" << std::endl;

        std::cout << "cell types=[";
        for(int i=0;i<cut_cell._types.size();i++)
        {
            std::cout << i << ": " << cell_type_to_str(cut_cell._types[i]); 
            if(i<cut_cell._types.size()-1)
                std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }

    /// Cut a cell of type cell_type with the level set values in its vertices and
    /// return a list of elements and their cell types
    void cut(const type cell_type, const std::span<const double> vertex_coordinates, const int gdim,
             const std::span<const double> ls_values, const std::string& cut_type_str,
             CutCell& cut_cell, bool triangulate)
    {
        switch(cell_type)
        {
            case type::interval: interval::cut(vertex_coordinates, gdim, ls_values, cut_type_str, cut_cell);
                                 break;
            case type::triangle: triangle::cut(vertex_coordinates, gdim, ls_values, cut_type_str, cut_cell, triangulate);
                                 break;
            case type::tetrahedron: tetrahedron::cut(vertex_coordinates, gdim, ls_values, cut_type_str, cut_cell, triangulate);
                                 break;
            default: throw std::invalid_argument("Only intervals, triangles and tetrahedra are implemented for cutting so far.");
                                break;
        }
    }

    void cut(const type cell_type, const std::span<const double> vertex_coordinates, const int gdim, 
             const std::span<const double> ls_values, const std::vector<std::string>& cut_type_str,
             std::vector<CutCell>& cut_cell, bool triangulate)
    {
        switch(cell_type)
        {
            case type::triangle: triangle::cut(vertex_coordinates, gdim, ls_values, cut_type_str, cut_cell, triangulate);
                                 break;
        }
    }

    //Cutting of 2nd order triangles (6-node) and tetrahedra (10-node)
    CutCell higher_order_cut(const type cell_type, const std::span<const double> vertex_coordinates, const int gdim, 
             const std::span<const double> ls_values, const std::string& cut_type_str,
             bool triangulate)
    {
      switch(cell_type)
      {
        case cutcells::cell::type::triangle:
        { if(ls_values.size()!=6)
              throw std::invalid_argument("more than 6 level set values, only 2nd order triangles are supported");
          break;
        }
        case cutcells::cell::type::tetrahedron:
        { if(ls_values.size()!=10)
             throw std::invalid_argument("more than 10 level set values, only 2nd order tetrahedra are supported");
          break;
        }
        default: throw std::invalid_argument("Only intervals, triangles and tetrahedra are implemented for cutting so far.");
                break;
      }

      cutcells::cell::domain cell_domain = cutcells::cell::classify_cell_domain(ls_values);

      if(cell_domain == cutcells::cell::domain::intersected)
      {
        std::size_t num_sub_cells = 0;
        switch(cell_type)
        {
          case cutcells::cell::type::triangle: {num_sub_cells = cutcells::cell::triangle_subdivision_table.size();
                                        break;}
          case cutcells::cell::type::tetrahedron: {num_sub_cells = cutcells::cell::tetrahedron_subdivision_table.size();
                                        break;}
        }

        std::vector<cutcells::cell::CutCell> sub_cut_cells;

        int sub_cut_cell_id = 0;

        //Iterate over sub cells
        for(std::size_t i=0;i<num_sub_cells;i++)
        {
          const std::size_t num_vertices = cutcells::cell::get_num_vertices(cell_type);
          std::span<int> sub_tet;

          switch(cell_type)
          {
            case cutcells::cell::type::triangle: {sub_tet = cutcells::cell::triangle_subdivision_table[i];
                                          break;}
            case cutcells::cell::type::tetrahedron: {sub_tet = cutcells::cell::tetrahedron_subdivision_table[i];
                                          break;}
          }

          std::vector<double> sub_ls_values(num_vertices);
          std::vector<double> sub_vertex_coordinates(num_vertices*gdim);

          for(std::size_t j=0;j<num_vertices;j++)
          {
            std::size_t vertex_id = sub_tet[j];
            sub_ls_values[j] = ls_values[vertex_id];

            for(std::size_t k=0;k<gdim;k++)
            {
              sub_vertex_coordinates[j*gdim+k] = vertex_coordinates[vertex_id*gdim+k];
            }
          }

          //Determine if sub_cell is intersected or inside or outside
          cutcells::cell::domain cell_domain = cutcells::cell::classify_cell_domain(sub_ls_values);

          switch(cell_domain)
          {
            case cutcells::cell::domain::intersected:
            {
              cutcells::cell::CutCell tmp_cut_cell;
              cutcells::cell::cut(cell_type, sub_vertex_coordinates, gdim, sub_ls_values, cut_type_str, tmp_cut_cell, triangulate);
              sub_cut_cells.push_back(tmp_cut_cell);
              sub_cut_cells[sub_cut_cell_id]._parent_cell_index.push_back(-1);
              sub_cut_cell_id++;
              break;
            }
            case cutcells::cell::domain::inside:
            {
              if(cut_type_str=="phi<0")
              {
                cutcells::cell::CutCell tmp_cut_cell = cutcells::cell::create_cut_cell(cell_type, sub_vertex_coordinates, gdim);
                sub_cut_cells.push_back(tmp_cut_cell);
                sub_cut_cells[sub_cut_cell_id]._parent_cell_index.push_back(-1);
                sub_cut_cell_id++;
              }
              break;
            }
            case cutcells::cell::domain::outside:
            {
              if(cut_type_str=="phi>0")
              {
                cutcells::cell::CutCell tmp_cut_cell = cutcells::cell::create_cut_cell(cell_type, sub_vertex_coordinates, gdim);
                sub_cut_cells.push_back(tmp_cut_cell);
                sub_cut_cells[sub_cut_cell_id]._parent_cell_index.push_back(-1);
                sub_cut_cell_id++;
              }
              break;
            }
          }
        }
        cutcells::cell::CutCell merged_cut_cell = cutcells::cell::merge(sub_cut_cells);
        return merged_cut_cell;
      }
      else
      {
        throw std::invalid_argument("cell is not intersected and therefore cannot be cut");
      }
    }

    /// Merge vector of CutCell objects into one CutCell
    /// This merging operation is intended to merge several cutcells with the same parent cell
    /// of the same geometrical and topological dimension
    /// this is needed for linearly approximated high order cuts
    CutCell merge(std::vector<CutCell> cut_cell_vec)
    {
      //nothing to merge if only one cell in vector
      if(cut_cell_vec.size()==1)
      {
        return cut_cell_vec[0];
      }

      CutCell merged_cut_cell;
      std::size_t gdim = cut_cell_vec[0]._gdim;
      std::size_t tdim = cut_cell_vec[0]._tdim;

      merged_cut_cell._gdim=gdim;
      merged_cut_cell._tdim=tdim;
      merged_cut_cell._parent_cell_index = cut_cell_vec[0]._parent_cell_index;

      //Count the total number of cells in vector
      int num_cells =0;
      for(auto & cut_cell : cut_cell_vec)
      {
        num_cells += cut_cell._connectivity.size();
      }

      merged_cut_cell._connectivity.resize(num_cells);
      merged_cut_cell._types.resize(num_cells);

      int merged_vertex_id = 0;
      int sub_cell_offset=0;
      int vertex_counter = 0;

      //all cutcells in vector above should have the same gdim and tdim
      for(auto & cut_cell : cut_cell_vec)
      {
        if(cut_cell._vertex_coords.size()==0)
        {
          continue;
        }
        //check that current cut_cell has same dimensions and parent
        if((cut_cell._gdim!=gdim)||(cut_cell._tdim!=tdim)||(cut_cell._parent_cell_index[0]!=merged_cut_cell._parent_cell_index[0]))
        {
          std::cout << "gdim: (" << gdim << ", " << cut_cell._gdim << ")" << std::endl;
          std::cout << "tdim: (" << tdim << ", " << cut_cell._tdim << ")" << std::endl;
          std::cout << "parent_cell: (" << merged_cut_cell._parent_cell_index[0] << ", " << cut_cell._parent_cell_index[0] << ")" << std::endl;
          throw std::runtime_error ("Error in merging cutcells have differing dimensions or parent cell");
        }

        int num_cut_cell_vertices = cut_cell._vertex_coords.size()/gdim;

        int local_num_cells = cut_cell._connectivity.size();

        //Map from vertex id in current cutcell to merged cutcell
        std::map<int,int> local_merged_vertex_ids;

        for(int local_id=0;local_id<num_cut_cell_vertices;local_id++)
        {
          //check if vertex already exists to avoid doubling of vertices
          int id = cutcells::utils::vertex_exists(merged_cut_cell._vertex_coords, cut_cell._vertex_coords, local_id, gdim);

          if(id==-1) //not found
          {
            //add vertex
            merged_vertex_id=vertex_counter;
            vertex_counter++;

            for(int j=0;j<gdim;j++)
            {
              merged_cut_cell._vertex_coords.push_back(cut_cell._vertex_coords[local_id*gdim+j]);
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

        for(int i=0;i<cut_cell._connectivity.size();i++)
        {
          for(int j=0;j<cut_cell._connectivity[i].size();j++)
          {
              int64_t index = cut_cell._connectivity[i][j];
              merged_cut_cell._connectivity[sub_cell_offset+i].push_back(local_merged_vertex_ids[index]);
          }
          merged_cut_cell._types[sub_cell_offset+i]=cut_cell._types[i];
        }

        sub_cell_offset+=local_num_cells;
      }

      return merged_cut_cell;
    }

    //create a cut cell from vertex coordinates and type and gdim
    CutCell create_cut_cell(const type& cell_type, std::span<const double> vertex_coords, const int& gdim)
    {
      CutCell cut_cell;
      cut_cell._vertex_coords.resize(vertex_coords.size());
      for(std::size_t i=0;i<vertex_coords.size();i++)
      {
        cut_cell._vertex_coords[i] = vertex_coords[i];
      }

      cut_cell._gdim = gdim;
      cut_cell._tdim = get_tdim(cell_type);

      cut_cell._connectivity.resize(1);
      int num_vertices = vertex_coords.size()/gdim;
      cut_cell._types.push_back(cell_type);

      cut_cell._connectivity[0].resize(num_vertices);
      for(std::size_t i=0;i<num_vertices;i++)
      {
        cut_cell._connectivity[0][i] = i;
      }

      return cut_cell;
    }

  //cut a cut_cell and return result in cut_cell
  //used for recursive cutting
  void cut_cut_cell(cutcells::cell::CutCell &cut_cell,
                    std::span<const double> ls_vals_all,
                    const int& parent_cell_index,
                    const std::string& cut_type_str,
                    bool triangulate)
  {
    int num_cells = cut_cell._connectivity.size();
    int gdim = cut_cell._gdim;

    std::vector<cutcells::cell::CutCell> cut_cells(num_cells);

    // std::cout << "ls_vals calculated " << ls_vals_all.size() << std::endl;
    // std::cout << "num_cells= " << num_cells << std::endl;
    //std::cout << "level_set_i=" << level_set_i << std::endl;
    int cut_cell_id = 0;

    //iterate over cells in cut_cell and cut each one
    for(std::size_t j=0;j<num_cells;j++)
    {
      int num_vertices = cut_cell._connectivity[j].size();
      auto cut_cell_type = cut_cell._types[j];
      //std::cout << "num_vertices" << num_vertices << std::endl;
      std::vector<double> vertex_coords(num_vertices*gdim);
      std::vector<double> ls_vals(num_vertices);

      for(std::size_t k=0;k<num_vertices;k++)
      {
        int vertex_id = cut_cell._connectivity[j][k];
        for(std::size_t l=0;l<gdim;l++)
        {
          vertex_coords[k*gdim+l] = cut_cell._vertex_coords[vertex_id*gdim+l];
          //std::cout << "vertex_coords=" << vertex_coords[k*gdim+l] << std::endl;
        }

        ls_vals[k] = ls_vals_all[vertex_id];
      }

      // std::cout << "ready to cut cell" << std::endl;
      if(cutcells::cell::is_intersected(ls_vals))
      {
        cutcells::cell::cut(cut_cell_type, vertex_coords, gdim, ls_vals, cut_type_str, cut_cells[cut_cell_id], triangulate);
        cut_cells[cut_cell_id]._parent_cell_index.push_back(parent_cell_index);
        cut_cell_id++;
      }
      else{
        //std::cout << "Warning subcell is not cut!!!!!!!!!!!!!!!!!!!" << std::endl;
        //sub-cell completely inside
        if(cut_type_str=="phi<0")
        {
          cut_cells[cut_cell_id] = create_cut_cell(cut_cell_type, vertex_coords, gdim);
          cut_cells[cut_cell_id]._parent_cell_index.push_back(parent_cell_index);
          cut_cell_id++;
        }
        //sub-cell completely outside
        else
        {
          cut_cells.pop_back();
        }
      }
    }

    //std::cout << "merging cells " << cut_cells.size() << std::endl;
    if(cut_cells.size()>1)
    {
      //Merge cut_cells into one cut cell to cut them again
      cut_cell = cutcells::cell::merge(cut_cells);
      //std::cout << "Cells merged" << std::endl;
    }
    else
    {
      cut_cell = cut_cells[0];
    }
  }

}