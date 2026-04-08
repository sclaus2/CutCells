// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <bit>
#include <numeric>
#include <vector>
#include <stdexcept>

#include "write_vtk.h"
#include "cell_types.h"

namespace
{

constexpr int VTK_LAGRANGE_TRIANGLE = 69;
constexpr int VTK_LAGRANGE_TETRAHEDRON = 71;

constexpr const char* vtk_byte_order()
{
  if constexpr (std::endian::native == std::endian::little)
    return "LittleEndian";
  else
    return "BigEndian";
}

enum class HOCellFamily
{
  triangle,
  tetrahedron
};

int vec_pop(std::vector<int>& v, int i)
{
  const auto pos = (i < 0) ? v.end() + i : v.begin() + i;
  const int value = *pos;
  v.erase(pos);
  return value;
}

int triangle_degree_from_num_nodes(const int num_nodes)
{
  const int n = static_cast<int>((std::sqrt(1 + 8 * num_nodes) - 1) / 2);
  if (2 * num_nodes != n * (n + 1))
    throw std::runtime_error("write_level_set_vtu: invalid triangle node count");
  return n - 1;
}

int tetrahedron_degree_from_num_nodes(const int num_nodes)
{
  int n = 0;
  while (n * (n + 1) * (n + 2) < 6 * num_nodes)
    ++n;
  if (n * (n + 1) * (n + 2) != 6 * num_nodes)
    throw std::runtime_error("write_level_set_vtu: invalid tetrahedron node count");
  return n - 1;
}

std::vector<int> vtk_triangle_remainders(std::vector<int> remainders)
{
  std::vector<int> map;
  map.reserve(remainders.size());

  while (!remainders.empty())
  {
    if (remainders.size() == 1)
    {
      map.push_back(vec_pop(remainders, 0));
      break;
    }

    const int degree = triangle_degree_from_num_nodes(
        static_cast<int>(remainders.size()));

    map.push_back(vec_pop(remainders, 0));
    map.push_back(vec_pop(remainders, degree - 1));
    map.push_back(vec_pop(remainders, -1));

    for (int i = 0; i < degree - 1; ++i)
      map.push_back(vec_pop(remainders, 0));

    for (int i = 1, k = degree * (degree - 1) / 2; i < degree;
         k -= degree - i, ++i)
    {
      map.push_back(vec_pop(remainders, -k));
    }

    for (int i = 1, k = 1; i < degree; k += i, ++i)
      map.push_back(vec_pop(remainders, -k));
  }

  return map;
}

std::vector<int> basix_to_vtk_triangle(const int num_nodes)
{
  std::vector<int> map;
  map.reserve(num_nodes);
  map.insert(map.end(), {0, 1, 2});

  const int degree = triangle_degree_from_num_nodes(num_nodes);
  for (int k = 1; k < degree; ++k)
    map.push_back(3 + 2 * (degree - 1) + k - 1);
  for (int k = 1; k < degree; ++k)
    map.push_back(3 + k - 1);
  for (int k = 1; k < degree; ++k)
    map.push_back(2 * degree - (k - 1));

  if (degree < 3)
    return map;

  std::vector<int> rem(num_nodes - static_cast<int>(map.size()));
  std::iota(rem.begin(), rem.end(), 3 * degree);
  std::vector<int> rem_map = vtk_triangle_remainders(std::move(rem));
  map.insert(map.end(), rem_map.begin(), rem_map.end());
  return map;
}

std::vector<int> vtk_tetrahedron_remainders(std::vector<int> remainders)
{
  std::vector<int> map;
  map.reserve(remainders.size());

  while (!remainders.empty())
  {
    if (remainders.size() == 1)
    {
      map.push_back(vec_pop(remainders, 0));
      break;
    }

    const int deg = tetrahedron_degree_from_num_nodes(
                        static_cast<int>(remainders.size()))
                    + 1;
    map.push_back(vec_pop(remainders, 0));
    map.push_back(vec_pop(remainders, deg - 2));
    map.push_back(vec_pop(remainders, deg * (deg + 1) / 2 - 3));
    map.push_back(vec_pop(remainders, -1));

    if (deg > 2)
    {
      for (int i = 0; i < deg - 2; ++i)
        map.push_back(vec_pop(remainders, 0));

      {
        int d = deg - 2;
        for (int i = 0; i < deg - 2; ++i)
        {
          map.push_back(vec_pop(remainders, d));
          d += (deg - 3 - i);
        }
      }

      {
        int d = (deg - 2) * (deg - 1) / 2 - 1;
        for (int i = 0; i < deg - 2; ++i)
        {
          map.push_back(vec_pop(remainders, d));
          d -= (2 + i);
        }
      }

      {
        int d = (deg - 3) * (deg - 2) / 2;
        for (int i = 0; i < deg - 2; ++i)
        {
          map.push_back(vec_pop(remainders, d));
          d += ((deg - i) * (deg - i - 1) / 2 - 1);
        }
      }

      {
        int d = (deg - 3) * (deg - 2) / 2 + deg - 3;
        for (int i = 0; i < deg - 2; ++i)
        {
          map.push_back(vec_pop(remainders, d));
          d += ((deg - 2 - i) * (deg - 1 - i) / 2 + deg - 4 - i);
        }
      }

      {
        int d = (deg - 3) * (deg - 2) / 2 + deg - 3
              + (deg - 2) * (deg - 1) / 2 - 1;
        for (int i = 0; i < deg - 2; ++i)
        {
          map.push_back(vec_pop(remainders, d));
          d += ((deg - 3 - i) * (deg - 2 - i) / 2 + deg - i - 5);
        }
      }
    }

    if (deg > 3)
    {
      {
        std::vector<int> dofs;
        int d = (deg - 3) * (deg - 2) / 2;
        for (int i = 0; i < deg - 3; ++i)
        {
          for (int ii = 0; ii < deg - 3 - i; ++ii)
            dofs.push_back(vec_pop(remainders, d));
          d += ((deg - 2 - i) * (deg - 1 - i) / 2 - 1);
        }

        std::vector<int> tri_map = vtk_triangle_remainders(std::move(dofs));
        map.insert(map.end(), tri_map.begin(), tri_map.end());
      }

      {
        std::vector<int> dofs;
        int start = deg * deg - 4 * deg + 2;
        int sub_i_start = deg - 3;
        for (int i = 0; i < deg - 3; ++i)
        {
          int d = start;
          int sub_i = sub_i_start;
          for (int ii = 0; ii < deg - 3 - i; ++ii)
          {
            dofs.push_back(vec_pop(remainders, d));
            d += sub_i * (sub_i + 1) / 2 - 2 - i;
            sub_i -= 1;
          }

          start -= (2 + i);
        }

        std::vector<int> tri_map = vtk_triangle_remainders(std::move(dofs));
        map.insert(map.end(), tri_map.begin(), tri_map.end());
      }

      {
        std::vector<int> dofs;
        int start = (deg - 3) * (deg - 2) / 2;
        int sub_i_start = deg - 3;
        for (int i = 0; i < deg - 3; ++i)
        {
          int d = start;
          int sub_i = sub_i_start;
          for (int ii = 0; ii < deg - 3 - i; ++ii)
          {
            dofs.push_back(vec_pop(remainders, d));
            d += sub_i * (sub_i + 1) / 2 - 1 - 2 * i;
            sub_i -= 1;
          }

          start += (deg - 4 - i);
        }

        std::vector<int> tri_map = vtk_triangle_remainders(std::move(dofs));
        map.insert(map.end(), tri_map.begin(), tri_map.end());
      }

      {
        std::vector<int> dofs;
        int add_start = deg - 4;
        for (int i = 0; i < deg - 3; ++i)
        {
          int d = 0;
          int add = add_start;
          for (int ii = 0; ii < deg - 3 - i; ++ii)
          {
            dofs.push_back(vec_pop(remainders, d));
            d += add;
            add -= 1;
          }

          add_start -= 1;
        }

        std::vector<int> tri_map = vtk_triangle_remainders(std::move(dofs));
        map.insert(map.end(), tri_map.begin(), tri_map.end());
      }
    }
  }

  return map;
}

std::vector<int> basix_to_vtk_tetrahedron(const int num_nodes)
{
  const int degree = tetrahedron_degree_from_num_nodes(num_nodes);
  std::vector<int> map;
  map.reserve(num_nodes);
  map.insert(map.end(), {0, 1, 2, 3});

  if (degree < 2)
    return map;

  int base = 4;
  const int edge_dofs = degree - 1;
  for (int edge : {5, 2, 4, 3, 1, 0})
  {
    if (edge == 4)
    {
      for (int i = 0; i < edge_dofs; ++i)
        map.push_back(base + edge_dofs * (edge + 1) - 1 - i);
    }
    else
    {
      for (int i = 0; i < edge_dofs; ++i)
        map.push_back(base + edge_dofs * edge + i);
    }
  }

  if (degree < 3)
    return map;

  base += 6 * edge_dofs;
  const int n_face_dofs = (degree - 1) * (degree - 2) / 2;
  for (int face : {2, 0, 1, 3})
  {
    std::vector<int> face_dofs;
    face_dofs.reserve(n_face_dofs);
    if (face == 2)
    {
      for (int i = 0; i < n_face_dofs; ++i)
        face_dofs.push_back(base + n_face_dofs * face + i);
    }
    else if (face == 0)
    {
      for (int i = degree - 3; i >= 0; --i)
      {
        int d = i;
        for (int ii = 0; ii <= i; ++ii)
        {
          face_dofs.push_back(base + n_face_dofs * face + d);
          d += degree - 3 - ii;
        }
      }
    }
    else
    {
      for (int i = 0; i < degree - 2; ++i)
      {
        int d = i;
        for (int ii = 0; ii < degree - 2 - i; ++ii)
        {
          face_dofs.push_back(base + n_face_dofs * face + d);
          d += degree - 2 - ii;
        }
      }
    }

    std::vector<int> face_map = vtk_triangle_remainders(std::move(face_dofs));
    map.insert(map.end(), face_map.begin(), face_map.end());
  }

  if (degree < 4)
    return map;

  base += 4 * n_face_dofs;
  std::vector<int> remainders((degree - 1) * (degree - 2) * (degree - 3) / 6);
  std::iota(remainders.begin(), remainders.end(), base);
  std::vector<int> rem_map = vtk_tetrahedron_remainders(std::move(remainders));
  map.insert(map.end(), rem_map.begin(), rem_map.end());
  return map;
}

HOCellFamily infer_cell_family(const cutcells::LevelSetMeshData<double, int>& mesh_data,
                               const int cell_id)
{
  if (!mesh_data.cell_types.empty())
  {
    const int vtk_type = mesh_data.cell_types[static_cast<std::size_t>(cell_id)];
    if (vtk_type == static_cast<int>(cutcells::cell::vtk_types::VTK_TRIANGLE)
        || vtk_type == VTK_LAGRANGE_TRIANGLE)
    {
      return HOCellFamily::triangle;
    }
    if (vtk_type == static_cast<int>(cutcells::cell::vtk_types::VTK_TETRA)
        || vtk_type == VTK_LAGRANGE_TETRAHEDRON)
    {
      return HOCellFamily::tetrahedron;
    }

    throw std::runtime_error(
        "write_level_set_vtu: unsupported VTK cell type " + std::to_string(vtk_type)
        + " (supported: triangle, tetrahedron)");
  }

  if (mesh_data.tdim == 2)
    return HOCellFamily::triangle;
  if (mesh_data.tdim == 3)
    return HOCellFamily::tetrahedron;

  throw std::runtime_error(
      "write_level_set_vtu: unsupported tdim without cell_types (supported tdim: 2, 3)");
}

int expected_local_dofs(const HOCellFamily family, const int degree)
{
  if (family == HOCellFamily::triangle)
    return (degree + 1) * (degree + 2) / 2;
  return (degree + 1) * (degree + 2) * (degree + 3) / 6;
}

std::vector<int> basix_to_vtk_lagrange_permutation(const HOCellFamily family,
                                                   const int local_dofs,
                                                   const int degree)
{
  if (degree < 2 || degree > 4)
  {
    throw std::runtime_error(
        "write_level_set_vtu: unsupported polynomial degree " + std::to_string(degree)
        + " (supported: 2, 3, 4)");
  }

  if (local_dofs != expected_local_dofs(family, degree))
  {
    throw std::runtime_error(
        "write_level_set_vtu: local dof count does not match cell family/degree");
  }

  if (family == HOCellFamily::triangle)
    return basix_to_vtk_triangle(local_dofs);
  return basix_to_vtk_tetrahedron(local_dofs);
}

} // namespace

namespace cutcells::io
{
    void write_vtk(std::string filename, const std::span<const double> element_vertex_coords,
                     const std::span<const int> connectivity,
                     const std::span<const int> offsets,
                     const std::span<cell::type> element_types,
                     const int gdim)
    {
        std::ofstream ofs;
        ofs.open(filename.c_str(), std::ios::out );

        if(ofs)
        {
            int num_points = element_vertex_coords.size()/gdim;
            int num_cells = static_cast<int>(offsets.size()) - 1;

            ofs << "<?xml version=\"1.0\"?>\n" 
                << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\""
                << vtk_byte_order() << "\">\n"
                << "\t<UnstructuredGrid>\n"
                << "\t\t<Piece NumberOfPoints=\"" << num_points
                << "\" NumberOfCells=\"" << num_cells << "\">\n";

            ofs << "\t\t\t<Points>\n"
                <<  "\t\t\t  <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"; 

            //coordinates
            for(int i=0;i<num_points;i++)
            {
                double x = element_vertex_coords[i*gdim];
                double y = element_vertex_coords[i*gdim+1];
                double z = 0.0;

                if(gdim==3)
                {
                   z = element_vertex_coords[i*gdim+2];
                }

                ofs << x << " " << y << " " << z << " ";
            }

            ofs <<  "</DataArray>\n";
            ofs << "\t\t\t</Points>\n";
            ofs << "\t\t\t<Cells>\n";
            ofs << "\t\t\t  <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
            for(auto& vertex_id: connectivity)
                ofs << vertex_id << " ";
            ofs <<  "</DataArray>\n"; 
            ofs << "\t\t\t  <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
            for(std::size_t i=1;i<offsets.size();i++)
                ofs << offsets[i] << " ";
            ofs <<  "</DataArray>\n"; 
            ofs << "\t\t\t  <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
            for(auto& type: element_types)
            {
                ofs << static_cast<int>(cell::map_cell_type_to_vtk(type)) << " ";
            }
            ofs <<  "</DataArray>\n";
            ofs << "\t\t\t</Cells>\n";

            ofs << "\t\t</Piece>\n"
                << "\t</UnstructuredGrid>\n"
                << "</VTKFile>\n";

            ofs.close();
        }
        else std:: cout << "Unable to open file " << filename << " to write" << std::endl;
    }

    void write_vtk(std::string filename, cell::CutCell<double>& cut_cell)
    {
        //std::as_const()
        write_vtk(filename, cut_cell._vertex_coords, cut_cell._connectivity, cut_cell._offset,
                     cut_cell._types,
                     cut_cell._gdim);
    }

    void write_level_set_vtu(std::string filename,
                             const cutcells::LevelSetFunction<double>& ls,
                             std::string field_name)
    {
        if (!ls.has_mesh_data())
            throw std::runtime_error("write_level_set_vtu: level set has no mesh_data");
        if (!ls.has_dof_values())
            throw std::runtime_error("write_level_set_vtu: level set has no dof_values");

        const auto& mesh_data = *ls.mesh_data;
        if (mesh_data.gdim < 2 || mesh_data.gdim > 3)
        {
            throw std::runtime_error(
                "write_level_set_vtu: unsupported geometric dimension "
                + std::to_string(mesh_data.gdim) + " (supported: 2 or 3)");
        }

        if (ls.dof_values.size() != static_cast<std::size_t>(mesh_data.num_dofs()))
        {
            throw std::runtime_error(
                "write_level_set_vtu: dof_values size does not match mesh_data.num_dofs()");
        }

        const int num_cells = mesh_data.num_cells();
        std::vector<int> connectivity;
        connectivity.reserve(mesh_data.cell_dofs.size());
        std::vector<int> offsets;
        offsets.reserve(static_cast<std::size_t>(num_cells));
        std::vector<int> types;
        types.reserve(static_cast<std::size_t>(num_cells));

        for (int cell_id = 0; cell_id < num_cells; ++cell_id)
        {
            const std::span<const int> cell_dofs = mesh_data.cell_dofs_span(cell_id);
            const HOCellFamily family = infer_cell_family(mesh_data, cell_id);
            const std::vector<int> perm = basix_to_vtk_lagrange_permutation(
                family,
                static_cast<int>(cell_dofs.size()),
                mesh_data.degree);

            for (const int p : perm)
                connectivity.push_back(cell_dofs[static_cast<std::size_t>(p)]);

            offsets.push_back(static_cast<int>(connectivity.size()));
            types.push_back(
                (family == HOCellFamily::triangle)
                ? VTK_LAGRANGE_TRIANGLE
                : VTK_LAGRANGE_TETRAHEDRON);
        }

        std::ofstream ofs(filename.c_str(), std::ios::out);
        if (!ofs)
            throw std::runtime_error("write_level_set_vtu: unable to open file " + filename);

        const int num_points = mesh_data.num_dofs();
        ofs << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\""
            << vtk_byte_order() << "\">\n"
            << "\t<UnstructuredGrid>\n"
            << "\t\t<Piece NumberOfPoints=\"" << num_points
            << "\" NumberOfCells=\"" << num_cells << "\">\n";

        ofs << "\t\t\t<Points>\n"
            << "\t\t\t  <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
        for (int i = 0; i < num_points; ++i)
        {
            const std::size_t base = static_cast<std::size_t>(i) * static_cast<std::size_t>(mesh_data.gdim);
            const double x = mesh_data.dof_coordinates[base];
            const double y = mesh_data.dof_coordinates[base + 1];
            const double z = (mesh_data.gdim == 3) ? mesh_data.dof_coordinates[base + 2] : 0.0;
            ofs << x << " " << y << " " << z << " ";
        }
        ofs << "</DataArray>\n"
            << "\t\t\t</Points>\n";

        ofs << "\t\t\t<PointData Scalars=\"" << field_name << "\">\n"
            << "\t\t\t  <DataArray type=\"Float64\" Name=\"" << field_name
            << "\" format=\"ascii\">";
        for (const double value : ls.dof_values)
            ofs << value << " ";
        ofs << "</DataArray>\n"
            << "\t\t\t</PointData>\n";

        ofs << "\t\t\t<Cells>\n";
        ofs << "\t\t\t  <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
        for (const int dof : connectivity)
            ofs << dof << " ";
        ofs << "</DataArray>\n";

        ofs << "\t\t\t  <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
        for (const int off : offsets)
            ofs << off << " ";
        ofs << "</DataArray>\n";

        ofs << "\t\t\t  <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
        for (const int type : types)
            ofs << type << " ";
        ofs << "</DataArray>\n";
        ofs << "\t\t\t</Cells>\n";

        ofs << "\t\t</Piece>\n"
            << "\t</UnstructuredGrid>\n"
            << "</VTKFile>\n";
    }

}
