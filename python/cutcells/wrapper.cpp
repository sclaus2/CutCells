// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT


#include <iostream>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <span>

#include <cutcells/cell_flags.h>
#include <cutcells/cell_types.h>
#include <cutcells/cut_cell.h>
#include <cutcells/cut_mesh.h>
#include <cutcells/write_vtk.h>

namespace nb = nanobind;

using namespace cutcells;

namespace
{
const std::string& cell_domain_to_str(cell::domain domain_id)
{
  static const std::map<cell::domain, std::string> type_to_name
      = {{cell::domain::inside, "inside"},
         {cell::domain::intersected, "intersected"},
         {cell::domain::outside, "outside"}};

  auto it = type_to_name.find(domain_id);
  if (it == type_to_name.end())
    throw std::runtime_error("Can't find type");

  return it->second;
}

template <typename V>
auto as_nbarray(V&& x, std::size_t ndim, const std::size_t* shape)
{
  using _V = std::decay_t<V>;
  _V* ptr = new _V(std::move(x));
  return nb::ndarray<typename _V::value_type, nb::numpy>(
      ptr->data(), ndim, shape,
      nb::capsule(ptr, [](void* p) noexcept { delete (_V*)p; }));
}

template <typename V>
auto as_nbarray(V&& x, const std::initializer_list<std::size_t> shape)
{
  return as_nbarray(x, shape.size(), shape.begin());
}

template <typename V>
auto as_nbarray(V&& x)
{
  return as_nbarray(std::move(x), {x.size()});
}

template <typename V, std::size_t U>
auto as_nbarrayp(std::pair<V, std::array<std::size_t, U>>&& x)
{
  return as_nbarray(std::move(x.first), x.second.size(), x.second.data());
}

template <typename T>
void declare_float(nb::module_& m, std::string type)
{
    std::string classify_cell_domain_name = "classify_cell_domain_" + type;
    m.def("classify_cell_domain", [](nb::ndarray<const T>& ls_values){
          cell::domain domain_id = cell::classify_cell_domain<T>(std::span{ls_values.data(),static_cast<unsigned long>(ls_values.size())});
          auto domain_str = cell_domain_to_str(domain_id);
          return domain_str;
        }
        , "classify a cell domain");

    std::string name = "CutCell_" + type;
    nb::class_<cell::CutCell<T>>(m, name.c_str(), "Cut Cell")
        .def(nb::init<>())
        //shape and classes for cutcell vertex coords, connectivity and types are for visualization with pyvista
        .def_prop_ro(
          "vertex_coords",
          [](const cell::CutCell<T>& self) {
            const unsigned long num_vertices = self._vertex_coords.size()/self._gdim;
            std::vector<T> vertex_coords(3*num_vertices);

            int idx = 0;

            for(unsigned long i=0;i<num_vertices;i++)
            {
              T z = 0;

              if(self._gdim==3)
              {
                  z = self._vertex_coords[i*self._gdim+2];
              }
              vertex_coords[idx] = self._vertex_coords[i*self._gdim];
              idx++;
              vertex_coords[idx] = self._vertex_coords[i*self._gdim+1];
              idx++;
              vertex_coords[idx] = z;
              idx++;
            }

            return vertex_coords;
        })
        .def_prop_ro(
          "connectivity",
          [](const cell::CutCell<T>& self) {
            nb::list inner;

            for(std::size_t i=0;i<self._connectivity.size();i++)
            {
              inner.append(self._connectivity[i].size());
              for(std::size_t j=0;j<self._connectivity[i].size();j++)
              {
                inner.append(self._connectivity[i][j]);
              }
            }
            return inner;
          })
        .def_prop_ro(
          "types",
          [](const cell::CutCell<T>& self) {
            //allocate memory
            nb::list types;

              for(std::size_t i=0;i<self._types.size();i++)
              {
                types.append(static_cast<int>(map_cell_type_to_vtk(self._types[i])));
              }

              return types;
          })
        .def("str", [](const cell::CutCell<T>& self) {cell::str(self); return ;})
        .def("volume", [](const cell::CutCell<T>& self) {return cell::volume(self);})
        .def("write_vtk", [](cell::CutCell<double>& self, std::string fname) {io::write_vtk(fname,self); return ;});

  name = "CutMesh_" + type;
  nb::class_<mesh::CutMesh<T>>(m, name.c_str(), "Cut Mesh")
        .def(nb::init<>())
        //shape and classes for cutcell vertex coords, connectivity and types are for visualization with pyvista
        .def_prop_ro(
          "vertex_coords",
          [](const mesh::CutMesh<T>& self) {
            const unsigned long num_vertices = self._vertex_coords.size()/self._gdim;
            std::vector<T> vertex_coords(3*num_vertices);

            int idx = 0;

            for(unsigned long i=0;i<num_vertices;i++)
            {
              T z = 0;

              if(self._gdim==3)
              {
                  z = self._vertex_coords[i*self._gdim+2];
              }
              vertex_coords[idx] = self._vertex_coords[i*self._gdim];
              idx++;
              vertex_coords[idx] = self._vertex_coords[i*self._gdim+1];
              idx++;
              vertex_coords[idx] = z;
              idx++;
            }

            return vertex_coords;
        })
        .def_prop_ro(
          "connectivity",
          [](const mesh::CutMesh<T>& self) {
            nb::list inner;

            for(std::size_t i=0;i<self._connectivity.size();i++)
            {
              std::size_t num_vertices = self._connectivity[i].size();
              inner.append(num_vertices);

              for(std::size_t j=0;j<num_vertices;j++)
              {
                inner.append(self._connectivity[i][j]);
              }
            }
            return inner;
          })
        .def_prop_ro(
          "types",
          [](const mesh::CutMesh<T>& self) {
            //allocate memory
            nb::list types;

              for(std::size_t i=0;i<self._types.size();i++)
              {
                types.append(static_cast<int>(map_cell_type_to_vtk(self._types[i])));
              }

              return types;
          });

  name = "create_cut_mesh_" + type;
  m.def("create_cut_mesh", [](std::vector<cell::CutCell<T>>& cut_cells){
              return mesh::create_cut_mesh(cut_cells);
             }
             , "Creating a cut mesh");

  name = "cut_" + type;
  m.def("cut", [](cell::type cell_type, const nb::ndarray<T>& vertex_coordinates, const int gdim,
             const nb::ndarray<T>& ls_values, const std::string& cut_type_str, bool triangulate){
              cell::CutCell<T> cut_cell;
              cell::cut<T>(cell_type, std::span{vertex_coordinates.data(),static_cast<unsigned long>(vertex_coordinates.size())}, gdim, std::span{ls_values.data(),static_cast<unsigned long>(ls_values.size())}, cut_type_str, cut_cell, triangulate);
              return cut_cell;
             }
             , "cut a cell");

  name = "higher_order_cut_" + type;
  m.def("higher_order_cut", [](cell::type cell_type, const nb::ndarray<T>& vertex_coordinates, const int gdim,
             const nb::ndarray<const T>& ls_values, const std::string& cut_type_str, bool triangulate){
              cell::CutCell<T> cut_cell = cell::higher_order_cut<T>(cell_type, std::span{vertex_coordinates.data(),static_cast<unsigned long>(vertex_coordinates.size())}, gdim, std::span{ls_values.data(),static_cast<unsigned long>(ls_values.size())}, cut_type_str, triangulate);
              return cut_cell;
             }
             , "cut a second order cell");

}

} // namespace

NB_MODULE(_cutcellscpp, m)
{
  // Create module for C++ wrappers
  m.doc() = "CutCells Python interface";

  nb::enum_<cell::type>(m, "CellType")
    .value("point", cell::type::point)
    .value("interval", cell::type::interval)
    .value("triangle", cell::type::triangle)
    .value("tetrahedron", cell::type::tetrahedron)
    .value("quadrilateral", cell::type::quadrilateral)
    .value("hexahedron", cell::type::hexahedron)
    .value("prism", cell::type::prism)
    .value("pyramid", cell::type::pyramid);

  declare_float<float>(m, "float32");
  declare_float<double>(m, "float64");
}
