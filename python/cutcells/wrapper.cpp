// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT


#include <iostream>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <span>

#include <cutcells/cell_flags.h>
#include <cutcells/cell_types.h>
#include <cutcells/cut_cell.h>
#include <cutcells/write_vtk.h>

namespace py = pybind11;

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

} // namespace


PYBIND11_MODULE(_cutcellscpp, m)
{
  // Create module for C++ wrappers
  m.doc() = "CutCells Python interface";

  m.def("classify_cell_domain", [](const py::array_t<double>& ls_values){
          cell::domain domain_id = cell::classify_cell_domain(std::span{ls_values.data(),static_cast<unsigned long>(ls_values.size())});
          auto domain_str = cell_domain_to_str(domain_id);
          return domain_str;
        }
        , "classify a cell domain");

  py::enum_<cell::type>(m, "CellType")
    .value("point", cell::type::point)
    .value("interval", cell::type::interval)
    .value("triangle", cell::type::triangle)
    .value("tetrahedron", cell::type::tetrahedron)
    .value("quadrilateral", cell::type::quadrilateral)
    .value("hexahedron", cell::type::hexahedron)
    .value("prism", cell::type::prism)
    .value("pyramid", cell::type::pyramid);

  py::class_<cell::CutCell>(m, "CutCell", "Cut Cell")
        .def(py::init<>())
        .def_property_readonly(
          "vertex_coords",
          [](const cell::CutCell& self) {
            int num_points = self._vertex_coords.size()/self._gdim;
            py::array_t<double> vertex_coords(3*num_points);
            py::buffer_info buf1 = vertex_coords.request();
            double *ptr1 = static_cast<double *>(buf1.ptr);

            int idx = 0;

            for(int i=0;i<num_points;i++)
            {
              double x = self._vertex_coords[i*self._gdim];
              double y = self._vertex_coords[i*self._gdim+1];
              double z = 0;
              if(self._gdim==3)
              {
                  z = self._vertex_coords[i*self._gdim+2];
              }
              ptr1[idx] = x;
              idx++;
              ptr1[idx] = y;
              idx++;
              ptr1[idx] = z;
              idx++;
            }

            return vertex_coords; 
          })
        .def_property_readonly(
          "connectivity",
          [](const cell::CutCell& self) {
              int size = 0;
              for(int i=0;i<self._connectivity.size();i++)
              {
                size+= self._connectivity[i].size() + 1;
              }

              py::array_t<int> connectivity(size); 
              py::buffer_info buf1 = connectivity.request();
              int *ptr1 = static_cast<int *>(buf1.ptr);

              int idx = 0;
              for(int i=0;i<self._connectivity.size();i++)
              {
                ptr1[idx] = self._connectivity[i].size();
                idx++;
                for(int j=0;j<self._connectivity[i].size();j++)
                {
                  ptr1[idx] = self._connectivity[i][j];
                  idx++;
                }
              }
              return connectivity;
          })
        .def_property_readonly(
          "types", 
          [](const cell::CutCell& self) {
              py::array_t<int> types(self._types.size()); 
              py::buffer_info buf1 = types.request();
              int *ptr1 = static_cast<int *>(buf1.ptr);

              for(int i=0;i<self._types.size();i++)
              {
                ptr1[i] = static_cast<int>(map_cell_type_to_vtk(self._types[i]));
              }
              return types;
          })
        .def("str", [](const cell::CutCell& self) {cell::str(self); return ;})
        .def("write_vtk", [](cell::CutCell& self, std::string fname) {io::write_vtk(fname,self); return ;});

  m.def("cut", [](cell::type cell_type, const py::array_t<double>& vertex_coordinates, const int gdim, 
             const py::array_t<double>& ls_values, const std::string& cut_type_str, bool triangulate){
              cell::CutCell cut_cell;
              cell::cut(cell_type, std::span{vertex_coordinates.data(),static_cast<unsigned long>(vertex_coordinates.size())}, gdim, std::span{ls_values.data(),static_cast<unsigned long>(ls_values.size())}, cut_type_str, cut_cell, triangulate);
              return cut_cell;
             }
             , "cut a cell");

  m.def("higher_order_cut", [](cell::type cell_type, const py::array_t<double>& vertex_coordinates, const int gdim,
             const py::array_t<const double>& ls_values, const std::string& cut_type_str, bool triangulate){
              cell::CutCell cut_cell = cell::higher_order_cut(cell_type, std::span{vertex_coordinates.data(),static_cast<unsigned long>(vertex_coordinates.size())}, gdim, std::span{ls_values.data(),static_cast<unsigned long>(ls_values.size())}, cut_type_str, triangulate);
              return cut_cell;
             }
             , "cut a second order cell");
}
