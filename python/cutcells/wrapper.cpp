// Copyright (C) 2017 Chris N. Richardson and Garth N. Wells
//
// This file is part of DOLFINx (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include <iostream>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <span>

#include <cutcells/cell_flags.h>
#include <cutcells/cut_cell.h>

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
          cell::domain domain_id = cell::classify_cell_domain(std::span{ls_values.data(),ls_values.size()});
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
            return py::array_t<double>(self._vertex_coords.size(), self._vertex_coords.data(), py::cast(self));
          })
        .def("str", [](const cell::CutCell& self) {cell::str(self); return ;});

//FIXME: does this copy the cut cell? 
  m.def("cut", [](cell::type cell_type, const py::array_t<double>& vertex_coordinates, const int gdim, 
             const py::array_t<double>& ls_values, const std::string& cut_type_str, bool triangulate){
              cell::CutCell cut_cell;
              cell::cut(cell_type, std::span{vertex_coordinates.data(),vertex_coordinates.size()}, gdim, std::span{ls_values.data(),ls_values.size()}, cut_type_str, cut_cell);
              return cut_cell;
             }
             , "cut a cell");
}
