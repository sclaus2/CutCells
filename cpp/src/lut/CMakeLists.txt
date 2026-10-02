# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT

# lut: the lookup tables. cut_<cell> cut one cell by the linear interpolant
# of a level set (tables in generated/), cut_cell and cut_mesh cut cells and
# meshes, triangulation and the midpoint splits divide the pieces into
# simplices, and iso_refine holds the Pk-iso-P1 templates. cell_pieces and
# piece_rules are the lookup-table backend of the front end (namespace
# cutcells::lut): cut cells subdivided by the templates, each sub-cell cut by
# the tables, and rules on the straight pieces. Sources join the cutcells
# target; the headers go back to cpp/src/CMakeLists.txt and are installed
# under include/cutcells/lut/ (the tables under lut/generated/).

set(LUT_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_interval.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_triangle.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_quadrilateral.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_tetrahedron.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_hexahedron.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_prism.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_pyramid.h
  ${CMAKE_CURRENT_SOURCE_DIR}/quad_midpoint_split.h
  ${CMAKE_CURRENT_SOURCE_DIR}/prism_midpoint_split.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_subdivision.h
  ${CMAKE_CURRENT_SOURCE_DIR}/triangulation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_cell.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_mesh.h
  ${CMAKE_CURRENT_SOURCE_DIR}/iso_refine.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_pieces.h
  ${CMAKE_CURRENT_SOURCE_DIR}/piece_rules.h
  PARENT_SCOPE)

target_sources(cutcells PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_interval.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_triangle.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_quadrilateral.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_tetrahedron.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_hexahedron.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_prism.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_pyramid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/quad_midpoint_split.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/prism_midpoint_split.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_cell.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_mesh.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/iso_refine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_pieces.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/piece_rules.cpp
)
