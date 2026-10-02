# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT

# part: the front end without AdaptCell (namespace cutcells::part): cells
# classified by level sets, mesh parts by selection expressions, their
# quadrature and visualisation by backend. Sources join the cutcells target;
# the headers go back to cpp/src/CMakeLists.txt and are installed under
# include/cutcells/part/.

set(PART_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_source.h
  ${CMAKE_CURRENT_SOURCE_DIR}/classify.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_result.h
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh_part.h
  ${CMAKE_CURRENT_SOURCE_DIR}/output.h
  PARENT_SCOPE)

target_sources(cutcells PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_source.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/classify.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_result.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh_part.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/output.cpp
)
