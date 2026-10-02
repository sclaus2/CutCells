# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT

# quadrays: height-function quadrature on cut cells (namespace
# cutcells::quadrays). Sources join the cutcells target; the headers go back
# to cpp/src/CMakeLists.txt and are installed under include/cutcells/quadrays/.

set(QUADRAYS_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/analytic.h
  ${CMAKE_CURRENT_SOURCE_DIR}/box_bernstein.h
  ${CMAKE_CURRENT_SOURCE_DIR}/clipped_box.h
  ${CMAKE_CURRENT_SOURCE_DIR}/engine.h
  ${CMAKE_CURRENT_SOURCE_DIR}/leaves.h
  ${CMAKE_CURRENT_SOURCE_DIR}/rules.h
  ${CMAKE_CURRENT_SOURCE_DIR}/source.h
  ${CMAKE_CURRENT_SOURCE_DIR}/taylor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/adapters/shapeforest_tape.h
  PARENT_SCOPE)

target_sources(cutcells PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/analytic.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/box_bernstein.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/clipped_box.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/leaves.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/rules.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/source.cpp
)
