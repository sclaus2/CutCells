# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT

# lut: the lookup-table backend (namespace cutcells::lut): cut cells subdivided
# by Pk-iso-P1 templates, each sub-cell cut by the lookup tables with the P1
# interpolants of the level sets, and rules on the straight pieces. Sources
# join the cutcells target; the headers go back to cpp/src/CMakeLists.txt and
# are installed under include/cutcells/lut/.

set(LUT_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_pieces.h
  ${CMAKE_CURRENT_SOURCE_DIR}/piece_rules.h
  PARENT_SCOPE)

target_sources(cutcells PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_pieces.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/piece_rules.cpp
)
