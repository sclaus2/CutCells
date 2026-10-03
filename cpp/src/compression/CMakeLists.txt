# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT

# compression: positive quadrature rules reduced to at most as many points as
# a polynomial space has moments (namespace cutcells::compression). Sources
# join the cutcells target; the headers go back to cpp/src/CMakeLists.txt and
# are installed under include/cutcells/compression/.

set(COMPRESSION_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/compress.h
  PARENT_SCOPE)

target_sources(cutcells PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/compress.cpp
)
