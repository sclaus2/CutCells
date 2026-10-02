# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT

# quadrays: height-function quadrature on cut cells (namespace
# cutcells::quadrays). Sources join the cutcells target with
# target_sources(cutcells PRIVATE ...); the headers go back to
# cpp/src/CMakeLists.txt and are installed under include/cutcells/quadrays/.

set(QUADRAYS_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/rules.h
  PARENT_SCOPE)
