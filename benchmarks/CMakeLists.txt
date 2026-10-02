# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT

# Comparison drivers for CutCells' quadrature backends. This is a CMake
# project of its own: it builds against an installed CutCells, found through
# CMAKE_PREFIX_PATH. Drivers that compare against algoim need
# -DCUTCELLS_WITH_ALGOIM=ON; the library itself never includes algoim.
#
#   cmake -S benchmarks -B build-benchmarks
#         -DCMAKE_PREFIX_PATH="<CutCells prefix>;<conda prefix>"
#         [-DCUTCELLS_WITH_ALGOIM=ON]
#
# With CUTCELLS_WITH_ALGOIM, configure with CONDA_PREFIX set: otherwise CMake
# picks Apple's Accelerate, which lacks the LAPACKE symbols algoim needs.

cmake_minimum_required(VERSION 3.21)
project(cutcells_benchmarks LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(CUTCELLS_WITH_ALGOIM
  "Build the comparisons with algoim (vendored in third_party/algoim)." OFF)

set(CUTCELLS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/..")

# The installed CutCells target may depend on OpenMP::OpenMP_CXX.
find_package(OpenMP QUIET COMPONENTS CXX)
find_package(CutCells REQUIRED CONFIG)
message(STATUS "Found CutCells in: ${CutCells_DIR}")

# Several CutCells installs may be on the search path; refuse one without
# the modules the drivers use.
get_target_property(_cutcells_include_dirs CUTCELLS::cutcells
                    INTERFACE_INCLUDE_DIRECTORIES)
set(_cutcells_has_quadrays OFF)
foreach(_dir IN LISTS _cutcells_include_dirs)
  if(IS_DIRECTORY "${_dir}/cutcells/quadrays")
    set(_cutcells_has_quadrays ON)
  endif()
endforeach()
if(NOT _cutcells_has_quadrays)
  message(FATAL_ERROR
    "The CutCells found in ${CutCells_DIR} has no quadrays module. "
    "Put the prefix of a CutCells install with quadrays first in "
    "-DCMAKE_PREFIX_PATH.")
endif()

# ============================================================================
# algoim (optional): header-only, but its multi-polynomial engine needs
# LAPACKE. Drivers link benchmarks_algoim and are compiled with
# CUTCELLS_WITH_ALGOIM=1.
# ============================================================================
if(CUTCELLS_WITH_ALGOIM)
  set(CUTCELLS_ALGOIM_DIR "${CUTCELLS_SOURCE_DIR}/third_party/algoim/algoim")
  if(NOT EXISTS "${CUTCELLS_ALGOIM_DIR}/quadrature_multipoly.hpp")
    message(FATAL_ERROR
      "CUTCELLS_WITH_ALGOIM=ON requires the vendored algoim headers at "
      "${CUTCELLS_ALGOIM_DIR}")
  endif()

  # Prefer the active conda environment's liblapack/libblas over Accelerate.
  if(NOT DEFINED BLA_VENDOR AND DEFINED ENV{CONDA_PREFIX})
    set(BLA_VENDOR Generic)
  endif()
  find_package(LAPACK REQUIRED)

  add_library(benchmarks_algoim INTERFACE)
  target_include_directories(benchmarks_algoim SYSTEM INTERFACE
    "${CUTCELLS_ALGOIM_DIR}"
    "${CUTCELLS_SOURCE_DIR}/third_party/lapacke_compat")
  target_link_libraries(benchmarks_algoim INTERFACE ${LAPACK_LIBRARIES})
  target_compile_definitions(benchmarks_algoim INTERFACE CUTCELLS_WITH_ALGOIM=1)
endif()

add_subdirectory(quadrays)
