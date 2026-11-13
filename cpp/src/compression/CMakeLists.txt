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

# The compression kernels for the build machine's CPU: 5 to 45 % faster single
# threaded with AVX2 and FMA (dense Givens, QR and elimination loops). The
# library then runs only on CPUs with the same instructions, so it is off by
# default. Only compress.cpp is compiled so; the cut cells' rounding is unchanged.
if(CUTCELLS_NATIVE_COMPRESSION)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-march=native" CUTCELLS_HAS_MARCH_NATIVE)
  if(CUTCELLS_HAS_MARCH_NATIVE)
    set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/compress.cpp
      TARGET_DIRECTORY cutcells PROPERTIES COMPILE_OPTIONS "-march=native")
  else()
    message(STATUS "CUTCELLS_NATIVE_COMPRESSION: the compiler does not take -march=native")
  endif()
endif()
