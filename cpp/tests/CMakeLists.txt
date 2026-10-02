# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT

# C++ tests: plain executables that return non-zero on failure, run by ctest.
# Configure with -DBUILD_TESTING=ON.

# cutcells_add_test(<name> <source>...): an executable linked to the cutcells
# library, registered with ctest under <name>.
function(cutcells_add_test name)
  add_executable(${name} ${ARGN})
  target_link_libraries(${name} PRIVATE cutcells)
  add_test(NAME ${name} COMMAND ${name})
endfunction()

add_subdirectory(quadrays)
add_subdirectory(lut)
add_subdirectory(part)
