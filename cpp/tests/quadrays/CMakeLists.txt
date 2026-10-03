# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT

# Tests of the quadrays backend. support/ holds exact references, test meshes
# and the robustness cases, shared with benchmarks/quadrays.

cutcells_add_test(quadrays_setup test_setup.cpp)
cutcells_add_test(quadrays_exact_reference test_exact_reference.cpp)
cutcells_add_test(quadrays_box_bernstein test_box_bernstein.cpp)
cutcells_add_test(quadrays_plane_exactness test_plane_exactness.cpp)
cutcells_add_test(quadrays_sphere test_sphere.cpp)
cutcells_add_test(quadrays_robustness test_robustness.cpp)
cutcells_add_test(quadrays_pyramid_apex test_pyramid_apex.cpp)
cutcells_add_test(quadrays_leaves test_leaves.cpp)
cutcells_add_test(quadrays_analytic test_analytic.cpp)
cutcells_add_test(quadrays_cells_2d test_cells_2d.cpp)
cutcells_add_test(quadrays_two_roots test_two_roots.cpp)
# Bernstein forms and analytic level sets as two tests each, which run side by side
add_executable(quadrays_prism_pyramid test_prism_pyramid.cpp)
target_link_libraries(quadrays_prism_pyramid PRIVATE cutcells)
add_test(NAME quadrays_prism_pyramid_bernstein COMMAND quadrays_prism_pyramid bernstein)
add_test(NAME quadrays_prism_pyramid_analytic COMMAND quadrays_prism_pyramid analytic)
add_executable(quadrays_level_sets test_level_sets.cpp)
target_link_libraries(quadrays_level_sets PRIVATE cutcells)
add_test(NAME quadrays_level_sets_bernstein COMMAND quadrays_level_sets bernstein)
add_test(NAME quadrays_level_sets_analytic COMMAND quadrays_level_sets analytic)
