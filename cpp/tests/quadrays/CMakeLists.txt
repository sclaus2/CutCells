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
cutcells_add_test(quadrays_leaves test_leaves.cpp)
cutcells_add_test(quadrays_analytic test_analytic.cpp)
