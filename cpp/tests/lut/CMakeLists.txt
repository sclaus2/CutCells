# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT

# Tests of the lookup-table backend in cpp/src/lut and of the front end's
# backend "lut"; the meshes and level sets come from the part and quadrays
# test supports.

cutcells_add_test(lut_backend test_lut.cpp)
