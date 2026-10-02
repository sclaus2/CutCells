# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT

# Tests of the front end in cpp/src/part. support/ builds MeshViews of the
# test meshes of cpp/tests/quadrays/support, whose exact references it uses.

cutcells_add_test(part_front_end test_part.cpp)
