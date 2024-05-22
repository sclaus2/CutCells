# Copyright (c) 2022 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier:    MIT
import cutcells
import numpy as np
import pytest

level_set_values = [(np.array([0.1,-0.1,0.2,0.2]),4./27.),
                    (np.array([-0.1,0.1,0.2,0.2]),4./27.),
                    (np.array([0.1,0.1,-0.2,0.2]),16./27.),
                    (np.array([-0.1,-0.1,0.2,0.2]),56./81.),
                    (np.array([0.1,-0.1,-0.2,0.2]),4./3.),
                    (np.array([-0.1,0.1,-0.2,0.2]),4./3.),
                    (np.array([0.1,-0.1,0.2,-0.2]),4./3.),
                    (np.array([-0.1,0.1,0.2,-0.2]),4./3.),
                    (np.array([0.1,0.1,-0.2,-0.2]),160./81.),
                    (np.array([-0.1,-0.1,0.2,-0.2]),56./27.),
                    (np.array([0.1,-0.1,-0.2,-0.2]),68./27.),
                    (np.array([-0.1,0.1,-0.2,-0.2]),68./27.),
                    (np.array([0.1,0.1,0.2,-0.2]),16./27.),
                    (np.array([-0.1,-0.1,-0.2,0.2]),56./27.)]

@pytest.mark.parametrize("ls_values, vol_ex", level_set_values)
def test_tetrahedron_interior(ls_values, vol_ex):
  vertex_coordinates = np.array([1.,1.,1., 1.,-1., -1., -1, 1., -1., -1., -1, 1])

  cell_type = cutcells.CellType.tetrahedron
  triangulate = True
  gdim = 3

  cut_cell = cutcells.cut(cell_type, vertex_coordinates,  gdim, ls_values, "phi<0", triangulate)
  vol = cut_cell.volume()

  assert np.isclose(vol,vol_ex)

level_set_values = [(np.array([0.1,-0.1,0.2,0.2]),68./27.),
                    (np.array([-0.1,0.1,0.2,0.2]),68./27.),
                    (np.array([0.1,0.1,-0.2,0.2]),56./27.),
                    (np.array([-0.1,-0.1,0.2,0.2]),160./81.),
                    (np.array([0.1,-0.1,-0.2,0.2]),4./3.),
                    (np.array([-0.1,0.1,-0.2,0.2]),4./3.),
                    (np.array([0.1,-0.1,0.2,-0.2]),4./3.),
                    (np.array([-0.1,0.1,0.2,-0.2]),4./3.),
                    (np.array([0.1,0.1,-0.2,-0.2]),56./81.),
                    (np.array([-0.1,-0.1,0.2,-0.2]),16./27.),
                    (np.array([0.1,-0.1,-0.2,-0.2]),4./27.),
                    (np.array([-0.1,0.1,-0.2,-0.2]),4./27.),
                    (np.array([0.1,0.1,0.2,-0.2]),56./27.),
                    (np.array([-0.1,-0.1,-0.2,0.2]),16./27.)]

@pytest.mark.parametrize("ls_values, vol_ex", level_set_values)
def test_tetrahedron_exterior(ls_values, vol_ex):
  vertex_coordinates = np.array([1.,1.,1., 1.,-1., -1., -1, 1., -1., -1., -1, 1])

  cell_type = cutcells.CellType.tetrahedron
  triangulate = True
  gdim = 3

  cut_cell = cutcells.cut(cell_type, vertex_coordinates,  gdim, ls_values, "phi>0", triangulate)
  vol = cut_cell.volume()

  assert np.isclose(vol,vol_ex)