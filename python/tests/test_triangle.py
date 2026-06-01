# Copyright (c) 2022 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier:    MIT
import cutcells
import numpy as np
import pytest

TRIANGLE = int(cutcells.CellType.triangle.value)

level_set_values = [(np.array([0.1,-0.1,0.2]),1./12.),
                    (np.array([-0.1,0.1,0.2]),1./12.),
                    (np.array([0.1,0.1,-0.2]),2./9.),
                    (np.array([-0.1,-0.1,0.2]),5./18.),
                    (np.array([0.1,-0.1,-0.2]),5./12.),
                    (np.array([-0.1,0.1,-0.2]),5./12.)]

@pytest.mark.parametrize("ls_values, vol_ex", level_set_values)
def test_triangle_interior(ls_values, vol_ex):
  vertex_coordinates = np.array([0.,0.,1.,0.,1.,1.])

  cell_type = cutcells.CellType.triangle
  triangulate = True
  gdim = 2

  cut_cell = cutcells.cut(cell_type, vertex_coordinates,  gdim, ls_values, "phi<0", triangulate)
  vol = cut_cell.volume()

  assert np.isclose(vol,vol_ex)

level_set_values = [(np.array([0.1,-0.1,0.2]),5./12.),
                    (np.array([-0.1,0.1,0.2]),5./12.),
                    (np.array([0.1,0.1,-0.2]),5./18.),
                    (np.array([-0.1,-0.1,0.2]),2./9.),
                    (np.array([0.1,-0.1,-0.2]),1./12.),
                    (np.array([-0.1,0.1,-0.2]),1./12.)]

@pytest.mark.parametrize("ls_values, vol_ex", level_set_values)
def test_triangle_exterior(ls_values, vol_ex):
  vertex_coordinates = np.array([0.,0.,1.,0.,1.,1.])

  cell_type = cutcells.CellType.triangle
  triangulate = True
  gdim = 2

  cut_cell = cutcells.cut(cell_type, vertex_coordinates,  gdim, ls_values, "phi>0", triangulate)
  vol = cut_cell.volume()

  assert np.isclose(vol,vol_ex)


def test_triangle_quad_midpoint_split():
  vertex_coordinates = np.array([0.,0.,1.,0.,1.,1.], dtype=np.float64)
  ls_values = np.array([0.1,-0.1,-0.2], dtype=np.float64)

  cut_cell = cutcells.cut(
      cutcells.CellType.triangle,
      vertex_coordinates,
      2,
      ls_values,
      "phi<0",
      True,
      triangulation="midpoint",
  )

  coords = np.asarray(cut_cell.vertex_coords, dtype=np.float64)

  assert list(cut_cell.types) == [TRIANGLE, TRIANGLE, TRIANGLE]
  assert coords.shape == (5, 2)
  assert np.any(np.all(np.isclose(coords, np.array([1.0, 0.5])), axis=1))


def test_triangle_quad_classical_triangulation_is_default():
  vertex_coordinates = np.array([0.,0.,1.,0.,1.,1.], dtype=np.float64)
  ls_values = np.array([0.1,-0.1,-0.2], dtype=np.float64)

  default_cut = cutcells.cut(
      cutcells.CellType.triangle,
      vertex_coordinates,
      2,
      ls_values,
      "phi<0",
      True,
  )
  classical_cut = cutcells.cut(
      cutcells.CellType.triangle,
      vertex_coordinates,
      2,
      ls_values,
      "phi<0",
      True,
      triangulation="classical",
  )

  assert list(default_cut.types) == [TRIANGLE, TRIANGLE]
  assert np.asarray(default_cut.vertex_coords, dtype=np.float64).shape == (4, 2)
  np.testing.assert_allclose(default_cut.vertex_coords, classical_cut.vertex_coords)
  np.testing.assert_array_equal(default_cut.connectivity, classical_cut.connectivity)
  np.testing.assert_array_equal(default_cut.types, classical_cut.types)
