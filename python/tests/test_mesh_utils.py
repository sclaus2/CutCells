import pytest
import cutcells


def test_rectangle_triangle():
    grid = cutcells.rectangle_triangle_mesh(-1, -1, 1, 1, 10, 10)
    assert grid.n_points == 100
    assert grid.n_cells > 0
    assert int(grid.celltypes[0]) == 5  # VTK_TRIANGLE


def test_mesh_from_pyvista_infers_tdim():
    grid = cutcells.rectangle_triangle_mesh(-1, -1, 1, 1, 5, 5)
    mesh = cutcells.mesh_from_pyvista(grid)
    assert mesh.tdim == 2
    assert mesh.gdim == 3
    assert mesh.num_cells() == grid.n_cells
