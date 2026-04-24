import pytest
import cutcells


def test_structured_triangle_mesh_view():
    mesh = cutcells.structured_triangle_mesh_view(-1, -1, 1, 1, 5, 4)
    assert mesh.tdim == 2
    assert mesh.gdim == 2
    assert mesh.num_nodes() == 20
    assert mesh.num_cells() == 24


def test_safe_part_name():
    assert cutcells.safe_part_name("phi_left = 0 and phi_right < 0") == "phi_lefteq0andphi_rightlt0"


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
