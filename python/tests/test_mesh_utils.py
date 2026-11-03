import numpy as np
import pytest
import cutcells


def test_rectangle_mesh():
    mesh = cutcells.rectangle_mesh(-1, -1, 1, 1, 5, 4, cell_type="triangle")
    assert mesh.tdim == 2
    assert mesh.gdim == 2
    assert mesh.num_nodes() == 20
    assert mesh.num_cells() == 24


def test_safe_part_name():
    assert cutcells.safe_part_name("phi_left = 0 and phi_right < 0") == "phi_lefteq0andphi_rightlt0"


def test_rectangle_triangle():
    pytest.importorskip("pyvista")
    grid = cutcells.rectangle_triangle_mesh(-1, -1, 1, 1, 10, 10)
    assert grid.n_points == 100
    assert grid.n_cells > 0
    assert int(grid.celltypes[0]) == 5  # VTK_TRIANGLE


def test_mesh_from_pyvista_infers_tdim():
    pytest.importorskip("pyvista")
    grid = cutcells.rectangle_triangle_mesh(-1, -1, 1, 1, 5, 5)
    mesh = cutcells.mesh_from_pyvista(grid)
    assert mesh.tdim == 2
    assert mesh.gdim == 3
    assert mesh.num_cells() == grid.n_cells


@pytest.mark.parametrize(
    "make, num_cells, volume",
    [
        (lambda: cutcells.rectangle_mesh(-1, -1, 1, 1, 5, 4, cell_type="quadrilateral"), 12, 4.0),
        (lambda: cutcells.box_mesh(0, 0, 0, 1, 2, 1, 3, 4, 2, cell_type="tetrahedron"), 36, 2.0),
        (lambda: cutcells.box_mesh(0, 0, 0, 1, 2, 1, 3, 4, 2, cell_type="hexahedron"), 6, 2.0),
    ],
)
def test_structured_mesh_views(make, num_cells, volume):
    mesh = make()
    assert mesh.num_cells() == num_cells
    assert mesh.gdim == mesh.tdim
    # A plane away from the box leaves every cell inside: the box's volume.
    plane = cutcells.create_level_set(mesh, lambda X: X[0] - 5.0, degree=1, name="plane")
    rules = cutcells.cut(mesh, plane)[plane < 0].quadrature(order=1)
    assert rules.weights.sum() == pytest.approx(volume)
    assert (rules.weights > 0).all()


def test_mesh_builder_names():
    with pytest.raises(ValueError):
        cutcells.rectangle_mesh(0, 0, 1, 1, 3, 3, cell_type="hexahedron")
    with pytest.warns(DeprecationWarning, match="rectangle_mesh"):
        old = cutcells.structured_triangle_mesh_view(0, 0, 1, 1, 3, 3)
    new = cutcells.rectangle_mesh(0, 0, 1, 1, 3, 3)
    assert np.array_equal(old.connectivity, new.connectivity)
    assert np.array_equal(old.coordinates, new.coordinates)


def test_tetrahedra_are_conforming():
    mesh = cutcells.box_mesh(0, 0, 0, 1, 1, 1, 3, 3, 3, cell_type="tetrahedron")
    tets = mesh.connectivity.reshape(-1, 4)
    x = mesh.coordinates[tets]
    det = np.linalg.det(x[:, 1:] - x[:, :1])
    assert (det > 0).all()
    # Every interior face is shared by exactly two tetrahedra.
    faces = np.sort(tets[:, [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]].reshape(-1, 3), axis=1)
    _, counts = np.unique(faces, axis=0, return_counts=True)
    assert set(counts.tolist()) <= {1, 2}
    boundary = (counts == 1).sum()
    assert boundary == 6 * 2 * 2 * 2  # two triangles per boundary grid square


def test_to_pyvista():
    pytest.importorskip("pyvista")
    mesh = cutcells.rectangle_mesh(-1, -1, 1, 1, 9, 9, cell_type="triangle")
    ball = cutcells.analytic_sphere([0.0, 0.0, 0.0], 0.6, name="ball")
    result = cutcells.cut(mesh, ball)
    grid = cutcells.to_pyvista(mesh)
    assert grid.n_cells == mesh.num_cells()
    assert int(grid.celltypes[0]) == 5
    curve = cutcells.to_pyvista(result[ball == 0].visualization_mesh(mode="cut_only", degree=2))
    assert curve.n_cells > 0
    assert "parent" in curve.cell_data
