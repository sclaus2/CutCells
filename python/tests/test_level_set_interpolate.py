import numpy as np
import pytest

import cutcells


def _make_mesh(cell_family: str):
    if cell_family == "triangle":
        coordinates = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
            dtype=np.float64,
        )
        connectivity = np.array([0, 1, 2], dtype=np.int32)
        offsets = np.array([0, 3], dtype=np.int32)
        cell_types = np.array([5], dtype=np.int32)  # VTK_TRIANGLE
        tdim = 2
    elif cell_family == "tetra":
        coordinates = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=np.float64,
        )
        connectivity = np.array([0, 1, 2, 3], dtype=np.int32)
        offsets = np.array([0, 4], dtype=np.int32)
        cell_types = np.array([10], dtype=np.int32)  # VTK_TETRA
        tdim = 3
    else:
        raise ValueError(f"Unsupported cell_family={cell_family}")

    return cutcells.MeshView(
        coordinates=coordinates,
        connectivity=connectivity,
        offsets=offsets,
        cell_types=cell_types,
        tdim=tdim,
    )


@pytest.mark.parametrize("cell_family", ["triangle", "tetra"])
@pytest.mark.parametrize("degree", [2, 3, 4])
def test_interpolate_level_set_calls_callback_once_and_sets_discrete_data(cell_family, degree):
    mesh = _make_mesh(cell_family)
    callback_count = 0
    callback_shapes = []

    def phi(X):
        nonlocal callback_count
        callback_count += 1
        callback_shapes.append(X.shape)
        return X[:, 0] + 2.0 * X[:, 1] - 0.5 * X[:, 2]

    ls = cutcells.interpolate_level_set(mesh, phi, degree=degree)

    assert callback_count == 1
    assert callback_shapes[0] == (ls.mesh_data.num_dofs(), mesh.gdim)

    assert ls.has_mesh_data() is True
    assert ls.has_dof_values() is True
    assert ls.mesh_data.degree == degree
    assert len(ls.dof_values) == ls.mesh_data.num_dofs()

    coords = np.asarray(ls.mesh_data.dof_coordinates)
    np.testing.assert_allclose(np.asarray(ls.dof_values), phi(coords))


def test_create_level_set_function_from_mesh_data_and_values():
    dof_coordinates = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.0, 0.5, 0.0],
            [0.5, 0.0, 0.0],
        ],
        dtype=np.float64,
    )
    cell_dofs = np.array([0, 1, 2, 3, 4, 5], dtype=np.int32)
    cell_offsets = np.array([0, 6], dtype=np.int32)
    cell_types = np.array([5], dtype=np.int32)
    dof_values = np.linspace(-1.0, 1.0, 6, dtype=np.float64)

    mesh_data = cutcells.create_level_set_mesh_data(
        dof_coordinates=dof_coordinates,
        cell_dofs=cell_dofs,
        cell_offsets=cell_offsets,
        degree=2,
        tdim=2,
        cell_types=cell_types,
    )
    ls = cutcells.create_level_set_function(mesh_data, dof_values)

    assert ls.has_mesh_data() is True
    assert ls.has_dof_values() is True
    np.testing.assert_allclose(np.asarray(ls.dof_values), dof_values)


def test_legacy_level_set_constructor_paths_still_work():
    ls_callable = cutcells.LevelSetFunction(
        value=lambda x: float(x[0] - 0.25),
        gdim=3,
    )
    assert ls_callable.has_value() is True
    assert ls_callable.has_nodal_values() is False
    assert ls_callable.has_mesh_data() is False
    assert ls_callable.has_dof_values() is False
    value = ls_callable.value(np.array([0.5, 0.0, 0.0], dtype=np.float64))
    assert value == pytest.approx(0.25)

    nodal = np.array([0.1, -0.2, 0.3], dtype=np.float64)
    ls_nodal = cutcells.LevelSetFunction(nodal_values=nodal, gdim=3)
    assert ls_nodal.has_nodal_values() is True
    assert ls_nodal.value_at_node(1) == pytest.approx(-0.2)
    assert ls_nodal.has_mesh_data() is False
    assert ls_nodal.has_dof_values() is False
