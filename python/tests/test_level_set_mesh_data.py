import numpy as np

import cutcells


def test_create_level_set_mesh_data_triangle_p2_shared_dofs():
    coords = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ],
        dtype=np.float64,
    )
    connectivity = np.array([0, 1, 2, 1, 3, 2], dtype=np.int32)
    offsets = np.array([0, 3, 6], dtype=np.int32)
    cell_types = np.array([5, 5], dtype=np.int32)  # VTK_TRIANGLE

    mesh = cutcells.MeshView(coords, connectivity, offsets, cell_types, 2)
    data = cutcells.create_level_set_mesh_data(mesh, 2)

    assert data.gdim == 2
    assert data.tdim == 2
    assert data.degree == 2
    assert data.num_cells() == 2
    assert data.num_dofs() == 9

    cell_offsets = np.asarray(data.cell_offsets)
    cell_dofs = np.asarray(data.cell_dofs)
    dof_coords = np.asarray(data.dof_coordinates)
    dof_parent_dim = np.asarray(data.dof_parent_dim)
    dof_parent_id = np.asarray(data.dof_parent_id)
    dof_parent_param = np.asarray(data.dof_parent_param)
    dof_parent_param_offset = np.asarray(data.dof_parent_param_offset)

    first_cell = cell_dofs[cell_offsets[0] : cell_offsets[1]]
    second_cell = cell_dofs[cell_offsets[1] : cell_offsets[2]]

    # Basix Triangle6 ordering:
    # vertices (0,1,2), then edge dofs on (1,2), (0,2), (0,1).
    expected_first = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [0.5, 0.5],
            [0.0, 0.5],
            [0.5, 0.0],
        ]
    )
    np.testing.assert_allclose(dof_coords[first_cell], expected_first)

    shared = set(first_cell).intersection(set(second_cell))
    assert len(shared) == 3

    # First 4 dofs are parent vertices; remaining 5 are edge-born.
    np.testing.assert_array_equal(dof_parent_dim[:4], np.zeros(4, dtype=np.int8))
    np.testing.assert_array_equal(dof_parent_id[:4], np.arange(4, dtype=np.int32))
    np.testing.assert_array_equal(dof_parent_dim[4:], np.ones(5, dtype=np.int8))
    assert dof_parent_param_offset.shape[0] == data.num_dofs() + 1

    # The shared midpoint on the diagonal edge should have one edge parameter at t=0.5.
    shared_edge_dof = first_cell[3]
    pb = dof_parent_param_offset[shared_edge_dof]
    pe = dof_parent_param_offset[shared_edge_dof + 1]
    np.testing.assert_allclose(dof_parent_param[pb:pe], np.array([0.5]))


def test_create_level_set_mesh_data_from_arrays():
    dof_coordinates = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ],
        dtype=np.float64,
    )
    cell_dofs = np.array([0, 1, 2], dtype=np.int32)
    cell_offsets = np.array([0, 3], dtype=np.int32)
    cell_types = np.array([5], dtype=np.int32)

    data = cutcells.create_level_set_mesh_data(
        dof_coordinates, cell_dofs, cell_offsets, degree=1, tdim=2, cell_types=cell_types
    )

    assert data.gdim == 2
    assert data.tdim == 2
    assert data.degree == 1
    np.testing.assert_array_equal(np.asarray(data.cell_dofs), cell_dofs)
    np.testing.assert_array_equal(np.asarray(data.cell_offsets), cell_offsets)
    np.testing.assert_allclose(np.asarray(data.dof_coordinates), dof_coordinates)
    assert np.asarray(data.dof_parent_dim).size == 0


def test_make_cell_level_set_matches_generated_mesh_data_without_provenance():
    coords = np.array(
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
    mesh = cutcells.MeshView(coords, connectivity, offsets, cell_types, 3)

    ls = cutcells.create_level_set(
        mesh,
        lambda X: X[0] * X[0] + 0.5 * X[1] * X[2] - 0.25 * X[2] + 0.1,
        degree=2,
        name="phi",
    )
    ls_cell_generated = cutcells.make_cell_level_set(ls, 0)

    mesh_data = ls.mesh_data
    mesh_data_from_arrays = cutcells.create_level_set_mesh_data(
        np.array(mesh_data.dof_coordinates),
        np.array(mesh_data.cell_dofs),
        np.array(mesh_data.cell_offsets),
        degree=mesh_data.degree,
        tdim=mesh_data.tdim,
        cell_types=np.array([10], dtype=np.int32),
    )
    ls_from_arrays = cutcells.create_level_set_function(
        mesh_data_from_arrays,
        np.array(ls.dof_values),
        name="phi",
    )
    ls_cell_from_arrays = cutcells.make_cell_level_set(ls_from_arrays, 0)

    np.testing.assert_allclose(
        np.asarray(ls_cell_from_arrays.bernstein_coeffs),
        np.asarray(ls_cell_generated.bernstein_coeffs),
    )
