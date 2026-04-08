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
