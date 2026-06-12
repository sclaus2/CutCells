import numpy as np

import cutcells


def test_direct_quadrilateral_cut_uses_basix_ordering():
    points = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
        dtype=np.float64,
    )
    ls_values = points[:, 0] - 0.3

    cut_cell = cutcells.cut(
        cutcells.CellType.quadrilateral,
        points.ravel(),
        2,
        ls_values,
        "phi<0",
        False,
        triangulation="none",
    )

    coords = np.asarray(cut_cell.vertex_coords)
    conn = np.asarray(cut_cell.connectivity, dtype=np.int32)
    offsets = np.asarray(cut_cell.offsets, dtype=np.int32)

    assert offsets.tolist() == [0, 4]
    cell_coords = coords[conn[offsets[0] : offsets[1]]]
    np.testing.assert_allclose(
        cell_coords,
        np.array([[0.0, 0.0], [0.3, 0.0], [0.0, 1.0], [0.3, 1.0]]),
        atol=1.0e-14,
    )

    tokens = np.asarray(cut_cell.vertex_parent_entity, dtype=np.int32)
    assert tokens[conn].tolist() == [100, 0, 102, 3]


def test_direct_hexahedron_cut_uses_basix_ordering():
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ],
        dtype=np.float64,
    )
    ls_values = points[:, 0] - 0.3

    cut_cell = cutcells.cut(
        cutcells.CellType.hexahedron,
        points.ravel(),
        3,
        ls_values,
        "phi<0",
        False,
        triangulation="none",
    )

    coords = np.asarray(cut_cell.vertex_coords)
    conn = np.asarray(cut_cell.connectivity, dtype=np.int32)
    offsets = np.asarray(cut_cell.offsets, dtype=np.int32)

    assert offsets.tolist() == [0, 8]
    cell_coords = coords[conn[offsets[0] : offsets[1]]]
    jac = np.column_stack(
        [cell_coords[1] - cell_coords[0],
         cell_coords[2] - cell_coords[0],
         cell_coords[4] - cell_coords[0]]
    )
    assert np.isclose(abs(np.linalg.det(jac)), 0.3, atol=1.0e-14)
    assert np.isclose(cell_coords[:, 0].min(), 0.0)
    assert np.isclose(cell_coords[:, 0].max(), 0.3)

    tokens = np.asarray(cut_cell.vertex_parent_entity, dtype=np.int32)
    assert set(tokens.tolist()) == {100, 102, 104, 106, 0, 5, 8, 11}
