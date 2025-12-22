import numpy as np
import cutcells


def test_quadrilateral_vertex_parent_entity_tokens_are_valid():
    # Unit square quad, intersected.
    vertex_coordinates = np.array(
        [0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0], dtype=np.float64
    )
    ls_values = np.array([-0.1, 0.2, -0.2, 0.3], dtype=np.float64)

    cut_cell = cutcells.cut(
        cutcells.CellType.quadrilateral,
        vertex_coordinates,
        2,
        ls_values,
        "phi<0",
        False,
    )

    parents = np.asarray(cut_cell.vertex_parent_entity, dtype=np.int32)
    assert parents.ndim == 1
    coords = np.asarray(cut_cell.vertex_coords)
    assert coords.ndim == 2
    assert parents.size == coords.shape[0]

    # For quads we only expect edge tokens 0..3 and original-vertex tokens 100..103.
    assert np.all(
        ((0 <= parents) & (parents < 4)) | ((100 <= parents) & (parents < 104))
    )


def test_hexahedron_vertex_parent_entity_includes_special_point_token_200():
    # Same case as test_hexahedron_interface_uses_special_point_n0, but assert the token appears.
    ls_values = np.array([-1.0, -1.1, 1.2, 1.3, 1.4, 1.5, -1.6, 1.7], dtype=np.float64)
    vertex_coordinates = np.array(
        [
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            1.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            1.0,
            0.0,
            1.0,
            1.0,
            1.0,
            1.0,
            0.0,
            1.0,
            1.0,
        ],
        dtype=np.float64,
    )

    cut_cell = cutcells.cut(
        cutcells.CellType.hexahedron, vertex_coordinates, 3, ls_values, "phi=0", False
    )
    parents = np.asarray(cut_cell.vertex_parent_entity, dtype=np.int32)
    assert 200 in set(map(int, parents))


def test_triangle_vertex_parent_entity_tokens_are_valid():
    # Unit right-triangle, intersected.
    vertex_coordinates = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
    ls_values = np.array([-0.1, 0.2, -0.3], dtype=np.float64)

    cut_cell = cutcells.cut(
        cutcells.CellType.triangle,
        vertex_coordinates,
        2,
        ls_values,
        "phi<0",
        False,
    )

    parents = np.asarray(cut_cell.vertex_parent_entity, dtype=np.int32)
    assert parents.ndim == 1
    coords = np.asarray(cut_cell.vertex_coords)
    assert coords.ndim == 2
    assert parents.size == coords.shape[0]

    # Triangles have 3 edges and 3 vertices; no special points.
    assert np.all(
        ((0 <= parents) & (parents < 3)) | ((100 <= parents) & (parents < 103))
    )


def test_tetrahedron_vertex_parent_entity_tokens_are_valid():
    # Unit tetrahedron, intersected.
    vertex_coordinates = np.array(
        [
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
        ],
        dtype=np.float64,
    )
    ls_values = np.array([-0.1, 0.2, -0.3, 0.4], dtype=np.float64)

    cut_cell = cutcells.cut(
        cutcells.CellType.tetrahedron,
        vertex_coordinates,
        3,
        ls_values,
        "phi<0",
        False,
    )

    parents = np.asarray(cut_cell.vertex_parent_entity, dtype=np.int32)
    assert parents.ndim == 1
    coords = np.asarray(cut_cell.vertex_coords)
    assert coords.ndim == 2
    assert parents.size == coords.shape[0]

    # Tets have 6 edges and 4 vertices; no special points.
    assert np.all(
        ((0 <= parents) & (parents < 6)) | ((100 <= parents) & (parents < 104))
    )
