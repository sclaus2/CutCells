from pathlib import Path

import numpy as np
import pytest

import cutcells


def _single_tetra_mesh():
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
    cell_types = np.array([10], dtype=np.int32)
    return cutcells.MeshView(coords, connectivity, offsets, cell_types, tdim=3)


def _single_triangle_mesh():
    coords = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
        ],
        dtype=np.float64,
    )
    connectivity = np.array([0, 1, 2], dtype=np.int32)
    offsets = np.array([0, 3], dtype=np.int32)
    cell_types = np.array([5], dtype=np.int32)
    return cutcells.MeshView(coords, connectivity, offsets, cell_types, tdim=2)


def _single_quad_mesh():
    coords = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ],
        dtype=np.float64,
    )
    connectivity = np.array([0, 1, 3, 2], dtype=np.int32)
    offsets = np.array([0, 4], dtype=np.int32)
    cell_types = np.array([9], dtype=np.int32)
    return cutcells.MeshView(coords, connectivity, offsets, cell_types, tdim=2)


def _single_hex_mesh():
    coords = np.array(
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
    connectivity = np.array([0, 1, 3, 2, 4, 5, 7, 6], dtype=np.int32)
    offsets = np.array([0, 8], dtype=np.int32)
    cell_types = np.array([12], dtype=np.int32)
    return cutcells.MeshView(coords, connectivity, offsets, cell_types, tdim=3)


def _structured_hex_mesh(n: int):
    grid = np.linspace(-1.0, 1.0, n + 1)
    coords = np.array(
        [[x, y, z] for z in grid for y in grid for x in grid],
        dtype=np.float64,
    )

    def node(i: int, j: int, k: int) -> int:
        return k * (n + 1) * (n + 1) + j * (n + 1) + i

    cells = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                cells.extend(
                    [
                        node(i, j, k),
                        node(i + 1, j, k),
                        node(i + 1, j + 1, k),
                        node(i, j + 1, k),
                        node(i, j, k + 1),
                        node(i + 1, j, k + 1),
                        node(i + 1, j + 1, k + 1),
                        node(i, j + 1, k + 1),
                    ]
                )

    connectivity = np.asarray(cells, dtype=np.int32)
    offsets = np.arange(0, connectivity.size + 1, 8, dtype=np.int32)
    cell_types = np.full(n**3, 12, dtype=np.int32)
    return cutcells.MeshView(coords, connectivity, offsets, cell_types, tdim=3)


def _two_tetra_mesh():
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.2, 0.0, 0.0],
            [0.0, 0.2, 0.0],
            [0.0, 0.0, 0.2],
            [1.0, 0.0, 0.0],
            [1.2, 0.0, 0.0],
            [1.0, 0.2, 0.0],
            [1.0, 0.0, 0.2],
        ],
        dtype=np.float64,
    )
    connectivity = np.array([0, 1, 2, 3, 4, 5, 6, 7], dtype=np.int32)
    offsets = np.array([0, 4, 8], dtype=np.int32)
    cell_types = np.array([10, 10], dtype=np.int32)
    return cutcells.MeshView(coords, connectivity, offsets, cell_types, tdim=3)


def _edge_key(xa: np.ndarray, xb: np.ndarray, digits: int = 12):
    pa = tuple(np.round(np.asarray(xa, dtype=np.float64), digits))
    pb = tuple(np.round(np.asarray(xb, dtype=np.float64), digits))
    return tuple(sorted((pa, pb)))


def _vtk_edges(vtk_type: int):
    if vtk_type == 5:  # triangle
        return ((0, 1), (1, 2), (2, 0))
    if vtk_type == 9:  # quad
        return ((0, 1), (1, 2), (2, 3), (3, 0))
    raise ValueError(f"Unsupported vtk_type {vtk_type} for boundary-edge extraction")


def _edge_counts(cut_mesh):
    points = np.asarray(cut_mesh.vertex_coords, dtype=np.float64)
    cells = np.asarray(cut_mesh.cells, dtype=np.int64)
    vtk_types = np.asarray(cut_mesh.vtk_types, dtype=np.int64)

    counts = {}
    cursor = 0
    for vtk_type in vtk_types:
        num_verts = int(cells[cursor])
        local_ids = cells[cursor + 1 : cursor + 1 + num_verts]
        cursor += num_verts + 1

        for local_a, local_b in _vtk_edges(int(vtk_type)):
            xa = points[int(local_ids[local_a])]
            xb = points[int(local_ids[local_b])]
            key = _edge_key(xa, xb)
            counts[key] = counts.get(key, 0) + 1

    return counts


def _assert_zero_edges_have_only_zero_vertices(adapt_cell):
    zero_masks = np.asarray(adapt_cell.zero_mask_per_vertex, dtype=np.uint64)
    edge_connectivity = np.asarray(adapt_cell.edge_connectivity, dtype=np.int32)
    edge_offsets = np.asarray(adapt_cell.edge_offsets, dtype=np.int32)
    zero_entity_dim = np.asarray(adapt_cell.zero_entity_dim, dtype=np.uint8)
    zero_entity_id = np.asarray(adapt_cell.zero_entity_id, dtype=np.int32)

    for dim, entity_id in zip(zero_entity_dim, zero_entity_id):
        if int(dim) != 1:
            continue
        verts = edge_connectivity[edge_offsets[entity_id] : edge_offsets[entity_id + 1]]
        assert all(zero_masks[vertex_id] & 1 for vertex_id in verts)


def _assert_vertex_provenance_is_populated(adapt_cell, parent_cell_id: int):
    parent_dim = np.asarray(adapt_cell.vertex_parent_dim, dtype=np.int8)
    parent_id = np.asarray(adapt_cell.vertex_parent_id, dtype=np.int32)
    parent_offsets = np.asarray(adapt_cell.vertex_parent_param_offset, dtype=np.int32)
    parent_params = np.asarray(adapt_cell.vertex_parent_param, dtype=np.float64)

    assert len(parent_dim) == adapt_cell.num_vertices()
    assert len(parent_id) == adapt_cell.num_vertices()
    assert len(parent_offsets) == adapt_cell.num_vertices() + 1

    for vertex_id, dim in enumerate(parent_dim):
        dim = int(dim)
        begin = int(parent_offsets[vertex_id])
        end = int(parent_offsets[vertex_id + 1])
        params = parent_params[begin:end]
        assert dim >= 0
        assert int(parent_id[vertex_id]) >= 0
        assert len(params) == dim
        if dim == adapt_cell.tdim:
            assert int(parent_id[vertex_id]) == parent_cell_id
        if dim == 1:
            assert np.all(params >= -1.0e-12)
            assert np.all(params <= 1.0 + 1.0e-12)


def test_triangle_lut_triangulated_quad_connects_uncut_edge_vertices_to_adjacent_roots():
    vertex_coordinates = np.array(
        [
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
        ],
        dtype=np.float64,
    )
    level_set_values = np.array([-0.7, 0.3, 0.3], dtype=np.float64)

    cut = cutcells.cut(
        cutcells.CellType.triangle,
        vertex_coordinates,
        2,
        level_set_values,
        "phi>0",
        True,
        triangulation="midpoint",
    )

    parent_tokens = np.asarray(cut.vertex_parent_entity, dtype=np.int32)
    connectivity = np.asarray(cut.connectivity, dtype=np.int32).reshape(-1, 3)
    token_tris = {tuple(sorted(int(parent_tokens[v]) for v in tri)) for tri in connectivity}

    assert token_tris == {
        tuple(sorted((101, 2, 0))),
        tuple(sorted((0, 2, 1))),
        tuple(sorted((102, 0, 1))),
    }


def test_homeshpart_straight_output_bridge(tmp_path: Path):
    mesh = _single_tetra_mesh()
    ls = cutcells.create_level_set(
        mesh,
        lambda X: X[0] + X[1] - 0.6,
        degree=1,
        name="phi",
    )

    result = cutcells.cut(mesh, ls)
    negative = result["phi < 0"]
    interface = result["phi = 0"]

    vis_cut = negative.visualization_mesh(mode="cut_only")
    vis_full = negative.visualization_mesh(mode="full")
    vis_interface = interface.visualization_mesh(mode="cut_only")

    assert len(np.asarray(vis_cut.types)) > 0
    assert len(np.asarray(vis_full.types)) >= len(np.asarray(vis_cut.types))
    assert len(np.asarray(vis_interface.types)) > 0
    assert np.unique(np.asarray(vis_cut.vtk_types)).tolist() == [13]
    assert np.unique(np.asarray(vis_interface.vtk_types)).tolist() == [9]

    q_cut = negative.quadrature(order=3, mode="cut_only")
    q_full = negative.quadrature(order=3, mode="full")
    q_interface = interface.quadrature(order=3, mode="cut_only")

    assert q_cut.weights.sum() > 0.0
    assert q_full.weights.sum() >= q_cut.weights.sum()
    assert q_interface.weights.sum() > 0.0

    negative_path = tmp_path / "negative_full.vtu"
    interface_path = tmp_path / "interface.vtu"
    negative.write_vtu(str(negative_path), mode="full")
    interface.write_vtu(str(interface_path), mode="cut_only")

    assert negative_path.exists()
    assert interface_path.exists()
    assert "Name=\"types\" format=\"ascii\">13 " in negative_path.read_text()
    assert "Name=\"types\" format=\"ascii\">9 " in interface_path.read_text()


def test_mesh_part_string_selection_supports_volume_or():
    mesh = _two_tetra_mesh()
    left = cutcells.create_level_set(
        mesh,
        lambda X: X[0] - 0.5,
        degree=1,
        name="left",
    )
    right = cutcells.create_level_set(
        mesh,
        lambda X: 0.5 - X[0],
        degree=1,
        name="right",
    )

    result = cutcells.cut(mesh, [left, right])

    left_part = result["left < 0"]
    right_part = result["right < 0"]
    union_part = result["left < 0 or right < 0"]

    np.testing.assert_array_equal(np.asarray(left_part.uncut_cell_ids), np.array([0]))
    np.testing.assert_array_equal(np.asarray(right_part.uncut_cell_ids), np.array([1]))
    np.testing.assert_array_equal(
        np.asarray(union_part.uncut_cell_ids),
        np.array([0, 1]),
    )
    assert union_part.num_cut_cells == 0


def test_mesh_part_string_selection_rejects_mixed_dimension_or():
    mesh = _single_tetra_mesh()
    ls = cutcells.create_level_set(
        mesh,
        lambda X: X[0] + X[1] - 0.6,
        degree=1,
        name="phi",
    )

    result = cutcells.cut(mesh, ls)
    with pytest.raises(RuntimeError, match="same entity dimension"):
        result["phi < 0 or phi = 0"]


def test_mesh_part_string_selection_supports_surface_or():
    mesh = _single_tetra_mesh()
    phi = cutcells.create_level_set(
        mesh,
        lambda X: X[0] + X[1] - 0.6,
        degree=1,
        name="phi",
    )
    psi = cutcells.create_level_set(
        mesh,
        lambda X: X[0] + X[2] - 0.6,
        degree=1,
        name="psi",
    )

    result = cutcells.cut(mesh, [phi, psi])
    union_part = result["phi = 0 or psi = 0"]
    union_mesh = union_part.visualization_mesh(mode="cut_only")

    assert union_part.dim == 2
    np.testing.assert_array_equal(np.asarray(union_part.cut_cell_ids), np.array([0]))
    assert len(np.asarray(union_mesh.types)) > 0


def test_triangulated_tetra_prism_midpoints_keep_masks_and_source_edges():
    mesh = _single_tetra_mesh()
    ls = cutcells.create_level_set(
        mesh,
        lambda X: X[0] + X[1] - 0.6,
        degree=1,
        name="phi",
    )

    result = cutcells.cut(mesh, ls, triangulate=True, triangulation="midpoint")
    adapt_cell = result.adapt_cell(0)
    parent_cell_id = int(np.asarray(result.parent_cell_ids, dtype=np.int32)[0])

    source_edges = np.asarray(adapt_cell.vertex_source_edge_id, dtype=np.int32)
    zero_masks = np.asarray(adapt_cell.zero_mask_per_vertex, dtype=np.uint64)
    negative_masks = np.asarray(adapt_cell.negative_mask_per_vertex, dtype=np.uint64)
    vertex_coords = np.asarray(adapt_cell.vertex_coords, dtype=np.float64)

    # The first four vertices are the original tetra vertices. All vertices
    # inserted by one-root localization or prism triangulation must retain an
    # originating leaf edge so subsequent topology updates can preserve masks.
    assert np.all(source_edges[4:] >= 0)
    _assert_vertex_provenance_is_populated(adapt_cell, parent_cell_id)

    phi = vertex_coords[:, 0] + vertex_coords[:, 1] - 0.6
    for vertex_id in range(vertex_coords.shape[0]):
        is_zero = bool(zero_masks[vertex_id] & 1)
        is_negative = bool(negative_masks[vertex_id] & 1)
        assert not (is_zero and is_negative)
        if not is_zero:
            assert is_negative == bool(phi[vertex_id] < 0.0)

    _assert_zero_edges_have_only_zero_vertices(adapt_cell)


def test_triangulated_triangle_quad_midpoint_keeps_masks_and_parent_edge():
    mesh = _single_triangle_mesh()
    ls = cutcells.create_level_set(
        mesh,
        lambda X: 0.1 - 0.2 * X[0] - 0.1 * X[1],
        degree=1,
        name="phi",
    )

    result = cutcells.cut(mesh, ls, triangulate=True, triangulation="midpoint")
    adapt_cell = result.adapt_cell(0)
    parent_cell_id = int(np.asarray(result.parent_cell_ids, dtype=np.int32)[0])

    source_edges = np.asarray(adapt_cell.vertex_source_edge_id, dtype=np.int32)
    parent_dim = np.asarray(adapt_cell.vertex_parent_dim, dtype=np.int8)
    parent_id = np.asarray(adapt_cell.vertex_parent_id, dtype=np.int32)
    parent_offsets = np.asarray(adapt_cell.vertex_parent_param_offset, dtype=np.int32)
    parent_params = np.asarray(adapt_cell.vertex_parent_param, dtype=np.float64)
    zero_masks = np.asarray(adapt_cell.zero_mask_per_vertex, dtype=np.uint64)
    negative_masks = np.asarray(adapt_cell.negative_mask_per_vertex, dtype=np.uint64)
    vertex_coords = np.asarray(adapt_cell.vertex_coords, dtype=np.float64)

    # Vertices 3 and 4 are one-root vertices; vertex 5 is the triangulation
    # midpoint on the uncut parent edge of the triangle-derived quadrilateral.
    assert np.all(source_edges[3:] >= 0)
    assert int(parent_dim[5]) == 1
    assert int(parent_id[5]) == int(source_edges[5])
    begin = int(parent_offsets[5])
    end = int(parent_offsets[6])
    np.testing.assert_allclose(parent_params[begin:end], [0.5])
    _assert_vertex_provenance_is_populated(adapt_cell, parent_cell_id)

    ls_cell = cutcells.make_cell_level_set(ls, parent_cell_id)
    phi = np.array(
        [
            cutcells.evaluate_bernstein(
                ls_cell.cell_type,
                ls_cell.bernstein_order,
                np.asarray(ls_cell.bernstein_coeffs),
                vertex_coords[vertex_id],
            )
            for vertex_id in range(vertex_coords.shape[0])
        ],
        dtype=np.float64,
    )
    for vertex_id in range(vertex_coords.shape[0]):
        is_zero = bool(zero_masks[vertex_id] & 1)
        is_negative = bool(negative_masks[vertex_id] & 1)
        assert not (is_zero and is_negative)
        if not is_zero:
            assert is_negative == bool(phi[vertex_id] < 0.0)

    _assert_zero_edges_have_only_zero_vertices(adapt_cell)


def test_interface_output_reflects_adaptcell_quad_leaf():
    mesh = _single_tetra_mesh()
    ls = cutcells.create_level_set(
        mesh,
        lambda X: X[0] + X[1] - 0.6,
        degree=1,
        name="phi",
    )

    result = cutcells.cut(mesh, ls)
    interface = result["phi = 0"]

    vis_interface = interface.visualization_mesh(mode="cut_only")

    quad_edges = _edge_counts(vis_interface)

    quad_boundary = {edge for edge, count in quad_edges.items() if count == 1}

    assert np.asarray(vis_interface.vtk_types).tolist() == [9]
    assert len(quad_boundary) == 4


def test_triangle_volume_output_reflects_adaptcell_quad_leaf():
    mesh = _single_triangle_mesh()
    ls = cutcells.create_level_set(
        mesh,
        lambda X: 0.1 - 0.2 * X[0] - 0.1 * X[1],
        degree=1,
        name="phi",
    )

    result = cutcells.cut(mesh, ls)
    negative = result["phi < 0"]

    vis_base = negative.visualization_mesh(mode="cut_only")
    q_base = negative.quadrature(order=3, mode="cut_only")

    assert np.asarray(vis_base.vtk_types).tolist() == [9]
    assert np.asarray(vis_base.vertex_coords).shape[0] == 4
    assert q_base.weights.sum() > 0.0


@pytest.mark.parametrize(
    "mesh_factory,expected_volume_type,expected_interface_type,expected_counts",
    [
        (_single_quad_mesh, 9, 3, {"negative": 2, "positive": 4, "interface": 2}),
        (_single_hex_mesh, 12, 9, {"negative": 4, "positive": 8, "interface": 4}),
    ],
)
def test_iso_p1_tensor_product_quad_hex_output(
    mesh_factory,
    expected_volume_type,
    expected_interface_type,
    expected_counts,
):
    mesh = mesh_factory()
    ls = cutcells.create_level_set(
        mesh,
        lambda X: X[0] - 0.3,
        degree=1,
        name="phi",
    )

    result = cutcells.cut(
        mesh,
        ls,
        triangulate=False,
        cut_approximation="iso_p1",
        cut_approximation_order=2,
    )

    parts = {
        "negative": result["phi < 0"],
        "positive": result["phi > 0"],
        "interface": result["phi = 0"],
    }
    expected = {
        "negative": (expected_volume_type, 0.3),
        "positive": (expected_volume_type, 0.7),
        "interface": (expected_interface_type, 1.0),
    }

    assert result.num_cut_cells == 1
    for name, part in parts.items():
        expected_type, expected_measure = expected[name]
        vis = part.visualization_mesh(mode="cut_only")
        q = part.quadrature(order=3, mode="cut_only")

        assert np.asarray(vis.vtk_types, dtype=np.int32).tolist() == [
            expected_type
        ] * expected_counts[name]
        assert q.tdim == mesh.tdim
        assert np.asarray(q.points).size == np.asarray(q.weights).size * mesh.tdim
        np.testing.assert_allclose(q.weights.sum(), expected_measure, atol=1.0e-12)


def test_hexahedron_diagonal_plane_interface_uses_all_six_cut_vertices():
    mesh = _single_hex_mesh()
    ls = cutcells.create_level_set(
        mesh,
        lambda X: X[0] + X[1] + X[2] - 1.5,
        degree=1,
        name="phi",
    )

    result = cutcells.cut(mesh, ls, triangulate=True)
    interface = result["phi = 0"]

    vis = interface.visualization_mesh(mode="cut_only")
    q = interface.quadrature(order=2, mode="cut_only")
    coords = np.asarray(vis.vertex_coords, dtype=np.float64)
    expected_points = np.array(
        [
            [1.0, 0.5, 0.0],
            [1.0, 0.0, 0.5],
            [0.5, 1.0, 0.0],
            [0.0, 1.0, 0.5],
            [0.5, 0.0, 1.0],
            [0.0, 0.5, 1.0],
        ],
        dtype=np.float64,
    )

    assert result.num_cut_cells == 1
    for point in expected_points:
        assert np.min(np.linalg.norm(coords - point, axis=1)) < 1.0e-12
    np.testing.assert_allclose(
        q.weights.sum(),
        3.0 * np.sqrt(3.0) / 4.0,
        atol=1.0e-12,
    )


def test_hexahedron_sphere_cell_interface_is_not_counted_from_both_sides():
    radius = 0.63
    coords = np.array(
        [
            [0.25, -0.5, -0.5],
            [0.5, -0.5, -0.5],
            [0.25, -0.25, -0.5],
            [0.5, -0.25, -0.5],
            [0.25, -0.5, -0.25],
            [0.5, -0.5, -0.25],
            [0.25, -0.25, -0.25],
            [0.5, -0.25, -0.25],
        ],
        dtype=np.float64,
    )
    connectivity = np.array([0, 1, 3, 2, 4, 5, 7, 6], dtype=np.int32)
    offsets = np.array([0, 8], dtype=np.int32)
    cell_types = np.array([12], dtype=np.int32)
    mesh = cutcells.MeshView(coords, connectivity, offsets, cell_types, tdim=3)
    level_set = lambda X: np.sum(X * X, axis=0) - radius**2

    ls = cutcells.create_level_set(mesh, level_set, degree=1, name="phi")
    result = cutcells.cut(mesh, ls, triangulate=True)
    interface = result["phi = 0"]
    q = interface.quadrature(order=4, mode="cut_only")

    values = level_set(coords.T)
    direct = cutcells.cut(
        cutcells.CellType.hexahedron,
        coords.ravel(),
        3,
        values,
        "phi=0",
        True,
    )
    direct_coords = np.asarray(direct.vertex_coords, dtype=np.float64)
    direct_cells = np.asarray(direct.connectivity, dtype=np.int32).reshape(-1, 3)
    direct_area = 0.0
    for tri in direct_cells:
        a, b, c = direct_coords[tri]
        direct_area += 0.5 * np.linalg.norm(np.cross(b - a, c - a))

    np.testing.assert_allclose(q.weights.sum(), direct_area, atol=1.0e-12)


def test_iso_p1_hexahedron_sphere_exports_generated_zero_faces():
    mesh = _structured_hex_mesh(2)
    radius = 0.8
    ls = cutcells.create_level_set(
        mesh,
        lambda X: X[0] ** 2 + X[1] ** 2 + X[2] ** 2 - radius**2,
        degree=2,
        name="phi",
    )

    result = cutcells.cut(
        mesh,
        ls,
        triangulate=False,
        cut_approximation="iso_p1",
        cut_approximation_order=2,
    )

    interface = result["phi = 0"]
    q = interface.quadrature(order=3, mode="cut_only")
    weights = np.asarray(q.weights, dtype=np.float64)

    assert result.num_cut_cells == 8
    assert np.asarray(q.parent_map).size > result.num_cut_cells
    assert weights.sum() > 1.0
    assert np.asarray(q.points).size == weights.size * mesh.tdim
