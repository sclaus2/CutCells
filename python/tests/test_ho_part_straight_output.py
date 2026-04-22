from pathlib import Path

import numpy as np

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
    vis_cut_triangulated = negative.visualization_mesh(mode="cut_only", triangulate=True)
    vis_interface_triangulated = interface.visualization_mesh(
        mode="cut_only", triangulate=True
    )

    assert len(np.asarray(vis_cut.types)) > 0
    assert len(np.asarray(vis_full.types)) >= len(np.asarray(vis_cut.types))
    assert len(np.asarray(vis_interface.types)) > 0
    assert np.unique(np.asarray(vis_cut.vtk_types)).tolist() == [13]
    assert np.unique(np.asarray(vis_interface.vtk_types)).tolist() == [9]
    assert np.unique(np.asarray(vis_cut_triangulated.vtk_types)).tolist() == [10]
    assert np.unique(np.asarray(vis_interface_triangulated.vtk_types)).tolist() == [5]

    q_cut = negative.quadrature(order=3, mode="cut_only")
    q_full = negative.quadrature(order=3, mode="full")
    q_interface = interface.quadrature(order=3, mode="cut_only")
    q_cut_triangulated = negative.quadrature(
        order=3, mode="cut_only", triangulate=True
    )

    assert q_cut.weights.sum() > 0.0
    assert q_full.weights.sum() >= q_cut.weights.sum()
    assert q_interface.weights.sum() > 0.0
    assert np.isclose(q_cut_triangulated.weights.sum(), q_cut.weights.sum())

    negative_path = tmp_path / "negative_full.vtu"
    interface_path = tmp_path / "interface.vtu"
    negative.write_vtu(str(negative_path), mode="full")
    interface.write_vtu(str(interface_path), mode="cut_only")

    assert negative_path.exists()
    assert interface_path.exists()
    assert "Name=\"types\" format=\"ascii\">13 " in negative_path.read_text()
    assert "Name=\"types\" format=\"ascii\">9 " in interface_path.read_text()


def test_triangulated_interface_preserves_quad_boundary_without_bowtie():
    mesh = _single_tetra_mesh()
    ls = cutcells.create_level_set(
        mesh,
        lambda X: X[0] + X[1] - 0.6,
        degree=1,
        name="phi",
    )

    result = cutcells.cut(mesh, ls)
    interface = result["phi = 0"]

    vis_interface = interface.visualization_mesh(mode="cut_only", triangulate=False)
    vis_interface_triangulated = interface.visualization_mesh(
        mode="cut_only", triangulate=True
    )

    quad_edges = _edge_counts(vis_interface)
    tri_edges = _edge_counts(vis_interface_triangulated)

    quad_boundary = {edge for edge, count in quad_edges.items() if count == 1}
    tri_boundary = {edge for edge, count in tri_edges.items() if count == 1}
    tri_interior = {edge for edge, count in tri_edges.items() if count == 2}

    assert np.asarray(vis_interface.vtk_types).tolist() == [9]
    assert np.asarray(vis_interface_triangulated.vtk_types).tolist() == [5, 5]
    assert tri_boundary == quad_boundary
    assert len(tri_interior) == 1


def test_triangle_volume_output_uses_quad_midpoint_split():
    mesh = _single_triangle_mesh()
    ls = cutcells.create_level_set(
        mesh,
        lambda X: 0.1 - 0.2 * X[0] - 0.1 * X[1],
        degree=1,
        name="phi",
    )

    result = cutcells.cut(mesh, ls)
    negative = result["phi < 0"]

    vis_base = negative.visualization_mesh(mode="cut_only", triangulate=False)
    vis_triangulated = negative.visualization_mesh(mode="cut_only", triangulate=True)
    q_base = negative.quadrature(order=3, mode="cut_only", triangulate=False)
    q_triangulated = negative.quadrature(order=3, mode="cut_only", triangulate=True)

    tri_points = np.asarray(vis_triangulated.vertex_coords, dtype=np.float64)

    assert np.asarray(vis_base.vtk_types).tolist() == [9]
    assert np.asarray(vis_triangulated.vtk_types).tolist() == [5, 5, 5]
    assert np.asarray(vis_base.vertex_coords).shape[0] == 4
    assert tri_points.shape[0] == 5
    assert np.any(np.all(np.isclose(tri_points, np.array([1.0, 0.5])), axis=1))
    assert np.isclose(q_triangulated.weights.sum(), q_base.weights.sum())
