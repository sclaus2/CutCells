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
