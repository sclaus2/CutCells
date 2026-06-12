import math

import numpy as np
import pytest

import cutcells


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


def _single_embedded_interval_mesh():
    coords = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
        ],
        dtype=np.float64,
    )
    connectivity = np.array([0, 1], dtype=np.int32)
    offsets = np.array([0, 2], dtype=np.int32)
    cell_types = np.array([3], dtype=np.int32)
    return cutcells.MeshView(coords, connectivity, offsets, cell_types, tdim=1)


def _single_embedded_quad_mesh():
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 0.0, 3.0],
            [2.0, 0.0, 3.0],
        ],
        dtype=np.float64,
    )
    connectivity = np.array([0, 1, 3, 2], dtype=np.int32)
    offsets = np.array([0, 4], dtype=np.int32)
    cell_types = np.array([9], dtype=np.int32)
    return cutcells.MeshView(coords, connectivity, offsets, cell_types, tdim=2)


def _quadrature_or_skip(part, *, order: int, mode: str, backend: str = "algoim"):
    try:
        return part.quadrature(order=order, mode=mode, backend=backend)
    except RuntimeError as exc:
        if "without Algoim" in str(exc):
            pytest.skip("CutCells was built without Algoim support")
        raise


def test_algoim_quadrature_single_quad_linear_cut():
    mesh = _single_quad_mesh()
    level_set = cutcells.create_level_set(
        mesh, lambda x: x[0] + x[1] - 0.8, degree=1, name="phi"
    )
    result = cutcells.cut(mesh, level_set)

    volume = _quadrature_or_skip(
        result["phi < 0"], order=4, mode="cut_only"
    )
    interface = _quadrature_or_skip(
        result["phi = 0"], order=4, mode="cut_only"
    )

    assert np.sum(volume.weights) == pytest.approx(0.5 * 0.8**2, abs=1.0e-12)
    assert np.sum(interface.weights) == pytest.approx(
        math.sqrt(2.0) * 0.8, abs=1.0e-12
    )


def test_algoim_general_backend_single_quad_linear_cut():
    mesh = _single_quad_mesh()
    level_set = cutcells.create_level_set(
        mesh, lambda x: x[0] + x[1] - 0.8, degree=1, name="phi"
    )
    result = cutcells.cut(mesh, level_set)

    volume = _quadrature_or_skip(
        result["phi < 0"], order=4, mode="cut_only", backend="algoim_general"
    )

    assert np.sum(volume.weights) == pytest.approx(0.5 * 0.8**2, abs=1.0e-12)


def test_algoim_quadrature_single_hex_linear_cut():
    mesh = _single_hex_mesh()
    level_set = cutcells.create_level_set(
        mesh, lambda x: x[0] + x[1] + x[2] - 0.8, degree=1, name="phi"
    )
    result = cutcells.cut(mesh, level_set)

    volume = _quadrature_or_skip(
        result["phi < 0"], order=4, mode="cut_only"
    )
    interface = _quadrature_or_skip(
        result["phi = 0"], order=4, mode="cut_only"
    )

    assert np.sum(volume.weights) == pytest.approx(0.8**3 / 6.0, abs=1.0e-12)
    assert np.sum(interface.weights) == pytest.approx(
        math.sqrt(3.0) * 0.8**2 / 2.0, abs=1.0e-12
    )


@pytest.mark.parametrize("backend", ["algoim", "algoim_general"])
def test_algoim_interval_interface_uses_root_point_quadrature(backend):
    mesh = _single_embedded_interval_mesh()
    mesh_data = cutcells.create_level_set_mesh_data(
        np.array([[0.0, 0.0], [1.0, 0.0]], dtype=np.float64),
        np.array([0, 1], dtype=np.int32),
        np.array([0, 2], dtype=np.int32),
        degree=1,
        tdim=1,
        cell_types=np.array([3], dtype=np.int32),
    )
    level_set = cutcells.create_level_set_function(
        mesh_data, np.array([-0.3, 0.7], dtype=np.float64), name="phi"
    )
    result = cutcells.cut(mesh, level_set)

    interface = _quadrature_or_skip(
        result["phi = 0"], order=4, mode="cut_only", backend=backend
    )

    points = np.asarray(interface.points, dtype=np.float64).reshape(-1, mesh.tdim)
    assert points.shape == (1, 1)
    assert np.sum(interface.weights) == pytest.approx(1.0, abs=1.0e-14)
    assert points[0, 0] == pytest.approx(0.3, abs=1.0e-12)


@pytest.mark.parametrize("backend", ["algoim", "algoim_general"])
def test_algoim_embedded_quad_interface_scales_with_induced_metric(backend):
    mesh = _single_embedded_quad_mesh()
    level_set = cutcells.create_level_set(
        mesh, lambda x: x[0] - 0.8, degree=1, name="phi"
    )
    result = cutcells.cut(mesh, level_set)

    interface = _quadrature_or_skip(
        result["phi = 0"], order=4, mode="cut_only", backend=backend
    )

    points = np.asarray(interface.points, dtype=np.float64).reshape(-1, mesh.tdim)
    assert points.shape[0] > 0
    assert np.sum(interface.weights) == pytest.approx(3.0, abs=1.0e-12)
    np.testing.assert_allclose(points[:, 0], 0.4, rtol=1.0e-12, atol=1.0e-12)
