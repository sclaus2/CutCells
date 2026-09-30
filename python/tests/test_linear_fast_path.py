"""The linear LUT fast path for one degree-1 level set on simplices."""

import itertools

import numpy as np
import pytest

import cutcells


def _simplex_box_mesh(n: int, tdim: int):
    """Structured simplex mesh of [-1, 1]^tdim (Kuhn subdivision)."""
    grid = np.linspace(-1.0, 1.0, n + 1)
    if tdim == 2:
        coords = np.array([[x, y] for y in grid for x in grid], dtype=np.float64)

        def node(i, j):
            return j * (n + 1) + i

        cells = []
        for j in range(n):
            for i in range(n):
                cells += [node(i, j), node(i + 1, j), node(i + 1, j + 1)]
                cells += [node(i, j), node(i + 1, j + 1), node(i, j + 1)]
        nv, vtk_type = 3, 5
    else:
        coords = np.array(
            [[x, y, z] for z in grid for y in grid for x in grid], dtype=np.float64
        )

        def node(i, j, k):
            return (k * (n + 1) + j) * (n + 1) + i

        cells = []
        for k, j, i in itertools.product(range(n), repeat=3):
            for perm in itertools.permutations(range(3)):
                c = [i, j, k]
                cells.append(node(*c))
                for axis in perm:
                    c[axis] += 1
                    cells.append(node(*c))
        nv, vtk_type = 4, 10
    connectivity = np.asarray(cells, dtype=np.int32)
    offsets = np.arange(0, connectivity.size + 1, nv, dtype=np.int32)
    cell_types = np.full(connectivity.size // nv, vtk_type, dtype=np.int32)
    return cutcells.MeshView(coords, connectivity, offsets, cell_types, tdim=tdim), coords


def _sphere(X):
    return np.sqrt(np.sum(X * X, axis=0)) - 0.7


@pytest.mark.parametrize("tdim,n", [(2, 12), (3, 6)])
def test_linear_fast_path_matches_generic_path(tdim, n):
    mesh, coords = _simplex_box_mesh(n, tdim)
    # No zero vertex: here the generic path cuts every cell directly as well.
    assert np.min(np.abs(_sphere(coords.T))) > 1.0e-6
    ls = cutcells.create_level_set(mesh, _sphere, degree=1, name="phi")

    fast = cutcells.cut(mesh, ls, triangulate=True)
    generic = cutcells.cut(mesh, ls, triangulate=True, linear_fast_path=False)
    for selector in ("phi < 0", "phi > 0", "phi = 0"):
        q_fast = fast[selector].quadrature(order=2, mode="cut_only")
        q_generic = generic[selector].quadrature(order=2, mode="cut_only")
        assert q_fast.weights.size > 0
        np.testing.assert_array_equal(q_fast.points, q_generic.points)
        np.testing.assert_array_equal(q_fast.weights, q_generic.weights)


@pytest.mark.parametrize(
    "tdim,n,level_set,volume,interface",
    [
        # Zero sets through mesh vertices, edges and faces.
        (3, 4, lambda X: X[0], 4.0, 4.0),
        (3, 4, lambda X: X[0] + X[1] - 0.5, 5.75, 3.0 * np.sqrt(2.0)),
        (2, 8, lambda X: X[0] - 0.5, 3.0, 2.0),
        (2, 8, lambda X: X[1] - X[0], 2.0, 2.0 * np.sqrt(2.0)),
    ],
)
def test_linear_fast_path_is_exact_for_mesh_aligned_zero_sets(
    tdim, n, level_set, volume, interface
):
    mesh, _ = _simplex_box_mesh(n, tdim)
    ls = cutcells.create_level_set(mesh, level_set, degree=1, name="phi")
    result = cutcells.cut(mesh, ls, triangulate=True)

    q_volume = result["phi < 0"].quadrature(order=1, mode="full")
    q_interface = result["phi = 0"].quadrature(order=1, mode="cut_only")
    np.testing.assert_allclose(q_volume.weights.sum(), volume, rtol=0.0, atol=1.0e-12)
    np.testing.assert_allclose(
        q_interface.weights.sum(), interface, rtol=0.0, atol=1.0e-12
    )
