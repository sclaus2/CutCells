# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier:    MIT
"""Linear level sets whose zero set contains mesh vertices, edges or facets:
the parts are exact, cells only touching the zero set take no time, and an
interface lying in facets (also the inner facets of the triangles' and
tetrahedra's diagonal splits) is counted once.
"""

import itertools
import time

import numpy as np
import pytest

import cutcells

VTK_TYPES = {"triangle": 5, "quadrilateral": 9, "tetrahedron": 10, "hexahedron": 12}

def _structured_mesh(cell, n):
    """Mesh of [-1, 1]^d with n cells per direction (VTK vertex order)."""
    tdim = 2 if cell in ("triangle", "quadrilateral") else 3
    x = np.linspace(-1.0, 1.0, n + 1)
    coords = np.column_stack(
        [g.ravel() for g in np.meshgrid(*([x] * tdim), indexing="ij")]
    )

    cells = []
    for ijk in itertools.product(range(n), repeat=tdim):
        corner = {
            off: int(np.ravel_multi_index(
                tuple(a + b for a, b in zip(ijk, off)), (n + 1,) * tdim))
            for off in itertools.product((0, 1), repeat=tdim)
        }
        if cell == "triangle":
            cells += [[corner[0, 0], corner[1, 0], corner[1, 1]],
                      [corner[0, 0], corner[1, 1], corner[0, 1]]]
        elif cell == "quadrilateral":
            cells.append([corner[0, 0], corner[1, 0], corner[1, 1], corner[0, 1]])
        elif cell == "tetrahedron":
            # Conforming Kuhn split of each cube along its main diagonal.
            for perm in itertools.permutations(range(3)):
                p = [0, 0, 0]
                tet = [corner[0, 0, 0]]
                for axis in perm:
                    p[axis] = 1
                    tet.append(corner[tuple(p)])
                cells.append(tet)
        else:
            cells.append([corner[0, 0, 0], corner[1, 0, 0], corner[1, 1, 0],
                          corner[0, 1, 0], corner[0, 0, 1], corner[1, 0, 1],
                          corner[1, 1, 1], corner[0, 1, 1]])

    nv = len(cells[0])
    return cutcells.MeshView(
        coords.astype(np.float64),
        np.asarray(cells, dtype=np.int32).ravel(),
        np.arange(0, nv * len(cells) + 1, nv, dtype=np.int32),
        np.full(len(cells), VTK_TYPES[cell], dtype=np.int32),
        tdim=tdim,
    )


# (level set, exact 2D (|phi<0|, |phi>0|, |phi=0|)) on [-1, 1]^2 with n = 4:
#  - x = 0.5 is a layer of mesh vertices and facets,
#  - x + y = 0.5 runs through mesh vertices (and hexahedron edges) and cuts
#    the cells in between diagonally,
#  - y = x runs along the diagonals that split squares into triangles and
#    cubes into tetrahedra.
# In 3D all measures are multiplied by the extent 2 in z.
LEVEL_SETS = {
    "facet_aligned": (lambda X: X[0] - 0.5, (3.0, 1.0, 2.0)),
    "vertex_diagonal": (lambda X: X[0] + X[1] - 0.5,
                        (2.875, 1.125, 1.5 * np.sqrt(2.0))),
    "split_diagonal": (lambda X: X[1] - X[0], (2.0, 2.0, 2.0 * np.sqrt(2.0))),
}


@pytest.mark.parametrize("level_set", sorted(LEVEL_SETS))
@pytest.mark.parametrize("cell", sorted(VTK_TYPES))
def test_vertex_aligned_linear_zero_sets_are_exact(cell, level_set):
    mesh = _structured_mesh(cell, 4)
    phi, exact = LEVEL_SETS[level_set]
    ls = cutcells.create_level_set(mesh, phi, degree=1, name="phi")

    start = time.perf_counter()
    result = cutcells.cut(mesh, ls, triangulate=True)
    elapsed = time.perf_counter() - start
    assert elapsed < 5.0, f"cut took {elapsed:.3f}s"

    scale = 1.0 if mesh.tdim == 2 else 2.0
    for part, value in zip(("phi < 0", "phi > 0", "phi = 0"), exact):
        weights = result[part].quadrature(order=2, mode="full").weights
        assert float(np.sum(weights)) == pytest.approx(scale * value, rel=1e-12), part
