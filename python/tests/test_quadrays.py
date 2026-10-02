# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""The quadrays backend through HOMeshPart.quadrature and the per-cell bindings."""

import itertools
import math

import numpy as np
import pytest

import cutcells

CENTRE = np.array([0.0123, -0.0371, 0.0217])
RADIUS = 0.7


def box_mesh(kind: str, n: int) -> cutcells.MeshView:
    """[-1, 1]^3 with n^3 hexahedra or 6 n^3 Kuhn tetrahedra, VTK vertex order."""
    g = np.linspace(-1.0, 1.0, n + 1)
    coords = np.array([[x, y, z] for z in g for y in g for x in g], dtype=np.float64)

    def node(i, j, k):
        return i + (n + 1) * (j + (n + 1) * k)

    connectivity, types = [], []
    for k, j, i in itertools.product(range(n), repeat=3):
        if kind == "hex":
            connectivity += [node(i, j, k), node(i + 1, j, k), node(i + 1, j + 1, k), node(i, j + 1, k),
                             node(i, j, k + 1), node(i + 1, j, k + 1), node(i + 1, j + 1, k + 1),
                             node(i, j + 1, k + 1)]
            types.append(12)
            continue
        for perm in itertools.permutations(range(3)):
            vertex = [i, j, k]
            path = [tuple(vertex)]
            for axis in perm:
                vertex[axis] += 1
                path.append(tuple(vertex))
            connectivity += [node(*v) for v in path]
            types.append(10)
    width = 8 if kind == "hex" else 4
    offsets = np.arange(0, len(connectivity) + 1, width, dtype=np.int32)
    return cutcells.MeshView(coords, np.array(connectivity, dtype=np.int32), offsets,
                             np.array(types, dtype=np.int32), tdim=3)


def plane(x):
    return x[0] + 0.3 * x[1] - 0.2 * x[2]


def sphere(x):
    return (x[0] - CENTRE[0]) ** 2 + (x[1] - CENTRE[1]) ** 2 + (x[2] - CENTRE[2]) ** 2 - RADIUS**2


@pytest.mark.parametrize("kind,n", [("hex", 5), ("tet", 4)])
@pytest.mark.parametrize("degree", [1, 2])
def test_plane_is_exact(kind, n, degree):
    mesh = box_mesh(kind, n)
    result = cutcells.cut(mesh, cutcells.create_level_set(mesh, plane, degree=degree, name="phi"))
    for selection, exact in [("phi < 0", 4.0), ("phi > 0", 4.0), ("phi = 0", 4.0 * math.sqrt(1.13))]:
        rules = result[selection].quadrature(order=3, mode="full", backend="quadrays")
        assert np.sum(rules.weights) == pytest.approx(exact, rel=1e-13), selection
        points = np.asarray(rules.points).reshape(-1, 3)
        assert np.all(points >= -1e-12) and np.all(points <= 1 + 1e-12)


@pytest.mark.parametrize("kind", ["hex", "tet"])
def test_sphere_totals(kind):
    mesh = box_mesh(kind, 8)
    result = cutcells.cut(mesh, cutcells.create_level_set(mesh, sphere, degree=2, name="phi"))
    volume = result["phi < 0"].quadrature(order=5, mode="full", backend="quadrays")
    area = result["phi = 0"].quadrature(order=5, mode="cut_only", backend="quadrays")
    assert np.sum(volume.weights) == pytest.approx(4.0 / 3.0 * math.pi * RADIUS**3, rel=1e-7)
    assert np.sum(area.weights) == pytest.approx(4.0 * math.pi * RADIUS**2, rel=1e-6)
    assert np.all(np.asarray(volume.weights) > 0) and np.all(np.asarray(area.weights) > 0)
    # at most one rule per cut cell, none empty: cells the front end flags as cut
    # may hold no interface
    assert len(area.parent_map) <= len(result["phi = 0"].cut_cell_ids)
    assert np.all(np.diff(area.offset) > 0)


def test_options_reach_the_engine():
    mesh = box_mesh("tet", 4)
    result = cutcells.cut(mesh, cutcells.create_level_set(mesh, sphere, degree=2, name="phi"))
    part = result["phi = 0"]
    default = cutcells.quadrays_quadrature(part, order=3, mode="cut_only")
    options = cutcells.QuadraysOptions()
    options.margin = 0.5  # stricter certification bisects more
    strict = cutcells.quadrays_quadrature(part, order=3, mode="cut_only", options=options)
    assert len(strict.weights) > len(default.weights)
    assert np.sum(strict.weights) == pytest.approx(np.sum(default.weights), rel=1e-3)


def test_cell_rules_and_stats():
    tet = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    # degree-1 simplex Bernstein coefficients are vertex values: phi = x + y + z - 0.5
    coeffs = np.array([-0.5, 0.5, 0.5, 0.5])
    rules, stats = cutcells.quadrays_cell_rules(cutcells.CellType.tetrahedron, tet, 1, coeffs, "phi < 0", q=3)
    assert np.sum(rules.weights) == pytest.approx(0.5**3 / 6, rel=1e-14)
    assert stats.bisections == 0
    surface, _ = cutcells.quadrays_cell_rules(cutcells.CellType.tetrahedron, tet, 1, coeffs, "psi = 0",
                                              level_set_name="psi")
    assert np.sum(surface.weights) == pytest.approx(math.sqrt(3) / 2 * 0.5**2, rel=1e-14)

    rules32, _ = cutcells.quadrays_cell_rules_float32(cutcells.CellType.tetrahedron, tet.astype(np.float32), 1,
                                                       coeffs.astype(np.float32), "phi < 0")
    assert np.sum(rules32.weights) == pytest.approx(0.5**3 / 6, rel=1e-6)


def test_terms_on_one_level_set():
    mesh = box_mesh("hex", 2)
    result = cutcells.cut(mesh, cutcells.create_level_set(mesh, plane, degree=1, name="phi"))
    rules = result["phi < 0 or phi > 0"].quadrature(order=3, backend="quadrays")
    assert np.sum(rules.weights) == pytest.approx(8.0, rel=1e-13)
