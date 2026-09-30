# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier:    MIT
"""Interface measures of selections combining several linear level sets.

A later level set cuts the leaves on both sides of an earlier interface and
may triangulate the same piece of that interface with different diagonals on
the two sides, e.g. when it splits a quadrilateral interface piece. Each piece
must be registered once. Before, the ownership test relied on certification
tags that are reset for the children of later cuts: pieces were counted from
both sides, or dropped when the untouched side was positive.

The reference is exact for P1 level sets: per simplex, the planar zero set of
one level set is clipped by the sign conditions of the others.
"""

import itertools

import numpy as np
import pytest

import cutcells

VTK_TYPES = {"triangle": 5, "tetrahedron": 10}


def _structured_mesh(cell, n):
    """Mesh of [-1, 1]^d with n cells per direction (VTK vertex order)."""
    tdim = 2 if cell == "triangle" else 3
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
        else:
            # Conforming Kuhn split of each cube along its main diagonal.
            for perm in itertools.permutations(range(3)):
                p = [0, 0, 0]
                tet = [corner[0, 0, 0]]
                for axis in perm:
                    p[axis] = 1
                    tet.append(corner[tuple(p)])
                cells.append(tet)

    cells = np.asarray(cells, dtype=np.int32)
    nv = cells.shape[1]
    mesh = cutcells.MeshView(
        coords.astype(np.float64),
        cells.ravel(),
        np.arange(0, nv * len(cells) + 1, nv, dtype=np.int32),
        np.full(len(cells), VTK_TYPES[cell], dtype=np.int32),
        tdim=tdim,
    )
    return mesh, coords, cells


def _clip(points, values):
    """Part of a convex polygon or a segment where the linear function with
    the given point values is <= 0 (Sutherland-Hodgman)."""
    n = len(points)
    edges = [(i, (i + 1) % n) for i in range(n)] if n > 2 else [(0, 1)]
    out = []
    for i, j in edges:
        if values[i] <= 0:
            out.append(points[i])
        if (values[i] <= 0) != (values[j] <= 0):
            t = values[i] / (values[i] - values[j])
            out.append(points[i] + t * (points[j] - points[i]))
    if n == 2 and values[1] <= 0:
        out.append(points[1])
    return out


def _measure(points, tdim):
    if tdim == 2:
        return float(np.linalg.norm(points[1] - points[0])) if len(points) == 2 else 0.0
    area = np.zeros(3)
    for p, q in zip(points[1:-1], points[2:]):
        area += np.cross(p - points[0], q - points[0])
    return 0.5 * float(np.linalg.norm(area))


def _exact_interface_measure(coords, cells, values, zero, signs):
    """Measure of {phi_zero = 0} and {s * phi_k > 0 for (k, s) in signs}.

    values[k] holds the vertex values of the P1 level set phi_k. The zero
    level set must not vanish at mesh vertices.
    """
    gdim = coords.shape[1]
    total = 0.0
    for cell in cells:
        # Points carry their coordinates and all level-set values, which are
        # linear on the simplex and hence interpolated exactly.
        P = np.column_stack([coords[cell], values[:, cell].T])
        f = P[:, gdim + zero]
        points = [
            P[i] + f[i] / (f[i] - f[j]) * (P[j] - P[i])
            for i, j in itertools.combinations(range(len(P)), 2)
            if (f[i] < 0) != (f[j] < 0)
        ]
        if not points:
            continue
        if len(points) == 4:
            # Planar quadrilateral: order the points around their centroid.
            x = [p[:gdim] for p in points]
            c = np.mean(x, axis=0)
            e1 = x[0] - c
            e2 = np.cross(np.cross(x[1] - x[0], x[2] - x[0]), e1)
            order = np.argsort([np.arctan2((p - c) @ e2, (p - c) @ e1) for p in x])
            points = [points[k] for k in order]
        for k, s in signs.items():
            points = _clip(points, [-s * p[gdim + k] for p in points])
            if len(points) < 2:
                break
        else:
            total += _measure([p[:gdim] for p in points], gdim)
    return total


def _sphere(x):
    return np.sqrt(sum(xi**2 for xi in x)) - 0.7


def _plane(x):
    return x[0] + 0.3 * x[1] - (0.2 * x[2] if len(x) > 2 else 0.0) - 0.0513


def _vertex_plane(x):
    # Vanishes at mesh vertices and, for tetrahedra, along mesh edges.
    return x[0] + x[1] - 0.5


def _third_plane(x):
    return 0.4 * x[0] - x[1] + (0.25 * x[2] if len(x) > 2 else 0.0) + 0.1234


# (level sets, selections as (expression, zero level set, sign conditions)).
CASES = {
    "generic": (
        {"a": _sphere, "b": _plane},
        [
            ("a = 0", "a", {}),
            ("a = 0 and b < 0", "a", {"b": -1}),
            ("a = 0 and b > 0", "a", {"b": 1}),
            ("b = 0 and a < 0", "b", {"a": -1}),
            ("b = 0 and a > 0", "b", {"a": 1}),
        ],
    ),
    "through_vertices": (
        {"a": _sphere, "b": _vertex_plane},
        [
            ("a = 0", "a", {}),
            ("a = 0 and b < 0", "a", {"b": -1}),
            ("a = 0 and b > 0", "a", {"b": 1}),
        ],
    ),
    "three_level_sets": (
        {"a": _sphere, "b": _plane, "c": _third_plane},
        [
            ("a = 0 and b < 0 and c > 0", "a", {"b": -1, "c": 1}),
            ("b = 0 and c < 0", "b", {"c": -1}),
            ("c = 0 and a < 0 and b > 0", "c", {"a": -1, "b": 1}),
        ],
    ),
}


@pytest.mark.parametrize("case", sorted(CASES))
@pytest.mark.parametrize("cell, n", [("triangle", 8), ("tetrahedron", 4)])
def test_multi_level_set_interface_measure_is_exact(cell, n, case):
    functions, selections = CASES[case]
    mesh, coords, cells = _structured_mesh(cell, n)
    names = list(functions)
    level_sets = [
        cutcells.create_level_set(mesh, functions[name], degree=1, name=name)
        for name in names
    ]
    values = np.array([[functions[name](x) for x in coords] for name in names])
    result = cutcells.cut(mesh, level_sets, triangulate=True)

    for expr, zero, signs in selections:
        weights = result[expr].quadrature(order=2, mode="full").weights
        exact = _exact_interface_measure(
            coords, cells, values, names.index(zero),
            {names.index(k): s for k, s in signs.items()})
        assert exact > 0.0, expr
        assert float(np.sum(weights)) == pytest.approx(exact, rel=1e-12), expr


def test_second_level_set_keeps_first_interface():
    """Cutting with a second level set leaves the first interface unchanged."""
    mesh, _, _ = _structured_mesh("tetrahedron", 4)
    a = cutcells.create_level_set(mesh, _sphere, degree=1, name="a")
    b = cutcells.create_level_set(mesh, _plane, degree=1, name="b")

    def measure(result, expr):
        return float(np.sum(result[expr].quadrature(order=2, mode="full").weights))

    single = measure(cutcells.cut(mesh, a, triangulate=True), "a = 0")
    multi = cutcells.cut(mesh, [a, b], triangulate=True)
    assert measure(multi, "a = 0") == pytest.approx(single, rel=1e-12)
    assert measure(multi, "a = 0 and b < 0") + measure(multi, "a = 0 and b > 0") \
        == pytest.approx(single, rel=1e-12)
