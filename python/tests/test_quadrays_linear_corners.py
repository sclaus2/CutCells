# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Level sets that are planes or lines, meeting inside a cell, integrate exactly.

They become linear functions in quadrays, as the cell's clip planes are: their
corners split the base exactly instead of being searched for. The cases come
from the CutElementIntegration benchmark (Loibl et al. 2026), where these
corners were off by 1e-4 to 1e-2.
"""

import math

import numpy as np
import pytest

import cutcells


def _line(p1, p2):
    """geo_line_through_points: negative on the left of p1 -> p2."""
    d = np.array(p2, float) - np.array(p1, float)
    n = np.array([-d[1], d[0]]) / np.linalg.norm(d)
    return lambda X: n[0] * (p1[0] - X[0]) + n[1] * (p1[1] - X[1])


def _polygon_area(points, cells, q):
    mesh = cutcells.rectangle_mesh(0.0, 0.0, 1.0, 1.0, cells + 1, cells + 1, cell_type="quadrilateral")
    names = [f"l{i}" for i in range(len(points) - 1)]
    level_sets = [
        cutcells.create_level_set(mesh, _line(points[i], points[i + 1]), degree=1, name=names[i])
        for i in range(len(points) - 1)
    ]
    result = cutcells.cut(mesh, level_sets)
    options = cutcells.QuadraysOptions()
    options.points_per_segment = q
    rules = result[" and ".join(f"{n} < 0" for n in names)].quadrature(order=2 * q - 1, options=options)
    return float(np.asarray(rules.weights).sum())


@pytest.mark.parametrize("cells,q", [(8, 2), (16, 3)])
def test_sliver_triangle(cells, q):
    """Three lines, two of them almost parallel, meeting at a mesh vertex."""
    d = 1.42799666072263
    p1, p2, p3 = [0.0, 0.0], [0.5, 0.0], [0.5 * math.cos(math.pi / 2 - d), 0.5 * math.sin(math.pi / 2 - d)]
    exact = 0.5 * abs(p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1]))
    assert _polygon_area([p1, p2, p3, p1], cells, q) == pytest.approx(exact, rel=1e-13)


@pytest.mark.parametrize("cells", [1, 2, 4])
def test_cusp(cells):
    """A kite whose sharp corner lies inside a cell."""
    p1, p2, p3, p4 = [1.0, 1.0], [0.7, 1.0], [0.6, 0.6], [1.0, 0.7]
    b = math.dist(p2, p4)
    exact = 0.5 * (p1[1] - p4[1]) * (p1[0] - p2[0]) + 0.5 * b * math.sqrt(math.dist(p2, p3) ** 2 - (b / 2) ** 2)
    assert _polygon_area([p1, p2, p3, p4, p1], cells, 3) == pytest.approx(exact, rel=1e-13)


def _rotation(axis, a):
    c, s = math.cos(a), math.sin(a)
    if axis == 0:
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    if axis == 1:
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


@pytest.mark.parametrize("cell_type", ["hexahedron", "tetrahedron"])
def test_rotated_cube(cell_type):
    """Six planes; their edges and corners lie inside cells."""
    s = 0.5
    corners = np.array([[0, s, 0, s, 0, s, 0, s], [0, 0, s, s, 0, 0, s, s], [0, 0, 0, 0, s, s, s, s]], float)
    corners = _rotation(0, math.pi / 2) @ _rotation(1, math.pi / 3) @ _rotation(2, math.pi / 4) @ corners
    corners -= 0.5 * (corners[:, [0]] + corners[:, [7]])
    faces = [(1, 3, 2), (1, 2, 5), (1, 5, 3), (8, 7, 6), (8, 6, 4), (8, 4, 7)]
    mesh = cutcells.box_mesh(-0.51, -0.51, -0.51, 0.51, 0.51, 0.51, 3, 3, 3, cell_type=cell_type)
    level_sets = []
    for i, (a, b, c) in enumerate(faces):
        p = corners[:, a - 1]
        n = np.cross(corners[:, b - 1] - p, corners[:, c - 1] - p)
        n /= np.linalg.norm(n)
        level_sets.append(
            cutcells.create_level_set(
                mesh, lambda X, n=n, p=p: n[0] * (X[0] - p[0]) + n[1] * (X[1] - p[1]) + n[2] * (X[2] - p[2]),
                degree=1, name=f"p{i}"))
    result = cutcells.cut(mesh, level_sets)
    options = cutcells.QuadraysOptions()
    options.points_per_segment = 2
    rules = result[" and ".join(f"p{i} < 0" for i in range(6))].quadrature(order=3, options=options)
    assert float(np.asarray(rules.weights).sum()) == pytest.approx(s ** 3, rel=1e-13)


@pytest.mark.parametrize("cell_type", ["triangle", "tetrahedron", "prism"])
def test_reference_frame_agrees(cell_type):
    """The reference frame and the orthogonal one integrate the same part."""
    if cell_type == "triangle":
        mesh = cutcells.rectangle_mesh(-1.0, -1.0, 1.0, 1.0, 6, 6, cell_type="triangle")
        phi = lambda X: (X[0] - 0.1) ** 2 + (X[1] + 0.05) ** 2 - 0.36
        exact = math.pi * 0.36
    else:
        if cell_type == "prism":
            g = np.linspace(-1.0, 1.0, 5)
            coords = [[x, y, z] for z in g for y in g for x in g]
            node = lambda i, j, k: i + 5 * (j + 5 * k)
            conn = []
            for k in range(4):
                for j in range(4):
                    for i in range(4):
                        for tri in ([(i, j), (i + 1, j), (i, j + 1)], [(i + 1, j), (i + 1, j + 1), (i, j + 1)]):
                            conn += [node(a, b, k) for a, b in tri] + [node(a, b, k + 1) for a, b in tri]
            m = len(conn) // 6
            mesh = cutcells.MeshView(np.asarray(coords, float), np.asarray(conn, np.int32),
                                     np.arange(0, 6 * m + 1, 6, dtype=np.int32), np.full(m, 13, np.int32), tdim=3)
        else:
            mesh = cutcells.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 5, 5, 5, cell_type=cell_type)
        phi = lambda X: (X[0] - 0.1) ** 2 + (X[1] + 0.05) ** 2 + (X[2] - 0.02) ** 2 - 0.36
        exact = 4 / 3 * math.pi * 0.6 ** 3
    result = cutcells.cut(mesh, [cutcells.create_level_set(mesh, phi, degree=2, name="a")])
    values, moments = [], []
    for reference in (False, True):
        options = cutcells.QuadraysOptions()
        options.reference_frame = reference
        rules = result["a < 0"].quadrature(order=8, options=options)
        w = np.asarray(rules.weights)
        points = np.asarray(rules.points).reshape(len(w), -1)   # parent reference coordinates
        values.append(float(w.sum()))
        moments.append(float(w @ (points ** 2).sum(axis=1)))
    assert values[0] == pytest.approx(exact, rel=1e-10)
    assert values[1] == pytest.approx(exact, rel=1e-10)
    # the points, not only the weights: both frames report parent reference coordinates
    assert moments[1] == pytest.approx(moments[0], rel=1e-10)
