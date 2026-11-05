# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Curve parts in cells of dimension 3 follow the curve (QuadraysOptions.follow_curves).

The curves here are known exactly: circles and two ellipses that cross. Their lengths must come out on every
cell type and at any number of points, also where the ellipses cross inside a cell, with few points per cell.
"""

import itertools
import math

import numpy as np
import pytest

import cutcells as cc

C = np.array((0.0123, -0.0371, 0.0217))
F = np.array((-0.21, 0.07, 0.04))
G = F + 0.5 * np.array((1.0, 2.0, 2.0)) / 3.0
NORMAL = np.array((0.93, 0.33, -0.19))
R = 0.6


def ball(c):
    return lambda x: sum((x[i] - c[i]) ** 2 for i in range(3)) - R**2


def cylinder(axis):
    return lambda x: sum((x[i] - C[i]) ** 2 for i in range(3) if i != axis) - R**2


def plane(n, d):
    return lambda x: sum(n[i] * x[i] for i in range(3)) - d


def wedge_mesh(kind, n):
    """[-1, 1]^3 with 2 n^3 prisms or 6 n^3 pyramids, VTK vertex order."""
    g = np.linspace(-1.0, 1.0, n + 1)
    coords = [[x, y, z] for z in g for y in g for x in g]

    def node(i, j, k):
        return i + (n + 1) * (j + (n + 1) * k)

    connectivity, types = [], []
    for k, j, i in itertools.product(range(n), repeat=3):
        if kind == "prism":
            for tri in ([(i, j), (i + 1, j), (i, j + 1)], [(i + 1, j), (i + 1, j + 1), (i, j + 1)]):
                connectivity += [node(a, b, k) for a, b in tri] + [node(a, b, k + 1) for a, b in tri]
                types.append(13)
            continue
        apex = len(coords)
        coords.append([(g[i] + g[i + 1]) / 2, (g[j] + g[j + 1]) / 2, (g[k] + g[k + 1]) / 2])
        for axis, side in itertools.product(range(3), range(2)):
            corners = []
            for a, b in [(0, 0), (1, 0), (1, 1), (0, 1)]:
                v = [i, j, k]
                v[axis] += side
                v[(axis + 1) % 3] += a
                v[(axis + 2) % 3] += b
                corners.append(node(*v))
            connectivity += corners + [apex]
            types.append(14)
    width = 6 if kind == "prism" else 5
    offsets = np.arange(0, len(connectivity) + 1, width, dtype=np.int32)
    return cc.MeshView(np.array(coords, dtype=np.float64), np.array(connectivity, dtype=np.int32), offsets,
                       np.array(types, dtype=np.int32), tdim=3)


def mesh_of(kind, n):
    if kind in ("prism", "pyramid"):
        return wedge_mesh(kind, n)
    return cc.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n + 1, n + 1, n + 1, cell_type=kind)


def level_sets(case, mesh):
    if case == "lens":
        return [cc.create_level_set(mesh, ball(F), degree=2, name="a"), cc.create_level_set(mesh, ball(G), degree=2, name="b")]
    if case == "ball_plane":
        return [cc.create_level_set(mesh, ball(F), degree=2, name="a"),
                cc.create_level_set(mesh, plane(NORMAL, 0.05), degree=1, name="b")]
    return [cc.create_level_set(mesh, cylinder(2), degree=2, name="a"), cc.create_level_set(mesh, cylinder(1), degree=2, name="b")]


def exact_length(case):
    if case == "lens":
        return 2 * math.pi * math.sqrt(R**2 - (np.linalg.norm(G - F) / 2) ** 2)
    if case == "ball_plane":
        distance = (0.05 - NORMAL @ F) / np.linalg.norm(NORMAL)
        return 2 * math.pi * math.sqrt(R**2 - distance**2)
    # two ellipses of semi-axes r and r sqrt(2)
    x, w = np.polynomial.legendre.leggauss(200)
    t = math.pi * (x + 1)
    return 2 * math.pi * float(w @ (R * np.sqrt(1 + np.cos(t) ** 2)))


def rules(result, q, **attributes):
    options = cc.QuadraysOptions()
    options.points_per_segment = q
    options.max_points = -1
    for key, value in attributes.items():
        setattr(options, key, value)
    return result["a = 0 and b = 0"].quadrature(order=8, mode="cut_only", options=options)


@pytest.mark.parametrize("case", ["lens", "ball_plane", "two_cylinders"])
@pytest.mark.parametrize("kind", ["hexahedron", "tetrahedron", "prism", "pyramid"])
@pytest.mark.parametrize("n", [5, 9])
def test_length_on_every_cell_type(case, kind, n):
    """Exact at any q: rays dropped a point on a crossing at odd q, and lost pieces on prisms and pyramids."""
    mesh = mesh_of(kind, n)
    result = cc.cut(mesh, level_sets(case, mesh))
    exact = exact_length(case)
    for q in (2, 3, 5, 6):
        r = rules(result, q)
        weights = np.asarray(r.weights)
        assert weights.min() > 0
        assert weights.sum() == pytest.approx(exact, abs=1e-10)


@pytest.mark.parametrize("kind", ["hexahedron", "tetrahedron"])
def test_points_lie_on_both_zero_sets(kind):
    n = 5
    mesh = mesh_of(kind, n)
    result = cc.cut(mesh, level_sets("two_cylinders", mesh))
    r = rules(result, 6)
    width = 8 if kind == "hexahedron" else 4
    coords = np.asarray(mesh.coordinates).reshape(-1, 3)
    conn = np.asarray(mesh.connectivity).reshape(-1, width)
    points = np.asarray(r.points).reshape(-1, 3)
    cells = np.repeat(np.asarray(r.parent_map), np.diff(np.asarray(r.offset)))
    axes = (1, 3, 4) if kind == "hexahedron" else (1, 2, 3)
    x0 = coords[conn[cells, 0]]
    x = x0 + sum(points[:, j, None] * (coords[conn[cells, axes[j]]] - x0) for j in range(3))
    for fn in (cylinder(2), cylinder(1)):
        assert np.max(np.abs(fn(x.T))) < 1e-12


def test_few_points_where_the_ellipses_cross():
    """Four arcs at most meet in the cell of a crossing, each with a bounded rule: no bisections to it."""
    mesh = mesh_of("hexahedron", 5)
    result = cc.cut(mesh, level_sets("two_cylinders", mesh))
    q = 6
    sizes = np.diff(np.asarray(rules(result, q).offset))
    assert sizes.max() <= 4 * (2 * q + 4)
    assert sizes.mean() < 3 * q


@pytest.mark.parametrize("case", ["lens", "ball_plane"])
def test_accuracy_at_high_order(case):
    """A moment of the mass matrix per cell, order 10: the curves that converged late under rays."""
    n, degree = 9, 5
    mesh = mesh_of("hexahedron", n)
    result = cc.cut(mesh, level_sets(case, mesh))
    def moments(rule):
        p = np.asarray(rule.points).reshape(-1, 3)
        w = np.asarray(rule.weights)
        off = np.asarray(rule.offset)
        L = np.polynomial.legendre.Legendre.basis(degree)
        shape = sum(L(2 * p[:, j] - 1) ** 2 for j in range(3))
        return {int(c): w[off[i]:off[i + 1]] @ shape[off[i]:off[i + 1]] for i, c in enumerate(np.asarray(rule.parent_map))}
    reference = moments(rules(result, 24, moment_tolerance=0.0))
    # q as order 10 gives it, the check adds the rest
    computed = moments(rules(result, 6))
    assert set(computed) == set(reference)
    assert sum(abs(computed[c] - reference[c]) for c in reference) < 1e-6


def test_closed_curve_inside_one_cell():
    mesh = mesh_of("hexahedron", 1)
    result = cc.cut(mesh, level_sets("lens", mesh))
    r = rules(result, 4)
    assert np.asarray(r.weights).sum() == pytest.approx(exact_length("lens"), abs=1e-10)
    assert len(np.asarray(r.weights)) < 100


@pytest.mark.parametrize("kind", ["hexahedron", "tetrahedron"])
def test_further_level_set_splits_the_arcs(kind):
    """The curve below a plane, against the exact curve."""
    mesh = mesh_of(kind, 5)
    normal, offset = np.array((0.3, -0.5, 0.8)), 0.1
    result = cc.cut(mesh, level_sets("two_cylinders", mesh) + [cc.create_level_set(mesh, plane(normal, offset), degree=1, name="c")])
    options = cc.QuadraysOptions()
    options.points_per_segment = 6
    r = result["a = 0 and b = 0 and c < 0"].quadrature(order=8, mode="cut_only", options=options)
    # two ellipses y - C1 = +- (z - C2): the plane cuts each where A cos t + B sin t + D = 0
    x, w = np.polynomial.legendre.leggauss(200)
    total = 0.0
    for s in (1.0, -1.0):
        A, B, D = R * normal[0], R * (normal[1] + s * normal[2]), normal @ C - offset
        phase, half = math.atan2(B, A), math.acos(-D / math.hypot(A, B))
        t0, t1 = phase + half, phase - half + 2 * math.pi  # A cos t + B sin t + D < 0 between them
        t = 0.5 * (t0 + t1) + 0.5 * (t1 - t0) * x
        assert A * math.cos(0.5 * (t0 + t1)) + B * math.sin(0.5 * (t0 + t1)) + D < 0
        total += 0.5 * (t1 - t0) * float(w @ (R * np.sqrt(1 + np.cos(t) ** 2)))
    assert np.asarray(r.weights).sum() == pytest.approx(total, abs=1e-10)


@pytest.mark.parametrize("distance", [False, True])
def test_analytic_level_sets_take_the_same_path(distance):
    mesh = mesh_of("hexahedron", 5)
    a = cc.analytic_sphere(list(F), R, signed_distance=distance, name="a")
    b = cc.analytic_sphere(list(G), R, signed_distance=distance, name="b")
    result = cc.cut(mesh, [a, b])
    assert np.asarray(rules(result, 5).weights).sum() == pytest.approx(exact_length("lens"), abs=1e-10)


@pytest.mark.parametrize("gap", [1e-3, 1e-6, 1e-9, -1e-9, -1e-6])
def test_curve_next_to_a_mesh_face(gap):
    """A plane just off the mesh faces at x = 0.2: the circle stays whole, in the cells on its side."""
    mesh = mesh_of("hexahedron", 5)
    a = cc.create_level_set(mesh, ball(F), degree=2, name="a")
    b = cc.create_level_set(mesh, plane((1.0, 0.0, 0.0), 0.2 + gap), degree=1, name="b")
    result = cc.cut(mesh, [a, b])
    exact = 2 * math.pi * math.sqrt(R**2 - (0.2 + gap - F[0]) ** 2)
    assert np.asarray(rules(result, 6).weights).sum() == pytest.approx(exact, abs=1e-9)


@pytest.mark.parametrize("kind", ["hexahedron", "tetrahedron", "prism", "pyramid"])
def test_curve_in_mesh_faces_counts_once(kind):
    """A plane through mesh faces at x = 0.2: one cell per face takes the curve, moved just off the face."""
    mesh = mesh_of(kind, 5)
    a = cc.create_level_set(mesh, ball(F), degree=2, name="a")
    b = cc.create_level_set(mesh, plane((1.0, 0.0, 0.0), 0.2), degree=1, name="b")
    result = cc.cut(mesh, [a, b])
    exact = 2 * math.pi * math.sqrt(R**2 - (0.2 - F[0]) ** 2)
    for mode in ("cut_only", "full"):
        weights = np.asarray(result["a = 0 and b = 0"].quadrature(order=8, mode=mode).weights)
        assert weights.min() > 0
        assert weights.sum() == pytest.approx(exact, abs=2e-8)


def test_analytic_curve_in_mesh_faces_counts_once():
    mesh = mesh_of("hexahedron", 5)
    a = cc.analytic_sphere(list(F), R, signed_distance=False, name="a")
    b = cc.create_level_set(mesh, plane((1.0, 0.0, 0.0), 0.2), degree=1, name="b")
    weights = np.asarray(cc.cut(mesh, [a, b])["a = 0 and b = 0"].quadrature(order=8).weights)
    assert weights.sum() == pytest.approx(2 * math.pi * math.sqrt(R**2 - (0.2 - F[0]) ** 2), abs=2e-8)


def test_lookup_tables_do_not_find_a_curve_in_mesh_faces():
    mesh = mesh_of("hexahedron", 5)
    a = cc.create_level_set(mesh, ball(F), degree=2, name="a")
    b = cc.create_level_set(mesh, plane((1.0, 0.0, 0.0), 0.2), degree=1, name="b")
    rule = cc.cut(mesh, [a, b])["a = 0 and b = 0"].quadrature(order=2, backend="lut")
    assert len(np.asarray(rule.weights)) == 0


def test_height_lines_stay_available():
    mesh = mesh_of("hexahedron", 5)
    result = cc.cut(mesh, level_sets("lens", mesh))
    r = rules(result, 6, follow_curves=False)
    assert np.asarray(r.weights).sum() == pytest.approx(exact_length("lens"), abs=1e-6)
