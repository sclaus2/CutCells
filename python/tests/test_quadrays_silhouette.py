# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Splitting a cut cell at the silhouette of its level set instead of bisecting it:
simple shapes take few points and no bisections, at the accuracy of a bisected rule."""

import math

import numpy as np
import pytest

import cutcells as cc


def box_mesh(n):
    """A grid of n hexahedra per direction on [-1, 1]^3."""
    return cc.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n + 1, n + 1, n + 1, cell_type="hexahedron")


def sphere(centre, radius):
    def phi(x):
        return sum((x[i] - centre[i]) ** 2 for i in range(3)) - radius**2

    return phi


def legendre_squares(points):
    """Squares of Legendre polynomials of degree 4 along each axis of the cell, as a
    mass matrix sees shape functions of degree 4 (points in the cell's coordinates)."""
    p4 = np.polynomial.legendre.Legendre.basis(4)
    return sum(p4(2 * points[:, j] - 1) ** 2 for j in range(3))


def clipped_ball(radius, half):
    """Volume and area of a ball about the centre of a cube of half width half, when
    the faces cut caps that do not meet."""
    h = radius - half
    return (4 / 3 * math.pi * radius**3 - 6 * math.pi * h**2 * (3 * radius - h) / 3,
            4 * math.pi * radius**2 - 6 * 2 * math.pi * radius * h)


def statistics(part, **kw):
    rules = part.quadrature(mode="cut_only", **kw)
    return rules, np.diff(np.asarray(rules.offset))


@pytest.fixture
def ball_in_cell():
    """A ball slightly larger than the cell [0, 0.5]^3, its caps in the neighbours."""
    mesh = box_mesh(4)
    phi = cc.create_level_set(mesh, sphere((0.25, 0.25, 0.25), 0.26), degree=2, name="phi")
    return cc.cut(mesh, phi)


def cell_of(rules, sizes):
    """The points and weights of the cell with the most points."""
    i = int(np.argmax(sizes))
    off = np.asarray(rules.offset)
    return (np.asarray(rules.points).reshape(-1, 3)[off[i]:off[i + 1]], np.asarray(rules.weights)[off[i]:off[i + 1]])


def test_ball_in_a_cell_volume(ball_in_cell):
    rules, sizes = statistics(ball_in_cell["phi < 0"], order=8)
    points, weights = cell_of(rules, sizes)
    volume, _ = clipped_ball(0.26, 0.25)
    assert sizes.max() <= 10000
    assert weights.min() > 0
    assert abs(weights.sum() - volume) < 1e-12


def test_ball_in_a_cell_interface(ball_in_cell):
    rules, sizes = statistics(ball_in_cell["phi = 0"], order=8)
    _, area = clipped_ball(0.26, 0.25)
    points, weights = cell_of(rules, sizes)
    assert sizes.max() <= 4000
    assert weights.min() > 0
    assert abs(weights.sum() - area) < 1e-11


def test_ball_in_a_cell_integrates_shape_functions(ball_in_cell):
    """What assembly sees: against the rule without silhouette splits."""
    options = cc.QuadraysOptions()
    options.silhouette_splits = False
    options.max_points = -1
    for selection in ("phi < 0", "phi = 0"):
        rules, sizes = statistics(ball_in_cell[selection], order=8)
        reference, _ = statistics(ball_in_cell[selection], order=8, options=options)
        assert np.array_equal(np.asarray(rules.parent_map), np.asarray(reference.parent_map))
        for mine, theirs in zip(np.split(np.arange(len(rules.weights)), np.asarray(rules.offset)[1:-1]),
                                np.split(np.arange(len(reference.weights)), np.asarray(reference.offset)[1:-1])):
            x = np.asarray(rules.points).reshape(-1, 3)[mine]
            y = np.asarray(reference.points).reshape(-1, 3)[theirs]
            a = np.asarray(rules.weights)[mine] @ legendre_squares(x)
            b = np.asarray(reference.weights)[theirs] @ legendre_squares(y)
            assert abs(a - b) < 1e-10


@pytest.mark.parametrize("n", [5, 9])
def test_ball_on_a_mesh(n):
    """Few points on every cell, and the exact measures, at any mesh size."""
    mesh = box_mesh(n)
    phi = cc.create_level_set(mesh, sphere((-0.21, 0.07, 0.04), 0.5), degree=2, name="phi")
    result = cc.cut(mesh, phi)
    volume = result["phi < 0"].quadrature(order=8, mode="full")
    area = result["phi = 0"].quadrature(order=8, mode="cut_only")
    assert abs(np.sum(volume.weights) - 4 / 3 * math.pi * 0.5**3) < 1e-12
    assert abs(np.sum(area.weights) - 4 * math.pi * 0.5**2) < 1e-11
    assert np.asarray(volume.weights).min() > 0 and np.asarray(area.weights).min() > 0
    _, sizes = statistics(result["phi < 0"], order=8)
    assert sizes.max() <= {5: 4000, 9: 6000}[n]


def test_tilted_ellipsoid():
    """Its silhouette is no plane: it crosses the faces of a cell, where it touches
    the zero set's traces on them."""
    a = np.array([[1.3, 0.4, -0.2], [0.4, 1.0, 0.3], [-0.2, 0.3, 2.1]])
    centre = np.array([0.03, -0.02, 0.05])

    def phi(x):
        d = [x[i] - centre[i] for i in range(3)]
        return sum(a[i, j] * d[i] * d[j] for i in range(3) for j in range(3)) - 0.25

    exact = 4 / 3 * math.pi * 0.25**1.5 / math.sqrt(np.linalg.det(a))
    for n in (2, 3, 5):
        mesh = box_mesh(n)
        result = cc.cut(mesh, cc.create_level_set(mesh, phi, degree=2, name="phi"))
        rules = result["phi < 0"].quadrature(order=8, mode="full")
        assert abs(np.sum(rules.weights) - exact) < 1e-11
        assert np.asarray(rules.weights).min() > 0
        _, sizes = statistics(result["phi < 0"], order=8)
        assert sizes.max() <= 15000


def test_cells_where_the_equator_lies_in_a_face():
    """The silhouette in a cell face coincides with the zero set's trace on it."""
    mesh = box_mesh(4)
    phi = cc.create_level_set(mesh, sphere((0.1, 0.07, 0.0), 0.4), degree=2, name="phi")
    result = cc.cut(mesh, phi)
    assert abs(np.sum(result["phi < 0"].quadrature(order=8, mode="full").weights) - 4 / 3 * math.pi * 0.4**3) < 1e-12
    _, sizes = statistics(result["phi < 0"], order=8)
    assert sizes.max() <= 15000


def test_splits_can_be_turned_off(ball_in_cell):
    options = cc.QuadraysOptions()
    assert options.silhouette_splits
    options.silhouette_splits = False
    options.max_points = -1
    _, sizes = statistics(ball_in_cell["phi < 0"], order=8, options=options)
    assert sizes.max() > 10000
