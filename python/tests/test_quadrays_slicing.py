# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Rules summed over plane slices of a cell (MeshPart.sliced_quadrature)."""

import itertools
import math

import numpy as np
import pytest

import cutcells as cc

CENTRE = np.array([-0.21, 0.07, 0.04])
AXIS = np.array([1.0, 2.0, 2.0]) / 3.0
RADIUS = 0.6
APART = 0.5  # distance between the centres of the two balls of the lens
EXPONENTS = [e for e in itertools.product(range(5), repeat=3) if sum(e) <= 4]


def ball(centre, radius):
    return lambda X: sum((X[i] - centre[i]) ** 2 for i in range(3)) - radius * radius


def hexahedron():
    """One hexahedron [-1, 1]^3; its reference coordinates xi give x = 2 xi - 1."""
    return cc.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 2, 2, 2, cell_type="hexahedron")


def lens(mesh):
    a = cc.create_level_set(mesh, ball(CENTRE - 0.5 * APART * AXIS, RADIUS), degree=2, name="a")
    b = cc.create_level_set(mesh, ball(CENTRE + 0.5 * APART * AXIS, RADIUS), degree=2, name="b")
    return cc.cut(mesh, [a, b])["a < 0 and b < 0"]


def lens_volume():
    d, r = APART, RADIUS
    return math.pi * (2 * r - d) ** 2 * (d * d + 4 * d * r) / (12 * d)


def moments(points, weights):
    """Integrals of the monomials of total degree <= 4 in the reference coordinates."""
    return np.array([weights @ (points[:, 0] ** i * points[:, 1] ** j * points[:, 2] ** k) for i, j, k in EXPONENTS])


def rule(part, sweep, origin, direction, positions, weights, q=6):
    r = part.sliced_quadrature(0, sweep, np.asarray(origin, float), np.asarray(direction, float),
                               np.asarray(positions, float), np.asarray(weights, float), points_per_segment=q)
    w = np.asarray(r.weights)
    return np.asarray(r.points).reshape(len(w), 3), w


def gauss(a, b, n):
    t, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (a + b) + 0.5 * (b - a) * t, 0.5 * (b - a) * w


def reference(part):
    r = part.quadrature(order=10, mode="cut_only")
    w = np.asarray(r.weights)
    return np.asarray(r.points).reshape(len(w), 3), w


def test_one_ball_parallel_planes():
    """The slices of a ball are discs whose integrals are polynomials of the position."""
    mesh = hexahedron()
    part = cc.cut(mesh, [cc.create_level_set(mesh, ball(CENTRE, RADIUS), degree=2, name="a")])["a < 0"]
    t, w = gauss(-RADIUS, RADIUS, 6)
    X, W = rule(part, "parallel", CENTRE, AXIS, t, w)
    assert (W > 0).all()
    assert abs(W.sum() - 4 / 3 * math.pi * RADIUS**3) < 1e-12
    assert np.abs(moments(X, W) - moments(*reference(part))).max() < 1e-11


def test_lens_parallel_planes_split_at_the_ridge():
    """The ridge of the lens lies in the slice at 0: two panels, each a polynomial."""
    part = lens(hexahedron())
    half = RADIUS - 0.5 * APART
    t0, w0 = gauss(-half, 0.0, 6)
    t1, w1 = gauss(0.0, half, 6)
    X, W = rule(part, "parallel", CENTRE, AXIS, np.concatenate([t0, t1]), np.concatenate([w0, w1]))
    assert (W > 0).all()
    assert abs(W.sum() - lens_volume()) < 1e-12
    assert np.abs(moments(X, W) - moments(*reference(part))).max() < 1e-10
    # one panel across the ridge is not a rule for this part
    t, w = gauss(-half, half, 12)
    assert abs(rule(part, "parallel", CENTRE, AXIS, t, w)[1].sum() - lens_volume()) > 1e-7


def test_lens_turning_planes():
    """Half-planes about a line through the ridge's centre cross the ridge everywhere: no event."""
    part = lens(hexahedron())
    tilt = np.cross(AXIS, [0.0, 0.0, 1.0])
    tilt /= np.linalg.norm(tilt)
    direction = math.cos(math.radians(12)) * AXIS + math.sin(math.radians(12)) * tilt
    n = 24
    X, W = rule(part, "turning", CENTRE + 0.04 * tilt, direction, 2 * math.pi * (np.arange(n) + 0.37) / n,
                np.full(n, 2 * math.pi / n))
    assert (W > 0).all()
    assert abs(W.sum() - lens_volume()) < 1e-11
    assert np.abs(moments(X, W) - moments(*reference(part))).max() < 1e-10
    # far fewer points than the height directions take on this cell
    assert len(W) < len(reference(part)[1]) / 4


def test_ball_in_a_tetrahedron():
    """Slices of a tetrahedron are clipped by its faces."""
    v = np.array([[-1.0, -0.9, -0.8], [1.1, -0.7, -0.6], [-0.2, 1.2, -0.7], [0.1, 0.0, 0.9]])
    mesh = cc.MeshView(v, np.arange(4, dtype=np.int32), np.array([0, 4], np.int32), np.array([10], np.int32), tdim=3)
    centre, radius = v.mean(0), 0.25
    part = cc.cut(mesh, [cc.create_level_set(mesh, ball(centre, radius), degree=2, name="a")])["a < 0"]
    t, w = gauss(-radius, radius, 6)
    X, W = rule(part, "parallel", centre, [0.3, -0.5, 0.8], t, w)
    assert (W > 0).all()
    assert abs(W.sum() - 4 / 3 * math.pi * radius**3) < 1e-12
    x = v[0] + X @ (v[1:] - v[0])
    assert np.abs(W @ x / W.sum() - centre).max() < 1e-12


def test_a_ball_that_leaves_the_cell():
    """Where the cell's faces cut the part the slices are clipped: the rule is that of the part in the cell."""
    mesh = hexahedron()
    centre = np.array([0.7, 0.1, -0.2])
    part = cc.cut(mesh, [cc.create_level_set(mesh, ball(centre, RADIUS), degree=2, name="a")])["a < 0"]
    # planes orthogonal to the face x = 1: the part ends there, and the slice integral is a polynomial up to it
    t, w = gauss(-RADIUS, 1.0 - centre[0], 6)
    X, W = rule(part, "parallel", centre, [1.0, 0.0, 0.0], t, w)
    h = RADIUS - (1.0 - centre[0])
    exact = 4 / 3 * math.pi * RADIUS**3 - math.pi * h * h * (3 * RADIUS - h) / 3
    assert (W > 0).all()
    assert abs(W.sum() - exact) < 1e-12
    assert (2 * X[:, 0] - 1 < 1.0).all()


def test_sliced_quadrature_refuses_what_it_does_not_do():
    mesh = hexahedron()
    result = cc.cut(mesh, [cc.create_level_set(mesh, ball(CENTRE, RADIUS), degree=2, name="a")])
    t, w = gauss(-RADIUS, RADIUS, 4)
    with pytest.raises(Exception):
        rule(result["a = 0"], "parallel", CENTRE, AXIS, t, w)
    with pytest.raises(Exception):
        rule(result["a < 0"], "spiral", CENTRE, AXIS, t, w)
    with pytest.raises(Exception):
        rule(result["a < 0"], "parallel", CENTRE, AXIS, t, w[:-1])
