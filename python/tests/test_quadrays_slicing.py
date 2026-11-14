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
    a = cc.create_level_set(mesh, ball(CENTRE - 0.5 * APART * AXIS, RADIUS), degree=2, name="a")
    b = cc.create_level_set(mesh, ball(CENTRE + 0.5 * APART * AXIS, RADIUS), degree=2, name="b")
    result = cc.cut(mesh, [a, b])
    t, w = gauss(-RADIUS, RADIUS, 4)
    with pytest.raises(Exception):
        rule(result["a = 0 and b = 0"], "parallel", CENTRE, AXIS, t, w)  # a curve
    with pytest.raises(Exception):
        rule(result["a < 0 and b < 0"], "spiral", CENTRE, AXIS, t, w)
    with pytest.raises(Exception):
        rule(result["a < 0 and b < 0"], "parallel", CENTRE, AXIS, t, w[:-1])


def test_sweep_events_of_a_lens():
    """Along x the lens begins and ends where its ridge, a circle, is tangent to the slices."""
    part = lens(hexahedron())
    found, events = part.sweep_events(0, 0)
    ridge = math.sqrt(RADIUS**2 - (0.5 * APART) ** 2) * math.sqrt(1.0 - AXIS[0] ** 2)
    expected = (np.array([CENTRE[0] - ridge, CENTRE[0] + ridge]) + 1.0) / 2.0
    assert found
    assert np.asarray(events) == pytest.approx(expected, abs=2e-6)
    # along y and z a sphere's pole lies on the lens as well: more events
    for axis in (1, 2):
        found, events = part.sweep_events(0, axis)
        assert found and len(events) == 4


def test_slicing_replaces_a_costly_rule():
    """A lens inside one cell: the rule summed over slices is as accurate and has a fraction of the points."""
    part = lens(hexahedron())
    usual = cc.QuadraysOptions()
    usual.points_per_segment = 5
    sliced = cc.QuadraysOptions()
    sliced.points_per_segment = 5
    sliced.slicing = 4000
    r0 = part.quadrature(order=7, mode="cut_only", options=usual)
    r1 = part.quadrature(order=7, mode="cut_only", options=sliced)
    w0, w1 = np.asarray(r0.weights), np.asarray(r1.weights)
    assert len(w0) > 4000
    assert len(w1) < len(w0) / 2
    assert (w1 > 0).all()
    assert abs(w1.sum() - lens_volume()) < 5e-9 * lens_volume()
    X1 = np.asarray(r1.points).reshape(len(w1), 3)
    assert np.abs(moments(X1, w1) - moments(*reference(part))).max() < 5e-9


def test_slicing_leaves_cheap_cells_alone():
    """On a mesh that resolves the lens no cell reaches the threshold: the rules are the usual ones."""
    mesh = cc.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 5, 5, 5, cell_type="hexahedron")
    part = lens(mesh)
    sliced = cc.QuadraysOptions()
    sliced.slicing = 4000
    r0 = part.quadrature(order=4, mode="cut_only")
    r1 = part.quadrature(order=4, mode="cut_only", options=sliced)
    assert np.array_equal(np.asarray(r0.weights), np.asarray(r1.weights))
    assert np.array_equal(np.asarray(r0.points), np.asarray(r1.points))


def cap_area():
    """The piece of one sphere of the lens inside the other."""
    return 2 * math.pi * RADIUS * (RADIUS - 0.5 * APART)


def lens_parts(mesh):
    a = cc.create_level_set(mesh, ball(CENTRE - 0.5 * APART * AXIS, RADIUS), degree=2, name="a")
    b = cc.create_level_set(mesh, ball(CENTRE + 0.5 * APART * AXIS, RADIUS), degree=2, name="b")
    result = cc.cut(mesh, [a, b])
    return result["a < 0 and b < 0"], result["a = 0 and b < 0"]


def options(slicing, q=5):
    o = cc.QuadraysOptions()
    o.points_per_segment = q
    o.slicing = slicing
    return o


@pytest.mark.parametrize("cell_type", ["hexahedron", "tetrahedron", "prism", "pyramid"])
def test_slicing_on_every_cell_type(cell_type):
    """The lens in the cells of one grid cell: with every cut cell tried, the totals are exact to the
    check's tolerance and no rule has more points than the usual one."""
    mesh = cc.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 2, 2, 2, cell_type=cell_type)
    volume, surface = lens_parts(mesh)
    usual = np.asarray(volume.quadrature(order=7, mode="cut_only", options=options(0)).weights)
    sliced = np.asarray(volume.quadrature(order=7, mode="cut_only", options=options(1)).weights)
    assert (sliced > 0).all()
    assert len(sliced) < len(usual)
    assert abs(sliced.sum() - lens_volume()) < 2e-9 * lens_volume()
    if cell_type == "pyramid":
        return  # surfaces on pyramids keep the usual rule
    usual = np.asarray(surface.quadrature(order=7, mode="cut_only", options=options(0)).weights)
    sliced = np.asarray(surface.quadrature(order=7, mode="cut_only", options=options(1)).weights)
    assert (sliced > 0).all()
    assert len(sliced) < len(usual)
    assert abs(sliced.sum() - cap_area()) < 2e-9 * cap_area()


def test_surface_by_slices():
    """The cap of the lens in one hexahedron: curves in the slices, weighted to the surface's measure."""
    volume, surface = lens_parts(hexahedron())
    found, r = surface.quadrature_by_slices(0, 5)
    w = np.asarray(r.weights)
    assert found and (w > 0).all()
    assert abs(w.sum() - cap_area()) < 1e-9 * cap_area()
    X = np.asarray(r.points).reshape(len(w), 3)
    ref = surface.quadrature(order=10, mode="cut_only")
    wr = np.asarray(ref.weights)
    assert np.abs(moments(X, w) - moments(np.asarray(ref.points).reshape(len(wr), 3), wr)).max() < 1e-9
    assert len(w) < len(np.asarray(surface.quadrature(order=7, mode="cut_only", options=options(0)).weights)) / 2


def test_one_level_set_by_slices():
    """A single zero set has events too: where the slices begin and end on it."""
    mesh = hexahedron()
    part = cc.cut(mesh, [cc.create_level_set(mesh, ball(CENTRE, RADIUS), degree=2, name="a")])["a < 0"]
    found, r = part.quadrature_by_slices(0, 5)
    w = np.asarray(r.weights)
    assert found and (w > 0).all()
    assert abs(w.sum() - 4 / 3 * math.pi * RADIUS**3) < 1e-9


def test_sweep_events_see_the_faces_of_a_tetrahedron():
    """A ball that leaves a tetrahedron through a face: the curve on that face has events, and so do the
    tetrahedron's vertices that lie in the ball."""
    v = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    mesh = cc.MeshView(v, np.arange(4, dtype=np.int32), np.array([0, 4], np.int32), np.array([10], np.int32), tdim=3)
    part = cc.cut(mesh, [cc.create_level_set(mesh, ball(np.array([0.3, 0.3, -0.1]), 0.35), degree=2, name="a")])["a < 0"]
    for axis in range(3):
        found, events = part.sweep_events(0, axis)
        assert found and len(events) >= 1
    found, r = part.quadrature_by_slices(0, 5)
    w = np.asarray(r.weights)
    ref = np.asarray(part.quadrature(order=10, mode="cut_only").weights)
    assert found and (w > 0).all()
    assert abs(w.sum() - ref.sum()) < 1e-9 * ref.sum()


def analytic_lens(mesh, signed_distance):
    a = cc.analytic_sphere(list(CENTRE - 0.5 * APART * AXIS), RADIUS, signed_distance=signed_distance, name="a")
    b = cc.analytic_sphere(list(CENTRE + 0.5 * APART * AXIS), RADIUS, signed_distance=signed_distance, name="b")
    result = cc.cut(mesh, [a, b])
    return result["a < 0 and b < 0"], result["a = 0 and b < 0"]


@pytest.mark.parametrize("signed_distance", [False, True])
def test_analytic_level_sets_by_slices(signed_distance):
    """Analytic level sets are read on the slices, and their bounds over boxes find the events."""
    volume, surface = analytic_lens(hexahedron(), signed_distance)
    found, r = volume.quadrature_by_slices(0, 5)
    w = np.asarray(r.weights)
    assert found and (w > 0).all()
    assert abs(w.sum() - lens_volume()) < 2e-9 * lens_volume()
    found, r = surface.quadrature_by_slices(0, 5)
    w = np.asarray(r.weights)
    assert found and (w > 0).all()
    assert abs(w.sum() - cap_area()) < 2e-9 * cap_area()
    # the events are those of the same spheres as polynomials
    bernstein = lens(hexahedron())
    for axis in range(3):
        assert np.asarray(volume.sweep_events(0, axis)[1]) == pytest.approx(np.asarray(bernstein.sweep_events(0, axis)[1]), abs=2e-6)


@pytest.mark.parametrize("signed_distance", [False, True])
def test_slicing_option_with_analytic_level_sets(signed_distance):
    volume, surface = analytic_lens(cc.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 2, 2, 2, cell_type="tetrahedron"),
                                    signed_distance)
    usual = np.asarray(volume.quadrature(order=7, mode="cut_only", options=options(0)).weights)
    sliced = np.asarray(volume.quadrature(order=7, mode="cut_only", options=options(1)).weights)
    assert (sliced > 0).all()
    assert len(sliced) < len(usual)
    assert abs(sliced.sum() - lens_volume()) < 2e-9 * lens_volume()
    usual = np.asarray(surface.quadrature(order=7, mode="cut_only", options=options(0)).weights)
    sliced = np.asarray(surface.quadrature(order=7, mode="cut_only", options=options(1)).weights)
    assert (sliced > 0).all()
    assert len(sliced) <= len(usual)
    assert abs(sliced.sum() - cap_area()) < 1e-8 * cap_area()
