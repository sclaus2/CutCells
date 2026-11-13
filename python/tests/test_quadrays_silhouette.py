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
    options.anchor_folds = 2
    options.moment_tolerance = 0
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


def measures(result, selection, order=8):
    """The total weight of a part over the mesh, its smallest weight, and the most
    points in one cell."""
    rules = result[selection].quadrature(order=order, mode="full")
    weights = np.asarray(rules.weights)
    _, sizes = statistics(result[selection], order=order)
    return weights.sum(), weights.min(), sizes.max()


def cylinder(centre, radius, axis):
    def phi(x):
        return sum((x[i] - centre[i]) ** 2 for i in range(3) if i != axis) - radius**2

    return phi


CENTRE = (0.0123, -0.0371, 0.0217)


@pytest.mark.parametrize("n", [1, 3, 5, 9])
def test_lens(n):
    """Two balls: the curve where they meet, their silhouettes and traces in the base
    of a box, which needs no bisection however large the cell."""
    radius, distance = 0.6, 0.5
    axis = np.array([1.0, 2.0, 2.0]) / 3
    first = np.array([-0.21, 0.07, 0.04])
    mesh = box_mesh(n)
    a = cc.create_level_set(mesh, sphere(first, radius), degree=2, name="a")
    b = cc.create_level_set(mesh, sphere(first + distance * axis, radius), degree=2, name="b")
    result = cc.cut(mesh, [a, b])

    total, smallest, worst = measures(result, "a < 0 and b < 0")
    assert abs(total - math.pi * (4 * radius + distance) * (2 * radius - distance) ** 2 / 12) < 1e-11
    assert smallest > 0
    assert worst <= {1: 45000, 3: 30000, 5: 25000, 9: 15000}[n]

    total, smallest, worst = measures(result, "a = 0 and b < 0")
    assert abs(total - 2 * math.pi * radius * (radius - distance / 2)) < 1e-10
    assert smallest > 0
    assert worst <= 6000


@pytest.mark.parametrize("n", [1, 3, 5, 9])
def test_ball_and_plane(n):
    radius, offset = 0.6, 0.05
    normal = np.array([0.93, 0.33, -0.19])
    normal /= np.linalg.norm(normal)
    centre = np.array([-0.21, 0.07, 0.04])
    mesh = box_mesh(n)
    a = cc.create_level_set(mesh, sphere(centre, radius), degree=2, name="a")
    b = cc.create_level_set(mesh, lambda x: sum(normal[i] * x[i] for i in range(3)) - offset, degree=1, name="b")
    result = cc.cut(mesh, [a, b])
    # the height of the cap on the side where b < 0
    height = radius + offset - normal @ centre

    total, smallest, worst = measures(result, "a < 0 and b < 0")
    assert abs(total - math.pi * height**2 * (3 * radius - height) / 3) < 1e-11
    assert smallest > 0
    assert worst <= {1: 25000, 3: 10000, 5: 15000, 9: 12000}[n]

    total, smallest, worst = measures(result, "a = 0 and b < 0")
    assert abs(total - 2 * math.pi * radius * height) < 1e-11
    assert smallest > 0
    assert worst <= 8000


@pytest.mark.parametrize("n", [3, 5, 9])
def test_two_cylinders(n):
    """Their zero curves on faces are straight, some along the height of the base. On
    the interface, the cell where the cylinders touch is less accurate."""
    radius = 0.6
    mesh = box_mesh(n)
    a = cc.create_level_set(mesh, cylinder(CENTRE, radius, 2), degree=2, name="a")
    b = cc.create_level_set(mesh, cylinder(CENTRE, radius, 1), degree=2, name="b")
    result = cc.cut(mesh, [a, b])

    total, smallest, worst = measures(result, "a < 0 and b < 0")
    assert abs(total - 16 / 3 * radius**3) < 1e-11
    assert smallest > 0
    assert worst <= 20000

    total, smallest, worst = measures(result, "a = 0 and b < 0")
    assert abs(total - 8 * radius**2) < 1e-8
    assert smallest > 0
    assert worst <= 2000


@pytest.mark.parametrize("n", [3, 5, 9])
def test_three_cylinders(n):
    radius = 0.6
    mesh = box_mesh(n)
    a = cc.create_level_set(mesh, cylinder(CENTRE, radius, 2), degree=2, name="a")
    b = cc.create_level_set(mesh, cylinder(CENTRE, radius, 1), degree=2, name="b")
    c = cc.create_level_set(mesh, cylinder(CENTRE, radius, 0), degree=2, name="c")
    result = cc.cut(mesh, [a, b, c])

    total, smallest, worst = measures(result, "a < 0 and b < 0 and c < 0")
    assert abs(total - 8 * (2 - math.sqrt(2)) * radius**3) < 1e-11
    assert smallest > 0
    assert worst <= 20000

    total, smallest, worst = measures(result, "a = 0 and b < 0 and c < 0")
    assert abs(total - 8 * (2 - math.sqrt(2)) * radius**2) < 1e-7
    assert smallest > 0
    assert worst <= 2000


def test_lens_integrates_shape_functions():
    """What assembly sees: against the rule without silhouette splits, cell by cell."""
    radius, distance = 0.6, 0.5
    axis = np.array([1.0, 2.0, 2.0]) / 3
    first = np.array([-0.21, 0.07, 0.04])
    mesh = box_mesh(5)
    a = cc.create_level_set(mesh, sphere(first, radius), degree=2, name="a")
    b = cc.create_level_set(mesh, sphere(first + distance * axis, radius), degree=2, name="b")
    result = cc.cut(mesh, [a, b])
    options = cc.QuadraysOptions()
    options.silhouette_splits = False
    options.max_points = -1
    options.anchor_folds = 2
    options.moment_tolerance = 0
    for selection in ("a < 0 and b < 0", "a = 0 and b < 0"):
        rules, _ = statistics(result[selection], order=8)
        reference, _ = statistics(result[selection], order=8, options=options)
        assert np.array_equal(np.asarray(rules.parent_map), np.asarray(reference.parent_map))
        for mine, theirs in zip(np.split(np.arange(len(rules.weights)), np.asarray(rules.offset)[1:-1]),
                                np.split(np.arange(len(reference.weights)), np.asarray(reference.offset)[1:-1])):
            x = np.asarray(rules.points).reshape(-1, 3)[mine]
            y = np.asarray(reference.points).reshape(-1, 3)[theirs]
            assert abs(np.asarray(rules.weights)[mine] @ legendre_squares(x)
                       - np.asarray(reference.weights)[theirs] @ legendre_squares(y)) < 1e-7


@pytest.mark.parametrize("gap", [1e-3, 1e-5, 1e-7, 1e-9])
def test_small_features_are_not_lost(gap):
    """A closed curve smaller than the grid of a base's search shows in no sign change,
    and no point of a rule lands on it for the accuracy check to see."""
    mesh = box_mesh(4)
    distance = 2 * 0.3 - gap
    a = cc.create_level_set(mesh, sphere((0.13, -0.03 - distance / 2, 0.21), 0.3), degree=2, name="a")
    b = cc.create_level_set(mesh, sphere((0.13, -0.03 + distance / 2, 0.21), 0.3), degree=2, name="b")
    result = cc.cut(mesh, [a, b])
    total, _, _ = measures(result, "a = 0 and b < 0")
    assert abs(total - 2 * math.pi * 0.3 * (0.3 - distance / 2)) < 2e-12

    c = cc.create_level_set(mesh, sphere((0.1, 0.07, 0.02), 0.4), degree=2, name="c")
    d = cc.create_level_set(mesh, lambda x: 0.42 - gap - x[2], degree=1, name="d")
    result = cc.cut(mesh, [c, d])
    total, _, _ = measures(result, "c = 0 and d < 0")
    assert abs(total - 2 * math.pi * 0.4 * gap) < 2e-12


@pytest.mark.parametrize("gap", [0.0, 1e-13, 1e-9, 1e-5])
def test_events_on_the_borders_of_cells(gap):
    """Folds and crossings exactly on, or just off, faces and edges of the mesh."""
    radius, distance = 0.6, 0.5
    rho = math.sqrt(radius**2 - distance**2 / 4)
    exact_volume = math.pi * (4 * radius + distance) * (2 * radius - distance) ** 2 / 12
    exact_area = 2 * math.pi * radius * (radius - distance / 2)
    z = 0.11
    centres = {"in a face": (gap, 0.07, 0.04),
               "touching a face": (0.13, 0.5 - rho + gap, 0.04),
               "through an edge": (0.13, 0.5 - math.sqrt(rho**2 - (0.5 - z) ** 2) + gap, z)}
    mesh = box_mesh(4)
    for x, y, z0 in centres.values():
        a = cc.create_level_set(mesh, sphere((x - distance / 2, y, z0), radius), degree=2, name="a")
        b = cc.create_level_set(mesh, sphere((x + distance / 2, y, z0), radius), degree=2, name="b")
        result = cc.cut(mesh, [a, b])
        total, smallest, worst = measures(result, "a < 0 and b < 0")
        assert abs(total - exact_volume) < 5e-11
        assert smallest > 0
        assert worst <= 15000
        total, smallest, worst = measures(result, "a = 0 and b < 0")
        assert abs(total - exact_area) < 5e-11
        assert smallest > 0
        assert worst <= 2000


@pytest.mark.parametrize("n", [3, 5, 9])
def test_lens_on_tetrahedra(n):
    """A tetrahedron is a box clipped by a plane, which bounds the lines like a face."""
    radius, distance = 0.6, 0.5
    axis = np.array([1.0, 2.0, 2.0]) / 3
    first = np.array([-0.21, 0.07, 0.04])
    mesh = cc.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n + 1, n + 1, n + 1, cell_type="tetrahedron")
    a = cc.create_level_set(mesh, sphere(first, radius), degree=2, name="a")
    b = cc.create_level_set(mesh, sphere(first + distance * axis, radius), degree=2, name="b")
    result = cc.cut(mesh, [a, b])
    total, smallest, worst = measures(result, "a < 0 and b < 0")
    assert abs(total - math.pi * (4 * radius + distance) * (2 * radius - distance) ** 2 / 12) < 2e-11
    assert smallest > 0
    assert worst <= {3: 40000, 5: 30000, 9: 15000}[n]
    total, smallest, worst = measures(result, "a = 0 and b < 0")
    assert abs(total - 2 * math.pi * radius * (radius - distance / 2)) < 2e-11
    assert smallest > 0
    assert worst <= 2000


@pytest.mark.parametrize("signed_distance, limits", [(False, {3: 30000, 5: 25000, 9: 15000}),
                                                     (True, {3: 70000, 5: 30000, 9: 20000})])
@pytest.mark.parametrize("n", [3, 5, 9])
def test_analytic_lens(n, signed_distance, limits):
    """Analytic level sets split at silhouettes too, their derivatives read through
    their second-derivative bounds: |x - c|^2 - r^2 and the distance |x - c| - r."""
    radius, distance = 0.6, 0.5
    axis = np.array([1.0, 2.0, 2.0]) / 3
    first = np.array([-0.21, 0.07, 0.04])
    a = cc.analytic_sphere(list(first), radius, signed_distance=signed_distance, name="a")
    b = cc.analytic_sphere(list(first + distance * axis), radius, signed_distance=signed_distance, name="b")
    result = cc.cut(box_mesh(n), [a, b])
    total, smallest, worst = measures(result, "a < 0 and b < 0")
    assert abs(total - math.pi * (4 * radius + distance) * (2 * radius - distance) ** 2 / 12) < 2e-11
    assert smallest > 0
    assert worst <= limits[n]
    total, smallest, worst = measures(result, "a = 0 and b < 0")
    assert abs(total - 2 * math.pi * radius * (radius - distance / 2)) < 2e-11
    assert smallest > 0
    assert worst <= 3000


def test_anchored_points_follow_what_assembly_sees():
    """One anchored factor takes far fewer points, and the measure alone does not see
    what that costs a polynomial. The check on a moment of the mass matrix does, and
    gives the boxes that need it two factors."""
    radius, distance = 0.6, 0.5
    axis = np.array([1.0, 2.0, 2.0]) / 3
    first = np.array([-0.21, 0.07, 0.04])
    mesh = box_mesh(5)
    a = cc.create_level_set(mesh, sphere(first, radius), degree=2, name="a")
    b = cc.create_level_set(mesh, sphere(first + distance * axis, radius), degree=2, name="b")
    part = cc.cut(mesh, [a, b])["a < 0 and b < 0"]

    def per_cell(options):
        rules = part.quadrature(order=8, mode="cut_only", options=options)
        points = np.asarray(rules.points).reshape(-1, 3)
        weights = np.asarray(rules.weights)
        offset = np.asarray(rules.offset)
        moments = {int(c): weights[offset[i]:offset[i + 1]] @ legendre_squares(points[offset[i]:offset[i + 1]])
                   for i, c in enumerate(np.asarray(rules.parent_map))}
        return moments, np.diff(offset).mean(), weights.min()

    fine = cc.QuadraysOptions()
    fine.points_per_segment = 12
    fine.max_points = -1
    fine.anchor_folds = 2
    fine.moment_tolerance = 0
    reference, _, _ = per_cell(fine)

    def l1(moments):
        return sum(abs(moments.get(c, 0.0) - reference[c]) for c in reference)

    two_factors = cc.QuadraysOptions()
    two_factors.anchor_folds = 2
    two_factors.moment_tolerance = 0
    heavy, heavy_mean, _ = per_cell(two_factors)

    light = cc.QuadraysOptions()
    light.moment_tolerance = 0
    unchecked, _, _ = per_cell(light)

    checked, checked_mean, smallest = per_cell(cc.QuadraysOptions())

    assert cc.QuadraysOptions().anchor_folds == 1
    assert cc.QuadraysOptions().moment_tolerance == 1e-8
    assert smallest > 0
    assert checked_mean < 0.75 * heavy_mean
    assert l1(checked) < 1e-8
    assert l1(checked) < 0.1 * l1(unchecked)
    assert l1(heavy) < 1e-9


@pytest.mark.parametrize("kind", ["hexahedron", "tetrahedron"])
def test_fine_meshes_take_no_more_points_than_two_factors(kind):
    """On a resolved sphere most boxes pass with one factor or have no anchored segment:
    the default is about the rule with two factors at worst."""
    mesh = cc.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 18, 18, 18, cell_type=kind)
    phi = cc.create_level_set(mesh, sphere((0.1, -0.05, 0.03), 0.6), degree=2, name="phi")
    part = cc.cut(mesh, phi)["phi = 0"]
    two_factors = cc.QuadraysOptions()
    two_factors.anchor_folds = 2
    two_factors.moment_tolerance = 0
    default = part.quadrature(order=6, mode="cut_only")
    heavy = part.quadrature(order=6, mode="cut_only", options=two_factors)
    assert len(default.weights) <= 1.1 * len(heavy.weights)
    assert np.sum(default.weights) == pytest.approx(4 * math.pi * 0.6**2, abs=1e-9)


def test_slow_cells_are_written_for_plotting(tmp_path):
    """debug_seconds writes a cell for its time, whatever its points."""
    debug = pytest.importorskip("cutcells.quadrays_debug")
    mesh = box_mesh(3)
    phi = cc.create_level_set(mesh, sphere((-0.21, 0.07, 0.04), 0.5), degree=2, name="phi")
    options = cc.QuadraysOptions()
    options.debug_dir = str(tmp_path)
    options.debug_points = 10**9
    options.debug_seconds = 0.0
    options.debug_reports = 2
    cc.cut(mesh, phi)["phi < 0"].quadrature(order=8, mode="cut_only", options=options)
    files = sorted(tmp_path.glob("*.json"))
    assert len(files) == 2
    report = debug.read_report(str(files[0]))
    assert report["seconds"] > 0
    assert "ms" in debug.summary(report)


def test_splits_can_be_turned_off(ball_in_cell):
    options = cc.QuadraysOptions()
    assert options.silhouette_splits
    options.silhouette_splits = False
    options.max_points = -1
    _, sizes = statistics(ball_in_cell["phi < 0"], order=8, options=options)
    assert sizes.max() > 10000
