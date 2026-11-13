# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""The accuracy check's moment on interfaces (QuadraysOptions.moment_tolerance).

q points are exact to degree 2 q - 1 on a flat piece. Curvature takes a Legendre polynomial squared along the cell's
axes past that, on any box, also where no segment is anchored, and a box then takes more points per direction before
it is bisected. Errors are per cell against rules with 14 points per segment, since totals cancel across shared faces.
"""

import numpy as np
import pytest

import cutcells as cc

F = np.array((-0.21, 0.07, 0.04))
NORMAL = np.array((0.93, 0.33, -0.19))
R = 0.6
ORDER = 10
DEGREE = ORDER // 2


def ball(x):
    return sum((x[i] - F[i]) ** 2 for i in range(3)) - R**2


def plane(x):
    return sum(NORMAL[i] * x[i] for i in range(3)) - 0.05


def part_of(kind, n, phi=ball, degree=2):
    mesh = cc.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n + 1, n + 1, n + 1, cell_type=kind)
    return cc.cut(mesh, [cc.create_level_set(mesh, phi, degree=degree, name="a")])["a = 0"]


def moments(rules):
    points = np.asarray(rules.points).reshape(-1, 3)
    weights = np.asarray(rules.weights)
    offset = np.asarray(rules.offset)
    legendre = np.polynomial.legendre.Legendre.basis(DEGREE)
    shape = sum(legendre(2 * points[:, j] - 1) ** 2 for j in range(3))
    return {int(c): weights[offset[i] : offset[i + 1]] @ shape[offset[i] : offset[i + 1]] for i, c in enumerate(np.asarray(rules.parent_map))}


def reference(part):
    options = cc.QuadraysOptions()
    options.anchor_folds = 2
    options.moment_tolerance = 0
    options.points_per_segment = 14
    options.max_points = -1
    return moments(part.quadrature(order=8, mode="cut_only", options=options))


def measured(part, exact, **values):
    options = cc.QuadraysOptions()
    for name, value in values.items():
        setattr(options, name, value)
    rules = part.quadrature(order=ORDER, mode="cut_only", options=options)
    got = moments(rules)
    assert set(got) == set(exact)
    return sum(abs(got[c] - exact[c]) for c in exact), np.asarray(rules.weights).size / len(exact)


@pytest.mark.parametrize("kind", ["hexahedron", "tetrahedron"])
def test_a_flat_piece_needs_no_more_points(kind):
    part = part_of(kind, 9, plane, 1)
    error, mean = measured(part, reference(part))
    assert error < 1e-13
    assert mean <= 2 * 6**2


def test_the_check_reaches_boxes_without_anchored_segments():
    part = part_of("tetrahedron", 9)
    exact = reference(part)
    checked, _ = measured(part, exact)
    unchecked, _ = measured(part, exact, moment_tolerance=0.0)
    # at this size most boxes on tetrahedra have no anchored segment
    assert checked < 1e-7
    assert unchecked > 3 * checked


@pytest.mark.parametrize("kind, n, largest_error, largest_mean", [("hexahedron", 17, 2e-8, 170), ("tetrahedron", 17, 1e-7, 140)])
def test_more_points_before_bisection(kind, n, largest_error, largest_mean):
    # Bisecting took 183 points per cut cell on hexahedra at n = 17, for 3e-8, and 146 on tetrahedra, for 7e-8.
    part = part_of(kind, n)
    error, mean = measured(part, reference(part))
    assert error < largest_error
    assert mean < largest_mean


def test_a_tighter_tolerance_costs_few_points():
    part = part_of("hexahedron", 9)
    exact = reference(part)
    error, mean = measured(part, exact)
    tighter_error, tighter_mean = measured(part, exact, moment_tolerance=1e-12)
    assert tighter_error < error / 100
    assert tighter_mean < 1.4 * mean


def test_a_box_keeps_its_rules_for_the_next_cell():
    # the rules a box raises for itself are back to the run's own when the next cell starts, in any order of calls
    part = part_of("hexahedron", 5)
    options = cc.QuadraysOptions()
    first = part.quadrature(order=ORDER, mode="cut_only", options=options)
    again = part.quadrature(order=ORDER, mode="cut_only", options=options)
    lower = part.quadrature(order=2, mode="cut_only", options=options)
    plain = cc.QuadraysOptions()
    plain.moment_tolerance = 0
    plain_lower = part.quadrature(order=2, mode="cut_only", options=plain)
    assert np.array_equal(np.asarray(first.weights), np.asarray(again.weights))
    assert np.asarray(lower.weights).size >= np.asarray(plain_lower.weights).size
