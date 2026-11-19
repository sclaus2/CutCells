# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Level sets that are constant along an axis of the cell (cylinders and other extruded
shapes): along that axis the height lines take no root of them, and the volume rule runs
its lines there. Exact volumes, and few points."""

import math

import numpy as np
import pytest

import cutcells as cc

CENTRE = (0.0123, -0.0371, 0.0217)


def box_mesh(n, cell_type="hexahedron"):
    """A grid of n cells per direction on [-1, 1]^3 (six tetrahedra per cell)."""
    return cc.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n + 1, n + 1, n + 1, cell_type=cell_type)


def cylinder(centre, radius, axis):
    def phi(x):
        return sum((x[i] - centre[i]) ** 2 for i in range(3) if i != axis) - radius**2

    return phi


def sphere(centre, radius):
    def phi(x):
        return sum((x[i] - centre[i]) ** 2 for i in range(3)) - radius**2

    return phi


def volume_rule(mesh, a, b, **kw):
    part = cc.cut(mesh, [a, b])["a < 0 and b < 0"]
    return np.asarray(part.quadrature(order=6, mode="full", **kw).weights)


@pytest.mark.parametrize("first,second", [(2, 1), (1, 2), (0, 2)])
@pytest.mark.parametrize("cell_type,n,most", [("hexahedron", 1, 1000), ("hexahedron", 2, 5000), ("tetrahedron", 1, 8000)])
def test_two_crossing_cylinders(first, second, cell_type, n, most):
    """Two cylinders of equal radius whose axes meet: 16 r^3 / 3, whichever is listed first."""
    mesh = box_mesh(n, cell_type)
    a = cc.create_level_set(mesh, cylinder(CENTRE, 0.6, first), degree=2, name="a")
    b = cc.create_level_set(mesh, cylinder(CENTRE, 0.6, second), degree=2, name="b")
    weights = volume_rule(mesh, a, b)
    assert weights.min() > 0
    assert abs(weights.sum() - 16 * 0.6**3 / 3) < 1e-9
    assert len(weights) <= most


@pytest.mark.parametrize("cylinder_first", [True, False])
@pytest.mark.parametrize("cell_type,n,most", [("hexahedron", 1, 3000), ("hexahedron", 2, 6000), ("tetrahedron", 2, 20000)])
def test_cylinder_through_a_ball(cylinder_first, cell_type, n, most):
    """A cylinder through the centre of a ball: 4 pi / 3 (R^3 - (R^2 - r^2)^(3/2))."""
    mesh = box_mesh(n, cell_type)
    tube = cylinder(CENTRE, 0.25, 2)
    ball = sphere(CENTRE, 0.6)
    shapes = (tube, ball) if cylinder_first else (ball, tube)
    a = cc.create_level_set(mesh, shapes[0], degree=2, name="a")
    b = cc.create_level_set(mesh, shapes[1], degree=2, name="b")
    weights = volume_rule(mesh, a, b)
    assert weights.min() > 0
    assert abs(weights.sum() - 4 * math.pi / 3 * (0.6**3 - (0.6**2 - 0.25**2) ** 1.5)) < 1e-9
    assert len(weights) <= most


@pytest.mark.parametrize("height,base", [(1, 2), (2, 1)])
def test_lines_along_a_cylinder_settle_at_once(height, base):
    """Height lines along one cylinder's axis and base lines along the other's: the wall's
    outline splits the base, and nothing is bisected (375 points; 9,300 before along x)."""
    mesh = box_mesh(1)
    a = cc.create_level_set(mesh, cylinder(CENTRE, 0.6, 2), degree=2, name="a")
    b = cc.create_level_set(mesh, cylinder(CENTRE, 0.6, 1), degree=2, name="b")
    options = cc.QuadraysOptions()
    options.height_axis, options.base_axis = height, base
    weights = volume_rule(mesh, a, b, options=options)
    assert abs(weights.sum() - 16 * 0.6**3 / 3) < 1e-11
    assert len(weights) <= 500
