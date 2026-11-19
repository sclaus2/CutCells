# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Two level sets that touch each other where one of them turns back along the height lines
of a cell (on its silhouette): cylinders of one radius whose axes meet touch in two points,
and the curve in which they meet crosses itself there. The base lines must split at such a
point although no function changes its sign at it."""

import math

import numpy as np
import pytest

import cutcells as cc
from cutcells import mesh_utils

CENTRE = (0.0123, -0.0371, 0.0217)
RADIUS = 0.6


def two_axes(axis):
    """A cylinder as the sum of the squares over the two other axes."""
    return lambda x: sum((x[i] - CENTRE[i]) ** 2 for i in range(3) if i != axis) - RADIUS**2


def three_less_one(axis):
    """The same polynomial, written so that its nodal values differ by rounding."""
    return lambda x: sum((x[i] - CENTRE[i]) ** 2 for i in range(3)) - (x[axis] - CENTRE[axis]) ** 2 - RADIUS * RADIUS


def cylinders(mesh, form):
    return [cc.create_level_set(mesh, form(axis), degree=2, name=name) for name, axis in (("a", 2), ("b", 1), ("c", 0))]


SURFACE = "a = 0 and b < 0 and c < 0"


@pytest.mark.parametrize("form", [two_axes, three_less_one])
@pytest.mark.parametrize("order", [6, 8])
def test_surface_of_three_cylinders_on_six_tetrahedra(form, order):
    """The piece of one of three crossing cylinders inside the other two: (16 - 8 sqrt 2) r^2,
    in either rounding of the polynomials (which decides the base direction of one cell)."""
    mesh = cc.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 2, 2, 2, cell_type="tetrahedron")
    weights = np.asarray(cc.cut(mesh, cylinders(mesh, form))[SURFACE].quadrature(order=order, mode="full").weights)
    assert weights.min() > 0
    assert abs(weights.sum() - (16 - 8 * math.sqrt(2)) * RADIUS**2) < 1e-7


def test_every_order_of_directions_in_the_touching_cell():
    """The tetrahedron that holds a touching point 0.0047 of its size from a face: the same area
    with every order of the height directions that the surface allows. With the height lines
    along x and the base lines along y, the base lines' ends pass the touching point."""
    vertices = np.array([[-1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])
    mesh = mesh_utils._mesh_view(vertices, np.array([[0, 1, 2, 3]]), 10, 3)
    part = cc.cut(mesh, cylinders(mesh, two_axes))[SURFACE]
    areas = []
    for height, base in ((0, 1), (0, 2), (1, 0), (1, 2)):
        options = cc.QuadraysOptions()
        options.height_axis, options.base_axis = height, base
        weights = np.asarray(part.quadrature(order=6, mode="full", options=options).weights)
        assert weights.min() > 0
        areas.append(weights.sum())
    assert max(areas) - min(areas) < 1e-7
    assert abs(areas[0] - 0.39774676) < 1e-7
