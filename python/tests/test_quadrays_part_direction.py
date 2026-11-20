# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""The height direction of an interface part is chosen where the part lies. A ball may fill
a cell in which the part is a cap that a cylinder cuts out of it; where the cap's rim clips
a corner of a neighbouring cell, lines that cross the ball well on average there run along
the clipped piece, which is then too thin in the base to be found."""

import math

import numpy as np
import pytest

import cutcells as cc

CENTRE = np.array([0.0123, -0.0371, 0.0217])
BALL, TUBE = 0.6, 0.25


def rotation(axis, degrees):
    a = np.asarray(axis, dtype=float)
    a /= np.linalg.norm(a)
    t = math.radians(degrees)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + math.sin(t) * K + (1 - math.cos(t)) * (K @ K)


def cylinder(direction, radius):
    """A cylinder about the line through CENTRE along a unit vector: a quadratic polynomial."""
    return lambda x: (
        sum((x[i] - CENTRE[i]) ** 2 for i in range(3)) - sum(direction[i] * (x[i] - CENTRE[i]) for i in range(3)) ** 2 - radius**2
    )


def ball(radius):
    return lambda x: sum((x[i] - CENTRE[i]) ** 2 for i in range(3)) - radius**2


@pytest.mark.parametrize("cell_type", ["hexahedron", "tetrahedron"])
@pytest.mark.parametrize("order", [6, 8])
def test_caps_of_a_ball_inside_a_tilted_cylinder(cell_type, order):
    """The two caps that a cylinder through the centre of a ball cuts out of its surface:
    4 pi R (R - sqrt(R^2 - r^2)). On 2 cells per direction one rim clips a corner of a cell,
    a piece of area 6.5e-7 that is 0.00075 by 0.0017 of the cell."""
    mesh = cc.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 3, 3, 3, cell_type=cell_type)
    along = rotation((1, 2, 2), 30)[:, 2]
    a = cc.create_level_set(mesh, cylinder(along, TUBE), degree=2, name="a")
    b = cc.create_level_set(mesh, ball(BALL), degree=2, name="b")
    weights = np.asarray(cc.cut(mesh, [a, b])["b = 0 and a < 0"].quadrature(order=order, mode="full").weights)
    exact = 4 * math.pi * BALL * (BALL - math.sqrt(BALL**2 - TUBE**2))
    assert weights.min() > 0
    assert abs(weights.sum() - exact) < 1e-7
