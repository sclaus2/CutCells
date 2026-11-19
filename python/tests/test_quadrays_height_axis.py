# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""QuadraysOptions.height_axis and base_axis: a cell of dimension 3 integrated with given
height directions, for studying what the choice of directions costs."""

import math

import numpy as np
import pytest

import cutcells as cc

ORDERS = [(k, j) for k in range(3) for j in range(3) if j != k]


def box_mesh(n):
    """A grid of n hexahedra per direction on [-1, 1]^3."""
    return cc.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n + 1, n + 1, n + 1, cell_type="hexahedron")


def sphere(centre, radius):
    def phi(x):
        return sum((x[i] - centre[i]) ** 2 for i in range(3)) - radius**2

    return phi


def distinct(points):
    """Distinct coordinates of the points along each axis."""
    return [len(np.unique(np.round(points[:, j], 12))) for j in range(3)]


def test_default_leaves_the_choice_to_the_engine():
    options = cc.QuadraysOptions()
    assert options.height_axis == -1 and options.base_axis == -1


@pytest.mark.parametrize("height,base", ORDERS)
def test_ball_in_one_cell_along_given_axes(height, base):
    mesh = box_mesh(1)
    a = cc.create_level_set(mesh, sphere((-0.21, 0.07, 0.04), 0.6), degree=2, name="a")
    part = cc.cut(mesh, [a])["a < 0"]
    options = cc.QuadraysOptions()
    options.height_axis, options.base_axis = height, base
    rules = part.quadrature(order=6, mode="cut_only", options=options)
    weights = np.asarray(rules.weights)
    points = np.asarray(rules.points).reshape(-1, 3)
    assert weights.min() > 0
    assert abs(weights.sum() - 4 / 3 * math.pi * 0.6**3) < 1e-11
    # the height lines run along height_axis, the lines of the base along base_axis: the
    # third axis carries the outer nodes, the fewest distinct coordinates
    counts = distinct(points)
    assert counts[3 - height - base] < counts[base] < counts[height]


@pytest.mark.parametrize("height,base", ORDERS)
def test_lens_on_a_mesh_along_given_axes(height, base):
    """Two level sets: every order gives the volume of the engine's own choice."""
    mesh = box_mesh(4)
    a = cc.create_level_set(mesh, sphere((-0.21, 0.07, 0.04), 0.6), degree=2, name="a")
    b = cc.create_level_set(mesh, sphere((-0.21 + 0.5 / 3, 0.07 + 1.0 / 3, 0.04 + 1.0 / 3), 0.6), degree=2, name="b")
    part = cc.cut(mesh, [a, b])["a < 0 and b < 0"]
    free = np.asarray(part.quadrature(order=6, mode="cut_only").weights)
    options = cc.QuadraysOptions()
    options.height_axis, options.base_axis = height, base
    forced = np.asarray(part.quadrature(order=6, mode="cut_only", options=options).weights)
    assert forced.min() > 0
    assert abs(forced.sum() - free.sum()) < 1e-10
