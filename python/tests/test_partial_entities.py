# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Pieces on mesh facets and edges that another level set cuts.

A zero set can lie on facets or edges of the mesh, a grid plane say. Where a
second level set cuts every cell sharing such an entity, the part holds only
a piece of the entity, which is integrated on the entity itself.
"""

import numpy as np
import pytest

import cutcells

RADIUS = 0.7


def box(cell_type, n=13):
    return cutcells.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n, n, n, cell_type=cell_type)


@pytest.mark.parametrize("cell_type", ["tetrahedron", "hexahedron"])
@pytest.mark.parametrize("backend", ["quadrays", "lut"])
def test_disk_on_grid_facets(cell_type, backend):
    # the plane y = 0 holds grid facets; the ball cuts some of them
    mesh = box(cell_type)
    ball = cutcells.analytic_sphere([0.0, 0.0, 0.0], RADIUS, name="ball")
    plane = cutcells.create_level_set(mesh, lambda X: X[1], degree=1, name="plane")
    result = cutcells.cut(mesh, [ball, plane], backend=backend)
    disk = result[(plane == 0) & (ball < 0)]
    area = disk.quadrature(order=4).weights.sum()
    assert area == pytest.approx(np.pi * RADIUS**2, rel=1e-6)
    on = disk.entity_quadrature(order=4)
    assert on.facets.weights.sum() == pytest.approx(area, rel=1e-12)


@pytest.mark.parametrize("cell_type", ["tetrahedron", "hexahedron"])
def test_disk_on_grid_facets_interpolated(cell_type):
    # with a P2 ball the facets take its polynomial restricted to them
    mesh = box(cell_type)
    ball = cutcells.create_level_set(mesh, lambda X: np.sqrt(X[0] ** 2 + X[1] ** 2 + X[2] ** 2) - RADIUS,
                                     degree=2, name="ball")
    plane = cutcells.create_level_set(mesh, lambda X: X[1], degree=1, name="plane")
    result = cutcells.cut(mesh, [ball, plane])
    area = result[(plane == 0) & (ball < 0)].quadrature(order=4).weights.sum()
    assert area == pytest.approx(np.pi * RADIUS**2, rel=1e-4)
    # the same disk just off the grid planes, integrated in the cells
    shifted = cutcells.create_level_set(mesh, lambda X: X[1] - 1e-3, degree=1, name="plane")
    off = cutcells.cut(mesh, [ball, shifted])
    off_area = off[(shifted == 0) & (ball < 0)].quadrature(order=4).weights.sum()
    assert area == pytest.approx(off_area, rel=1e-4)


@pytest.mark.parametrize("cell_type", ["triangle", "quadrilateral"])
def test_segment_on_grid_edges(cell_type):
    mesh = cutcells.rectangle_mesh(-1.0, -1.0, 1.0, 1.0, 13, 13, cell_type=cell_type)
    disk = cutcells.analytic_sphere([0.0, 0.0], RADIUS, name="disk")
    line = cutcells.create_level_set(mesh, lambda X: X[1], degree=1, name="line")
    result = cutcells.cut(mesh, [disk, line])
    length = result[(line == 0) & (disk < 0)].quadrature(order=4).weights.sum()
    assert length == pytest.approx(2.0 * RADIUS, rel=1e-10)


@pytest.mark.parametrize("cell_type", ["tetrahedron", "hexahedron"])
def test_segment_on_grid_edges_3d(cell_type):
    # the planes x = 0 and y = 0 meet along grid edges; the ball cuts some of them
    mesh = box(cell_type)
    ball = cutcells.analytic_sphere([0.0, 0.0, 0.0], RADIUS, name="ball")
    px = cutcells.create_level_set(mesh, lambda X: X[0], degree=1, name="px")
    py = cutcells.create_level_set(mesh, lambda X: X[1], degree=1, name="py")
    result = cutcells.cut(mesh, [ball, px, py])
    length = result[(px == 0) & (py == 0) & (ball < 0)].quadrature(order=4).weights.sum()
    assert length == pytest.approx(2.0 * RADIUS, rel=1e-10)


def test_whole_facets_unchanged():
    # a plane on grid facets with no other level set holds them whole
    mesh = box("tetrahedron")
    plane = cutcells.create_level_set(mesh, lambda X: X[1], degree=1, name="plane")
    result = cutcells.cut(mesh, plane)
    assert result[plane == 0].quadrature(order=4).weights.sum() == pytest.approx(4.0, rel=1e-12)
