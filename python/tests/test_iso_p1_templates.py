# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Pk-iso-P1 templates, and the template order that cut() hands to the lookup
tables through cut_approximation and cut_approximation_order."""

import numpy as np
import pytest

import cutcells


def triangle_mesh():
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float64)
    return cutcells.MeshView(coords, np.array([0, 1, 2], dtype=np.int32), np.array([0, 3], dtype=np.int32),
                             np.array([5], dtype=np.int32), tdim=2)


@pytest.mark.parametrize(
    "cell_type,order,n_vertices,n_cells,child_type",
    [
        (cutcells.CellType.interval, 4, 5, 4, cutcells.CellType.interval),
        (cutcells.CellType.triangle, 3, 10, 9, cutcells.CellType.triangle),
        (cutcells.CellType.quadrilateral, 2, 9, 4, cutcells.CellType.quadrilateral),
        (cutcells.CellType.tetrahedron, 3, 20, 27, cutcells.CellType.tetrahedron),
        (cutcells.CellType.hexahedron, 2, 27, 8, cutcells.CellType.hexahedron),
    ],
)
def test_template_counts(cell_type, order, n_vertices, n_cells, child_type):
    tpl = cutcells.iso_p1_template(cell_type, order)
    assert tpl.n_vertices == n_vertices
    assert tpl.n_cells == n_cells
    assert tpl.child_cell_type == child_type
    assert len(np.asarray(tpl.vertex_parent_dim)) == n_vertices
    assert len(np.asarray(tpl.vertex_parent_id)) == n_vertices
    assert len(np.asarray(tpl.cell_connectivity)) == n_cells * tpl.vertices_per_cell


def pieces(result, expr="phi < 0"):
    return len(np.asarray(result[expr].visualization_mesh(mode="cut_only").types))


def test_cut_approximation_sets_the_template_order():
    mesh = triangle_mesh()
    ls = cutcells.create_level_set(mesh, lambda x: x[0] + x[1] - 0.9, degree=1, name="phi")
    default = cutcells.cut(mesh, ls)
    refined = cutcells.cut(mesh, ls, cut_approximation="iso_p1", cut_approximation_order=3)
    assert refined.num_cut_cells == 1
    assert refined.options.template_order == 3
    assert pieces(refined) > pieces(default)
    # the plane is exact on any template
    for result in (default, refined):
        assert np.sum(result["phi < 0"].quadrature(order=1).weights) == pytest.approx(0.5 * 0.9**2, rel=1e-14)


def test_higher_order_level_set_uses_its_degree():
    mesh = triangle_mesh()
    ls = cutcells.create_level_set(mesh, lambda x: x[0] + x[1] - 0.9, degree=3, name="phi")
    automatic = cutcells.cut(mesh, ls)
    explicit = cutcells.cut(mesh, ls, cut_approximation="iso_p1", cut_approximation_order=3)
    linear = cutcells.cut(mesh, ls, cut_approximation="linear")
    assert automatic.options.template_order == 0 and linear.options.template_order == 1
    assert pieces(automatic) == pieces(explicit) > pieces(linear)


def test_invalid_template_orders_are_rejected():
    mesh = triangle_mesh()
    ls = cutcells.create_level_set(mesh, lambda x: x[0] + x[1] - 0.9, degree=1, name="phi")
    with pytest.raises(ValueError, match="linear"):
        cutcells.cut(mesh, ls, cut_approximation="linear", cut_approximation_order=2)
    with pytest.raises(ValueError, match="cut_approximation"):
        cutcells.cut(mesh, ls, cut_approximation="quadratic")
    with pytest.raises(ValueError, match="1 to 4"):
        cutcells.cut(mesh, ls, cut_approximation="iso_p1", cut_approximation_order=5)
