# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""The front end without AdaptCell (cutcells.part) with the quadrays backend:
parts selected by expressions against exact values, for Pk and analytic level
sets, interfaces lying in mesh faces, and the API."""

import math

import numpy as np
import pytest

import cutcells

from test_quadrays import CENTRE, RADIUS, box_mesh, sphere

BALL = 4.0 / 3.0 * math.pi * RADIUS**3
SPHERE = 4.0 * math.pi * RADIUS**2


def total(part, order=5, mode="full"):
    return float(np.sum(part.quadrature(order=order, mode=mode).weights))


def level_set(kind, mesh):
    if kind == "P2":
        return cutcells.create_level_set(mesh, sphere, degree=2, name="phi")
    return cutcells.analytic_sphere(CENTRE, RADIUS)


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
@pytest.mark.parametrize("ls_kind", ["P2", "analytic"])
def test_sphere_parts(mesh_kind, ls_kind):
    mesh = box_mesh(mesh_kind, 8)
    result = cutcells.part.cut(mesh, level_set(ls_kind, mesh))
    assert total(result["phi < 0"]) == pytest.approx(BALL, rel=1e-7)
    assert total(result["phi <= 0"]) == pytest.approx(BALL, rel=1e-7)
    assert total(result["phi > 0"]) == pytest.approx(8.0 - BALL, rel=1e-7)
    assert total(result["phi < 0 or phi > 0"]) == pytest.approx(8.0, rel=1e-12)
    assert total(result["phi = 0"], mode="cut_only") == pytest.approx(SPHERE, rel=1e-6)
    rules = result["phi < 0"].quadrature(order=3)
    points = np.asarray(rules.points).reshape(-1, 3)
    assert np.all(np.asarray(rules.weights) > 0)
    assert np.all(points >= -1e-12) and np.all(points <= 1 + 1e-12)
    # one rule per cell, cells ascending
    assert np.all(np.diff(np.asarray(rules.parent_map)) > 0)


def test_cut_cells_hold_the_interface():
    """Bounds over each tetrahedron and its halves find exactly the cells the
    sphere crosses: each holds a piece of the interface."""
    mesh = box_mesh("tet", 8)
    result = cutcells.part.cut(mesh, level_set("P2", mesh))
    rules = result["phi = 0"].quadrature(order=3, mode="cut_only")
    assert result.num_cut_cells == len(rules.parent_map) == 627
    np.testing.assert_array_equal(result.cut_cells, rules.parent_map)


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
@pytest.mark.parametrize("ls_kind", ["P1", "analytic"])
def test_plane_in_mesh_faces(mesh_kind, ls_kind):
    """The plane x = 0.25 lies in mesh faces: no cell is cut, and every face in
    it is owned once, by the cell below, so the interface is integrated once."""
    mesh = box_mesh(mesh_kind, 8)
    if ls_kind == "P1":
        phi = cutcells.create_level_set(mesh, lambda x: x[0] - 0.25, degree=1, name="phi")
    else:
        normal, offset = np.array([1.0, 0.0, 0.0]), 0.25

        def box_bounds(lo, hi):
            return [lo[0] - offset, hi[0] - offset, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]

        phi = cutcells.AnalyticLevelSet(lambda x: float(x[0] - offset), lambda x: normal, box_bounds)
    result = cutcells.part.cut(mesh, phi)
    assert result.num_cut_cells == 0
    zero_faces = result.zero_faces
    assert len(zero_faces) == (64 if mesh_kind == "hex" else 128)
    domains = result.domains[0]
    assert np.all(domains[zero_faces[:, 1]] == 0)  # owners lie below the plane
    assert total(result["phi = 0"], order=2, mode="cut_only") == pytest.approx(4.0, rel=1e-12)
    assert total(result["phi < 0"], order=2) == pytest.approx(5.0, rel=1e-12)
    assert total(result["phi > 0"], order=2) == pytest.approx(3.0, rel=1e-12)


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
def test_ball_inside_one_cell(mesh_kind):
    """A ball smaller than a cell reaching 0.01 into its neighbours, away from
    all vertices: classification and integration both use the level set itself,
    no interpolant is involved."""
    centre, radius = [0.25, 0.25, 0.25], 0.26
    mesh = box_mesh(mesh_kind, 4)
    result = cutcells.part.cut(mesh, cutcells.analytic_sphere(centre, radius))
    # exactly the cells the ball crosses (dense sampling gives the same 18 tetrahedra)
    assert result.num_cut_cells == (7 if mesh_kind == "hex" else 18)
    assert total(result["phi < 0"]) == pytest.approx(4.0 / 3.0 * math.pi * radius**3, rel=1e-8)
    assert total(result["phi = 0"], mode="cut_only") == pytest.approx(4.0 * math.pi * radius**2, rel=1e-8)


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
def test_two_level_sets(mesh_kind):
    """A ball and a plane that never meet in a cell; terms on both."""
    mesh = box_mesh(mesh_kind, 8)
    ball = cutcells.analytic_sphere([-0.45, 0.02, -0.03], 0.5)

    def plane_bounds(lo, hi):
        return [lo[0] - 0.3, hi[0] - 0.3, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]

    plane = cutcells.AnalyticLevelSet(lambda x: float(x[0] - 0.3), lambda x: [1.0, 0.0, 0.0], plane_bounds)
    result = cutcells.part.cut(mesh, [ball, plane])
    assert result.level_set_names == ["phi1", "phi2"]
    volume, area = 4.0 / 3.0 * math.pi * 0.5**3, 4.0 * math.pi * 0.5**2
    assert total(result["phi1 < 0 and phi2 < 0"]) == pytest.approx(volume, rel=1e-7)
    assert total(result["phi1 > 0 and phi2 < 0"]) == pytest.approx(4 * 1.3 - volume, rel=1e-7)
    assert total(result["phi1 = 0 and phi2 < 0"], mode="cut_only") == pytest.approx(area, rel=1e-6)
    assert total(result["phi2 = 0 and phi1 > 0"], mode="cut_only") == pytest.approx(4.0, rel=1e-12)
    assert total(result["phi1 < 0 or phi2 > 0"]) == pytest.approx(volume + 4 * 0.7, rel=1e-7)


def test_two_level_sets_in_one_cell_are_refused():
    mesh = box_mesh("hex", 4)
    first = cutcells.create_level_set(mesh, lambda x: x[0] - 0.1, degree=1, name="a")
    second = cutcells.create_level_set(mesh, lambda x: x[1] - 0.1, degree=1, name="b")
    result = cutcells.part.cut(mesh, [first, second])
    with pytest.raises(RuntimeError, match="phase 6"):
        result["a < 0 and b < 0"].quadrature(order=3)
    # each level set alone is fine
    assert total(result["a < 0"], order=2) == pytest.approx(2 * 2 * 1.1, rel=1e-12)


def test_visualization(tmp_path):
    mesh = box_mesh("hex", 4)
    result = cutcells.part.cut(mesh, level_set("analytic", mesh))
    volume = result["phi < 0"].visualization_mesh(mode="full", degree=2)
    interface = result["phi = 0"].visualization_mesh(mode="cut_only", degree=2)
    assert set(np.unique(volume.vtk_types)) <= {12, 72}
    assert set(np.unique(interface.vtk_types)) == {70}
    assert interface.n_cells() > 0 and volume.n_cells() > interface.n_cells()
    path = tmp_path / "ball.vtu"
    result["phi < 0"].write_vtu(str(path), mode="full", degree=2)
    assert path.stat().st_size > 0
    # zero faces show as linear faces
    plane = cutcells.create_level_set(mesh, lambda x: x[0] - 0.5, degree=1, name="phi")
    faces = cutcells.part.cut(mesh, plane)["phi = 0"].visualization_mesh(mode="cut_only")
    assert set(np.unique(faces.vtk_types)) == {9} and faces.n_cells() == 16


def test_api_checks():
    mesh = box_mesh("hex", 2)
    phi = cutcells.analytic_sphere(CENTRE, RADIUS)
    result = cutcells.part.cut(mesh, phi, names=["ball"])
    assert result.level_set_names == ["ball"]
    assert result["ball < 0"].dim == 3 and result["ball = 0"].dim == 2
    with pytest.raises(ValueError, match="unknown backend"):
        result["ball < 0"].quadrature(order=3, backend="bogus")
    with pytest.raises(ValueError, match="named"):
        cutcells.part.cut(mesh, [phi, phi], names=["a", "a"])
    with pytest.raises(TypeError):
        cutcells.part.cut(mesh, "phi")
