# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""The front end without AdaptCell (cutcells.part) with the quadrays backend:
parts selected by expressions against exact values, for Pk and analytic level
sets, interfaces lying in mesh faces, triangles and quadrilaterals, prisms and
pyramids, several level sets meeting in cells, and the API."""

import itertools
import math

import numpy as np
import pytest

import cutcells

from test_part_lut import CENTRE_2D, circle, square_mesh
from test_quadrays import CENTRE, RADIUS, box_mesh, plane, sphere

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
    assert total(result["phi = 0"], mode="cut_only") == pytest.approx(4.0 * math.pi * radius**2, rel=1e-7)


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


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
def test_two_level_sets_in_one_cell(mesh_kind):
    """Two P1 planes crossing in cells: the corner, the union and the corner's
    faces are exact; the line where both vanish comes from the lookup tables."""
    mesh = box_mesh(mesh_kind, 4)
    first = cutcells.create_level_set(mesh, lambda x: x[0] - 0.1, degree=1, name="a")
    second = cutcells.create_level_set(mesh, lambda x: x[1] - 0.1, degree=1, name="b")
    result = cutcells.part.cut(mesh, [first, second])
    assert total(result["a < 0 and b < 0"], order=2) == pytest.approx(2 * 1.1 * 1.1, rel=1e-12)
    assert total(result["a < 0 or b < 0"], order=2) == pytest.approx(8 - 2 * 0.9 * 0.9, rel=1e-12)
    assert total(result["a < 0 and b > 0"], order=2) == pytest.approx(2 * 1.1 * 0.9, rel=1e-12)
    faces = result["a = 0 and b < 0 or b = 0 and a < 0"]
    assert total(result["a = 0 and b < 0"], order=2, mode="cut_only") == pytest.approx(2 * 1.1, rel=1e-12)
    assert total(faces, order=2, mode="cut_only") == pytest.approx(4 * 1.1, rel=1e-12)
    with pytest.raises(ValueError, match="lookup tables"):
        result["a = 0 and b = 0"].quadrature(order=2)
    line = result["a = 0 and b = 0"].quadrature(order=2, mode="cut_only", backend="lut")
    assert float(np.sum(line.weights)) == pytest.approx(2.0, rel=1e-12)


def wedge_mesh(kind, n):
    """[-1, 1]^3 with 2 n^3 prisms (the cubes' bottom triangles extruded) or
    6 n^3 pyramids (the cubes' faces as bases, their centres as apices), VTK
    vertex order."""
    g = np.linspace(-1.0, 1.0, n + 1)
    coords = [[x, y, z] for z in g for y in g for x in g]

    def node(i, j, k):
        return i + (n + 1) * (j + (n + 1) * k)

    connectivity, types = [], []
    for k, j, i in itertools.product(range(n), repeat=3):
        if kind == "prism":
            for tri in ([(i, j), (i + 1, j), (i, j + 1)], [(i + 1, j), (i + 1, j + 1), (i, j + 1)]):
                connectivity += [node(a, b, k) for a, b in tri] + [node(a, b, k + 1) for a, b in tri]
                types.append(13)
            continue
        apex = len(coords)
        coords.append([(g[i] + g[i + 1]) / 2, (g[j] + g[j + 1]) / 2, (g[k] + g[k + 1]) / 2])
        for axis, side in itertools.product(range(3), range(2)):
            # the face of the cube normal to axis, its vertices around it
            corners = []
            for a, b in [(0, 0), (1, 0), (1, 1), (0, 1)]:
                v = [i, j, k]
                v[axis] += side
                v[(axis + 1) % 3] += a
                v[(axis + 2) % 3] += b
                corners.append(node(*v))
            connectivity += corners + [apex]
            types.append(14)
    width = 6 if kind == "prism" else 5
    offsets = np.arange(0, len(connectivity) + 1, width, dtype=np.int32)
    return cutcells.MeshView(np.array(coords, dtype=np.float64), np.array(connectivity, dtype=np.int32), offsets,
                             np.array(types, dtype=np.int32), tdim=3)


@pytest.mark.parametrize("mesh_kind", ["quad", "tri"])
@pytest.mark.parametrize("ls_kind", ["P2", "analytic"])
def test_disk_parts(mesh_kind, ls_kind):
    """Triangles and quadrilaterals: the disk, its complement and the circle."""
    mesh = square_mesh(mesh_kind, 16)
    if ls_kind == "P2":
        phi = cutcells.create_level_set(mesh, circle, degree=2, name="phi")
    else:
        phi = cutcells.analytic_sphere(CENTRE_2D, RADIUS)
    result = cutcells.part.cut(mesh, phi)
    disk, length = math.pi * RADIUS**2, 2 * math.pi * RADIUS
    assert total(result["phi < 0"]) == pytest.approx(disk, rel=1e-9)
    assert total(result["phi > 0"]) == pytest.approx(4.0 - disk, rel=1e-9)
    assert total(result["phi = 0"], mode="cut_only") == pytest.approx(length, rel=1e-8)
    rules = result["phi = 0"].quadrature(order=3, mode="cut_only")
    points = np.asarray(rules.points).reshape(-1, 2)
    assert np.all(np.asarray(rules.weights) > 0)
    assert np.all(points >= -1e-12) and np.all(points <= 1 + 1e-12)
    np.testing.assert_array_equal(result.cut_cells, rules.parent_map)


@pytest.mark.parametrize("mesh_kind", ["prism", "pyramid"])
@pytest.mark.parametrize("ls_kind", ["P2", "analytic"])
def test_sphere_on_prisms_and_pyramids(mesh_kind, ls_kind):
    """Prisms and pyramids: a P2 level set (on pyramids a rational function, as
    Basix's pyramid elements span) or an analytic one."""
    mesh = wedge_mesh(mesh_kind, 6)
    result = cutcells.part.cut(mesh, level_set(ls_kind, mesh))
    assert total(result["phi < 0"]) == pytest.approx(BALL, rel=1e-7)
    assert total(result["phi > 0"]) == pytest.approx(8.0 - BALL, rel=1e-7)
    assert total(result["phi = 0"], mode="cut_only") == pytest.approx(SPHERE, rel=1e-6)


@pytest.mark.parametrize("mesh_kind", ["prism", "pyramid"])
@pytest.mark.parametrize("degree", [1, 2])
def test_plane_on_prisms_and_pyramids(mesh_kind, degree):
    """A plane in general position as a P1 or P2 level set is exact."""
    mesh = wedge_mesh(mesh_kind, 3)
    result = cutcells.part.cut(mesh, cutcells.create_level_set(mesh, plane, degree=degree, name="phi"))
    # x + 0.3 y - 0.2 z < 0 halves [-1, 1]^3; the cut is 4 sqrt(1.13)
    assert total(result["phi < 0"], order=2) == pytest.approx(4.0, rel=1e-13)
    assert total(result["phi = 0"], order=2, mode="cut_only") == pytest.approx(4.0 * math.sqrt(1.13), rel=1e-13)


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
@pytest.mark.parametrize("ls_kind", ["P2", "analytic"])
def test_ball_and_half_space(mesh_kind, ls_kind):
    """A ball and a plane meeting in cells: the cap, its sphere, the disk the
    plane cuts from the ball, and the union. The plane is P1; the ball P2 or
    analytic, so one cut mixes the two kinds."""
    radius, height = 0.8, 0.3  # the plane z = CENTRE[2] + height
    mesh = box_mesh(mesh_kind, 8)
    if ls_kind == "P2":
        ball = cutcells.create_level_set(
            mesh, lambda x: sum((x[i] - CENTRE[i]) ** 2 for i in range(3)) - radius**2, degree=2, name="ball")
    else:
        ball = cutcells.analytic_sphere(CENTRE, radius, signed_distance=False)
    plane = cutcells.create_level_set(mesh, lambda x: x[2] - CENTRE[2] - height, degree=1, name="plane")
    result = cutcells.part.cut(mesh, [ball, plane], names=["ball", "plane"])
    h = radius - height
    cap = math.pi * h**2 * (3 * radius - h) / 3
    above = 4.0 * (1.0 - CENTRE[2] - height)
    assert total(result["ball < 0 and plane > 0"]) == pytest.approx(cap, rel=1e-7)
    assert total(result["ball = 0 and plane > 0"], mode="cut_only") == pytest.approx(2 * math.pi * radius * h,
                                                                                     rel=1e-6)
    assert total(result["plane = 0 and ball < 0"], mode="cut_only") == pytest.approx(
        math.pi * (radius**2 - height**2), rel=1e-7)
    ball_volume = 4.0 / 3.0 * math.pi * radius**3
    assert total(result["ball < 0 or plane > 0"]) == pytest.approx(ball_volume + above - cap, rel=1e-7)


@pytest.mark.parametrize("mesh_kind", ["quad", "tri"])
def test_visualization_2d(mesh_kind):
    mesh = square_mesh(mesh_kind, 8)
    result = cutcells.part.cut(mesh, cutcells.analytic_sphere(CENTRE_2D, RADIUS))
    disk = result["phi < 0"].visualization_mesh(mode="full", degree=2)
    curve = result["phi = 0"].visualization_mesh(mode="cut_only", degree=2)
    linear = 9 if mesh_kind == "quad" else 5
    assert set(np.unique(disk.vtk_types)) <= {linear, 70}
    assert 70 in set(np.unique(disk.vtk_types))
    assert set(np.unique(curve.vtk_types)) == {68}


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
