# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""The lookup-table backend of cutcells.part (backend="lut"): rules on the
straight pieces of Pk-iso-P1 templates against today's straight backend
(AdaptCell) cell by cell, convergence in the template order, several level sets
in one cell, faces in zero sets, visualisation and the API."""

import math

import numpy as np
import pytest

import cutcells

from test_quadrays import CENTRE, RADIUS, box_mesh, sphere

BALL = 4.0 / 3.0 * math.pi * RADIUS**3
SPHERE = 4.0 * math.pi * RADIUS**2
CENTRE_2D = CENTRE[:2]


def square_mesh(kind, n):
    """[-1, 1]^2 with n^2 quadrilaterals or 2 n^2 triangles, VTK vertex order."""
    g = np.linspace(-1.0, 1.0, n + 1)
    coords = np.array([[x, y] for y in g for x in g], dtype=np.float64)

    def node(i, j):
        return i + (n + 1) * j

    connectivity, types = [], []
    for j in range(n):
        for i in range(n):
            if kind == "quad":
                connectivity += [node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)]
                types.append(9)
            else:
                connectivity += [node(i, j), node(i + 1, j), node(i + 1, j + 1)]
                connectivity += [node(i, j), node(i + 1, j + 1), node(i, j + 1)]
                types += [5, 5]
    width = 4 if kind == "quad" else 3
    offsets = np.arange(0, len(connectivity) + 1, width, dtype=np.int32)
    return cutcells.MeshView(coords, np.array(connectivity, dtype=np.int32), offsets,
                             np.array(types, dtype=np.int32), tdim=2)


def circle(x):
    return (x[0] - CENTRE_2D[0]) ** 2 + (x[1] - CENTRE_2D[1]) ** 2 - RADIUS**2


def mesh_of(kind, n):
    return box_mesh(kind, n) if kind in ("hex", "tet") else square_mesh(kind, n)


def round_level_set(kind):
    return sphere if kind in ("hex", "tet") else circle


def plane(x):
    return x[0] + 0.3 * x[1] - 0.1


def per_cell(rules, moment=None):
    """Weights summed per cell, or the moment sum w xi[moment] of the
    reference coordinates."""
    w = np.asarray(rules.weights)
    if moment is not None:
        w = w * np.asarray(rules.points).reshape(len(w), -1)[:, moment]
    offsets, parents = np.asarray(rules.offset), np.asarray(rules.parent_map)
    sums = {}
    for i, c in enumerate(parents):
        sums[int(c)] = sums.get(int(c), 0.0) + float(np.sum(w[offsets[i]:offsets[i + 1]]))
    return sums


def assert_same_cells(a, b, atol):
    cells = set(a) | set(b)
    worst = max(abs(a.get(c, 0.0) - b.get(c, 0.0)) for c in cells)
    assert worst < atol


def total(part, order=2, mode="full", **kwargs):
    return float(np.sum(part.quadrature(order=order, mode=mode, backend="lut", **kwargs).weights))


@pytest.mark.parametrize("kind", ["hex", "tet", "quad", "tri"])
@pytest.mark.parametrize("degree", [1, 2])
def test_matches_straight_backend(kind, degree):
    """One Pk level set: lut's template is today's iso-Pk refinement, so every
    cell gets the same measure and first moments as the straight backend."""
    mesh = mesh_of(kind, 6)
    phi = cutcells.create_level_set(mesh, round_level_set(kind), degree=degree, name="phi")
    new, old = cutcells.part.cut(mesh, phi), cutcells.cut(mesh, phi)
    tdim = 3 if kind in ("hex", "tet") else 2
    for expr, mode in [("phi < 0", "full"), ("phi > 0", "full"), ("phi = 0", "cut_only")]:
        lut = new[expr].quadrature(order=2, mode=mode, backend="lut")
        straight = old[expr].quadrature(order=3, mode=mode, backend="straight")
        assert_same_cells(per_cell(lut), per_cell(straight), 1e-14)
        for d in range(tdim):
            assert_same_cells(per_cell(lut, d), per_cell(straight, d), 1e-14)


@pytest.mark.parametrize("kind", ["hex", "tet"])
def test_template_order(kind):
    """The analytic sphere's values at the template's vertices: errors fall as
    (h / k)^2 with the template order k."""
    mesh = box_mesh(kind, 6)
    result = cutcells.part.cut(mesh, cutcells.analytic_sphere(CENTRE, RADIUS))
    volume, area = [], []
    for k in range(1, 5):
        options = cutcells.LutOptions(template_order=k)
        volume.append(abs(total(result["phi < 0"], options=options) / BALL - 1))
        area.append(abs(total(result["phi = 0"], mode="cut_only", options=options) / SPHERE - 1))
        assert total(result["phi < 0 or phi > 0"], options=options) == pytest.approx(8.0, rel=1e-13)
    assert math.log(volume[0] / volume[3]) / math.log(4) == pytest.approx(2.0, abs=0.1)
    assert math.log(area[0] / area[3]) / math.log(4) == pytest.approx(2.0, abs=0.1)


@pytest.mark.parametrize("kind", ["hex", "tet", "quad", "tri"])
def test_crossing_level_sets(kind):
    """A P2 sphere (circle) and a plane crossing in cells: the four parts fill
    the box, unions are counted once, and the sphere's surface splits into its
    parts on both sides of the plane. On simplices and quadrilaterals they
    match the straight backend cell by cell (which refuses this on hexahedra)."""
    mesh = mesh_of(kind, 6)
    a = cutcells.create_level_set(mesh, round_level_set(kind), degree=2, name="a")
    b = cutcells.create_level_set(mesh, plane, degree=1, name="b")
    result = cutcells.part.cut(mesh, [a, b])
    box = 8.0 if kind in ("hex", "tet") else 4.0
    quarters = [total(result[f"a {s} 0 and b {t} 0"]) for s in "<>" for t in "<>"]
    assert sum(quarters) == pytest.approx(box, rel=1e-13)
    assert total(result["a < 0 or b < 0"]) == pytest.approx(box - quarters[3], rel=1e-13)
    assert total(result["a < 0"]) == pytest.approx(quarters[0] + quarters[1], rel=1e-13)
    surface = total(result["a = 0"], mode="cut_only")
    split = total(result["a = 0 and b < 0"], mode="cut_only") + total(result["a = 0 and b > 0"], mode="cut_only")
    assert split == pytest.approx(surface, rel=1e-13)
    if kind == "hex":
        return
    old = cutcells.cut(mesh, [a, b])
    for expr, mode in [("a < 0 and b < 0", "full"), ("a < 0 or b < 0", "full"), ("a > 0 and b < 0", "full"),
                       ("a = 0 and b < 0", "cut_only"), ("b = 0 and a > 0", "cut_only")]:
        lut = result[expr].quadrature(order=2, mode=mode, backend="lut")
        straight = old[expr].quadrature(order=3, mode=mode, backend="straight")
        assert_same_cells(per_cell(lut), per_cell(straight), 1e-14)


@pytest.mark.parametrize("degree", [1, 2])
def test_hexahedra_with_multilinear_values(degree):
    """The Qk interpolant of the distance is multilinear on the template's
    hexahedra; they are cut as Kuhn tetrahedra, so both sides fill every cell."""

    def distance(x):
        return np.sqrt((x[0] - CENTRE[0]) ** 2 + (x[1] - CENTRE[1]) ** 2 + (x[2] - CENTRE[2]) ** 2) - RADIUS

    mesh = box_mesh("hex", 6)
    result = cutcells.part.cut(mesh, cutcells.create_level_set(mesh, distance, degree=degree, name="phi"))
    below = per_cell(result["phi < 0"].quadrature(order=1, mode="cut_only", backend="lut"))
    above = per_cell(result["phi > 0"].quadrature(order=1, mode="cut_only", backend="lut"))
    cell = (2.0 / 6) ** 3
    assert max(abs(below.get(c, 0.0) + above.get(c, 0.0) - cell) for c in result.cut_cells) < 1e-15


@pytest.mark.parametrize("kind", ["hex", "tet"])
def test_planes_in_faces(kind):
    """A plane through mesh faces (x = 0.25) is integrated on the owned zero
    faces; one through the faces of the template's sub-cells (x = 0.125, order
    2) is counted once by the lookup tables."""
    mesh = box_mesh(kind, 8)
    phi = cutcells.create_level_set(mesh, lambda x: x[0] - 0.25, degree=1, name="phi")
    result = cutcells.part.cut(mesh, phi)
    assert result.num_cut_cells == 0
    assert total(result["phi = 0"], mode="cut_only") == pytest.approx(4.0, rel=1e-13)
    assert total(result["phi < 0"]) == pytest.approx(5.0, rel=1e-13)
    assert total(result["phi > 0"]) == pytest.approx(3.0, rel=1e-13)

    phi = cutcells.create_level_set(mesh, lambda x: x[0] - 0.125, degree=1, name="phi")
    result = cutcells.part.cut(mesh, phi)
    options = cutcells.LutOptions(template_order=2)
    assert total(result["phi = 0"], mode="cut_only", options=options) == pytest.approx(4.0, rel=1e-13)
    assert total(result["phi < 0"], options=options) == pytest.approx(4.5, rel=1e-13)
    assert total(result["phi > 0"], options=options) == pytest.approx(3.5, rel=1e-13)


# Kuhn tetrahedra of the cells of a CutMesh (Basix vertex order), by cell type
TETRAHEDRA = {3: [(0, 1, 2, 3)], 5: [(0, 1, 3, 7), (0, 1, 5, 7), (0, 2, 3, 7), (0, 2, 6, 7), (0, 4, 5, 7),
                                     (0, 4, 6, 7)],
              6: [(0, 2, 1, 3), (1, 3, 5, 4), (1, 2, 5, 3)], 7: [(0, 1, 2, 4), (1, 3, 2, 4)]}


def cut_mesh_volume(mesh):
    x = np.asarray(mesh.vertex_coords)
    connectivity, offsets, types = np.asarray(mesh.connectivity), np.asarray(mesh.offset), np.asarray(mesh.types)
    volume = 0.0
    for c, t in enumerate(types):
        v = x[connectivity[offsets[c]:offsets[c + 1]]]
        for tet in TETRAHEDRA[int(t)]:
            a, b, d, e = v[list(tet)]
            volume += abs(np.linalg.det(np.array([b - a, d - a, e - a]))) / 6
    return volume


@pytest.mark.parametrize("kind", ["hex", "tet"])
def test_visualization(kind, tmp_path):
    mesh = box_mesh(kind, 4)
    result = cutcells.part.cut(mesh, cutcells.analytic_sphere(CENTRE, RADIUS))
    options = cutcells.LutOptions(template_order=2)
    volume = result["phi < 0"].visualization_mesh(mode="full", backend="lut", options=options)
    assert isinstance(volume, cutcells.CutMesh_float64)
    assert set(np.unique(volume.vtk_types)) <= {10, 12, 13, 14}
    assert cut_mesh_volume(volume) == pytest.approx(total(result["phi < 0"], options=options), rel=1e-12)
    # whole cells, and the cut cells whose template values reach below 0
    parents = set(np.asarray(volume.parent_map).tolist())
    part = result["phi < 0"]
    assert set(part.uncut_cells.tolist()) <= parents <= set(part.cut_cells.tolist()) | set(part.uncut_cells.tolist())
    interface = result["phi = 0"].visualization_mesh(mode="cut_only", backend="lut")
    assert set(np.unique(interface.vtk_types)) <= {5, 9}
    path = tmp_path / "ball.vtu"
    result["phi < 0"].write_vtu(str(path), mode="full", backend="lut")
    assert path.stat().st_size > 0


def test_api():
    mesh = box_mesh("hex", 2)
    result = cutcells.part.cut(mesh, cutcells.analytic_sphere(CENTRE, RADIUS))
    options = cutcells.LutOptions(template_order=3, triangulate=True)
    assert options.template_order == 3 and options.triangulate
    assert cutcells.LutOptions().template_order == 0
    with pytest.raises(TypeError, match="LutOptions"):
        result["phi < 0"].quadrature(backend="lut", options=cutcells.QuadraysOptions())
    with pytest.raises(TypeError, match="QuadraysOptions"):
        result["phi < 0"].quadrature(backend="quadrays", options=options)
    with pytest.raises(ValueError, match="template order"):
        result["phi < 0"].quadrature(backend="lut", options=cutcells.LutOptions(template_order=5))
    with pytest.raises(ValueError, match="unknown backend"):
        result["phi < 0"].visualization_mesh(backend="straight")
    planes = cutcells.part.cut(mesh, [cutcells.create_level_set(mesh, lambda x: x[0] - 0.1, degree=1, name="a"),
                                      cutcells.create_level_set(mesh, lambda x: x[1] - 0.1, degree=1, name="b")])
    with pytest.raises(ValueError, match="two level sets vanish"):
        planes["a = 0 and b = 0"].quadrature(backend="lut")
    # triangulated pieces give the same rules' measure
    assert total(result["phi < 0"], options=options) == pytest.approx(
        total(result["phi < 0"], options=cutcells.LutOptions(template_order=3)), rel=1e-13)
