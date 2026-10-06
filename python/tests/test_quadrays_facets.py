# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""quadrays on facets: meshes whose cells have a lower dimension than their
space, as the facets of a mesh do (triangles and quadrilaterals in 3D,
intervals in 2D and 3D), and intervals in 1D. Each facet integrates as the
cell it is: areas and lengths on the facet, with points in its reference
cell."""

import math

import numpy as np
import pytest

import cutcells

from test_part_lut import CENTRE_2D, square_mesh


def rotation():
    """A rotation about the axis (1, 2, 3) by 0.7."""
    a = np.array([1.0, 2.0, 3.0]) / math.sqrt(14.0)
    c, s = math.cos(0.7), math.sin(0.7)
    k = np.array([[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]])
    return c * np.eye(3) + s * k + (1 - c) * np.outer(a, a)


Q = rotation()
SHIFT = np.array([0.3, -0.2, 0.5])


def move(x):
    """Points (rows) rotated by Q and shifted."""
    return np.asarray(x, dtype=np.float64) @ Q.T + SHIFT


def cube_surface(kind, n, moved=True):
    """The boundary of [-1, 1]^3 with n x n quadrilaterals per face, or twice
    as many triangles, as a mesh in 3D (VTK vertex order)."""
    nodes, coords = {}, []

    def node(x):
        key = tuple(int(round(c * n)) for c in x)
        if key not in nodes:
            nodes[key] = len(coords)
            coords.append(x)
        return nodes[key]

    h = 2.0 / n
    connectivity, types = [], []
    for a in range(3):
        for side in (-1.0, 1.0):
            for j in range(n):
                for i in range(n):
                    def corner(di, dj):
                        x = [0.0, 0.0, 0.0]
                        x[a] = side
                        x[(a + 1) % 3] = -1 + h * (i + di)
                        x[(a + 2) % 3] = -1 + h * (j + dj)
                        return node(x)

                    v0, v1, v2, v3 = corner(0, 0), corner(1, 0), corner(0, 1), corner(1, 1)
                    if kind == "quad":
                        connectivity += [v0, v1, v3, v2]
                        types.append(9)
                    else:
                        connectivity += [v0, v1, v3, v0, v3, v2]
                        types += [5, 5]
    x = np.array(coords)
    if moved:
        x = move(x)
    width = 4 if kind == "quad" else 3
    offsets = np.arange(0, len(connectivity) + 1, width, dtype=np.int32)
    return cutcells.MeshView(np.ascontiguousarray(x), np.array(connectivity, dtype=np.int32), offsets,
                             np.array(types, dtype=np.int32), tdim=2)


def square_sides(n, gdim):
    """The boundary of [-1, 1]^2 with n intervals per side, in 2D or moved into 3D."""
    h = 2.0 / n
    points = []
    for side in range(4):
        for i in range(n):
            t = -1 + h * i
            points.append([(t, -1.0), (1.0, t), (-t, 1.0), (-1.0, -t)][side])
    x = np.array(points)
    if gdim == 3:
        x = move(np.column_stack([x, np.zeros(len(x))]))
    m = len(points)
    connectivity = np.array([[i, (i + 1) % m] for i in range(m)], dtype=np.int32).ravel()
    offsets = np.arange(0, 2 * m + 1, 2, dtype=np.int32)
    return cutcells.MeshView(np.ascontiguousarray(x), connectivity, offsets, np.full(m, 3, dtype=np.int32), tdim=1)


def interval_mesh(n):
    """[-1, 1] with n intervals, in 1D."""
    x = np.linspace(-1.0, 1.0, n + 1).reshape(-1, 1)
    connectivity = np.array([[i, i + 1] for i in range(n)], dtype=np.int32).ravel()
    offsets = np.arange(0, 2 * n + 1, 2, dtype=np.int32)
    return cutcells.MeshView(x, connectivity, offsets, np.full(n, 3, dtype=np.int32), tdim=1)


def facet_points(mesh, rules):
    """The rules' points in physical coordinates on affine cells: vertex 0 plus
    the reference coordinates along the edges to Basix vertices 1 (and 2)."""
    coords = np.asarray(mesh.coordinates)
    conn, offsets = np.asarray(mesh.connectivity), np.asarray(mesh.offsets)
    tdim = rules.tdim
    xi = np.asarray(rules.points).reshape(-1, tdim)
    offset = np.asarray(rules.offset)
    out = np.empty((len(xi), coords.shape[1]))
    for k, c in enumerate(np.asarray(rules.parent_map)):
        v = coords[conn[offsets[c]:offsets[c + 1]]]
        if len(v) == 4:  # a quadrilateral in VTK order: Basix vertex 2 is VTK vertex 3
            v = v[[0, 1, 3]]
        p = slice(offset[k], offset[k + 1])
        out[p] = v[0] + xi[p] @ (v[1:tdim + 1] - v[0])
    return out


def total(part, order=8, mode="full"):
    return float(np.sum(part.quadrature(order=order, mode=mode).weights))


# ============================================================================
# Triangles and quadrilaterals in 3D
# ============================================================================

BALL_CENTRE = np.array([0.1, -0.05, 0.08])
BALL_RADIUS = 1.2


def disks_and_circles():
    """The sphere cuts a disk from each face of the cube: their areas and circles."""
    areas, lengths = 0.0, 0.0
    for a in range(3):
        for side in (-1.0, 1.0):
            rho2 = BALL_RADIUS**2 - (side - BALL_CENTRE[a]) ** 2
            areas += math.pi * rho2
            lengths += 2 * math.pi * math.sqrt(rho2)
    return areas, lengths


@pytest.mark.parametrize("kind", ["quad", "tri"])
@pytest.mark.parametrize("ls_kind", ["P2", "analytic"])
@pytest.mark.parametrize("moved", [False, True])
def test_cube_surface(kind, ls_kind, moved):
    """The surface of the cube, in the coordinate planes or rotated: the disks
    a sphere cuts from its faces, their circles, and the faces' total area."""
    mesh = cube_surface(kind, 16, moved)
    assert mesh.gdim == 3 and mesh.tdim == 2
    centre = move(BALL_CENTRE) if moved else BALL_CENTRE
    if ls_kind == "P2":
        phi = cutcells.create_level_set(
            mesh, lambda x: sum((x[i] - centre[i]) ** 2 for i in range(3)) - BALL_RADIUS**2, degree=2, name="phi")
    else:
        phi = cutcells.analytic_sphere(centre, BALL_RADIUS)
    result = cutcells.cut(mesh, phi)
    areas, lengths = disks_and_circles()
    assert total(result["phi < 0"]) == pytest.approx(areas, rel=1e-9)
    assert total(result["phi > 0"]) == pytest.approx(24.0 - areas, rel=1e-9)
    assert total(result["phi = 0"], mode="cut_only") == pytest.approx(lengths, rel=1e-8)
    rules = result["phi = 0"].quadrature(order=4, mode="cut_only")
    assert rules.tdim == 2
    xi = np.asarray(rules.points).reshape(-1, 2)
    assert np.all(xi >= -1e-12) and np.all(xi <= 1 + 1e-12)
    if kind == "tri":
        assert np.all(xi.sum(axis=1) <= 1 + 1e-12)
    assert np.all(np.asarray(rules.weights) > 0)
    assert set(np.asarray(rules.parent_map)) <= set(np.asarray(result.cut_cells))
    x = facet_points(mesh, rules)
    np.testing.assert_allclose(np.linalg.norm(x - centre, axis=1), BALL_RADIUS, atol=1e-10)


@pytest.mark.parametrize("backend", ["quadrays", "lut"])
def test_plane_on_the_surface(backend):
    """A plane cuts the rotated cube's surface along a polygon: both backends
    integrate it exactly, from a P1 level set."""
    mesh = cube_surface("tri", 4)
    normal = Q @ np.array([0.0, 0.0, 1.0])
    phi = cutcells.create_level_set(mesh, lambda x: sum(normal[i] * (x[i] - SHIFT[i]) for i in range(3)) - 0.3,
                                    degree=1, name="phi")
    result = cutcells.cut(mesh, phi, backend=backend)
    # z' < 0.3: the bottom face, and the four sides up to z' = 0.3
    assert total(result["phi < 0"], order=2) == pytest.approx(4.0 + 4 * 2 * 1.3, rel=1e-13)
    assert total(result["phi = 0"], order=2, mode="cut_only") == pytest.approx(8.0, rel=1e-13)


@pytest.mark.parametrize("backend", ["quadrays", "lut"])
def test_plane_on_mesh_edges_of_the_surface(backend):
    """The plane z = 0.5 of the rotated cube on its surface mesh: its zero set,
    the square around the cube's sides, lies on 16 mesh edges, rules on the
    edges measured in 3D with the two triangles sharing each on either side."""
    mesh = cube_surface("tri", 4)
    normal = Q @ np.array([0.0, 0.0, 1.0])
    phi = cutcells.create_level_set(mesh, lambda x: sum(normal[i] * (x[i] - SHIFT[i]) for i in range(3)) - 0.5,
                                    degree=1, name="phi")
    result = cutcells.cut(mesh, phi, backend=backend)
    assert len(result.cut_cells) == 0
    part = result["phi = 0"]
    edges = part.entity_quadrature(order=2).facets
    assert len(edges.cell_offsets) - 1 == 16
    np.testing.assert_array_equal(np.diff(edges.cell_offsets), np.full(16, 2))
    assert np.sum(edges.weights) == pytest.approx(8.0, rel=1e-13)
    assert total(part, order=2) == pytest.approx(8.0, rel=1e-13)
    assert total(result["phi < 0"], order=2) == pytest.approx(16.0, rel=1e-13)


def test_off_plane_sphere_on_a_2d_mesh():
    """On a mesh in the plane z = 0, a sphere centred off the plane cuts the
    circle of radius sqrt(r^2 - z^2): its normal derivative does not count."""
    mesh = square_mesh("tri", 16)
    r, z = 0.7, 0.4
    result = cutcells.cut(mesh, cutcells.analytic_sphere([CENTRE_2D[0], CENTRE_2D[1], z], r))
    rho = math.sqrt(r * r - z * z)
    assert total(result["phi < 0"]) == pytest.approx(math.pi * rho**2, rel=1e-9)
    assert total(result["phi = 0"], mode="cut_only") == pytest.approx(2 * math.pi * rho, rel=1e-8)


def test_crossing_points_on_facets():
    """Where two level sets vanish on facets in 3D: points of weight 1."""
    mesh = cube_surface("tri", 8)
    e = Q.T  # rows: the cube's axes after rotation
    a = cutcells.create_level_set(mesh, lambda x: sum(e[0][i] * (x[i] - SHIFT[i]) for i in range(3)) - 0.37,
                                  degree=1, name="a")
    b = cutcells.create_level_set(mesh, lambda x: sum(e[1][i] * (x[i] - SHIFT[i]) for i in range(3)) - 0.61,
                                  degree=1, name="b")
    result = cutcells.cut(mesh, [a, b])
    part = result["a = 0 and b = 0"]
    assert part.dim == 0
    rules = part.quadrature(order=2)
    np.testing.assert_array_equal(np.asarray(rules.weights), [1.0, 1.0])
    exact = move(np.array([[0.37, 0.61, -1.0], [0.37, 0.61, 1.0]]))
    found = facet_points(mesh, rules)
    found = found[np.argsort(found @ (Q @ np.array([0.0, 0.0, 1.0])))]
    np.testing.assert_allclose(found, exact, atol=1e-12)


def test_cell_rules_in_3d():
    """The per-cell bindings take a facet's vertices in 3D."""
    tri = np.array([[0.1, 0.2], [0.9, 0.35], [0.3, 0.8]])
    tri3 = move(np.column_stack([tri, np.zeros(3)]))
    values = tri @ np.array([1.0, 0.4]) - 0.42  # a line, P1: the vertex values
    for selection in ("phi < 0", "phi = 0"):
        own, _ = cutcells.quadrays_cell_rules(cutcells.CellType.triangle, tri, 1, values, selection)
        moved, _ = cutcells.quadrays_cell_rules(cutcells.CellType.triangle, tri3, 1, values, selection)
        assert moved.tdim == 2
        assert np.sum(moved.weights) == pytest.approx(np.sum(own.weights), rel=1e-14)
    # a parallelogram keeps its reference box: lifted to z = 0.7, the same rule
    # for an analytic sphere; rotated, Taylor models bound the sphere a little
    # differently, and bisect elsewhere, but the measures agree
    quad = np.array([[0.1, 0.1], [0.8, 0.2], [0.3, 0.7], [1.0, 0.8]])
    lifted = np.column_stack([quad, np.full(4, 0.7)])
    rotated = move(np.column_stack([quad, np.zeros(4)]))
    centre = np.array([0.5, 0.45, 0.0])
    for selection in ("phi < 0", "phi = 0"):
        own, _ = cutcells.quadrays_cell_rules(cutcells.CellType.quadrilateral, quad,
                                              cutcells.analytic_sphere(centre, 0.3), selection)
        up, _ = cutcells.quadrays_cell_rules(cutcells.CellType.quadrilateral, lifted,
                                             cutcells.analytic_sphere(centre + [0.0, 0.0, 0.7], 0.3), selection)
        moved, _ = cutcells.quadrays_cell_rules(cutcells.CellType.quadrilateral, rotated,
                                                cutcells.analytic_sphere(move(centre), 0.3), selection)
        np.testing.assert_allclose(np.asarray(up.points), np.asarray(own.points), atol=1e-14)
        np.testing.assert_allclose(np.asarray(up.weights), np.asarray(own.weights), rtol=1e-14)
        assert np.sum(moved.weights) == pytest.approx(np.sum(own.weights), rel=1e-12)


# ============================================================================
# Intervals
# ============================================================================

CIRCLE_CENTRE = np.array([0.13, -0.07])
CIRCLE_RADIUS = 1.2


def chords():
    """The circle's chords on the four sides of the square, all inside them."""
    return sum(2 * math.sqrt(CIRCLE_RADIUS**2 - d * d)
               for d in (1 + CIRCLE_CENTRE[0], 1 - CIRCLE_CENTRE[0], 1 + CIRCLE_CENTRE[1], 1 - CIRCLE_CENTRE[1]))


@pytest.mark.parametrize("gdim", [2, 3])
@pytest.mark.parametrize("ls_kind", ["P2", "analytic"])
def test_square_sides(gdim, ls_kind):
    """The square's sides as intervals in 2D, and rotated into 3D, where a
    sphere centred off their plane cuts the same circle: its chords, and the 8
    points where it crosses the sides, each with weight 1."""
    mesh = square_sides(8, gdim)
    if gdim == 2:
        centre, radius = np.array([*CIRCLE_CENTRE, 0.0]), CIRCLE_RADIUS
    else:
        centre, radius = move(np.array([*CIRCLE_CENTRE, 0.5])), math.sqrt(CIRCLE_RADIUS**2 + 0.25)
    if ls_kind == "P2":
        phi = cutcells.create_level_set(
            mesh, lambda x: sum((x[i] - centre[i]) ** 2 for i in range(gdim)) - radius**2, degree=2, name="phi")
    else:
        phi = cutcells.analytic_sphere(centre, radius)
    result = cutcells.cut(mesh, phi)
    assert total(result["phi < 0"], order=4) == pytest.approx(chords(), rel=1e-12)
    assert total(result["phi > 0"], order=4) == pytest.approx(8.0 - chords(), rel=1e-12)
    rules = result["phi = 0"].quadrature(order=4)
    assert rules.tdim == 1
    np.testing.assert_array_equal(np.asarray(rules.weights), np.ones(8))
    x = facet_points(mesh, rules)
    np.testing.assert_allclose(np.linalg.norm(x - centre[:gdim], axis=1), radius, atol=1e-12)


def test_intervals_in_1d():
    """[-1, 1]: (x - 0.13)^2 < 0.6^2 on (-0.47, 0.73), and its two ends."""
    mesh = interval_mesh(8)
    phi = cutcells.create_level_set(mesh, lambda x: (x[0] - 0.13) ** 2 - 0.36, degree=2, name="phi")
    result = cutcells.cut(mesh, phi)
    assert total(result["phi < 0"], order=4) == pytest.approx(1.2, rel=1e-13)
    assert total(result["phi > 0"], order=4) == pytest.approx(0.8, rel=1e-13)
    rules = result["phi = 0"].quadrature(order=4)
    np.testing.assert_array_equal(np.asarray(rules.weights), [1.0, 1.0])
    np.testing.assert_allclose(np.sort(facet_points(mesh, rules)[:, 0]), [-0.47, 0.73], atol=1e-13)


def test_zero_at_vertices():
    """A level set vanishing at two mesh vertices, changing sign across them:
    no cell is cut there, and each point is the rule of its vertex, the facet
    of the two intervals sharing it, one on either side. By default each
    interval takes half of it; with entities="owner" the lower one takes it."""
    mesh = interval_mesh(8)
    phi = cutcells.create_level_set(mesh, lambda x: (x[0] - 0.25) * (x[0] + 0.5), degree=2, name="phi")
    result = cutcells.cut(mesh, phi)
    assert len(result.cut_cells) == 0
    assert total(result["phi < 0"], order=4) == pytest.approx(0.75, rel=1e-13)
    part = result["phi = 0"]
    on = part.entity_quadrature(order=4).facets
    np.testing.assert_array_equal(on.weights, [1.0, 1.0])
    np.testing.assert_array_equal(np.diff(on.cell_offsets), [2, 2])
    assert np.all((np.asarray(on.negative).reshape(2, 2) != 0).sum(axis=1) == 1)
    np.testing.assert_array_equal(np.asarray(part.quadrature(order=4).weights), np.full(4, 0.5))
    owner = part.quadrature(order=4, entities="owner")
    np.testing.assert_array_equal(np.asarray(owner.weights), [1.0, 1.0])
    np.testing.assert_allclose(np.sort(facet_points(mesh, owner)[:, 0]), [-0.5, 0.25], atol=1e-15)
    assert part.visualization_mesh(mode="cut_only").n_cells() == 2


def test_two_zero_sets_on_intervals():
    mesh = interval_mesh(4)
    a = cutcells.create_level_set(mesh, lambda x: x[0] - 0.3, degree=1, name="a")
    b = cutcells.create_level_set(mesh, lambda x: x[0] - 0.31, degree=1, name="b")
    with pytest.raises(ValueError, match="two level sets vanish"):
        cutcells.cut(mesh, [a, b])["a = 0 and b = 0"].quadrature(order=2)


# ============================================================================
# Leaves
# ============================================================================

def test_facet_leaves(tmp_path):
    """Leaves on facets: Lagrange quadrilaterals along the facets' normals,
    curves on the sphere; on intervals, curves and vertices."""
    mesh = cube_surface("tri", 8)
    centre = move(BALL_CENTRE)
    result = cutcells.cut(mesh, cutcells.analytic_sphere(centre, BALL_RADIUS))
    disks = result["phi < 0"].visualization_mesh(mode="full", degree=2)
    assert set(np.unique(disks.vtk_types)) == {5, 70}
    # each Lagrange quadrilateral faces as its facet does: the cross product of
    # its diagonals (leaves may collapse, or have no width where a segment
    # does not, which leaves keep)
    points = np.asarray(disks.points).reshape(-1, 3)
    coords, conn, offsets = np.asarray(mesh.coordinates), np.asarray(mesh.connectivity), np.asarray(mesh.offsets)
    cells_conn, cells_offsets = np.asarray(disks.connectivity), np.asarray(disks.offsets)
    parents = np.asarray(disks.parent)
    facing = []
    for k in np.flatnonzero(np.asarray(disks.vtk_types) == 70):
        x = points[cells_conn[cells_offsets[k]:cells_offsets[k] + 4]]
        v = coords[conn[offsets[parents[k]]:offsets[parents[k] + 1]]]
        facet_normal = np.cross(v[1] - v[0], v[2] - v[0])
        leaf_normal = np.cross(x[2] - x[0], x[3] - x[1])
        facing.append(np.dot(leaf_normal, facet_normal) / np.dot(facet_normal, facet_normal))
    assert min(facing) > -1e-12 and max(facing) > 0
    circles = result["phi = 0"].visualization_mesh(mode="cut_only", degree=2)
    assert set(np.unique(circles.vtk_types)) == {68}
    np.testing.assert_allclose(np.linalg.norm(np.asarray(circles.points).reshape(-1, 3) - centre, axis=1),
                               BALL_RADIUS, atol=1e-9)
    path = tmp_path / "disks.vtu"
    result["phi < 0"].write_vtu(str(path), mode="full", degree=2)
    assert path.stat().st_size > 0

    sides = square_sides(8, 3)
    result = cutcells.cut(sides, cutcells.analytic_sphere(move([*CIRCLE_CENTRE, 0.0]), CIRCLE_RADIUS))
    chords = result["phi < 0"].visualization_mesh(mode="full", degree=2)
    assert set(np.unique(chords.vtk_types)) == {3, 68}
    ends = result["phi = 0"].visualization_mesh(mode="cut_only")
    assert set(np.unique(ends.vtk_types)) == {1} and ends.n_cells() == 8
