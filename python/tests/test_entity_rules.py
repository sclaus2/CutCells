# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Pieces of a part lying on facets, edges and vertices of the mesh: the
entities near the zero sets in cut (their sharing cells, Basix local numbers,
the level sets vanishing on them and the cells' sides), one threshold for
vanishing and touching zero, rules made once per entity (entity_quadrature),
and how cells take them (quadrature(entities=...))."""

import itertools
import math

import numpy as np
import pytest

import cutcells

from test_part import wedge_mesh
from test_part_lut import square_mesh
from test_part_order import physical_points
from test_quadrays import box_mesh

basix = pytest.importorskip("basix")

BACKENDS = ("quadrays", "lut")
POLICIES = ("split", "negative", "positive", "owner")
BASIX_CELL = {10: "tetrahedron", 12: "hexahedron", 13: "prism", 14: "pyramid", 5: "triangle", 9: "quadrilateral"}
VTK_TO_BASIX = {10: [0, 1, 2, 3], 12: [0, 1, 3, 2, 4, 5, 7, 6], 13: [0, 1, 2, 3, 4, 5], 14: [0, 1, 3, 2, 4],
                5: [0, 1, 2], 9: [0, 1, 3, 2]}


def mesh_3d(kind, n):
    return box_mesh(kind, n) if kind in ("hex", "tet") else wedge_mesh(kind, n)


def vtk_type(kind):
    return {"hex": 12, "tet": 10, "prism": 13, "pyramid": 14, "quad": 9, "tri": 5}[kind]


def cell_nodes(mesh, kind, c):
    """The mesh nodes of cell c's vertices in Basix order (the test meshes are in VTK order)."""
    conn, offsets = np.asarray(mesh.connectivity), np.asarray(mesh.offsets)
    row = conn[offsets[c]:offsets[c + 1]]
    return [int(row[v]) for v in VTK_TO_BASIX[vtk_type(kind)]]


def axis_plane(offset, axis=0):
    """The analytic plane x_axis = offset."""
    normal = np.eye(3)[axis]

    def box_bounds(lo, hi):
        b = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        b[0], b[1] = lo[axis] - offset, hi[axis] - offset
        b[2 + 2 * axis] = b[3 + 2 * axis] = 1.0
        return b

    return cutcells.AnalyticLevelSet(lambda x: float(x[axis] - offset), lambda x: normal, box_bounds)


def plane_level_set(mesh, source, offset):
    if source == "P1":
        return cutcells.create_level_set(mesh, lambda x: x[0] - offset, degree=1, name="phi")
    return axis_plane(offset)


def total(rules):
    return float(np.sum(rules.weights))


def n_cells(visualization):
    """Cells of a QuadraysLeafMesh or a CutMesh."""
    return visualization.n_cells() if hasattr(visualization, "n_cells") else len(visualization.parent_map)


def tensor_tet_mesh(gx, gy, gz):
    """Kuhn tetrahedra of the tensor grid gx x gy x gz, VTK vertex order."""
    coords = np.array([[x, y, z] for z in gz for y in gy for x in gx], dtype=np.float64)
    nx, ny = len(gx), len(gy)

    def node(i, j, k):
        return i + nx * (j + ny * k)

    connectivity = []
    for k, j, i in itertools.product(range(len(gz) - 1), range(ny - 1), range(nx - 1)):
        for perm in itertools.permutations(range(3)):
            vertex = [i, j, k]
            path = [tuple(vertex)]
            for axis in perm:
                vertex[axis] += 1
                path.append(tuple(vertex))
            connectivity += [node(*v) for v in path]
    n = len(connectivity) // 4
    return cutcells.MeshView(coords, np.array(connectivity, dtype=np.int32),
                             np.arange(0, 4 * n + 1, 4, dtype=np.int32), np.full(n, 10, dtype=np.int32), tdim=3)


# ============================================================================
# Step 1: entities in cut
# ============================================================================


@pytest.mark.parametrize("kind", ["hex", "tet", "prism", "pyramid"])
@pytest.mark.parametrize("source", ["P1", "analytic"])
def test_plane_on_grid_faces_gives_each_facet_once(kind, source):
    """The plane x = 0.25 on grid faces: each facet in it once, with both cells
    sharing it and their sides (the cell below on the negative side)."""
    mesh = mesh_3d(kind, 8)
    result = cutcells.cut(mesh, plane_level_set(mesh, source, 0.25))
    facets = result.entities(2)
    zero = [e for e in range(facets.num_entities) if facets.zero[e] & 1]
    assert len(zero) == (128 if kind == "tet" else 64)
    x = np.asarray(mesh.coordinates).reshape(-1, 3)
    for e in zero:
        nodes = facets.nodes[facets.node_offsets[e]:facets.node_offsets[e + 1]]
        np.testing.assert_allclose(x[nodes, 0], 0.25, atol=1e-15)
        s0, s1 = facets.cell_offsets[e], facets.cell_offsets[e + 1]
        assert s1 - s0 == 2 and not facets.boundary[e]
        for s in range(s0, s1):
            below = x[cell_nodes(mesh, kind, facets.cells[s]), 0].mean() < 0.25
            assert (facets.negative[s], facets.positive[s]) == ((1, 0) if below else (0, 1))
    # every other facet holding the plane's nodes only is not in it
    assert int(np.count_nonzero(facets.vanishing)) == len(zero)
    # zero faces are derived, one per facet, owned by the cell below
    assert len(result.zero_faces) == len(zero)
    domains = result.domains[0]
    assert np.all(domains[result.zero_faces[:, 1]] == 0)


@pytest.mark.parametrize("kind", ["hex", "tet", "prism", "pyramid", "quad", "tri"])
def test_local_numbers_are_basix(kind):
    """The local number of every entity in every sharing cell names the
    entity's nodes in basix.topology; planes on grid faces in every direction
    (and the slanted faces of prisms, pyramids and Kuhn tetrahedra) reach
    every facet number."""
    two_d = kind in ("quad", "tri")
    mesh = square_mesh(kind, 4) if two_d else mesh_3d(kind, 4)
    tdim = 2 if two_d else 3
    planes = [lambda x, a=a: x[a] - 0.5 for a in range(tdim)]
    planes += {"prism": [lambda x: x[0] + x[1]], "pyramid": [lambda x: x[1] - x[2]],
               "tet": [lambda x: x[0] - x[1], lambda x: x[1] - x[2]]}.get(kind, [])
    level_sets = [cutcells.create_level_set(mesh, f, degree=1, name=f"p{i}") for i, f in enumerate(planes)]
    result = cutcells.cut(mesh, level_sets)
    topology = basix.topology(getattr(basix.CellType, BASIX_CELL[vtk_type(kind)]))
    seen = set()
    for dim in range(tdim):
        table = result.entities(dim)
        assert table.dim == dim and table.num_entities > 0
        for e in range(table.num_entities):
            nodes = set(table.nodes[table.node_offsets[e]:table.node_offsets[e + 1]].tolist())
            for s in range(table.cell_offsets[e], table.cell_offsets[e + 1]):
                cell = cell_nodes(mesh, kind, table.cells[s])
                assert {cell[v] for v in topology[dim][table.local_index[s]]} == nodes
                seen.add((dim, int(table.local_index[s])))
    assert {f for d, f in seen if d == tdim - 1} == set(range(len(topology[tdim - 1])))


def test_chart_of_quadrilaterals():
    """A quadrilateral's chart starts at its lowest node, then that node's
    lower and higher neighbours, then the opposite node: Basix's order."""
    mesh = box_mesh("hex", 4)
    result = cutcells.cut(mesh, cutcells.create_level_set(mesh, lambda x: x[0] - 0.5, degree=1, name="phi"))
    facets = result.entities(2)
    x = np.asarray(mesh.coordinates).reshape(-1, 3)
    for e in range(facets.num_entities):
        a, b, c, d = facets.nodes[facets.node_offsets[e]:facets.node_offsets[e + 1]]
        assert a == min(a, b, c, d) and b < c
        # (0,0), (1,0), (0,1), (1,1): the diagonals a-d and b-c cross at the centre
        np.testing.assert_allclose(x[a] + x[d], x[b] + x[c], atol=1e-14)


@pytest.mark.parametrize("kind", ["hex", "tet"])
@pytest.mark.parametrize("source", ["P1", "analytic"])
def test_planes_off_grid_faces(kind, source):
    """Planes offset by 0 to 1e-6 from grid faces: the interface measure is 4
    whatever the offset, in both backends and with every option, never
    doubled; the volume below is 5 + 4 d, continuous."""
    mesh = box_mesh(kind, 8)
    for d in [0.0, 1e-16, 1e-15, 3e-15, 1e-14, 2e-14, 5e-14, 1e-13, 1e-12, 1e-9, 1e-6, -1e-15, -1e-14, -1e-13]:
        result = cutcells.cut(mesh, plane_level_set(mesh, source, 0.25 + d))
        # one threshold: a facet is in the zero set, or a cell beside it is cut, not both
        assert result.num_cut_cells == 0 or len(result.zero_faces) == 0, d
        for backend in BACKENDS:
            for entities in POLICIES:
                rules = result["phi = 0"].quadrature(order=2, mode="cut_only", backend=backend, entities=entities)
                assert total(rules) == pytest.approx(4.0, rel=1e-13), (d, backend, entities)
            below = result["phi < 0"].quadrature(order=2, backend=backend)
            assert total(below) == pytest.approx(5.0 + 4.0 * d, rel=1e-13), (d, backend)


BOXES = (((-0.5, -0.5, -0.5), (0.0, 0.5, 0.5)), ((0.0, -0.5, -0.5), (0.5, 0.5, 0.5)))


def box_union(x):
    """min(box1, box2), each box the max of its six half-spaces."""
    return np.min([np.max([lo[0] - x[0], x[0] - hi[0], lo[1] - x[1], x[1] - hi[1], lo[2] - x[2], x[2] - hi[2]], axis=0)
                   for lo, hi in BOXES], axis=0)


def analytic_box_union():
    """box_union as an analytic level set, with interval bounds over boxes."""

    def gradient(x):
        e = np.eye(3)
        branches = []
        for lo, hi in BOXES:
            terms = [(lo[i] - x[i], -e[i]) for i in range(3)] + [(x[i] - hi[i], e[i]) for i in range(3)]
            branches.append(max(terms, key=lambda t: t[0]))
        return min(branches, key=lambda t: t[0])[1]

    def box_bounds(blo, bhi):
        lows, highs = [], []
        for lo, hi in BOXES:
            terms = [(lo[i] - bhi[i], lo[i] - blo[i]) for i in range(3)] + [(blo[i] - hi[i], bhi[i] - hi[i])
                                                                            for i in range(3)]
            lows.append(max(t[0] for t in terms))
            highs.append(max(t[1] for t in terms))
        return [min(lows), min(highs), -1.0, 1.0, -1.0, 1.0, -1.0, 1.0]

    return cutcells.AnalyticLevelSet(lambda x: float(box_union(np.asarray(x))), gradient, box_bounds)


@pytest.mark.parametrize("kind,source,n", [("hex", "P1", 8), ("hex", "analytic", 4), ("tet", "analytic", 4)])
def test_two_boxes_as_one_level_set(kind, source, n):
    """min(box1, box2) of two boxes side by side, faces on grid faces: it
    vanishes on the shared face, negative on both sides, so that face is no
    surface; phi = 0 is the union's surface, 6, in both backends. (A P1
    interpolant on tetrahedra vanishes on the tetrahedra whose vertices all
    lie on the boxes' faces, so tetrahedra take the level set itself.)"""
    mesh = box_mesh(kind, n)
    phi = (cutcells.create_level_set(mesh, box_union, degree=1, name="phi") if source == "P1"
           else analytic_box_union())
    result = cutcells.cut(mesh, phi)
    assert result.num_cut_cells == 0
    facets = result.entities(2)
    x = np.asarray(mesh.coordinates).reshape(-1, 3)
    shared = [e for e in range(facets.num_entities)
              if np.allclose(x[facets.nodes[facets.node_offsets[e]:facets.node_offsets[e + 1]], 0], 0.0)]
    per_face = (n // 2) ** 2 * (2 if kind == "tet" else 1)
    assert len(shared) == per_face
    for e in shared:
        assert facets.vanishing[e] == 1 and facets.zero[e] == 0
        assert np.all(facets.negative[facets.cell_offsets[e]:facets.cell_offsets[e + 1]] == 1)
    for backend in BACKENDS:
        assert total(result["phi = 0"].quadrature(order=2, backend=backend)) == pytest.approx(6.0, rel=1e-13)
        assert total(result["phi < 0"].quadrature(order=2, backend=backend)) == pytest.approx(1.0, rel=1e-13)
    assert result["phi = 0"].entity_quadrature(order=2).facets.num_rules == 6 * per_face


# ============================================================================
# Step 2: rules on entities, and how cells take them
# ============================================================================


@pytest.mark.parametrize("kind", ["hex", "tet", "prism", "pyramid"])
@pytest.mark.parametrize("source", ["P1", "analytic"])
def test_plane_area_once_with_every_option(kind, source):
    """The plane on grid faces gives its area once, in both backends and with
    every option; with exclude the cells keep none of it, and
    entity_quadrature has it all."""
    mesh = mesh_3d(kind, 8)
    result = cutcells.cut(mesh, plane_level_set(mesh, source, 0.25))
    part = result["phi = 0"]
    for backend in BACKENDS:
        for entities in POLICIES:
            assert total(part.quadrature(order=3, backend=backend, entities=entities)) == pytest.approx(4.0, rel=1e-13)
        assert total(part.quadrature(order=3, backend=backend, entities="exclude")) == pytest.approx(0.0, abs=1e-13)
        on = part.entity_quadrature(order=3, backend=backend)
        assert on.facets.num_rules == len(result.zero_faces)
        assert float(np.sum(on.facets.weights)) == pytest.approx(4.0, rel=1e-13)
        assert on.ridges.num_rules == 0 and on.vertices.num_rules == 0


@pytest.mark.parametrize("kind", ["hex", "tet"])
def test_split_gives_each_cell_half(kind):
    """With split each sharing cell takes half of each facet; with negative
    the cell below takes it all, as its zero faces did."""
    mesh = box_mesh(kind, 8)
    result = cutcells.cut(mesh, cutcells.create_level_set(mesh, lambda x: x[0] - 0.25, degree=1, name="phi"))
    part = result["phi = 0"]
    on = part.entity_quadrature(order=2).facets
    per_cell = {}
    for r in range(on.num_rules):
        area = float(np.sum(on.weights[on.offsets[r]:on.offsets[r + 1]]))
        for c in on.cells[on.cell_offsets[r]:on.cell_offsets[r + 1]]:
            per_cell[int(c)] = per_cell.get(int(c), 0.0) + area
    split = part.quadrature(order=2)
    assert len(split.parent_map) == len(per_cell) == 2 * on.num_rules
    for k, c in enumerate(split.parent_map):
        w = float(np.sum(split.weights[split.offset[k]:split.offset[k + 1]]))
        assert w == pytest.approx(0.5 * per_cell[int(c)], rel=1e-14)
    assert total(split) == pytest.approx(4.0, rel=1e-13)
    negative = part.quadrature(order=2, entities="negative")
    owners = sorted(result.zero_faces[:, 1].tolist())
    assert sorted(negative.parent_map.tolist()) == owners
    assert np.all(result.domains[0][negative.parent_map] == 0)
    assert total(negative) == pytest.approx(4.0, rel=1e-13)
    positive = part.quadrature(order=2, entities="positive")
    assert np.all(result.domains[0][positive.parent_map] == 2)


@pytest.mark.parametrize("kind", ["hex", "tet"])
def test_entity_points_in_cells(kind):
    """An entity rule's points: in the facet's reference cell, and in each
    sharing cell's reference coordinates at the same physical points."""
    mesh = box_mesh(kind, 4)
    result = cutcells.cut(mesh, cutcells.create_level_set(mesh, lambda x: x[0] - 0.5, degree=1, name="phi"))
    on = result["phi = 0"].entity_quadrature(order=4).facets
    assert on.entity_dim == 2 and on.tdim == 3
    assert on.points.shape == (len(on.weights), 2) and on.cell_points.shape[1] == 3
    assert np.all(on.weights > 0)
    x = np.asarray(mesh.coordinates).reshape(-1, 3)
    element = basix.create_element(basix.ElementFamily.P,
                                   getattr(basix.CellType, BASIX_CELL[vtk_type(kind)]), 1)
    for r in range(on.num_rules):
        physical = []
        for s in range(on.cell_offsets[r], on.cell_offsets[r + 1]):
            xi = on.cell_points[on.cell_point_offsets[s]:on.cell_point_offsets[s + 1]]
            assert len(xi) == on.offsets[r + 1] - on.offsets[r]
            vertices = x[cell_nodes(mesh, kind, on.cells[s])]
            physical.append(element.tabulate(0, xi)[0, :, :, 0] @ vertices)
        np.testing.assert_allclose(physical[0], physical[1], atol=1e-14)
        np.testing.assert_allclose(physical[0][:, 0], 0.5, atol=1e-14)


@pytest.mark.parametrize("backend", BACKENDS)
def test_line_along_mesh_edges(backend):
    """x = 0.1 and y = 0.1 on Kuhn tetrahedra whose grid has planes at 0.1:
    the line runs along mesh edges where both vanish, and gives length 2
    once, with every option; each edge rule lists all tetrahedra around it."""
    g = np.array([-1.0, -0.5, 0.0, 0.1, 0.5, 1.0])
    mesh = tensor_tet_mesh(g, g, np.linspace(-1.0, 1.0, 5))
    a = cutcells.create_level_set(mesh, lambda x: x[0] - 0.1, degree=1, name="a")
    b = cutcells.create_level_set(mesh, lambda x: x[1] - 0.1, degree=1, name="b")
    result = cutcells.cut(mesh, [a, b])
    part = result["a = 0 and b = 0"]
    for entities in POLICIES:
        assert total(part.quadrature(order=2, backend=backend, entities=entities)) == pytest.approx(2.0, rel=1e-14)
    on = part.entity_quadrature(order=2, backend=backend).ridges
    assert on.num_rules == 4 and float(np.sum(on.weights)) == pytest.approx(2.0, rel=1e-14)
    edges = result.entities(1)
    for r in range(on.num_rules):
        cells = on.cells[on.cell_offsets[r]:on.cell_offsets[r + 1]]
        x = np.asarray(mesh.coordinates).reshape(-1, 3)
        nodes = on.entity_nodes[on.node_offsets[r]:on.node_offsets[r + 1]]
        np.testing.assert_allclose(x[nodes, :2], 0.1, atol=1e-15)
        # every tetrahedron holding both nodes shares the edge
        conn = np.asarray(mesh.connectivity).reshape(-1, 4)
        holding = np.flatnonzero(np.isin(conn, nodes).sum(axis=1) == 2)
        np.testing.assert_array_equal(cells, holding)
        assert on.zero[r] == 3
    assert int(np.count_nonzero(edges.zero == 3)) == 4


@pytest.mark.parametrize("kind", ["quad", "tri"])
@pytest.mark.parametrize("backend", BACKENDS)
def test_points_at_vertices_and_on_edges(kind, backend):
    """Two lines crossing at a grid vertex, and on a grid edge: each point
    once, with weight 1 and all the cells sharing its entity, drawn once."""
    mesh = square_mesh(kind, 8)
    for where, point in (("vertex", (0.25, 0.5)), ("edge", (0.25, 0.4))):
        lines = [cutcells.create_level_set(mesh, lambda x, s=s: s * (x[0] - point[0]) - (x[1] - point[1]),
                                           degree=1, name=name) for s, name in ((0.37, "a"), (-3.45, "b"))]
        part = cutcells.cut(mesh, lines)["a = 0 and b = 0"]
        on = part.entity_quadrature(order=2, backend=backend)
        rules = on.ridges if where == "vertex" else on.facets
        other = on.facets if where == "vertex" else on.ridges
        assert rules.num_rules == 1 and other.num_rules == 0
        np.testing.assert_array_equal(rules.weights, [1.0])
        n = {"vertex": {"quad": 4, "tri": 6}[kind], "edge": 2}[where]
        assert len(rules.cells) == n
        for entities in POLICIES:
            assert total(part.quadrature(order=2, backend=backend, entities=entities)) == pytest.approx(1.0, rel=1e-14)
        split = part.quadrature(order=2, backend=backend)
        np.testing.assert_allclose(physical_points(mesh, split), np.tile(point, (n, 1)), atol=1e-12)
        assert n_cells(part.visualization_mesh(backend=backend)) == 1


def test_lut_count_facets_of_cut_cells_once():
    """A P2 level set vanishing on grid faces and cutting the cells beside
    them elsewhere: (x - 0.25)(x - 0.45). The lookup tables gave the facets
    twice (from the sub-cells below them and as zero faces); now once."""
    for kind in ("hex", "tet"):
        mesh = box_mesh(kind, 8)
        phi = cutcells.create_level_set(mesh, lambda x: (x[0] - 0.25) * (x[0] - 0.45), degree=2, name="phi")
        result = cutcells.cut(mesh, phi)
        assert len(result.zero_faces) == (128 if kind == "tet" else 64)
        for backend in BACKENDS:
            for entities in POLICIES:
                rules = result["phi = 0"].quadrature(order=2, mode="cut_only", backend=backend, entities=entities)
                assert total(rules) == pytest.approx(8.0, rel=1e-12), (kind, backend, entities)


@pytest.mark.parametrize("backend", BACKENDS)
def test_visualization_draws_each_facet_once(backend):
    mesh = box_mesh("hex", 4)
    result = cutcells.cut(mesh, cutcells.create_level_set(mesh, lambda x: x[0] - 0.5, degree=1, name="phi"))
    assert n_cells(result["phi = 0"].visualization_mesh(mode="cut_only", backend=backend)) == 16


def test_entities_of_a_part_and_bad_options():
    mesh = box_mesh("hex", 4)
    result = cutcells.cut(mesh, cutcells.create_level_set(mesh, lambda x: x[0] - 0.5, degree=1, name="phi"))
    part = result["phi = 0"]
    held = part.entities(2)
    assert len(held) == 16 and np.all(result.entities(2).zero[held] == 1)
    assert len(part.entities(1)) == 0 and len(result["phi < 0"].entities(2)) == 0
    with pytest.raises(ValueError, match="entities takes"):
        part.quadrature(order=2, entities="half")
    with pytest.raises(IndexError):
        result.entities(3)
    on = part.entity_quadrature(order=2)
    assert len(on) == 3 and on[2].num_rules == on.facets.num_rules == 16
    two_d = cutcells.cut(square_mesh("quad", 4),
                         cutcells.create_level_set(square_mesh("quad", 4), lambda x: x[0] - 0.5, degree=1, name="p"))
    assert two_d["p = 0"].entity_quadrature(order=2).vertices is None


def test_sides_of_a_cut_cell():
    """A cell that the level set cuts and that shares a facet it vanishes on
    takes its side from the derivative into it."""
    mesh = box_mesh("hex", 8)
    phi = cutcells.create_level_set(mesh, lambda x: (x[0] - 0.25) * (x[0] - 0.45), degree=2, name="phi")
    result = cutcells.cut(mesh, phi)
    facets = result.entities(2)
    x = np.asarray(mesh.coordinates).reshape(-1, 3)
    for e in np.flatnonzero(facets.zero):
        for s in range(facets.cell_offsets[e], facets.cell_offsets[e + 1]):
            c = facets.cells[s]
            above = x[cell_nodes(mesh, "hex", c), 0].mean() > 0.25
            # above: cut, negative next to the facet; below: positive
            assert result.domains[0][c] == (1 if above else 2)
            assert (facets.negative[s], facets.positive[s]) == ((1, 0) if above else (0, 1))
