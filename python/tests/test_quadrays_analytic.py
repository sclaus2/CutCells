# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Analytic level sets in the quadrays backend: compiled (analytic_sphere),
Python callables, capsules and ShapeForest tapes (skipped without ShapeForest)."""

import math

import numpy as np
import pytest

import cutcells

from test_quadrays import CENTRE, RADIUS, box_mesh

BALL = 4.0 / 3.0 * math.pi * RADIUS**3
SPHERE = 4.0 * math.pi * RADIUS**2


def unit_tet():
    return np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])


def plane_callables(a, b):
    """phi = a . x - b as Python callables; its bounds are exact."""
    a = np.asarray(a, dtype=float)

    def value(x):
        return float(a @ x - b)

    def gradient(x):
        return a

    def box_bounds(lo, hi):
        low = float(np.sum(np.minimum(a * lo, a * hi)) - b)
        high = float(np.sum(np.maximum(a * lo, a * hi)) - b)
        return [low, high, a[0], a[0], a[1], a[1], a[2], a[2]]

    return cutcells.AnalyticLevelSet(value, gradient, box_bounds)


def distance_callables(centre, radius):
    """phi = |x - c| - r as Python callables, with interval bounds over boxes."""
    c = np.asarray(centre, dtype=float)

    def value(x):
        return float(np.linalg.norm(x - c) - radius)

    def gradient(x):
        return (x - c) / np.linalg.norm(x - c)

    def box_bounds(lo, hi):
        d_lo, d_hi = lo - c, hi - c
        nearest = np.where(d_lo > 0, d_lo, np.where(d_hi < 0, d_hi, 0.0))
        farthest = np.maximum(np.abs(d_lo), np.abs(d_hi))
        dmin, dmax = np.linalg.norm(nearest), np.linalg.norm(farthest)
        b = [dmin - radius, dmax - radius]
        if dmin <= 0.0:
            return b  # the gradient has no bound at the centre
        for i in range(3):
            # (x_i - c_i) / d with d in [dmin, dmax]
            q = [d_lo[i] / dmin, d_lo[i] / dmax, d_hi[i] / dmin, d_hi[i] / dmax]
            b += [min(q), max(q)]
        return b

    return cutcells.AnalyticLevelSet(value, gradient, box_bounds)


def totals(result, order=7):
    volume = result["phi < 0"].quadrature(order=order, mode="full", backend="quadrays")
    area = result["phi = 0"].quadrature(order=order, mode="cut_only", backend="quadrays")
    return volume, area


@pytest.mark.parametrize("kind", ["hex", "tet"])
@pytest.mark.parametrize("signed_distance", [True, False])
def test_analytic_sphere_totals(kind, signed_distance):
    mesh = box_mesh(kind, 8)
    phi = cutcells.analytic_sphere(CENTRE, RADIUS, signed_distance=signed_distance)
    result = cutcells.cut(mesh, phi)
    volume, area = totals(result)
    # as for the polynomial sphere (test_quadrays.py)
    assert np.sum(volume.weights) == pytest.approx(BALL, rel=1e-7)
    assert np.sum(area.weights) == pytest.approx(SPHERE, rel=1e-6)
    assert np.all(np.asarray(volume.weights) > 0) and np.all(np.asarray(area.weights) > 0)
    points = np.asarray(area.points).reshape(-1, 3)
    assert np.all(points >= -1e-12) and np.all(points <= 1 + 1e-12)


def test_analytic_beats_its_interpolant():
    """quadrays integrates the analytic sphere, the lookup tables its P2 interpolant."""
    mesh = box_mesh("hex", 6)
    result = cutcells.cut(mesh, cutcells.analytic_sphere(CENTRE, RADIUS))
    exact_error = abs(np.sum(totals(result)[0].weights) - BALL)
    straight = result["phi < 0"].quadrature(order=4, mode="full", backend="lut",
                                            options=cutcells.LutOptions(template_order=2))
    assert exact_error < 1e-8
    assert abs(np.sum(straight.weights) - BALL) > 100 * exact_error


@pytest.mark.parametrize("kind", ["hex", "tet"])
def test_cells_are_classified_by_the_analytic_level_set(kind):
    """A ball about the centre of the cell [0, 0.5]^3 that reaches 0.01 beyond
    its faces: all vertices of all cells lie outside, so a P1 interpolant sees
    no ball at all, and the caps lie away from every vertex. The cells are
    classified by the level set's own bounds, so quadrays finds them; the
    lookup tables see the level set at the template's vertices, P1 ones
    nothing. Dense sampling finds the same 18 tetrahedra crossed."""
    centre, radius = [0.25, 0.25, 0.25], 0.26
    mesh = box_mesh(kind, 4)
    result = cutcells.cut(mesh, cutcells.analytic_sphere(centre, radius))
    assert result.num_cut_cells == (7 if kind == "hex" else 18)
    volume = result["phi < 0"].quadrature(order=5, mode="full", backend="quadrays")
    area = result["phi = 0"].quadrature(order=5, mode="cut_only", backend="quadrays")
    assert np.sum(volume.weights) == pytest.approx(4.0 / 3.0 * math.pi * radius**3, rel=1e-6)
    assert np.sum(area.weights) == pytest.approx(4.0 * math.pi * radius**2, rel=1e-5)
    straight = result["phi < 0"].quadrature(order=2, mode="full", backend="lut",
                                            options=cutcells.LutOptions(template_order=1))
    assert np.sum(straight.weights) == 0.0


@pytest.mark.parametrize("kind,n", [("hex", 3), ("tet", 2)])
def test_plane_from_callables_is_exact(kind, n):
    mesh = box_mesh(kind, n)
    phi = plane_callables([1.0, 0.3, -0.2], 0.0)
    result = cutcells.cut(mesh, phi)
    volume, area = totals(result, order=3)
    assert np.sum(volume.weights) == pytest.approx(4.0, rel=1e-13)
    assert np.sum(area.weights) == pytest.approx(4.0 * math.sqrt(1.13), rel=1e-13)


def test_callables_match_the_compiled_sphere():
    """The callables bound the sphere over boxes, the compiled sphere by Taylor
    models over sub-boxes, which certify larger boxes: the rules differ, their
    values agree to the accuracy of q."""
    vertices = np.array([0.5, 0.25, 0.0, 0.75, 0.25, 0.0, 0.5, 0.5, 0.0, 0.5, 0.25, 0.25])
    compiled = cutcells.analytic_sphere(CENTRE, RADIUS)
    python = distance_callables(CENTRE, RADIUS)
    assert not python.has_taylor_bounds and compiled.has_taylor_bounds
    for selection in ["phi < 0", "phi = 0"]:
        for q, rel in [(5, 1e-6), (12, 1e-12)]:
            a, _ = cutcells.quadrays_cell_rules(cutcells.CellType.tetrahedron, vertices, compiled, selection, q=q)
            b, _ = cutcells.quadrays_cell_rules(cutcells.CellType.tetrahedron, vertices, python, selection, q=q)
            assert np.sum(a.weights) > 0
            assert np.sum(b.weights) == pytest.approx(np.sum(a.weights), rel=rel), selection


def test_capsule_round_trip():
    phi = cutcells.analytic_sphere(CENTRE, RADIUS)
    copy = cutcells.AnalyticLevelSet.from_capsule(phi.capsule)
    del phi  # the capsule keeps the level set alive
    original = cutcells.analytic_sphere(CENTRE, RADIUS)
    for selection in ["phi < 0", "phi = 0"]:
        a, _ = cutcells.quadrays_cell_rules(cutcells.CellType.tetrahedron, unit_tet(), original, selection)
        b, _ = cutcells.quadrays_cell_rules(cutcells.CellType.tetrahedron, unit_tet(), copy, selection)
        np.testing.assert_array_equal(a.weights, b.weights)
        np.testing.assert_array_equal(a.points, b.points)


def test_interface_queries():
    phi = cutcells.analytic_sphere([0.0, 0.0, 0.0], 0.5)
    assert phi.value([0.5, 0.0, 0.0]) == pytest.approx(0.0, abs=1e-15)
    np.testing.assert_allclose(phi.gradient([0.0, 2.0, 0.0]), [0.0, 1.0, 0.0])
    b = phi.box_bounds([0.3, 0.3, 0.3], [0.4, 0.4, 0.4])
    assert b[0] <= math.sqrt(0.27) - 0.5 and b[1] >= math.sqrt(0.48) - 0.5
    # around the centre the value has bounds, the gradient none
    b = phi.box_bounds([-0.1, -0.1, -0.1], [0.1, 0.1, 0.1])
    assert len(b) == 2 and b[0] <= -0.5 and math.sqrt(0.03) - 0.5 <= b[1] < 0
    models = phi.taylor_bounds(np.array([0.35, 0.35, 0.35]), np.ascontiguousarray(0.05 * np.eye(3)))
    assert models.shape == (4, 5)
    assert phi.taylor_bounds(np.zeros(3), np.ascontiguousarray(0.05 * np.eye(3))).shape == (1, 5)


def test_level_set_function_keeps_the_analytic_level_set():
    mesh = box_mesh("hex", 2)
    phi = cutcells.analytic_sphere(CENTRE, RADIUS)
    ls = cutcells.create_level_set(mesh, phi, 2)
    assert ls.has_dof_values() and ls.has_value()
    x = np.array([0.3, 0.2, 0.1])
    assert ls.value(x) == pytest.approx(phi.value(x))


def test_callable_errors_reach_python():
    def value(x):
        raise ValueError("no value here")

    phi = cutcells.AnalyticLevelSet(value, lambda x: [1.0, 0.0, 0.0], lambda lo, hi: [0.0] * 8)
    with pytest.raises(RuntimeError, match="no value here"):
        phi.value([0.0, 0.0, 0.0])
    short = cutcells.AnalyticLevelSet(lambda x: 0.0, lambda x: [1.0], lambda lo, hi: [0.0] * 8)
    with pytest.raises(RuntimeError, match="3 numbers"):
        short.gradient([0.0, 0.0, 0.0])


def test_capsule_checks():
    phi = cutcells.analytic_sphere(CENTRE, RADIUS)
    import ctypes

    PyCapsule_New = ctypes.pythonapi.PyCapsule_New
    PyCapsule_New.restype = ctypes.py_object
    PyCapsule_New.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p]
    name = ctypes.create_string_buffer(b"something.else")  # outlives the capsule
    wrong = PyCapsule_New(ctypes.c_void_p(1), ctypes.cast(name, ctypes.c_char_p), None)
    with pytest.raises(ValueError, match="cutcells.AnalyticLevelSet"):
        cutcells.AnalyticLevelSet.from_capsule(wrong)
    del wrong
    assert phi.capsule is not None


def test_tape_arrays_refuse_unsupported_opcodes():
    arrays = dict(
        op=np.array([1, 16], dtype=np.uint8),  # INPUT_X, TAN
        a=np.array([-1, 0], dtype=np.int32),
        b=np.array([-1, -1], dtype=np.int32),
        c=np.array([-1, -1], dtype=np.int32),
        out=np.array([0, 1], dtype=np.int32),
        imm=np.zeros(2),
        n_registers=2,
        output=1,
        inputs=np.array([0, -1, -1], dtype=np.int32),
        extra_registers=np.zeros(0, dtype=np.int32),
        extra_values=np.zeros(0),
    )
    with pytest.raises(ValueError, match="unsupported opcode 16"):
        cutcells.analytic_level_set_from_tape(**arrays)


def test_leaves_of_an_analytic_level_set():
    phi = cutcells.analytic_sphere([0.2, 0.2, 0.2], 0.5)
    leaves = cutcells.quadrays_cell_leaves(cutcells.CellType.hexahedron,
                                           np.array([0.0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0,
                                                     0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1]),
                                           phi, "phi = 0", leaf_degree=2)
    assert leaves.n_cells() > 0
    points = np.asarray(leaves.points)
    distance = np.linalg.norm(points - 0.2, axis=1)
    assert np.max(np.abs(distance - 0.5)) < 1e-4  # nodes are pulled 1e-5 inside their segments


# ============================================================================
# ShapeForest (optional)
# ============================================================================


def test_shapeforest_sphere_matches_the_compiled_one():
    sf = pytest.importorskip("shapeforest")
    from cutcells import shapeforest as csf

    shape = sf.sphere(RADIUS).translate(*CENTRE)
    tape_phi = csf.analytic_level_set(shape)
    compiled = cutcells.analytic_sphere(CENTRE, RADIUS)
    mesh = box_mesh("tet", 6)
    a = totals(cutcells.cut(mesh, tape_phi))
    b = totals(cutcells.cut(mesh, compiled))
    for x, y in zip(a, b):
        assert np.sum(x.weights) == pytest.approx(np.sum(y.weights), rel=1e-12)
    assert np.sum(a[0].weights) == pytest.approx(BALL, rel=1e-8)


def test_shapeforest_union_with_kinks():
    sf = pytest.importorskip("shapeforest")
    from cutcells import shapeforest as csf

    # two overlapping balls: min of two distances, a kink where they meet
    shape = sf.union(sf.sphere(0.45).translate(-0.2, 0.0, 0.0), sf.sphere(0.45).translate(0.2, 0.0, 0.0))
    phi = csf.analytic_level_set(shape)
    mesh = box_mesh("hex", 8)
    volume, _ = totals(cutcells.cut(mesh, phi))
    # union volume: two balls minus their lens (distance d = 0.4)
    r, d = 0.45, 0.4
    lens = math.pi * (4 * r + d) * (2 * r - d) ** 2 / 12
    assert np.sum(volume.weights) == pytest.approx(2 * 4 / 3 * math.pi * r**3 - lens, rel=1e-6)


def test_shapeforest_write_tape(tmp_path):
    sf = pytest.importorskip("shapeforest")
    from cutcells import shapeforest as csf

    path = tmp_path / "sphere.tape"
    csf.write_tape(sf.sphere(RADIUS).translate(*CENTRE), path)
    lines = path.read_text().splitlines()
    assert lines[0] == "shapeforest-tape 1"
    n_regs, out, x, y, z, n = (int(v) for v in lines[1].split())
    assert len(lines) == 2 + n and 0 <= out < n_regs
