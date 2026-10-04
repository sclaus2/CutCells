# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Conditions (cutcells.conditions): selections built from the level sets
themselves, next to the strings they stand for. Every selection of
test_part.py, on the same cuts (one test here per test there, of the same
name), selects the same part as a condition and as a string; conditions write
themselves as strings that parse back to them; level sets compare with 0
only, and among themselves and in hashes by identity; conditions refuse
bool(); strings take == as well as =."""

import functools
import itertools
import operator

import numpy as np
import pytest

import cutcells

from test_part import (CURVE, TRICYLINDER_FACES, ball_level_set, curve_shapes, cylinder_level_set, disk_level_set,
                       level_set, wedge_mesh)
from test_part_lut import CENTRE_2D, circle, square_mesh
from test_quadrays import CENTRE, RADIUS, box_mesh, plane

# the relations of selection strings as operators; <= and >= select as < and >
RELATIONS = {"<": operator.lt, "<=": operator.lt, ">": operator.gt, ">=": operator.gt, "=": operator.eq,
             "==": operator.eq}


def condition(text, level_sets):
    """The condition that a selection string names, built with the operators
    from the level sets, given by name."""

    def clause(words):
        name, relation, zero = words.split()
        assert zero == "0", words
        return RELATIONS[relation](level_sets[name], 0)

    terms = (functools.reduce(operator.and_, map(clause, term.split(" and "))) for term in text.split(" or "))
    return functools.reduce(operator.or_, terms)


def assert_same_part(part, other, what):
    assert part.dim == other.dim, what
    for cells in ("uncut_cells", "cut_cells", "zero_faces"):
        np.testing.assert_array_equal(getattr(part, cells), getattr(other, cells), err_msg=f"{what}: {cells}")


def check(result, level_sets, selections):
    """Each selection string selects the same part, not empty, as its condition."""
    for text in selections:
        part = result[text]
        assert part.num_uncut_cells + part.num_cut_cells + len(part.zero_faces) > 0, text
        assert_same_part(result[condition(text, level_sets)], part, text)


def axis_plane(axis, offset):
    """The analytic plane x_axis = offset."""
    normal = np.eye(3)[axis]

    def box_bounds(lo, hi):
        bounds = [lo[axis] - offset, hi[axis] - offset, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        bounds[2 + 2 * axis] = bounds[3 + 2 * axis] = 1.0
        return bounds

    return cutcells.AnalyticLevelSet(lambda x: float(x[axis] - offset), lambda x: normal, box_bounds)


# ============================================================================
# The selections of test_part.py, as conditions and as strings
# ============================================================================


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
@pytest.mark.parametrize("ls_kind", ["P2", "analytic"])
def test_sphere_parts(mesh_kind, ls_kind):
    """Also the cut of test_cut_cells_hold_the_interface."""
    mesh = box_mesh(mesh_kind, 8)
    phi = level_set(ls_kind, mesh)
    check(cutcells.cut(mesh, phi), {"phi": phi},
          ["phi < 0", "phi <= 0", "phi > 0", "phi < 0 or phi > 0", "phi = 0"])


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
@pytest.mark.parametrize("ls_kind", ["P1", "analytic"])
def test_plane_in_mesh_faces(mesh_kind, ls_kind):
    mesh = box_mesh(mesh_kind, 8)
    if ls_kind == "P1":
        phi = cutcells.create_level_set(mesh, lambda x: x[0] - 0.25, degree=1, name="phi")
    else:
        phi = axis_plane(0, 0.25)
    check(cutcells.cut(mesh, phi), {"phi": phi}, ["phi = 0", "phi < 0", "phi > 0"])


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
def test_ball_inside_one_cell(mesh_kind):
    phi = cutcells.analytic_sphere([0.25, 0.25, 0.25], 0.26)
    check(cutcells.cut(box_mesh(mesh_kind, 4), phi), {"phi": phi}, ["phi < 0", "phi = 0"])


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
def test_two_level_sets(mesh_kind):
    ball, wall = cutcells.analytic_sphere([-0.45, 0.02, -0.03], 0.5), axis_plane(0, 0.3)
    check(cutcells.cut(box_mesh(mesh_kind, 8), [ball, wall]), {"phi1": ball, "phi2": wall},
          ["phi1 < 0 and phi2 < 0", "phi1 > 0 and phi2 < 0", "phi1 = 0 and phi2 < 0", "phi2 = 0 and phi1 > 0",
           "phi1 < 0 or phi2 > 0"])


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
def test_two_level_sets_in_one_cell(mesh_kind):
    mesh = box_mesh(mesh_kind, 4)
    a, b, off, c = (cutcells.create_level_set(mesh, lambda x, i=i, s=s: x[i] - s, degree=1, name=name)
                    for i, s, name in ((0, 0.1, "a"), (1, 0.1, "b"), (1, 0.13, "b"), (2, 0.1, "c")))
    check(cutcells.cut(mesh, [a, b]), {"a": a, "b": b},
          ["a < 0 and b < 0", "a < 0 or b < 0", "a < 0 and b > 0", "a = 0 and b < 0 or b = 0 and a < 0",
           "a = 0 and b < 0"])
    check(cutcells.cut(mesh, [a, off]), {"a": a, "b": off}, [CURVE])
    check(cutcells.cut(mesh, [a, b, c]), {"a": a, "b": b, "c": c}, ["a = 0 and b = 0 and c = 0"])


@pytest.mark.parametrize("mesh_kind", ["quad", "tri"])
@pytest.mark.parametrize("ls_kind", ["P2", "analytic"])
def test_disk_parts(mesh_kind, ls_kind):
    mesh = square_mesh(mesh_kind, 16)
    if ls_kind == "P2":
        phi = cutcells.create_level_set(mesh, circle, degree=2, name="phi")
    else:
        phi = cutcells.analytic_sphere(CENTRE_2D, RADIUS)
    check(cutcells.cut(mesh, phi), {"phi": phi}, ["phi < 0", "phi > 0", "phi = 0"])


@pytest.mark.parametrize("mesh_kind", ["prism", "pyramid"])
@pytest.mark.parametrize("ls_kind", ["P2", "analytic"])
def test_sphere_on_prisms_and_pyramids(mesh_kind, ls_kind):
    mesh = wedge_mesh(mesh_kind, 6)
    phi = level_set(ls_kind, mesh)
    check(cutcells.cut(mesh, phi), {"phi": phi}, ["phi < 0", "phi > 0", "phi = 0"])


@pytest.mark.parametrize("mesh_kind", ["prism", "pyramid"])
@pytest.mark.parametrize("degree", [1, 2])
def test_plane_on_prisms_and_pyramids(mesh_kind, degree):
    mesh = wedge_mesh(mesh_kind, 3)
    phi = cutcells.create_level_set(mesh, plane, degree=degree, name="phi")
    check(cutcells.cut(mesh, phi), {"phi": phi}, ["phi < 0", "phi = 0"])


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
@pytest.mark.parametrize("ls_kind", ["P2", "analytic"])
def test_ball_and_half_space(mesh_kind, ls_kind):
    mesh = box_mesh(mesh_kind, 8)
    ball = ball_level_set(mesh, ls_kind, CENTRE, 0.8, "ball")
    half_space = cutcells.create_level_set(mesh, lambda x: x[2] - CENTRE[2] - 0.3, degree=1, name="plane")
    check(cutcells.cut(mesh, [ball, half_space], names=["ball", "plane"]), {"ball": ball, "plane": half_space},
          ["ball < 0 and plane > 0", "ball = 0 and plane > 0", "plane = 0 and ball < 0", "ball < 0 or plane > 0"])


@pytest.mark.parametrize("mesh_kind", ["quad", "tri"])
@pytest.mark.parametrize("ls_kind", ["P2", "analytic"])
def test_two_disks(mesh_kind, ls_kind):
    """Also the cut of test_crossing_points."""
    mesh = square_mesh(mesh_kind, 16)
    a = disk_level_set(mesh, ls_kind, CENTRE_2D + np.array([-0.2, 0.0]), 0.62, "a")
    b = disk_level_set(mesh, ls_kind, CENTRE_2D + np.array([0.3, 0.1]), 0.5, "b")
    check(cutcells.cut(mesh, [a, b], names=["a", "b"]), {"a": a, "b": b},
          ["a < 0 and b < 0", "a < 0 and b > 0", "a < 0 or b < 0", "a = 0 and b < 0", "b = 0 and a < 0",
           "a = 0 and b > 0", CURVE])


@pytest.mark.parametrize("mesh_kind", ["quad", "tri"])
def test_three_disks(mesh_kind):
    mesh = square_mesh(mesh_kind, 16)
    centres = [CENTRE_2D + np.array(offset) for offset in ([-0.25, -0.1], [0.25, -0.12], [0.02, 0.3])]
    a, b, c = (disk_level_set(mesh, "P2", centre, r, name) for centre, r, name in zip(centres, [0.5, 0.45, 0.48], "abc"))
    signs = list(itertools.product("<>", repeat=3))
    selections = [f"a {sa} 0 and b {sb} 0 and c {sc} 0" for sa, sb, sc in signs]
    selections += [f"a = 0 and b {sb} 0 and c {sc} 0" for _, sb, sc in signs[:4]]
    check(cutcells.cut(mesh, [a, b, c], names=["a", "b", "c"]), {"a": a, "b": b, "c": c}, selections)


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
def test_napkin_ring(mesh_kind):
    mesh = box_mesh(mesh_kind, 8)
    ball, cyl = ball_level_set(mesh, "P2", CENTRE, 0.8, "ball"), cylinder_level_set(mesh, 0.45, 2, "cyl")
    check(cutcells.cut(mesh, [ball, cyl], names=["ball", "cyl"]), {"ball": ball, "cyl": cyl},
          ["ball < 0 and cyl > 0", "ball = 0 and cyl > 0", "cyl = 0 and ball < 0"])


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
def test_tricylinder(mesh_kind):
    """Also the cuts of test_face_leaves_stay_in_their_part."""
    mesh = box_mesh(mesh_kind, 8)
    a, b, c = (cylinder_level_set(mesh, 0.6, axis, name) for axis, name in zip((2, 1, 0), "abc"))
    check(cutcells.cut(mesh, [a, b, c], names=["a", "b", "c"]), {"a": a, "b": b, "c": c},
          ["a < 0 and b < 0 and c < 0"] + TRICYLINDER_FACES)


@pytest.mark.parametrize("n", [4, 8])
@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
@pytest.mark.parametrize("ls_kind", ["P2", "analytic"])
def test_curves(n, mesh_kind, ls_kind):
    """Also the cuts of test_lut_curves_converge (n = 4) and the lens of
    test_curve_leaves."""
    mesh = box_mesh(mesh_kind, n)
    for _, (a, b), _ in curve_shapes(mesh, ls_kind):
        check(cutcells.cut(mesh, [a, b], names=["a", "b"]), {"a": a, "b": b}, [CURVE])


@pytest.mark.parametrize("mesh_kind", ["quad", "tri"])
@pytest.mark.parametrize("where", ["vertex", "edge"])
def test_crossing_point_on_cells_shared(mesh_kind, where):
    mesh = square_mesh(mesh_kind, 16)
    point = np.array([0.125, -0.25]) if where == "vertex" else np.array([0.125, -0.21])
    a, b = (cutcells.create_level_set(mesh, lambda x, s=s: s * (x[0] - point[0]) - (x[1] - point[1]), degree=1,
                                      name=name) for s, name in ((0.37, "a"), (0.8, "b")))
    check(cutcells.cut(mesh, [a, b]), {"a": a, "b": b}, [CURVE])


def test_curve_leaves():
    square = square_mesh("quad", 8)
    a = disk_level_set(square, "analytic", CENTRE_2D + np.array([-0.2, 0.0]), 0.62, "a")
    b = disk_level_set(square, "analytic", CENTRE_2D + np.array([0.3, 0.1]), 0.5, "b")
    check(cutcells.cut(square, [a, b], names=["a", "b"]), {"a": a, "b": b}, [CURVE])


@pytest.mark.parametrize("mesh_kind", ["quad", "tri"])
def test_visualization_2d(mesh_kind):
    phi = cutcells.analytic_sphere(CENTRE_2D, RADIUS)
    check(cutcells.cut(square_mesh(mesh_kind, 8), phi), {"phi": phi}, ["phi < 0", "phi = 0"])


def test_visualization():
    mesh = box_mesh("hex", 4)
    phi = level_set("analytic", mesh)
    check(cutcells.cut(mesh, phi), {"phi": phi}, ["phi < 0", "phi = 0"])
    flat = cutcells.create_level_set(mesh, lambda x: x[0] - 0.5, degree=1, name="phi")
    check(cutcells.cut(mesh, flat), {"phi": flat}, ["phi = 0"])


def test_api_checks():
    phi = cutcells.analytic_sphere(CENTRE, RADIUS)
    check(cutcells.cut(box_mesh("hex", 2), phi, names=["ball"]), {"ball": phi}, ["ball < 0", "ball = 0"])


# ============================================================================
# Conditions
# ============================================================================


def planes(mesh, names="abc"):
    """P1 planes x_i = 0.1, named."""
    return [cutcells.create_level_set(mesh, lambda x, i=i: x[i] - 0.1, degree=1, name=name)
            for i, name in enumerate(names)]


def test_conditions_read_as_their_strings():
    """& binds more tightly than | and distributes over it; 0 > a is a < 0."""
    mesh = box_mesh("hex", 4)
    a, b, c = planes(mesh)
    result = cutcells.cut(mesh, [a, b, c])
    for selection, text in [
            ((a > 0) & (b < 0), "a > 0 and b < 0"),
            ((a == 0) & (b < 0) | (b == 0) & (a < 0), "a = 0 and b < 0 or b = 0 and a < 0"),
            ((a < 0) & ((b > 0) | (c > 0)), "a < 0 and b > 0 or a < 0 and c > 0"),
            (((a < 0) | (b < 0)) & ((b > 0) | (c > 0)),
             "a < 0 and b > 0 or a < 0 and c > 0 or b < 0 and b > 0 or b < 0 and c > 0"),
            ((a == 0) & (b == 0), "a = 0 and b = 0"),
            (0 > a, "a < 0")]:
        assert str(selection) == text
        assert_same_part(result[selection], result[text], text)


def test_text_parses_back():
    """str(condition) parses back to the same condition, with the level sets'
    own names; analytic level sets have none yet, and the text the result
    writes, with its names, parses back as well."""
    mesh = box_mesh("hex", 2)
    a, b, c = planes(mesh)
    assert [a.name, b.name, c.name] == ["a", "b", "c"]
    own = {"a": a, "b": b, "c": c}
    for selection in [a < 0, b > 0, c == 0, (a < 0) & (b > 0), (a == 0) & (b < 0) | (b == 0) & (a < 0),
                      (a < 0) & ((b > 0) | (c == 0)), ((a < 0) | (b > 0)) & ((c < 0) | (a == 0))]:
        assert condition(str(selection), own) == selection
    assert repr((a < 0) & (c == 0)) == "Condition('a < 0 and c = 0')"
    # conditions are equal when their clauses are, the level sets compared by identity
    twin = planes(mesh, "a")[0]
    assert (a < 0) == (a < 0) and hash(a < 0) == hash(a < 0)
    assert (a < 0) != (a > 0) and (a < 0) != (twin < 0) and str(a < 0) == str(twin < 0)

    ball = cutcells.analytic_sphere(CENTRE, RADIUS)
    assert str(ball < 0).startswith("<AnalyticLevelSet at 0x")
    result = cutcells.cut(mesh, [ball, a], names=["ball", "x"])
    selection = (ball < 0) & (a > 0) | (ball > 0) & (a < 0)
    text = cutcells.conditions.selection_text(selection, [ball, a], result.level_set_names)
    assert text == "ball < 0 and x > 0 or ball > 0 and x < 0"
    assert condition(text, {"ball": ball, "x": a}) == selection
    assert_same_part(result[selection], result[text], text)


def test_comparisons_take_only_zero():
    mesh = box_mesh("hex", 2)
    for phi in (cutcells.analytic_sphere(CENTRE, RADIUS), cutcells.create_level_set(mesh, plane, degree=1, name="phi"),
                cutcells.LevelSetFunction_float32(value=lambda x: x[0], gdim=3)):
        for zero in (0, 0.0, -0.0, np.float64(0), np.int32(0)):
            assert (phi < zero) == cutcells.Condition([[(phi, "<")]])
            assert (phi > zero) == cutcells.Condition([[(phi, ">")]])
            assert (phi == zero) == cutcells.Condition([[(phi, "=")]])
        assert (0 > phi) == (phi < 0) and (0 < phi) == (phi > 0) and (0 == phi) == (phi == 0)
        for value in (1, -0.5, 1e-300, float("nan"), np.float64(2)):
            for compare in (operator.lt, operator.gt, operator.eq):
                with pytest.raises(ValueError, match="with 0, not"):
                    compare(phi, value)
        for value in (0, 1):
            with pytest.raises(ValueError, match="not a selection"):
                phi != value
        for other in ("0", None, [0], phi):
            with pytest.raises(TypeError):
                phi < other
            with pytest.raises(TypeError):
                phi > other
        with pytest.raises(TypeError):
            phi <= 0


def test_level_sets_keep_their_identity():
    mesh = box_mesh("hex", 2)
    body = cutcells.analytic_sphere(CENTRE, RADIUS)
    fluid = cutcells.create_level_set(mesh, plane, degree=1, name="fluid")
    twin = cutcells.create_level_set(mesh, plane, degree=1, name="fluid")
    assert (body == fluid) is False and (fluid == twin) is False and (fluid == fluid) is True
    assert (body != fluid) is True and (fluid != fluid) is False
    assert body in [fluid, body] and fluid in [body, fluid] and twin not in [body, fluid]
    assert [body, fluid, twin].index(twin) == 2
    assert len({body, fluid, twin, body}) == 3
    names = {body: "body", fluid: "fluid"}
    assert names[body] == "body" and names[fluid] == "fluid" and twin not in names
    assert hash(body) == object.__hash__(body) and hash(fluid) == object.__hash__(fluid)


def test_conditions_have_no_truth_value():
    body, fluid = cutcells.analytic_sphere(CENTRE, RADIUS), cutcells.analytic_sphere(CENTRE, 0.5)
    with pytest.raises(ValueError, match=r"combine conditions with & and \|"):
        (body < 0) and (fluid > 0)
    with pytest.raises(ValueError, match=r"combine conditions with & and \|"):
        (body < 0) or (fluid > 0)
    with pytest.raises(ValueError, match="no truth value"):
        bool(body == 0)
    with pytest.raises(ValueError, match="no truth value"):
        0 < body < 0
    with pytest.raises(TypeError):
        ~(body < 0)


def test_double_equals_in_strings():
    """The plan's example: "body == 0" selects as "body = 0" and as body == 0."""
    mesh = box_mesh("hex", 8)
    body, fluid = cutcells.analytic_sphere(CENTRE, 0.5), axis_plane(2, 0.2)
    result = cutcells.cut(mesh, [body, fluid], names=["body", "fluid"])
    assert_same_part(result["body == 0"], result["body = 0"], "body == 0")
    assert_same_part(result["body==0 and fluid < 0"], result["body = 0 and fluid < 0"], "body==0 and fluid < 0")
    assert_same_part(result[body == 0], result["body == 0"], "body == 0")
    assert_same_part(result[(body > 0) & (fluid < 0)], result["body > 0 and fluid < 0"], "around")
    assert_same_part(result[body < 0], result["body < 0"], "solid")


def test_selection_errors():
    mesh = box_mesh("hex", 2)
    body, other = cutcells.analytic_sphere(CENTRE, RADIUS), cutcells.analytic_sphere(CENTRE, 0.5)
    result = cutcells.cut(mesh, body, names=["body"])
    with pytest.raises(ValueError, match="not one of the cut's level sets: body"):
        result[other < 0]
    twice = cutcells.cut(mesh, [body, body], names=["a", "b"])
    with pytest.raises(ValueError, match="cut twice, as a and b"):
        twice[body < 0]
    assert_same_part(twice["a < 0"], twice["b < 0"], "one level set under two names")
    with pytest.raises(TypeError, match="condition"):
        result[3]


def test_float32_results():
    mesh = box_mesh("hex", 4)
    single = cutcells.MeshView_float32(np.array(mesh.coordinates, dtype=np.float32),
                                       np.array(mesh.connectivity, dtype=np.int32),
                                       np.array(mesh.offsets, dtype=np.int32),
                                       np.full(len(mesh.offsets) - 1, 12, dtype=np.int32), tdim=3)
    phi = cutcells.analytic_sphere(CENTRE, RADIUS)
    result = cutcells.cut(single, phi)
    assert type(result) is cutcells.part.CutResult_float32
    check(result, {"phi": phi}, ["phi < 0", "phi = 0"])
