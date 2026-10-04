# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Level sets made of level sets: LevelSetFunction(level_sets=..., inside=...)
is cut as its level sets, named after it ("ring.ball"), and selections of it
expand into selections of them.

- The lens, the bicylinder, the tricylinder and the napkin ring of test_part.py,
  built as level sets made of level sets, select the parts of their
  hand-written selections, with their totals: inside, outside and faces.
- The napkin ring's faces: .on("cyl") gives 2 pi r h and .on("ball") 2 pi R h.
- Two boxes side by side, whose shared plane is one level set used with both
  signs, have no surface on that plane; a smaller second box leaves the ring
  around it; both backends.
- Nested insides (parentheses) and level sets made of level sets made of level
  sets give the parts of their hand-written expansions.
- Conditions and strings, short and qualified names, names of analytic level
  sets, float32, and the errors."""

import math
import re

import numpy as np
import pytest

import cutcells

from test_conditions import assert_same_part
from test_part import TRICYLINDER_FACES, ball_level_set, cylinder_level_set, total
from test_quadrays import CENTRE, box_mesh


def made_of(level_sets, inside, name):
    return cutcells.LevelSetFunction(level_sets=level_sets, inside=inside, name=name)


def qualified(text, owner):
    """A selection string over a, b, c with the full names of owner's level sets."""
    return re.sub(r"\b([abc])\b", rf"{owner}.\1", text)


def lens(mesh, ls_kind):
    return [ball_level_set(mesh, ls_kind, CENTRE + np.array([-0.25, 0.0, 0.0]), 0.6, "a"),
            ball_level_set(mesh, ls_kind, CENTRE + np.array([0.3, 0.05, 0.0]), 0.5, "b")]


def bicylinder(mesh, ls_kind):
    return [cylinder_level_set(mesh, 0.6, 2, "a"), cylinder_level_set(mesh, 0.6, 1, "b")]


def tricylinder(mesh, ls_kind):
    return [cylinder_level_set(mesh, 0.6, axis, name) for axis, name in zip((2, 1, 0), "abc")]


def napkin_ring(mesh, ls_kind):
    return [ball_level_set(mesh, ls_kind, CENTRE, 0.8, "a"), cylinder_level_set(mesh, 0.45, 2, "b")]


# the shape, its inside, and the hand-written selections of its outside and faces
SHAPES = {
    "lens": (lens, "a < 0 and b < 0", "a > 0 or b > 0", "a = 0 and b < 0 or b = 0 and a < 0"),
    "bicylinder": (bicylinder, "a < 0 and b < 0", "a > 0 or b > 0", "a = 0 and b < 0 or b = 0 and a < 0"),
    "tricylinder": (tricylinder, "a < 0 and b < 0 and c < 0", "a > 0 or b > 0 or c > 0",
                    " or ".join(TRICYLINDER_FACES)),
    "napkin ring": (napkin_ring, "a < 0 and b > 0", "a > 0 or b < 0", "a = 0 and b > 0 or b = 0 and a < 0"),
}


# balls P2 or analytic; cylinders P2
CASES = [(shape, mesh_kind, ls_kind) for shape in SHAPES
         for mesh_kind, ls_kind in ([("hex", "P2"), ("hex", "analytic"), ("tet", "analytic")]
                                    if shape in ("lens", "napkin ring") else [("hex", "P2"), ("tet", "P2")])]


@pytest.mark.parametrize("shape,mesh_kind,ls_kind", CASES)
def test_shapes_as_level_sets_made_of_level_sets(shape, mesh_kind, ls_kind):
    """body < 0, body > 0 and body == 0 select the parts of the hand-written
    selections of the inside, the outside and the faces, with their totals."""
    make, inside, outside, faces = SHAPES[shape]
    mesh = box_mesh(mesh_kind, 8)
    level_sets = make(mesh, ls_kind)
    names = list("abc"[:len(level_sets)])
    plain = cutcells.cut(mesh, level_sets, names=names)
    # analytic balls have no names: as the values of named LevelSetFunctions
    named = [cutcells.LevelSetFunction(value=ls, name=name) if isinstance(ls, cutcells.AnalyticLevelSet) else ls
             for ls, name in zip(level_sets, names)]
    body = made_of(named, inside, "body")
    result = cutcells.cut(mesh, body)
    assert result.level_set_names == [f"body.{name}" for name in names]
    for selection, text, mode in ((body < 0, inside, "full"), (body > 0, outside, "full"),
                                  (body == 0, faces, "cut_only")):
        part = result[selection]
        assert_same_part(part, result[qualified(text, "body")], text)
        assert total(part, mode=mode) == pytest.approx(total(plain[text], mode=mode), rel=1e-12), text


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
def test_napkin_ring_faces(mesh_kind):
    """The napkin ring made of a ball and a cylinder: its volume, and its faces
    on either, which depend only on its height h."""
    R, r = 0.8, 0.45
    h = 2 * math.sqrt(R**2 - r**2)
    mesh = box_mesh(mesh_kind, 8)
    ball = cutcells.analytic_sphere(CENTRE, R, signed_distance=False, name="ball")
    cyl = cylinder_level_set(mesh, r, 2, "cyl")
    ring = made_of([ball, cyl], (ball < 0) & (cyl > 0), "ring")
    result = cutcells.cut(mesh, ring)
    assert total(result[ring < 0]) == pytest.approx(math.pi * h**3 / 6, rel=1e-7)
    surface = result[ring == 0]
    assert total(surface.on("cyl"), mode="cut_only") == pytest.approx(2 * math.pi * r * h, rel=1e-6)
    assert total(surface.on("ball"), mode="cut_only") == pytest.approx(2 * math.pi * R * h, rel=1e-6)
    assert total(surface.on("ring.cyl"), mode="cut_only") == total(surface.on("cyl"), mode="cut_only")
    assert total(surface.on("ball", "cyl"), mode="cut_only") == pytest.approx(total(surface, mode="cut_only"),
                                                                               rel=1e-14)
    assert total(surface.on("ring"), mode="cut_only") == total(surface, mode="cut_only")
    assert_same_part(surface.on("cyl"), result["cyl = 0 and ball < 0"], "the cylinder's face")
    assert surface.on("cyl").dim == 2


def planes(mesh, specs):
    """P1 planes x_axis = offset, named."""
    return [cutcells.create_level_set(mesh, lambda x, i=axis, o=offset: x[i] - o, degree=1, name=name)
            for axis, offset, name in specs]


BOX = "ylo > 0 and yhi < 0 and zlo > 0 and zhi < 0"


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
@pytest.mark.parametrize("shared", [0.1, 0.0])
def test_two_boxes_sharing_a_plane(mesh_kind, shared):
    """Two boxes side by side, [-0.7, p] x Y x Z and [p, 0.7] x Y x Z with
    Y = Z = [-0.6, 0.45], their union made of planes, the shared one p used
    with both signs (at x = 0 it lies in mesh faces): no surface on p, and the
    union's volume and surface, with quadrays and the lookup tables."""
    mesh = box_mesh(mesh_kind, 8)
    level_sets = planes(mesh, [(0, -0.7, "xlo"), (0, shared, "p"), (0, 0.7, "xhi"), (1, -0.6, "ylo"),
                               (1, 0.45, "yhi"), (2, -0.6, "zlo"), (2, 0.45, "zhi")])
    union = made_of(level_sets, f"(xlo > 0 and p < 0 and {BOX}) or (p > 0 and xhi < 0 and {BOX})", "union")
    result = cutcells.cut(mesh, union)
    w = 1.05
    surface = result[union == 0]
    assert len(surface.on("p").quadrature(order=2, mode="cut_only").weights) == 0
    for backend in ("quadrays", "lut"):
        volume = float(np.sum(result[union < 0].quadrature(order=2, backend=backend).weights))
        area = float(np.sum(surface.quadrature(order=2, mode="cut_only", backend=backend).weights))
        assert volume == pytest.approx(1.4 * w * w, rel=1e-12), backend
        assert area == pytest.approx(2 * (2 * 1.4 * w + w * w), rel=1e-12), backend


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
def test_smaller_second_box(mesh_kind):
    """A smaller second box, Y = Z = [-0.35, 0.3]: the shared plane keeps the
    ring of the first box's face around it, with the first box's normal."""
    mesh = box_mesh(mesh_kind, 8)
    level_sets = planes(mesh, [(0, -0.7, "xlo"), (0, 0.1, "p"), (0, 0.7, "xhi"), (1, -0.6, "ylo"), (1, 0.45, "yhi"),
                               (2, -0.6, "zlo"), (2, 0.45, "zhi"), (1, -0.35, "ylo2"), (1, 0.3, "yhi2"),
                               (2, -0.35, "zlo2"), (2, 0.3, "zhi2")])
    second = "p > 0 and xhi < 0 and ylo2 > 0 and yhi2 < 0 and zlo2 > 0 and zhi2 < 0"
    union = made_of(level_sets, f"(xlo > 0 and p < 0 and {BOX}) or ({second})", "union")
    result = cutcells.cut(mesh, union)
    w, a2 = 1.05, 0.65 * 0.65
    assert total(result[union < 0], order=2) == pytest.approx(0.8 * w * w + 0.6 * a2, rel=1e-12)
    surface = result[union == 0]
    assert total(surface.on("p"), order=2, mode="cut_only") == pytest.approx(w * w - a2, rel=1e-12)
    assert total(surface, order=2, mode="cut_only") == pytest.approx(
        2 * (2 * 0.8 * w + w * w) + 2 * (2 * 0.6 * 0.65 + a2) - 2 * a2, rel=1e-12)


@pytest.mark.parametrize("mesh_kind", ["hex", "tet"])
def test_nested_insides(mesh_kind):
    """An inside with parentheses, (a < 0 or b < 0) and c > 0, two balls above
    a plane, selects the parts of its expansions written out."""
    mesh = box_mesh(mesh_kind, 8)
    a, b = lens(mesh, "P2")
    c = planes(mesh, [(2, 0.13, "c")])[0]
    body = made_of([a, b, c], "(a < 0 or b < 0) and c > 0", "body")
    result = cutcells.cut(mesh, body)
    for selection, text, mode in (
            (body < 0, "a < 0 and c > 0 or b < 0 and c > 0", "full"),
            (body > 0, "a > 0 and b > 0 or c < 0", "full"),
            (body == 0, "a = 0 and b > 0 and c > 0 or b = 0 and a > 0 and c > 0 or c = 0 and a < 0 or "
                        "c = 0 and b < 0", "cut_only")):
        hand = result[qualified(text, "body")]
        assert_same_part(result[selection], hand, text)
        assert total(result[selection], mode=mode) == pytest.approx(total(hand, mode=mode), rel=1e-12), text
    # the same inside as a condition, and strings with parentheses
    assert_same_part(result[body < 0], result["(body.a < 0 or body.b < 0) and body.c > 0"], "parentheses")
    twin = made_of([a, b, c], ((a < 0) | (b < 0)) & (c > 0), "body")
    assert twin.inside == "a < 0 and c > 0 or b < 0 and c > 0"
    assert_same_part(cutcells.cut(mesh, twin)[twin == 0], result[body == 0], "an inside written by a condition")


def test_level_sets_made_of_level_sets_made_of_level_sets():
    """outer, made of a lens (made of two balls), a third ball c and a plane,
    with inside lens < 0 or (c < 0 and plane > 0): its parts are those written
    out with the names of its level sets, nested."""
    mesh = box_mesh("hex", 8)
    lens_ = made_of(lens(mesh, "P2"), "a < 0 and b < 0", "lens")
    c = ball_level_set(mesh, "P2", CENTRE + np.array([0.0, 0.4, 0.1]), 0.35, "c")
    plane = planes(mesh, [(2, -0.2, "plane")])[0]
    outer = made_of([lens_, c, plane], (lens_ < 0) | (c < 0) & (plane > 0), "outer")
    assert outer.inside == "lens < 0 or c < 0 and plane > 0"
    result = cutcells.cut(mesh, outer)
    assert result.level_set_names == ["outer.lens.a", "outer.lens.b", "outer.c", "outer.plane"]
    assert result.insides == {"outer": "lens < 0 or c < 0 and plane > 0", "outer.lens": "a < 0 and b < 0"}
    a, b, c, p = result.level_set_names
    for selection, text, mode in (
            ("outer < 0", f"{a} < 0 and {b} < 0 or {c} < 0 and {p} > 0", "full"),
            ("outer > 0", f"{a} > 0 and {c} > 0 or {a} > 0 and {p} < 0 or {b} > 0 and {c} > 0 or {b} > 0 and {p} < 0",
             "full"),
            ("lens = 0", f"{a} = 0 and {b} < 0 or {b} = 0 and {a} < 0", "cut_only"),
            ("outer = 0", f"{a} = 0 and {b} < 0 and {c} > 0 or {a} = 0 and {b} < 0 and {p} < 0 or "
                          f"{b} = 0 and {a} < 0 and {c} > 0 or {b} = 0 and {a} < 0 and {p} < 0 or "
                          f"{c} = 0 and {a} > 0 and {p} > 0 or {c} = 0 and {b} > 0 and {p} > 0 or "
                          f"{p} = 0 and {a} > 0 and {c} < 0 or {p} = 0 and {b} > 0 and {c} < 0", "cut_only")):
        assert_same_part(result[selection], result[text], selection)
        assert total(result[selection], mode=mode) == pytest.approx(total(result[text], mode=mode), rel=1e-12)
    # faces of a level set made of level sets within: on("lens") keeps both of its faces
    surface = result["outer = 0"]
    assert_same_part(surface.on("lens"), surface.on("a", "b"), "on a level set made of level sets")
    assert_same_part(surface.on("lens.a"), surface.on("outer.lens.a"), "on a name within a name")


def test_conditions_and_names():
    """Conditions name a level set made of level sets as cut names it, strings
    take its name and those of its level sets, in full or short; the result and
    the level set say what they are made of."""
    mesh = box_mesh("hex", 4)
    ball = cutcells.analytic_sphere(CENTRE, 0.8, signed_distance=False, name="ball")
    cyl = cylinder_level_set(mesh, 0.45, 2, "cyl")
    fluid = planes(mesh, [(2, 0.2, "fluid")])[0]
    ring = made_of([ball, cyl], (ball < 0) & (cyl > 0), "ring")
    assert ring.name == "ring" and ring.inside == "ball < 0 and cyl > 0" and ring.has_level_sets()
    assert [ls.name for ls in ring.level_sets] == ["ball", "cyl"]
    assert fluid.level_sets == [] and fluid.inside is None and not fluid.has_level_sets()
    result = cutcells.cut(mesh, [ring, fluid], names=["body", "fluid"])
    assert result.level_set_names == ["body.ball", "body.cyl", "fluid"]
    assert result.insides == {"body": "ball < 0 and cyl > 0"}
    assert result.num_level_sets == 3 and result.domains.shape == (3, mesh.num_cells())
    for selection, text in ((ring < 0, "body < 0"), ((ring > 0) & (fluid < 0), "body > 0 and fluid < 0"),
                            (ring == 0, "body == 0"), ((ring == 0) & (fluid > 0), "body = 0 and fluid > 0"),
                            ((ring == 0) & (fluid == 0), "body = 0 and fluid = 0")):
        assert_same_part(result[selection], result[text], text)
    assert_same_part(result["cyl < 0"], result["body.cyl < 0"], "a short name")
    assert_same_part(result["body = 0"].on("cyl"), result["body = 0"].on("body.cyl"), "on a short name")
    assert result[(ring == 0) & (fluid == 0)].dim == 1
    # two of the same short name
    holes = [made_of([cutcells.analytic_sphere(CENTRE, radius, name="hole"), fluid], "hole > 0 and fluid < 0", name)
             for radius, name in ((0.3, "s"), (0.5, "t"))]
    both = cutcells.cut(mesh, holes)
    with pytest.raises(RuntimeError, match="could be any of s.hole, t.hole"):
        both["hole < 0"]
    assert_same_part(both["s = 0"].on("hole"), both["s = 0"].on("s.hole"), "a short name among the part's faces")
    with pytest.raises(ValueError, match="could be any of s.hole, t.hole"):
        both["s = 0 or t = 0"].on("hole")


def test_analytic_level_sets_have_names():
    mesh = box_mesh("hex", 2)
    ball = cutcells.analytic_sphere(CENTRE, 0.5, name="ball")
    assert ball.name == "ball" and str(ball < 0) == "ball < 0"
    assert cutcells.cut(mesh, ball).level_set_names == ["ball"]
    assert cutcells.cut(mesh, ball, names=["other"]).level_set_names == ["other"]
    unnamed = cutcells.analytic_sphere(CENTRE, 0.5)
    assert unnamed.name is None and str(unnamed < 0).startswith("<AnalyticLevelSet at 0x")
    assert cutcells.cut(mesh, [ball, unnamed]).level_set_names == ["ball", "phi2"]
    assert cutcells.AnalyticLevelSet.from_capsule(ball.capsule, name="copy").name == "copy"
    plane = cutcells.AnalyticLevelSet(lambda x: float(x[0]), lambda x: [1.0, 0.0, 0.0],
                                      lambda lo, hi: [lo[0], hi[0], 1.0, 1.0, 0.0, 0.0, 0.0, 0.0], name="plane")
    assert plane.name == "plane"
    # a LevelSetFunction with an analytic value is that analytic level set
    wrapped = cutcells.LevelSetFunction(value=ball)
    assert wrapped.name == "ball" and wrapped.value(np.array([*CENTRE])) == pytest.approx(-0.5)
    fine = box_mesh("hex", 8)
    assert total(cutcells.cut(fine, wrapped)["ball < 0"]) == pytest.approx(total(cutcells.cut(fine, ball)["ball < 0"]),
                                                                          rel=1e-14)


def test_float32():
    mesh = box_mesh("hex", 4)
    single = cutcells.MeshView_float32(np.array(mesh.coordinates, dtype=np.float32),
                                       np.array(mesh.connectivity, dtype=np.int32),
                                       np.array(mesh.offsets, dtype=np.int32),
                                       np.full(len(mesh.offsets) - 1, 12, dtype=np.int32), tdim=3)
    a = cutcells.analytic_sphere(CENTRE, 0.6, signed_distance=False, name="a")
    b = cutcells.analytic_sphere(CENTRE + np.array([0.3, 0.0, 0.0]), 0.5, signed_distance=False, name="b")
    lens32 = cutcells.LevelSetFunction_float32(level_sets=[a, b], inside=(a < 0) & (b < 0), name="lens")
    result = cutcells.cut(single, lens32)
    assert type(result) is cutcells.part.CutResult_float32
    assert_same_part(result[lens32 < 0], result["lens.a < 0 and lens.b < 0"], "float32")
    with pytest.raises(TypeError, match="LevelSetFunction_float64s or AnalyticLevelSets"):
        cutcells.LevelSetFunction(level_sets=[lens32], inside="lens < 0", name="x")


def test_errors():
    mesh = box_mesh("hex", 2)
    ball = cutcells.analytic_sphere(CENTRE, 0.5, name="ball")
    cyl = cylinder_level_set(mesh, 0.3, 2, "cyl")
    with pytest.raises(ValueError, match="need names"):
        made_of([cutcells.analytic_sphere(CENTRE, 0.5), cyl], "cyl < 0", "x")
    with pytest.raises(ValueError, match="two level sets of 'x' are named 'cyl'"):
        made_of([cyl, cyl], "cyl < 0", "x")
    with pytest.raises(ValueError, match="names 'hole', none of its level sets"):
        made_of([ball, cyl], "ball < 0 and hole > 0", "x")
    with pytest.raises(ValueError, match="an inside is a region"):
        made_of([ball, cyl], "ball = 0", "x")
    with pytest.raises(ValueError, match="right-hand side"):
        made_of([ball, cyl], "(ball < 1)", "x")
    with pytest.raises(ValueError, match="give both"):
        cutcells.LevelSetFunction(value=lambda x: x[0], inside="ball < 0")
    with pytest.raises(ValueError, match="give the inside"):
        cutcells.LevelSetFunction(level_sets=[ball, cyl], name="x")
    with pytest.raises(TypeError, match="AnalyticLevelSets"):
        made_of([ball, "cyl"], "ball < 0", "x")
    with pytest.raises(ValueError, match="is not one of the level sets of x: ball, cyl"):
        made_of([ball, cyl], cutcells.analytic_sphere(CENTRE, 0.2) < 0, "x")
    with pytest.raises(RuntimeError, match="requires at least one of value, nodal_values or level_sets"):
        cutcells.LevelSetFunction()
    ring = made_of([ball, cyl], (ball < 0) & (cyl > 0), "ring")
    with pytest.raises(ValueError, match="needs a name that selections can use"):
        cutcells.cut(mesh, ring, names=["my ring"])
    with pytest.raises(ValueError, match="would read as one of 'ring''s own"):
        cutcells.cut(mesh, [ring, planes(mesh, [(0, 0.1, "ring.x")])[0]])
    result = cutcells.cut(mesh, ring)
    with pytest.raises(ValueError, match="a volume part has no faces"):
        result[ring < 0].on("cyl")
    with pytest.raises(ValueError, match="no level set is named 'hole'"):
        result[ring == 0].on("hole")
    with pytest.raises(ValueError, match="none of the part's faces"):
        result[ring == 0].on("ball").on("cyl")
    with pytest.raises(ValueError, match="name one or more faces"):
        result[ring == 0].on()
    with pytest.raises(TypeError, match="strings"):
        result[ring == 0].on(cyl)
    many = planes(mesh, [(0, 0.01 * i, f"p{i}") for i in range(65)])
    with pytest.raises(ValueError, match="give 1 to 64 level sets, counting"):
        cutcells.cut(mesh, made_of(many, "p0 < 0", "x"))
    assert cutcells.cut(mesh, made_of(many[:64], "p0 < 0", "x")).num_level_sets == 64


def test_visualization():
    mesh = box_mesh("hex", 4)
    ball = cutcells.analytic_sphere(CENTRE, 0.8, signed_distance=False, name="ball")
    cyl = cylinder_level_set(mesh, 0.45, 2, "cyl")
    ring = made_of([ball, cyl], (ball < 0) & (cyl > 0), "ring")
    result = cutcells.cut(mesh, ring)
    faces = result[ring == 0].visualization_mesh(mode="cut_only", degree=2)
    hole = result[ring == 0].on("cyl").visualization_mesh(mode="cut_only", degree=2)
    assert faces.n_cells() > hole.n_cells() > 0
    assert len(result[ring < 0].visualization_mesh(mode="cut_only", backend="lut").vtk_types) > 0
