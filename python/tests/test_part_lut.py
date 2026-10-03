# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""The lookup-table backend of cutcells.part (backend="lut"): rules on the
straight pieces of Pk-iso-P1 templates against the AdaptCell straight backend
they replace, convergence in the template order, several level sets in one
cell, faces in zero sets, visualisation and the API."""

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


# The AdaptCell straight backend (commit 091f211, before its removal) on these
# meshes and level sets: per selection the total weight, the weights summed
# with their cell index + 1, and the first moments of the reference
# coordinates. The lookup tables gave every cell the same numbers to 1e-15.
STRAIGHT = {
    ('hex', 1, 'phi < 0'): [1.198854233378003, 132.12896196701985, 0.6019069128623977, 0.5927686344785889, 0.6036663936229846],
    ('hex', 1, 'phi > 0'): [6.801145766621997, 735.8710380329802, 3.3980930871376027, 3.4072313655214117, 3.396333606377016],
    ('hex', 1, 'phi = 0'): [5.571074586758331, 614.5000436513922, 2.78503466423923, 2.781070549699648, 2.7856497897684105],
    ('hex', 2, 'phi < 0'): [1.3767665066133155, 151.66557293289415, 0.6900736130909748, 0.6831844931409687, 0.6914493065794004],
    ('hex', 2, 'phi > 0'): [6.623233493386685, 716.334427067106, 3.3099263869090256, 3.3168155068590317, 3.3085506934206],
    ('hex', 2, 'phi = 0'): [6.0136179979813065, 665.1867724962202, 2.961419662127704, 3.092951877121945, 2.927897280742848],
    ('tet', 1, 'phi < 0'): [1.198854233378003, 789.7829315059929, 0.30016750481567117, 0.2999998909011662, 0.2993915514486556],
    ('tet', 1, 'phi > 0'): [6.801145766621997, 4398.217068494007, 1.6998324951843282, 1.7000001090988333, 1.7006084485513437],
    ('tet', 1, 'phi = 0'): [5.571074586758331, 3673.0079484174494, 1.4203159838228223, 1.4005534066610066, 1.3767774021874841],
    ('tet', 2, 'phi < 0'): [1.3767665066133155, 906.5542672318867, 0.34235506079554645, 0.34176344194989294, 0.3462751560386704],
    ('tet', 2, 'phi > 0'): [6.623233493386684, 4281.445732768114, 1.657644939204453, 1.6582365580501066, 1.653724843961329],
    ('tet', 2, 'phi = 0'): [6.013617997981307, 3975.911583925143, 1.5263953316369645, 1.5384013959157317, 1.45969023217469],
    ('tet', 'crossing', 'a < 0 and b < 0'): [0.8291966635903943, 537.8656669730722, 0.20681568424864608, 0.2058298064795407, 0.2063826447726179],
    ('tet', 'crossing', 'a < 0 or b < 0'): [4.947569843022921, 3146.5587319947585, 1.2489814042090133, 1.2210734375011798, 1.1981581401976356],
    ('tet', 'crossing', 'a > 0 and b < 0'): [3.5708033364096057, 2240.0044647628715, 0.9066263434134668, 0.8793099955512867, 0.8518829841589651],
    ('tet', 'crossing', 'a = 0 and b < 0'): [3.419508185584402, 2216.5274178724694, 0.8940594420392184, 0.8489312207687633, 0.8122237481759584],
    ('tet', 'crossing', 'b = 0 and a > 0'): [2.7089435320169413, 1741.8509372073752, 0.6579659603743446, 0.7089640755267912, 0.699166466622227],
    ('quad', 1, 'phi < 0'): [1.4365102556308538, 25.669145408893293, 0.7206422473389128, 0.7103872274770804],
    ('quad', 1, 'phi > 0'): [2.563489744369146, 48.330854591106714, 1.279357752661087, 1.2896127725229196],
    ('quad', 1, 'phi = 0'): [4.291826896487456, 74.66396852556163, 2.014196321631695, 2.481221050153342],
    ('quad', 2, 'phi < 0'): [1.511722135995354, 27.028283999581657, 0.7554091293689454, 0.7540166319776491],
    ('quad', 2, 'phi > 0'): [2.4882778640046457, 46.97171600041835, 1.2445908706310544, 1.2459833680223507],
    ('quad', 2, 'phi = 0'): [4.368805284933446, 75.90820674089052, 2.092143717602884, 2.5600391777359497],
    ('quad', 'crossing', 'a < 0 and b < 0'): [0.8865588966118776, 14.149550956479796, 0.43435895140789654, 0.4370728980650539],
    ('quad', 'crossing', 'a < 0 or b < 0'): [2.825163239383476, 47.252807117175934, 1.3669761038869748, 1.4002770672459284],
    ('quad', 'crossing', 'a > 0 and b < 0'): [1.3134411033881221, 20.224523117594277, 0.6115669745180292, 0.6462604352682794],
    ('quad', 'crossing', 'a = 0 and b < 0'): [2.374002274282702, 33.272668356411316, 1.179345933984927, 1.4327420679375988],
    ('quad', 'crossing', 'b = 0 and a > 0'): [0.7114953956431181, 13.624654094572309, 0.3349462804151267, 0.3629507874291539],
    ('tri', 1, 'phi < 0'): [1.4365102556308538, 50.618437692663534, 0.4843136501787987, 0.4706975134924034],
    ('tri', 1, 'phi > 0'): [2.5634897443691465, 95.38156230733648, 0.8490196831545348, 0.86263581984093],
    ('tri', 1, 'phi = 0'): [4.291826896487456, 147.90325159786346, 1.5030575577539038, 1.509278416161476],
    ('tri', 2, 'phi < 0'): [1.5117221359953543, 53.3082754361153, 0.507677318135142, 0.4983732346113081],
    ('tri', 2, 'phi > 0'): [2.4882778640046457, 92.6917245638847, 0.8256560151981913, 0.8349600987220251],
    ('tri', 2, 'phi = 0'): [4.368805284933446, 150.01129378480672, 1.3917313473540807, 1.7099890721146227],
    ('tri', 'crossing', 'a < 0 and b < 0'): [0.8865588966118777, 27.861772362760156, 0.2948511043755897, 0.2885331972329532],
    ('tri', 'crossing', 'a < 0 or b < 0'): [2.8251632393834765, 93.13240050925258, 0.9197122779960439, 0.9340350114451238],
    ('tri', 'crossing', 'a > 0 and b < 0'): [1.3134411033881224, 39.82412507313728, 0.4120349598609019, 0.43566177683381574],
    ('tri', 'crossing', 'a = 0 and b < 0'): [2.374002274282702, 65.55611547821717, 0.7986001719081947, 0.972216962848346],
    ('tri', 'crossing', 'b = 0 and a > 0'): [0.7114953956431184, 26.91914464997395, 0.20944752696043847, 0.23931603619010286],
}


def fingerprint(rules, tdim):
    weights = per_cell(rules)
    cells = sorted(weights)
    moments = [sum(per_cell(rules, d).values()) for d in range(tdim)]
    return [sum(weights[c] for c in cells), sum((c + 1) * weights[c] for c in cells)] + moments


def total(part, order=2, mode="full", **kwargs):
    return float(np.sum(part.quadrature(order=order, mode=mode, backend="lut", **kwargs).weights))


@pytest.mark.parametrize("kind", ["hex", "tet", "quad", "tri"])
@pytest.mark.parametrize("degree", [1, 2])
def test_matches_straight_backend(kind, degree):
    """One Pk level set: lut's template is the straight backend's iso-Pk
    refinement, so the parts have its measures and first moments."""
    mesh = mesh_of(kind, 6)
    phi = cutcells.create_level_set(mesh, round_level_set(kind), degree=degree, name="phi")
    result = cutcells.part.cut(mesh, phi)
    tdim = 3 if kind in ("hex", "tet") else 2
    for expr, mode in [("phi < 0", "full"), ("phi > 0", "full"), ("phi = 0", "cut_only")]:
        lut = result[expr].quadrature(order=2, mode=mode, backend="lut")
        np.testing.assert_allclose(fingerprint(lut, tdim), STRAIGHT[(kind, degree, expr)], rtol=1e-13)


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


@pytest.mark.parametrize("kind,top", [("prism", 4), ("pyramid", 2)])
def test_template_order_on_prisms_and_pyramids(kind, top):
    """Prisms take templates of order 1 to 4 (k^2 triangles in k layers),
    pyramids 1 and 2 (six pyramids and four tetrahedra): errors fall as
    (h / k)^2."""
    from test_part import wedge_mesh

    mesh = wedge_mesh(kind, 6)
    result = cutcells.part.cut(mesh, cutcells.analytic_sphere(CENTRE, RADIUS))
    volume, area = [], []
    for k in range(1, top + 1):
        options = cutcells.LutOptions(template_order=k)
        volume.append(abs(total(result["phi < 0"], options=options) / BALL - 1))
        area.append(abs(total(result["phi = 0"], mode="cut_only", options=options) / SPHERE - 1))
        assert total(result["phi < 0 or phi > 0"], options=options) == pytest.approx(8.0, rel=1e-13)
    assert math.log(volume[0] / volume[-1]) / math.log(top) == pytest.approx(2.0, abs=0.2)
    assert math.log(area[0] / area[-1]) / math.log(top) == pytest.approx(2.0, abs=0.2)


@pytest.mark.parametrize("kind", ["prism", "pyramid"])
@pytest.mark.parametrize("degree", [1, 2])
def test_plane_on_prisms_and_pyramids(kind, degree):
    """A plane as a P1 or P2 level set is exact on the templates of prisms and
    pyramids, and the parts fill the box."""
    from test_part import wedge_mesh

    mesh = wedge_mesh(kind, 3)
    result = cutcells.part.cut(mesh, cutcells.create_level_set(mesh, plane, degree=degree, name="phi"))
    # x + 0.3 y - 0.1 < 0 in [-1, 1]^3: 4 - 0.4 (the plane moved by 0.1 along x)
    assert total(result["phi < 0"]) == pytest.approx(4.0 + 4 * 0.1, rel=1e-13)
    assert total(result["phi < 0 or phi > 0"]) == pytest.approx(8.0, rel=1e-13)
    assert total(result["phi = 0"], mode="cut_only") == pytest.approx(4.0 * math.sqrt(1.09), rel=1e-13)


@pytest.mark.parametrize("kind", ["hex", "tet", "quad", "tri"])
def test_crossing_level_sets(kind):
    """A P2 sphere (circle) and a plane crossing in cells: the four parts fill
    the box, unions are counted once, and the sphere's surface splits into its
    parts on both sides of the plane. On simplices and quadrilaterals they
    match the straight backend (which refused this on hexahedra)."""
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
    tdim = 3 if kind == "tet" else 2
    for expr, mode in [("a < 0 and b < 0", "full"), ("a < 0 or b < 0", "full"), ("a > 0 and b < 0", "full"),
                       ("a = 0 and b < 0", "cut_only"), ("b = 0 and a > 0", "cut_only")]:
        lut = result[expr].quadrature(order=2, mode=mode, backend="lut")
        np.testing.assert_allclose(fingerprint(lut, tdim), STRAIGHT[(kind, "crossing", expr)], rtol=1e-13)


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
        result["phi < 0"].visualization_mesh(backend="bogus")
    planes = cutcells.part.cut(mesh, [cutcells.create_level_set(mesh, lambda x: x[0] - 0.1, degree=1, name="a"),
                                      cutcells.create_level_set(mesh, lambda x: x[1] - 0.1, degree=1, name="b"),
                                      cutcells.create_level_set(mesh, lambda x: x[2] - 0.1, degree=1, name="c")])
    with pytest.raises(ValueError, match="three do"):
        planes["a = 0 and b = 0 and c = 0"].quadrature(backend="lut")
    # triangulated pieces give the same rules' measure
    assert total(result["phi < 0"], options=options) == pytest.approx(
        total(result["phi < 0"], options=cutcells.LutOptions(template_order=3)), rel=1e-13)


@pytest.mark.parametrize("kind", ["hex", "tet", "quad", "tri"])
def test_curves_where_two_level_sets_vanish(kind):
    """The line x = 0.1, y = -0.2 through the box (length 2), or the point
    where the two lines cross in 2D (measure 1), once."""
    mesh = mesh_of(kind, 4)
    a = cutcells.create_level_set(mesh, lambda x: x[0] - 0.1, degree=1, name="a")
    b = cutcells.create_level_set(mesh, lambda x: x[1] + 0.2, degree=1, name="b")
    part = cutcells.cut(mesh, [a, b])["a = 0 and b = 0"]
    assert part.dim == (1 if kind in ("hex", "tet") else 0)
    rules = part.quadrature(order=2)
    assert np.sum(rules.weights) == pytest.approx(2.0 if kind in ("hex", "tet") else 1.0, rel=1e-13)
    cells = part.visualization_mesh(mode="cut_only")
    assert set(np.unique(cells.vtk_types)) == ({3} if kind in ("hex", "tet") else {1})
