# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Accuracy of quadrays on the parts two curved level sets select.

Two balls of radius 0.5 meet in a lens. Its volume, the caps of each sphere
inside the other ball, and the circle (2D: the two points) where the spheres
meet are compared cell by cell with a reference of 14 Gauss points per
segment: whole-mesh totals hide errors that cancel between cells sharing a
face. With exact level sets the cells' errors are rounding; with P3
interpolants the quadrature error stays far below the interpolation error,
which falls at the interpolation rate h^4.
"""

import math

import numpy as np
import pytest

import cutcells

R = 0.5
C1 = (-0.21, 0.07, 0.04)
C2 = (0.24, 0.03, -0.02)


def exact_values(gdim):
    d = math.dist(C1[:gdim], C2[:gdim])
    if gdim == 2:
        th = math.acos(d / (2 * R))
        return {"lens": 2 * R * R * th - d / 2 * math.sqrt(4 * R * R - d * d), "cap": 2 * R * th, "curve": 2.0}
    h, a = R - d / 2, math.sqrt(R * R - d * d / 4)
    return {"lens": math.pi * (4 * R + d) * (2 * R - d) ** 2 / 12, "cap": 2 * math.pi * R * h,
            "curve": 2 * math.pi * a}


def mesh_of(kind, n):
    if kind in ("triangle", "quadrilateral"):
        return cutcells.rectangle_mesh(-1.0, -1.0, 1.0, 1.0, n, n, cell_type=kind)
    return cutcells.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n, n, n, cell_type=kind)


def ball(centre, gdim, squared):
    def phi(X):
        r2 = sum((X[i] - centre[i]) ** 2 for i in range(gdim))
        return r2 - R * R if squared else np.sqrt(r2) - R

    return phi


def parts(a, b):
    return {"lens": (a < 0) & (b < 0), "cap": (a == 0) & (b < 0), "curve": (a == 0) & (b == 0)}


def cell_totals(rules):
    """The sum of the weights of each cell."""
    offsets, weights = np.asarray(rules.offset), np.asarray(rules.weights)
    totals = {}
    for i, c in enumerate(np.asarray(rules.parent_map)):
        totals[int(c)] = totals.get(int(c), 0.0) + weights[offsets[i]:offsets[i + 1]].sum()
    return totals


def cell_error(totals, reference):
    """The sum over the cells of the error of their measure."""
    return sum(abs(totals.get(c, 0.0) - reference.get(c, 0.0)) for c in set(totals) | set(reference))


def reference_options():
    options = cutcells.QuadraysOptions()
    options.points_per_segment = 14
    return options


@pytest.mark.parametrize("kind, n, order", [("triangle", 17, 8), ("quadrilateral", 17, 8), ("tetrahedron", 5, 6),
                                            ("hexahedron", 5, 6)])
def test_exact_spheres_are_accurate_in_every_cell(kind, n, order):
    """Quadratic level sets |x - c|^2 - r^2 are exact as P2 interpolants. In
    2D the cells' errors are rounding. On the coarse 3D meshes the accuracy
    check bisects until the rule's own rate is met, (h kappa / 2)^(2 q): the
    lens is rounding, the caps and the circle within 1e-10 and 1e-9 (without
    the check and the anchoring at folds, the caps err by 6e-9 to 2e-7). The
    3D lens, whose reference costs most, is checked by its total."""
    gdim = 2 if kind in ("triangle", "quadrilateral") else 3
    bounds = {"lens": 1e-13, "cap": 1e-13, "curve": 1e-13} if gdim == 2 else {"lens": 5e-13, "cap": 1e-10,
                                                                              "curve": 1e-9}
    mesh = mesh_of(kind, n)
    a = cutcells.create_level_set(mesh, ball(C1, gdim, True), degree=2, name="a")
    b = cutcells.create_level_set(mesh, ball(C2, gdim, True), degree=2, name="b")
    result = cutcells.cut(mesh, [a, b])
    exact = exact_values(gdim)
    for name, selection in parts(a, b).items():
        part = result[selection]
        rules = part.quadrature(order=order)
        assert np.all(np.asarray(rules.weights) > 0)
        assert abs(np.sum(rules.weights) - exact[name]) < bounds[name], name
        if gdim == 3 and name == "lens":
            continue
        reference = cell_totals(part.quadrature(order=order, options=reference_options()))
        assert cell_error(cell_totals(rules), reference) < bounds[name], name


@pytest.mark.parametrize("kind", ["triangle", "quadrilateral"])
def test_p3_spheres_converge_at_the_interpolation_rate(kind):
    """P3 interpolants of the distances |x - c| - r: the quadrature error of
    each cell stays far below its interpolation error (the exact spheres'
    cell values), and the cells' interpolation errors of the lens fall as h^4."""
    exact_errors = []
    for n in (9, 17, 33):
        mesh = mesh_of(kind, n)
        a = cutcells.create_level_set(mesh, ball(C1, 2, False), degree=3, name="a")
        b = cutcells.create_level_set(mesh, ball(C2, 2, False), degree=3, name="b")
        ea = cutcells.analytic_sphere(list(C1[:2]) + [0.0], R, name="a")
        eb = cutcells.analytic_sphere(list(C2[:2]) + [0.0], R, name="b")
        result, exact_result = cutcells.cut(mesh, [a, b]), cutcells.cut(mesh, [ea, eb])
        for name in ("lens", "cap"):
            part = result[parts(a, b)[name]]
            totals = cell_totals(part.quadrature(order=6))
            reference = cell_totals(part.quadrature(order=6, options=reference_options()))
            exact = cell_totals(exact_result[parts(ea, eb)[name]].quadrature(order=6, options=reference_options()))
            quadrature, interpolation = cell_error(totals, reference), cell_error(reference, exact)
            assert quadrature < 1e-3 * interpolation, (name, n, quadrature, interpolation)
            if name == "lens":
                exact_errors.append(interpolation)
    rates = [math.log2(e0 / e1) for e0, e1 in zip(exact_errors, exact_errors[1:])]
    assert min(rates) > 3.3, (exact_errors, rates)
