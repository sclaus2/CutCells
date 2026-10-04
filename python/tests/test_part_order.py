# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""The quadrature order is the polynomial degree integrated exactly on flat
pieces, in both backends: the moments of the cells of [-1, 1]^d below a plane
and on it, whole cells included, and of the line where two planes vanish."""

import itertools
import math

import numpy as np
import pytest

import cutcells

from test_part_lut import mesh_of

# phi = x0 + a . (x1, x2) - c: the plane is the graph x0 = c - a . (x1, x2),
# which stays between -0.4 and 0.6 over the box
A = np.array([0.3, -0.2])
C = 0.1


def plane_level_set(mesh):
    d = mesh.gdim
    return cutcells.create_level_set(mesh, lambda x: x[0] + sum(A[i] * x[i + 1] for i in range(d - 1)) - C, degree=1,
                                     name="phi")


def moments(weights, x, order):
    """sum_n w_n x_n^e for all exponents e up to order in each coordinate."""
    d = x.shape[1]
    powers = [np.vander(x[:, i], order + 1, increasing=True) for i in range(d)]
    return np.einsum("n," + ",".join(f"n{c}" for c in "abc"[:d]) + "->" + "abc"[:d], weights, *powers)


def exact_moments(d, order, surface):
    """The moments of {x0 < c - a . x'} in [-1, 1]^d, or of the plane in it,
    for all exponents up to order in each coordinate: the integral over x0
    in closed form, Gauss-Legendre over x' (exact)."""
    t, w = np.polynomial.legendre.leggauss(12)
    grid = np.array(list(itertools.product(t, repeat=d - 1)))
    gw = np.prod(np.array(list(itertools.product(w, repeat=d - 1))), axis=1)
    u = C - grid @ A[: d - 1]
    k = np.arange(order + 1)
    if surface:
        first = u[:, None] ** k * math.sqrt(1.0 + A[: d - 1] @ A[: d - 1])
    else:
        first = (u[:, None] ** (k + 1) - (-1.0) ** (k + 1)) / (k + 1)
    powers = [first] + [np.vander(grid[:, i], order + 1, increasing=True) for i in range(d - 1)]
    return np.einsum("n," + ",".join(f"n{c}" for c in "abc"[:d]) + "->" + "abc"[:d], gw, *powers)


def physical_points(mesh, rules):
    """The rules' points in physical coordinates: affine simplices, and
    axis-aligned quadrilaterals and hexahedra."""
    d = mesh.gdim
    coords = np.asarray(mesh.coordinates)[:, :d]
    conn, offsets = np.asarray(mesh.connectivity), np.asarray(mesh.offsets)
    xi = np.asarray(rules.points).reshape(-1, d)
    x = np.empty_like(xi)
    offset = np.asarray(rules.offset)
    for k, c in enumerate(np.asarray(rules.parent_map)):
        v = coords[conn[offsets[c]:offsets[c + 1]]]
        p = slice(offset[k], offset[k + 1])
        if len(v) == d + 1:
            x[p] = v[0] + xi[p] @ (v[1:] - v[0])
        else:
            x[p] = v.min(axis=0) + xi[p] * (v.max(axis=0) - v.min(axis=0))
    return x


@pytest.mark.parametrize("kind", ["tet", "hex", "tri", "quad"])
@pytest.mark.parametrize("backend", ["quadrays", "lut"])
@pytest.mark.parametrize("surface", [False, True])
def test_order_is_the_exact_degree(kind, backend, surface):
    mesh = mesh_of(kind, 3)
    d = mesh.gdim
    result = cutcells.cut(mesh, plane_level_set(mesh))
    part = result["phi = 0" if surface else "phi < 0"]
    for order in range(1, 11):
        rules = part.quadrature(order=order, backend=backend)
        x = physical_points(mesh, rules)
        weights = np.asarray(rules.weights)
        # tetrahedra skip the reference rules of degree 3, 7 and 8, which have a negative weight
        assert np.all(weights > 0), order
        degree = np.sum(np.indices((order + 1,) * d), axis=0)
        np.testing.assert_allclose(moments(weights, x, order)[degree <= order],
                                   exact_moments(d, order, surface)[degree <= order], rtol=1e-12, atol=1e-12,
                                   err_msg=f"order {order}")


# two planes in general position, NA . x = DA and NB . x = DB; the curve where
# both vanish is a line through [-1, 1]^3
NA, DA = np.array([1.0, 0.3, -0.2]), 0.1
NB, DB = np.array([-0.25, 1.0, 0.4]), -0.05


def line_moments(order):
    """The moments of the planes' line inside [-1, 1]^3, for all exponents up
    to order in each coordinate: Gauss-Legendre in its parameter (exact for
    total degrees up to 23)."""
    t = np.cross(NA, NB)
    p0 = (DA * np.cross(NB, t) + DB * np.cross(t, NA)) / (t @ t)
    lo, hi = -np.inf, np.inf
    for i in range(3):
        s1, s2 = sorted(((-1.0 - p0[i]) / t[i], (1.0 - p0[i]) / t[i]))
        lo, hi = max(lo, s1), min(hi, s2)
    g, w = np.polynomial.legendre.leggauss(12)
    s = lo + (hi - lo) * (g + 1) / 2
    return moments(w * (hi - lo) / 2 * np.linalg.norm(t), p0 + s[:, None] * t, order)


@pytest.mark.parametrize("kind", ["tet", "hex"])
@pytest.mark.parametrize("backend", ["quadrays", "lut"])
def test_order_on_curves(kind, backend):
    """On the line where two planes vanish, order is the polynomial degree
    integrated exactly too: quadrays' points per segment, ceil((order + 1) / 2),
    on a line along which the parameter is affine."""
    mesh = mesh_of(kind, 3)
    level_sets = [cutcells.create_level_set(mesh, lambda x, n=n, d=d: n[0] * x[0] + n[1] * x[1] + n[2] * x[2] - d,
                                            degree=1, name=name) for n, d, name in ((NA, DA, "a"), (NB, DB, "b"))]
    part = cutcells.cut(mesh, level_sets)["a = 0 and b = 0"]
    degree = np.sum(np.indices((11,) * 3), axis=0)
    for order in range(1, 11):
        rules = part.quadrature(order=order, backend=backend)
        weights = np.asarray(rules.weights)
        assert np.all(weights > 0), order
        mask = degree[: order + 1, : order + 1, : order + 1] <= order
        np.testing.assert_allclose(moments(weights, physical_points(mesh, rules), order)[mask],
                                   line_moments(order)[mask], rtol=1e-12, atol=1e-12, err_msg=f"order {order}")


def test_points_per_segment():
    """order 1 integrates the volume of a planar cut exactly; one Gauss point
    per segment, as QuadraysOptions.points_per_segment can ask, does not."""
    mesh = mesh_of("tet", 1)
    result = cutcells.cut(mesh, plane_level_set(mesh))
    exact = exact_moments(3, 0, False)[0, 0, 0]
    assert np.sum(result["phi < 0"].quadrature(order=1).weights) == pytest.approx(exact, rel=1e-13)
    options = cutcells.QuadraysOptions()
    options.points_per_segment = 1
    assert abs(np.sum(result["phi < 0"].quadrature(order=1, options=options).weights) - exact) > 1e-3


@pytest.mark.parametrize("backend", ["quadrays", "lut"])
def test_order_range(backend):
    mesh = mesh_of("tri", 2)
    part = cutcells.cut(mesh, plane_level_set(mesh))["phi < 0"]
    for order in (0, 11):
        with pytest.raises(ValueError, match="1 to 10"):
            part.quadrature(order=order, backend=backend)
