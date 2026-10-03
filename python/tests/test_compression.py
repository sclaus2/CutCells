# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Compression of quadrature rules (compress_rules): quadrays rules of a ball
on hexahedra and tetrahedra keep their polynomial moments on at most as many
points as the space has moments, with positive weights."""

import math

import numpy as np
import pytest

import cutcells

from test_quadrays import CENTRE, RADIUS, box_mesh


def moments(rules, degree, space):
    """Monomial moments of every rule in reference coordinates, shape (rules, monomials)."""
    tdim = rules.tdim
    points = np.asarray(rules.points).reshape(-1, tdim)
    weights = np.asarray(rules.weights)
    offset = np.asarray(rules.offset)
    exps = [e for e in np.ndindex(*(tdim * [degree + 1])) if space == "tensor" or sum(e) <= degree]
    values = np.stack([np.prod(points ** np.array(e), axis=1) for e in exps], axis=1)
    return np.add.reduceat(values * weights[:, None], offset[:-1], axis=0), np.add.reduceat(np.abs(weights), offset[:-1])


@pytest.mark.parametrize(
    "kind, space, degree, n_moments",
    [("hex", "tensor", 4, 125), ("hex", "tensor", 2, 27), ("tet", "total", 4, 35), ("tet", "total", 2, 10)],
)
def test_compressed_quadrays_rules(kind, space, degree, n_moments):
    mesh = box_mesh(kind, 4)
    result = cutcells.cut(mesh, cutcells.analytic_sphere(CENTRE, RADIUS))
    rules = result["phi < 0"].quadrature(order=5, mode="cut_only", backend="quadrays")
    compressed, stats = cutcells.compress_rules(rules, degree, space)

    counts = np.diff(np.asarray(compressed.offset))
    assert np.array_equal(np.asarray(compressed.parent_map), np.asarray(rules.parent_map))
    original = np.diff(np.asarray(rules.offset))
    assert np.all((counts <= n_moments) | (counts == original))
    assert np.all(np.asarray(compressed.weights) > 0)
    assert stats.n_rules == len(counts) and stats.n_skipped == 0 and stats.n_compressed > 0
    assert stats.points_after < stats.points_before
    assert stats.max_residual < 1e-13

    m0, scale = moments(rules, degree, space)
    m1, _ = moments(compressed, degree, space)
    assert np.max(np.abs(m1 - m0) / scale[:, None]) < 1e-12
    volume = 4.0 / 3.0 * math.pi * RADIUS**3
    whole = result["phi < 0"].quadrature(order=5, mode="full", backend="quadrays")
    uncut = np.sum(np.asarray(whole.weights)) - np.sum(np.asarray(rules.weights))
    assert abs(np.sum(np.asarray(compressed.weights)) + uncut - volume) < 1e-6


def test_small_rules_are_copied():
    rules, _ = cutcells.quadrays_cell_rules(
        cutcells.CellType.hexahedron,
        np.array([0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1], dtype=float),
        cutcells.analytic_sphere([0.0, 0.0, 0.0], 0.5),
        "phi < 0",
        q=2,
    )
    assert len(rules.weights) <= 11**3
    compressed, stats = cutcells.compress_rules(rules, 10)
    assert len(compressed.weights) == len(rules.weights)
    assert stats.n_compressed == 0


def test_float64_name_and_unknown_space():
    mesh = box_mesh("hex", 3)
    result = cutcells.cut(mesh, cutcells.analytic_sphere(CENTRE, RADIUS))
    rules = result["phi < 0"].quadrature(order=5, mode="cut_only", backend="quadrays")
    compressed, _ = cutcells.compress_rules_float64(rules, 4, "tensor")
    assert np.diff(np.asarray(compressed.offset)).max() <= 125
    with pytest.raises(ValueError):
        cutcells.compress_rules(rules, 4, "simplex")
