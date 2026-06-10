# Copyright (c) 2025 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier:    MIT
import cutcells
import numpy as np
import pytest


def marching_squares_flag(ls_values: np.ndarray) -> int:
    flag = 0
    for i, v in enumerate(ls_values):
        if v < 0.0:
            flag |= 1 << i
    return flag


def asymptotic_decider(ls_values: np.ndarray) -> bool:
    f0, f1, f2, f3 = ls_values
    return f0 * f3 - f1 * f2 >= 0.0


def test_opposite_corner_masks_are_ambiguous():
    ls_a = np.array([-1.0, 1.0, 1.0, -1.0])
    ls_b = np.array([1.0, -1.0, -1.0, 1.0])

    assert marching_squares_flag(ls_a) == 0b1001
    assert marching_squares_flag(ls_b) == 0b0110


def test_asymptotic_decider_scale_invariant_and_different_diagonals():
    base = np.array([-1.0, 0.5, 1.2, -0.25])
    scaled = 7.5 * base

    decision_base = asymptotic_decider(base)
    decision_scaled = asymptotic_decider(scaled)

    assert decision_base == decision_scaled

    # Flip sign on one corner to flip the diagonal choice
    flipped = base.copy()
    flipped[2] *= -1.0
    decision_flipped = asymptotic_decider(flipped)

    assert decision_base != decision_flipped


def _segments_from_cutcell(cut_cell):
    coords = np.asarray(cut_cell.vertex_coords)
    conn = np.asarray(cut_cell.connectivity)
    offsets = np.asarray(cut_cell.offsets)
    segments = []
    for i in range(len(offsets) - 1):
        begin = int(offsets[i])
        end = int(offsets[i + 1])
        n = end - begin
        assert n == 2
        a = int(conn[begin])
        b = int(conn[begin + 1])
        segments.append((coords[a], coords[b]))
    return segments


def _edge_id_for_point(p, tol=1e-12):
    x, y = float(p[0]), float(p[1])
    if abs(y - 0.0) < tol:
        return 0  # edge 0: v0-v1
    if abs(x - 0.0) < tol:
        return 1  # edge 1: v0-v2
    if abs(x - 1.0) < tol:
        return 2  # edge 2: v1-v3
    if abs(y - 1.0) < tol:
        return 3  # edge 3: v2-v3
    raise AssertionError(f"Point not on boundary: {p}")


@pytest.mark.parametrize(
    "ls_values, expected_pairs",
    [
        # d = f0*f3 - f1*f2 > 0 -> connect diagonal (0,3).
        (np.array([-1.0, 0.1, 0.1, -1.0]), {frozenset({0, 2}), frozenset({1, 3})}),
        # d < 0 -> connect diagonal (1,2).
        (np.array([-0.1, 1.0, 1.0, -0.1]), {frozenset({0, 1}), frozenset({2, 3})}),
    ],
)
def test_quad_interface_ambiguous_pairs_follow_decider(ls_values, expected_pairs):
    vertex_coordinates = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0])
    cell_type = cutcells.CellType.quadrilateral
    gdim = 2

    cut_cell = cutcells.cut(
        cell_type, vertex_coordinates, gdim, ls_values, "phi=0", True
    )
    segments = _segments_from_cutcell(cut_cell)
    assert len(segments) == 2

    pairs = set()
    for a, b in segments:
        ea = _edge_id_for_point(a)
        eb = _edge_id_for_point(b)
        pairs.add(frozenset({ea, eb}))

    assert pairs == expected_pairs


@pytest.mark.parametrize(
    "ls_values, expected_mask",
    [
        # mask 1001: connected variant is d>=0 (variant 0)
        (np.array([-2.0, 1.0, 1.0, -2.0]), 0b1001),
        # mask 0110: connected variant is d<0 (variant 1)
        (np.array([1.0, -2.0, -2.0, 1.0]), 0b0110),
    ],
)
def test_quad_ambiguous_connected_variant_keeps_quads(ls_values, expected_mask):
    vertex_coordinates = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0])
    cell_type = cutcells.CellType.quadrilateral
    gdim = 2

    assert marching_squares_flag(ls_values) == expected_mask

    inside = cutcells.cut(
        cell_type, vertex_coordinates, gdim, ls_values, "phi<0", False
    )
    types = list(inside.types)

    assert types == [
        int(cutcells.CellType.quadrilateral.value),
        int(cutcells.CellType.quadrilateral.value),
    ]

    conn = np.asarray(inside.connectivity)
    offsets = np.asarray(inside.offsets)
    assert len(offsets) == 3
    for i in range(2):
        assert int(offsets[i + 1] - offsets[i]) == 4
    assert int(offsets[-1]) == len(conn)
