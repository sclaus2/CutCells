# Copyright (c) 2025 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier:    MIT
import numpy as np
import pytest

import cutcells

TRIANGLE = int(cutcells.CellType.triangle.value)
QUADRILATERAL = int(cutcells.CellType.quadrilateral.value)


def _square_vertices():
    return np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0], dtype=float)


def test_quadrilateral_opposite_corners_disconnected():
    vertices = _square_vertices()
    ls_values = np.array([-0.5, 0.5, 0.5, -0.5], dtype=float)

    cut_inside = cutcells.cut(
        cutcells.CellType.quadrilateral, vertices, 2, ls_values, "phi<0", False
    )
    assert list(cut_inside.types) == [TRIANGLE, TRIANGLE]
    assert np.isclose(cut_inside.volume(), 0.25)

    cut_outside = cutcells.cut(
        cutcells.CellType.quadrilateral, vertices, 2, ls_values, "phi>0", False
    )
    assert list(cut_outside.types) == [QUADRILATERAL, QUADRILATERAL]
    assert np.isclose(cut_outside.volume(), 0.75)


def test_quadrilateral_pentagon_splits_triangle_quad():
    vertices = _square_vertices()
    ls_values = np.array([0.5, -0.5, -0.5, -0.5], dtype=float)

    cut_inside = cutcells.cut(
        cutcells.CellType.quadrilateral, vertices, 2, ls_values, "phi<0", False
    )
    assert list(cut_inside.types) == [TRIANGLE, QUADRILATERAL]
    assert np.isclose(cut_inside.volume(), 0.875)

    cut_outside = cutcells.cut(
        cutcells.CellType.quadrilateral, vertices, 2, ls_values, "phi>0", False
    )
    assert list(cut_outside.types) == [TRIANGLE]
    assert np.isclose(cut_outside.volume(), 0.125)
