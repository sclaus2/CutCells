# Copyright (c) 2024 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier:    MIT
"""Tests for mapping.h / compute_physical_cut_vertices & complete_from_physical.

Design:
  - CutCells normally cuts in PHYSICAL space:
      _vertex_coords from cut() = physical coordinates.
    To obtain the reference coordinates: call complete_from_physical().
      -> copies _vertex_coords -> _vertex_coords_phys (physical)
      -> pulls back to fill _vertex_coords with reference coords

  - complete_from_physical (the standard workflow):
      After: _vertex_coords = reference,  _vertex_coords_phys = physical.

  - compute_physical_cut_vertices (reference-space cut workflow):
      If _vertex_coords were already in reference coords, push them to physical.
      After: _vertex_coords = unchanged (reference), _vertex_coords_phys = physical.

Coverage:
  - complete_from_physical: pull-back produces correct reference coordinates.
  - compute_physical_cut_vertices: push-forward from reference to physical.
  - round-trip consistency.
"""

import numpy as np
import pytest
import cutcells


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _unit_triangle():
    """Unit right triangle (0,0),(1,0),(0,1) in physical/reference space."""
    vertex_coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
    return vertex_coords, 2


def _scaled_triangle(sx, sy):
    """Triangle (0,0),(sx,0),(0,sy) — physical coords differ from reference."""
    vertex_coords = np.array([0.0, 0.0, sx, 0.0, 0.0, sy], dtype=np.float64)
    return vertex_coords, 2


def _standard_tet():
    """Standard tetrahedron: (0,0,0),(1,0,0),(0,1,0),(0,0,1)."""
    vertex_coords = np.array(
        [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], dtype=np.float64
    )
    return vertex_coords, 3


# ---------------------------------------------------------------------------
# 1.  complete_from_physical  (standard workflow: cut in physical space)
# ---------------------------------------------------------------------------


class TestCompleteFromPhysical:
    """cut() is performed with physical vertex coordinates.
    complete_from_physical() sets _vertex_coords_phys = copy of _vertex_coords
    (i.e. the physical cut coords), then pulls back to fill _vertex_coords with
    reference coords."""

    def test_unit_triangle_identity(self):
        """For the unit triangle (identity map) ref == phys after call."""
        vertex_coords, gdim = _unit_triangle()
        ls = np.array([0.5, -0.5, 0.5], dtype=np.float64)

        cut_cell = cutcells.cut(
            cutcells.CellType.triangle, vertex_coords, gdim, ls, "phi<0", True
        )
        assert cut_cell.vertex_coords_phys.size == 0, (
            "vertex_coords_phys must be empty before complete_from_physical"
        )
        phys_from_cutter = cut_cell.vertex_coords.copy()

        cutcells.complete_from_physical(cut_cell)

        vc_ref = cut_cell.vertex_coords
        vc_phys = cut_cell.vertex_coords_phys

        # physical should be what came out of the cutter
        np.testing.assert_allclose(
            vc_phys, phys_from_cutter.reshape(-1, gdim), atol=1e-14
        )
        # reference == physical for identity map
        np.testing.assert_allclose(
            vc_ref, phys_from_cutter.reshape(-1, gdim), atol=1e-14
        )

    def test_scaled_triangle_pull_back(self):
        """Scaled triangle (3×2): pull-back yields ref = phys / scale."""
        sx, sy = 3.0, 2.0
        vertex_coords, gdim = _scaled_triangle(sx, sy)
        ls = np.array([0.1, -0.1, 0.2], dtype=np.float64)

        cut_cell = cutcells.cut(
            cutcells.CellType.triangle, vertex_coords, gdim, ls, "phi<0", True
        )
        phys_from_cutter = cut_cell.vertex_coords.copy()

        cutcells.complete_from_physical(cut_cell)

        vc_ref = cut_cell.vertex_coords
        vc_phys = cut_cell.vertex_coords_phys

        np.testing.assert_allclose(
            vc_phys, phys_from_cutter.reshape(-1, gdim), atol=1e-13
        )
        # For affine map v0=0, J=diag(sx,sy):  ref = phys / [sx, sy]
        expected_ref = phys_from_cutter.reshape(-1, gdim) / np.array([sx, sy])
        np.testing.assert_allclose(vc_ref, expected_ref, atol=1e-13)

    def test_tetrahedron_pull_back(self):
        """Standard tet (identity map): ref == phys after complete_from_physical."""
        vertex_coords, gdim = _standard_tet()
        ls = np.array([0.1, -0.1, 0.2, 0.2], dtype=np.float64)

        cut_cell = cutcells.cut(
            cutcells.CellType.tetrahedron, vertex_coords, gdim, ls, "phi<0", True
        )
        phys_from_cutter = cut_cell.vertex_coords.copy()

        cutcells.complete_from_physical(cut_cell)

        np.testing.assert_allclose(
            cut_cell.vertex_coords_phys, phys_from_cutter.reshape(-1, gdim), atol=1e-14
        )
        np.testing.assert_allclose(
            cut_cell.vertex_coords, phys_from_cutter.reshape(-1, gdim), atol=1e-14
        )

    def test_phys_preserved_after_first_call(self):
        """complete_from_physical preserves the original physical cut coords in
        vertex_coords_phys and converts vertex_coords to reference coordinates.
        Calling it a second time would further pull back (not idempotent) — verify
        that only the first call gives meaningful physical-space vertex_coords_phys."""
        vertex_coords, gdim = _scaled_triangle(2.5, 1.5)
        ls = np.array([0.5, -0.5, 0.5], dtype=np.float64)

        cut_cell = cutcells.cut(
            cutcells.CellType.triangle, vertex_coords, gdim, ls, "phi<0", True
        )
        phys_original = cut_cell.vertex_coords.copy()  # physical from cut()

        cutcells.complete_from_physical(cut_cell)

        # After first call: _vertex_coords_phys == original physical cut coords
        np.testing.assert_allclose(
            cut_cell.vertex_coords_phys,
            phys_original.reshape(-1, gdim),
            atol=1e-13,
            err_msg="vertex_coords_phys must equal original physical cut coords",
        )
        # _vertex_coords must have been pulled back
        expected_ref = phys_original.reshape(-1, gdim) / np.array([2.5, 1.5])
        np.testing.assert_allclose(cut_cell.vertex_coords, expected_ref, atol=1e-13)


# ---------------------------------------------------------------------------
# 2.  compute_physical_cut_vertices  (reference-frame cut workflow)
# ---------------------------------------------------------------------------


class TestComputePhysicalCutVertices:
    """compute_physical_cut_vertices pushes _vertex_coords (treated as reference)
    to physical space using the affine map defined by _parent_vertex_coords.

    For the standard physical-space cutting workflow, _vertex_coords are already
    physical, so compute_physical_cut_vertices would push them *again* (incorrect
    usage). The correct function in that case is complete_from_physical.

    Here we test the correct usage: cutting with reference (unit) vertex coordinates
    and then pushing forward to physical."""

    def test_unit_triangle_identity(self):
        """For identity map push-forward gives phys == ref."""
        vertex_coords, gdim = _unit_triangle()
        ls = np.array([0.5, -0.5, 0.5], dtype=np.float64)

        cut_cell = cutcells.cut(
            cutcells.CellType.triangle, vertex_coords, gdim, ls, "phi<0", True
        )
        assert cut_cell.vertex_coords_phys.size == 0

        cutcells.compute_physical_cut_vertices(cut_cell)

        vc_ref = cut_cell.vertex_coords.reshape(-1, gdim)
        vc_phys = cut_cell.vertex_coords_phys
        np.testing.assert_allclose(vc_phys, vc_ref, atol=1e-14)

    def test_tetrahedron_identity(self):
        """Standard tet (identity map): phys == ref."""
        vertex_coords, gdim = _standard_tet()
        ls = np.array([0.1, -0.1, 0.2, 0.2], dtype=np.float64)

        cut_cell = cutcells.cut(
            cutcells.CellType.tetrahedron, vertex_coords, gdim, ls, "phi<0", True
        )
        cutcells.compute_physical_cut_vertices(cut_cell)

        vc_ref = cut_cell.vertex_coords.reshape(-1, gdim)
        vc_phys = cut_cell.vertex_coords_phys
        np.testing.assert_allclose(vc_phys, vc_ref, atol=1e-14)

    def test_idempotent(self):
        """Calling compute_physical_cut_vertices twice gives the same physical coords."""
        vertex_coords, gdim = _unit_triangle()
        ls = np.array([0.5, -0.5, 0.5], dtype=np.float64)

        cut_cell = cutcells.cut(
            cutcells.CellType.triangle, vertex_coords, gdim, ls, "phi<0", True
        )
        cutcells.compute_physical_cut_vertices(cut_cell)
        phys1 = cut_cell.vertex_coords_phys.copy()
        cutcells.compute_physical_cut_vertices(cut_cell)
        phys2 = cut_cell.vertex_coords_phys
        np.testing.assert_array_equal(phys1, phys2)
