# Copyright (c) 2024 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier:    MIT
"""Tests for quadrature.h / make_quadrature.

Coverage:
  - QuadratureRules struct (tdim, points, weights, offset, parent_map).
  - make_quadrature: full-cell case integrates to correct area/volume.
  - make_quadrature: cut-cell case integrates subdomain area.
  - Batch interface: multiple cut cells, offset/parent_map metadata.
  - Cell types: triangle, tetrahedron, hexahedron.
"""

import numpy as np
import pytest
import cutcells


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_and_enrich(
    cell_type, vertex_coords, gdim, ls_values, cut_type, triangulate=True
):
    """Cut (physical coords), complete frame (complete_from_physical), return enriched CutCell."""
    cut_cell = cutcells.cut(
        cell_type, vertex_coords, gdim, ls_values, cut_type, triangulate
    )
    cutcells.complete_from_physical(cut_cell)
    return cut_cell


def _integrate(cut_cells_list, order=3):
    """Return QuadratureRules and summed weight for cut_cells_list."""
    qr = cutcells.make_quadrature(cut_cells_list, order)
    return qr, np.sum(qr.weights)


# ---------------------------------------------------------------------------
# 1. QuadratureRules struct
# ---------------------------------------------------------------------------


class TestQuadratureRulesStruct:
    def test_default_empty(self):
        qr = cutcells._cutcellscpp.QuadratureRules_float64()
        assert qr.tdim == 0
        assert qr.points.size == 0
        assert qr.weights.size == 0
        assert qr.offset.size == 0
        assert qr.parent_map.size == 0

    def test_types(self):
        """After make_quadrature the arrays have expected dtypes."""
        vertex_coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
        ls = np.array([0.1, -0.1, 0.2], dtype=np.float64)
        cut_cell = _make_and_enrich(
            cutcells.CellType.triangle, vertex_coords, 2, ls, "phi<0"
        )
        qr = cutcells.make_quadrature([cut_cell], 2)
        assert qr.points.dtype == np.float64
        assert qr.weights.dtype == np.float64
        assert qr.offset.dtype == np.int32
        assert qr.parent_map.dtype == np.int32


# ---------------------------------------------------------------------------
# 2. Full-cell integration (no cut)
# ---------------------------------------------------------------------------


class TestFullCellIntegration:
    """Barely-intersected cells (epsilon level-set) must integrate to ~full area/volume.
    A level-set like [-1,-1,...,+eps] creates a tiny cut but leaves essentially the full
    interior intact, so the integral must be close to the full cell measure."""

    def test_triangle_full_cell_area(self):
        """Unit right triangle has area 0.5."""
        vertex_coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
        # Tiny positive at one vertex: interior ~= full triangle
        eps = 1e-4
        ls = np.array([-1.0, -1.0, eps], dtype=np.float64)
        cut_cell = _make_and_enrich(
            cutcells.CellType.triangle, vertex_coords, 2, ls, "phi<0"
        )
        _, total = _integrate([cut_cell], order=4)
        np.testing.assert_allclose(total, 0.5, rtol=1e-3)

    def test_triangle_scaled_area(self):
        """Scaled triangle (3×2) has area 3.0."""
        vertex_coords = np.array([0.0, 0.0, 3.0, 0.0, 0.0, 2.0], dtype=np.float64)
        eps = 1e-4
        ls = np.array([-1.0, -1.0, eps], dtype=np.float64)
        cut_cell = _make_and_enrich(
            cutcells.CellType.triangle, vertex_coords, 2, ls, "phi<0"
        )
        _, total = _integrate([cut_cell], order=4)
        np.testing.assert_allclose(total, 3.0, rtol=1e-3)

    def test_tetrahedron_full_cell_volume(self):
        """Standard tet has volume 1/6."""
        vertex_coords = np.array(
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            dtype=np.float64,
        )
        eps = 1e-4
        ls = np.array([-1.0, -1.0, -1.0, eps], dtype=np.float64)
        cut_cell = _make_and_enrich(
            cutcells.CellType.tetrahedron, vertex_coords, 3, ls, "phi<0"
        )
        _, total = _integrate([cut_cell], order=3)
        np.testing.assert_allclose(total, 1.0 / 6.0, rtol=1e-3)

    def test_hexahedron_full_cell_volume(self):
        """Unit cube has volume 1.0.

        Basix hexahedron vertex ordering:
          v0=(0,0,0)  v1=(1,0,0)  v2=(0,1,0)  v3=(1,1,0)
          v4=(0,0,1)  v5=(1,0,1)  v6=(0,1,1)  v7=(1,1,1)
        """
        vertex_coords = np.array(
            [
                0.0,
                0.0,
                0.0,  # v0
                1.0,
                0.0,
                0.0,  # v1
                0.0,
                1.0,
                0.0,  # v2
                1.0,
                1.0,
                0.0,  # v3
                0.0,
                0.0,
                1.0,  # v4
                1.0,
                0.0,
                1.0,  # v5
                0.0,
                1.0,
                1.0,  # v6
                1.0,
                1.0,
                1.0,  # v7
            ],
            dtype=np.float64,
        )
        eps = 1e-4
        # v7=(1,1,1) has eps>0; all others negative -> phi<0 is almost the full cube
        ls = np.array([-1.0] * 7 + [eps], dtype=np.float64)
        cut_cell = _make_and_enrich(
            cutcells.CellType.hexahedron, vertex_coords, 3, ls, "phi<0"
        )
        _, total = _integrate([cut_cell], order=3)
        np.testing.assert_allclose(total, 1.0, rtol=1e-3)


# ---------------------------------------------------------------------------
# 3. Cut-cell integration
# ---------------------------------------------------------------------------


class TestCutCellIntegration:
    """The interior cut of a known fraction must integrate to that fraction."""

    def test_triangle_half_cut(self):
        """Level set x - 0.5 = 0 cuts the unit triangle.
        The interior phi<0 part has area 0.5 * 0.5 * (1 - 0.5^2) ... compute
        analytically: triangle with vertices (0,0),(0.5,0),(0,1),(... tricky).
        Use the known vol from test_triangle.py instead: ls=[0.1,-0.1,0.2] -> 1/12."""
        vertex_coords = np.array([0.0, 0.0, 1.0, 0.0, 1.0, 1.0], dtype=np.float64)
        ls = np.array([0.1, -0.1, 0.2], dtype=np.float64)
        cut_cell = _make_and_enrich(
            cutcells.CellType.triangle, vertex_coords, 2, ls, "phi<0"
        )
        _, total = _integrate([cut_cell], order=5)
        np.testing.assert_allclose(total, 1.0 / 12.0, rtol=1e-9)

    def test_tetrahedron_cut(self):
        """Known tet cut: ls=[0.1,-0.1,0.2,0.2] interior vol = 4/27."""
        vertex_coords = np.array(
            [1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0],
            dtype=np.float64,
        )
        ls = np.array([0.1, -0.1, 0.2, 0.2], dtype=np.float64)
        cut_cell = _make_and_enrich(
            cutcells.CellType.tetrahedron, vertex_coords, 3, ls, "phi<0"
        )
        _, total = _integrate([cut_cell], order=5)
        np.testing.assert_allclose(total, 4.0 / 27.0, rtol=1e-8)

    @pytest.mark.parametrize("order", [1, 2, 3, 4, 5])
    def test_triangle_order_convergence(self, order):
        """All orders should give the exact answer for a linear level-set cut
        because the sub-cells are simplices (exact for polynomials of any order)."""
        vertex_coords = np.array([0.0, 0.0, 1.0, 0.0, 1.0, 1.0], dtype=np.float64)
        ls = np.array([0.1, -0.1, 0.2], dtype=np.float64)
        cut_cell = _make_and_enrich(
            cutcells.CellType.triangle, vertex_coords, 2, ls, "phi<0"
        )
        _, total = _integrate([cut_cell], order=order)
        np.testing.assert_allclose(
            total, 1.0 / 12.0, rtol=1e-9, err_msg=f"order={order} failed"
        )


# ---------------------------------------------------------------------------
# 4. Batch interface
# ---------------------------------------------------------------------------


class TestBatchIntegration:
    """make_quadrature([c1, c2, ...]) must produce correct offset/parent_map."""

    def _make_triangle_cells(self, n_cells):
        """n_cells copies of a half-cut unit triangle."""
        vertex_coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
        ls = np.array([0.5, -0.5, 0.5], dtype=np.float64)
        cells = []
        for _ in range(n_cells):
            c = _make_and_enrich(
                cutcells.CellType.triangle, vertex_coords, 2, ls, "phi<0"
            )
            cells.append(c)
        return cells

    def test_offset_length(self):
        cells = self._make_triangle_cells(4)
        qr = cutcells.make_quadrature(cells, 3)
        assert len(qr.offset) == len(cells) + 1, "offset length must be num_rules + 1"
        assert qr.offset[0] == 0, "offset[0] must be 0"
        assert qr.offset[-1] == len(qr.weights), (
            "last offset must equal total number of quadrature points"
        )

    def test_parent_map_length(self):
        cells = self._make_triangle_cells(3)
        qr = cutcells.make_quadrature(cells, 3)
        assert len(qr.parent_map) == len(cells)

    def test_parent_map_values(self):
        cells = self._make_triangle_cells(3)
        qr = cutcells.make_quadrature(cells, 3)
        np.testing.assert_array_equal(qr.parent_map, [0, 1, 2])

    def test_batch_sum_equals_individual_sum(self):
        cells = self._make_triangle_cells(5)
        qr_batch = cutcells.make_quadrature(cells, 3)
        batch_total = np.sum(qr_batch.weights)

        individual_total = 0.0
        for c in cells:
            qr_i = cutcells.make_quadrature([c], 3)
            individual_total += np.sum(qr_i.weights)

        np.testing.assert_allclose(batch_total, individual_total, rtol=1e-14)

    def test_points_shape(self):
        cells = self._make_triangle_cells(2)
        qr = cutcells.make_quadrature(cells, 3)
        total_pts = qr.offset[-1]
        assert qr.points.size == total_pts * qr.tdim, (
            "flat points array size must be total_pts * tdim"
        )

    def test_tdim_triangle(self):
        cells = self._make_triangle_cells(1)
        qr = cutcells.make_quadrature(cells, 2)
        assert qr.tdim == 2, "tdim for triangle parent should be 2"

    def test_tdim_tetrahedron(self):
        vertex_coords = np.array(
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            dtype=np.float64,
        )
        ls = np.array([0.1, -0.1, 0.2, 0.2], dtype=np.float64)
        cut_cell = _make_and_enrich(
            cutcells.CellType.tetrahedron, vertex_coords, 3, ls, "phi<0"
        )
        qr = cutcells.make_quadrature([cut_cell], 2)
        assert qr.tdim == 3


# ---------------------------------------------------------------------------
# 5. Float32 variant
# ---------------------------------------------------------------------------


class TestFloat32:
    """Smoke-test the float32 code path."""

    def test_triangle_full_cell_float32(self):
        vertex_coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float32)
        ls = np.array([0.1, -0.1, 0.2], dtype=np.float32)
        cut_cell = cutcells.cut(
            cutcells.CellType.triangle, vertex_coords, 2, ls, "phi<0", True
        )
        cutcells.complete_from_physical(cut_cell)
        qr = cutcells.make_quadrature([cut_cell], 3)
        total = np.sum(qr.weights)
        np.testing.assert_allclose(total, 1.0 / 12.0, rtol=1e-4)
        assert qr.points.dtype == np.float32
        assert qr.weights.dtype == np.float32
