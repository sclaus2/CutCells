# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier:    MIT

import time

import numpy as np

import cutcells

VTK_HEXAHEDRON = 12


def _single_unit_hex_vtk():
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ],
        dtype=float,
    )

    connectivity = np.array([0, 1, 2, 3, 4, 5, 6, 7], dtype=np.int32)
    offset = np.array([0], dtype=np.int32)
    vtk_type = np.array([VTK_HEXAHEDRON], dtype=np.int32)

    return points, connectivity, offset, vtk_type


def test_cut_vtk_mesh_terminates_quickly_on_small_mesh():
    points, connectivity, offset, vtk_type = _single_unit_hex_vtk()
    ls_vals = points[:, 0] - 0.3

    iterations = 200
    start = time.perf_counter()

    for _ in range(iterations):
        cut_mesh = cutcells.cut_vtk_mesh(
            ls_vals,
            points.ravel(),
            connectivity,
            offset,
            vtk_type,
            "phi<0",
            triangulate=True,
        )

    elapsed = time.perf_counter() - start

    assert cut_mesh.types.size > 0
    assert np.asarray(cut_mesh.vertex_coords).shape[1] == 3
    assert elapsed < 5.0, (
        f"cut_vtk_mesh appears too slow/non-terminating for a trivial mesh: "
        f"{iterations} calls took {elapsed:.3f}s"
    )
