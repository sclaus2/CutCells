# Copyright (c) 2025 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier:    MIT

import numpy as np

import cutcells

VTK_TRIANGLE = 5
VTK_QUAD = 9


def _quad_mesh_2x2_single_cell():
    """One VTK_QUAD cell on [0,1]x[0,1] in z=0 plane."""
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        dtype=float,
    )

    connectivity = np.array([0, 1, 2, 3], dtype=np.int32)
    offset = np.array([0], dtype=np.int32)
    vtk_type = np.array([VTK_QUAD], dtype=np.int32)

    return points, connectivity, offset, vtk_type


def test_cut_vtk_mesh_triangulate_flag_controls_output_types():
    points, connectivity, offset, vtk_type = _quad_mesh_2x2_single_cell()

    # Level set: phi = x - 0.3 produces a vertical cut line inside the quad.
    # This is a 2-inside/2-outside case that can naturally produce a quad when not triangulated.
    ls_vals = points[:, 0] - 0.3

    cut_mesh_quads = cutcells.cut_vtk_mesh(
        ls_vals, points, connectivity, offset, vtk_type, "phi<0", triangulate=False
    )
    assert VTK_QUAD in list(cut_mesh_quads.types)

    cut_mesh_tris = cutcells.cut_vtk_mesh(
        ls_vals, points, connectivity, offset, vtk_type, "phi<0", triangulate=True
    )
    assert VTK_QUAD not in list(cut_mesh_tris.types)
    assert VTK_TRIANGLE in list(cut_mesh_tris.types)
