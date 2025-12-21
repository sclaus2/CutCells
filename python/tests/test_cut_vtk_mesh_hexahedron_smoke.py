# Copyright (c) 2025 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier:    MIT

import numpy as np

import cutcells

VTK_TETRA = 10
VTK_HEXAHEDRON = 12


def _single_unit_hex_vtk():
    """One VTK_HEXAHEDRON cell on the unit cube in VTK vertex order."""
    points = np.array(
        [
            [0.0, 0.0, 0.0],  # 0
            [1.0, 0.0, 0.0],  # 1
            [1.0, 1.0, 0.0],  # 2
            [0.0, 1.0, 0.0],  # 3
            [0.0, 0.0, 1.0],  # 4
            [1.0, 0.0, 1.0],  # 5
            [1.0, 1.0, 1.0],  # 6
            [0.0, 1.0, 1.0],  # 7
        ],
        dtype=float,
    )

    connectivity = np.array([0, 1, 2, 3, 4, 5, 6, 7], dtype=np.int32)
    offset = np.array([0], dtype=np.int32)
    vtk_type = np.array([VTK_HEXAHEDRON], dtype=np.int32)

    return points, connectivity, offset, vtk_type


def test_cut_vtk_mesh_single_hexahedron_smoke():
    points, connectivity, offset, vtk_type = _single_unit_hex_vtk()

    # Intersecting cut: split by plane x=0.3
    ls_vals = points[:, 0] - 0.3

    cut_mesh = cutcells.cut_vtk_mesh(
        ls_vals,
        points,
        connectivity,
        offset,
        vtk_type,
        "phi<0",
        triangulate=True,
    )

    assert len(cut_mesh.types) > 0
    # Fallback hex implementation cuts via tetrahedra; with triangulate=True it should be tetra-dominant.
    assert VTK_TETRA in list(cut_mesh.types)
