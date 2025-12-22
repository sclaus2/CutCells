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

    coords = np.asarray(cut_mesh.vertex_coords, dtype=float).reshape(-1, 3)
    # All generated vertices must stay within the parent unit cube (tolerant for FP noise).
    tol = 1e-10
    assert coords[:, 0].min() >= -tol and coords[:, 0].max() <= 1.0 + tol
    assert coords[:, 1].min() >= -tol and coords[:, 1].max() <= 1.0 + tol
    assert coords[:, 2].min() >= -tol and coords[:, 2].max() <= 1.0 + tol

    assert len(cut_mesh.types) > 0
    # Do not constrain emitted cell types here; the table-driven hex cutter may emit
    # hex/prism/pyramid/tet depending on the case.
