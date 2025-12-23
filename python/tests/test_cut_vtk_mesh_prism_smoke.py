import numpy as np

import cutcells

VTK_WEDGE = 13


def _single_unit_prism_vtk():
    """One VTK_WEDGE cell on a unit right triangular prism in VTK vertex order."""
    points = np.array(
        [
            [0.0, 0.0, 0.0],  # 0
            [1.0, 0.0, 0.0],  # 1
            [0.0, 1.0, 0.0],  # 2
            [0.0, 0.0, 1.0],  # 3
            [1.0, 0.0, 1.0],  # 4
            [0.0, 1.0, 1.0],  # 5
        ],
        dtype=float,
    )

    connectivity = np.array([0, 1, 2, 3, 4, 5], dtype=np.int32)
    offset = np.array([0], dtype=np.int32)
    vtk_type = np.array([VTK_WEDGE], dtype=np.int32)

    return points, connectivity, offset, vtk_type


def test_cut_vtk_mesh_single_prism_smoke():
    points, connectivity, offset, vtk_type = _single_unit_prism_vtk()

    # Intersecting cut: split by plane z=0.3
    ls_vals = points[:, 2] - 0.3

    cut_mesh = cutcells.cut_vtk_mesh(
        ls_vals,
        points,
        connectivity,
        offset,
        vtk_type,
        "phi<0",
        triangulate=True,
    )

    coords = np.asarray(cut_mesh.vertex_coords)
    tol = 1e-10
    assert coords[:, 0].min() >= -tol and coords[:, 0].max() <= 1.0 + tol
    assert coords[:, 1].min() >= -tol and coords[:, 1].max() <= 1.0 + tol
    assert coords[:, 2].min() >= -tol and coords[:, 2].max() <= 1.0 + tol

    assert len(cut_mesh.types) > 0
