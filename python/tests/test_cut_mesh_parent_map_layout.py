import numpy as np

import cutcells

VTK_QUAD = 9


def _two_quad_strip_vtk():
    points = np.array(
        [
            [0.0, 0.0, 0.0],  # 0
            [1.0, 0.0, 0.0],  # 1
            [2.0, 0.0, 0.0],  # 2
            [0.0, 1.0, 0.0],  # 3
            [1.0, 1.0, 0.0],  # 4
            [2.0, 1.0, 0.0],  # 5
        ],
        dtype=float,
    )

    connectivity = np.array(
        [
            0,
            1,
            4,
            3,
            1,
            2,
            5,
            4,
        ],
        dtype=np.int32,
    )
    offset = np.array([0, 4], dtype=np.int32)
    vtk_type = np.array([VTK_QUAD, VTK_QUAD], dtype=np.int32)

    return points, connectivity, offset, vtk_type


def test_cut_mesh_parent_map_row_layout_smoke():
    points, connectivity, offset, vtk_type = _two_quad_strip_vtk()

    ls_vals = points[:, 1] - 0.5

    cut_mesh = cutcells.cut_vtk_mesh(
        ls_vals,
        points.ravel(),
        connectivity,
        offset,
        vtk_type,
        "phi<0",
        triangulate=False,
    )

    types = np.asarray(cut_mesh.types)
    parent_map = np.asarray(cut_mesh.parent_map)

    assert types.size > 0
    assert parent_map.size == types.size
    assert set(parent_map.tolist()).issubset({0, 1})
    assert 0 in parent_map and 1 in parent_map
