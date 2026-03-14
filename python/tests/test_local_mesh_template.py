import numpy as np

import cutcells


def test_init_local_mesh_from_template_triangle_p2():
    tpl = cutcells.iso_p1_template(cutcells.CellType.triangle, 2)

    # Basix P2 triangle node coordinates (flat, gdim=2)
    parent_cell_coords = np.array(
        [
            0.0, 0.0,
            1.0, 0.0,
            0.0, 1.0,
            0.5, 0.5,
            0.0, 0.5,
            0.5, 0.0,
        ],
        dtype=np.float64,
    )

    lm = cutcells.init_local_mesh_from_template(
        tpl,
        parent_cell_coords,
        cutcells.CellType.triangle,
        0,
        1,
    )

    assert lm.gdim == 2
    assert lm.tdim == 2
    assert lm.n_vertices() == 6
    assert lm.n_cells() == 4
    assert lm.n_edges() == 9

    np.testing.assert_allclose(np.asarray(lm.vertex_x), parent_cell_coords)
