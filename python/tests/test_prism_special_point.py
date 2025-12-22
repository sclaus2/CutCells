import numpy as np
import cutcells


def test_prism_interface_uses_special_point_n0_case10():
    # Pick a concrete case that (per generated prism interface tables) uses token 200 (N0).
    # Wedge/prism has 6 vertices => 64 cases; case id 10 corresponds to mask bits set at vertices {1,3}.
    ls_values = np.array([1.0, -1.0, 1.0, -1.0, 1.0, 1.0], dtype=np.float64)

    # VTK_WEDGE / prism vertex ordering:
    # bottom tri: 0:(0,0,0) 1:(1,0,0) 2:(0,1,0)
    # top tri:    3:(0,0,1) 4:(1,0,1) 5:(0,1,1)
    vertex_coordinates = np.array(
        [
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            1.0,
            0.0,
            1.0,
            0.0,
            1.0,
            1.0,
        ],
        dtype=np.float64,
    )

    cut_cell = cutcells.cut(
        cutcells.CellType.prism,
        vertex_coordinates,
        3,
        ls_values,
        "phi=0",
        False,
    )

    coords = np.asarray(cut_cell.vertex_coords, dtype=np.float64).reshape(-1, 3)

    tol = 1e-10
    assert coords[:, 0].min() >= -tol and coords[:, 0].max() <= 1.0 + tol
    assert coords[:, 1].min() >= -tol and coords[:, 1].max() <= 1.0 + tol
    assert coords[:, 2].min() >= -tol and coords[:, 2].max() <= 1.0 + tol

    # Edge intersections live on the prism boundary. N0 should be strictly interior for this case.
    inside = (
        (coords[:, 0] > tol)
        & (coords[:, 1] > tol)
        & (coords[:, 0] + coords[:, 1] < 1.0 - tol)
        & (coords[:, 2] > tol)
        & (coords[:, 2] < 1.0 - tol)
    )
    assert bool(np.any(inside))
