import numpy as np
import cutcells


def test_pyramid_uses_special_point_n0_case19():
    # Pyramid has 5 vertices => 32 cases; case id 19 corresponds to mask bits set at vertices {0,1,4}.
    # (CutCells convention: bit i set when phi_i < 0.)
    ls_values = np.array([-1.0, -1.0, 1.0, 1.0, -1.0], dtype=np.float64)

    # VTK_PYRAMID vertex ordering:
    # base quad: 0:(0,0,0) 1:(1,0,0) 2:(1,1,0) 3:(0,1,0)
    # apex:      4:(0.5,0.5,1)
    vertex_coordinates = np.array(
        [
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            1.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.5,
            0.5,
            1.0,
        ],
        dtype=np.float64,
    )

    cut_cell = cutcells.cut(
        cutcells.CellType.pyramid,
        vertex_coordinates,
        3,
        ls_values,
        "phi<0",
        False,
    )

    coords = np.asarray(cut_cell.vertex_coords)

    # Basic bounds
    tol = 1e-10
    assert coords[:, 0].min() >= -tol and coords[:, 0].max() <= 1.0 + tol
    assert coords[:, 1].min() >= -tol and coords[:, 1].max() <= 1.0 + tol
    assert coords[:, 2].min() >= -tol and coords[:, 2].max() <= 1.0 + tol

    # Edge intersections lie on edges/faces. N0 should be strictly interior for this case.
    # For this unit pyramid (apex centered), the interior condition at height z is:
    # |x-0.5| < 0.5*(1-z) and |y-0.5| < 0.5*(1-z).
    x = coords[:, 0]
    y = coords[:, 1]
    z = coords[:, 2]

    half = 0.5 * (1.0 - z)
    inside = (
        (z > tol)
        & (z < 1.0 - tol)
        & (np.abs(x - 0.5) < half - tol)
        & (np.abs(y - 0.5) < half - tol)
    )
    assert bool(np.any(inside))
