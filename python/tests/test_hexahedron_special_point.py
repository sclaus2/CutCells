import numpy as np
import cutcells


def _is_on_unit_cube_edge(p: np.ndarray, tol: float = 1e-12) -> bool:
    # On a unit-cube edge => at least two coordinates are on {0,1} (within tol)
    on01 = [abs(p[i] - 0.0) < tol or abs(p[i] - 1.0) < tol for i in range(3)]
    return sum(on01) >= 2


def test_hexahedron_interface_uses_special_point_n0():
    # Pick a concrete VTK case that is known (from generated tables) to reference token 200 (N0).
    # Case id 67 corresponds to mask bits set at vertices {0,1,6}.
    ls_values = np.array([-1.0, -1.1, 1.2, 1.3, 1.4, 1.5, -1.6, 1.7], dtype=np.float64)

    # VTK hexahedron vertex ordering:
    # 0:(0,0,0) 1:(1,0,0) 2:(1,1,0) 3:(0,1,0) 4:(0,0,1) 5:(1,0,1) 6:(1,1,1) 7:(0,1,1)
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
            0.0,
            0.0,
            1.0,
            1.0,
            0.0,
            1.0,
            1.0,
            1.0,
            1.0,
            0.0,
            1.0,
            1.0,
        ],
        dtype=np.float64,
    )

    cut_cell = cutcells.cut(
        cutcells.CellType.hexahedron,
        vertex_coordinates,
        3,
        ls_values,
        "phi=0",
        False,
    )

    coords = np.asarray(cut_cell.vertex_coords)

    tol = 1e-10
    assert coords[:, 0].min() >= -tol and coords[:, 0].max() <= 1.0 + tol
    assert coords[:, 1].min() >= -tol and coords[:, 1].max() <= 1.0 + tol
    assert coords[:, 2].min() >= -tol and coords[:, 2].max() <= 1.0 + tol

    # All edge-intersection points lie on unit-cube edges; the special point should not.
    assert any(not _is_on_unit_cube_edge(p) for p in coords)
