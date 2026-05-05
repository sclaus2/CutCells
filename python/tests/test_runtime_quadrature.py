"""Pytest suite for runtime_quadrature and physical_points."""

import numpy as np
import pytest
import cutcells

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _full_cell_rules(pts, conn_indices, vtk_type_id, triangulate, order=4):
    """Return QuadratureRules for a fully-inside cell."""
    pts_arr = np.array(pts, dtype=np.float64).flatten()
    conn = np.array(conn_indices, dtype=np.int32)
    offs = np.array([0, len(conn_indices)], dtype=np.int32)
    vtype = np.array([vtk_type_id], dtype=np.int32)
    ls = np.full(len(pts), -1.0, dtype=np.float64)  # all inside
    rules = cutcells.runtime_quadrature(
        ls, pts_arr, conn, offs, vtype, "phi<0", triangulate, order
    )
    phys = cutcells.physical_points(rules, pts_arr, conn, offs, vtype)
    return rules, phys.reshape(-1, 3)


# ---------------------------------------------------------------------------
# Test data: (name, pts, conn, vtk_id, expected_vol)
# ---------------------------------------------------------------------------

INTERVAL_PTS = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
TRI_PTS = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
TET_PTS = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
QUAD_PTS = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
HEX_PTS = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
]
PRISM_PTS = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [0.0, 1.0, 1.0],
]
PYR_PTS = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.5, 0.5, 1.0],
]

CELL_CASES = [
    # (name, pts, conn, vtk_id, expected_volume)
    ("interval", INTERVAL_PTS, [0, 1], 3, 2.0),
    ("triangle", TRI_PTS, [0, 1, 2], 5, 0.5),
    ("tetrahedron", TET_PTS, [0, 1, 2, 3], 10, 1.0 / 6.0),
    ("quad", QUAD_PTS, [0, 1, 2, 3], 9, 1.0),
    ("hexahedron", HEX_PTS, list(range(8)), 12, 1.0),
    ("prism", PRISM_PTS, list(range(6)), 13, 0.5),
    ("pyramid", PYR_PTS, list(range(5)), 14, 1.0 / 3.0),
]


# ---------------------------------------------------------------------------
# Parametrized volume tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name,pts,conn,vtk_id,expected_vol", CELL_CASES)
@pytest.mark.parametrize("triangulate", [False, True])
def test_full_cell_volume(name, pts, conn, vtk_id, expected_vol, triangulate):
    """Sum of quadrature weights equals the physical volume of the cell."""
    rules, _ = _full_cell_rules(pts, conn, vtk_id, triangulate)
    assert abs(rules.weights.sum() - expected_vol) < 1e-9, (
        f"{name} triangulate={triangulate}: "
        f"got {rules.weights.sum()}, expected {expected_vol}"
    )


@pytest.mark.parametrize("name,pts,conn,vtk_id,expected_vol", CELL_CASES)
@pytest.mark.parametrize("triangulate", [False, True])
def test_full_cell_physical_points_shape(
    name, pts, conn, vtk_id, expected_vol, triangulate
):
    """physical_points returns array of shape (nq, 3)."""
    rules, phys = _full_cell_rules(pts, conn, vtk_id, triangulate)
    assert phys.shape == (len(rules.weights), 3)


# ---------------------------------------------------------------------------
# Cut cell tests (tet only — requires tdim == gdim == 3)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "cut_type,frac",
    [
        ("phi<0", 1.0 / 6.0),  # cut from one corner; vol < full tet
        ("phi>0", 1.0 / 6.0),  # complement
    ],
)
def test_cut_tet_volume_sign(cut_type, frac):
    """Cut tet quadrature produces sensible weights (sum < full volume)."""
    pts_arr = np.array(TET_PTS, dtype=np.float64).flatten()
    ls = np.array([-1.0, 1.0, 1.0, 1.0], dtype=np.float64)  # 1 vertex inside
    conn = np.array([0, 1, 2, 3], dtype=np.int32)
    offs = np.array([0, 4], dtype=np.int32)
    vtype = np.array([10], dtype=np.int32)
    rules = cutcells.runtime_quadrature(
        ls, pts_arr, conn, offs, vtype, cut_type, True, 3
    )
    s = rules.weights.sum()
    assert 0 < s < frac + 1e-10


def test_cut_tet_physical_points():
    """physical_points for a cut tet has correct shape."""
    pts_arr = np.array(TET_PTS, dtype=np.float64).flatten()
    ls = np.array([-1.0, 1.0, 1.0, 1.0], dtype=np.float64)
    conn = np.array([0, 1, 2, 3], dtype=np.int32)
    offs = np.array([0, 4], dtype=np.int32)
    vtype = np.array([10], dtype=np.int32)
    rules = cutcells.runtime_quadrature(
        ls, pts_arr, conn, offs, vtype, "phi<0", True, 3
    )
    phys = cutcells.physical_points(rules, pts_arr, conn, offs, vtype)
    assert phys.reshape(-1, 3).shape[1] == 3


def test_physical_points_rejects_non_vtk_parent_type_ids():
    """physical_points expects VTK parent ids, not internal CutCells type ids."""
    pts_arr = np.array(TET_PTS, dtype=np.float64).flatten()
    ls = np.full(4, -1.0, dtype=np.float64)
    conn = np.array([0, 1, 2, 3], dtype=np.int32)
    offs = np.array([0, 4], dtype=np.int32)
    vtk_tetra = np.array([10], dtype=np.int32)
    internal_tetra = np.array([3], dtype=np.int32)
    rules = cutcells.runtime_quadrature(
        ls, pts_arr, conn, offs, vtk_tetra, "phi<0", False, 2
    )

    with pytest.raises(ValueError, match="dimension does not match"):
        cutcells.physical_points(rules, pts_arr, conn, offs, internal_tetra)


# ---------------------------------------------------------------------------
# Fully outside → empty result
# ---------------------------------------------------------------------------


def test_fully_outside_empty():
    """A cell with all ls > 0 and cut_type='phi<0' produces zero quadrature points."""
    pts_arr = np.array(TET_PTS, dtype=np.float64).flatten()
    ls = np.full(4, 1.0, dtype=np.float64)  # all outside
    conn = np.array([0, 1, 2, 3], dtype=np.int32)
    offs = np.array([0, 4], dtype=np.int32)
    vtype = np.array([10], dtype=np.int32)
    rules = cutcells.runtime_quadrature(
        ls, pts_arr, conn, offs, vtype, "phi<0", True, 3
    )
    assert len(rules.weights) == 0


def test_fully_inside_all_ls_negative():
    """A cell with all ls < 0 and cut_type='phi<0' gives the full volume."""
    pts_arr = np.array(TET_PTS, dtype=np.float64).flatten()
    ls = np.full(4, -1.0, dtype=np.float64)
    conn = np.array([0, 1, 2, 3], dtype=np.int32)
    offs = np.array([0, 4], dtype=np.int32)
    vtype = np.array([10], dtype=np.int32)
    rules = cutcells.runtime_quadrature(
        ls, pts_arr, conn, offs, vtype, "phi<0", True, 4
    )
    assert abs(rules.weights.sum() - 1.0 / 6.0) < 1e-9
