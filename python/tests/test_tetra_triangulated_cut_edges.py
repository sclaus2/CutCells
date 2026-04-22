import numpy as np

import cutcells


EDGE_PATTERNS = {
    5: ((0, 1), (1, 2), (2, 0)),  # triangle
    9: ((0, 1), (1, 2), (2, 3), (3, 0)),  # quad
    10: ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)),  # tet
    13: (
        (0, 1),
        (1, 2),
        (2, 0),
        (3, 4),
        (4, 5),
        (5, 3),
        (0, 3),
        (1, 4),
        (2, 5),
    ),  # prism / wedge
}


def _single_tetra_mesh():
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )
    connectivity = np.array([0, 1, 2, 3], dtype=np.int32)
    offsets = np.array([0, 4], dtype=np.int32)
    cell_types = np.array([10], dtype=np.int32)
    return cutcells.MeshView(coords, connectivity, offsets, cell_types, tdim=3), coords


def _edge_key(xa: np.ndarray, xb: np.ndarray, digits: int = 12):
    pa = tuple(np.round(np.asarray(xa, dtype=np.float64), digits))
    pb = tuple(np.round(np.asarray(xb, dtype=np.float64), digits))
    return tuple(sorted((pa, pb)))


def _collect_unique_edges(cut_cell):
    points = np.asarray(cut_cell.vertex_coords, dtype=np.float64)
    cells = np.asarray(cut_cell.cells, dtype=np.int64)
    vtk_types = np.asarray(cut_cell.vtk_types, dtype=np.int64)

    edges = {}
    cursor = 0
    for vtk_type in vtk_types:
        num_verts = int(cells[cursor])
        local_ids = cells[cursor + 1 : cursor + 1 + num_verts]
        cursor += num_verts + 1

        for local_a, local_b in EDGE_PATTERNS[int(vtk_type)]:
            xa = points[int(local_ids[local_a])]
            xb = points[int(local_ids[local_b])]
            edges[_edge_key(xa, xb)] = (
                np.array(xa, dtype=np.float64, copy=True),
                np.array(xb, dtype=np.float64, copy=True),
            )

    return edges


def _classify_new_edges(ls_cell, base_cut, triangulated_cut):
    base_edges = _collect_unique_edges(base_cut)
    triangulated_edges = _collect_unique_edges(triangulated_cut)

    records = []
    for key in sorted(set(triangulated_edges) - set(base_edges)):
        xa, xb = triangulated_edges[key]
        coeffs = np.asarray(
            cutcells.restrict_edge_bernstein_exact(
                ls_cell.cell_type,
                ls_cell.bernstein_order,
                np.asarray(ls_cell.bernstein_coeffs),
                np.asarray(xa, dtype=np.float64),
                np.asarray(xb, dtype=np.float64),
            ),
            dtype=np.float64,
        )
        tag, _ = cutcells.classify_edge_roots(
            coeffs,
            zero_tol=1.0e-12,
            sign_tol=1.0e-12,
            max_depth=20,
        )
        phi_a = float(
            cutcells.evaluate_bernstein(
                ls_cell.cell_type,
                ls_cell.bernstein_order,
                np.asarray(ls_cell.bernstein_coeffs),
                np.asarray(xa, dtype=np.float64),
            )
        )
        phi_b = float(
            cutcells.evaluate_bernstein(
                ls_cell.cell_type,
                ls_cell.bernstein_order,
                np.asarray(ls_cell.bernstein_coeffs),
                np.asarray(xb, dtype=np.float64),
            )
        )
        has_zero_endpoint = abs(phi_a) <= 1.0e-12 or abs(phi_b) <= 1.0e-12
        records.append(
            {
                "tag": tag,
                "has_zero_endpoint": has_zero_endpoint,
                "has_interior_root": tag == cutcells.EdgeRootTag.multiple_roots
                or (tag == cutcells.EdgeRootTag.one_root and not has_zero_endpoint),
            }
        )

    return records


def test_triangulated_tetra_cut_edges_are_classified_against_zero_level_set():
    mesh, coords = _single_tetra_mesh()
    nodal_phi = coords[:, 0] + coords[:, 1] - 0.6
    ls = cutcells.create_level_set(mesh, lambda X: X[0] + X[1] - 0.6, degree=1, name="phi")
    ls_cell = cutcells.make_cell_level_set(ls, 0)

    expected = {
        "phi<0": {
            "base_types": [13],
            "triangulated_types": [10, 10, 10, 10],
            "new_edge_count": 7,
            "zero_tag_count": 1,
            "zero_endpoint_count": 5,
            "interior_root_count": 0,
        },
        "phi=0": {
            "base_types": [9],
            "triangulated_types": [5, 5],
            "new_edge_count": 1,
            "zero_tag_count": 1,
            "zero_endpoint_count": 1,
            "interior_root_count": 0,
        },
        "phi>0": {
            "base_types": [13],
            "triangulated_types": [10, 10, 10, 10],
            "new_edge_count": 7,
            "zero_tag_count": 1,
            "zero_endpoint_count": 5,
            "interior_root_count": 0,
        },
    }

    for expr, exp in expected.items():
        base_cut = cutcells.cut(
            cutcells.CellType.tetrahedron,
            coords.ravel(),
            3,
            nodal_phi,
            expr,
            False,
        )
        triangulated_cut = cutcells.cut(
            cutcells.CellType.tetrahedron,
            coords.ravel(),
            3,
            nodal_phi,
            expr,
            True,
        )

        assert np.asarray(base_cut.vtk_types).tolist() == exp["base_types"]
        assert np.asarray(triangulated_cut.vtk_types).tolist() == exp["triangulated_types"]

        edge_records = _classify_new_edges(ls_cell, base_cut, triangulated_cut)
        assert len(edge_records) == exp["new_edge_count"]
        assert sum(record["tag"] == cutcells.EdgeRootTag.zero for record in edge_records) == exp[
            "zero_tag_count"
        ]
        assert sum(record["has_zero_endpoint"] for record in edge_records) == exp[
            "zero_endpoint_count"
        ]
        assert sum(record["has_interior_root"] for record in edge_records) == exp[
            "interior_root_count"
        ]
