import unittest

import numpy as np

import cutcells


def _mesh(coords, connectivity, vtk_type, tdim):
    coords = np.asarray(coords, dtype=np.float64)
    connectivity = np.asarray(connectivity, dtype=np.int32)
    offsets = np.array([0, len(connectivity)], dtype=np.int32)
    cell_types = np.array([vtk_type], dtype=np.int32)
    return cutcells.MeshView(coords, connectivity, offsets, cell_types, tdim=tdim)


def _interval_mesh():
    return _mesh([[0.0], [1.0]], [0, 1], 3, 1)


def _triangle_mesh():
    return _mesh([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], [0, 1, 2], 5, 2)


def _quad_mesh():
    return _mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
        [0, 1, 3, 2],
        9,
        2,
    )


def _tetra_mesh():
    return _mesh(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        [0, 1, 2, 3],
        10,
        3,
    )


def _hex_mesh():
    return _mesh(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ],
        [0, 1, 3, 2, 4, 5, 7, 6],
        12,
        3,
    )


class AdaptCellIsoRefineTests(unittest.TestCase):
    def test_template_counts_for_supported_v1_cells(self):
        cases = [
            (cutcells.CellType.interval, 4, 5, 4, cutcells.CellType.interval),
            (cutcells.CellType.triangle, 3, 10, 9, cutcells.CellType.triangle),
            (cutcells.CellType.quadrilateral, 2, 9, 4, cutcells.CellType.quadrilateral),
            (cutcells.CellType.tetrahedron, 3, 20, 27, cutcells.CellType.tetrahedron),
            (cutcells.CellType.hexahedron, 2, 27, 8, cutcells.CellType.hexahedron),
        ]

        for cell_type, order, n_vertices, n_cells, child_type in cases:
            with self.subTest(cell_type=cell_type, order=order):
                tpl = cutcells.iso_p1_template(cell_type, order)
                self.assertEqual(tpl.n_vertices, n_vertices)
                self.assertEqual(tpl.n_cells, n_cells)
                self.assertEqual(tpl.child_cell_type, child_type)
                self.assertEqual(len(np.asarray(tpl.vertex_parent_dim)), n_vertices)
                self.assertEqual(len(np.asarray(tpl.vertex_parent_id)), n_vertices)
                self.assertEqual(
                    len(np.asarray(tpl.cell_connectivity)),
                    n_cells * tpl.vertices_per_cell,
                )

    def test_make_adapt_cell_applies_iso_refine_once(self):
        cases = [
            (_interval_mesh(), 3, 4, 3),
            (_triangle_mesh(), 4, 15, 16),
            (_quad_mesh(), 3, 16, 9),
            (_tetra_mesh(), 2, 10, 8),
            (_hex_mesh(), 2, 27, 8),
        ]

        for mesh, order, n_vertices, n_cells in cases:
            with self.subTest(order=order, tdim=mesh.tdim):
                adapt = cutcells.make_adapt_cell(
                    mesh,
                    0,
                    cut_approximation="iso_p1",
                    cut_approximation_order=order,
                )
                self.assertEqual(adapt.num_vertices(), n_vertices)
                self.assertEqual(adapt.num_cells(), n_cells)
                np.testing.assert_array_equal(
                    np.asarray(adapt.cell_refinement_generation, dtype=np.int32),
                    np.ones(n_cells, dtype=np.int32),
                )
                np.testing.assert_array_equal(
                    np.asarray(adapt.cell_refinement_reason, dtype=np.uint8),
                    np.full(n_cells, 6, dtype=np.uint8),
                )

    def test_cut_option_uses_iso_refined_starting_topology(self):
        mesh = _triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] + X[1] - 0.9,
            degree=1,
            name="phi",
        )

        default = cutcells.cut(mesh, ls)
        refined = cutcells.cut(
            mesh,
            ls,
            cut_approximation="iso_p1",
            cut_approximation_order=3,
        )

        self.assertEqual(refined.num_cut_cells, 1)
        self.assertGreater(
            refined.adapt_cell(0).num_vertices(),
            default.adapt_cell(0).num_vertices(),
        )

    def test_higher_order_level_set_automatically_uses_iso_p1(self):
        mesh = _triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] + X[1] - 0.9,
            degree=3,
            name="phi",
        )

        automatic = cutcells.cut(mesh, ls)
        explicit = cutcells.cut(
            mesh,
            ls,
            cut_approximation="iso_p1",
            cut_approximation_order=3,
        )
        forced_linear = cutcells.cut(mesh, ls, cut_approximation="linear")

        self.assertEqual(automatic.num_cut_cells, 1)
        self.assertEqual(
            automatic.adapt_cell(0).num_vertices(),
            explicit.adapt_cell(0).num_vertices(),
        )
        self.assertGreater(
            automatic.adapt_cell(0).num_vertices(),
            forced_linear.adapt_cell(0).num_vertices(),
        )

    def test_invalid_linear_order_is_rejected(self):
        mesh = _triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] + X[1] - 0.9,
            degree=1,
            name="phi",
        )

        with self.assertRaisesRegex(Exception, "linear"):
            cutcells.cut(
                mesh,
                ls,
                cut_approximation="linear",
                cut_approximation_order=2,
            )


if __name__ == "__main__":
    unittest.main()
