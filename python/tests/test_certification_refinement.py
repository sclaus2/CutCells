import unittest

import numpy as np

import cutcells


def _single_triangle_mesh():
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float64)
    connectivity = np.array([0, 1, 2], dtype=np.int32)
    offsets = np.array([0, 3], dtype=np.int32)
    cell_types = np.array([5], dtype=np.int32)  # VTK_TRIANGLE
    return cutcells.MeshView(coords, connectivity, offsets, cell_types, tdim=2)


def _single_tetra_mesh():
    coords = np.array(
        [[0.0, 0.0, 0.0],
         [1.0, 0.0, 0.0],
         [0.0, 1.0, 0.0],
         [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    connectivity = np.array([0, 1, 2, 3], dtype=np.int32)
    offsets = np.array([0, 4], dtype=np.int32)
    cell_types = np.array([10], dtype=np.int32)  # VTK_TETRA
    return cutcells.MeshView(coords, connectivity, offsets, cell_types, tdim=3)




class CertificationRefinementTests(unittest.TestCase):
    def test_restrict_edge_bernstein_exact_matches_parent_evaluation(self):
        mesh = _single_triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] * X[0] + X[0] * X[1] + 0.25 * X[1] * X[1] - 0.1,
            degree=2,
        )
        ls_cell = cutcells.make_cell_level_set(ls, 0)

        xi_a = np.array([0.1, 0.15], dtype=np.float64)
        xi_b = np.array([0.7, 0.2], dtype=np.float64)
        edge_coeffs = cutcells.restrict_edge_bernstein_exact(
            ls_cell.cell_type,
            ls_cell.bernstein_order,
            np.asarray(ls_cell.bernstein_coeffs),
            xi_a,
            xi_b,
        )

        for t in np.linspace(0.0, 1.0, 9):
            xi = (1.0 - t) * xi_a + t * xi_b
            parent_val = cutcells.evaluate_bernstein(
                ls_cell.cell_type,
                ls_cell.bernstein_order,
                np.asarray(ls_cell.bernstein_coeffs),
                xi,
            )
            edge_val = cutcells.evaluate_bernstein(
                cutcells.CellType.interval,
                ls_cell.bernstein_order,
                np.asarray(edge_coeffs),
                np.array([t], dtype=np.float64),
            )
            self.assertAlmostEqual(edge_val, parent_val, places=10)

    def test_edge_endpoint_zero_semantics(self):
        one_root_tag, split = cutcells.classify_edge_roots(
            np.array([0.0, 0.0, 1.0], dtype=np.float64),
            max_depth=18,
        )
        self.assertEqual(one_root_tag, cutcells.EdgeRootTag.one_root)
        self.assertIsNone(split)

        multi_tag, split = cutcells.classify_edge_roots(
            np.array([0.0, -0.3, 0.4], dtype=np.float64),
            max_depth=18,
        )
        self.assertEqual(multi_tag, cutcells.EdgeRootTag.multiple_roots)
        self.assertGreater(split, 0.0)
        self.assertLess(split, 1.0)

    def test_green_refinement_on_multiple_root_edge(self):
        mesh = _single_triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: (X[0] - 0.2) * (X[0] - 0.8),
            degree=2,
        )
        adapt = cutcells.make_adapt_cell(mesh, 0)
        cutcells.build_edges(adapt)
        ls_cell = cutcells.make_cell_level_set(ls, 0)

        cutcells.classify_new_edges(adapt, ls_cell, 0)
        tags = np.asarray(adapt.edge_root_tags(0))
        self.assertGreaterEqual(
            np.count_nonzero(tags == cutcells.EdgeRootTag.multiple_roots.value),
            1,
        )

        changed = cutcells.refine_green_on_multiple_root_edges(adapt, 0)
        self.assertTrue(changed)
        self.assertEqual(adapt.num_vertices(), 4)
        self.assertEqual(adapt.num_cells(), 2)

    def test_green_refinement_only_invalidates_new_edges(self):
        mesh = _single_triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: (X[0] - 0.2) * (X[0] - 0.8),
            degree=2,
        )
        adapt = cutcells.make_adapt_cell(mesh, 0)
        cutcells.build_edges(adapt)
        ls_cell = cutcells.make_cell_level_set(ls, 0)

        cutcells.classify_new_edges(adapt, ls_cell, 0)
        before = np.asarray(adapt.edge_root_tags(0))
        self.assertEqual(
            np.count_nonzero(before == cutcells.EdgeRootTag.not_classified.value),
            0,
        )

        changed = cutcells.refine_green_on_multiple_root_edges(adapt, 0)
        self.assertTrue(changed)

        after_refine = np.asarray(adapt.edge_root_tags(0))
        self.assertGreater(
            np.count_nonzero(after_refine == cutcells.EdgeRootTag.not_classified.value),
            0,
        )
        self.assertGreater(
            np.count_nonzero(after_refine != cutcells.EdgeRootTag.not_classified.value),
            0,
        )

        cutcells.classify_new_edges(adapt, ls_cell, 0)
        after_recertify = np.asarray(adapt.edge_root_tags(0))
        self.assertEqual(
            np.count_nonzero(after_recertify == cutcells.EdgeRootTag.not_classified.value),
            0,
        )

    def test_red_refinement_on_ambiguous_triangle(self):
        mesh = _single_triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] * X[1] * (1.0 - X[0] - X[1]) - 1e-3,
            degree=3,
        )
        adapt = cutcells.make_adapt_cell(mesh, 0)
        cutcells.build_edges(adapt)
        ls_cell = cutcells.make_cell_level_set(ls, 0)

        cutcells.classify_new_edges(adapt, ls_cell, 0)
        cutcells.classify_leaf_cells(adapt, ls_cell, 0)
        self.assertEqual(
            np.asarray(adapt.cell_cert_tags(0)).tolist(),
            [cutcells.CellCertTag.ambiguous.value],
        )

        changed = cutcells.refine_red_on_ambiguous_cells(adapt, 0)
        self.assertTrue(changed)
        self.assertEqual(adapt.num_vertices(), 6)
        self.assertEqual(adapt.num_cells(), 4)
        self.assertEqual(adapt.num_edges(), 9)

    def test_triangle_cell_is_marked_ready_to_cut(self):
        mesh = _single_triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] + X[1] - 0.3,
            degree=1,
        )
        adapt = cutcells.make_adapt_cell(mesh, 0)
        cutcells.build_edges(adapt)
        ls_cell = cutcells.make_cell_level_set(ls, 0)

        cutcells.classify_new_edges(adapt, ls_cell, 0)
        cutcells.classify_leaf_cells(adapt, ls_cell, 0)
        self.assertEqual(
            np.asarray(adapt.cell_cert_tags(0)).tolist(),
            [cutcells.CellCertTag.ready_to_cut.value],
        )

        edge_tags = np.asarray(adapt.edge_root_tags(0))
        self.assertEqual(
            np.count_nonzero(edge_tags == cutcells.EdgeRootTag.one_root.value),
            2,
        )
        self.assertEqual(
            np.count_nonzero(edge_tags == cutcells.EdgeRootTag.multiple_roots.value),
            0,
        )
        self.assertEqual(
            np.count_nonzero(edge_tags == cutcells.EdgeRootTag.zero.value),
            0,
        )

    def test_triangle_ready_to_cut_cell_is_replaced_by_polygonal_lut_cells(self):
        mesh = _single_triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] + X[1] - 0.3,
            degree=1,
        )
        adapt = cutcells.make_adapt_cell(mesh, 0)
        cutcells.build_edges(adapt)
        ls_cell = cutcells.make_cell_level_set(ls, 0)

        cutcells.certify_refine_and_process_ready_cells(adapt, ls_cell, 0)
        self.assertEqual(adapt.num_cells(), 2)

        cell_types = np.asarray(adapt.cell_types)
        self.assertIn(cutcells.CellType.triangle.value, cell_types.tolist())
        self.assertIn(cutcells.CellType.quadrilateral.value, cell_types.tolist())

        tags = np.asarray(adapt.cell_cert_tags(0))
        self.assertEqual(
            sorted(np.unique(tags).tolist()),
            sorted([
                cutcells.CellCertTag.negative.value,
                cutcells.CellCertTag.positive.value,
            ]),
        )

    def test_ready_to_cut_triangle_uses_linear_one_root_vertices(self):
        mesh = _single_triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] * X[0] + X[1] - 0.25,
            degree=2,
        )
        adapt = cutcells.make_adapt_cell(mesh, 0)
        cutcells.build_edges(adapt)
        ls_cell = cutcells.make_cell_level_set(ls, 0)

        cutcells.certify_refine_and_process_ready_cells(adapt, ls_cell, 0)

        coords = np.asarray(adapt.vertex_coords).reshape(-1, 2)
        source_edges = np.asarray(adapt.vertex_source_edge_id)

        expected_linear_root = np.array([0.25, 0.0])
        expected_second_root = np.array([0.0, 0.25])

        d_linear = np.linalg.norm(coords - expected_linear_root, axis=1)
        d_second = np.linalg.norm(coords - expected_second_root, axis=1)

        linear_id = int(np.argmin(d_linear))
        second_id = int(np.argmin(d_second))

        self.assertLess(d_linear[linear_id], 1e-10)
        self.assertLess(d_second[second_id], 1e-10)
        self.assertGreaterEqual(int(source_edges[linear_id]), 0)
        self.assertGreaterEqual(int(source_edges[second_id]), 0)

    def test_ho_cut_smoke_uses_certified_triangle_pipeline(self):
        mesh = _single_triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] + X[1] - 0.3,
            degree=1,
            name="phi",
        )

        result = cutcells.ho_cut(mesh, ls)

        self.assertEqual(result.num_cut_cells, 1)
        np.testing.assert_array_equal(np.asarray(result.parent_cell_ids), np.array([0]))
        self.assertEqual(np.asarray(result.cell_domains).shape, (1, 1))

    def test_cut_overload_supports_mesh_part_selection(self):
        mesh = _single_triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] + X[1] - 0.3,
            degree=1,
            name="phi",
        )

        result = cutcells.cut(mesh, ls)
        negative = result["phi < 0"]
        interface = result["phi = 0"]

        self.assertEqual(result.num_cut_cells, 1)
        np.testing.assert_array_equal(np.asarray(negative.cut_cell_ids), np.array([0]))
        np.testing.assert_array_equal(np.asarray(interface.cut_cell_ids), np.array([0]))

        adapt = result.adapt_cell(0)
        self.assertGreater(adapt.num_cells(), 0)
        self.assertGreater(adapt.num_vertices(), 3)

    def test_zero_face_final_incidence_uses_face_to_cell_connectivity(self):
        mesh = _single_tetra_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] + X[1] - 0.5,
            degree=1,
            name="phi",
        )

        result = cutcells.ho_cut(
            mesh,
            ls,
        )
        adapt = result.adapt_cell(0)

        zero_dims = np.asarray(adapt.zero_entity_dim, dtype=np.int32)
        zero_ids = np.asarray(adapt.zero_entity_id, dtype=np.int32)
        zero_face_ids = zero_ids[zero_dims == 2]
        self.assertEqual(zero_face_ids.size, 1)

        face_id = int(zero_face_ids[0])
        face_to_cell = np.asarray(adapt.face_to_cell_connectivity, dtype=np.int32)
        face_offsets = np.asarray(adapt.face_to_cell_offsets, dtype=np.int32)
        incident_cells = face_to_cell[face_offsets[face_id]:face_offsets[face_id + 1]]
        self.assertEqual(incident_cells.size, 2)

        cell_types = np.asarray(adapt.cell_types, dtype=np.int32)
        self.assertEqual(cell_types[incident_cells].tolist(), [6, 6])

    def test_cut_detects_quadratic_interior_tetra_intersection(self):
        mesh = _single_tetra_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: (
                (X[0] - 0.2) * (X[0] - 0.2)
                + (X[1] - 0.2) * (X[1] - 0.2)
                + (X[2] - 0.2) * (X[2] - 0.2)
                - 0.09
            ),
            degree=2,
            name="phi",
        )

        vertex_vals = np.array(
            [
                (node[0] - 0.2) * (node[0] - 0.2)
                + (node[1] - 0.2) * (node[1] - 0.2)
                + (node[2] - 0.2) * (node[2] - 0.2)
                - 0.09
                for node in np.asarray(mesh.coordinates)
            ],
            dtype=np.float64,
        )
        self.assertTrue(np.all(vertex_vals > 0.0))

        ls_cell = cutcells.make_cell_level_set(ls, 0)
        self.assertLess(np.min(np.asarray(ls_cell.bernstein_coeffs)), 0.0)

        result = cutcells.cut(mesh, ls)
        interface = result["phi = 0"]

        self.assertEqual(result.num_cut_cells, 1)
        np.testing.assert_array_equal(np.asarray(interface.cut_cell_ids), np.array([0]))

    def test_certify_and_refine_runs_triangle_workflow(self):
        mesh = _single_triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] * X[1] * (1.0 - X[0] - X[1]) - 1e-3,
            degree=3,
        )
        adapt = cutcells.make_adapt_cell(mesh, 0)
        cutcells.build_edges(adapt)
        ls_cell = cutcells.make_cell_level_set(ls, 0)

        cutcells.certify_and_refine(adapt, ls_cell, 0, max_iterations=2)
        self.assertGreater(adapt.num_cells(), 1)
        self.assertGreater(adapt.num_vertices(), 3)

    def test_green_refinement_on_multiple_root_tetra_edge(self):
        mesh = _single_tetra_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: (X[0] - 0.1) * (X[0] - 0.6),
            degree=2,
        )
        adapt = cutcells.make_adapt_cell(mesh, 0)
        cutcells.build_edges(adapt)
        ls_cell = cutcells.make_cell_level_set(ls, 0)

        cutcells.classify_new_edges(adapt, ls_cell, 0)
        tags = np.asarray(adapt.edge_root_tags(0))
        splits = np.asarray(adapt.edge_green_split_params(0))
        edge_connectivity = np.asarray(adapt.edge_connectivity)
        edge_offsets = np.asarray(adapt.edge_offsets)
        self.assertGreaterEqual(
            np.count_nonzero(tags == cutcells.EdgeRootTag.multiple_roots.value),
            3,
        )

        split_edge = int(np.flatnonzero(tags == cutcells.EdgeRootTag.multiple_roots.value)[0])
        a, b = edge_connectivity[edge_offsets[split_edge]:edge_offsets[split_edge + 1]]
        split_t = float(splits[split_edge])
        self.assertNotAlmostEqual(split_t, 0.5, places=6)

        vertex_coords = np.asarray(adapt.vertex_coords).reshape(-1, 3)
        expected = (1.0 - split_t) * vertex_coords[a] + split_t * vertex_coords[b]

        changed = cutcells.refine_green_on_multiple_root_edges(adapt, 0)
        self.assertTrue(changed)
        self.assertEqual(adapt.num_vertices(), 5)
        self.assertEqual(adapt.num_cells(), 2)
        np.testing.assert_allclose(
            np.asarray(adapt.vertex_coords).reshape(-1, 3)[-1],
            expected,
            atol=1e-12,
        )

    def test_tetra_cell_is_marked_ready_to_cut_for_triangular_interface(self):
        mesh = _single_tetra_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] + X[1] + X[2] - 0.2,
            degree=1,
        )
        adapt = cutcells.make_adapt_cell(mesh, 0)
        cutcells.build_edges(adapt)
        ls_cell = cutcells.make_cell_level_set(ls, 0)

        cutcells.classify_new_edges(adapt, ls_cell, 0)
        cutcells.classify_leaf_cells(adapt, ls_cell, 0)
        self.assertEqual(
            np.asarray(adapt.cell_cert_tags(0)).tolist(),
            [cutcells.CellCertTag.ready_to_cut.value],
        )

        edge_tags = np.asarray(adapt.edge_root_tags(0))
        self.assertEqual(
            np.count_nonzero(edge_tags == cutcells.EdgeRootTag.one_root.value),
            3,
        )
        self.assertEqual(
            np.count_nonzero(edge_tags == cutcells.EdgeRootTag.multiple_roots.value),
            0,
        )
        self.assertEqual(
            np.count_nonzero(edge_tags == cutcells.EdgeRootTag.zero.value),
            0,
        )

    def test_tetra_cell_is_marked_ready_to_cut_for_quadrilateral_interface(self):
        mesh = _single_tetra_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] + X[1] - X[2] - 0.2,
            degree=1,
        )
        adapt = cutcells.make_adapt_cell(mesh, 0)
        cutcells.build_edges(adapt)
        ls_cell = cutcells.make_cell_level_set(ls, 0)

        cutcells.classify_new_edges(adapt, ls_cell, 0)
        cutcells.classify_leaf_cells(adapt, ls_cell, 0)
        self.assertEqual(
            np.asarray(adapt.cell_cert_tags(0)).tolist(),
            [cutcells.CellCertTag.ready_to_cut.value],
        )

        edge_tags = np.asarray(adapt.edge_root_tags(0))
        self.assertEqual(
            np.count_nonzero(edge_tags == cutcells.EdgeRootTag.one_root.value),
            4,
        )
        self.assertEqual(
            np.count_nonzero(edge_tags == cutcells.EdgeRootTag.multiple_roots.value),
            0,
        )
        self.assertEqual(
            np.count_nonzero(edge_tags == cutcells.EdgeRootTag.zero.value),
            0,
        )

    def test_tetra_ready_to_cut_cell_is_replaced_without_tetra_triangulation(self):
        mesh = _single_tetra_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] + X[1] + X[2] - 0.2,
            degree=1,
        )
        adapt = cutcells.make_adapt_cell(mesh, 0)
        cutcells.build_edges(adapt)
        ls_cell = cutcells.make_cell_level_set(ls, 0)

        cutcells.certify_refine_and_process_ready_cells(adapt, ls_cell, 0)
        self.assertGreaterEqual(adapt.num_cells(), 2)

        cell_types = np.asarray(adapt.cell_types)
        self.assertIn(cutcells.CellType.tetrahedron.value, cell_types.tolist())
        self.assertIn(cutcells.CellType.prism.value, cell_types.tolist())

        tags = np.asarray(adapt.cell_cert_tags(0))
        self.assertEqual(
            sorted(np.unique(tags).tolist()),
            sorted([
                cutcells.CellCertTag.negative.value,
                cutcells.CellCertTag.positive.value,
            ]),
        )

    def test_certify_and_refine_recursively_splits_tetra_multiple_root_edges(self):
        mesh = _single_tetra_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: (X[0] - 0.2) * (X[0] - 0.8),
            degree=2,
        )
        adapt = cutcells.make_adapt_cell(mesh, 0)
        cutcells.build_edges(adapt)
        ls_cell = cutcells.make_cell_level_set(ls, 0)

        cutcells.certify_and_refine(adapt, ls_cell, 0, max_iterations=8)
        self.assertGreater(adapt.num_cells(), 2)
        self.assertGreater(adapt.num_vertices(), 5)

        cutcells.classify_new_edges(adapt, ls_cell, 0)
        tags = np.asarray(adapt.edge_root_tags(0))
        self.assertEqual(
            np.count_nonzero(tags == cutcells.EdgeRootTag.multiple_roots.value),
            0,
        )


if __name__ == "__main__":
    unittest.main()
