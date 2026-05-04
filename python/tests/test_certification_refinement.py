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


def _tetra_line_interval(seed, direction):
    lo = -np.inf
    hi = np.inf
    constraints = [
        (np.array([1.0, 0.0, 0.0]), 0.0),
        (np.array([0.0, 1.0, 0.0]), 0.0),
        (np.array([0.0, 0.0, 1.0]), 0.0),
        (np.array([-1.0, -1.0, -1.0]), -1.0),
    ]
    for normal, lower in constraints:
        value = float(np.dot(normal, seed) - lower)
        slope = float(np.dot(normal, direction))
        if abs(slope) <= 1.0e-14:
            if value < -1.0e-14:
                return np.nan, np.nan
            continue
        bound = -value / slope
        if slope > 0.0:
            lo = max(lo, bound)
        else:
            hi = min(hi, bound)
    return lo, hi


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

    def test_ready_to_cut_triangle_uses_exact_one_root_vertices(self):
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

        expected_curved_root = np.array([0.5, 0.0])
        expected_linear_root = np.array([0.25, 0.0])
        expected_second_root = np.array([0.0, 0.25])

        d_curved = np.linalg.norm(coords - expected_curved_root, axis=1)
        d_second = np.linalg.norm(coords - expected_second_root, axis=1)
        d_linear = np.linalg.norm(coords - expected_linear_root, axis=1)

        curved_id = int(np.argmin(d_curved))
        second_id = int(np.argmin(d_second))

        self.assertLess(d_curved[curved_id], 1e-10)
        self.assertLess(d_second[second_id], 1e-10)
        self.assertGreater(np.min(d_linear), 1e-2)
        self.assertGreaterEqual(int(source_edges[curved_id]), 0)
        self.assertGreaterEqual(int(source_edges[second_id]), 0)

    def test_graph_check_failure_refines_original_uncut_ready_triangle(self):
        mesh = _single_triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] * X[0] + X[1] - 0.25,
            degree=2,
        )
        adapt = cutcells.make_adapt_cell(mesh, 0)
        cutcells.build_edges(adapt)
        ls_cell = cutcells.make_cell_level_set(ls, 0)

        diagnostics = cutcells.certify_refine_graph_check_and_process_ready_cells(
            adapt,
            ls_cell,
            0,
            projection_direction="straight_zero_entity_normal",
            graph_max_refinements=2,
            max_relative_correction_distance=0.01,
        )

        self.assertFalse(bool(diagnostics["accepted"]))
        self.assertEqual(int(diagnostics["graph_refinements"]), 2)
        self.assertEqual(
            str(diagnostics["first_failure_reason"]),
            "excessive_correction_distance",
        )
        self.assertGreater(adapt.num_cells(), 2)

        tags = np.asarray(adapt.cell_cert_tags(0))
        self.assertEqual(
            sorted(np.unique(tags).tolist()),
            sorted([
                cutcells.CellCertTag.negative.value,
                cutcells.CellCertTag.positive.value,
                cutcells.CellCertTag.cut.value,
            ]),
        )

    def test_ho_cut_smoke_uses_certified_triangle_pipeline(self):
        mesh = _single_triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] + X[1] - 0.3,
            degree=1,
            name="phi",
        )

        result = cutcells.ho_cut(mesh, ls, graph_enabled=False)

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

    def test_graph_check_data_is_attached_to_zero_interface_cells(self):
        mesh = _single_triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] * X[0] + X[1] - 0.25,
            degree=2,
            name="phi",
        )

        result = cutcells.cut(
            mesh,
            ls,
            graph_max_refinements=0,
            min_level_set_gradient_host_alignment=0.0,
        )
        interface = result["phi = 0"]
        zero_mesh = interface.visualization_mesh(mode="cut_only", geometry_order=1)
        data = interface.graph_check_zero_entity_data()
        summary = result.graph_check_summary()

        num_zero_cells = np.asarray(zero_mesh.offset).size - 1
        self.assertGreater(num_zero_cells, 0)
        self.assertTrue(np.all(np.asarray(summary["graph_refinements"]) == 0))
        self.assertEqual(np.asarray(data["local_zero_entity_id"]).shape, (num_zero_cells,))
        self.assertEqual(np.asarray(data["graph_accepted"]).shape, (num_zero_cells,))
        self.assertTrue(np.all(np.asarray(data["zero_entity_dim"]) == interface.dim))
        self.assertTrue(np.all(np.asarray(data["graph_accepted"]) >= 0))
        self.assertTrue(np.all(np.isfinite(np.asarray(data["graph_max_correction"]))))
        for key in (
            "graph_failed_projection_seed",
            "graph_failed_projection_direction",
            "graph_failed_projection_clip_lo",
            "graph_failed_projection_clip_hi",
            "graph_failed_projection_root_t",
        ):
            self.assertIn(key, data)
        self.assertEqual(
            np.asarray(data["graph_failed_projection_seed"]).shape,
            (num_zero_cells, 2),
        )
        self.assertEqual(
            np.asarray(data["graph_failed_projection_direction"]).shape,
            (num_zero_cells, 2),
        )
        self.assertEqual(
            np.asarray(data["graph_failed_projection_clip_lo"]).shape,
            (num_zero_cells,),
        )

    def test_sphere_graph_check_reports_surface_jacobian_and_node_data(self):
        grid = cutcells.box_tetrahedron_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 3, 3, 3)
        mesh = cutcells.mesh_from_pyvista(grid)
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] * X[0] + X[1] * X[1] + X[2] * X[2] - 0.36,
            degree=2,
            name="phi",
        )

        result = cutcells.cut(mesh, ls, graph_max_refinements=0)
        interface = result["phi = 0"]
        data = interface.graph_check_zero_entity_data()
        node_data = interface.graph_check_node_data()

        self.assertIn("graph_min_level_set_gradient_host_alignment", data)
        self.assertIn("graph_min_surface_jacobian_ratio", data)
        self.assertIn("graph_failed_surface_jacobian_ratio", data)

        node_accepted = np.asarray(node_data["node_accepted"], dtype=np.int32)
        self.assertGreater(node_accepted.size, 0)
        gradient_alignment = np.asarray(
            node_data["level_set_gradient_host_alignment"], dtype=np.float64)
        finite = np.isfinite(gradient_alignment)
        self.assertTrue(np.any(finite))
        self.assertTrue(np.all(gradient_alignment[finite] >= 0.0))

    def test_cut_graph_check_accepts_level_set_gradient_mode(self):
        mesh = _single_triangle_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] * X[0] + X[1] - 0.25,
            degree=2,
            name="phi",
        )

        result = cutcells.cut(
            mesh,
            ls,
            graph_max_refinements=0,
            graph_projection_direction="level_set_gradient",
        )
        data = result["phi = 0"].graph_check_zero_entity_data()

        self.assertTrue(np.all(np.asarray(result.graph_check_summary()["graph_refinements"]) == 0))
        self.assertGreater(np.asarray(data["graph_accepted"]).size, 0)

    def test_cut_graph_check_accepts_red_failed_cell_refinement_mode(self):
        grid = cutcells.box_tetrahedron_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 3, 3, 3)
        mesh = cutcells.mesh_from_pyvista(grid)
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] * X[0] + X[1] * X[1] + X[2] * X[2] - 0.36,
            degree=2,
            name="phi",
        )

        result = cutcells.cut(
            mesh,
            ls,
            graph_max_refinements=1,
            graph_refinement_mode="red_failed_cell",
        )
        summary = result.graph_check_summary()

        self.assertGreater(np.asarray(summary["accepted"]).size, 0)
        self.assertTrue(np.all(np.asarray(summary["graph_refinements"]) >= 0))

    def test_cut_graph_check_accepts_orthogonal_surface_edge_refinement_mode(self):
        grid = cutcells.box_tetrahedron_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 3, 3, 3)
        mesh = cutcells.mesh_from_pyvista(grid)
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] * X[0] + X[1] * X[1] + X[2] * X[2] - 0.36,
            degree=2,
            name="phi",
        )

        result = cutcells.cut(
            mesh,
            ls,
            graph_max_refinements=1,
            graph_refinement_mode="green_orthogonal_surface_edge",
        )
        summary = result.graph_check_summary()

        self.assertGreater(np.asarray(summary["accepted"]).size, 0)
        self.assertTrue(np.all(np.asarray(summary["graph_refinements"]) >= 0))

    def test_cut_graph_check_accepts_surface_error_refinement_modes(self):
        grid = cutcells.box_tetrahedron_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 3, 3, 3)
        mesh = cutcells.mesh_from_pyvista(grid)
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] * X[0] + X[1] * X[1] + X[2] * X[2] - 0.36,
            degree=2,
            name="phi",
        )

        for mode in ("green_midpoint_residual", "green_normal_variation"):
            result = cutcells.cut(
                mesh,
                ls,
                graph_max_refinements=1,
                graph_refinement_mode=mode,
            )
            summary = result.graph_check_summary()
            self.assertGreater(np.asarray(summary["accepted"]).size, 0)
            self.assertTrue(np.all(np.asarray(summary["graph_refinements"]) >= 0))

    def test_scalar_projection_records_parent_clipped_tetra_bracket(self):
        mesh = _single_tetra_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] * X[0] + X[1] * X[1] + X[2] * X[2] - 0.36,
            degree=2,
            name="phi",
        )

        result = cutcells.ho_cut(
            mesh,
            ls,
            min_level_set_gradient_host_alignment=0.0,
        )
        curved = result.curved_zero_nodes(geometry_order=2, node_family="gll")

        seeds = np.asarray(curved["node_seed"], dtype=np.float64)
        directions = np.asarray(curved["node_direction"], dtype=np.float64)
        clip_lo = np.asarray(curved["node_clip_lo"], dtype=np.float64)
        clip_hi = np.asarray(curved["node_clip_hi"], dtype=np.float64)
        root_t = np.asarray(curved["node_root_t"], dtype=np.float64)

        projected = np.flatnonzero(np.isfinite(root_t))
        self.assertGreater(projected.size, 0)
        for idx in projected:
            seed = seeds[idx]
            direction = directions[idx]
            lo, hi = _tetra_line_interval(seed, direction)
            np.testing.assert_allclose([clip_lo[idx], clip_hi[idx]], [lo, hi], atol=1e-12)
            self.assertFalse(np.isclose(clip_lo[idx], -2.0))
            self.assertFalse(np.isclose(clip_hi[idx], 2.0))
            self.assertLessEqual(clip_lo[idx] - 1e-12, root_t[idx])
            self.assertLessEqual(root_t[idx], clip_hi[idx] + 1e-12)

            phi_seed = float(np.dot(seed, seed) - 0.36)
            grad_seed = 2.0 * seed
            self.assertGreaterEqual(float(np.dot(grad_seed, direction)), -1e-12)
            if phi_seed > 0.0:
                self.assertLessEqual(root_t[idx], 1e-12)
            elif phi_seed < 0.0:
                self.assertGreaterEqual(root_t[idx], -1e-12)

    def test_straight_normal_projection_is_gradient_oriented(self):
        mesh = _single_tetra_mesh()
        ls = cutcells.create_level_set(
            mesh,
            lambda X: X[0] * X[0] + X[1] * X[1] + X[2] * X[2] - 0.36,
            degree=2,
            name="phi",
        )

        result = cutcells.ho_cut(
            mesh,
            ls,
            min_level_set_gradient_host_alignment=0.0,
        )
        curved = result.curved_zero_nodes(
            geometry_order=2,
            node_family="gll",
            projection_direction="straight_zero_entity_normal",
        )

        seeds = np.asarray(curved["node_seed"], dtype=np.float64)
        directions = np.asarray(curved["node_direction"], dtype=np.float64)
        root_t = np.asarray(curved["node_root_t"], dtype=np.float64)
        projected = np.flatnonzero(np.isfinite(root_t))
        self.assertGreater(projected.size, 0)
        for idx in projected:
            self.assertGreaterEqual(float(np.dot(2.0 * seeds[idx], directions[idx])), -1e-12)

    def test_triangulated_zero_quad_faces_curve_across_artificial_diagonal(self):
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
            triangulate=True,
            graph_max_refinements=0,
            min_level_set_gradient_host_alignment=0.0,
        )
        curved = result.curved_zero_nodes(geometry_order=2, node_family="gll")
        graph_nodes = result["phi = 0"].graph_check_node_data()

        dim = np.asarray(curved["dim"], dtype=np.int32)
        status = np.asarray(curved["status"], dtype=np.uint8)
        edge_states = np.flatnonzero(dim == 1)
        face_states = np.flatnonzero(dim == 2)

        self.assertGreaterEqual(edge_states.size, 5)
        self.assertTrue(np.all(status[edge_states] == 2))  # CurvingStatus::curved
        self.assertGreaterEqual(face_states.size, 2)
        self.assertTrue(np.all(status[face_states] == 2))  # CurvingStatus::curved
        face_interior_nodes = (
            np.asarray(graph_nodes["node_kind"], dtype=np.int32) == 3
        )
        self.assertTrue(np.any(
            np.asarray(graph_nodes["node_index"], dtype=np.int32)[face_interior_nodes] > 0
        ))

    def test_curving_accepts_near_zero_multi_level_set_node_without_projection(self):
        mesh = _single_triangle_mesh()

        def phi(X):
            return X[0] + X[1] - 0.25 + 5.12e-11 * X[0] * X[1]

        ls0 = cutcells.create_level_set(mesh, phi, degree=2, name="phi")
        ls1 = cutcells.create_level_set(mesh, phi, degree=2, name="psi")

        result = cutcells.ho_cut(mesh, [ls0, ls1])
        curved = result.curved_zero_nodes(geometry_order=2, node_family="gll")

        status = np.asarray(curved["status"], dtype=np.uint8)
        dim = np.asarray(curved["dim"], dtype=np.int32)
        zero_mask = np.asarray(curved["zero_mask"], dtype=np.uint64)
        offsets = np.asarray(curved["offsets"], dtype=np.int32)
        stats_offsets = np.asarray(curved["stats_offsets"], dtype=np.int32)
        node_iterations = np.asarray(curved["node_iterations"], dtype=np.int32)
        node_residual = np.asarray(curved["node_residual"], dtype=np.float64)
        node_projection_mode = np.asarray(curved["node_projection_mode"], dtype=np.uint8)
        points = np.asarray(curved["points"], dtype=np.float64)

        self.assertTrue(np.all(status == 2))  # CurvingStatus::curved
        edge_state = int(np.flatnonzero((dim == 1) & (zero_mask == 3))[0])
        node_begin = int(offsets[edge_state])
        node_end = int(offsets[edge_state + 1])
        stat_begin = int(stats_offsets[edge_state])
        stat_end = int(stats_offsets[edge_state + 1])

        self.assertEqual(node_end - node_begin, 3)
        self.assertEqual(stat_end - stat_begin, 3)
        mid_node = node_begin + 1
        mid_stat = stat_begin + 1

        np.testing.assert_allclose(points[mid_node], [0.125, 0.125], atol=1e-15)
        self.assertLess(node_residual[mid_stat], 1.0e-12)
        self.assertEqual(int(node_iterations[mid_stat]), 0)
        self.assertEqual(int(node_projection_mode[mid_stat]), 0)  # CurvingProjectionMode::none

    def test_curving_keeps_only_short_edge_faces_straight(self):
        mesh = _single_tetra_mesh()
        eps = 1.0e-4

        def phi(X):
            return X[0] + X[1] - eps

        ls = cutcells.create_level_set(mesh, phi, degree=1, name="phi")
        result = cutcells.ho_cut(mesh, ls, graph_enabled=False)
        curved = result.curved_zero_nodes(
            geometry_order=2,
            node_family="gll",
            small_entity_tol=2.0e-2,
        )

        dim = np.asarray(curved["dim"], dtype=np.int32)
        stats_offsets = np.asarray(curved["stats_offsets"], dtype=np.int32)
        node_failure_code = np.asarray(curved["node_failure_code"], dtype=np.uint8)
        small_code = list(curved["failure_code_names"]).index("small_entity_kept_straight")

        face_states = np.flatnonzero(dim == 2)
        self.assertGreater(len(face_states), 0)
        for face_state in face_states:
            begin = int(stats_offsets[face_state])
            end = int(stats_offsets[face_state + 1])
            self.assertFalse(np.all(node_failure_code[begin:end] == small_code))

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
