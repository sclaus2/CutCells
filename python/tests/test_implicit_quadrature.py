import math

import numpy as np
import pytest

import cutcells


def _single_triangle_mesh():
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        dtype=np.float64,
    )
    connectivity = np.array([0, 1, 2], dtype=np.int32)
    offsets = np.array([0, 3], dtype=np.int32)
    cell_types = np.array([5], dtype=np.int32)
    return cutcells.MeshView(coords, connectivity, offsets, cell_types, tdim=2)


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
    return cutcells.MeshView(coords, connectivity, offsets, cell_types, tdim=3)


def _quadratic_sphere_mesh(n: int = 5):
    center = np.array([0.07, -0.04, 0.03], dtype=np.float64)
    radius = 0.55
    grid = cutcells.box_tetrahedron_mesh(
        -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n, n, n
    )
    mesh = cutcells.mesh_from_pyvista(grid, tdim=3)
    ls = cutcells.interpolate_level_set(
        mesh,
        lambda X: np.sum((X - center) ** 2, axis=1) - radius**2,
        degree=2,
        name="phi",
    )
    return mesh, ls


def _parent_chart_signatures(rules):
    signatures = {}
    for parent, local_cell, path, depth in zip(
        np.asarray(rules.parent_map, dtype=np.int32),
        np.asarray(rules.debug_local_cell_id, dtype=np.int32),
        np.asarray(rules.debug_chart_path, dtype=np.int32),
        np.asarray(rules.debug_refinement_depth, dtype=np.int32),
    ):
        signatures.setdefault(int(parent), []).append(
            (int(local_cell), int(path), int(depth))
        )
    return {key: tuple(sorted(value)) for key, value in signatures.items()}


def _parent_plan_hash_signatures(rules):
    signatures = {}
    for parent, plan_hash in zip(
        np.asarray(rules.parent_map, dtype=np.int32),
        np.asarray(rules.debug_chart_plan_hash, dtype=np.int64),
    ):
        signatures.setdefault(int(parent), []).append(int(plan_hash))
    return {key: tuple(sorted(value)) for key, value in signatures.items()}


def test_implicit_options_are_bound_and_used():
    opts = cutcells.ImplicitQuadratureOptions()
    opts.max_refinement_iterations = 10
    opts.root_method = cutcells.EdgeRootMethod.itp
    opts.enable_height_direction_recovery = False
    opts.max_surface_q_rel_error = 0.12
    opts.min_surface_area_ratio = 0.03
    opts.max_surface_area_ratio = 30.0
    opts.max_diagnostic_refinement_depth = 5
    opts.min_height_area_ratio = 0.2
    opts.max_height_area_ratio = 5.0

    assert opts.order == 6
    assert opts.max_refinement_iterations == 10
    assert opts.root_method == cutcells.EdgeRootMethod.itp
    assert opts.enable_height_direction_recovery is False
    assert opts.max_surface_q_rel_error == pytest.approx(0.12)
    assert opts.min_surface_area_ratio == pytest.approx(0.03)
    assert opts.max_surface_area_ratio == pytest.approx(30.0)
    assert opts.max_diagnostic_refinement_depth == 5
    assert opts.min_height_area_ratio == pytest.approx(0.2)
    assert opts.max_height_area_ratio == pytest.approx(5.0)


def test_tetra_sphere_chart_path_is_independent_of_quadrature_order():
    mesh, ls = _quadratic_sphere_mesh(n=5)
    result = cutcells.cut(mesh, ls)

    signatures_by_order = {}
    plan_hashes_by_order = {}
    for order in (4, 6, 8, 10):
        opts = cutcells.ImplicitQuadratureOptions()
        opts.root_tol = 1.0e-13
        rules = result["phi = 0"].implicit_quadrature(
            order=order,
            mode="cut_only",
            options=opts,
        )
        signatures_by_order[order] = _parent_chart_signatures(rules)
        plan_hashes_by_order[order] = _parent_plan_hash_signatures(rules)
        num_rules = len(np.asarray(rules.parent_map, dtype=np.int32))
        assert len(np.asarray(rules.debug_chart_plan_hash, dtype=np.int64)) == num_rules
        assert len(np.asarray(rules.debug_candidate_mask, dtype=np.int32)) == num_rules
        assert len(np.asarray(rules.debug_rejection_reason, dtype=np.int32)) == num_rules
        assert len(np.asarray(rules.debug_measure_probe, dtype=np.float64)) == num_rules
        assert (
            len(np.asarray(rules.debug_validation_weight_sum, dtype=np.float64))
            == num_rules
        )

    assert signatures_by_order[4] == signatures_by_order[6]
    assert signatures_by_order[6] == signatures_by_order[8]
    assert signatures_by_order[8] == signatures_by_order[10]
    assert plan_hashes_by_order[4] == plan_hashes_by_order[6]
    assert plan_hashes_by_order[6] == plan_hashes_by_order[8]
    assert plan_hashes_by_order[8] == plan_hashes_by_order[10]


def test_triangle_linear_cut_matches_exact_area_and_length():
    mesh = _single_triangle_mesh()
    cut_at = 0.6
    ls = cutcells.interpolate_level_set(
        mesh,
        lambda X: X[:, 0] + X[:, 1] - cut_at,
        degree=1,
        name="phi",
    )

    result = cutcells.cut(mesh, ls)
    negative = result["phi < 0"].implicit_quadrature(order=5, mode="cut_only")
    positive = result["phi > 0"].implicit_quadrature(order=5, mode="cut_only")
    interface = result["phi = 0"].implicit_quadrature(order=5, mode="cut_only")

    assert np.asarray(negative.weights).sum() == pytest.approx(cut_at**2 / 2)
    assert np.asarray(positive.weights).sum() == pytest.approx(0.5 - cut_at**2 / 2)
    assert np.asarray(interface.weights).sum() == pytest.approx(math.sqrt(2) * cut_at)


def test_triangle_height_recovery_handles_bad_sign_chart_direction():
    mesh = _single_triangle_mesh()
    cut_at = 0.2
    ls = cutcells.interpolate_level_set(
        mesh,
        lambda X: X[:, 0] - cut_at,
        degree=1,
        name="phi",
    )

    opts = cutcells.ImplicitQuadratureOptions()
    opts.min_transversality = 1.0
    opts.min_height_transversality = 0.05
    opts.enable_height_direction_recovery = True

    result = cutcells.cut(mesh, ls)
    negative = result["phi < 0"].implicit_quadrature(
        order=5, mode="cut_only", options=opts
    )
    positive = result["phi > 0"].implicit_quadrature(
        order=5, mode="cut_only", options=opts
    )
    interface = result["phi = 0"].implicit_quadrature(
        order=5, mode="cut_only", options=opts
    )

    assert np.asarray(negative.weights).sum() == pytest.approx(0.18)
    assert np.asarray(positive.weights).sum() == pytest.approx(0.32)
    assert np.asarray(interface.weights).sum() == pytest.approx(0.8)


def test_triangle_soft_transversality_does_not_force_recovery():
    mesh = _single_triangle_mesh()
    ls = cutcells.interpolate_level_set(
        mesh,
        lambda X: X[:, 0] - 0.2,
        degree=1,
        name="phi",
    )

    opts = cutcells.ImplicitQuadratureOptions()
    opts.min_transversality = 1.0
    opts.enable_height_direction_recovery = False

    result = cutcells.cut(mesh, ls)
    rules = result["phi < 0"].implicit_quadrature(
        order=5, mode="cut_only", options=opts
    )
    assert np.asarray(rules.weights).sum() == pytest.approx(0.18)


def test_tetra_linear_cut_matches_exact_volume_and_area():
    mesh = _single_tetra_mesh()
    cut_at = 0.6
    ls = cutcells.interpolate_level_set(
        mesh,
        lambda X: X[:, 0] + X[:, 1] + X[:, 2] - cut_at,
        degree=1,
        name="phi",
    )

    result = cutcells.cut(mesh, ls)
    negative = result["phi < 0"].implicit_quadrature(order=5, mode="cut_only")
    interface = result["phi = 0"].implicit_quadrature(order=5, mode="cut_only")

    assert np.asarray(negative.weights).sum() == pytest.approx(cut_at**3 / 6)
    assert np.asarray(interface.weights).sum() == pytest.approx(
        math.sqrt(3) * cut_at**2 / 2
    )


def test_tetra_height_recovery_handles_bad_sign_chart_direction():
    mesh = _single_tetra_mesh()
    cut_at = 0.2
    ls = cutcells.interpolate_level_set(
        mesh,
        lambda X: X[:, 0] - cut_at,
        degree=1,
        name="phi",
    )

    opts = cutcells.ImplicitQuadratureOptions()
    opts.min_transversality = 1.0
    opts.min_height_transversality = 0.05
    opts.enable_height_direction_recovery = True

    result = cutcells.cut(mesh, ls)
    negative = result["phi < 0"].implicit_quadrature(
        order=5, mode="cut_only", options=opts
    )
    positive = result["phi > 0"].implicit_quadrature(
        order=5, mode="cut_only", options=opts
    )
    interface = result["phi = 0"].implicit_quadrature(
        order=5, mode="cut_only", options=opts
    )

    exact_negative = (1.0 - (1.0 - cut_at) ** 3) / 6.0
    exact_positive = 1.0 / 6.0 - exact_negative
    exact_interface = (1.0 - cut_at) ** 2 / 2.0

    assert np.asarray(negative.weights).sum() == pytest.approx(exact_negative)
    assert np.asarray(positive.weights).sum() == pytest.approx(exact_positive)
    assert np.asarray(interface.weights).sum() == pytest.approx(exact_interface)


def test_tetra_soft_transversality_does_not_force_recovery():
    mesh = _single_tetra_mesh()
    ls = cutcells.interpolate_level_set(
        mesh,
        lambda X: X[:, 0] - 0.2,
        degree=1,
        name="phi",
    )

    opts = cutcells.ImplicitQuadratureOptions()
    opts.min_transversality = 1.0
    opts.enable_height_direction_recovery = False

    result = cutcells.cut(mesh, ls)
    rules = result["phi < 0"].implicit_quadrature(
        order=5, mode="cut_only", options=opts
    )
    assert np.asarray(rules.weights).sum() == pytest.approx(0.08133333333333333)


def test_triangle_circle_cut_matches_analytic_reference_and_roots():
    mesh = _single_triangle_mesh()
    radius = 0.45
    ls = cutcells.interpolate_level_set(
        mesh,
        lambda X: X[:, 0] ** 2 + X[:, 1] ** 2 - radius**2,
        degree=2,
        name="phi",
    )

    result = cutcells.cut(mesh, ls)
    negative = result["phi < 0"].implicit_quadrature(order=8, mode="cut_only")
    interface = result["phi = 0"].implicit_quadrature(order=8, mode="cut_only")

    points = np.asarray(interface.points).reshape(-1, interface.tdim)
    residual = points[:, 0] ** 2 + points[:, 1] ** 2 - radius**2

    assert np.asarray(negative.weights).sum() == pytest.approx(
        math.pi * radius**2 / 4, rel=4.0e-4
    )
    assert np.asarray(interface.weights).sum() == pytest.approx(
        math.pi * radius / 2, rel=4.0e-4
    )
    assert np.max(np.abs(residual)) < 1.0e-11


def test_tetra_sphere_cut_matches_analytic_reference_and_roots():
    mesh = _single_tetra_mesh()
    radius = 0.45
    ls = cutcells.interpolate_level_set(
        mesh,
        lambda X: X[:, 0] ** 2 + X[:, 1] ** 2 + X[:, 2] ** 2 - radius**2,
        degree=2,
        name="phi",
    )

    result = cutcells.cut(mesh, ls)
    negative = result["phi < 0"].implicit_quadrature(order=8, mode="cut_only")
    interface = result["phi = 0"].implicit_quadrature(order=8, mode="cut_only")

    points = np.asarray(interface.points).reshape(-1, interface.tdim)
    residual = np.sum(points**2, axis=1) - radius**2

    assert np.asarray(negative.weights).sum() == pytest.approx(
        math.pi * radius**3 / 6, rel=2.0e-3
    )
    assert np.asarray(interface.weights).sum() == pytest.approx(
        math.pi * radius**2 / 2, rel=2.0e-3
    )
    assert np.max(np.abs(residual)) < 1.0e-11
