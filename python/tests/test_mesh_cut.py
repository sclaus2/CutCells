# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""cutcells.cut on a mesh: quadrays by default, the lookup tables with
backend="lut" and the keywords of the former cut(), and the retired names
cutcells.part.cut, cutcells.ho_cut, HOCutResult and HOMeshPart."""

import numpy as np
import pytest

import cutcells


def single_triangle_mesh():
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float64)
    return cutcells.MeshView(coords, np.array([0, 1, 2], dtype=np.int32), np.array([0, 3], dtype=np.int32),
                             np.array([5], dtype=np.int32), tdim=2)


def single_tetra_mesh():
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)
    return cutcells.MeshView(coords, np.array([0, 1, 2, 3], dtype=np.int32), np.array([0, 4], dtype=np.int32),
                             np.array([10], dtype=np.int32), tdim=3)


def test_names_and_defaults():
    mesh = single_triangle_mesh()
    ls = cutcells.create_level_set(mesh, lambda x: x[0] + x[1] - 0.3, degree=1, name="phi")
    for result, backend in ((cutcells.cut(mesh, ls), "quadrays"), (cutcells.cut(mesh, ls, backend="lut"), "lut")):
        assert isinstance(result, cutcells.part.CutResult) and result.backend == backend
        if backend == "lut":
            assert isinstance(result.options, cutcells.LutOptions)
        else:
            assert result.options is None
        assert result.num_cut_cells == 1 and result.num_level_sets == 1
        np.testing.assert_array_equal(np.asarray(result.parent_cell_ids), np.array([0]))
        assert np.asarray(result.cell_domains).shape == (1, 1)
        negative, interface = result["phi < 0"], result["phi = 0"]
        assert isinstance(negative, cutcells.part.MeshPart) and negative.backend == backend
        np.testing.assert_array_equal(np.asarray(negative.cut_cell_ids), np.array([0]))
        np.testing.assert_array_equal(np.asarray(interface.cut_cell_ids), np.array([0]))
        assert np.sum(negative.quadrature(order=1).weights) == pytest.approx(0.5 * 0.3**2, rel=1e-14)
        assert np.sum(interface.quadrature(order=1).weights) == pytest.approx(0.3 * np.sqrt(2.0), rel=1e-14)


def same_rules(a, b):
    for name in ("points", "weights", "offset", "parent_map"):
        np.testing.assert_array_equal(np.asarray(getattr(a, name)), np.asarray(getattr(b, name)))


def test_retired_cut_functions_warn_and_forward():
    mesh = single_triangle_mesh()
    ls = cutcells.create_level_set(mesh, lambda x: (x[0] - 0.1) ** 2 + x[1] ** 2 - 0.25, degree=2, name="phi")
    new = cutcells.cut(mesh, ls)
    for name in ("cut", "cut_float64"):
        with pytest.warns(DeprecationWarning, match=f"cutcells.part.{name} is deprecated; use cutcells.cut"):
            old = getattr(cutcells.part, name)(mesh, ls)
        assert old.backend == new.backend == "quadrays"
        np.testing.assert_array_equal(np.asarray(old.cell_domains), np.asarray(new.cell_domains))
        for selection in ("phi < 0", "phi = 0"):
            same_rules(old[selection].quadrature(order=3), new[selection].quadrature(order=3))
    with pytest.warns(DeprecationWarning, match="cutcells.ho_cut is deprecated; use cutcells.cut"):
        old = cutcells.ho_cut(mesh, ls)
    assert old.backend == "quadrays"
    same_rules(old["phi < 0"].quadrature(order=3), new["phi < 0"].quadrature(order=3))
    with pytest.warns(DeprecationWarning, match="cutcells.part.cut_float32 is deprecated"):
        with pytest.raises(TypeError):
            cutcells.part.cut_float32(mesh, ls)


def test_retired_type_names_warn():
    for old, new in [("HOCutResult", "CutResult"), ("HOMeshPart", "MeshPart")]:
        for suffix in ("", "_float32", "_float64"):
            with pytest.warns(DeprecationWarning, match=f"cutcells.{old}{suffix} is deprecated"):
                assert getattr(cutcells, old + suffix) is getattr(cutcells.part, new + suffix)
            assert old + suffix not in dir(cutcells)
    with pytest.warns(DeprecationWarning):
        from cutcells import HOMeshPart  # noqa: F401
    with pytest.raises(AttributeError, match="no_such_name"):
        cutcells.no_such_name


def test_backend_per_call_and_per_result():
    mesh = single_tetra_mesh()
    # a ball strictly inside the tetrahedron
    ls = cutcells.create_level_set(
        mesh, lambda x: (x[0] - 0.22) ** 2 + (x[1] - 0.22) ** 2 + (x[2] - 0.22) ** 2 - 0.12**2, degree=2, name="phi")
    result = cutcells.cut(mesh, ls, backend="lut")
    part = result["phi < 0"]
    exact = 4.0 / 3.0 * np.pi * 0.12**3
    straight = np.sum(part.quadrature(order=2).weights)
    assert np.sum(part.quadrature(order=2, backend="straight").weights) == straight
    curved = np.sum(part.quadrature(order=5, backend="quadrays").weights)
    assert curved == pytest.approx(exact, rel=1e-5)
    assert abs(straight / exact - 1) > 1e-3
    # a result's backend applies to the parts selected afterwards
    result.backend = "quadrays"
    assert result.options is None
    assert np.sum(result["phi < 0"].quadrature(order=5).weights) == curved
    assert type(result["phi < 0"].visualization_mesh()).__name__.startswith("QuadraysLeafMesh")
    with pytest.raises(ValueError, match="benchmarks"):
        part.quadrature(backend="algoim")
    with pytest.raises(TypeError, match="LutOptions"):
        cutcells.cut(mesh, ls, backend="lut", options=cutcells.QuadraysOptions())
    with pytest.raises(TypeError, match="QuadraysOptions"):
        cutcells.cut(mesh, ls, options=cutcells.LutOptions())


def test_quadratic_interior_intersection_is_found():
    """A ball inside a tetrahedron whose vertices are all outside it: the
    bounds of the P2 level set find the cut, and quadrays integrates it."""
    mesh = single_tetra_mesh()

    def ball(x):
        return (x[0] - 0.2) ** 2 + (x[1] - 0.2) ** 2 + (x[2] - 0.2) ** 2 - 0.09

    ls = cutcells.create_level_set(mesh, ball, degree=2, name="phi")
    assert np.all(np.array([ball(v) for v in np.asarray(mesh.coordinates)]) > 0.0)
    assert np.min(np.asarray(cutcells.make_cell_level_set(ls, 0).bernstein_coeffs)) < 0.0
    result = cutcells.cut(mesh, ls)
    assert result.num_cut_cells == 1
    np.testing.assert_array_equal(np.asarray(result["phi = 0"].cut_cell_ids), np.array([0]))
    assert np.sum(result["phi < 0"].quadrature(order=4, backend="quadrays").weights) > 0.0


def test_legacy_keywords():
    mesh = single_tetra_mesh()
    ls = cutcells.create_level_set(mesh, lambda x: x[0] + x[1] - 0.6, degree=1, name="phi")
    result = cutcells.cut(mesh, ls, backend="lut", triangulate=True, triangulation="midpoint",
                          cut_approximation="iso_p1", cut_approximation_order=2)
    options = result.options
    assert (options.template_order, options.triangulate, options.triangulation) == (2, True, "midpoint")
    cells = result["phi < 0"].visualization_mesh(mode="cut_only")
    assert set(np.unique(cells.vtk_types)) == {10}
    total = np.sum(result["phi < 0"].quadrature(order=1).weights)
    plain = cutcells.cut(mesh, ls, backend="lut")
    assert total == pytest.approx(np.sum(plain["phi < 0"].quadrature(order=1).weights), rel=1e-14)
    phi = cutcells.analytic_sphere([0.2, 0.2, 0.2], 0.25)
    named = cutcells.cut(mesh, phi, name="ball", backend="lut", degree=3)
    assert named.level_set_names == ["ball"] and named.options.template_order == 3
    with pytest.raises(ValueError, match="name or names"):
        cutcells.cut(mesh, phi, name="a", names=["b"])


@pytest.mark.parametrize("keywords", [dict(triangulate=True), dict(triangulation="midpoint"),
                                      dict(cut_approximation="iso_p1"), dict(cut_approximation_order=2),
                                      dict(degree=2)])
def test_legacy_keywords_need_the_lookup_tables(keywords):
    mesh = single_tetra_mesh()
    ls = cutcells.create_level_set(mesh, lambda x: x[0] + x[1] - 0.6, degree=1, name="phi")
    for backend in ({}, dict(backend="quadrays")):
        with pytest.raises(ValueError, match='backend="lut"'):
            cutcells.cut(mesh, ls, **backend, **keywords)
    assert cutcells.cut(mesh, ls, backend="lut", **keywords).backend == "lut"
