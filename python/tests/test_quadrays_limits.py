# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""The per-cell limits of quadrays, their warning, and the debug path that
writes the cells that take many points."""

import glob
import json
import os

import numpy as np
import pytest

import cutcells as cc

R = 0.5
C1, C2 = (-0.21, 0.07, 0.04), (0.24, 0.03, -0.02)


def lens(n=2, distance=True):
    mesh = cc.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n, n, n, cell_type="hexahedron")
    a = cc.analytic_sphere(list(C1), R, signed_distance=distance, name="a")
    b = cc.analytic_sphere(list(C2), R, signed_distance=distance, name="b")
    return cc.cut(mesh, [a, b]), a, b


def worst_cell(rules):
    return int(np.diff(np.asarray(rules.offset)).max())


def test_point_limit_bounds_cells_and_warns(capfd):
    result, a, b = lens()
    options = cc.QuadraysOptions()
    options.max_points = 20000
    rules = result[(a < 0) & (b < 0)].quadrature(order=8, options=options)
    err = capfd.readouterr().err
    assert worst_cell(rules) <= 20000
    assert "reached a limit" in err and "max_points = 20000" in err
    assert np.asarray(rules.weights).min() > 0
    # less accurate, but still the lens
    exact = np.pi * (4 * R + np.linalg.norm(np.subtract(C1, C2))) * (2 * R - np.linalg.norm(np.subtract(C1, C2))) ** 2 / 12
    assert abs(np.asarray(rules.weights).sum() - exact) < 1e-3


def test_no_warning_below_the_limits(capfd):
    result, a, b = lens(n=5, distance=False)
    rules = result[(a < 0) & (b < 0)].quadrature(order=4)
    assert "warning" not in capfd.readouterr().err
    assert worst_cell(rules) > 0


def test_bisection_limit_warns(capfd):
    result, a, b = lens()
    options = cc.QuadraysOptions()
    options.max_bisections = 8
    result[(a < 0) & (b < 0)].quadrature(order=8, options=options)
    assert "max_bisections = 8" in capfd.readouterr().err


def test_nothing_is_written_by_default(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    result, a, b = lens()
    options = cc.QuadraysOptions()
    options.max_points = 20000
    result[(a < 0) & (b < 0)].quadrature(order=8, options=options)
    assert os.listdir(tmp_path) == []


def test_debug_dir_writes_the_cell(tmp_path):
    result, a, b = lens(n=2, distance=False)
    options = cc.QuadraysOptions()
    options.debug_dir = str(tmp_path / "cells")
    options.debug_points = 1000
    rules = result[(a < 0) & (b < 0)].quadrature(order=8, options=options)
    reports = glob.glob(str(tmp_path / "cells" / "*.json"))
    assert 1 <= len(reports) <= 16

    report = json.loads(open(reports[0]).read())
    assert report["selection"] == "a < 0 and b < 0"
    assert report["level_sets"] == ["a", "b"]
    assert report["q"] == 6 and report["order"] == 8
    assert report["n_points"] >= 1000
    assert report["n_points"] in np.diff(np.asarray(rules.offset)).tolist()
    assert len(report["values"]) == 2
    assert len(report["values"][0]) == report["grid"] ** 3


def test_debug_report_regions_partition_the_points(tmp_path):
    """Boxes in rotated frames integrate the region they replace: every point
    lies in its top box, and the top boxes do not overlap."""
    debug = pytest.importorskip("cutcells.quadrays_debug")
    result, a, b = lens(n=2, distance=True)
    options = cc.QuadraysOptions()
    options.debug_dir = str(tmp_path)
    options.max_points = -1
    result[(a < 0) & (b < 0)].quadrature(order=8, options=options)
    report = debug.read_report(glob.glob(str(tmp_path / "*.json"))[0])
    boxes, leaves = report["boxes"], debug.top_leaves(report)
    assert np.any(boxes["fate"][leaves] == debug.FATES.index("rotated"))
    assert boxes["points"][leaves].sum() == report["n_points"]

    lo, hi = boxes["lo"][leaves], boxes["hi"][leaves]
    overlap = 0.0
    for i in range(len(leaves)):
        extent = np.minimum(hi[i], hi[i + 1 :]) - np.maximum(lo[i], lo[i + 1 :])
        overlap += np.prod(np.clip(extent, 0, None), axis=1).sum()
    assert overlap == 0.0

    start, stride, points = 0, report["point_stride"], report["points"]
    for i in leaves:
        n = boxes["points"][i]
        kept = np.arange(start, start + n)
        kept = kept[kept % stride == 0] // stride
        kept = kept[kept < len(points)]
        assert np.all(points[kept] >= boxes["lo"][i] - 1e-12) and np.all(points[kept] <= boxes["hi"][i] + 1e-12)
        start += n


def test_debug_plot(tmp_path):
    pytest.importorskip("matplotlib")
    debug = pytest.importorskip("cutcells.quadrays_debug")
    result, a, b = lens(n=2, distance=False)
    options = cc.QuadraysOptions()
    options.debug_dir = str(tmp_path)
    result[(a < 0) & (b < 0)].quadrature(order=8, options=options)
    path = glob.glob(str(tmp_path / "*.json"))[0]
    image = tmp_path / "cell.png"
    debug.plot_report(path, image)
    assert image.stat().st_size > 0
    assert "a < 0 and b < 0" in debug.summary(debug.read_report(path))


def test_debug_report_in_2d(tmp_path):
    debug = pytest.importorskip("cutcells.quadrays_debug")
    mesh = cc.rectangle_mesh(-1.0, -1.0, 1.0, 1.0, 2, 2, cell_type="quadrilateral")
    a = cc.analytic_sphere([-0.21, 0.07, 0.0], R, name="a")
    b = cc.analytic_sphere([0.24, 0.03, 0.0], R, name="b")
    options = cc.QuadraysOptions()
    options.debug_dir = str(tmp_path)
    options.debug_points = 1
    cc.cut(mesh, [a, b])[(a < 0) & (b < 0)].quadrature(order=4, options=options)
    reports = glob.glob(str(tmp_path / "*.json"))
    assert reports
    report = debug.read_report(reports[0])
    assert report["tdim"] == 2 and report["values"][0].ndim == 2
    pytest.importorskip("matplotlib")
    debug.plot_report(report, tmp_path / "cell.png")
    assert (tmp_path / "cell.png").stat().st_size > 0
