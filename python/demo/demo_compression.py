#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""Compression: the same integrals of a polynomial space on far fewer points.

`compress_rules` replaces every rule by a subset of its points with new positive
weights, which integrates a polynomial space exactly as the rule did. A space of
M functions needs at most M points, whatever the rule had. The demo measures

* how many points remain, cell shape by cell shape, on volumes and surfaces;
* that mass matrices of Pk and Qk functions do not change, and how much
  assembling them saves against what compressing costs;
* what happens beyond the space: a coefficient that is not polynomial, and the
  margin of degree that protects it.

Run it with ``python demo_compression.py``; ``--no-plot`` skips the PyVista window.
"""

from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

import numpy as np

import cutcells

CENTRE = (0.1, -0.05, 0.03)
RADIUS = 0.6
# Total degree on simplices, the degree in each coordinate on the other shapes
SPACES = {"triangle": "total", "quadrilateral": "tensor", "tetrahedron": "total",
          "hexahedron": "tensor", "prism": "tensor", "pyramid": "tensor"}


def round_shape(X: np.ndarray) -> np.ndarray:
    return np.sqrt(sum((X[d] - CENTRE[d]) ** 2 for d in range(X.shape[0]))) - RADIUS


def cut_mesh(shape: str, n: int):
    if shape in ("triangle", "quadrilateral"):
        mesh = cutcells.rectangle_mesh(-1, -1, 1, 1, n, n, cell_type=shape)
    else:
        mesh = cutcells.box_mesh(-1, -1, -1, 1, 1, 1, n, n, n, cell_type=shape)
    phi = cutcells.create_level_set(mesh, round_shape, degree=2, name="phi")
    return mesh, phi, cutcells.cut(mesh, phi)


def n_moments(space: str, tdim: int, degree: int) -> int:
    return (degree + 1) ** tdim if space == "tensor" else math.comb(degree + tdim, tdim)


def monomials(points: np.ndarray, degree: int, space: str) -> np.ndarray:
    exponents = [e for e in np.ndindex(*(points.shape[1] * [degree + 1])) if space == "tensor" or sum(e) <= degree]
    return np.stack([np.prod(points ** np.array(e), axis=1) for e in exponents], axis=1)


def assemble(rules, k: int, space: str, coefficient=None) -> list[np.ndarray]:
    """The mass matrix of every cut cell, for the monomials of degree k."""
    points = rules.points.reshape(-1, rules.tdim)
    matrices = []
    for a, b in zip(rules.offset[:-1], rules.offset[1:]):
        basis = monomials(points[a:b], k, space)
        w = rules.weights[a:b]
        if coefficient is not None:
            w = w * coefficient(points[a:b])
        matrices.append((basis.T * w) @ basis)
    return matrices


def difference(full: list[np.ndarray], small: list[np.ndarray]) -> float:
    """The largest entry of a difference relative to the largest entry of the matrix, over the cells."""
    return max(np.abs(a - b).max() / np.abs(a).max() for a, b in zip(full, small))


def timed(function):
    start = time.perf_counter()
    value = function()
    return value, time.perf_counter() - start


def point_counts(n: int, order: int, degree: int) -> list[dict]:
    rows = []
    for shape, space in SPACES.items():
        mesh, phi, result = cut_mesh(shape, n)
        for name, part, mode in (("volume", phi < 0, "cut_only"), ("surface", phi == 0, "full")):
            rules = result[part].quadrature(order=order, mode=mode)
            (_, stats), seconds = timed(lambda: cutcells.compress_rules(rules, degree=degree, space=space))
            rows.append({
                "shape": shape, "part": name, "space": space,
                "moments": n_moments(space, mesh.tdim, degree),
                "before": stats.points_before / stats.n_rules,
                "after": stats.points_after / stats.n_rules,
                "ms": 1e3 * seconds / stats.n_rules,
                "residual": stats.max_residual,
            })
    return rows


def mass_matrices(n: int, order: int) -> list[dict]:
    rows = []
    for shape in ("tetrahedron", "hexahedron"):
        space = SPACES[shape]
        _, phi, result = cut_mesh(shape, n)
        rules = result[phi < 0].quadrature(order=order, mode="cut_only")
        for k in (1, 2, 3):
            (small, stats), compress = timed(lambda: cutcells.compress_rules(rules, degree=2 * k, space=space))
            full_matrices, full_time = timed(lambda: assemble(rules, k, space))
            small_matrices, small_time = timed(lambda: assemble(small, k, space))
            saved = full_time - small_time
            rows.append({
                "shape": shape, "k": k, "degree": 2 * k,
                "before": stats.points_before / stats.n_rules,
                "after": stats.points_after / stats.n_rules,
                "difference": difference(full_matrices, small_matrices),
                "compress ms": 1e3 * compress / stats.n_rules,
                "full ms": 1e3 * full_time / stats.n_rules,
                "small ms": 1e3 * small_time / stats.n_rules,
                "pays after": compress / saved if saved > 0 else math.inf,
            })
    return rows


def margins(n: int, order: int, k: int = 2) -> list[dict]:
    """A smooth coefficient, which no finite space holds, against a margin of degree."""
    coefficients = {"exp(x + y + z)": lambda p: np.exp(p.sum(axis=1)),
                    "1 / (1 + 4 r^2)": lambda p: 1.0 / (1.0 + 4.0 * np.sum(p**2, axis=1))}
    rows = []
    for shape in ("tetrahedron", "hexahedron"):
        space = SPACES[shape]
        _, phi, result = cut_mesh(shape, n)
        rules = result[phi < 0].quadrature(order=order, mode="cut_only")
        full = {name: assemble(rules, k, space, f) for name, f in coefficients.items()}
        for margin in (0, 2, 4):
            small, stats = cutcells.compress_rules(rules, degree=2 * k + margin, space=space)
            row = {"shape": shape, "degree": 2 * k + margin, "after": stats.points_after / stats.n_rules}
            for name, f in coefficients.items():
                row[name] = difference(full[name], assemble(small, k, space, f))
            rows.append(row)
    return rows


def build(n: int = 9, order: int = 8) -> dict:
    return {"n": n, "order": order, "counts": point_counts(n, order, 4), "mass": mass_matrices(n, order),
            "margins": margins(n, order)}


def report(state: dict) -> None:
    print(f"Ball, P2 level set, {state['n']} grid points per side, rules of order {state['order']}")
    print("\nCompressed onto degree 4, points per rule")
    print(f"{'shape':>14}{'part':>9}{'space':>8}{'moments':>9}{'before':>8}{'after':>8}{'ms/rule':>9}{'residual':>10}")
    for r in state["counts"]:
        print(f"{r['shape']:>14}{r['part']:>9}{r['space']:>8}{r['moments']:9d}{r['before']:8.0f}{r['after']:8.0f}"
              f"{r['ms']:9.2f}{r['residual']:10.1e}")
    print("\nMass matrices of degree k on the cut cells, compressed onto degree 2 k, times per rule")
    print(f"{'shape':>14}{'k':>3}{'points':>16}{'difference':>12}{'compress':>11}{'assemble, full':>16}{'compressed':>12}{'pays after':>12}")
    for r in state["mass"]:
        points = f"{r['before']:.0f} -> {r['after']:.0f}"
        print(f"{r['shape']:>14}{r['k']:3d}{points:>16}{r['difference']:12.1e}{r['compress ms']:9.2f}ms"
              f"{r['full ms']:14.2f}ms{r['small ms']:10.2f}ms{r['pays after']:12.1f}")
    print("\nMass matrices of degree 2 with a coefficient, against the full rule")
    names = [key for key in state["margins"][0] if key not in ("shape", "degree", "after")]
    print(f"{'shape':>14}{'degree':>8}{'points':>8}" + "".join(f"{name:>18}" for name in names))
    for r in state["margins"]:
        print(f"{r['shape']:>14}{r['degree']:8d}{r['after']:8.0f}" + "".join(f"{r[name]:18.1e}" for name in names))


def reference_cell(shape: str):
    """The reference cell the rules live on, for drawing."""
    import pyvista as pv

    if shape == "tetrahedron":
        vertices, vtk = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], pv.CellType.TETRA
    else:
        vertices, vtk = [[x, y, z] for z in (0, 1) for x, y in ((0, 0), (1, 0), (1, 1), (0, 1))], pv.CellType.HEXAHEDRON
    return pv.UnstructuredGrid({vtk: np.arange(len(vertices)).reshape(1, -1)}, np.asarray(vertices, dtype=float))


def plot(state: dict, screenshot: Path | None = None) -> None:
    import pyvista as pv

    plotter = pv.Plotter(shape=(1, 4), off_screen=screenshot is not None, window_size=(1600, 460), border=False)
    column = 0
    for shape, grid_n in (("tetrahedron", 5), ("hexahedron", 5)):
        mesh, phi, result = cut_mesh(shape, grid_n)
        rules = result[phi < 0].quadrature(order=6, mode="cut_only")
        small, _ = cutcells.compress_rules(rules, degree=4, space=SPACES[shape])
        # the cell whose rule has the median number of points
        sizes = np.diff(rules.offset)
        rule = int(np.argsort(sizes)[sizes.size // 2])
        cell = reference_cell(shape)
        for name, r in (("full rule", rules), ("compressed", small)):
            plotter.subplot(0, column)
            a, b = r.offset[rule], r.offset[rule + 1]
            reference = r.points.reshape(-1, 3)[a:b]
            cloud = pv.PolyData(reference)
            cloud["weight"] = r.weights[a:b]
            plotter.add_mesh(cell, style="wireframe", color="#64748b", line_width=2)
            plotter.add_points(cloud, scalars="weight", cmap="viridis", point_size=7 if name == "full rule" else 14,
                               render_points_as_spheres=True, show_scalar_bar=False)
            plotter.add_text(f"{shape}, {name}: {b - a} points", font_size=10, color="#1e293b", position="upper_left")
            plotter.view_isometric()
            column += 1
    if screenshot is not None:
        plotter.screenshot(str(screenshot))
        plotter.close()
    else:
        plotter.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--n", type=int, default=9, help="grid points per side")
    parser.add_argument("--order", type=int, default=8, help="order of the rules that are compressed")
    parser.add_argument("--no-plot", action="store_true", help="skip the PyVista window")
    parser.add_argument("--screenshot", type=Path, help="render the view to a PNG instead")
    args = parser.parse_args()
    state = build(args.n, args.order)
    report(state)
    if args.screenshot is not None or not args.no_plot:
        plot(state, args.screenshot)


if __name__ == "__main__":
    main()
