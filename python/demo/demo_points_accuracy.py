#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""How many points does an accuracy take? Lookup tables, quadrays, plain points.

The level set is a quadratic sphere as a P2 function, so that it holds the exact
geometry on every cell and every error is the quadrature's. Three rules are
compared on one hexahedral mesh, against a rule of 14 points per segment:

* ``lut``, the straight pieces of the lookup tables on a Pk-iso-P1 template;
* ``quadrays``, the default: curved pieces whose base segments end at the folds
  of the height functions, and a check of the accuracy of every box;
* ``plain``, quadrays without both, the construction of Saye (2015).

Errors are summed over the cells, because the errors of neighbouring cells
cancel in a total over the mesh. Both are printed.

With ``--sweep file.csv``, the output of the driver ``quadrays_sweep`` in
``benchmarks/``, the demo plots the rules against Algoim's over many placements
of the sphere.

Run it with ``python demo_points_accuracy.py``; ``--no-plot`` skips the plot.
"""

from __future__ import annotations

import argparse
import csv
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

import cutcells

CENTRE = (0.0123, -0.0371, 0.0217)
RADIUS = 0.7
POINTS_PER_SEGMENT = (1, 2, 3, 4, 5, 6, 8)


def quadratic_sphere(X: np.ndarray) -> np.ndarray:
    return (X[0] - CENTRE[0]) ** 2 + (X[1] - CENTRE[1]) ** 2 + (X[2] - CENTRE[2]) ** 2 - RADIUS**2


def per_cell(rules, n_cells: int) -> np.ndarray:
    cells = np.repeat(rules.parent_map, np.diff(rules.offset))
    return np.bincount(cells, weights=rules.weights, minlength=n_cells)


def quadrays_options(q: int, plain: bool):
    options = cutcells.QuadraysOptions()
    options.points_per_segment = q
    if plain:
        options.anchor_folds = 0
        options.tolerance = -1.0
    return options


def row(rules, seconds: float, reference: np.ndarray, n_cells: int, n_cut: int) -> dict:
    cells = per_cell(rules, n_cells)
    return {
        "points": rules.weights.size / n_cut,
        "us": 1e6 * seconds / n_cut,
        "per cell": np.abs(cells - reference).sum(),
        "total": abs(cells.sum() - reference.sum()),
    }


def build(cell_type: str = "hexahedron", n: int = 17) -> dict:
    mesh = cutcells.box_mesh(-1, -1, -1, 1, 1, 1, n, n, n, cell_type=cell_type)
    phi = cutcells.create_level_set(mesh, quadratic_sphere, degree=2, name="phi")
    result = cutcells.cut(mesh, phi)
    n_cells, n_cut = mesh.num_cells(), result.num_cut_cells
    fine = quadrays_options(14, False)
    parts = {"volume": phi < 0, "surface": phi == 0}
    reference = {name: per_cell(result[part].quadrature(order=10, mode="cut_only", options=fine), n_cells)
                 for name, part in parts.items()}
    rows = defaultdict(list)
    for name, part in parts.items():
        selected = result[part]
        for method in ("quadrays", "plain"):
            for q in POINTS_PER_SEGMENT:
                options = quadrays_options(q, method == "plain")
                start = time.perf_counter()
                rules = selected.quadrature(order=2 * q, mode="cut_only", options=options) if q <= 5 else \
                    selected.quadrature(order=10, mode="cut_only", options=options)
                rows[name, method].append({"q": q} | row(rules, time.perf_counter() - start, reference[name], n_cells, n_cut))
        for k in (1, 2, 4):
            options = cutcells.LutOptions(template_order=k)
            for order in (1, 4):
                start = time.perf_counter()
                rules = selected.quadrature(order=order, mode="cut_only", backend="lut", options=options)
                rows[name, f"lut, k = {k}"].append(
                    {"q": order} | row(rules, time.perf_counter() - start, reference[name], n_cells, n_cut))
    return {"cell_type": cell_type, "n_cut": n_cut, "rows": rows}


def report(state: dict) -> None:
    print(f"{state['cell_type']}, cut cells: {state['n_cut']}; q is the points per segment (the order for lut)")
    print(f"{'part':>8}{'rule':>11}{'q':>3}{'points/cut':>12}{'us/cut':>8}{'per-cell err':>14}{'total err':>11}")
    for (part, method), rows in state["rows"].items():
        for r in rows:
            print(f"{part:>8}{method:>11}{r['q']:3d}{r['points']:12.1f}{r['us']:8.0f}{r['per cell']:14.1e}{r['total']:11.1e}")


def read_sweep(path: Path) -> dict:
    """Median and worst error, and mean points, per generator, part and q."""
    samples = defaultdict(list)
    with open(path) as f:
        for r in csv.DictReader(f):
            samples[r["generator"], r["part"], int(r["q"])].append(
                (float(r["points"]) / float(r["cut_cells"]), float(r["l1"]), float(r["total_error"]), float(r["seconds"]) / float(r["cut_cells"])))
    curves = defaultdict(list)
    for (generator, part, q), rows in sorted(samples.items()):
        a = np.array(rows)
        curves[generator, part].append({
            "q": q, "points": a[:, 0].mean(), "median": np.median(a[:, 1]), "worst": a[:, 1].max(),
            "total": np.median(a[:, 2]), "us": 1e6 * np.median(a[:, 3]),
        })
    return curves


def report_sweep(curves: dict) -> None:
    print("\nOver the placements of the sphere: per-cell error, median and worst; total, median")
    print(f"{'part':>8}{'rule':>10}{'q':>3}{'points/cut':>12}{'median':>10}{'worst':>10}{'total':>10}{'us/cut':>9}")
    for (generator, part), rows in curves.items():
        if generator in ("quadrays", "algoim"):
            for r in rows:
                print(f"{part:>8}{generator:>10}{r['q']:3d}{r['points']:12.1f}{r['median']:10.1e}{r['worst']:10.1e}{r['total']:10.1e}{r['us']:9.0f}")


def plot(curves: dict, screenshot: Path | None = None) -> None:
    import matplotlib

    if screenshot is not None:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    styles = {
        "lut:1": ("lut, k = 1", "#2a78d6", "o", ":"),
        "lut:2": ("lut, k = 2", "#2a78d6", "o", "--"),
        "lut:4": ("lut, k = 4", "#2a78d6", "o", "-"),
        "algoim": ("Algoim", "#1baf7a", "^", "-"),
        "quadrays": ("quadrays", "#eb6834", "s", "-"),
    }
    panels = (("phi < 0", "median", "volume, per-cell error"), ("phi = 0", "median", "surface, per-cell error"),
              ("phi = 0", "total", "surface, total over the mesh"))
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.3), sharey=True)
    for ax, (part, key, title) in zip(axes, panels):
        for generator, (label, color, marker, line) in styles.items():
            rows = [r for r in curves[generator, part]]
            x = np.array([r["points"] for r in rows])
            y = np.maximum([r[key] for r in rows], 1e-16)
            ax.loglog(x, y, marker=marker, markersize=6, linewidth=2, linestyle=line, color=color, label=label)
            if key == "median":
                worst = np.maximum([r["worst"] for r in rows], 1e-16)
                ax.fill_between(x, y, worst, color=color, alpha=0.15, linewidth=0)
        ax.set_title(title, fontsize=11, color="#1e293b")
        ax.set_xlabel("points per cut cell")
        ax.grid(True, which="major", color="#e2e8f0", linewidth=0.8)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("error relative to the measure")
    axes[0].legend(frameon=False, loc="lower left")
    fig.tight_layout()
    if screenshot is not None:
        fig.savefig(screenshot, dpi=150)
        plt.close(fig)
    else:
        plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cell-type", default="hexahedron", choices=("hexahedron", "tetrahedron"))
    parser.add_argument("--n", type=int, default=17, help="grid points per side")
    parser.add_argument("--sweep", type=Path, help="output of quadrays_sweep, to plot against Algoim")
    parser.add_argument("--no-plot", action="store_true", help="skip the matplotlib window")
    parser.add_argument("--screenshot", type=Path, help="save the plot to a PNG instead")
    args = parser.parse_args()
    report(build(args.cell_type, args.n))
    if args.sweep is not None:
        curves = read_sweep(args.sweep)
        report_sweep(curves)
        if args.screenshot is not None or not args.no_plot:
            plot(curves, args.screenshot)


if __name__ == "__main__":
    main()
