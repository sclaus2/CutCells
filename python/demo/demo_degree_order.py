#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""Two dials: the degree of the level set and the order of the quadrature.

The degree decides how well the level set describes the geometry, the order
how well the rule integrates over what the level set describes. This demo
turns both on one ball in a hexahedral mesh and measures, cell by cell, the
error against the exact sphere, so that errors of neighbouring cells cannot
cancel. It also prints what each degree costs: the time of `cut`, the time of
the rules and the points per cut cell.

Run it with ``python demo_degree_order.py``; ``--no-plot`` skips the plot and
``--python-level-set`` also times an analytic level set written in Python.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

import cutcells

CENTRE = (0.1, -0.05, 0.03)
RADIUS = 0.6
ORDERS = (1, 2, 4, 6, 8, 10)
DEGREES = (1, 2, 3, 4, 5, 6)


def sphere(X: np.ndarray) -> np.ndarray:
    return np.sqrt((X[0] - CENTRE[0]) ** 2 + (X[1] - CENTRE[1]) ** 2 + (X[2] - CENTRE[2]) ** 2) - RADIUS


def per_cell(rules, n_cells: int) -> np.ndarray:
    """The sum of the weights of every cell."""
    cells = np.repeat(rules.parent_map, np.diff(rules.offset))
    return np.bincount(cells, weights=rules.weights, minlength=n_cells)


def measure(mesh, phi, reference: dict, order: int) -> dict:
    """Per-cell errors, points and time of one level set and order."""
    result = cutcells.cut(mesh, phi)
    n = mesh.num_cells()
    start = time.perf_counter()
    cut_volume = result[phi < 0].quadrature(order=order, mode="cut_only")
    surface = result[phi == 0].quadrature(order=order)
    seconds = time.perf_counter() - start
    n_cut = result.num_cut_cells
    return {
        "volume": np.abs(per_cell(result[phi < 0].quadrature(order=order), n) - reference["volume"]).sum(),
        "area": np.abs(per_cell(surface, n) - reference["area"]).sum(),
        "volume points": cut_volume.weights.size / n_cut,
        "surface points": surface.weights.size / n_cut,
        "us per cell": 1e6 * seconds / n_cut,
    }


def build(n: int = 17) -> dict:
    mesh = cutcells.box_mesh(-1, -1, -1, 1, 1, 1, n, n, n, cell_type="hexahedron")
    exact = cutcells.analytic_sphere(list(CENTRE), RADIUS, name="phi")
    result = cutcells.cut(mesh, exact)
    fine = cutcells.QuadraysOptions()
    fine.points_per_segment = 14
    reference = {
        "volume": per_cell(result[exact < 0].quadrature(order=10, options=fine), mesh.num_cells()),
        "area": per_cell(result[exact == 0].quadrature(order=10, options=fine), mesh.num_cells()),
    }
    rows, setup = {}, {}
    for degree in DEGREES:
        start = time.perf_counter()
        phi = cutcells.create_level_set(mesh, sphere, degree=degree, name="phi")
        interpolated = time.perf_counter()
        cutcells.cut(mesh, phi)
        setup[f"P{degree}"] = (interpolated - start, time.perf_counter() - interpolated)
        rows[f"P{degree}"] = {order: measure(mesh, phi, reference, order) for order in ORDERS}
    start = time.perf_counter()
    cutcells.cut(mesh, exact)
    setup["exact"] = (0.0, time.perf_counter() - start)
    rows["exact"] = {order: measure(mesh, exact, reference, order) for order in ORDERS}
    return {"mesh": mesh, "rows": rows, "setup": setup}


def report(state: dict) -> None:
    rows = state["rows"]
    print(f"hexahedra: {state['mesh'].num_cells()}, errors summed over the cells, against the exact sphere")
    for key, title in (("volume", "volume error"), ("area", "area error")):
        print(f"\n{title}\n{'order':>6}" + "".join(f"{order:>9}" for order in ORDERS))
        for name, row in rows.items():
            print(f"{name:>6}" + "".join(f"{row[order][key]:9.1e}" for order in ORDERS))
    print(f"\npoints per cut cell, volume / surface\n{'order':>6}" + "".join(f"{order:>10}" for order in ORDERS))
    for name, row in rows.items():
        print(f"{name:>6}" + "".join(
            f"{row[order]['volume points']:6.0f}/{row[order]['surface points']:<3.0f}" for order in ORDERS))
    print(f"\ncost at order 6\n{'':>6}{'interpolate':>13}{'cut':>10}{'rules, us per cut cell':>24}")
    for name, row in rows.items():
        interpolate, cut = state["setup"][name]
        print(f"{name:>6}{interpolate:12.3f}s{cut:9.3f}s{row[6]['us per cell']:24.0f}")


def python_level_set(n: int = 9, order: int = 4) -> None:
    """The same sphere as an AnalyticLevelSet with Python callbacks, against the C++ one."""
    c = np.asarray(CENTRE)

    def value(x):
        return float(np.linalg.norm(x - c) - RADIUS)

    def gradient(x):
        return (x - c) / np.linalg.norm(x - c)

    def box_bounds(lower, upper):
        # phi first, then every component of the gradient: quadrays needs the
        # gradient's bounds to certify a box
        below, above = lower - c, upper - c
        nearest = np.where(below > 0, below, np.where(above < 0, above, 0.0))
        d_min = np.linalg.norm(nearest)
        d_max = np.linalg.norm(np.maximum(np.abs(below), np.abs(above)))
        bounds = [d_min - RADIUS, d_max - RADIUS]
        if d_min > 0.0:
            for i in range(3):
                q = [below[i] / d_min, below[i] / d_max, above[i] / d_min, above[i] / d_max]
                bounds += [min(q), max(q)]
        return bounds

    mesh = cutcells.box_mesh(-1, -1, -1, 1, 1, 1, n, n, n, cell_type="hexahedron")
    print(f"\nanalytic level sets, {mesh.num_cells()} hexahedra, order {order}")
    in_python = cutcells.AnalyticLevelSet(value, gradient, box_bounds, name="phi")
    in_cpp = cutcells.analytic_sphere(list(CENTRE), RADIUS, name="phi")
    for name, phi in (("C++", in_cpp), ("Python", in_python)):
        start = time.perf_counter()
        result = cutcells.cut(mesh, phi)
        rules = result[phi < 0].quadrature(order=order)
        print(f"{name:>8}: {time.perf_counter() - start:8.3f} s, volume {rules.weights.sum():.12f}")


def plot(state: dict, screenshot: Path | None = None) -> None:
    import matplotlib

    if screenshot is not None:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = state["rows"]
    # one hue from light to dark for the degree, ink for the exact sphere
    ramp = ["#b3cdf0", "#86aee6", "#5a8fdc", "#2a78d6", "#1d5ca8", "#123f75"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), sharey=True)
    for ax, key in zip(axes, ("volume", "area")):
        for name, color in zip(list(rows), ramp + ["#1e293b"]):
            e = np.maximum([rows[name][order][key] for order in ORDERS], 1e-16)
            ax.semilogy(ORDERS, e, marker="o", markersize=5, linewidth=2, color=color, label=name)
            ax.annotate(name, (ORDERS[-1], e[-1]), textcoords="offset points", xytext=(7, -3),
                        fontsize=9, color="#334155")
        ax.set_title(f"{key} error, summed over the cells", fontsize=11, color="#1e293b")
        ax.set_xlabel("quadrature order")
        ax.set_xticks(ORDERS)
        ax.grid(True, which="major", color="#e2e8f0", linewidth=0.8)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("error")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="lower center", ncol=len(labels))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    if screenshot is not None:
        fig.savefig(screenshot, dpi=150)
        plt.close(fig)
    else:
        plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--n", type=int, default=17, help="grid points per side")
    parser.add_argument("--python-level-set", action="store_true", help="also time an analytic level set in Python")
    parser.add_argument("--no-plot", action="store_true", help="skip the matplotlib window")
    parser.add_argument("--screenshot", type=Path, help="save the plot to a PNG instead")
    args = parser.parse_args()
    state = build(args.n)
    report(state)
    if args.python_level_set:
        python_level_set()
    if args.screenshot is not None or not args.no_plot:
        plot(state, args.screenshot)


if __name__ == "__main__":
    main()
