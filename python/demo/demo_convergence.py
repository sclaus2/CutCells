#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""Convergence of the area and perimeter of a disk under mesh refinement.

The disk {|x - c| < R} is cut out of triangular meshes of the square
[-1, 1]^2 with mesh sizes h = 2 / (n - 1). Three ways to describe and
integrate it are compared:

* ``P1 + lut``: a P1 level set, integrated on the straight pieces of the
  lookup tables. The pieces are polygons, so the error decays as h^2.
* ``P3 + quadrays``: a P3 level set, integrated by quadrays on the curved
  pieces. The error is the interpolation error of the level set, h^4.
* ``exact + quadrays``: the analytic circle, integrated by quadrays. Only the
  quadrature error remains, and it falls to round-off.

The demo prints the errors with their observed rates and, with matplotlib,
draws them in a log-log plot.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import cutcells

CENTER = (0.1, -0.05)
RADIUS = 0.6
AREA = np.pi * RADIUS**2
PERIMETER = 2.0 * np.pi * RADIUS


def circle(X: np.ndarray) -> np.ndarray:
    return np.sqrt((X[0] - CENTER[0]) ** 2 + (X[1] - CENTER[1]) ** 2) - RADIUS


def level_set(mesh, method: str):
    if method == "P1 + lut":
        return cutcells.create_level_set(mesh, circle, degree=1, name="phi"), "lut"
    if method == "P3 + quadrays":
        return cutcells.create_level_set(mesh, circle, degree=3, name="phi"), "quadrays"
    exact = cutcells.analytic_sphere([CENTER[0], CENTER[1], 0.0], RADIUS, name="phi")
    return exact, "quadrays"


METHODS = ("P1 + lut", "P3 + quadrays", "exact + quadrays")


def errors(sizes, order: int = 6) -> dict:
    """Area and perimeter errors per method and mesh size."""
    table = {method: {"h": [], "area": [], "perimeter": []} for method in METHODS}
    for n in sizes:
        mesh = cutcells.rectangle_mesh(-1.0, -1.0, 1.0, 1.0, n, n, cell_type="triangle")
        for method in METHODS:
            phi, backend = level_set(mesh, method)
            result = cutcells.cut(mesh, phi, backend=backend)
            area = result[phi < 0].quadrature(order=order).weights.sum()
            perimeter = result[phi == 0].quadrature(order=order).weights.sum()
            row = table[method]
            row["h"].append(2.0 / (n - 1))
            row["area"].append(abs(area - AREA))
            row["perimeter"].append(abs(perimeter - PERIMETER))
    return table


def report(table: dict) -> None:
    for quantity in ("area", "perimeter"):
        print(f"\n{quantity} error (rate)")
        print(f"{'h':>9}" + "".join(f"{m:>26}" for m in METHODS))
        h = table[METHODS[0]]["h"]
        for i in range(len(h)):
            line = f"{h[i]:9.4f}"
            for method in METHODS:
                e = table[method][quantity]
                rate = ""
                if i > 0 and e[i] > 0.0:
                    rate = f"({np.log(e[i - 1] / e[i]) / np.log(h[i - 1] / h[i]):5.2f})"
                line += f"{e[i]:18.3e} {rate:>7}"
            print(line)


def plot(table: dict, screenshot: Path | None = None) -> None:
    import matplotlib

    if screenshot is not None:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"P1 + lut": "#2a78d6", "P3 + quadrays": "#eb6834", "exact + quadrays": "#1baf7a"}
    markers = {"P1 + lut": "o", "P3 + quadrays": "s", "exact + quadrays": "^"}
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for ax, quantity in zip(axes, ("area", "perimeter")):
        for method in METHODS:
            h = np.asarray(table[method]["h"])
            e = np.maximum(np.asarray(table[method][quantity]), 1e-16)
            ax.loglog(h, e, marker=markers[method], markersize=7, linewidth=2,
                      color=colors[method], label=method)
            ax.annotate(method, (h[0], e[0]), textcoords="offset points", xytext=(8, -3),
                        fontsize=9, color="#334155")
        ax.set_title(f"{quantity} error", fontsize=11, color="#1e293b")
        ax.set_xlabel("mesh size h")
        ax.grid(True, which="major", color="#e2e8f0", linewidth=0.8)
        ax.set_xlim(right=ax.get_xlim()[1] * 2.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("absolute error")
    axes[1].legend(frameon=False, loc="lower right")
    fig.tight_layout()
    if screenshot is not None:
        fig.savefig(screenshot, dpi=150)
        plt.close(fig)
    else:
        plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sizes", type=int, nargs="+", default=[9, 17, 33, 65],
                        help="grid points per side of the meshes")
    parser.add_argument("--order", type=int, default=6, help="quadrature order")
    parser.add_argument("--no-plot", action="store_true", help="skip the matplotlib window")
    parser.add_argument("--screenshot", type=Path, help="save the plot to a PNG instead")
    args = parser.parse_args()

    table = errors(args.sizes, args.order)
    report(table)
    if args.screenshot is not None or not args.no_plot:
        plot(table, args.screenshot)


if __name__ == "__main__":
    main()
