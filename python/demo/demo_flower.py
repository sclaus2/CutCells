#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""A flower cut out of a quadrilateral mesh: a non-convex curve of high curvature.

The curve is r(theta) = R (1 + a cos(k theta)) in polar coordinates, the zero
set of

    phi(x) = |x| - R (1 + a cos(k theta)).

Its area and its polar moment have closed forms,

    |Omega| = pi R^2 (1 + a^2 / 2),   int_Omega |x|^2 dx = int_0^2pi r^4 / 4 dtheta,

and its length is the integral of sqrt(r^2 + r'^2) over theta, which the
trapezoidal rule gives to round-off for this periodic function. The demo
integrates the three for level sets of several degrees on one mesh, with
the lookup tables' straight pieces for comparison.

Run it with ``python demo_flower.py``; ``--no-plot`` skips the PyVista window.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import cutcells

RADIUS = 0.6
AMPLITUDE = 0.3
PETALS = 5


def flower(X: np.ndarray) -> np.ndarray:
    theta = np.arctan2(X[1], X[0])
    return np.hypot(X[0], X[1]) - RADIUS * (1.0 + AMPLITUDE * np.cos(PETALS * theta))


def exact_values() -> dict:
    """Area, length and polar moment; the sums converge geometrically."""
    theta = np.linspace(0.0, 2.0 * np.pi, 4096, endpoint=False)
    r = RADIUS * (1.0 + AMPLITUDE * np.cos(PETALS * theta))
    dr = -RADIUS * AMPLITUDE * PETALS * np.sin(PETALS * theta)
    h = 2.0 * np.pi / theta.size
    return {
        "area": np.pi * RADIUS**2 * (1.0 + 0.5 * AMPLITUDE**2),
        "length": h * np.sum(np.hypot(r, dr)),
        "int |x|^2": h * np.sum(r**4) / 4.0,
    }


def physical_points(mesh, rules) -> np.ndarray:
    """Map the reference points of rules on a mesh of parallelograms to physical space."""
    reference = rules.points.reshape(-1, 2)
    cells = np.repeat(rules.parent_map, np.diff(rules.offset))
    nodes = mesh.connectivity.reshape(-1, 4)[cells]  # VTK order: around the cell
    x = mesh.coordinates[nodes]
    return x[:, 0] + reference[:, :1] * (x[:, 1] - x[:, 0]) + reference[:, 1:] * (x[:, 3] - x[:, 0])


def measures(mesh, level_set, backend: str, order: int) -> dict:
    result = cutcells.cut(mesh, level_set, backend=backend)
    inside = result[level_set < 0].quadrature(order=order, mode="full")
    x = physical_points(mesh, inside)
    return {
        "area": inside.weights.sum(),
        "length": result[level_set == 0].quadrature(order=order).weights.sum(),
        "int |x|^2": np.dot(inside.weights, np.sum(x**2, axis=1)),
    }


def build(n: int = 41, degrees: tuple[int, ...] = (1, 2, 3, 4), order: int = 6) -> dict:
    mesh = cutcells.rectangle_mesh(-1.0, -1.0, 1.0, 1.0, n, n, cell_type="quadrilateral")
    rows = {}
    for degree in degrees:
        phi = cutcells.create_level_set(mesh, flower, degree=degree, name="phi")
        rows[f"P{degree} + quadrays"] = measures(mesh, phi, "quadrays", order)
    phi = cutcells.create_level_set(mesh, flower, degree=3, name="phi")
    rows["P3 + lut"] = measures(mesh, phi, "lut", order)
    return {"mesh": mesh, "phi": phi, "result": cutcells.cut(mesh, phi), "rows": rows, "exact": exact_values()}


def report(state: dict) -> None:
    mesh, exact = state["mesh"], state["exact"]
    print(f"quadrilaterals: {mesh.num_cells()}, cut cells: {state['result'].num_cut_cells}")
    print("exact: " + ", ".join(f"{key} = {value:.12f}" for key, value in exact.items()))
    print(f"{'level set':>16}" + "".join(f"{key + ' error':>18}" for key in exact))
    for name, row in state["rows"].items():
        print(f"{name:>16}" + "".join(f"{abs(row[key] - exact[key]):18.2e}" for key in exact))


def write(state: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    phi, result = state["phi"], state["result"]
    result[phi < 0].write_vtu(str(output_dir / "flower.vtu"), mode="full", degree=3)
    result[phi == 0].write_vtu(str(output_dir / "curve.vtu"), mode="cut_only", degree=3)
    print(f"wrote flower.vtu and curve.vtu (Lagrange cells of degree 3) to {output_dir}")


def plot(state: dict, screenshot: Path | None = None) -> None:
    import pyvista as pv

    phi, result = state["phi"], state["result"]
    background = cutcells.to_pyvista(state["mesh"])
    background.cell_data["class"] = result.domains[0]
    curve = cutcells.to_pyvista(result[phi == 0].visualization_mesh(mode="cut_only", degree=3))
    curve = curve.extract_surface(nonlinear_subdivision=4, algorithm="dataset_surface")
    curve.points[:, 2] = 0.01

    plotter = pv.Plotter(off_screen=screenshot is not None, window_size=(900, 900))
    plotter.add_mesh(background, scalars="class", show_edges=True, edge_color="#94a3b8",
                     cmap=["#68b984", "#f59e42", "#e2e8f0"], clim=[-0.5, 2.5], show_scalar_bar=False,
                     lighting=False)
    plotter.add_mesh(curve, color="#0f766e", line_width=4)
    plotter.view_xy()
    plotter.enable_parallel_projection()
    plotter.camera.zoom(1.25)
    if screenshot is not None:
        plotter.screenshot(str(screenshot))
        plotter.close()
    else:
        plotter.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--n", type=int, default=41, help="grid points per side")
    parser.add_argument("--order", type=int, default=6, help="quadrature order")
    parser.add_argument("--output-dir", type=Path, default=Path("demo_flower_output"))
    parser.add_argument("--no-plot", action="store_true", help="skip the PyVista window")
    parser.add_argument("--screenshot", type=Path, help="render the view to a PNG instead")
    args = parser.parse_args()

    state = build(args.n, order=args.order)
    report(state)
    write(state, args.output_dir)
    if args.screenshot is not None or not args.no_plot:
        plot(state, args.screenshot)


if __name__ == "__main__":
    main()
