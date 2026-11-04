#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""A torus in a hexahedral mesh: a level set of your own, as Python functions.

The torus with radii R > r about the z axis is the zero set of

    phi(x) = sqrt((rho - R)^2 + z^2) - r,   rho = sqrt(x^2 + y^2),

with the exact volume and area

    |T| = 2 pi^2 R r^2,   |dT| = 4 pi^2 R r.

The demo interpolates phi with several degrees. With ``--analytic`` it also
gives phi as an analytic level set written here in Python: its value, its
gradient and bounds of both over a box. The bounds are what lets quadrays
decide, for a whole box at once, that a cell is not cut and which direction
the surface is a height function along.

Python callbacks are slow, and this run takes minutes; for production use,
compile them and pass a capsule (``cutcells.AnalyticLevelSet.from_capsule``).

Run it with ``python demo_torus.py``; ``--no-plot`` skips the PyVista window.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

import cutcells

MAJOR = 0.65
MINOR = 0.25


def torus(X: np.ndarray) -> np.ndarray:
    return np.hypot(np.hypot(X[0], X[1]) - MAJOR, X[2]) - MINOR


def value(x: np.ndarray) -> float:
    return float(np.hypot(np.hypot(x[0], x[1]) - MAJOR, x[2]) - MINOR)


def gradient(x: np.ndarray) -> np.ndarray:
    rho = np.hypot(x[0], x[1])
    d = np.hypot(rho - MAJOR, x[2])
    return np.array([(rho - MAJOR) / d * x[0] / rho, (rho - MAJOR) / d * x[1] / rho, x[2] / d])


def absolute_range(lo: float, hi: float) -> tuple[float, float]:
    """The range of |t| for t in [lo, hi]."""
    if lo > 0.0:
        return lo, hi
    if hi < 0.0:
        return -hi, -lo
    return 0.0, max(-lo, hi)


def quotient_range(a: tuple[float, float], lo: float, hi: float) -> tuple[float, float]:
    """The range of a / t for a in the interval ``a`` and t in [lo, hi], lo > 0."""
    q = (a[0] / lo, a[0] / hi, a[1] / lo, a[1] / hi)
    return min(q), max(q)


def product_range(a: tuple[float, float], b: tuple[float, float]) -> tuple[float, float]:
    p = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return min(p), max(p)


def box_bounds(lo: np.ndarray, hi: np.ndarray) -> list[float]:
    """Interval bounds of phi and, where they exist, of its gradient over a box."""
    ax, ay, az = (absolute_range(lo[i], hi[i]) for i in range(3))
    rho = (np.hypot(ax[0], ay[0]), np.hypot(ax[1], ay[1]))
    t = (rho[0] - MAJOR, rho[1] - MAJOR)
    at = absolute_range(*t)
    d = (np.hypot(at[0], az[0]), np.hypot(at[1], az[1]))
    bounds = [d[0] - MINOR, d[1] - MINOR]
    if d[0] <= 0.0 or rho[0] <= 0.0:
        return bounds  # the gradient has no bound on the circle rho = R, z = 0, nor on the axis
    radial = quotient_range(t, d[0], d[1])
    for i in range(2):
        bounds += product_range(radial, quotient_range((lo[i], hi[i]), rho[0], rho[1]))
    bounds += quotient_range((lo[2], hi[2]), d[0], d[1])
    return bounds


def build(n: int = 17, degrees: tuple[int, ...] = (2, 3, 4), order: int = 4, analytic: bool = False) -> dict:
    mesh = cutcells.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n, n, n, cell_type="hexahedron")
    state = {"mesh": mesh, "order": order, "parts": {}}
    level_sets = {f"P{k}": cutcells.create_level_set(mesh, torus, degree=k, name="phi") for k in degrees}
    if analytic:
        level_sets["analytic"] = cutcells.AnalyticLevelSet(value, gradient, box_bounds, name="phi")
    for key, phi in level_sets.items():
        start = time.perf_counter()
        result = cutcells.cut(mesh, phi)
        inside, surface = result[phi < 0], result[phi == 0]
        volume = inside.quadrature(order=order, mode="full").weights.sum()
        area = surface.quadrature(order=order).weights.sum()
        state["parts"][key] = {"phi": phi, "result": result, "inside": inside, "surface": surface,
                               "volume": volume, "area": area, "seconds": time.perf_counter() - start}
    return state


def report(state: dict) -> None:
    volume = 2.0 * np.pi**2 * MAJOR * MINOR**2
    area = 4.0 * np.pi**2 * MAJOR * MINOR
    first = next(iter(state["parts"].values()))
    print(f"hexahedra: {state['mesh'].num_cells()}, cut cells: {first['result'].num_cut_cells}")
    print(f"{'level set':>10}  {'volume error':>13}  {'area error':>11}  {'seconds':>8}")
    for key, part in state["parts"].items():
        print(f"{key:>10}  {abs(part['volume'] - volume):13.2e}  {abs(part['area'] - area):11.2e}"
              f"  {part['seconds']:8.2f}")


def write(state: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    part = list(state["parts"].values())[-1]
    part["inside"].write_vtu(str(output_dir / "torus.vtu"), mode="cut_only", degree=3)
    part["surface"].write_vtu(str(output_dir / "surface.vtu"), mode="cut_only", degree=3)
    print(f"wrote torus.vtu and surface.vtu (Lagrange cells of degree 3) to {output_dir}")


def plot(state: dict, screenshot: Path | None = None) -> None:
    import pyvista as pv

    part = list(state["parts"].values())[-1]
    grid = cutcells.to_pyvista(part["surface"].visualization_mesh(mode="cut_only", degree=3))
    surface = grid.extract_surface(nonlinear_subdivision=2, algorithm="dataset_surface")
    surface.cell_data.clear()

    plotter = pv.Plotter(off_screen=screenshot is not None, window_size=(1000, 800))
    plotter.add_mesh(cutcells.to_pyvista(state["mesh"]).outline(), color="#334155")
    plotter.add_mesh(surface, color="#f59e42", smooth_shading=True)
    plotter.camera_position = "iso"
    plotter.camera.Elevation(15)
    if screenshot is not None:
        plotter.screenshot(str(screenshot))
        plotter.close()
    else:
        plotter.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--n", type=int, default=17, help="grid points per side")
    parser.add_argument("--analytic", action="store_true", help="also run the analytic level set in Python (minutes)")
    parser.add_argument("--order", type=int, default=4, help="quadrature order")
    parser.add_argument("--output-dir", type=Path, default=Path("demo_torus_output"))
    parser.add_argument("--no-plot", action="store_true", help="skip the PyVista window")
    parser.add_argument("--screenshot", type=Path, help="render the view to a PNG instead")
    args = parser.parse_args()

    state = build(args.n, order=args.order, analytic=args.analytic)
    report(state)
    write(state, args.output_dir)
    if args.screenshot is not None or not args.no_plot:
        plot(state, args.screenshot)


if __name__ == "__main__":
    main()
