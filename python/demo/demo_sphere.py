#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""A ball in a tetrahedral mesh: exact and interpolated geometry.

The background mesh covers the cube [-1, 1]^3. The ball {|x - c| < R} is given
twice: as the analytic sphere ``cutcells.analytic_sphere``, which CutCells
evaluates exactly, and as a P2 interpolant of the same function on the mesh.
For both the demo integrates the volume and the surface area,

    |B| = 4/3 pi R^3,   |S| = 4 pi R^2,

and shows the curved pieces of the cut cells with PyVista: the sphere, and
the inside pieces of the cut cells behind a clipping plane.

Run it with ``python demo_sphere.py``; ``--no-plot`` skips the PyVista window.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import cutcells

CENTER = (0.05, -0.02, 0.03)
RADIUS = 0.7


def sphere(X: np.ndarray) -> np.ndarray:
    return np.sqrt(sum((X[i] - CENTER[i]) ** 2 for i in range(3))) - RADIUS


def build(n: int = 13, order: int = 4) -> dict:
    mesh = cutcells.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n, n, n, cell_type="tetrahedron")

    exact = cutcells.analytic_sphere(CENTER, RADIUS, name="ball")
    interpolated = cutcells.create_level_set(mesh, sphere, degree=2, name="ball")

    state = {"mesh": mesh, "order": order}
    for key, ball in (("exact", exact), ("P2", interpolated)):
        result = cutcells.cut(mesh, ball)
        state[key] = {"ball": ball, "result": result, "inside": result[ball < 0],
                      "surface": result[ball == 0]}
    return state


def report(state: dict) -> None:
    mesh, order = state["mesh"], state["order"]
    volume = 4.0 / 3.0 * np.pi * RADIUS**3
    area = 4.0 * np.pi * RADIUS**2
    print(f"background cells: {mesh.num_cells()}, cut cells: {state['exact']['result'].num_cut_cells}")
    print(f"{'level set':>10}  {'volume error':>13}  {'area error':>11}")
    for key in ("exact", "P2"):
        v = state[key]["inside"].quadrature(order=order).weights.sum()
        a = state[key]["surface"].quadrature(order=order).weights.sum()
        print(f"{key:>10}  {abs(v - volume):13.2e}  {abs(a - area):11.2e}")


def write(state: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    state["exact"]["inside"].write_vtu(str(output_dir / "ball.vtu"), mode="full", degree=3)
    state["exact"]["surface"].write_vtu(str(output_dir / "sphere.vtu"), mode="cut_only", degree=3)
    print(f"wrote ball.vtu and sphere.vtu (Lagrange cells of degree 3) to {output_dir}")


def plot(state: dict, screenshot: Path | None = None) -> None:
    import pyvista as pv

    parts = state["exact"]
    background = cutcells.to_pyvista(state["mesh"])
    sphere_surface = cutcells.to_pyvista(parts["surface"].visualization_mesh(mode="cut_only", degree=3))
    ball = cutcells.to_pyvista(parts["inside"].visualization_mesh(mode="full", degree=3))

    # Left: the curved sphere. Right: the ball behind the plane y = c_y, whole
    # cells in green and the curved pieces of the cut cells in orange.
    plotter = pv.Plotter(shape=(1, 2), off_screen=screenshot is not None, window_size=(1400, 700))
    plotter.subplot(0, 0)
    plotter.add_mesh(background.outline(), color="#334155")
    plotter.add_mesh(sphere_surface.extract_surface(nonlinear_subdivision=3, algorithm="dataset_surface"),
                     color="#f59e42", smooth_shading=True)
    plotter.add_text("ball == 0", font_size=11, color="#334155")
    plotter.subplot(0, 1)
    plotter.add_mesh(background.outline(), color="#334155")
    ball.cell_data["cut"] = np.isin(ball.cell_data["parent"], parts["result"].cut_cells).astype(int)
    half = ball.extract_cells(ball.cell_centers().points[:, 1] < CENTER[1])
    plotter.add_mesh(half.extract_surface(nonlinear_subdivision=2, algorithm="dataset_surface"),
                     scalars="cut", cmap=["#68b984", "#f59e42"], show_scalar_bar=False)
    plotter.add_mesh(half.extract_surface(algorithm="dataset_surface").extract_all_edges(),
                     color="#334155", line_width=0.5)
    plotter.add_text("ball < 0, mode=\"full\"", font_size=11, color="#334155")
    plotter.link_views()
    plotter.camera_position = "iso"
    plotter.camera.Azimuth(-20)
    if screenshot is not None:
        plotter.screenshot(str(screenshot))
        plotter.close()
    else:
        plotter.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--n", type=int, default=13, help="grid points per side")
    parser.add_argument("--order", type=int, default=4, help="quadrature order")
    parser.add_argument("--output-dir", type=Path, default=Path("demo_sphere_output"))
    parser.add_argument("--no-plot", action="store_true", help="skip the PyVista window")
    parser.add_argument("--screenshot", type=Path, help="render the view to a PNG instead")
    args = parser.parse_args()

    state = build(args.n, args.order)
    report(state)
    write(state, args.output_dir)
    if args.screenshot is not None or not args.no_plot:
        plot(state, args.screenshot)


if __name__ == "__main__":
    main()
