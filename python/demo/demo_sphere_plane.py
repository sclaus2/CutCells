#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""Two level sets: a ball and a plane, and the parts they make together.

A ball {|x| < R} and the half space {z < z0} are cut at once. Conditions on
both select parts of three dimensions:

* ``(ball < 0) & (plane < 0)``: the ball below the plane, a volume;
* ``(plane == 0) & (ball < 0)``: the disk where the plane crosses the ball;
* ``(ball == 0) & (plane == 0)``: the circle where the sphere meets the plane.

The exact measures are, with the disk radius a = sqrt(R^2 - z0^2),

    |ball below| = 4/3 pi R^3 - pi (R - z0)^2 (2 R + z0) / 3,
    |disk| = pi a^2,   |circle| = 2 pi a,

and the sphere below the plane, ``(ball == 0) & (plane < 0)``, has the area
2 pi R (R + z0).

Run it with ``python demo_sphere_plane.py``; ``--no-plot`` skips the PyVista
window.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import cutcells

RADIUS = 0.8
HEIGHT = 0.3


def build(n: int = 13, order: int = 4) -> dict:
    mesh = cutcells.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n, n, n, cell_type="tetrahedron")
    ball = cutcells.analytic_sphere([0.0, 0.0, 0.0], RADIUS, name="ball")
    plane = cutcells.create_level_set(mesh, lambda X: X[2] - HEIGHT, degree=1, name="plane")
    result = cutcells.cut(mesh, [ball, plane])

    a = np.sqrt(RADIUS**2 - HEIGHT**2)
    selections = {
        "ball below the plane": ((ball < 0) & (plane < 0),
                                 4.0 / 3.0 * np.pi * RADIUS**3
                                 - np.pi * (RADIUS - HEIGHT) ** 2 * (2.0 * RADIUS + HEIGHT) / 3.0),
        "sphere below the plane": ((ball == 0) & (plane < 0), 2.0 * np.pi * RADIUS * (RADIUS + HEIGHT)),
        "disk": ((plane == 0) & (ball < 0), np.pi * a**2),
        "circle": ((ball == 0) & (plane == 0), 2.0 * np.pi * a),
    }
    parts = {name: (result[condition], condition, exact)
             for name, (condition, exact) in selections.items()}
    return {"mesh": mesh, "result": result, "parts": parts, "order": order}


def report(state: dict) -> None:
    result = state["result"]
    print(f"level sets: {result.level_set_names}, cut cells: {result.num_cut_cells}")
    print(f"{'part':>24}  {'selection':>28}  {'dim':>3}  {'measure':>16}  {'error':>9}")
    for name, (part, condition, exact) in state["parts"].items():
        measure = part.quadrature(order=state["order"]).weights.sum()
        print(f"{name:>24}  {str(condition):>28}  {part.dim:>3}  {measure:16.12f}  {abs(measure - exact):9.2e}")


def write(state: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, (part, _, _) in state["parts"].items():
        filename = name.replace(" ", "_") + ".vtu"
        part.write_vtu(str(output_dir / filename), mode="full" if part.dim == 3 else "cut_only", degree=3)
    print(f"wrote {len(state['parts'])} VTU files to {output_dir}")


def plot(state: dict, screenshot: Path | None = None) -> None:
    import pyvista as pv

    def surface(name: str, subdivision: int = 3):
        part = state["parts"][name][0]
        mode = "full" if part.dim == 3 else "cut_only"
        grid = cutcells.to_pyvista(part.visualization_mesh(mode=mode, degree=3))
        return grid.extract_surface(nonlinear_subdivision=subdivision, algorithm="dataset_surface")

    plotter = pv.Plotter(off_screen=screenshot is not None, window_size=(1000, 800))
    plotter.add_mesh(cutcells.to_pyvista(state["mesh"]).outline(), color="#334155")
    plotter.add_mesh(surface("sphere below the plane"), color="#f59e42", smooth_shading=True)
    plotter.add_mesh(surface("disk"), color="#2563eb", opacity=0.85)
    plotter.add_mesh(surface("circle", 4).tube(radius=0.012), color="#be185d", smooth_shading=True)
    plotter.camera_position = "iso"
    plotter.camera.Elevation(10)
    plotter.camera.Zoom(1.05)
    if screenshot is not None:
        plotter.screenshot(str(screenshot))
        plotter.close()
    else:
        plotter.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--n", type=int, default=13, help="grid points per side")
    parser.add_argument("--order", type=int, default=4, help="quadrature order")
    parser.add_argument("--output-dir", type=Path, default=Path("demo_sphere_plane_output"))
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
