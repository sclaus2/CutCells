#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""A napkin ring: a level set made of level sets, with exact sharp edges.

A cylindrical hole of radius a is drilled through a ball of radius R. What is
left is a napkin ring of height h = 2 sqrt(R^2 - a^2), with the famous volume

    |ring| = pi h^3 / 6,

which depends on h alone. Its surface has two faces, the sphere zone of area
2 pi R h and the wall of the hole of area 2 pi a h, which meet in two circles
of radius a: sharp edges.

One level set for the whole shape, such as max(ball, -cylinder), has a kink
along those edges, and its interpolant smears them. Instead the ring is a
level set made of two smooth level sets, the exact ball and a P3 cylinder, with the inside
``(ball < 0) & (cyl > 0)``. CutCells cuts with both, so the edges are exact,
and ``ring == 0`` gives the surface, ``.on("cyl")`` the face on the cylinder.

Run it with ``python demo_napkin_ring.py``; ``--no-plot`` skips the PyVista
window.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import cutcells

BALL_RADIUS = 0.8
HOLE_RADIUS = 0.45
SECTION_Y = 0.02  # the plane that cuts the ring open in the view


def cylinder(X: np.ndarray) -> np.ndarray:
    return np.sqrt(X[0] ** 2 + X[1] ** 2) - HOLE_RADIUS


def build(n: int = 13, degree: int = 3, order: int = 4) -> dict:
    mesh = cutcells.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n, n, n, cell_type="tetrahedron")

    ball = cutcells.analytic_sphere([0.0, 0.0, 0.0], BALL_RADIUS, name="ball")
    cyl = cutcells.create_level_set(mesh, cylinder, degree=degree, name="cyl")
    ring = cutcells.LevelSetFunction(level_sets=[ball, cyl], inside=(ball < 0) & (cyl > 0), name="ring")

    result = cutcells.cut(mesh, ring)
    surface = result[ring == 0]
    return {
        "mesh": mesh,
        "level_set": ring,
        "result": result,
        "ring": result[ring < 0],
        "surface": surface,
        "zone": surface.on("ball"),
        "hole": surface.on("cyl"),
        "rims": result[(ball == 0) & (cyl == 0)],
        "order": order,
    }


def report(state: dict) -> None:
    result, order = state["result"], state["order"]
    R, a = BALL_RADIUS, HOLE_RADIUS
    h = 2.0 * np.sqrt(R**2 - a**2)
    print(f"level sets: {result.level_set_names}, insides: {result.insides}")
    print(f"cut cells: {result.num_cut_cells}, ring height h = {h:.6f}")
    exact = {
        "ring": np.pi * h**3 / 6.0,
        "surface": 2.0 * np.pi * (R + a) * h,
        "zone": 2.0 * np.pi * R * h,
        "hole": 2.0 * np.pi * a * h,
        "rims": 2.0 * 2.0 * np.pi * a,
    }
    print(f"{'part':>8}  {'dim':>3}  {'measure':>14}  {'exact':>14}  {'error':>9}")
    for key, value in exact.items():
        part = state[key]
        measure = part.quadrature(order=order).weights.sum()
        print(f"{key:>8}  {part.dim:>3}  {measure:14.10f}  {value:14.10f}  {abs(measure - value):9.2e}")


def write(state: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    state["ring"].write_vtu(str(output_dir / "ring.vtu"), mode="full", degree=3)
    state["zone"].write_vtu(str(output_dir / "zone.vtu"), mode="cut_only", degree=3)
    state["hole"].write_vtu(str(output_dir / "hole.vtu"), mode="cut_only", degree=3)
    print(f"wrote ring.vtu, zone.vtu and hole.vtu to {output_dir}")


def plot(state: dict, screenshot: Path | None = None) -> None:
    import pyvista as pv

    def surface(part, mode="cut_only", subdivision=3):
        grid = cutcells.to_pyvista(part.visualization_mesh(mode=mode, degree=3))
        return grid.extract_surface(nonlinear_subdivision=subdivision, algorithm="dataset_surface")

    # Left: the two faces and the rims. Right: the solid ring cut open at y = 0.
    plotter = pv.Plotter(shape=(1, 2), off_screen=screenshot is not None, window_size=(1400, 700))
    plotter.subplot(0, 0)
    plotter.add_mesh(surface(state["zone"]), color="#f59e42", smooth_shading=True)
    plotter.add_mesh(surface(state["hole"]), color="#2563eb", smooth_shading=True)
    plotter.add_mesh(surface(state["rims"], subdivision=4).tube(radius=0.01), color="#be185d")
    plotter.add_text("ring == 0: .on(\"ball\") and .on(\"cyl\")", font_size=11, color="#334155")
    plotter.subplot(0, 1)
    # Cut the ring open with the plane y = y0, a second level set: its surface
    # behind the plane and the section in it, coloured by whether the
    # background cell is cut by the ring.
    ring = state["level_set"]
    plane = cutcells.create_level_set(state["mesh"], lambda X: X[1] - SECTION_Y, degree=1, name="plane")
    opened = cutcells.cut(state["mesh"], [ring, plane])
    section = surface(opened[(plane == 0) & (ring < 0)], subdivision=2)
    section.cell_data["cut"] = np.isin(section.cell_data["parent"], state["result"].cut_cells).astype(int)
    plotter.add_mesh(surface(opened[(ring == 0) & (plane > 0)]), color="#f59e42", smooth_shading=True)
    plotter.add_mesh(section, scalars="cut", cmap=["#68b984", "#f59e42"], show_scalar_bar=False, lighting=False)
    plotter.add_text("ring < 0, cut open", font_size=11, color="#334155")
    plotter.link_views()
    plotter.camera_position = [(1.6, -3.4, 2.4), (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]
    if screenshot is not None:
        plotter.screenshot(str(screenshot))
        plotter.close()
    else:
        plotter.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--n", type=int, default=13, help="grid points per side")
    parser.add_argument("--degree", type=int, default=3, help="degree of the cylinder level set")
    parser.add_argument("--order", type=int, default=4, help="quadrature order")
    parser.add_argument("--output-dir", type=Path, default=Path("demo_napkin_ring_output"))
    parser.add_argument("--no-plot", action="store_true", help="skip the PyVista window")
    parser.add_argument("--screenshot", type=Path, help="render the view to a PNG instead")
    args = parser.parse_args()

    state = build(args.n, args.degree, args.order)
    report(state)
    write(state, args.output_dir)
    if args.screenshot is not None or not args.no_plot:
        plot(state, args.screenshot)


if __name__ == "__main__":
    main()
