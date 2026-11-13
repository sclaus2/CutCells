#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""A gyroid in a hexahedral mesh: a triply periodic minimal surface.

The gyroid is approximated by the zero set of

    g(x) = sin(kx) cos(ky) + sin(ky) cos(kz) + sin(kz) cos(kx),  k = 2 pi / L,

on the unit cell [0, L]^3. It splits the cell into two congruent labyrinths,
so each of {g < 0} and {g > 0} fills exactly half of the cell. Its area in
one unit cell is about 3.09 L^2.

The demo cuts a hexahedral mesh by a P3 interpolant of g and

* checks the volume fraction 1/2 of {g < 0};
* integrates the area of {g = 0}, with quadrays and with the straight pieces
  of the lookup tables;
* selects a sheet of thickness about 2 t / |grad g| with two level sets,
  ``(lower > 0) & (upper < 0)`` for lower = g + t and upper = g - t, and
  integrates its volume.

Run it with ``python demo_gyroid.py``; ``--no-plot`` skips the PyVista window.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import cutcells

LENGTH = 1.0
SHEET = 0.4  # t, the half width of the sheet in values of g


def gyroid(X: np.ndarray) -> np.ndarray:
    k = 2.0 * np.pi / LENGTH
    x, y, z = k * X[0], k * X[1], k * X[2]
    return np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)


def build(n: int = 17, degree: int = 3, order: int = 4) -> dict:
    mesh = cutcells.box_mesh(0.0, 0.0, 0.0, LENGTH, LENGTH, LENGTH, n, n, n, cell_type="hexahedron")
    g = cutcells.create_level_set(mesh, gyroid, degree=degree, name="g")
    result = cutcells.cut(mesh, g)

    lower = cutcells.create_level_set(mesh, lambda X: gyroid(X) + SHEET, degree=degree, name="lower")
    upper = cutcells.create_level_set(mesh, lambda X: gyroid(X) - SHEET, degree=degree, name="upper")
    sheet_result = cutcells.cut(mesh, [lower, upper])

    return {
        "mesh": mesh,
        "result": result,
        "solid": result[g < 0],
        "surface": result[g == 0],
        "sheet": sheet_result[(lower > 0) & (upper < 0)],
        "sheet_faces": (sheet_result[lower == 0], sheet_result[upper == 0]),
        "order": order,
    }


def report(state: dict) -> None:
    mesh, result, order = state["mesh"], state["result"], state["order"]
    print(f"hexahedra: {mesh.num_cells()}, cut by g: {result.num_cut_cells}")
    fraction = state["solid"].quadrature(order=order).weights.sum() / LENGTH**3
    area = state["surface"].quadrature(order=order).weights.sum() / LENGTH**2
    straight = state["surface"].quadrature(order=order, backend="lut").weights.sum() / LENGTH**2
    sheet = state["sheet"].quadrature(order=order).weights.sum() / LENGTH**3
    print(f"volume fraction of g < 0:  {fraction:.12f}  (exact 1/2, error {abs(fraction - 0.5):.1e})")
    print(f"area of g = 0, quadrays:   {area:.8f} L^2")
    print(f"area of g = 0, lut:        {straight:.8f} L^2")
    print(f"volume fraction of sheet:  {sheet:.8f}")


def write(state: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    state["surface"].write_vtu(str(output_dir / "gyroid.vtu"), mode="cut_only", degree=3)
    state["sheet"].write_vtu(str(output_dir / "sheet.vtu"), mode="full", degree=2)
    print(f"wrote gyroid.vtu and sheet.vtu to {output_dir}")


def plot(state: dict, screenshot: Path | None = None) -> None:
    import pyvista as pv

    background = cutcells.to_pyvista(state["mesh"])
    surface = cutcells.to_pyvista(state["surface"].visualization_mesh(mode="cut_only", degree=3))
    surface = surface.extract_surface(nonlinear_subdivision=3, algorithm="dataset_surface")
    surface.point_data["z"] = surface.points[:, 2]
    # The sheet is drawn by its two faces, lower == 0 and upper == 0.
    faces = [cutcells.to_pyvista(face.visualization_mesh(mode="cut_only", degree=3))
             .extract_surface(nonlinear_subdivision=3, algorithm="dataset_surface")
             for face in state["sheet_faces"]]

    plotter = pv.Plotter(shape=(1, 2), off_screen=screenshot is not None, window_size=(1400, 700))
    plotter.subplot(0, 0)
    plotter.add_mesh(background.outline(), color="#334155")
    plotter.add_mesh(surface, scalars="z", cmap="viridis", smooth_shading=True, show_scalar_bar=False)
    plotter.add_text("g == 0", font_size=11, color="#334155")
    plotter.subplot(0, 1)
    plotter.add_mesh(background.outline(), color="#334155")
    plotter.add_mesh(faces[0], color="#2563eb", smooth_shading=True)
    plotter.add_mesh(faces[1], color="#f59e42", smooth_shading=True)
    plotter.add_text("(lower > 0) & (upper < 0)", font_size=11, color="#334155")
    plotter.link_views()
    plotter.camera_position = "iso"
    if screenshot is not None:
        plotter.screenshot(str(screenshot))
        plotter.close()
    else:
        plotter.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--n", type=int, default=17, help="grid points per side")
    parser.add_argument("--degree", type=int, default=3, help="degree of the level sets")
    parser.add_argument("--order", type=int, default=4, help="quadrature order")
    parser.add_argument("--output-dir", type=Path, default=Path("demo_gyroid_output"))
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
