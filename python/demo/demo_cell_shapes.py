#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""One level set on every cell shape, with both backends.

The same disk (2D) or ball (3D) is cut out of meshes of triangles,
quadrilaterals, tetrahedra, hexahedra, prisms and pyramids that share one grid.
The demo reports what each shape costs and gains, and how a P1 level set
differs: on simplices its zero set is flat and the two backends agree exactly,
on the other shapes it is curved.

Run it with ``python demo_cell_shapes.py``; ``--no-plot`` skips the PyVista window.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

import cutcells

CENTRE = (0.1, -0.05, 0.03)
RADIUS = 0.6
SHAPES = ("triangle", "quadrilateral", "tetrahedron", "hexahedron", "prism", "pyramid")

# The reference cells, whose physical coordinates are their reference coordinates
REFERENCE = {
    "triangle": (5, 2, [[0, 0], [1, 0], [0, 1]]),
    "quadrilateral": (9, 2, [[0, 0], [1, 0], [1, 1], [0, 1]]),
    "tetrahedron": (10, 3, [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    "hexahedron": (12, 3, [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]),
    "prism": (13, 3, [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]]),
    "pyramid": (14, 3, [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1]]),
}


def round_shape(X: np.ndarray) -> np.ndarray:
    return np.sqrt(sum((X[d] - CENTRE[d]) ** 2 for d in range(X.shape[0]))) - RADIUS


def background(shape: str, n: int):
    if shape in ("triangle", "quadrilateral"):
        return cutcells.rectangle_mesh(-1, -1, 1, 1, n, n, cell_type=shape)
    return cutcells.box_mesh(-1, -1, -1, 1, 1, 1, n, n, n, cell_type=shape)


def exact_measures(tdim: int) -> tuple[float, float]:
    if tdim == 2:
        return np.pi * RADIUS**2, 2.0 * np.pi * RADIUS
    return 4.0 / 3.0 * np.pi * RADIUS**3, 4.0 * np.pi * RADIUS**2


def survey(shape: str, n: int, degree: int, backend: str, order: int) -> dict:
    mesh = background(shape, n)
    phi = cutcells.create_level_set(mesh, round_shape, degree=degree, name="phi")
    result = cutcells.cut(mesh, phi, backend=backend)
    start = time.perf_counter()
    volume = result[phi < 0].quadrature(order=order, mode="cut_only")
    surface = result[phi == 0].quadrature(order=order)
    seconds = time.perf_counter() - start
    exact_volume, exact_area = exact_measures(mesh.tdim)
    whole = result[phi < 0].quadrature(order=order).weights.sum()
    holding = np.unique(surface.parent_map[np.diff(surface.offset) > 0]).size
    return {
        "cells": mesh.num_cells(),
        "cut": result.num_cut_cells,
        "holding": holding,
        "volume points": volume.weights.size / holding,
        "surface points": surface.weights.size / holding,
        "total points": volume.weights.size + surface.weights.size,
        "seconds": seconds,
        "volume": whole,
        "area": surface.weights.sum(),
        "volume error": abs(whole - exact_volume),
        "area error": abs(surface.weights.sum() - exact_area),
    }


def build(n: int = 13, order: int = 6) -> dict:
    curved = {}
    for shape in SHAPES:
        degree = 2 if shape == "pyramid" else 3  # pyramids take degree 2 at most
        for backend in ("quadrays", "lut"):
            curved[shape, backend] = (degree, survey(shape, n, degree, backend, order))
    flat = {(shape, backend): survey(shape, n, 1, backend, 2) for shape in SHAPES for backend in ("quadrays", "lut")}
    return {"n": n, "order": order, "curved": curved, "flat": flat}


def report(state: dict) -> None:
    print(f"{state['n']} grid points per side, quadrature order {state['order']}, P3 level set (P2 on pyramids)")
    print(f"{'shape':>14}{'backend':>10}{'cells':>7}{'marked cut':>12}{'on surface':>12}{'all points':>12}"
          f"{'per surface cell':>18}{'seconds':>9}{'volume err':>12}{'area err':>10}")
    for (shape, backend), (_, row) in state["curved"].items():
        points = f"{row['volume points']:.0f} + {row['surface points']:.0f}"
        print(f"{shape:>14}{backend:>10}{row['cells']:7d}{row['cut']:12d}{row['holding']:12d}{row['total points']:12d}"
              f"{points:>18}{row['seconds']:9.3f}{row['volume error']:12.1e}{row['area error']:10.1e}")
    print("\nP1 level set, order 2")
    print(f"{'shape':>14}{'backend':>10}{'per surface cell':>18}{'seconds':>9}{'volume':>18}{'area':>18}")
    for (shape, backend), row in state["flat"].items():
        points = f"{row['volume points']:.0f} + {row['surface points']:.0f}"
        print(f"{shape:>14}{backend:>10}{points:>18}{row['seconds']:9.3f}{row['volume']:18.12f}{row['area']:18.12f}")


def one_cell(shape: str, backend: str, order: int):
    """Points of one reference cell cut by a ball around a point near its origin."""
    vtk, tdim, vertices = REFERENCE[shape]
    vertices = np.asarray(vertices, dtype=np.float64)
    cells = np.arange(len(vertices), dtype=np.int32)
    offsets = np.array([0, len(vertices)], dtype=np.int32)
    mesh = cutcells.MeshView(vertices, cells, offsets, np.array([vtk], dtype=np.int32), tdim=tdim)
    centre = np.array([-0.1, -0.05, -0.05])[:tdim]
    radius = 0.62

    def ball(X):
        return np.sqrt(sum((X[d] - centre[d]) ** 2 for d in range(tdim))) - radius

    phi = cutcells.create_level_set(mesh, ball, degree=2, name="phi")
    result = cutcells.cut(mesh, phi, backend=backend)
    pad = lambda p: np.column_stack([p, np.zeros((p.shape[0], 3 - tdim))])
    volume = result[phi < 0].quadrature(order=order, mode="cut_only")
    surface = result[phi == 0].quadrature(order=order)
    return mesh, pad(volume.points.reshape(-1, tdim)), pad(surface.points.reshape(-1, tdim))


def plot(state: dict, screenshot: Path | None = None) -> None:
    import pyvista as pv

    order = 4
    plotter = pv.Plotter(shape=(2, len(SHAPES)), off_screen=screenshot is not None,
                         window_size=(1800, 640), border=False)
    for row, backend in enumerate(("quadrays", "lut")):
        for column, shape in enumerate(SHAPES):
            plotter.subplot(row, column)
            mesh, volume, surface = one_cell(shape, backend, order)
            grid = cutcells.to_pyvista(mesh)
            plotter.add_mesh(grid, style="wireframe", color="#64748b", line_width=2)
            plotter.add_points(pv.PolyData(volume), color="#64748b", point_size=5, render_points_as_spheres=True)
            plotter.add_points(pv.PolyData(surface), color="#eb6834", point_size=6, render_points_as_spheres=True)
            label = f"{backend}, {shape}\n{volume.shape[0]} + {surface.shape[0]} points"
            plotter.add_text(label, font_size=9, color="#1e293b", position="upper_left")
            if mesh.tdim == 2:
                plotter.view_xy()
            else:
                plotter.view_isometric()
            plotter.camera.zoom(0.9)
    if screenshot is not None:
        plotter.screenshot(str(screenshot))
        plotter.close()
    else:
        plotter.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--n", type=int, default=13, help="grid points per side")
    parser.add_argument("--order", type=int, default=6, help="quadrature order")
    parser.add_argument("--no-plot", action="store_true", help="skip the PyVista window")
    parser.add_argument("--screenshot", type=Path, help="render the view to a PNG instead")
    args = parser.parse_args()
    state = build(args.n, args.order)
    report(state)
    if args.screenshot is not None or not args.no_plot:
        plot(state, args.screenshot)


if __name__ == "__main__":
    main()
