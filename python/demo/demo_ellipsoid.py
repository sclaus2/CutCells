#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""Integrating functions on a rotated ellipsoid: moments, normals and the divergence theorem.

A quadrature rule is points and weights, so any function can be integrated
with it. The ellipsoid has semi-axes (a, b, c) and centre x0, in axes turned
by a rotation Q, with y = Q (x - x0):

    phi(x) = sqrt((y_1/a)^2 + (y_2/b)^2 + (y_3/c)^2) - 1.

The demo maps the points of the rules to physical space and checks

    |Omega| = 4/3 pi a b c,
    int_Omega x dx = |Omega| x0,
    int_Omega y_i^2 dx = |Omega| a_i^2 / 5,

and, with the unit normal grad phi / |grad phi| at the surface points, the
divergence theorem for F(x) = x - x0, whose divergence is 3:

    int_Gamma F . n ds = 3 |Omega|.

Run it with ``python demo_ellipsoid.py``; ``--no-plot`` skips the PyVista window.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import cutcells

CENTER = np.array([0.04, -0.03, 0.02])
AXES = np.array([0.8, 0.55, 0.35])


def rotation() -> np.ndarray:
    """A rotation that keeps the ellipsoid's axes off the mesh's."""
    cx, sx = np.cos(0.5), np.sin(0.5)
    cz, sz = np.cos(0.7), np.sin(0.7)
    return np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]]) @ np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])


Q = rotation()


def ellipsoid(X: np.ndarray) -> np.ndarray:
    y = Q @ (X.T - CENTER).T
    return np.sqrt(np.sum((y.T / AXES) ** 2, axis=1)) - 1.0


def normals(x: np.ndarray) -> np.ndarray:
    """Unit normals of phi at the points ``x``, shape (num_points, 3)."""
    y = (x - CENTER) @ Q.T
    grad = (y / AXES**2) @ Q
    return grad / np.linalg.norm(grad, axis=1, keepdims=True)


def physical_points(mesh, rules) -> np.ndarray:
    """Map the reference points of rules on a tetrahedral mesh to physical space."""
    reference = rules.points.reshape(-1, 3)
    cells = np.repeat(rules.parent_map, np.diff(rules.offset))
    x = mesh.coordinates[mesh.connectivity.reshape(-1, 4)[cells]]
    return x[:, 0] + np.einsum("pi,pig->pg", reference, x[:, 1:] - x[:, :1])


def build(n: int = 17, degree: int = 3, order: int = 6) -> dict:
    mesh = cutcells.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n, n, n, cell_type="tetrahedron")
    phi = cutcells.create_level_set(mesh, ellipsoid, degree=degree, name="phi")
    result = cutcells.cut(mesh, phi)
    inside, surface = result[phi < 0], result[phi == 0]
    return {
        "mesh": mesh,
        "phi": phi,
        "result": result,
        "inside": inside,
        "surface": surface,
        "inside_rules": inside.quadrature(order=order, mode="full"),
        "surface_rules": surface.quadrature(order=order),
    }


def report(state: dict) -> None:
    mesh = state["mesh"]
    print(f"tetrahedra: {mesh.num_cells()}, cut cells: {state['result'].num_cut_cells}")

    rules = state["inside_rules"]
    x = physical_points(mesh, rules)
    w = rules.weights
    volume = w.sum()
    exact_volume = 4.0 / 3.0 * np.pi * np.prod(AXES)
    y = (x - CENTER) @ Q.T

    surface = state["surface_rules"]
    xs = physical_points(mesh, surface)
    flux = np.dot(surface.weights, np.sum((xs - CENTER) * normals(xs), axis=1))

    rows = [("volume", volume, exact_volume), ("flux of x - x0", flux, 3.0 * exact_volume)]
    rows += [(f"centroid {'xyz'[i]}", np.dot(w, x[:, i]) / volume, CENTER[i]) for i in range(3)]
    rows += [(f"int y_{i + 1}^2", np.dot(w, y[:, i] ** 2), exact_volume * AXES[i] ** 2 / 5.0) for i in range(3)]
    print(f"{'quantity':>16}  {'computed':>16}  {'exact':>16}  {'error':>9}")
    for name, value, exact in rows:
        print(f"{name:>16}  {value:16.12f}  {exact:16.12f}  {abs(value - exact):9.2e}")
    print(f"surface area: {surface.weights.sum():.10f} (no closed form)")


def write(state: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    state["inside"].write_vtu(str(output_dir / "ellipsoid.vtu"), mode="full", degree=3)
    state["surface"].write_vtu(str(output_dir / "surface.vtu"), mode="cut_only", degree=3)
    print(f"wrote ellipsoid.vtu and surface.vtu (Lagrange cells of degree 3) to {output_dir}")


def plot(state: dict, screenshot: Path | None = None) -> None:
    import pyvista as pv

    grid = cutcells.to_pyvista(state["surface"].visualization_mesh(mode="cut_only", degree=3))
    surface = grid.extract_surface(nonlinear_subdivision=2, algorithm="dataset_surface")

    rules = state["surface_rules"]
    x = physical_points(state["mesh"], rules)
    points = pv.PolyData(x)
    points["F . n"] = np.sum((x - CENTER) * normals(x), axis=1)

    plotter = pv.Plotter(shape=(1, 2), off_screen=screenshot is not None, window_size=(1400, 700))
    plotter.subplot(0, 0)
    plotter.add_mesh(surface, color="#f59e42", smooth_shading=True)
    plotter.add_text("phi == 0", font_size=11, color="#334155")
    plotter.subplot(0, 1)
    plotter.add_mesh(points, scalars="F . n", cmap="viridis", point_size=3, render_points_as_spheres=True,
                     scalar_bar_args={"title": "F . n", "color": "#334155"})
    plotter.add_text("the surface's quadrature points", font_size=11, color="#334155")
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
    parser.add_argument("--degree", type=int, default=3, help="degree of the level set")
    parser.add_argument("--order", type=int, default=6, help="quadrature order")
    parser.add_argument("--output-dir", type=Path, default=Path("demo_ellipsoid_output"))
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
