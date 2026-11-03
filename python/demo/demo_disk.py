#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""A disk cut out of a triangular mesh: area, perimeter and an integral.

The background mesh covers the square [-1, 1]^2. A circular level set

    phi(x) = |x - c| - R

splits it into the disk Omega = {phi < 0}, its boundary Gamma = {phi = 0} and
the rest {phi > 0}. The demo selects the three parts, integrates 1 over the
disk and its boundary, integrates f(x) = |x - c|^2 over the disk, and compares
with the exact values

    |Omega| = pi R^2,   |Gamma| = 2 pi R,   int_Omega f dx = pi R^4 / 2.

Run it with ``python demo_disk.py``; ``--no-plot`` skips the PyVista window.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import cutcells

CENTER = (0.1, -0.05)
RADIUS = 0.6


def circle(X: np.ndarray) -> np.ndarray:
    """Signed distance to the circle; ``X`` has shape (gdim, num_points)."""
    return np.sqrt((X[0] - CENTER[0]) ** 2 + (X[1] - CENTER[1]) ** 2) - RADIUS


def physical_points(mesh, rules) -> np.ndarray:
    """Map the reference points of rules on a simplex mesh to physical space."""
    tdim = rules.tdim
    reference = rules.points.reshape(-1, tdim)
    cells = np.repeat(rules.parent_map, np.diff(rules.offset))
    nodes = mesh.connectivity.reshape(-1, tdim + 1)[cells]
    x = mesh.coordinates[nodes]  # (num_points, tdim + 1, gdim)
    return x[:, 0] + np.einsum("pi,pig->pg", reference, x[:, 1:] - x[:, :1])


def build(n: int = 25, degree: int = 2, order: int = 4) -> dict:
    """Cut the mesh and build the quadrature rules of the demo."""
    mesh = cutcells.rectangle_mesh(-1.0, -1.0, 1.0, 1.0, n, n, cell_type="triangle")
    phi = cutcells.create_level_set(mesh, circle, degree=degree, name="phi")
    result = cutcells.cut(mesh, phi)

    inside = result[phi < 0]
    interface = result[phi == 0]
    outside = result[phi > 0]

    return {
        "mesh": mesh,
        "phi": phi,
        "result": result,
        "inside": inside,
        "interface": interface,
        "outside": outside,
        "inside_rules": inside.quadrature(order=order, mode="full"),
        "interface_rules": interface.quadrature(order=order),
    }


def report(state: dict) -> None:
    mesh, result = state["mesh"], state["result"]
    inside, interface, outside = state["inside"], state["interface"], state["outside"]
    print(f"background cells: {mesh.num_cells()}, cut cells: {result.num_cut_cells}")
    print(f"  phi < 0: {inside.num_uncut_cells} whole cells + {inside.num_cut_cells} cut cells")
    print(f"  phi > 0: {outside.num_uncut_cells} whole cells + {outside.num_cut_cells} cut cells")
    print(f"  phi = 0: {interface.num_cut_cells} cut cells, dimension {interface.dim}")

    area_rules = state["inside_rules"]
    perimeter_rules = state["interface_rules"]
    area = area_rules.weights.sum()
    perimeter = perimeter_rules.weights.sum()

    # Integrate f(x) = |x - c|^2 with the points mapped to physical space.
    x = physical_points(mesh, area_rules)
    f = (x[:, 0] - CENTER[0]) ** 2 + (x[:, 1] - CENTER[1]) ** 2
    moment = np.dot(area_rules.weights, f)

    exact = {
        "area": np.pi * RADIUS**2,
        "perimeter": 2.0 * np.pi * RADIUS,
        "int |x-c|^2": 0.5 * np.pi * RADIUS**4,
    }
    computed = {"area": area, "perimeter": perimeter, "int |x-c|^2": moment}
    print(f"{'quantity':>12}  {'computed':>18}  {'exact':>18}  {'error':>9}")
    for key, value in computed.items():
        print(f"{key:>12}  {value:18.14f}  {exact[key]:18.14f}  {abs(value - exact[key]):9.2e}")


def write(state: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    state["inside"].write_vtu(str(output_dir / "disk.vtu"), mode="full", degree=2)
    state["interface"].write_vtu(str(output_dir / "circle.vtu"), mode="cut_only", degree=2)
    state["outside"].write_vtu(str(output_dir / "outside.vtu"), mode="full", degree=2)
    print(f"wrote disk.vtu, circle.vtu and outside.vtu to {output_dir}")


def plot(state: dict, screenshot: Path | None = None) -> None:
    import pyvista as pv

    # Background cells by class: 0 inside, 1 cut, 2 outside.
    background = cutcells.to_pyvista(state["mesh"])
    background.cell_data["class"] = state["result"].domains[0]
    circle_curve = cutcells.to_pyvista(state["interface"].visualization_mesh(mode="cut_only", degree=3))
    circle_curve = circle_curve.extract_surface(nonlinear_subdivision=4, algorithm="dataset_surface")
    circle_curve.points[:, 2] = 0.04  # above the points

    # Quadrature points of the disk in the cut cells.
    rules = state["inside_rules"]
    cells = np.repeat(rules.parent_map, np.diff(rules.offset))
    in_cut = np.isin(cells, state["result"].cut_cells)
    points = physical_points(state["mesh"], rules)[in_cut]
    points = np.column_stack([points, np.full(len(points), 0.03)])

    plotter = pv.Plotter(off_screen=screenshot is not None, window_size=(900, 900))
    plotter.add_mesh(background, scalars="class", show_edges=True, edge_color="#64748b",
                     cmap=["#68b984", "#f59e42", "#e2e8f0"], clim=[-0.5, 2.5], show_scalar_bar=False,
                     lighting=False)
    plotter.add_mesh(circle_curve, color="#0f766e", line_width=5)
    plotter.add_points(points, color="#1d4ed8", point_size=3, render_points_as_spheres=True)
    plotter.view_xy()
    plotter.enable_parallel_projection()  # no perspective shift of the lifted points
    plotter.camera.zoom(1.25)
    if screenshot is not None:
        plotter.screenshot(str(screenshot))
        plotter.close()
    else:
        plotter.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--n", type=int, default=25, help="grid points per side")
    parser.add_argument("--degree", type=int, default=2, help="degree of the level set")
    parser.add_argument("--order", type=int, default=4, help="quadrature order")
    parser.add_argument("--output-dir", type=Path, default=Path("demo_disk_output"))
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
