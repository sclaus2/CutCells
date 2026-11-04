#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""A plane that lies on mesh facets: pieces shared by two cells are integrated once.

Fitted geometry is common: a level set can vanish on whole facets of the mesh.
No cell is cut there, and the surface belongs to both cells that share each
facet. The demo uses two planes in a tetrahedral mesh of [-1, 1]^3, one level
set at a time,

    x = 0.25 on grid faces,   x = y on the diagonal faces of the tetrahedra,

integrates the surface once on the facets (``entity_quadrature``), and folds
those rules into the cell rules with each policy of ``quadrature(entities=...)``:
every policy integrates the same surface, and sends it to different cells. The
exact areas are 4 and 4 sqrt(2).

Run it with ``python demo_mesh_aligned.py``; ``--no-plot`` skips the PyVista window.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import cutcells

POLICIES = ("split", "negative", "positive", "owner", "exclude")


def build(n: int = 9, order: int = 2) -> dict:
    mesh = cutcells.box_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n, n, n, cell_type="tetrahedron")
    planes = {
        "x = 0.25": (lambda X: X[0] - 0.25, 4.0, 4.0 * 1.25),
        "x = y": (lambda X: X[0] - X[1], 4.0 * np.sqrt(2.0), 4.0),
    }
    state = {"mesh": mesh, "order": order, "planes": {}}
    for name, (function, area, volume) in planes.items():
        phi = cutcells.create_level_set(mesh, function, degree=1, name="phi")
        result = cutcells.cut(mesh, phi)
        state["planes"][name] = {"phi": phi, "result": result, "face": result[phi == 0],
                                 "negative": result[phi < 0], "area": area, "volume": volume}
    return state


def report(state: dict) -> None:
    order = state["order"]
    print(f"tetrahedra: {state['mesh'].num_cells()}")
    for name, plane in state["planes"].items():
        face = plane["face"]
        on = face.entity_quadrature(order=order).facets
        print(f"{name}: cut cells {plane['result'].num_cut_cells}, facets on the plane {on.cell_offsets.size - 1}")
        print(f"  on the facets, once each: area {on.weights.sum():.14f}, exact {plane['area']:.14f}")
        print(f"  volume phi < 0: {plane['negative'].quadrature(order=order, mode='full').weights.sum():.14f},"
              f" exact {plane['volume']:.14f}")
        print(f"  {'entities=':>12}  {'area in cell rules':>20}  {'cells with points':>18}")
        for policy in POLICIES:
            rules = face.quadrature(order=order, entities=policy)
            print(f"  {policy:>12}  {rules.weights.sum():20.14f}  {np.unique(rules.parent_map).size:18d}")


def write(state: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for k, plane in enumerate(state["planes"].values()):
        plane["face"].write_vtu(str(output_dir / f"plane_{k}.vtu"), mode="cut_only", degree=1)
    print(f"wrote plane_0.vtu and plane_1.vtu to {output_dir}")


def plot(state: dict, screenshot: Path | None = None) -> None:
    import pyvista as pv

    background = cutcells.to_pyvista(state["mesh"])
    plotter = pv.Plotter(shape=(1, 2), off_screen=screenshot is not None, window_size=(1400, 700))
    for k, (name, plane) in enumerate(state["planes"].items()):
        plotter.subplot(0, k)
        grid = cutcells.to_pyvista(plane["face"].visualization_mesh(mode="cut_only", degree=1))
        plotter.add_mesh(background.outline(), color="#334155")
        plotter.add_mesh(grid, color="#f59e42", show_edges=True, edge_color="#7c2d12", lighting=False)
        plotter.add_text(name, font_size=11, color="#334155")
        plotter.camera_position = "iso"
        plotter.camera.Azimuth(-30)  # the iso view looks along the plane x = y
    if screenshot is not None:
        plotter.screenshot(str(screenshot))
        plotter.close()
    else:
        plotter.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--n", type=int, default=9, help="grid points per side, 8 k + 1 so that x = 0.25 is on the grid")
    parser.add_argument("--order", type=int, default=2, help="quadrature order")
    parser.add_argument("--output-dir", type=Path, default=Path("demo_mesh_aligned_output"))
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
