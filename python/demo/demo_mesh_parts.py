#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""
Demo: MeshView + create_level_set + cut(mesh, ls) straight output bridge.

This is the first user-facing step for the HO pipeline:
1. build a pyvista tetrahedral mesh,
2. convert it to a MeshView,
3. build one polynomial level set with create_level_set(...),
4. call cut(mesh, ls),
5. select HOMeshPart objects,
6. obtain straight hybrid visualization meshes and straight quadrature from HOMeshPart,
7. write VTU files.

Curving is intentionally out of scope here. This demo only exercises the
straight decomposition / straight output bridge.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

try:
    import pyvista as pv
except ImportError as exc:
    raise SystemExit(
        "pyvista is required for this demo. "
        "Install with: python -m pip install pyvista\n"
        f"Import error: {exc}"
    )

import cutcells
from cutcells import box_tetrahedron_mesh, mesh_from_pyvista


CENTER = np.array([0.1, -0.05, 0.0], dtype=np.float64)
RADIUS = 0.55


def phi_batch(X: np.ndarray) -> np.ndarray:
    return np.sqrt(
        (X[0] - CENTER[0]) ** 2
        + (X[1] - CENTER[1]) ** 2
        + (X[2] - CENTER[2]) ** 2
    ) - RADIUS


def cutmesh_to_pyvista(cut_mesh) -> pv.UnstructuredGrid:
    return pv.UnstructuredGrid(
        np.asarray(cut_mesh.cells, dtype=np.int64),
        np.asarray(cut_mesh.vtk_types, dtype=np.uint8),
        np.asarray(cut_mesh.vertex_coords, dtype=np.float64),
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Demo: use cut(mesh, ls) on a MeshView from a pyvista tetra mesh."
    )
    parser.add_argument(
        "--n",
        type=int,
        default=8,
        help="Structured tetrahedral mesh resolution parameter.",
    )
    parser.add_argument(
        "--degree",
        type=int,
        default=2,
        help="Polynomial degree for create_level_set(...).",
    )
    parser.add_argument(
        "--quadrature-order",
        type=int,
        default=4,
        help="Quadrature order for HOMeshPart.quadrature(...).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("demo_meshview_ho_cut_output"),
        help="Directory for VTU files.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip interactive pyvista plotting.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Creating tetrahedral pyvista mesh with n={args.n} ...")
    grid = box_tetrahedron_mesh(
        -1.0, -1.0, -1.0,
         1.0,  1.0,  1.0,
        args.n, args.n, args.n,
    )
    mesh = mesh_from_pyvista(grid, tdim=3)
    print(f"  pyvista cells={grid.n_cells}, points={grid.n_points}")
    print(f"  MeshView gdim={mesh.gdim}, tdim={mesh.tdim}")

    print(f"Building polynomial level set with degree={args.degree} ...")
    ls = cutcells.create_level_set(mesh, phi_batch, degree=args.degree, name="phi")
    print(f"  level-set dofs={ls.mesh_data.num_dofs()}")

    print("Running cut(mesh, ls) ...")
    result = cutcells.cut(mesh, ls)
    negative = result["phi < 0"]
    interface = result["phi = 0"]
    positive = result["phi > 0"]

    print(f"  num_cut_cells={result.num_cut_cells}")
    print(f"  phi < 0  : cut={negative.num_cut_cells}, uncut={negative.num_uncut_cells}")
    print(f"  phi = 0  : cut={interface.num_cut_cells}, uncut={interface.num_uncut_cells}")
    print(f"  phi > 0  : cut={positive.num_cut_cells}, uncut={positive.num_uncut_cells}")

    print("Building straight hybrid visualization meshes from HOMeshPart ...")
    inside_full = negative.visualization_mesh(mode="full")
    inside_cut = negative.visualization_mesh(mode="cut_only")
    interface_mesh = interface.visualization_mesh(mode="cut_only")

    print(f"  inside full cells={len(np.asarray(inside_full.types))}")
    print(f"  inside cut-only cells={len(np.asarray(inside_cut.types))}")
    print(f"  interface cells={len(np.asarray(interface_mesh.types))}")

    print("Computing straight quadrature from HOMeshPart ...")
    q_inside_full = negative.quadrature(
        order=args.quadrature_order,
        mode="full",
    )
    q_inside_cut = negative.quadrature(
        order=args.quadrature_order,
        mode="cut_only",
    )
    q_interface = interface.quadrature(
        order=args.quadrature_order,
        mode="cut_only",
    )

    print(f"  inside full quadrature sum    = {q_inside_full.weights.sum():.12g}")
    print(f"  inside cut-only quadrature sum= {q_inside_cut.weights.sum():.12g}")
    print(f"  interface quadrature sum      = {q_interface.weights.sum():.12g}")

    inside_full_path = args.output_dir / "phi_negative_full_straight.vtu"
    inside_cut_path = args.output_dir / "phi_negative_cut_only_straight.vtu"
    interface_path = args.output_dir / "phi_interface_straight.vtu"

    negative.write_vtu(str(inside_full_path), mode="full")
    negative.write_vtu(str(inside_cut_path), mode="cut_only")
    interface.write_vtu(str(interface_path), mode="cut_only")

    print(f"  wrote {inside_full_path}")
    print(f"  wrote {inside_cut_path}")
    print(f"  wrote {interface_path}")

    if args.no_plot:
        return

    pv_inside_full = cutmesh_to_pyvista(inside_full)
    pv_interface = cutmesh_to_pyvista(interface_mesh)

    plotter = pv.Plotter(shape=(1, 2), title="CutCells HO straight-output demo")
    plotter.subplot(0, 0)
    plotter.add_title("phi < 0, mode='full'")
    plotter.add_mesh(grid, style="wireframe", color="lightgrey", opacity=0.25)
    plotter.add_mesh(pv_inside_full, color="steelblue", opacity=0.75, show_edges=True)

    plotter.subplot(0, 1)
    plotter.add_title("phi = 0")
    plotter.add_mesh(grid, style="wireframe", color="lightgrey", opacity=0.15)
    plotter.add_mesh(pv_interface, color="crimson", show_edges=True, line_width=2)

    plotter.link_views()
    plotter.show()


if __name__ == "__main__":
    main()
