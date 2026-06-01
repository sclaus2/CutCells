#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""
Demo: cut a tetrahedral mesh with two level sets.

The example uses a sphere and a plane:
  sphere(x) = |x - c| - r
  plane(x)  = n . x - offset

It writes straight VTU output for several mesh-part selections, including
single-interface parts, a codimension-two intersection curve, and volume
regions selected by both level sets.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import cutcells
from cutcells import box_tetrahedron_mesh, mesh_from_pyvista


DEFAULT_CENTER = np.array([0.1, -0.08, 0.05], dtype=np.float64)
DEFAULT_RADIUS = 0.58
DEFAULT_PLANE_NORMAL = np.array([1.0, -0.35, 0.2], dtype=np.float64)
DEFAULT_PLANE_OFFSET = 0.0
TRIANGULATE = True
OUTPUT_DIR = Path("demo_two_level_sets_output")


def sphere_phi(X: np.ndarray, center: np.ndarray, radius: float) -> np.ndarray:
    return (
        np.sqrt(
            (X[0] - center[0]) ** 2 + (X[1] - center[1]) ** 2 + (X[2] - center[2]) ** 2
        )
        - radius
    )


def plane_phi(X: np.ndarray, normal: np.ndarray, offset: float) -> np.ndarray:
    return normal[0] * X[0] + normal[1] * X[1] + normal[2] * X[2] - offset


def write_part(part, path: Path, mode: str = "cut_only") -> None:
    part.write_vtu(str(path), mode=mode)
    print(
        f"  {path.name:<38} dim={part.dim}, "
        f"cut={part.num_cut_cells}, uncut={part.num_uncut_cells}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cut a tetrahedral mesh with a sphere and a plane level set."
    )
    parser.add_argument(
        "--n",
        type=int,
        default=7,
        help="Number of grid points per coordinate direction.",
    )
    args = parser.parse_args()

    if args.n < 2:
        raise SystemExit("--n must be at least 2")

    center = DEFAULT_CENTER
    normal = DEFAULT_PLANE_NORMAL.copy()
    norm = float(np.linalg.norm(normal))
    if norm == 0.0:
        raise RuntimeError("DEFAULT_PLANE_NORMAL must be nonzero")
    normal /= norm

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Building tetrahedral box mesh with n={args.n} ...")
    grid = box_tetrahedron_mesh(
        -1.0,
        -1.0,
        -1.0,
        1.0,
        1.0,
        1.0,
        args.n,
        args.n,
        args.n,
    )
    mesh = mesh_from_pyvista(grid, tdim=3)
    print(f"  cells={mesh.num_cells()}, points={mesh.num_nodes()}")

    background_path = OUTPUT_DIR / "background_mesh.vtu"
    grid.save(background_path)
    print(f"  wrote {background_path}")

    print("Interpolating level sets ...")
    sphere = cutcells.create_level_set(
        mesh,
        lambda X: sphere_phi(X, center, DEFAULT_RADIUS),
        degree=1,
        name="sphere",
    )
    plane = cutcells.create_level_set(
        mesh,
        lambda X: plane_phi(X, normal, DEFAULT_PLANE_OFFSET),
        degree=1,
        name="plane",
    )
    print(f"  sphere dofs={sphere.mesh_data.num_dofs()}")
    print(f"  plane dofs={plane.mesh_data.num_dofs()}")

    print(f"Running cut(mesh, [sphere, plane], triangulate={TRIANGULATE}) ...")
    result = cutcells.cut(mesh, [sphere, plane], triangulate=TRIANGULATE)
    print(f"  level sets={list(result.level_set_names)}")
    print(f"  num_cut_cells={result.num_cut_cells}")

    outputs = [
        ("sphere = 0 and plane < 0", "sphere_interface_plane_negative.vtu", "cut_only"),
        ("sphere = 0 and plane > 0", "sphere_interface_plane_positive.vtu", "cut_only"),
        ("plane = 0 and sphere < 0", "plane_disk_inside_sphere.vtu", "cut_only"),
        ("sphere = 0 and plane = 0", "sphere_plane_intersection_curve.vtu", "cut_only"),
        ("sphere < 0 and plane < 0", "inside_sphere_plane_negative.vtu", "cut_only"),
        ("sphere < 0 and plane > 0", "inside_sphere_plane_positive.vtu", "cut_only"),
        ("sphere < 0 or plane < 0", "inside_sphere_or_plane_negative.vtu", "cut_only"),
    ]

    print("Writing selected mesh parts ...")
    for expr, filename, mode in outputs:
        part = result[expr]
        write_part(part, OUTPUT_DIR / filename, mode=mode)


if __name__ == "__main__":
    main()
