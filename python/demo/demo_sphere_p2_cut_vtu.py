#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""
Demo: cut one sphere level set and write a P2 Lagrange VTU.

The demo builds a structured tetrahedral MeshView, interpolates one sphere
level set with P=2 nodes, cuts the mesh with CutCells, and writes the curved
P2 zero interface as a higher-order VTU file. The cut path uses
triangulation=true by default.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import cutcells


P = 2
VTK_TETRA = 10
DEFAULT_CENTER = np.array([0.15, -0.1, 0.05], dtype=np.float64)
DEFAULT_RADIUS = 0.55


def sphere_phi(X: np.ndarray, center: np.ndarray, radius: float) -> np.ndarray:
    return np.sqrt(
        (X[0] - center[0]) ** 2
        + (X[1] - center[1]) ** 2
        + (X[2] - center[2]) ** 2
    ) - radius


def structured_tetra_mesh(
    x0: float,
    y0: float,
    z0: float,
    x1: float,
    y1: float,
    z1: float,
    nx: int,
    ny: int,
    nz: int,
) -> cutcells.MeshView:
    xs = np.linspace(x0, x1, num=nx)
    ys = np.linspace(y0, y1, num=ny)
    zs = np.linspace(z0, z1, num=nz)
    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing="ij")
    points = np.c_[xx.ravel(), yy.ravel(), zz.ravel()].astype(np.float64, copy=False)

    def vertex_id(i: int, j: int, k: int) -> int:
        return i + nx * (j + ny * k)

    tets: list[list[int]] = []
    for k in range(nz - 1):
        for j in range(ny - 1):
            for i in range(nx - 1):
                v000 = vertex_id(i, j, k)
                v100 = vertex_id(i + 1, j, k)
                v010 = vertex_id(i, j + 1, k)
                v110 = vertex_id(i + 1, j + 1, k)
                v001 = vertex_id(i, j, k + 1)
                v101 = vertex_id(i + 1, j, k + 1)
                v011 = vertex_id(i, j + 1, k + 1)
                v111 = vertex_id(i + 1, j + 1, k + 1)

                tets.extend(
                    (
                        [v000, v100, v110, v111],
                        [v000, v110, v010, v111],
                        [v000, v010, v011, v111],
                        [v000, v011, v001, v111],
                        [v000, v001, v101, v111],
                        [v000, v101, v100, v111],
                    )
                )

    connectivity = np.asarray(tets, dtype=np.int32).reshape(-1)
    offsets = np.arange(0, connectivity.size + 1, 4, dtype=np.int32)
    cell_types = np.full(offsets.size - 1, VTK_TETRA, dtype=np.int32)

    for cell_id in range(cell_types.size):
        start = int(offsets[cell_id])
        tet = connectivity[start : start + 4].copy()
        vertices = points[tet]
        jacobian = np.column_stack(
            (
                vertices[1] - vertices[0],
                vertices[2] - vertices[0],
                vertices[3] - vertices[0],
            )
        )
        if np.linalg.det(jacobian) < 0.0:
            connectivity[start + 1], connectivity[start + 2] = (
                connectivity[start + 2],
                connectivity[start + 1],
            )

    return cutcells.MeshView(points, connectivity, offsets, cell_types, tdim=3)


def write_background_mesh(path: Path, mesh: cutcells.MeshView) -> None:
    coordinates = np.ascontiguousarray(np.asarray(mesh.coordinates).reshape(-1))
    connectivity = np.ascontiguousarray(np.asarray(mesh.connectivity), dtype=np.int32)
    offsets = np.ascontiguousarray(np.asarray(mesh.offsets), dtype=np.int32)
    cell_types = np.ascontiguousarray(np.asarray(mesh.cell_types), dtype=np.int32)

    cutcells.write_vtk(
        str(path),
        coordinates,
        connectivity,
        offsets,
        cell_types,
        gdim=mesh.gdim,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cut one sphere level set and write a P2 Lagrange VTU."
    )
    parser.add_argument(
        "--n",
        type=int,
        default=7,
        help="Number of grid points per coordinate direction.",
    )
    parser.add_argument(
        "--center",
        type=float,
        nargs=3,
        default=DEFAULT_CENTER.tolist(),
        metavar=("X", "Y", "Z"),
        help="Sphere center.",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=DEFAULT_RADIUS,
        help="Sphere radius.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("demo_sphere_p2_cut_vtu_output"),
        help="Directory for generated VTU files.",
    )
    parser.add_argument(
        "--node-family",
        choices=["lagrange", "gll"],
        default="gll",
        help="Construction/output nodes for the P2 geometry.",
    )
    parser.add_argument(
        "--write-domains",
        action="store_true",
        help="Also write P2 VTUs for phi < 0 and phi > 0 on cut cells.",
    )
    parser.add_argument(
        "--triangulate",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Triangulate cut parts before extracting VTU output.",
    )
    args = parser.parse_args()

    if args.n < 2:
        raise SystemExit("--n must be at least 2")
    if args.radius <= 0.0:
        raise SystemExit("--radius must be positive")

    center = np.asarray(args.center, dtype=np.float64)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Building structured tetrahedral mesh with n={args.n} ...")
    mesh = structured_tetra_mesh(
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
    print(f"  cells={mesh.num_cells()}, points={mesh.num_nodes()}")

    background_path = args.output_dir / "sphere_background_mesh.vtu"
    write_background_mesh(background_path, mesh)
    print(f"  wrote {background_path}")

    print("Interpolating one P2 sphere level set ...")
    ls = cutcells.create_level_set(
        mesh,
        lambda X: sphere_phi(X, center, args.radius),
        degree=P,
        name="phi",
    )
    print(f"  level-set dofs={ls.mesh_data.num_dofs()}")

    print(f"Running CutCells cut(mesh, ls, triangulate={args.triangulate}) ...")
    result = cutcells.cut(mesh, ls, triangulate=args.triangulate)
    negative = result["phi < 0"]
    interface = result["phi = 0"]
    positive = result["phi > 0"]

    print(f"  num_cut_cells={result.num_cut_cells}")
    print(f"  phi < 0 : cut={negative.num_cut_cells}, uncut={negative.num_uncut_cells}")
    print(f"  phi = 0 : cut={interface.num_cut_cells}, uncut={interface.num_uncut_cells}")
    print(f"  phi > 0 : cut={positive.num_cut_cells}, uncut={positive.num_uncut_cells}")

    if interface.num_cut_cells == 0:
        raise RuntimeError(
            "The sphere does not cut the mesh. Adjust --center, --radius, or --n."
        )

    suffix = "_triangulated" if args.triangulate else ""
    interface_path = args.output_dir / f"sphere_interface{suffix}_p2.vtu"
    interface.write_vtu(
        str(interface_path),
        mode="cut_only",
        geometry_order=P,
        node_family=args.node_family,
    )
    print(f"  wrote {interface_path}")

    if args.write_domains:
        inside_path = args.output_dir / f"sphere_inside_cut_only{suffix}_p2.vtu"
        outside_path = args.output_dir / f"sphere_outside_cut_only{suffix}_p2.vtu"
        negative.write_vtu(
            str(inside_path),
            mode="cut_only",
            geometry_order=P,
            node_family=args.node_family,
        )
        positive.write_vtu(
            str(outside_path),
            mode="cut_only",
            geometry_order=P,
            node_family=args.node_family,
        )
        print(f"  wrote {inside_path}")
        print(f"  wrote {outside_path}")


if __name__ == "__main__":
    main()
