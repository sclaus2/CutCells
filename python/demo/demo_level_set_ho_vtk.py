# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier:    MIT

import argparse
from pathlib import Path

import numpy as np

import cutcells
from cutcells import box_tetrahedron_mesh, mesh_from_pyvista, rectangle_triangle_mesh


def _parse_degrees(arg: str) -> list[int]:
    if arg == "all":
        return [2, 3, 4]
    return [int(arg)]


def _phi(X: np.ndarray) -> np.ndarray:
    x = X[:, 0]
    y = X[:, 1]
    z = X[:, 2] if X.shape[1] > 2 else 0.0
    return np.sqrt((x - 0.15) ** 2 + (y + 0.1) ** 2 + (z - 0.05) ** 2) - 0.55


def _run_case(mesh, family: str, degree: int, output_dir: Path) -> None:
    callback_info: dict[str, tuple[int, int]] = {}

    def phi(X: np.ndarray) -> np.ndarray:
        callback_info["shape"] = tuple(X.shape)
        return _phi(X)

    ls = cutcells.interpolate_level_set(mesh, phi, degree=degree)
    output_path = output_dir / f"{family}_p{degree}.vtu"
    cutcells.write_level_set_vtu(str(output_path), ls, field_name="phi")

    print(f"[{family} p{degree}] cells={mesh.num_cells()}")
    print(f"[{family} p{degree}] dofs={ls.mesh_data.num_dofs()}")
    print(f"[{family} p{degree}] callback batch shape={callback_info.get('shape')}")
    print(f"[{family} p{degree}] wrote {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Interpolate a higher-order level set and write Lagrange VTU files."
    )
    parser.add_argument(
        "--cell-family",
        choices=["triangle", "tetra", "both"],
        default="both",
        help="Cell family to run.",
    )
    parser.add_argument(
        "--degree",
        choices=["2", "3", "4", "all"],
        default="all",
        help="Polynomial degree.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("demo_level_set_ho_vtu"),
        help="Output directory for VTU files.",
    )
    parser.add_argument(
        "--n",
        type=int,
        default=8,
        help="Structured mesh resolution parameter.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    degrees = _parse_degrees(args.degree)
    families = (
        ["triangle", "tetra"] if args.cell_family == "both" else [args.cell_family]
    )

    for family in families:
        if family == "triangle":
            grid = rectangle_triangle_mesh(-1.0, -1.0, 1.0, 1.0, args.n, args.n)
            mesh = mesh_from_pyvista(grid, tdim=2)
        else:
            grid = box_tetrahedron_mesh(
                -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, args.n, args.n, args.n
            )
            mesh = mesh_from_pyvista(grid, tdim=3)

        for degree in degrees:
            _run_case(mesh, family, degree, args.output_dir)


if __name__ == "__main__":
    main()
