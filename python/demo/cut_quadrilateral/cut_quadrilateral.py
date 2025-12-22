"""Demo: cut a level set through quadrilateral cells using CutCells.

Run (from repo root, fenicsx-dev env):
  conda run -n fenicsx-dev python python/demo/cut_quadrilateral/cut_quadrilateral.py

Optional:
  conda run -n fenicsx-dev python python/demo/cut_quadrilateral/cut_quadrilateral.py --write-vtk out_prefix

This uses the installed CutCells Python bindings, which call the C++ implementation
(and therefore the generated quadrilateral clip tables).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import cutcells


def cut_one_quad(ls_values: np.ndarray, write_vtk_prefix: str | None = None) -> None:
    # Unit square quad: (0,0)->(1,0)->(1,1)->(0,1)
    vertex_coordinates = np.array([0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0], dtype=float)
    gdim = 2

    inside = cutcells.cut(
        cutcells.CellType.quadrilateral,
        vertex_coordinates,
        gdim,
        ls_values,
        "phi<0",
        False,
    )
    outside = cutcells.cut(
        cutcells.CellType.quadrilateral,
        vertex_coordinates,
        gdim,
        ls_values,
        "phi>0",
        False,
    )
    interface = cutcells.cut(
        cutcells.CellType.quadrilateral,
        vertex_coordinates,
        gdim,
        ls_values,
        "phi=0",
        False,
    )

    print(f"ls={ls_values.tolist()}")
    print(f"  inside:  volume={inside.volume():.6f}, vtk_types={inside.types}")
    print(f"  outside: volume={outside.volume():.6f}, vtk_types={outside.types}")
    print(f"  iface:   cells={len(interface.types)}, vtk_types={interface.types}")

    if write_vtk_prefix:
        inside.write_vtk(f"{write_vtk_prefix}_inside.vtu")
        outside.write_vtk(f"{write_vtk_prefix}_outside.vtu")
        interface.write_vtk(f"{write_vtk_prefix}_interface.vtu")
        print(f"  wrote: {write_vtk_prefix}_inside.vtu / _outside.vtu / _interface.vtu")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cut a level set through a quadrilateral"
    )
    parser.add_argument(
        "--write-vtk",
        type=str,
        default=None,
        help="Prefix for VTK output files (writes *_inside.vtu, *_outside.vtu, *_interface.vtu)",
    )
    args = parser.parse_args()

    # Two canonical ambiguity checks + a generic slanted cut
    cases = [
        (np.array([-0.5, 0.5, -0.5, 0.5], dtype=float), "opp_corners"),
        (np.array([0.5, -0.5, -0.5, -0.5], dtype=float), "single_outside_corner"),
        (np.array([-1.0, -0.2, 0.7, 1.3], dtype=float), "slanted"),
    ]

    for ls, name in cases:
        prefix = None
        if args.write_vtk:
            prefix = str(
                Path(args.write_vtk).with_name(f"{Path(args.write_vtk).name}_{name}")
            )
        cut_one_quad(ls, prefix)


if __name__ == "__main__":
    main()
