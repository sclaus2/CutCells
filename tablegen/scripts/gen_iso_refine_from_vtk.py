#!/usr/bin/env python3
"""Generate iso-refine table blocks from VTK tessellation, remapped to Basix ordering.

Pipeline:
1) Call C++ helper (tablegen/tools/vtk_iso_refine_dump) to tessellate a VTK higher-order cell.
2) Query Basix equispaced P Lagrange interpolation points for the same cell/order.
3) Map VTK point ids to Basix point ids by coordinate matching in reference space.
4) Emit C++ arrays + RefinementTemplate block.

This script does not edit source files automatically; it writes generated code to a C++
header (default: cpp/src/generated/iso_refine_vtk_blocks.h, configurable via --out)
for manual integration.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Iterable

CELL_TO_CUTCELLS = {
    "interval": "cell::type::interval",
    "triangle": "cell::type::triangle",
    "quadrilateral": "cell::type::quadrilateral",
    "tetrahedron": "cell::type::tetrahedron",
    "hexahedron": "cell::type::hexahedron",
    "prism": "cell::type::prism",
    "pyramid": "cell::type::pyramid",
}


def child_simplex_type(cell: str) -> str:
    if cell in ("interval",):
        return "interval"
    if cell in ("triangle", "quadrilateral"):
        return "triangle"
    return "tetrahedron"


def flatten_rows(rows: Iterable[Iterable[int]]) -> list[int]:
    out: list[int] = []
    for r in rows:
        out.extend(r)
    return out


def fmt_int_array(name: str, values: list[int]) -> str:
    vals = ", ".join(str(v) for v in values)
    return f"static const int {name}[] = {{ {vals} }};"


def fmt_double_array(name: str, values: list[float]) -> str:
    vals = ", ".join(f"{v:.17g}" for v in values)
    return f"static const double {name}[] = {{ {vals} }};"


def map_points(vtk_points: list[list[float]], basix_points, tdim: int, tol: float) -> list[int]:
    used = set()
    mapping = [-1] * len(vtk_points)
    for i, p in enumerate(vtk_points):
        best = -1
        best_d = 1.0e100
        for j, q in enumerate(basix_points):
            if j in used:
                continue
            d2 = 0.0
            for d in range(tdim):
                dd = float(p[d]) - float(q[d])
                d2 += dd * dd
            if d2 < best_d:
                best_d = d2
                best = j
        if best < 0 or best_d > tol * tol:
            raise RuntimeError(
                f"Could not map VTK point {i}={p[:tdim]} to Basix points (best d2={best_d})"
            )
        mapping[i] = best
        used.add(best)
    return mapping


def parent_masks_from_basix(element) -> tuple[list[int], list[int]]:
    n = element.dim
    pd = [0] * n
    pid = [0] * n
    for dim, entities in enumerate(element.entity_dofs):
        for eid, dofs in enumerate(entities):
            for dof in dofs:
                pd[dof] = dim
                pid[dof] = eid
    return pd, pid


def run_helper(helper: Path, cell: str, order: int) -> dict:
    cmd = [str(helper), "--cell", cell, "--order", str(order)]
    out = subprocess.check_output(cmd, text=True)
    return json.loads(out)


def emit_block(cell: str, order: int, dump: dict, tol: float) -> str:
    import basix

    if not any(hasattr(basix, name) for name in ("CellType", "_basixcpp", "create_element", "finite_element")):
        raise RuntimeError(
            "Imported 'basix' does not look like FEniCS Basix. "
            f"Found module at {getattr(basix, '__file__', '<unknown>')} "
            f"with version {getattr(basix, '__version__', '<unknown>')}. "
            "Install FEniCS Basix (python package 'fenics-basix') and remove the unrelated 'basix' package."
        )

    def _resolve(name: str):
        if hasattr(basix, name):
            return getattr(basix, name)
        if hasattr(basix, "_basixcpp") and hasattr(basix._basixcpp, name):
            return getattr(basix._basixcpp, name)
        raise AttributeError(f"Basix symbol '{name}' not found")

    CellType = _resolve("CellType")
    ElementFamily = _resolve("ElementFamily")
    LagrangeVariant = _resolve("LagrangeVariant")
    create_element = getattr(basix, "create_element", None)
    if create_element is None and hasattr(basix, "finite_element"):
        create_element = getattr(basix.finite_element, "create_element", None)
    if create_element is None:
        raise AttributeError("Basix 'create_element' API not found")

    cell_to_basix = {
        "interval": CellType.interval,
        "triangle": CellType.triangle,
        "quadrilateral": CellType.quadrilateral,
        "tetrahedron": CellType.tetrahedron,
        "hexahedron": CellType.hexahedron,
        "prism": CellType.prism,
        "pyramid": CellType.pyramid,
    }

    tdim = int(dump["tdim"])
    vtk_points = dump["points"]
    simplices = dump["simplices_vtk"]
    bs_cell = cell_to_basix[cell]
    elem = create_element(
        ElementFamily.P,
        bs_cell,
        order,
        LagrangeVariant.equispaced,
    )

    basix_pts = elem.points
    vtk_to_basix = map_points(vtk_points, basix_pts, tdim, tol)
    remapped = [[vtk_to_basix[i] for i in s] for s in simplices]

    parent_dim, parent_id = parent_masks_from_basix(elem)
    basix_flat = [float(x) for row in basix_pts for x in row[:tdim]]
    conn_flat = flatten_rows(remapped)
    child = child_simplex_type(cell)
    child_vpc = len(remapped[0]) if remapped else (tdim + 1)

    tag = f"{cell}_p{order}_vtk_tess_basix"
    lines = []
    lines.append(f"// ---- {cell} P{order} from VTK tessellation (Basix remapped) ----")
    lines.append(fmt_double_array(f"{tag}_ref_coords", basix_flat))
    lines.append(fmt_int_array(f"{tag}_parent_dim", parent_dim))
    lines.append(fmt_int_array(f"{tag}_parent_id", parent_id))
    lines.append(fmt_int_array(f"{tag}_cells", conn_flat))
    lines.append(
        "static const RefinementTemplate "
        f"{tag}_tpl = {{"
        f"{len(parent_dim)}, {len(remapped)}, {tdim}, {child_vpc}, "
        f"{CELL_TO_CUTCELLS[cell]}, {CELL_TO_CUTCELLS[child]}, "
        f"std::vector<double>({tag}_ref_coords, {tag}_ref_coords + {len(basix_flat)}), "
        f"std::vector<int>({tag}_parent_dim, {tag}_parent_dim + {len(parent_dim)}), "
        f"std::vector<int>({tag}_parent_id, {tag}_parent_id + {len(parent_id)}), "
        f"std::vector<int>({tag}_cells, {tag}_cells + {len(conn_flat)})}};"
    )
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--helper", type=Path, required=True, help="Path to vtk_iso_refine_dump executable")
    p.add_argument(
        "--cell",
        nargs="+",
        default=["triangle", "quadrilateral", "tetrahedron", "hexahedron", "prism"],
        choices=list(CELL_TO_CUTCELLS.keys()),
    )
    p.add_argument("--orders", nargs="+", type=int, default=[2, 3, 4])
    p.add_argument("--tol", type=float, default=1.0e-12, help="Coordinate matching tolerance")
    p.add_argument(
        "--out",
        type=Path,
        default=Path("cpp/src/generated/iso_refine_vtk_blocks.h"),
        help="Output file (default: cpp/src/generated/iso_refine_vtk_blocks.h)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    chunks: list[str] = []
    for cell in args.cell:
        for order in args.orders:
            dump = run_helper(args.helper, cell, order)
            chunks.append(emit_block(cell, order, dump, args.tol))
    body = "\n".join(chunks)
    text = (
        "// AUTO-GENERATED by gen_iso_refine_from_vtk.py — do not edit.\n"
        "// Regenerate with tablegen/scripts/gen_iso_refine_from_vtk.py.\n"
        "// SPDX-License-Identifier: MIT\n"
        "#pragma once\n\n"
        "#include \"iso_refine.h\"\n"
        "#include <vector>\n\n"
        "namespace cutcells::iso_refine_generated\n"
        "{\n"
        f"{body}\n"
        "} // namespace cutcells::iso_refine_generated\n"
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
