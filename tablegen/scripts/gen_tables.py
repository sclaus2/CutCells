"""Generate CutCells clip tables (tetrahedron-like format).

Usage:
    conda activate fenicsx-dev
    python tablegen/scripts/gen_tables.py --cell-type quadrilateral

Generates tables in tet-like format (fixed-width 2D arrays with -1 padding)
instead of VTK's variable-length format.
Default output: cpp/src/generated/
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

# Allow running this script without installing the tablegen package
_THIS_DIR = Path(__file__).resolve().parent
_TABLEGEN_ROOT = _THIS_DIR.parent
sys.path.insert(0, str(_TABLEGEN_ROOT))


_VTK_DERIVED_CELL_TYPES = {"hexahedron", "prism", "pyramid"}


def _import_tablegen():
    # Import after we potentially set CUT_CELLS_VTK_REF/CUT_CELLS_VTK_HEADER
    # so vtk_tables_import loads tables from the intended VTK source.
    from cutcells_tablegen.emit_c_headers import emit_all_tet_like
    from cutcells_tablegen.normalize import normalize_cases
    from cutcells_tablegen.vtk_tables_import import get_cases, registered_cell_types

    return emit_all_tet_like, normalize_cases, get_cases, registered_cell_types


def generate_for_cell(
    cell_type: str,
    out_dir: Path,
    vtk_ref: str | None,
    vtk_header_path: str | None,
) -> None:
    emit_all_tet_like, normalize_cases, get_cases, _ = _import_tablegen()

    tables = {}
    for kind in ("inside", "outside", "interface"):
        cases = get_cases(cell_type, kind)
        tables[kind] = normalize_cases(cell_type, cases)
    provenance = "vtk" if cell_type in _VTK_DERIVED_CELL_TYPES else None
    emit_all_tet_like(
        cell_type,
        tables,
        out_dir,
        provenance=provenance,
        vtk_ref=vtk_ref if provenance == "vtk" else None,
        vtk_header_path=vtk_header_path if provenance == "vtk" else None,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate CutCells clip tables (tet-like format)"
    )
    parser.add_argument(
        "--cell-type",
        nargs="*",
        default=None,
        help="Cell types to generate (default: all registered types)",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("cpp/src/generated"),
        help="Output directory for generated headers",
    )
    parser.add_argument(
        "--vtk-ref",
        default=None,
        help="VTK git ref/tag/commit to download vtkTableBasedClipCases.h from (default: env CUT_CELLS_VTK_REF or 'master')",
    )
    parser.add_argument(
        "--vtk-header-path",
        default=None,
        help="Local path to vtkTableBasedClipCases.h (default: env CUT_CELLS_VTK_HEADER)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.vtk_ref is not None:
        os.environ["CUT_CELLS_VTK_REF"] = str(args.vtk_ref)
    if args.vtk_header_path is not None:
        os.environ["CUT_CELLS_VTK_HEADER"] = str(args.vtk_header_path)

    emit_all_tet_like, _, _, registered_cell_types = _import_tablegen()
    vtk_ref = os.environ.get("CUT_CELLS_VTK_REF", "master")
    vtk_header_path = os.environ.get("CUT_CELLS_VTK_HEADER")

    cell_types = args.cell_type
    if not cell_types:
        cell_types = registered_cell_types()

    for cell_type in cell_types:
        generate_for_cell(cell_type, args.out, vtk_ref=vtk_ref, vtk_header_path=vtk_header_path)


if __name__ == "__main__":
    main()
