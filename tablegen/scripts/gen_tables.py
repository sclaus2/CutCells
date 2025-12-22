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
from pathlib import Path
import sys

# Allow running this script without installing the tablegen package
_THIS_DIR = Path(__file__).resolve().parent
_TABLEGEN_ROOT = _THIS_DIR.parent
sys.path.insert(0, str(_TABLEGEN_ROOT))

from cutcells_tablegen.emit_c_headers import emit_all_tet_like
from cutcells_tablegen.normalize import normalize_cases
from cutcells_tablegen.vtk_tables_import import get_cases, registered_cell_types


def generate_for_cell(cell_type: str, out_dir: Path) -> None:
    tables = {}
    for kind in ("inside", "outside", "interface"):
        cases = get_cases(cell_type, kind)
        tables[kind] = normalize_cases(cell_type, cases)
    emit_all_tet_like(cell_type, tables, out_dir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate CutCells clip tables (tet-like format)"
    )
    available = registered_cell_types()
    parser.add_argument(
        "--cell-type",
        choices=available,
        nargs="*",
        default=available,
        help="Cell types to generate",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("cpp/src/generated"),
        help="Output directory for generated headers",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for cell_type in args.cell_type:
        generate_for_cell(cell_type, args.out)


if __name__ == "__main__":
    main()
