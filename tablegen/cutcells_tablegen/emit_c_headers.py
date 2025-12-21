"""Emit CutCells-style C++ headers from normalized tables."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

from .topology import num_edges

CELL_ENUM_LITERALS = [
    "cell::type::point",
    "cell::type::interval",
    "cell::type::triangle",
    "cell::type::tetrahedron",
    "cell::type::quadrilateral",
    "cell::type::hexahedron",
    "cell::type::prism",
    "cell::type::pyramid",
]


def _brace_join(rows: List[str]) -> str:
    return "{ " + ", ".join(rows) + " }"


def _format_int_list(values: List[int]) -> str:
    return _brace_join([str(v) for v in values])


def _format_cell_type_list(values: List[int]) -> str:
    mapped = []
    for v in values:
        if v < 0 or v >= len(CELL_ENUM_LITERALS):
            raise ValueError(f"Cell type enum value {v} out of range")
        mapped.append(CELL_ENUM_LITERALS[v])
    return _brace_join(mapped)


def _format_matrix(matrix: List[List[int]]) -> str:
    inner = ["{ " + ", ".join(str(v) for v in row) + " }" for row in matrix]
    return "{\n" + ",\n".join(f"    {row}" for row in inner) + "\n}"


def _header_guard(cell_type: str, table_kind: str) -> str:
    return f"CUT_CELLS_{cell_type.upper()}_{table_kind.upper()}_TABLES_H"


def emit_header(
    cell_type: str, table_kind: str, arrays: Dict[str, List], output_dir: Path
) -> Path:
    """Write a single header for one cell type + table kind (inside/outside/interface)."""
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"cut_{cell_type}_{table_kind}_tables.h"

    n_cases = len(arrays["intersected_edges"])
    n_edges = num_edges(cell_type)
    guard = _header_guard(cell_type, table_kind)

    lines: List[str] = []
    lines.append(f"#ifndef {guard}")
    lines.append(f"#define {guard}")
    lines.append("")
    lines.append('#include "cell_types.h"')
    lines.append("")
    lines.append("namespace cutcells::cell::generated {")
    lines.append("")
    lines.append(
        f"constexpr int cut_{cell_type}_{table_kind}_intersected_edges[{n_cases}][{n_edges}] = "
        + _format_matrix(arrays["intersected_edges"])
        + ";"
    )
    lines.append("")
    lines.append(
        f"constexpr int cut_{cell_type}_{table_kind}_sub_element_offset[{n_cases + 1}] = "
        + _format_int_list(arrays["sub_element_offset"])
        + ";"
    )
    lines.append("")
    lines.append(
        f"constexpr cell::type cut_{cell_type}_{table_kind}_sub_element_cell_types[{len(arrays['sub_element_cell_types'])}] = "
        + _format_cell_type_list(arrays["sub_element_cell_types"])
        + ";"
    )
    lines.append("")
    lines.append(
        f"constexpr int cut_{cell_type}_{table_kind}_sub_element[{len(arrays['sub_element'])}] = "
        + _format_int_list(arrays["sub_element"])
        + ";"
    )

    extra_coords = arrays.get("extra_coords", [])
    if extra_coords:
        lines.append("")
        lines.append(
            f"constexpr double cut_{cell_type}_{table_kind}_extra_coords[{len(extra_coords)}] = "
            + _format_int_list(extra_coords)
            + ";"
        )

    lines.append("")
    lines.append("} // namespace cutcells::cell::generated")
    lines.append("")
    lines.append(f"#endif // {guard}")

    path.write_text("\n".join(lines))
    return path


def emit_all(
    cell_type: str, tables_by_kind: Dict[str, Dict[str, List]], output_dir: Path
) -> List[Path]:
    written: List[Path] = []
    for kind, arrays in tables_by_kind.items():
        written.append(emit_header(cell_type, kind, arrays, output_dir))
    return written


def _emit_tet_like_header(
    cell_type: str, table_kind: str, arrays: Dict[str, List], output_dir: Path
) -> Path:
    """Write header in tetrahedron-like format with fixed-width 2D arrays."""
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"cut_{cell_type}_{table_kind}_tables.h"

    n_cases = len(arrays["intersected_edges"])
    n_edges = num_edges(cell_type)
    guard = _header_guard(cell_type, table_kind)

    # Convert VTK flat format to tet-like format
    # Extract subcells from flat stream
    subcells_data = _extract_subcells_from_flat(arrays)

    lines: List[str] = []
    lines.append(f"#ifndef {guard}")
    lines.append(f"#define {guard}")
    lines.append("")
    lines.append('#include "cell_types.h"')
    lines.append("")
    lines.append(f"namespace cutcells::cell::{cell_type} {{")
    lines.append("")

    # Topology and intersected_edges only in the "inside" header to avoid redefinition
    if table_kind == "inside":
        # Topology (if quad/hex)
        if cell_type == "quadrilateral":
            lines.append("// Topology")
            lines.append("constexpr int edges[4][2] = {{0,1}, {1,2}, {2,3}, {3,0}};")
            lines.append("")
        elif cell_type == "hexahedron":
            lines.append("// Topology")
            lines.append("constexpr int edges[12][2] = {")
            lines.append("    {0,1}, {1,2}, {2,3}, {3,0},  // bottom face")
            lines.append("    {4,5}, {5,6}, {6,7}, {7,4},  // top face")
            lines.append("    {0,4}, {1,5}, {2,6}, {3,7}   // vertical edges")
            lines.append("};")
            lines.append("")
            lines.append("constexpr int faces[6][4] = {")
            lines.append("    {0,3,2,1}, {4,5,6,7},  // bottom, top")
            lines.append("    {0,1,5,4}, {1,2,6,5},  // front, right")
            lines.append("    {2,3,7,6}, {3,0,4,7}   // back, left")
            lines.append("};")
            lines.append("")

        # Intersected edges
        lines.append(
            f"// Intersected edges per case (1 = intersected, 0 = not intersected)"
        )
        lines.append(
            f"constexpr int intersected_edges[{n_cases}][{n_edges}] = "
            + _format_matrix(arrays["intersected_edges"])
            + ";"
        )
        lines.append("")

    # Number of subcells per case
    num_subcells_per_case = []
    for i in range(n_cases):
        start = arrays["sub_element_offset"][i]
        end = arrays["sub_element_offset"][i + 1]
        num_subcells_per_case.append(end - start)

    lines.append(f"// Number of subcells produced for each case ({table_kind} volume)")
    lines.append(
        f"constexpr int num_subcells_{table_kind}[{n_cases}] = "
        + _format_int_list(num_subcells_per_case)
        + ";"
    )
    lines.append("")

    # Offset into subcell array
    lines.append(f"// Offset into subcell array for each case")
    lines.append(
        f"constexpr int case_subcell_offset_{table_kind}[{n_cases + 1}] = "
        + _format_int_list(arrays["sub_element_offset"])
        + ";"
    )
    lines.append("")

    # Cell types
    total_subcells = len(arrays["sub_element_cell_types"])
    lines.append(f"// Cell types for {table_kind} subcells")
    lines.append(
        f"constexpr type subcell_type_{table_kind}[{total_subcells}] = "
        + _format_cell_type_list(arrays["sub_element_cell_types"])
        + ";"
    )
    lines.append("")

    # Subcell vertices (fixed-width 2D array with -1 padding)
    max_verts = subcells_data["max_verts_per_subcell"]
    padded_verts = subcells_data["padded_vertices"]

    lines.append(
        f"// Subcell vertices (max {max_verts} vertices per subcell, -1 padding)"
    )
    lines.append(f"// Tokens: <100 = edge id, >=100 = 100+vertex_id")
    lines.append(
        f"constexpr int subcell_verts_{table_kind}[{total_subcells}][{max_verts}] = "
        + _format_matrix(padded_verts)
        + ";"
    )

    lines.append("")
    lines.append(f"}} // namespace cutcells::cell::{cell_type}")
    lines.append("")
    lines.append(f"#endif // {guard}")

    path.write_text("\n".join(lines))
    return path


def _extract_subcells_from_flat(arrays: Dict[str, List]) -> Dict:
    """Extract subcell data from VTK flat format."""
    sub_element = arrays["sub_element"]
    num_subcells = len(arrays["sub_element_cell_types"])

    subcells = []
    idx = 0
    for _ in range(num_subcells):
        nverts = sub_element[idx]
        verts = sub_element[idx + 1 : idx + 1 + nverts]
        subcells.append(verts)
        idx += 1 + nverts

    # Find max vertices
    max_verts = max(len(v) for v in subcells) if subcells else 0

    # Pad with -1
    padded = []
    for verts in subcells:
        row = list(verts) + [-1] * (max_verts - len(verts))
        padded.append(row)

    return {
        "max_verts_per_subcell": max_verts,
        "padded_vertices": padded,
    }


def emit_all_tet_like(
    cell_type: str, tables_by_kind: Dict[str, Dict[str, List]], output_dir: Path
) -> List[Path]:
    """Emit tables in tetrahedron-like format (fixed-width 2D arrays)."""
    written: List[Path] = []
    for kind, arrays in tables_by_kind.items():
        written.append(_emit_tet_like_header(cell_type, kind, arrays, output_dir))
    return written
