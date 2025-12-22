"""Normalize VTK clip tables to CutCells token streams."""

from __future__ import annotations

from typing import Dict, List, Sequence

from .topology import num_edges
from .vtk_tables_import import ClipCase, PointRef

# Mapping matches enum order in cpp/src/cell_types.h
CELL_TYPE_TO_ENUM: Dict[str, int] = {
    "point": 0,
    "interval": 1,
    "triangle": 2,
    "tetrahedron": 3,
    "quadrilateral": 4,
    "hexahedron": 5,
    "prism": 6,
    "pyramid": 7,
}


def token_from_ref(ref: PointRef) -> int:
    kind, idx = ref
    if kind == "V":
        return 100 + idx
    if kind == "E":
        return idx
    if kind == "N":
        return 200 + idx
    raise ValueError(f"Unsupported point reference: {ref}")


def build_intersected_edges(case: ClipCase, num_edges_in_cell: int) -> List[int]:
    flags = [0] * num_edges_in_cell
    for e in case.intersected_edges:
        if e < 0 or e >= num_edges_in_cell:
            raise ValueError(f"Edge id {e} outside range 0..{num_edges_in_cell - 1}")
        flags[e] = 1
    return flags


def normalize_cases(cell_type: str, cases: Sequence[ClipCase]) -> Dict[str, List]:
    n_edges = num_edges(cell_type)

    intersected_edges: List[List[int]] = []
    sub_element_offset: List[int] = [0]
    sub_element_cell_types: List[int] = []
    sub_element: List[int] = []

    special_point_count: List[int] = []
    special_point_offset: List[int] = [0]
    special_point_data: List[int] = []

    for case in cases:
        intersected_edges.append(build_intersected_edges(case, n_edges))

        # Special points (e.g. VTK centroid N0)
        # Encoding per case: repeated blocks [nrefs, ref0, ref1, ...]
        sp_list = list(case.special_points)
        sp_list.sort(key=lambda sp: sp.point_id)
        special_point_count.append(len(sp_list))
        for sp in sp_list:
            # Token for the special point itself is implicit: 200 + sp.point_id
            # References must not include other special points.
            for ref in sp.refs:
                if ref[0] == "N":
                    raise ValueError("Special point refs must not include 'N' refs")
            special_point_data.append(len(sp.refs))
            special_point_data.extend(token_from_ref(r) for r in sp.refs)
        special_point_offset.append(len(special_point_data))

        for emitted in case.cells:
            if emitted.cell_type not in CELL_TYPE_TO_ENUM:
                raise ValueError(f"Unknown emitted cell type: {emitted.cell_type}")
            sub_element_cell_types.append(CELL_TYPE_TO_ENUM[emitted.cell_type])
            sub_element.append(len(emitted.vertices))
            sub_element.extend(token_from_ref(ref) for ref in emitted.vertices)

        sub_element_offset.append(len(sub_element_cell_types))

    out = {
        "intersected_edges": intersected_edges,
        "sub_element_offset": sub_element_offset,
        "sub_element_cell_types": sub_element_cell_types,
        "sub_element": sub_element,
        "extra_coords": [],  # placeholder; fill if Steiner points are added later
    }

    if any(c > 0 for c in special_point_count):
        out["special_point_count"] = special_point_count
        out["special_point_offset"] = special_point_offset
        out["special_point_data"] = special_point_data

    return out
