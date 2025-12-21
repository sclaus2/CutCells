"""Adapters for VTK TableBasedClip tables.

Populate the stubs with canonical VTK connectivity to regenerate CutCells tables.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Literal, Tuple

from .topology import edge_endpoints

PointRef = Tuple[Literal["V", "E"], int]


@dataclass(frozen=True)
class EmittedCell:
    cell_type: str
    vertices: Tuple[PointRef, ...]


@dataclass(frozen=True)
class ClipCase:
    intersected_edges: Tuple[int, ...]
    cells: Tuple[EmittedCell, ...]


TableKind = Literal["inside", "outside", "interface"]


@dataclass(frozen=True)
class ClipTables:
    inside: List[ClipCase]
    outside: List[ClipCase]
    interface: List[ClipCase]


# TODO: replace stubs with actual VTK TableBasedClip data for each cell type.
VTK_TABLES: dict[str, ClipTables] = {}


def register_tables(cell_type: str, tables: ClipTables) -> None:
    VTK_TABLES[cell_type] = tables


def get_tables(cell_type: str) -> ClipTables:
    tables = VTK_TABLES.get(cell_type)
    if tables is None:
        raise NotImplementedError(
            f"VTK TableBasedClip data for '{cell_type}' not registered yet"
        )
    return tables


def get_cases(cell_type: str, table_kind: TableKind) -> List[ClipCase]:
    tables = get_tables(cell_type)
    if table_kind == "inside":
        return tables.inside
    if table_kind == "outside":
        return tables.outside
    if table_kind == "interface":
        return tables.interface
    raise ValueError(f"Unknown table kind: {table_kind}")


def registered_cell_types() -> List[str]:
    return list(VTK_TABLES.keys())


# --- Quadrilateral TableBasedClip (marching squares style) -----------------


def _mask_bits(mask: int) -> List[bool]:
    return [(mask >> i) & 1 == 1 for i in range(4)]


def _intersected_edges(mask: int) -> Tuple[int, ...]:
    inside = _mask_bits(mask)
    edges = []
    for e, (v0, v1) in enumerate(edge_endpoints("quadrilateral")):
        if inside[v0] != inside[v1]:
            edges.append(e)
    return tuple(edges)


def _polygon_for_mask(mask: int) -> List[PointRef]:
    inside = _mask_bits(mask)
    poly: List[PointRef] = []
    for i in range(4):
        if inside[i]:
            poly.append(("V", i))
        j = (i + 1) % 4
        if inside[i] != inside[j]:
            poly.append(("E", i))
    return poly


def _triangulate(poly: List[PointRef]) -> List[EmittedCell]:
    if len(poly) < 3:
        return []
    if len(poly) == 3:
        return [EmittedCell(cell_type="triangle", vertices=tuple(poly))]
    if len(poly) == 4:
        return [EmittedCell(cell_type="quadrilateral", vertices=tuple(poly))]
    # fan triangulation for n>4
    tris: List[EmittedCell] = []
    anchor = poly[0]
    for i in range(1, len(poly) - 1):
        tris.append(
            EmittedCell(cell_type="triangle", vertices=(anchor, poly[i], poly[i + 1]))
        )
    return tris


def _build_inside_case(mask: int) -> ClipCase:
    inside = _mask_bits(mask)
    inside_count = sum(inside)
    edges = _intersected_edges(mask)

    cells: List[EmittedCell] = []

    if inside_count == 0:
        return ClipCase(intersected_edges=edges, cells=tuple(cells))

    if inside_count == 4:
        cells.append(
            EmittedCell(
                cell_type="quadrilateral",
                vertices=(("V", 0), ("V", 1), ("V", 2), ("V", 3)),
            )
        )
        return ClipCase(intersected_edges=edges, cells=tuple(cells))

    if inside_count == 3:
        # Split pentagon into triangle + quad using diagonal from outside corner to opposite
        outside_vertex = next(i for i, flag in enumerate(inside) if not flag)
        prev_vertex = (outside_vertex + 3) % 4
        next_vertex = (outside_vertex + 1) % 4
        opp_vertex = (outside_vertex + 2) % 4
        edge_prev = prev_vertex
        edge_next = outside_vertex

        cells.append(
            EmittedCell(
                cell_type="triangle",
                vertices=(("E", edge_next), ("V", next_vertex), ("V", opp_vertex)),
            )
        )
        cells.append(
            EmittedCell(
                cell_type="quadrilateral",
                vertices=(
                    ("E", edge_next),
                    ("V", opp_vertex),
                    ("V", prev_vertex),
                    ("E", edge_prev),
                ),
            )
        )
        return ClipCase(intersected_edges=edges, cells=tuple(cells))

    # Opposite corners: emit two disjoint triangles to avoid ambiguity
    if mask == 0b0101:
        cells.append(
            EmittedCell(cell_type="triangle", vertices=(("V", 0), ("E", 0), ("E", 3)))
        )
        cells.append(
            EmittedCell(cell_type="triangle", vertices=(("V", 2), ("E", 1), ("E", 2)))
        )
        return ClipCase(intersected_edges=edges, cells=tuple(cells))

    if mask == 0b1010:
        cells.append(
            EmittedCell(cell_type="triangle", vertices=(("V", 1), ("E", 1), ("E", 0)))
        )
        cells.append(
            EmittedCell(cell_type="triangle", vertices=(("V", 3), ("E", 3), ("E", 2)))
        )
        return ClipCase(intersected_edges=edges, cells=tuple(cells))

    poly = _polygon_for_mask(mask)
    cells.extend(_triangulate(poly))
    return ClipCase(intersected_edges=edges, cells=tuple(cells))


def _build_outside_case(mask: int) -> ClipCase:
    # Ambiguous opposite-corner cases: emit outside as two quads (hex split)
    if mask == 0b0101:
        # inside corners: v0, v2; outside polygon (hex) order: V1, E1, E2, V3, E3, E0
        cells = (
            EmittedCell(
                cell_type="quadrilateral",
                vertices=(("V", 1), ("E", 1), ("E", 2), ("V", 3)),
            ),
            EmittedCell(
                cell_type="quadrilateral",
                vertices=(("V", 3), ("E", 3), ("E", 0), ("V", 1)),
            ),
        )
        return ClipCase(intersected_edges=_intersected_edges(mask), cells=cells)

    if mask == 0b1010:
        # inside corners: v1, v3; outside hex order: V0, E0, E1, V2, E2, E3
        cells = (
            EmittedCell(
                cell_type="quadrilateral",
                vertices=(("V", 0), ("E", 0), ("E", 1), ("V", 2)),
            ),
            EmittedCell(
                cell_type="quadrilateral",
                vertices=(("V", 2), ("E", 2), ("E", 3), ("V", 0)),
            ),
        )
        return ClipCase(intersected_edges=_intersected_edges(mask), cells=cells)

    # For general cases, use phi>0 polygon (complement of inside mask).
    # Invert the mask to get the outside region
    return _build_inside_case(mask ^ 0b1111)


def _build_interface_case(mask: int) -> ClipCase:
    edges = _intersected_edges(mask)
    cells: List[EmittedCell] = []

    if len(edges) == 0:
        return ClipCase(intersected_edges=edges, cells=tuple(cells))

    # Opposite corners -> two disjoint segments
    if mask == 0b0101:
        cells.append(EmittedCell(cell_type="interval", vertices=(("E", 0), ("E", 3))))
        cells.append(EmittedCell(cell_type="interval", vertices=(("E", 1), ("E", 2))))
        return ClipCase(intersected_edges=edges, cells=tuple(cells))

    if mask == 0b1010:
        cells.append(EmittedCell(cell_type="interval", vertices=(("E", 0), ("E", 1))))
        cells.append(EmittedCell(cell_type="interval", vertices=(("E", 2), ("E", 3))))
        return ClipCase(intersected_edges=edges, cells=tuple(cells))

    # All other cases: single segment connecting the two intersection edges in order
    if len(edges) == 2:
        cells.append(
            EmittedCell(
                cell_type="interval", vertices=(("E", edges[0]), ("E", edges[1]))
            )
        )
    else:
        # Should not happen for marching squares, fallback to empty
        cells = []

    return ClipCase(intersected_edges=edges, cells=tuple(cells))


def _build_quadrilateral_tables() -> ClipTables:
    inside_cases: List[ClipCase] = []
    outside_cases: List[ClipCase] = []
    interface_cases: List[ClipCase] = []

    for mask in range(16):
        inside_cases.append(_build_inside_case(mask))
        outside_cases.append(_build_outside_case(mask))
        interface_cases.append(_build_interface_case(mask))

    return ClipTables(
        inside=inside_cases, outside=outside_cases, interface=interface_cases
    )


# Register quad tables eagerly
register_tables("quadrilateral", _build_quadrilateral_tables())
