"""Adapters for VTK TableBasedClip tables.

Populate the stubs with canonical VTK connectivity to regenerate CutCells tables.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Literal, Sequence, Tuple

import re

from .topology import edge_endpoints

PointRef = Tuple[Literal["V", "E", "N"], int]


@dataclass(frozen=True)
class SpecialPoint:
    # Currently only VTK's centroid point N0 is expected.
    point_id: int
    refs: Tuple[PointRef, ...]


@dataclass(frozen=True)
class EmittedCell:
    cell_type: str
    vertices: Tuple[PointRef, ...]


@dataclass(frozen=True)
class ClipCase:
    intersected_edges: Tuple[int, ...]
    cells: Tuple[EmittedCell, ...]
    special_points: Tuple[SpecialPoint, ...] = ()


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


# --- Hexahedron TableBasedClip (VTK header parse) -------------------------


def _intersected_edges_from_mask(cell_type: str, mask: int) -> Tuple[int, ...]:
    edges = edge_endpoints(cell_type)
    out: List[int] = []
    for eid, (a, b) in enumerate(edges):
        sa = (mask >> a) & 1
        sb = (mask >> b) & 1
        if sa != sb:
            out.append(eid)
    return tuple(out)


def _vtk_name_to_point_ref(name: str) -> PointRef:
    # Vertices
    if re.fullmatch(r"P[0-7]", name):
        return ("V", int(name[1]))

    # Edges for hex: EA..EL correspond to 12 edges, contiguous.
    if re.fullmatch(r"E[A-L]", name):
        return ("E", ord(name[1]) - ord("A"))

    # Centroid
    if name == "N0":
        return ("N", 0)

    raise ValueError(f"Unsupported VTK point name: {name}")


def _vtk_st_to_cell_type(st: str) -> str:
    # We only need 3D volume cell types here.
    if st == "ST_TET":
        return "tetrahedron"
    if st == "ST_PYR":
        return "pyramid"
    if st == "ST_WDG":
        return "prism"
    if st == "ST_HEX":
        return "hexahedron"
    if st == "ST_PNT":
        return "__point__"
    # VTK tables sometimes contain lower-dimensional primitives (e.g. vertex-only
    # artifacts) that are not representable in CutCells' volume tables.
    if st in {"ST_VTX", "ST_LN", "ST_LIN", "ST_TRI", "ST_QUA", "ST_POLY"}:
        return "__skip__"
    raise ValueError(f"Unsupported VTK shape token: {st}")


def _faces_for_cell(cell_type: str) -> Tuple[Tuple[int, ...], ...]:
    # Face templates in VTK point ordering.
    if cell_type == "tetrahedron":
        return ((0, 1, 2), (0, 3, 1), (1, 3, 2), (0, 2, 3))
    if cell_type == "pyramid":
        return ((0, 1, 2, 3), (0, 4, 1), (1, 4, 2), (2, 4, 3), (3, 4, 0))
    if cell_type == "prism":
        return ((0, 1, 2), (3, 5, 4), (0, 3, 4, 1), (1, 4, 5, 2), (2, 5, 3, 0))
    if cell_type == "hexahedron":
        return (
            (0, 1, 2, 3),
            (4, 7, 6, 5),
            (0, 4, 5, 1),
            (1, 5, 6, 2),
            (2, 6, 7, 3),
            (3, 7, 4, 0),
        )
    raise ValueError(f"Unsupported cell type for face extraction: {cell_type}")


def _derive_interface_case_from_inside(inside: ClipCase) -> ClipCase:
    # Collect boundary faces from inside volume cells.
    face_counts: dict[Tuple[int, ...], int] = {}
    face_oriented: dict[Tuple[int, ...], Tuple[PointRef, ...]] = {}

    def key_for_face(face: Sequence[PointRef]) -> Tuple[int, ...]:
        # Stable key: sort by a numeric encoding of PointRef.
        def enc(ref: PointRef) -> int:
            kind, idx = ref
            if kind == "E":
                return idx
            if kind == "V":
                return 100 + idx
            if kind == "N":
                return 200 + idx
            raise ValueError(ref)

        return tuple(sorted(enc(r) for r in face))

    for cell in inside.cells:
        templ = _faces_for_cell(cell.cell_type)
        for f in templ:
            face = tuple(cell.vertices[i] for i in f)
            k = key_for_face(face)
            face_counts[k] = face_counts.get(k, 0) + 1
            face_oriented.setdefault(k, face)

    iface_cells: List[EmittedCell] = []
    for k, count in face_counts.items():
        if count != 1:
            continue
        face = face_oriented[k]

        # In the generic (phi!=0) case, interface faces should contain no original vertices.
        if any(kind == "V" for (kind, _) in face):
            continue

        if len(face) == 3:
            iface_cells.append(EmittedCell(cell_type="triangle", vertices=face))
        elif len(face) == 4:
            iface_cells.append(EmittedCell(cell_type="quadrilateral", vertices=face))

    # Only keep special points if referenced by interface faces.
    used_special: set[PointRef] = set()
    for c in iface_cells:
        for r in c.vertices:
            if r[0] == "N":
                used_special.add(r)
    special_points = tuple(
        sp for sp in inside.special_points if ("N", sp.point_id) in used_special
    )

    return ClipCase(
        intersected_edges=inside.intersected_edges,
        cells=tuple(iface_cells),
        special_points=special_points,
    )


def _build_hexahedron_tables() -> ClipTables:
    """Build hexahedron tables by parsing VTK's published case tables.

    This supports VTK's centroid point (ST_PNT -> N0) by storing its defining
    references in `ClipCase.special_points`.
    """

    from .vtk_tablebasedclip_parser import load_hexahedron_cases

    inside_cases: List[ClipCase] = []

    vtk_cases = load_hexahedron_cases()
    if len(vtk_cases) != 256:
        raise RuntimeError(f"Expected 256 hexahedron cases, got {len(vtk_cases)}")

    for case_id, vtk_case in enumerate(vtk_cases):
        if case_id != vtk_case.case_id:
            raise RuntimeError("Case ids not sequential")

        special_points: List[SpecialPoint] = []
        emitted: List[EmittedCell] = []

        for item in vtk_case.items:
            st = item[0]
            n = int(item[1])
            pts = item[2:]
            if len(pts) != n:
                raise RuntimeError("Malformed case item")

            cell_type = _vtk_st_to_cell_type(st)
            if cell_type == "__skip__":
                continue

            if cell_type == "__point__":
                # Define N0 as an unweighted average (duplicates => weights).
                refs = tuple(_vtk_name_to_point_ref(p) for p in pts)
                special_points.append(SpecialPoint(point_id=0, refs=refs))
                continue

            verts = tuple(_vtk_name_to_point_ref(p) for p in pts)
            emitted.append(EmittedCell(cell_type=cell_type, vertices=verts))

        # Derive intersected edges from the actual VTK case stream.
        #
        # Rationale: VTK's hexahedron case index does not appear to map 1:1 to a
        # simple vertex-bitmask in P0..P7 order (e.g. case 188 references EL, but
        # a naive mask-based edge test would miss it). By collecting all edge
        # references used by the emitted volume cells (and by special-point
        # definitions), we guarantee runtime will create every intersection token
        # that the tables may reference.
        used_edges: set[int] = set()
        for cell in emitted:
            for kind, idx in cell.vertices:
                if kind == "E":
                    used_edges.add(idx)
        for sp in special_points:
            for kind, idx in sp.refs:
                if kind == "E":
                    used_edges.add(idx)

        inside_cases.append(
            ClipCase(
                intersected_edges=tuple(sorted(used_edges)),
                cells=tuple(emitted),
                special_points=tuple(special_points),
            )
        )

    # Outside is just the inside of the complemented mask.
    outside_cases = [inside_cases[mask ^ 0xFF] for mask in range(256)]

    # Interface faces derived from inside volume boundary.
    interface_cases = [_derive_interface_case_from_inside(c) for c in inside_cases]

    return ClipTables(
        inside=inside_cases,
        outside=outside_cases,
        interface=interface_cases,
    )


def _build_prism_tables() -> ClipTables:
    """Build prism (VTK_WEDGE) tables by parsing VTK's published case tables.

    Like the hexahedron tables, this supports VTK's centroid point (ST_PNT -> N0)
    by storing its defining references in `ClipCase.special_points`.
    """

    from .vtk_tablebasedclip_parser import load_wedge_cases

    inside_cases: List[ClipCase] = []

    vtk_cases = load_wedge_cases()
    if len(vtk_cases) != 64:
        raise RuntimeError(f"Expected 64 wedge cases, got {len(vtk_cases)}")

    for case_id, vtk_case in enumerate(vtk_cases):
        if case_id != vtk_case.case_id:
            raise RuntimeError("Case ids not sequential")

        special_points: List[SpecialPoint] = []
        emitted: List[EmittedCell] = []

        for item in vtk_case.items:
            st = item[0]
            n = int(item[1])
            pts = item[2:]
            if len(pts) != n:
                raise RuntimeError("Malformed case item")

            cell_type = _vtk_st_to_cell_type(st)
            if cell_type == "__skip__":
                continue

            if cell_type == "__point__":
                refs = tuple(_vtk_name_to_point_ref(p) for p in pts)
                special_points.append(SpecialPoint(point_id=0, refs=refs))
                continue

            verts = tuple(_vtk_name_to_point_ref(p) for p in pts)
            emitted.append(EmittedCell(cell_type=cell_type, vertices=verts))

        used_edges: set[int] = set()
        for cell in emitted:
            for kind, idx in cell.vertices:
                if kind == "E":
                    used_edges.add(idx)
        for sp in special_points:
            for kind, idx in sp.refs:
                if kind == "E":
                    used_edges.add(idx)

        # Sanity: wedge has 9 edges (EA..EI). If we ever see larger edge ids,
        # it likely indicates we parsed the wrong table block.
        if any(eid >= 9 for eid in used_edges):
            raise RuntimeError(
                f"Wedge table references invalid edge id(s): {sorted(used_edges)}"
            )

        inside_cases.append(
            ClipCase(
                intersected_edges=tuple(sorted(used_edges)),
                cells=tuple(emitted),
                special_points=tuple(special_points),
            )
        )

    # Outside is just the inside of the complemented mask.
    outside_cases = [inside_cases[mask ^ 0x3F] for mask in range(64)]

    # Interface faces derived from inside volume boundary.
    interface_cases = [_derive_interface_case_from_inside(c) for c in inside_cases]

    return ClipTables(
        inside=inside_cases,
        outside=outside_cases,
        interface=interface_cases,
    )


def _build_pyramid_tables() -> ClipTables:
    """Build pyramid (VTK_PYRAMID) tables by parsing VTK's published case tables.

    Like the hexahedron/prism tables, this supports VTK's centroid point
    (ST_PNT -> N0) by storing its defining references in `ClipCase.special_points`.
    """

    from .vtk_tablebasedclip_parser import load_pyramid_cases

    inside_cases: List[ClipCase] = []

    vtk_cases = load_pyramid_cases()
    if len(vtk_cases) != 32:
        raise RuntimeError(f"Expected 32 pyramid cases, got {len(vtk_cases)}")

    for case_id, vtk_case in enumerate(vtk_cases):
        if case_id != vtk_case.case_id:
            raise RuntimeError("Case ids not sequential")

        special_points: List[SpecialPoint] = []
        emitted: List[EmittedCell] = []

        for item in vtk_case.items:
            st = item[0]
            n = int(item[1])
            pts = item[2:]
            if len(pts) != n:
                raise RuntimeError("Malformed case item")

            cell_type = _vtk_st_to_cell_type(st)
            if cell_type == "__skip__":
                continue

            if cell_type == "__point__":
                refs = tuple(_vtk_name_to_point_ref(p) for p in pts)
                special_points.append(SpecialPoint(point_id=0, refs=refs))
                continue

            verts = tuple(_vtk_name_to_point_ref(p) for p in pts)
            emitted.append(EmittedCell(cell_type=cell_type, vertices=verts))

        used_edges: set[int] = set()
        for cell in emitted:
            for kind, idx in cell.vertices:
                if kind == "E":
                    used_edges.add(idx)
        for sp in special_points:
            for kind, idx in sp.refs:
                if kind == "E":
                    used_edges.add(idx)

        # Sanity: pyramid has 8 edges (EA..EH)
        if any(eid >= 8 for eid in used_edges):
            raise RuntimeError(
                f"Pyramid table references invalid edge id(s): {sorted(used_edges)}"
            )

        inside_cases.append(
            ClipCase(
                intersected_edges=tuple(sorted(used_edges)),
                cells=tuple(emitted),
                special_points=tuple(special_points),
            )
        )

    # Outside is just the inside of the complemented mask.
    outside_cases = [inside_cases[mask ^ 0x1F] for mask in range(32)]

    # Interface faces derived from inside volume boundary.
    interface_cases = [_derive_interface_case_from_inside(c) for c in inside_cases]

    return ClipTables(
        inside=inside_cases,
        outside=outside_cases,
        interface=interface_cases,
    )


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

# Register hex tables eagerly (requires VTK at runtime)
register_tables("hexahedron", _build_hexahedron_tables())

# Register prism (wedge) tables eagerly (requires VTK at runtime)
register_tables("prism", _build_prism_tables())

# Register pyramid tables eagerly (requires VTK at runtime)
register_tables("pyramid", _build_pyramid_tables())
