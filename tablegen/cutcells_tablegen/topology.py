"""VTK-consistent topology helpers for CutCells table generation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple


@dataclass(frozen=True)
class CellTopology:
    name: str
    num_vertices: int
    edges: Tuple[Tuple[int, int], ...]


# VTK edge numbering (matches cut_tetrahedron conventions already used elsewhere).
CELL_TOPOLOGY = {
    "quadrilateral": CellTopology(
        name="quadrilateral",
        num_vertices=4,
        edges=((0, 1), (1, 2), (2, 3), (3, 0)),
    ),
    "hexahedron": CellTopology(
        name="hexahedron",
        num_vertices=8,
        edges=(
            (0, 1),  # 0
            (1, 2),  # 1
            (2, 3),  # 2
            (3, 0),  # 3
            (4, 5),  # 4
            (5, 6),  # 5
            (6, 7),  # 6
            (7, 4),  # 7
            (0, 4),  # 8
            (1, 5),  # 9
            # NOTE: VTK's TableBasedClipCases.h uses (3,7) as edge 10 and (2,6) as edge 11
            # for VTK_HEXAHEDRON (see vtkTableBasedClipCasesBase::CellEdges).
            (3, 7),  # 10
            (2, 6),  # 11
        ),
    ),
    "prism": CellTopology(
        name="prism",
        num_vertices=6,
        edges=(
            (0, 1),  # 0
            (1, 2),  # 1
            (2, 0),  # 2
            (3, 4),  # 3
            (4, 5),  # 4
            (5, 3),  # 5
            (0, 3),  # 6
            (1, 4),  # 7
            (2, 5),  # 8
        ),
    ),
    "pyramid": CellTopology(
        name="pyramid",
        num_vertices=5,
        edges=(
            (0, 1),  # 0
            (1, 2),  # 1
            (2, 3),  # 2
            (3, 0),  # 3
            (0, 4),  # 4
            (1, 4),  # 5
            (2, 4),  # 6
            (3, 4),  # 7
        ),
    ),
}


def get_topology(cell_type: str) -> CellTopology:
    topo = CELL_TOPOLOGY.get(cell_type)
    if topo is None:
        raise ValueError(f"Unsupported cell type: {cell_type}")
    return topo


def edge_endpoints(cell_type: str) -> Tuple[Tuple[int, int], ...]:
    return get_topology(cell_type).edges


def num_vertices(cell_type: str) -> int:
    return get_topology(cell_type).num_vertices


def num_edges(cell_type: str) -> int:
    return len(get_topology(cell_type).edges)


def edge_count(cell_type: str) -> int:
    return num_edges(cell_type)


def supported_cell_types() -> List[str]:
    return list(CELL_TOPOLOGY.keys())
