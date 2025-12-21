"""Plot all 16 marching-squares cases for quadrilateral clipping.

Run from repo root (fenicsx-dev env):
  PYTHONPATH=tablegen python python/demo/plot_quadrilateral_cases.py [--save out.png]
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import matplotlib.pyplot as plt

# Use the tablegen package directly
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tablegen"))

from cutcells_tablegen.topology import edge_endpoints  # noqa: E402
from cutcells_tablegen.vtk_tables_import import ClipCase, EmittedCell, get_cases  # noqa: E402

Point = Tuple[float, float]


def build_vertices(mask: int) -> Tuple[List[Point], List[Point]]:
    """Return base vertices and edge intersection points for a given mask.

    Mask bit i=1 means vertex i is inside (phi<0). We assign ls=-1 inside, +1 outside,
    so intersections land at the midpoint of sign-changing edges.
    """
    verts: List[Point] = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    edges = edge_endpoints("quadrilateral")
    ls = [-1.0 if (mask >> i) & 1 else 1.0 for i in range(4)]

    ip: List[Point] = [(float("nan"), float("nan")) for _ in range(len(edges))]
    for e, (v0, v1) in enumerate(edges):
        if (ls[v0] < 0) != (ls[v1] < 0):
            t = abs(ls[v0]) / (abs(ls[v0]) + abs(ls[v1]))  # with +-1 -> 0.5
            x = verts[v0][0] * (1 - t) + verts[v1][0] * t
            y = verts[v0][1] * (1 - t) + verts[v1][1] * t
            ip[e] = (x, y)
    return verts, ip


def ref_to_point(
    ref: Tuple[str, int], verts: Sequence[Point], ip: Sequence[Point]
) -> Point:
    kind, idx = ref
    if kind == "V":
        return verts[idx]
    if kind == "E":
        return ip[idx]
    raise ValueError(f"Unknown ref {ref}")


def draw_cells(
    ax,
    cells: Iterable[EmittedCell],
    verts: Sequence[Point],
    ip: Sequence[Point],
    color: str,
    label: str,
):
    for cell in cells:
        coords = [ref_to_point(r, verts, ip) for r in cell.vertices]
        xs, ys = zip(*coords)
        ax.fill(
            xs,
            ys,
            facecolor=color,
            alpha=0.5,
            edgecolor="k",
            linewidth=0.5,
            label=label,
        )
        label = "_nolegend_"  # only label first patch


def draw_segments(
    ax,
    cells: Iterable[EmittedCell],
    verts: Sequence[Point],
    ip: Sequence[Point],
    color: str,
    label: str,
):
    for cell in cells:
        coords = [ref_to_point(r, verts, ip) for r in cell.vertices]
        xs, ys = zip(*coords)
        ax.plot(xs, ys, color=color, linewidth=2.0, label=label)
        label = "_nolegend_"


def plot_cases(save: Path | None = None) -> None:
    inside_cases: List[ClipCase] = get_cases("quadrilateral", "inside")
    outside_cases: List[ClipCase] = get_cases("quadrilateral", "outside")
    interface_cases: List[ClipCase] = get_cases("quadrilateral", "interface")

    fig, axes = plt.subplots(4, 4, figsize=(12, 12), constrained_layout=True)

    for mask in range(16):
        ax = axes[mask // 4][mask % 4]
        verts, ip = build_vertices(mask)

        # Draw original quad boundary
        quad_x = [verts[0][0], verts[1][0], verts[2][0], verts[3][0], verts[0][0]]
        quad_y = [verts[0][1], verts[1][1], verts[2][1], verts[3][1], verts[0][1]]
        ax.plot(
            quad_x,
            quad_y,
            "k--",
            linewidth=1.5,
            alpha=0.3,
            label="quad" if mask == 0 else "_nolegend_",
        )

        # Draw outside first (so inside can overlay)
        draw_cells(
            ax, outside_cases[mask].cells, verts, ip, color="red", label="outside"
        )
        draw_cells(
            ax, inside_cases[mask].cells, verts, ip, color="blue", label="inside"
        )
        draw_segments(
            ax, interface_cases[mask].cells, verts, ip, color="black", label="interface"
        )

        ax.set_title(f"mask {mask:04b}")
        ax.set_aspect("equal")
        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-0.1, 1.1)
        ax.axis("off")

    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=4)

    if save:
        fig.savefig(save, dpi=200)
    else:
        plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot all quadrilateral clip cases")
    parser.add_argument(
        "--save",
        type=Path,
        default=None,
        help="Path to save the figure instead of showing",
    )
    args = parser.parse_args()
    plot_cases(args.save)


if __name__ == "__main__":
    main()
