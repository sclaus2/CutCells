"""Plot marching-squares cases for quadrilateral clipping.

The two ambiguous opposite-corner masks (0101 and 1010) are shown as *two*
separate subplots (one per asymptotic-decider variant), to avoid overplotting.

Run from repo root (fenicsx-dev env):
    conda run -n fenicsx-dev python python/demo/plot_quadrilateral_cases.py [--save out.png]
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

# Use the tablegen package directly
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tablegen"))

from cutcells_tablegen.topology import edge_endpoints  # noqa: E402
from cutcells_tablegen.vtk_tables_import import ClipCase, EmittedCell, get_cases  # noqa: E402

try:
    import cutcells  # type: ignore

    _HAS_CUTCELLS = True
except Exception:  # pragma: no cover
    _HAS_CUTCELLS = False

Point = Tuple[float, float]


VTK_LINE = 3
VTK_TRIANGLE = 5
VTK_QUAD = 9


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


def _iter_connectivity(connectivity: Sequence[int]) -> Iterable[List[int]]:
    i = 0
    n = len(connectivity)
    while i < n:
        nv = int(connectivity[i])
        verts = [int(connectivity[i + j]) for j in range(1, 1 + nv)]
        yield verts
        i += 1 + nv


def _plot_runtime_cutcell(
    ax,
    cut_cell,
    facecolor: str,
    edgecolor: str,
    alpha: float,
    linestyle: str,
    label: str,
):
    coords = list(cut_cell.vertex_coords)
    pts = [(coords[i], coords[i + 1]) for i in range(0, len(coords), 3)]
    connectivity = list(cut_cell.connectivity)
    types = list(cut_cell.types)

    for cell_type, verts in zip(types, _iter_connectivity(connectivity)):
        poly = [pts[v] for v in verts]
        xs, ys = zip(*poly)

        if cell_type == VTK_LINE:
            ax.plot(
                xs, ys, color=edgecolor, linewidth=2.0, linestyle=linestyle, label=label
            )
        else:
            ax.fill(
                xs,
                ys,
                facecolor=facecolor,
                alpha=alpha,
                edgecolor=edgecolor,
                linewidth=0.8,
                linestyle=linestyle,
                label=label,
            )
        label = "_nolegend_"


def _ambiguous_ls(mask: int, variant: int) -> List[float]:
    """Return ls values for an ambiguous mask with a chosen decider variant.

    Variant 0: d = f0*f2 - f1*f3 >= 0
    Variant 1: d < 0
    """
    if mask not in (0b0101, 0b1010):
        raise ValueError("mask must be 0101 or 1010")
    if variant not in (0, 1):
        raise ValueError("variant must be 0 or 1")

    # Set signs per mask
    if mask == 0b0101:
        # (-, +, -, +)
        if variant == 0:
            return [-2.0, 1.0, -2.0, 1.0]
        return [-1.0, 2.0, -1.0, 2.0]

    # mask == 0b1010 -> (+, -, +, -)
    if variant == 0:
        return [2.0, -1.0, 2.0, -1.0]
    return [1.0, -2.0, 1.0, -2.0]


def plot_cases(save: Path | None = None) -> None:
    inside_cases: List[ClipCase] = get_cases("quadrilateral", "inside")
    outside_cases: List[ClipCase] = get_cases("quadrilateral", "outside")
    interface_cases: List[ClipCase] = get_cases("quadrilateral", "interface")

    @dataclass(frozen=True)
    class Panel:
        kind: str  # "tablegen" | "runtime"
        mask: int
        variant: int | None = None

    panels: List[Panel] = []
    for mask in range(16):
        if mask in (0b0101, 0b1010):
            panels.append(Panel(kind="runtime", mask=mask, variant=0))
            panels.append(Panel(kind="runtime", mask=mask, variant=1))
        else:
            panels.append(Panel(kind="tablegen", mask=mask))

    n = len(panels)  # 14 + 4 = 18
    ncols = 6
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(3.0 * ncols, 3.0 * nrows),
        constrained_layout=True,
    )
    axes_flat = list(axes.ravel())

    for i, panel in enumerate(panels):
        ax = axes_flat[i]
        mask = panel.mask
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

        show_legend_labels = i == 0
        outside_label = "outside" if show_legend_labels else "_nolegend_"
        inside_label = "inside" if show_legend_labels else "_nolegend_"
        iface_label = "interface" if show_legend_labels else "_nolegend_"

        if panel.kind == "tablegen":
            # Draw outside first (so inside can overlay)
            draw_cells(
                ax,
                outside_cases[mask].cells,
                verts,
                ip,
                color="red",
                label=outside_label,
            )
            draw_cells(
                ax,
                inside_cases[mask].cells,
                verts,
                ip,
                color="blue",
                label=inside_label,
            )
            draw_segments(
                ax,
                interface_cases[mask].cells,
                verts,
                ip,
                color="black",
                label=iface_label,
            )
            ax.set_title(f"mask {mask:04b}")
        else:
            if not _HAS_CUTCELLS:
                ax.set_title(f"mask {mask:04b} (cutcells unavailable)")
            else:
                assert panel.variant is not None
                vertex_coordinates = np.asarray(
                    [0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0], dtype=np.float64
                )
                cell_type = cutcells.CellType.quadrilateral
                gdim = 2
                ls = np.asarray(_ambiguous_ls(mask, panel.variant), dtype=np.float64)

                # Keep quads as quads for easier interpretation
                inside = cutcells.cut(
                    cell_type, vertex_coordinates, gdim, ls, "phi<0", False
                )
                outside = cutcells.cut(
                    cell_type, vertex_coordinates, gdim, ls, "phi>0", False
                )
                iface = cutcells.cut(
                    cell_type, vertex_coordinates, gdim, ls, "phi=0", False
                )

                _plot_runtime_cutcell(
                    ax,
                    outside,
                    facecolor="red",
                    edgecolor="k",
                    alpha=0.35,
                    linestyle="-",
                    label=outside_label,
                )
                _plot_runtime_cutcell(
                    ax,
                    inside,
                    facecolor="blue",
                    edgecolor="k",
                    alpha=0.35,
                    linestyle="-",
                    label=inside_label,
                )
                _plot_runtime_cutcell(
                    ax,
                    iface,
                    facecolor="none",
                    edgecolor="black",
                    alpha=0.0,
                    linestyle="-",
                    label=iface_label,
                )

                ax.set_title(f"mask {mask:04b} (v{panel.variant})")
        ax.set_aspect("equal")
        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-0.1, 1.1)
        ax.axis("off")

    # Hide any unused axes
    for j in range(n, len(axes_flat)):
        axes_flat[j].axis("off")

    handles, labels = axes_flat[0].get_legend_handles_labels()
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
