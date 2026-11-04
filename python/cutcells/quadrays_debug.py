# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Plots of the cells that quadrays writes to ``QuadraysOptions.debug_dir``.

A cell is written when it reaches a limit of the engine or exceeds
``QuadraysOptions.debug_points``::

    options = cutcells.QuadraysOptions()
    options.debug_dir = "bad_cells"
    part.quadrature(order=8, options=options)

and plotted, one image next to each file, with::

    python -m cutcells.quadrays_debug bad_cells/*.json

Needs matplotlib, which the module imports only to plot.
"""

import argparse
import json
from pathlib import Path

import numpy as np

# BoxFate in cpp/src/quadrays/engine.h
FATES = ("certified", "split", "bisected", "refined", "rotated", "uncertified", "abandoned")
_CERTIFIED, _SPLIT, _BISECTED, _REFINED, _ROTATED, _UNCERTIFIED, _ABANDONED = range(7)

_INK, _MUTED, _GRID, _SURFACE = "#0b0b0b", "#52514e", "#e7e6e2", "#fcfcfb"
_AXIS_ALIGNED, _IN_ROTATED_FRAME = "#2a78d6", "#eb6834"
_RAMP = ("#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#104281", "#0d366b")
_LINE_STYLES = ("-", "--", ":", "-.")


def read_report(path):
    """Read a cell report, with its arrays as NumPy arrays."""
    report = json.loads(Path(path).read_text())
    boxes = report["boxes"]
    for key in ("parent", "level", "depth", "height", "fate", "rotated", "points"):
        boxes[key] = np.asarray(boxes[key], dtype=int)
    for key in ("lo", "hi"):
        boxes[key] = np.asarray(boxes[key], dtype=float).reshape(-1, 3)
    boxes["frame"] = np.asarray(boxes["frame"], dtype=float).reshape(-1, 12)
    tdim, n = report["tdim"], report["grid"]
    report["values"] = [np.asarray(v, dtype=float).reshape((n,) * tdim) for v in report["values"]]
    report["points"] = np.asarray(report["points"], dtype=float).reshape(-1, 3)
    report["clips"] = np.asarray(report["clips"], dtype=float).reshape(-1, 4)
    return report


def top_leaves(report):
    """The indices of the cell's top boxes that nothing below replaces.

    A box integrated in a rotated frame counts as itself: the frame's own box
    only encloses it, with the box's faces as clip planes, and is left out.
    """
    b = report["boxes"]
    n = len(b["parent"])
    in_rotated = np.zeros(n, dtype=bool)
    for i in range(n):
        p = b["parent"][i]
        in_rotated[i] = p >= 0 and (in_rotated[p] or b["fate"][p] == _ROTATED)
    keep = (b["level"] == report["tdim"]) & ~in_rotated
    keep &= (b["fate"] != _BISECTED) & (b["fate"] != _REFINED) & (b["fate"] != _ABANDONED)
    return np.flatnonzero(keep)


def summary(report):
    """What the engine did in the cell, as text."""
    b = report["boxes"]
    lines = [f"cell {report['cell']}, {report['cell_type']}: {report['selection']}"]
    lines.append(f"order {report['order']}, q = {report['q']}: {report['n_points']:,} points")
    s = report["stats"]
    lines.append(f"{s['bisections']} bisections, {s['rotations']} rotations, {s['fold_splits']} fold splits, "
                 f"{s['uncertified']} uncertified")
    if s["points_limited"]:
        lines.append(f"over max_points = {report['options']['max_points']:,}: the boxes are those of the run stopped there")
    if s["bisections_limited"]:
        lines.append(f"used all of max_bisections = {report['options']['max_bisections']:,}")
    # a rotated box counts once: its frame's box below it is no box of its own
    frames = np.array([p >= 0 and b["fate"][p] == _ROTATED for p in b["parent"]], dtype=bool)
    for level in sorted(set(b["level"]), reverse=True):
        at_level = (b["level"] == level) & ~frames
        counts = [f"{np.count_nonzero(at_level & (b['fate'] == i))} {name}"
                  for i, name in enumerate(FATES) if np.any(at_level & (b["fate"] == i))]
        lines.append(f"level {level}: " + ", ".join(counts))
    leaves = top_leaves(report)
    if len(leaves):
        rotated = b["fate"][leaves] == _ROTATED
        total = max(1, b["points"][leaves].sum())
        lines.append(f"top boxes: {len(leaves)}, {np.count_nonzero(rotated)} in rotated frames with "
                     f"{100 * b['points'][leaves[rotated]].sum() / total:.0f} % of the points")
    return "\n".join(lines)


def _style(ax, title, labels):
    ax.set_title(title, fontsize=10, color=_INK, loc="left")
    ax.tick_params(colors=_MUTED, labelsize=8, length=3)
    for s in ax.spines.values():
        s.set_color(_GRID)
    ax.set_xlabel(labels[0], color=_MUTED, fontsize=9)
    ax.set_ylabel(labels[1], color=_MUTED, fontsize=9)
    ax.set_facecolor(_SURFACE)


def _zero_sets(ax, report, axes, normal, at):
    """The zero sets of the level sets in the slice through @p at, on the grid."""
    n = report["grid"]
    grid = np.linspace(0, 1, n)
    index = [slice(None)] * report["tdim"]
    if normal is not None:
        index[normal] = int(round(at[normal] * (n - 1)))
    for v, style in zip(report["values"], _LINE_STYLES):
        values = v[tuple(index)]
        if np.nanmin(values) < 0 < np.nanmax(values):
            ax.contour(grid, grid, values.T, levels=[0], colors=_INK, linewidths=1.2, linestyles=style)


def _clip_lines(ax, report, axes, normal, at):
    """The clip planes of the cell (the faces of a tet, say) in the slice."""
    t = np.linspace(0, 1, 2)
    for row in report["clips"]:
        c, d = row[:3], row[3]
        rhs = d - (c[normal] * at[normal] if normal is not None else 0.0)
        a, b = c[axes[0]], c[axes[1]]
        if abs(b) > 1e-12:
            ax.plot(t, (rhs - a * t) / b, color=_MUTED, lw=0.8)
        elif abs(a) > 1e-12:
            ax.plot([rhs / a, rhs / a], [0, 1], color=_MUTED, lw=0.8)


def _boxes_panel(ax, report, leaves, norm, cmap, axes, normal, at):
    """The top boxes in a slice, by their points; hatched in a rotated frame."""
    from matplotlib.patches import Rectangle

    b = report["boxes"]
    for i in leaves:
        if normal is not None and not b["lo"][i, normal] <= at[normal] <= b["hi"][i, normal]:
            continue
        lo, hi = b["lo"][i, axes], b["hi"][i, axes]
        points = b["points"][i]
        ax.add_patch(Rectangle(lo, *(hi - lo), facecolor=cmap(norm(points)) if points > 0 else "#f0efec",
                               edgecolor=_SURFACE, lw=0.6, hatch="////" if b["fate"][i] == _ROTATED else None))
    _zero_sets(ax, report, axes, normal, at)
    _clip_lines(ax, report, axes, normal, at)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")


def plot_report(report, path=None, at=None):
    """Plot a cell report: its top boxes and zero sets in slices through the
    heaviest box (or the point @p at of the unit box), where the points land,
    and how the points divide by box depth.

    Args:
        report: a report from read_report, or the path of one.
        path: where to save the image; None leaves the figure to the caller.
        at: the point the slices go through, in box coordinates.

    Returns:
        The matplotlib figure.
    """
    import matplotlib

    if path is not None:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, LogNorm
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    if not isinstance(report, dict):
        report = read_report(report)
    tdim, b = report["tdim"], report["boxes"]
    if tdim < 2:
        raise ValueError("quadrays_debug: an interval has no boxes to plot")
    leaves = top_leaves(report)
    points = b["points"][leaves] if len(leaves) else np.zeros(1, dtype=int)
    heaviest = leaves[np.argmax(points)] if len(leaves) else None
    if at is None:
        at = 0.5 * (b["lo"][heaviest] + b["hi"][heaviest]) if heaviest is not None else np.full(3, 0.5)
    at = np.asarray(at, dtype=float)

    cmap = LinearSegmentedColormap.from_list("points", _RAMP)
    norm = LogNorm(vmin=max(1, points[points > 0].min() if np.any(points > 0) else 1), vmax=max(2, points.max()))
    if tdim == 3:
        fig, axs = plt.subplots(2, 3, figsize=(15, 10), facecolor=_SURFACE, constrained_layout=True)
        slices = [((0, 1), 2), ((0, 2), 1), ((1, 2), 0)]
        slice_axes = axs[0]
        density, depths, text = axs[1]
    else:
        fig, axs = plt.subplots(2, 2, figsize=(11, 10), facecolor=_SURFACE, constrained_layout=True)
        slices = [((0, 1), None)]
        slice_axes = [axs[0, 0]]
        density, depths, text = axs[0, 1], axs[1, 0], axs[1, 1]

    for ax, (axes, normal) in zip(slice_axes, slices):
        axes = list(axes)
        _boxes_panel(ax, report, leaves, norm, cmap, axes, normal, at)
        where = "" if normal is None else f", u{normal} = {at[normal]:.3f}"
        _style(ax, f"Top boxes by their points{where}", (f"u{axes[0]}", f"u{axes[1]}"))
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=list(slice_axes), fraction=0.02, pad=0.01,
                 label="points in the box")

    # where the points land: a slab about the first slice
    axes, normal = list(slices[0][0]), slices[0][1]
    x = report["points"]
    if normal is not None:
        x = x[np.abs(x[:, normal] - at[normal]) < 0.03]
    if len(x):
        h, xe, ye = np.histogram2d(x[:, axes[0]], x[:, axes[1]], bins=200, range=[[0, 1], [0, 1]])
        density.pcolormesh(xe, ye, np.ma.masked_equal(h.T, 0), cmap=cmap, norm=LogNorm(vmin=1, vmax=max(2, h.max())),
                           rasterized=True)
    _zero_sets(density, report, axes, normal, at)
    _clip_lines(density, report, axes, normal, at)
    density.set_xlim(0, 1)
    density.set_ylim(0, 1)
    density.set_aspect("equal")
    slab = "" if normal is None else f" within 0.03 of u{normal} = {at[normal]:.3f}"
    _style(density, f"Where the points land{slab}", (f"u{axes[0]}", f"u{axes[1]}"))

    # the points by the depth of their box
    if len(leaves):
        rotated = b["fate"][leaves] == _ROTATED
        levels = np.arange(b["depth"][leaves].min(), b["depth"][leaves].max() + 1)
        plain = np.array([b["points"][leaves[~rotated & (b["depth"][leaves] == d)]].sum() for d in levels])
        turned = np.array([b["points"][leaves[rotated & (b["depth"][leaves] == d)]].sum() for d in levels])
        depths.bar(levels, plain, width=0.8, color=_AXIS_ALIGNED)
        depths.bar(levels, turned, width=0.8, bottom=plain, color=_IN_ROTATED_FRAME)
        depths.grid(axis="y", color=_GRID, lw=0.6)
        depths.set_axisbelow(True)
        depths.legend(handles=[Patch(color=_AXIS_ALIGNED, label="axis-aligned boxes"),
                               Patch(color=_IN_ROTATED_FRAME, label="boxes in rotated frames")],
                      frameon=False, fontsize=8, labelcolor=_MUTED)
    _style(depths, "Points by the depth of their top box", ("bisections", "points"))

    text.axis("off")
    text.text(0, 1, summary(report), va="top", ha="left", fontsize=8.5, color=_INK, wrap=True, family="monospace",
              transform=text.transAxes)

    handles = [Line2D([], [], color=_INK, lw=1.2, ls=style, label=name)
               for name, style in zip(report["level_sets"], _LINE_STYLES)]
    handles.append(Patch(facecolor="white", edgecolor=_MUTED, hatch="////", label="integrated in a rotated frame"))
    fig.legend(handles=handles, loc="outside lower center", ncol=len(handles), frameon=False, fontsize=9,
               labelcolor=_MUTED)
    fig.suptitle(f"quadrays: cell {report['cell']}, {report['n_points']:,} points", fontsize=13, color=_INK,
                 x=0.01, ha="left")
    if path is not None:
        fig.savefig(path, dpi=150, facecolor=_SURFACE)
        plt.close(fig)
    return fig


def main():
    parser = argparse.ArgumentParser(description="Plot the cells quadrays wrote to QuadraysOptions.debug_dir.")
    parser.add_argument("reports", nargs="+", help="the cells' JSON files")
    parser.add_argument("--at", type=float, nargs=3, metavar=("U0", "U1", "U2"),
                        help="the point the slices go through (default: the middle of the heaviest box)")
    parser.add_argument("--no-plot", action="store_true", help="print the summaries only")
    args = parser.parse_args()
    for path in args.reports:
        report = read_report(path)
        print(f"{path}\n{summary(report)}\n")
        if not args.no_plot and report["tdim"] >= 2:
            image = Path(path).with_suffix(".png")
            plot_report(report, image, args.at)
            print(f"wrote {image}\n")


if __name__ == "__main__":
    main()
