# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT

"""Generate the figures and interactive views of the tutorial pages.

The figures are the views of the demos in ``python/demo``: each demo's
``plot(state, screenshot)`` renders its PyVista scene off screen. The
interactive views are small self-contained HTML pages that draw the
surfaces of the same parts on a canvas.

Run from the repository root, with CutCells and PyVista importable:

    python docs/tutorials/make_pyvista_scenes.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pyvista as pv

import cutcells

ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = ROOT / "docs" / "_static" / "tutorials"
sys.path.insert(0, str(ROOT / "python" / "demo"))

import demo_convergence  # noqa: E402
import demo_disk  # noqa: E402
import demo_gyroid  # noqa: E402
import demo_napkin_ring  # noqa: E402
import demo_sphere  # noqa: E402
import demo_sphere_plane  # noqa: E402

pv.OFF_SCREEN = True


# ============================================================================
# Interactive views
# ============================================================================


def _surface(part, mode: str = "cut_only", subdivision: int = 0) -> pv.PolyData:
    grid = cutcells.to_pyvista(part.visualization_mesh(mode=mode, degree=3))
    return grid.extract_surface(nonlinear_subdivision=subdivision, algorithm="dataset_surface").triangulate()


def _curve(part) -> list[list[list[float]]]:
    """The curve of a part of dimension 1 as polylines, one per piece."""
    grid = cutcells.to_pyvista(part.visualization_mesh(mode="cut_only", degree=3))
    lines = grid.extract_surface(nonlinear_subdivision=4, algorithm="dataset_surface")
    paths = []
    for cell in range(lines.n_cells):
        ids = lines.get_cell(cell).point_ids
        paths.append(lines.points[ids])
    return paths


def _view_data(surfaces: list[tuple[pv.PolyData, str]], paths: list[tuple[list, str]]) -> dict:
    """Points, triangles per color and paths of a scene, as flat rounded arrays."""
    points, groups = [], []
    for surface, color in surfaces:
        offset = len(points)
        points.extend(surface.points.tolist())
        faces = surface.faces.reshape(-1, 4)[:, 1:]
        # Slivers have no reliable normal and would shade as dark lines.
        x = surface.points[faces]
        area = 0.5 * np.linalg.norm(np.cross(x[:, 1] - x[:, 0], x[:, 2] - x[:, 0]), axis=1)
        faces = faces[area > 1e-3 * np.median(area)] + offset
        groups.append({"color": color, "triangles": faces.ravel().tolist()})
    allp = np.asarray(points)
    bounds = np.column_stack([allp.min(axis=0), allp.max(axis=0)]).ravel().tolist()
    return {
        "points": np.round(allp, 3).ravel().tolist(),
        "groups": groups,
        "paths": [{"points": np.round(p, 3).ravel().tolist(), "color": c, "width": 2.6}
                  for group, c in paths for p in group],
        "bounds": bounds,
    }


def _write_interactive_view(path: Path, title: str, data: dict) -> None:
    """A canvas view to turn and zoom, after the CutFEMx tutorial views."""
    data_json = json.dumps(data, separators=(",", ":"))
    template = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>__TITLE__</title>
  <style>
    :root { color-scheme: light; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }
    html, body { width: 100%; height: 100%; margin: 0; overflow: hidden; background: #ffffff; }
    canvas { width: 100%; height: 100%; display: block; cursor: grab; touch-action: none; background: #ffffff; }
    canvas.dragging { cursor: grabbing; }
    .toolbar { position: absolute; top: 10px; right: 12px; z-index: 2; display: flex; gap: 8px; }
    .hint { position: absolute; left: 12px; bottom: 10px; color: #64748b; font-size: 12px; }
    button { border: 1px solid #cbd5e1; background: rgba(255, 255, 255, 0.92); color: #334155;
             border-radius: 4px; padding: 6px 10px; font: inherit; cursor: pointer;
             box-shadow: 0 1px 3px rgba(15, 23, 42, 0.12); }
    button:hover { background: #f8fafc; }
  </style>
</head>
<body>
  <div class="toolbar"><button type="button" id="reset">Reset</button></div>
  <div class="hint">Drag to turn, scroll to zoom</div>
  <canvas id="plot" aria-label="__TITLE__"></canvas>
  <script>
    const data = __DATA__;
    const canvas = document.getElementById("plot");
    const ctx = canvas.getContext("2d");
    const initial = { yaw: -0.75, pitch: -1.15, zoom: 1.0 };  // z up, seen from above
    const view = { ...initial };
    const b = data.bounds;
    const center = [0.5 * (b[0] + b[1]), 0.5 * (b[2] + b[3]), 0.5 * (b[4] + b[5])];
    const span = Math.max(b[1] - b[0], b[3] - b[2], b[5] - b[4]);
    const rgb = {};
    let dragging = null;
    const points = [];
    for (let i = 0; i < data.points.length; i += 3) points.push(data.points.slice(i, i + 3));
    const cells = [];
    for (const group of data.groups) {
      const t = group.triangles;
      for (let i = 0; i < t.length; i += 3) cells.push({ ids: [t[i], t[i + 1], t[i + 2]], color: group.color });
    }
    const paths = data.paths.map((path) => {
      const q = [];
      for (let i = 0; i < path.points.length; i += 3) q.push(path.points.slice(i, i + 3));
      return { ...path, points: q };
    });

    function color(hex, shade) {
      if (!rgb[hex]) rgb[hex] = [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16));
      return `rgb(${rgb[hex].map((v) => Math.min(255, Math.round(v * shade))).join(",")})`;
    }

    function rotate(p) {
      const x0 = p[0] - center[0], y0 = p[1] - center[1], z0 = p[2] - center[2];
      const cy = Math.cos(view.yaw), sy = Math.sin(view.yaw);
      const cp = Math.cos(view.pitch), sp = Math.sin(view.pitch);
      const x1 = cy * x0 - sy * y0, y1 = sy * x0 + cy * y0;
      return [x1, cp * y1 - sp * z0, sp * y1 + cp * z0];
    }

    function resize() {
      const dpr = window.devicePixelRatio || 1;
      const rect = canvas.getBoundingClientRect();
      canvas.width = Math.max(1, Math.floor(rect.width * dpr));
      canvas.height = Math.max(1, Math.floor(rect.height * dpr));
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      render();
    }

    function render() {
      const rect = canvas.getBoundingClientRect();
      const w = rect.width, h = rect.height;
      ctx.fillStyle = "#ffffff";
      ctx.fillRect(0, 0, w, h);
      const scale = view.zoom * Math.min(w, h) / (span * 1.9);
      const r = points.map(rotate);
      const s = r.map((p) => [0.5 * w + scale * p[0], 0.5 * h - scale * p[1]]);
      // Triangles and path segments are drawn back to front together.
      const items = cells.map((c) => [(r[c.ids[0]][2] + r[c.ids[1]][2] + r[c.ids[2]][2]) / 3.0, c, null]);
      for (const path of paths) {
        const q = path.points.map(rotate);
        for (let i = 1; i < q.length; i++) items.push([0.5 * (q[i - 1][2] + q[i][2]) + 1e-3 * span, null, [q[i - 1], q[i], path]]);
      }
      items.sort((u, v) => u[0] - v[0]);
      ctx.lineJoin = "round";
      ctx.lineCap = "round";
      for (const [, cell, segment] of items) {
        if (segment) {
          const [p0, p1, path] = segment;
          ctx.beginPath();
          ctx.moveTo(0.5 * w + scale * p0[0], 0.5 * h - scale * p0[1]);
          ctx.lineTo(0.5 * w + scale * p1[0], 0.5 * h - scale * p1[1]);
          ctx.strokeStyle = path.color;
          ctx.lineWidth = path.width;
          ctx.stroke();
          continue;
        }
        const [a, bb, cc] = cell.ids;
        const u = [r[bb][0] - r[a][0], r[bb][1] - r[a][1], r[bb][2] - r[a][2]];
        const v = [r[cc][0] - r[a][0], r[cc][1] - r[a][1], r[cc][2] - r[a][2]];
        const n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        const light = Math.abs(-0.3 * n[0] + 0.45 * n[1] + 0.84 * n[2]) / (Math.hypot(...n) || 1.0);
        ctx.beginPath();
        ctx.moveTo(s[a][0], s[a][1]);
        ctx.lineTo(s[bb][0], s[bb][1]);
        ctx.lineTo(s[cc][0], s[cc][1]);
        ctx.closePath();
        ctx.fillStyle = color(cell.color, 0.55 + 0.5 * light);
        ctx.fill();
        // A stroke in the fill colour covers the hairline seams between pieces.
        ctx.strokeStyle = ctx.fillStyle;
        ctx.lineWidth = 1.1;
        ctx.stroke();
      }
    }

    canvas.addEventListener("pointerdown", (e) => {
      canvas.setPointerCapture(e.pointerId);
      canvas.classList.add("dragging");
      dragging = { x: e.clientX, y: e.clientY, yaw: view.yaw, pitch: view.pitch };
    });
    canvas.addEventListener("pointermove", (e) => {
      if (!dragging) return;
      view.yaw = dragging.yaw + 0.01 * (e.clientX - dragging.x);
      view.pitch = Math.max(-3.1, Math.min(0.2, dragging.pitch - 0.01 * (e.clientY - dragging.y)));
      render();
    });
    const stop = () => { dragging = null; canvas.classList.remove("dragging"); };
    canvas.addEventListener("pointerup", stop);
    canvas.addEventListener("pointercancel", stop);
    canvas.addEventListener("wheel", (e) => {
      e.preventDefault();
      view.zoom = Math.max(0.4, Math.min(4.0, view.zoom * (e.deltaY < 0 ? 1.1 : 0.9)));
      render();
    }, { passive: false });
    document.getElementById("reset").addEventListener("click", () => { Object.assign(view, initial); render(); });
    window.addEventListener("resize", resize);
    resize();
  </script>
</body>
</html>
"""
    path.write_text(template.replace("__TITLE__", title).replace("__DATA__", data_json), encoding="utf-8")


# ============================================================================
# Scenes
# ============================================================================

ORANGE = "#f59e42"
BLUE = "#2563eb"
MAGENTA = "#be185d"
GREEN = "#68b984"


def disk_scenes() -> None:
    demo_disk.plot(demo_disk.build(n=17), OUTPUT_DIR / "disk-scene.png")


def convergence_scenes() -> None:
    table = demo_convergence.errors([9, 17, 33, 65])
    demo_convergence.plot(table, OUTPUT_DIR / "convergence-plot.png")


def sphere_scenes() -> None:
    demo_sphere.plot(demo_sphere.build(), OUTPUT_DIR / "sphere-scene.png")


def sphere_plane_scenes() -> None:
    state = demo_sphere_plane.build()
    demo_sphere_plane.plot(state, OUTPUT_DIR / "sphere-plane-scene.png")
    parts = state["parts"]
    data = _view_data(
        [(_surface(parts["sphere below the plane"][0]), ORANGE), (_surface(parts["disk"][0]), BLUE)],
        [(_curve(parts["circle"][0]), MAGENTA)],
    )
    _write_interactive_view(OUTPUT_DIR / "sphere-plane-view.html", "Ball and plane", data)


def napkin_ring_scenes() -> None:
    state = demo_napkin_ring.build()
    demo_napkin_ring.plot(state, OUTPUT_DIR / "napkin-ring-scene.png")
    data = _view_data(
        [(_surface(state["zone"]), ORANGE), (_surface(state["hole"]), BLUE)],
        [(_curve(state["rims"]), MAGENTA)],
    )
    _write_interactive_view(OUTPUT_DIR / "napkin-ring-view.html", "Napkin ring", data)


def gyroid_scenes() -> None:
    state = demo_gyroid.build()
    demo_gyroid.plot(state, OUTPUT_DIR / "gyroid-scene.png")
    lower, upper = state["sheet_faces"]
    data = _view_data([(_surface(lower), BLUE), (_surface(upper), ORANGE)], [])
    _write_interactive_view(OUTPUT_DIR / "gyroid-view.html", "Gyroid sheet", data)


SCENES = {
    "disk": disk_scenes,
    "convergence": convergence_scenes,
    "sphere": sphere_scenes,
    "sphere-plane": sphere_plane_scenes,
    "napkin-ring": napkin_ring_scenes,
    "gyroid": gyroid_scenes,
}


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    names = sys.argv[1:] or list(SCENES)
    for name in names:
        print(f"rendering {name} ...", flush=True)
        SCENES[name]()


if __name__ == "__main__":
    main()
