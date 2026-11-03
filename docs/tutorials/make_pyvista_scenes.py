# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT

"""Generate the figures and interactive views of the tutorial pages.

The figures are the views of the demos in ``python/demo``: each demo's
``plot(state, screenshot)`` renders its PyVista scene off screen. The
interactive views are small self-contained HTML pages that draw the
surfaces of the same parts with WebGL.

Run from the repository root, with CutCells and PyVista importable:

    python docs/tutorials/make_pyvista_scenes.py
"""

from __future__ import annotations

import base64
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


def _surface(part, mode: str = "cut_only", subdivision: int = 1) -> pv.PolyData:
    grid = cutcells.to_pyvista(part.visualization_mesh(mode=mode, degree=3))
    return grid.extract_surface(nonlinear_subdivision=subdivision, algorithm="dataset_surface").triangulate()


def _curve(part, radius: float) -> pv.PolyData:
    """The curve of a part of dimension 1 as a tube, one piece per leaf."""
    grid = cutcells.to_pyvista(part.visualization_mesh(mode="cut_only", degree=3))
    lines = grid.extract_surface(nonlinear_subdivision=2, algorithm="dataset_surface")
    return lines.tube(radius=radius, n_sides=8).triangulate()


def _smooth(surface: pv.PolyData, grid: float = 1e-4):
    """Shared vertices and normals of a triangulated surface.

    Leaves do not share nodes, so points that agree to ``grid`` are merged.
    Normals are area-weighted sums of the faces' normals, each turned to agree
    with the sum so far: leaves need not be oriented alike, and the view lights
    both sides.
    """
    faces = surface.faces.reshape(-1, 4)[:, 1:]
    keys, index = np.unique(np.round(surface.points / grid).astype(np.int64), axis=0, return_inverse=True)
    index = index.ravel()
    points = np.zeros((len(keys), 3))
    points[index] = surface.points
    tri = index[faces]
    x = points[tri]
    n = np.cross(x[:, 1] - x[:, 0], x[:, 2] - x[:, 0])
    area = np.linalg.norm(n, axis=1)
    keep = (area > 1e-6 * np.median(area)) & (tri[:, 0] != tri[:, 1]) & (tri[:, 1] != tri[:, 2]) & (tri[:, 0] != tri[:, 2])
    tri, n = tri[keep], n[keep]
    normals = np.zeros_like(points)
    for corner in range(3):
        for _ in range(2):  # a second pass lets late faces follow the settled sum
            sign = np.sign(np.einsum("ij,ij->i", normals[tri[:, corner]], n))
            sign[sign == 0] = 1.0
            np.add.at(normals, tri[:, corner], sign[:, None] * n * 0.5)
    length = np.linalg.norm(normals, axis=1)
    normals /= np.where(length > 0, length, 1.0)[:, None]
    return points, normals, tri


def _view_data(groups: list[tuple[pv.PolyData, str]]) -> dict:
    """Vertices, normals and triangles per colour, packed as base64 arrays.

    Points are 16-bit fractions of the bounding box, normals 8-bit, triangle
    indices 16-bit where they fit and 32-bit otherwise, all little-endian.
    """
    smooth = [(_smooth(surface), color) for surface, color in groups]
    allp = np.vstack([points for (points, _, _), _ in smooth])
    lo, hi = allp.min(axis=0), allp.max(axis=0)
    size = np.where(hi > lo, hi - lo, 1.0)
    pack = lambda a: base64.b64encode(np.ascontiguousarray(a).tobytes()).decode("ascii")
    out = []
    for (points, normals, tri), color in smooth:
        out.append({
            "color": color,
            "points": pack(np.round((points - lo) / size * 65535).astype("<u2")),
            "normals": pack(np.round(normals * 127).astype("i1")),
            "triangles": pack(tri.astype("<u2" if len(points) < 65536 else "<u4")),
            "wide": len(points) >= 65536,
        })
    return {"groups": out, "lo": lo.tolist(), "size": size.tolist()}


def _write_interactive_view(path: Path, title: str, data: dict) -> None:
    """A WebGL view to turn and zoom, with a depth buffer and smooth shading."""
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
    .error { position: absolute; inset: 0; display: none; align-items: center; justify-content: center; color: #334155; }
    button { border: 1px solid #cbd5e1; background: rgba(255, 255, 255, 0.92); color: #334155;
             border-radius: 4px; padding: 6px 10px; font: inherit; cursor: pointer;
             box-shadow: 0 1px 3px rgba(15, 23, 42, 0.12); }
    button:hover { background: #f8fafc; }
  </style>
</head>
<body>
  <div class="toolbar"><button type="button" id="reset">Reset</button></div>
  <div class="hint">Drag to turn, scroll to zoom</div>
  <div class="error" id="error">This view needs WebGL.</div>
  <canvas id="plot" aria-label="__TITLE__"></canvas>
  <script>
    const data = __DATA__;
    const canvas = document.getElementById("plot");
    const gl = canvas.getContext("webgl", { antialias: true });
    const initial = { yaw: -0.75, pitch: -1.15, zoom: 1.0 };  // z up, seen from above
    const view = { ...initial };
    if (!gl) {
      document.getElementById("error").style.display = "flex";
    } else {
      gl.getExtension("OES_element_index_uint");
      const vs = `attribute vec3 p; attribute vec3 n; uniform mat3 R; uniform vec3 c; uniform vec2 s; uniform float dz;
        varying vec3 vn; void main() { vec3 q = R * (p - c); vn = R * n;
        gl_Position = vec4(q.x * s.x, q.y * s.y, -q.z * dz, 1.0); }`;
      const fs = `precision mediump float; varying vec3 vn; uniform vec3 color;
        void main() { vec3 m = normalize(vn); vec3 l = normalize(vec3(-0.3, 0.45, 0.84));
        float d = abs(dot(m, l)); float h = pow(abs(dot(m, normalize(l + vec3(0.0, 0.0, 1.0)))), 40.0);
        gl_FragColor = vec4(color * (0.42 + 0.62 * d) + 0.12 * h, 1.0); }`;
      function shader(type, src) { const x = gl.createShader(type); gl.shaderSource(x, src); gl.compileShader(x); return x; }
      const prog = gl.createProgram();
      gl.attachShader(prog, shader(gl.VERTEX_SHADER, vs));
      gl.attachShader(prog, shader(gl.FRAGMENT_SHADER, fs));
      gl.linkProgram(prog);
      gl.useProgram(prog);
      const loc = { p: gl.getAttribLocation(prog, "p"), n: gl.getAttribLocation(prog, "n"),
                    R: gl.getUniformLocation(prog, "R"), c: gl.getUniformLocation(prog, "c"),
                    s: gl.getUniformLocation(prog, "s"), dz: gl.getUniformLocation(prog, "dz"), color: gl.getUniformLocation(prog, "color") };
      const bytes = (text) => Uint8Array.from(atob(text), (ch) => ch.charCodeAt(0)).buffer;
      const groups = data.groups.map((g) => {
        const buffer = (target, array) => { const b = gl.createBuffer(); gl.bindBuffer(target, b); gl.bufferData(target, array, gl.STATIC_DRAW); return b; };
        const q = new Uint16Array(bytes(g.points)), p = new Float32Array(q.length);
        for (let i = 0; i < q.length; i++) p[i] = data.lo[i % 3] + data.size[i % 3] * q[i] / 65535;
        const n = Float32Array.from(new Int8Array(bytes(g.normals)), (v) => v / 127);
        const t = g.wide ? new Uint32Array(bytes(g.triangles)) : new Uint16Array(bytes(g.triangles));
        return { p: buffer(gl.ARRAY_BUFFER, p), n: buffer(gl.ARRAY_BUFFER, n), t: buffer(gl.ELEMENT_ARRAY_BUFFER, t),
                 count: t.length, type: g.wide ? gl.UNSIGNED_INT : gl.UNSIGNED_SHORT, color: [1, 3, 5].map((i) => parseInt(g.color.slice(i, i + 2), 16) / 255) };
      });
      const center = [0, 1, 2].map((i) => data.lo[i] + 0.5 * data.size[i]);
      const span = Math.max(...data.size);
      let dragging = null;

      function render() {
        const dpr = window.devicePixelRatio || 1;
        const rect = canvas.getBoundingClientRect();
        canvas.width = Math.max(1, Math.floor(rect.width * dpr));
        canvas.height = Math.max(1, Math.floor(rect.height * dpr));
        gl.viewport(0, 0, canvas.width, canvas.height);
        gl.clearColor(1, 1, 1, 1);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
        gl.enable(gl.DEPTH_TEST);
        const cy = Math.cos(view.yaw), sy = Math.sin(view.yaw), cp = Math.cos(view.pitch), sp = Math.sin(view.pitch);
        // rows: x' = (cy, -sy, 0), y' = cp (sy, cy, 0) - sp (0, 0, 1), z' = sp (sy, cy, 0) + cp (0, 0, 1); column-major
        gl.uniformMatrix3fv(loc.R, false, [cy, cp * sy, sp * sy, -sy, cp * cy, sp * cy, 0, -sp, cp]);
        gl.uniform3fv(loc.c, center);
        gl.uniform1f(loc.dz, 0.9 / span);
        const k = view.zoom * 2.0 / (span * 1.9);
        gl.uniform2fv(loc.s, [k * Math.min(rect.width, rect.height) / rect.width, k * Math.min(rect.width, rect.height) / rect.height]);
        for (const g of groups) {
          gl.bindBuffer(gl.ARRAY_BUFFER, g.p); gl.enableVertexAttribArray(loc.p); gl.vertexAttribPointer(loc.p, 3, gl.FLOAT, false, 0, 0);
          gl.bindBuffer(gl.ARRAY_BUFFER, g.n); gl.enableVertexAttribArray(loc.n); gl.vertexAttribPointer(loc.n, 3, gl.FLOAT, false, 0, 0);
          gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, g.t);
          gl.uniform3fv(loc.color, g.color);
          gl.drawElements(gl.TRIANGLES, g.count, g.type, 0);
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
      window.addEventListener("resize", render);
      render();
    }
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
    data = _view_data([(_surface(parts["sphere below the plane"][0]), ORANGE), (_surface(parts["disk"][0]), BLUE),
                       (_curve(parts["circle"][0], 0.012), MAGENTA)])
    _write_interactive_view(OUTPUT_DIR / "sphere-plane-view.html", "Ball and plane", data)


def napkin_ring_scenes() -> None:
    state = demo_napkin_ring.build()
    demo_napkin_ring.plot(state, OUTPUT_DIR / "napkin-ring-scene.png")
    data = _view_data([(_surface(state["zone"]), ORANGE), (_surface(state["hole"]), BLUE),
                       (_curve(state["rims"], 0.01), MAGENTA)])
    _write_interactive_view(OUTPUT_DIR / "napkin-ring-view.html", "Napkin ring", data)


def gyroid_scenes() -> None:
    state = demo_gyroid.build()
    demo_gyroid.plot(state, OUTPUT_DIR / "gyroid-scene.png")
    lower, upper = state["sheet_faces"]
    data = _view_data([(_surface(lower), BLUE), (_surface(upper), ORANGE)])
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
