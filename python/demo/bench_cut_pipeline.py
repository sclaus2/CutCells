"""
bench_cut_pipeline.py
=====================

Pure-Python micro-benchmark for the CutCells cutting pipeline.

Measures wall time for:
  1. cut_vtk_mesh         – for triangle, quad, tetrahedron, hex meshes
  2. runtime_quadrature   – phi<0 and phi=0 domains
  3. physical_points      – mapping reference → physical space

Run with:
    conda run -n fenicsxv10 python bench_cut_pipeline.py [--n N] [--repeat R] [--order O]

Columns:
    kernel          – what is being timed
    mesh / domain   – short description
    n_cells         – number of background cells
    n_pts           – number of quadrature points (quadrature kernels only)
    total_ms        – total wall time for <repeat> repetitions [ms]
    per_call_us     – time per single call [µs]
    Mcells/s        – million cells per second  (cut kernels)
    Mpts/s          – million quad pts per second (quadrature mapping)
"""

import argparse
import time
from typing import Callable

import numpy as np

try:
    import pyvista as pv

    HAS_PV = True
except ImportError:
    HAS_PV = False

import cutcells

# ---------------------------------------------------------------------------
# Timing helper
# ---------------------------------------------------------------------------


def timeit(fn: Callable, repeat: int = 5) -> float:
    """Return minimum wall time [seconds] over <repeat> calls."""
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)
    return min(times)


# ---------------------------------------------------------------------------
# Mesh creation helpers
# ---------------------------------------------------------------------------


def make_triangle_mesh_vtk(n: int):
    """Return (points_flat, conn, offset, celltypes, ls_vals) for an n×n Delaunay
    triangle mesh on [-1,1]^2, level-set = circle of radius 0.7."""
    if not HAS_PV:
        raise RuntimeError("pyvista required for triangle/quad mesh generation")
    x = np.linspace(-1.0, 1.0, n)
    y = np.linspace(-1.0, 1.0, n)
    xx, yy, zz = np.meshgrid(x, y, [0.0])
    pts = np.c_[xx.ravel(), yy.ravel(), zz.ravel()]
    grid = pv.UnstructuredGrid(pv.PolyData(pts).delaunay_2d())
    pts3 = np.asarray(grid.points, dtype=np.float64)
    conn = np.asarray(grid.cell_connectivity, dtype=np.int32)
    off = np.asarray(grid.offset, dtype=np.int32)
    ctypes = np.asarray(grid.celltypes, dtype=np.int32)
    ls = np.sqrt(pts3[:, 0] ** 2 + pts3[:, 1] ** 2) - 0.7
    return pts3.ravel(), conn, off, ctypes, ls, grid.n_cells


def make_quad_mesh_vtk(n: int):
    """Structured quad mesh on [-1,1]^2 (vectorised)."""
    x = np.linspace(-1.0, 1.0, n)
    y = np.linspace(-1.0, 1.0, n)
    xx, yy = np.meshgrid(x, y)  # both (n, n)
    zz = np.zeros_like(xx)
    pts3 = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=1)  # (n*n, 3)

    # Node index: node_id[i,j] = i*n + j
    i = np.arange(n - 1)
    j = np.arange(n - 1)
    ii, jj = np.meshgrid(i, j, indexing="ij")  # (n-1, n-1)
    ii = ii.ravel()
    jj = jj.ravel()
    n0 = ii * n + jj
    n1 = ii * n + (jj + 1)
    n2 = (ii + 1) * n + (jj + 1)
    n3 = (ii + 1) * n + jj
    conn = np.stack([n0, n1, n2, n3], axis=1).astype(np.int32).ravel()
    n_cells = (n - 1) * (n - 1)
    off = np.arange(0, (n_cells + 1) * 4, 4, dtype=np.int32)
    ct = np.full(n_cells, 9, dtype=np.int32)  # VTK_QUAD = 9
    ls = np.sqrt(pts3[:, 0] ** 2 + pts3[:, 1] ** 2) - 0.7
    return pts3.ravel(), conn, off, ct, ls, n_cells


def make_tet_mesh_vtk(n: int):
    """Random tet mesh built via pyvista Delaunay 3D."""
    if not HAS_PV:
        raise RuntimeError("pyvista required")
    x = np.linspace(-1.0, 1.0, n)
    y = np.linspace(-1.0, 1.0, n)
    z = np.linspace(-1.0, 1.0, n)
    xx, yy, zz = np.meshgrid(x, y, z)
    pts = np.c_[xx.ravel(), yy.ravel(), zz.ravel()]
    grid = pv.UnstructuredGrid(pv.PolyData(pts).delaunay_3d())
    pts3 = np.asarray(grid.points, dtype=np.float64)
    conn = np.asarray(grid.cell_connectivity, dtype=np.int32)
    off = np.asarray(grid.offset, dtype=np.int32)
    ctypes = np.asarray(grid.celltypes, dtype=np.int32)
    ls = np.sqrt(pts3[:, 0] ** 2 + pts3[:, 1] ** 2 + pts3[:, 2] ** 2) - 0.7
    return pts3.ravel(), conn, off, ctypes, ls, grid.n_cells


def make_hex_mesh_vtk(n: int):
    """Structured hex mesh on [-1,1]^3 (vectorised)."""
    x = np.linspace(-1.0, 1.0, n)
    y = np.linspace(-1.0, 1.0, n)
    z = np.linspace(-1.0, 1.0, n)
    # node_id[i,j,k] = i*n*n + j*n + k  (i=z-axis, j=y-axis, k=x-axis)
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")  # each (n,n,n)
    pts3 = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=1)  # (n^3, 3)

    i = np.arange(n - 1)
    j = np.arange(n - 1)
    k = np.arange(n - 1)
    ii, jj, kk = np.meshgrid(i, j, k, indexing="ij")  # each (n-1)^3
    ii = ii.ravel()
    jj = jj.ravel()
    kk = kk.ravel()

    def nid(a, b, c):
        return (a * n + b) * n + c

    conn = (
        np.stack(
            [
                nid(ii, jj, kk),
                nid(ii, jj, kk + 1),
                nid(ii, jj + 1, kk + 1),
                nid(ii, jj + 1, kk),
                nid(ii + 1, jj, kk),
                nid(ii + 1, jj, kk + 1),
                nid(ii + 1, jj + 1, kk + 1),
                nid(ii + 1, jj + 1, kk),
            ],
            axis=1,
        )
        .astype(np.int32)
        .ravel()
    )
    n_cells = (n - 1) ** 3
    off = np.arange(0, (n_cells + 1) * 8, 8, dtype=np.int32)
    ct = np.full(n_cells, 12, dtype=np.int32)  # VTK_HEXAHEDRON = 12
    ls = np.sqrt(pts3[:, 0] ** 2 + pts3[:, 1] ** 2 + pts3[:, 2] ** 2) - 0.7
    return pts3.ravel(), conn, off, ct, ls, n_cells


# ---------------------------------------------------------------------------
# Warm-up  (avoids measuring JIT / module-load costs)
# ---------------------------------------------------------------------------


def warmup():
    pts = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0], dtype=np.float64)
    ls = np.array([-0.1, 0.1, 0.1], dtype=np.float64)
    conn = np.array([0, 1, 2], dtype=np.int32)
    off = np.array([0, 3], dtype=np.int32)
    ct = np.array([5], dtype=np.int32)  # VTK_TRIANGLE
    for _ in range(10):
        cutcells.cut_vtk_mesh(ls, pts, conn, off, ct, "phi<0", True)


# ---------------------------------------------------------------------------
# Benchmark runners
# ---------------------------------------------------------------------------


def bench_cut_vtk_mesh(pts_flat, conn, off, ct, ls, n_cells, label, repeat):
    dt = timeit(
        lambda: cutcells.cut_vtk_mesh(ls, pts_flat, conn, off, ct, "phi<0", True),
        repeat=repeat,
    )
    per_call_us = dt * 1e6
    mcells_s = n_cells / dt / 1e6
    return label, n_cells, None, per_call_us, mcells_s, None


def bench_runtime_quadrature(
    pts_flat, conn, off, ct, ls, n_cells, domain, label, order, repeat
):
    def _run():
        cutcells.runtime_quadrature(ls, pts_flat, conn, off, ct, domain, True, order)

    # collect once to count points
    q = cutcells.runtime_quadrature(ls, pts_flat, conn, off, ct, domain, True, order)
    n_pts = len(q.weights)
    dt = timeit(_run, repeat=repeat)
    per_call_us = dt * 1e6
    mcells_s = n_cells / dt / 1e6
    return label, n_cells, n_pts, per_call_us, mcells_s, None


def bench_physical_points(
    pts_flat, conn, off, ct, ls, n_cells, domain, label, order, repeat
):
    q = cutcells.runtime_quadrature(ls, pts_flat, conn, off, ct, domain, True, order)
    n_pts = len(q.weights)

    def _run():
        cutcells.physical_points(q, pts_flat, conn, off, ct)

    dt = timeit(_run, repeat=repeat)
    per_call_us = dt * 1e6
    mpts_s = n_pts / dt / 1e6 if n_pts > 0 else 0.0
    return label, n_cells, n_pts, per_call_us, None, mpts_s


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

HDR = f"{'kernel':<36} {'n_cells':>8} {'n_pts':>8} {'us/call':>10} {'Mcells/s':>10} {'Mpts/s':>8}"
SEP = "-" * len(HDR)


def print_row(kernel, n_cells, n_pts, per_call_us, mcells_s, mpts_s):
    npts_s = f"{n_pts:>8}" if n_pts is not None else f"{'':>8}"
    mc_s = f"{mcells_s:>10.2f}" if mcells_s is not None else f"{'':>10}"
    mp_s = f"{mpts_s:>8.2f}" if mpts_s is not None else f"{'':>8}"
    print(f"{kernel:<36} {n_cells:>8} {npts_s} {per_call_us:>10.1f} {mc_s} {mp_s}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="CutCells pipeline benchmark")
    parser.add_argument(
        "--n2d",
        type=int,
        default=40,
        help="Grid points per axis for 2-D meshes (default: 40)",
    )
    parser.add_argument(
        "--n3d",
        type=int,
        default=12,
        help="Grid points per axis for 3-D meshes (default: 12)",
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=7,
        help="Repetitions per kernel (min is reported, default: 7)",
    )
    parser.add_argument(
        "--order", type=int, default=3, help="Quadrature order (default: 3)"
    )
    args = parser.parse_args()

    print(
        f"\nCutCells pipeline benchmark  (repeat={args.repeat}, quad order={args.order})"
    )
    print(
        f"  2-D mesh resolution: {args.n2d} pts/axis  (~{2 * (args.n2d - 1) ** 2} tri cells, {(args.n2d - 1) ** 2} quad cells)"
    )
    print(
        f"  3-D mesh resolution: {args.n3d} pts/axis  (~{(args.n3d - 1) ** 3} hex cells)"
    )

    print("\nBuilding meshes …")
    t0 = time.perf_counter()
    tri_data = make_triangle_mesh_vtk(args.n2d)
    quad_data = make_quad_mesh_vtk(args.n2d)
    hex_data = make_hex_mesh_vtk(args.n3d)
    tet_data = make_tet_mesh_vtk(args.n3d) if HAS_PV else None
    print(f"  mesh generation: {(time.perf_counter() - t0) * 1e3:.1f} ms")
    if tet_data is not None:
        print(
            f"  tri  {tri_data[5]:>6} cells | quad {quad_data[5]:>6} cells | "
            f"tet {tet_data[5]:>6} cells | hex {hex_data[5]:>6} cells"
        )
    else:
        print(
            f"  tri  {tri_data[5]:>6} cells | quad {quad_data[5]:>6} cells | "
            f"hex {hex_data[5]:>6} cells"
        )

    print("\nWarming up …")
    warmup()

    print()
    print(HDR)
    print(SEP)

    # ------------------------------------------------------------------
    # 2-D Triangle mesh
    # ------------------------------------------------------------------
    pts, conn, off, ct, ls, nc = tri_data
    for row in [
        bench_cut_vtk_mesh(
            pts, conn, off, ct, ls, nc, "cut_vtk_mesh  tri phi<0", args.repeat
        ),
        bench_runtime_quadrature(
            pts,
            conn,
            off,
            ct,
            ls,
            nc,
            "phi<0",
            "runtime_quad  tri phi<0",
            args.order,
            args.repeat,
        ),
        bench_runtime_quadrature(
            pts,
            conn,
            off,
            ct,
            ls,
            nc,
            "phi=0",
            "runtime_quad  tri phi=0",
            args.order,
            args.repeat,
        ),
        bench_physical_points(
            pts,
            conn,
            off,
            ct,
            ls,
            nc,
            "phi<0",
            "physical_pts  tri phi<0",
            args.order,
            args.repeat,
        ),
        bench_physical_points(
            pts,
            conn,
            off,
            ct,
            ls,
            nc,
            "phi=0",
            "physical_pts  tri phi=0",
            args.order,
            args.repeat,
        ),
    ]:
        print_row(*row)

    print()

    # ------------------------------------------------------------------
    # 2-D Quad mesh
    # ------------------------------------------------------------------
    pts, conn, off, ct, ls, nc = quad_data
    for row in [
        bench_cut_vtk_mesh(
            pts, conn, off, ct, ls, nc, "cut_vtk_mesh  quad phi<0", args.repeat
        ),
        bench_runtime_quadrature(
            pts,
            conn,
            off,
            ct,
            ls,
            nc,
            "phi<0",
            "runtime_quad  quad phi<0",
            args.order,
            args.repeat,
        ),
        bench_runtime_quadrature(
            pts,
            conn,
            off,
            ct,
            ls,
            nc,
            "phi=0",
            "runtime_quad  quad phi=0",
            args.order,
            args.repeat,
        ),
        bench_physical_points(
            pts,
            conn,
            off,
            ct,
            ls,
            nc,
            "phi<0",
            "physical_pts  quad phi<0",
            args.order,
            args.repeat,
        ),
    ]:
        print_row(*row)

    print()

    # ------------------------------------------------------------------
    # 3-D Tetrahedron mesh
    # ------------------------------------------------------------------
    if HAS_PV and tet_data is not None:
        pts, conn, off, ct, ls, nc = tet_data
        for row in [
            bench_cut_vtk_mesh(
                pts, conn, off, ct, ls, nc, "cut_vtk_mesh  tet phi<0", args.repeat
            ),
            bench_runtime_quadrature(
                pts,
                conn,
                off,
                ct,
                ls,
                nc,
                "phi<0",
                "runtime_quad  tet phi<0",
                args.order,
                args.repeat,
            ),
            bench_runtime_quadrature(
                pts,
                conn,
                off,
                ct,
                ls,
                nc,
                "phi=0",
                "runtime_quad  tet phi=0",
                args.order,
                args.repeat,
            ),
            bench_physical_points(
                pts,
                conn,
                off,
                ct,
                ls,
                nc,
                "phi<0",
                "physical_pts  tet phi<0",
                args.order,
                args.repeat,
            ),
        ]:
            print_row(*row)
        print()

    # ------------------------------------------------------------------
    # 3-D Hex mesh
    # ------------------------------------------------------------------
    pts, conn, off, ct, ls, nc = hex_data
    for row in [
        bench_cut_vtk_mesh(
            pts, conn, off, ct, ls, nc, "cut_vtk_mesh  hex phi<0", args.repeat
        ),
        bench_runtime_quadrature(
            pts,
            conn,
            off,
            ct,
            ls,
            nc,
            "phi<0",
            "runtime_quad  hex phi<0",
            args.order,
            args.repeat,
        ),
        bench_runtime_quadrature(
            pts,
            conn,
            off,
            ct,
            ls,
            nc,
            "phi=0",
            "runtime_quad  hex phi=0",
            args.order,
            args.repeat,
        ),
        bench_physical_points(
            pts,
            conn,
            off,
            ct,
            ls,
            nc,
            "phi<0",
            "physical_pts  hex phi<0",
            args.order,
            args.repeat,
        ),
    ]:
        print_row(*row)

    print(SEP)
    print("Times are minimum over all repetitions (best-of-N).")
    print("Mcells/s = background cells processed per second.")
    print("Mpts/s   = quadrature points mapped per second.")


if __name__ == "__main__":
    main()
