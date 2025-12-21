"""Demo: cut a complex level set through a *quadrilateral* (VTK_QUAD) mesh.

Run (from repo root):
  conda run -n fenicsx-dev python python/demo/cut_vtk_mesh/cut_quad_mesh_2D.py

Optional:
  conda run -n fenicsx-dev python python/demo/cut_vtk_mesh/cut_quad_mesh_2D.py --n 60 --save cut_quad_mesh.vtu

Notes:
- Uses `cutcells.cut_vtk_mesh` (C++ implementation + generated quad tables).
- Requires `pyvista` for plotting/saving VTU in Python.
"""

from __future__ import annotations

import argparse

import numpy as np

import cutcells


def level_set_popcorn_2d(xyz: np.ndarray) -> float:
    """2D 'popcorn'-like implicit function on z=0 plane.

    Negative values are "inside".
    """
    x, y = float(xyz[0]), float(xyz[1])

    r0 = 0.65
    r = (x * x + y * y) ** 0.5

    value = r - r0

    # Add bumps around a ring (gaussian perturbations)
    sigma2 = 0.02
    for k in range(10):
        ang = 2.0 * np.pi * k / 10.0
        cx = 0.9 * r0 * np.cos(ang)
        cy = 0.9 * r0 * np.sin(ang)
        d2 = (x - cx) ** 2 + (y - cy) ** 2
        value -= 0.22 * np.exp(-d2 / sigma2)

    # Add mild sinusoidal ripples to make the interface less symmetric
    value += 0.05 * np.sin(6.0 * x) * np.cos(5.0 * y)

    return value


def create_quad_mesh(x0: float, y0: float, x1: float, y1: float, nx: int, ny: int):
    """Create a structured quad mesh in VTK-unstructured-array form.

    Returns:
      points: (npts, 3)
      connectivity: (ncells*4,)
      offset: (ncells,)
      celltypes: (ncells,) filled with VTK_QUAD (9)
    """
    xs = np.linspace(x0, x1, nx)
    ys = np.linspace(y0, y1, ny)

    points = np.zeros((nx * ny, 3), dtype=float)
    p = 0
    for j in range(ny):
        for i in range(nx):
            points[p, 0] = xs[i]
            points[p, 1] = ys[j]
            points[p, 2] = 0.0
            p += 1

    def vid(i: int, j: int) -> int:
        return j * nx + i

    ncx = nx - 1
    ncy = ny - 1
    ncells = ncx * ncy

    connectivity = np.zeros((ncells * 4,), dtype=np.int32)
    offset = np.zeros((ncells,), dtype=np.int32)
    celltypes = np.full((ncells,), 9, dtype=np.int32)  # VTK_QUAD

    c = 0
    k = 0
    for j in range(ncy):
        for i in range(ncx):
            offset[c] = k
            # (i,j)->(i+1,j)->(i+1,j+1)->(i,j+1)
            connectivity[k + 0] = vid(i, j)
            connectivity[k + 1] = vid(i + 1, j)
            connectivity[k + 2] = vid(i + 1, j + 1)
            connectivity[k + 3] = vid(i, j + 1)
            k += 4
            c += 1

    return points, connectivity, offset, celltypes


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cut a complex level set through a quadrilateral mesh"
    )
    parser.add_argument(
        "--n", type=int, default=45, help="Points per axis (mesh is (n x n) points)"
    )
    parser.add_argument(
        "--save", type=str, default=None, help="Optional output VTU filename"
    )
    parser.add_argument(
        "--triangulate",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Triangulate cut cells (default: True). Use --no-triangulate to preserve quads.",
    )
    args = parser.parse_args()

    try:
        import pyvista as pv
    except Exception as exc:  # pragma: no cover
        raise SystemExit(
            "pyvista is required for this demo. Install with: `python -m pip install pyvista`\n"
            f"Import error: {exc}"
        )

    n = int(args.n)
    points, conn, off, celltypes = create_quad_mesh(-1.0, -1.0, 1.0, 1.0, n, n)

    # level set at points
    ls = np.array([level_set_popcorn_2d(p) for p in points], dtype=float)

    # Cut mesh and also extract all inside (phi<0) background cells for context
    cut_mesh = cutcells.cut_vtk_mesh(
        ls, points, conn, off, celltypes, "phi<0", triangulate=bool(args.triangulate)
    )
    inside_cells = cutcells.locate_cells(ls, points, conn, off, celltypes, "phi<0")

    # Build a pyvista grid for visualization
    cells_with_counts = np.empty((len(off) * 5,), dtype=np.int32)
    cells_with_counts[0::5] = 4
    cells_with_counts[1::5] = conn[0::4]
    cells_with_counts[2::5] = conn[1::4]
    cells_with_counts[3::5] = conn[2::4]
    cells_with_counts[4::5] = conn[3::4]

    bg = pv.UnstructuredGrid(cells_with_counts, celltypes, points)

    cut_points = np.asarray(cut_mesh.vertex_coords, dtype=float).reshape((-1, 3))
    cut_cells = np.asarray(cut_mesh.cells, dtype=np.int32)
    cut_types = np.asarray(cut_mesh.types, dtype=np.int32)
    pv_cut = pv.UnstructuredGrid(cut_cells, cut_types, cut_points)
    extract = bg.extract_cells(inside_cells)
    merged = extract.merge(pv_cut)

    pl = pv.Plotter()
    pl.add_mesh(bg, color="white", show_edges=True, opacity=0.25)
    pl.add_mesh(merged, show_edges=True)
    pl.camera_position = "xy"

    if args.save:
        merged.save(args.save)
        print(f"Saved: {args.save}")

    pl.show()


if __name__ == "__main__":
    main()
