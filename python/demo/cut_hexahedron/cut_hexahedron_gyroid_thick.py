import argparse

import cutcells
import numpy as np
import pyvista as pv


def gyroid(x: np.ndarray, k: float) -> float:
    # Classic gyroid implicit: sin(kx)cos(ky) + sin(ky)cos(kz) + sin(kz)cos(kx)
    xx, yy, zz = float(x[0]), float(x[1]), float(x[2])
    return (
        np.sin(k * xx) * np.cos(k * yy)
        + np.sin(k * yy) * np.cos(k * zz)
        + np.sin(k * zz) * np.cos(k * xx)
    )


def level_set_thick_gyroid(x: np.ndarray, k: float, thickness: float) -> float:
    # A "thick" gyroid shell: |g(x)| < thickness  =>  phi(x) = |g(x)| - thickness < 0
    return abs(gyroid(x, k)) - thickness


def create_structured_hex_mesh(
    x0, y0, z0, x1, y1, z1, nx, ny, nz
) -> pv.UnstructuredGrid:
    # Use VTK structured grid -> cast to unstructured hex mesh.
    xs = np.linspace(x0, x1, num=nx)
    ys = np.linspace(y0, y1, num=ny)
    zs = np.linspace(z0, z1, num=nz)

    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing="ij")
    sg = pv.StructuredGrid(xx, yy, zz)
    return sg.cast_to_unstructured_grid()


def create_cut_mesh(
    grid: pv.UnstructuredGrid, k: float, thickness: float, triangulate: bool
) -> pv.UnstructuredGrid:
    points = grid.points

    ls_values = np.zeros(len(points), dtype=float)
    for j, point in enumerate(points):
        ls_values[j] = level_set_thick_gyroid(point, k=k, thickness=thickness)

    cut_mesh = cutcells.cut_vtk_mesh(
        ls_values,
        points,
        grid.cell_connectivity,
        grid.offset,
        grid.celltypes,
        "phi<0",
        triangulate,
    )

    inside_cells = cutcells.locate_cells(
        ls_values,
        points,
        grid.cell_connectivity,
        grid.offset,
        grid.celltypes,
        "phi<0",
    )

    cut_points = np.asarray(cut_mesh.vertex_coords, dtype=float).reshape(-1, 3)
    pv_cut = pv.UnstructuredGrid(cut_mesh.cells, cut_mesh.types, cut_points)
    extract = grid.extract_cells(inside_cells)

    return extract.merge(pv_cut)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cut a thick gyroid (|g(x)| < t) from a hexahedral (VTK_HEXAHEDRON) mesh."
    )
    parser.add_argument(
        "--n", type=int, default=41, help="Number of points per axis (>= 2)."
    )
    parser.add_argument(
        "--periods",
        type=float,
        default=2.0,
        help="Number of gyroid periods across the domain.",
    )
    parser.add_argument(
        "--thickness",
        type=float,
        default=0.30,
        help="Shell half-thickness in gyroid level-set units (phi = |g|-thickness).",
    )
    parser.add_argument(
        "--triangulate",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Triangulate cut output (more robust for visualization).",
    )
    parser.add_argument(
        "--show-grid",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Overlay the input hex grid (wireframe).",
    )
    parser.add_argument(
        "--show-edges",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Show mesh edges (can look very busy for larger N).",
    )
    parser.add_argument(
        "--screenshot",
        default="gyroid_hex_mesh3D.png",
        help="Output screenshot filename (set to empty string to disable).",
    )
    parser.add_argument(
        "--show",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Show an interactive window (disable for headless runs).",
    )
    parser.add_argument(
        "--write-prefix",
        default="",
        help="If set, write interface/inside/outside VTK files using this prefix (e.g. /tmp/gyroid).",
    )
    args = parser.parse_args()

    if args.n < 2:
        raise ValueError("--n must be >= 2")

    # Domain [-1,1]^3. Choose k so that `periods` periods fit into length 2:
    # period = 2*pi/k, so periods = 2 / (2*pi/k) = k/pi  => k = periods*pi
    k = float(args.periods) * np.pi

    grid = create_structured_hex_mesh(-1, -1, -1, 1, 1, 1, args.n, args.n, args.n)
    mesh = create_cut_mesh(
        grid, k=k, thickness=float(args.thickness), triangulate=args.triangulate
    )

    print(
        f"input: n_points={grid.n_points}, n_cells={grid.n_cells}, unique_celltypes={sorted(set(map(int, grid.celltypes)))}"
    )
    print(f"output: n_points={mesh.n_points}, n_cells={mesh.n_cells}")

    pl = pv.Plotter(off_screen=not args.show)
    if args.show_grid:
        pl.add_mesh(grid, show_edges=True, style="wireframe")
    pl.add_mesh(mesh, show_edges=args.show_edges)
    pl.camera_position = "iso"

    # Optionally write interface, inside and outside VTK files
    if args.write_prefix:
        pts = grid.points
        ls_values = np.zeros(len(pts), dtype=float)
        for j, point in enumerate(pts):
            ls_values[j] = level_set_thick_gyroid(
                point, k=k, thickness=float(args.thickness)
            )

        conn = grid.cell_connectivity
        off = grid.offset
        ctypes = grid.celltypes

        # Interface (phi=0)
        try:
            iface = cutcells.cut_vtk_mesh(
                ls_values, pts, conn, off, ctypes, "phi=0", args.triangulate
            )
            iface_pts = np.asarray(iface.vertex_coords, dtype=float).reshape(-1, 3)
            iface_pv = pv.UnstructuredGrid(iface.cells, iface.types, iface_pts)
            iface_file = f"{args.write_prefix}_interface.vtu"
            iface_pv.save(iface_file)
            print(f"Wrote interface VTK -> {iface_file}")
        except Exception as e:  # pragma: no cover
            print("Failed to write interface mesh:", e)

        # Inside (phi<0)
        try:
            inside = cutcells.cut_vtk_mesh(
                ls_values, pts, conn, off, ctypes, "phi<0", args.triangulate
            )
            inside_pts = np.asarray(inside.vertex_coords, dtype=float).reshape(-1, 3)
            inside_pv = pv.UnstructuredGrid(inside.cells, inside.types, inside_pts)
            inside_file = f"{args.write_prefix}_inside.vtu"
            inside_pv.save(inside_file)
            print(f"Wrote inside VTK -> {inside_file}")
        except Exception as e:  # pragma: no cover
            print("Failed to write inside mesh:", e)

        # Outside (phi>0)
        try:
            outside = cutcells.cut_vtk_mesh(
                ls_values, pts, conn, off, ctypes, "phi>0", args.triangulate
            )
            outside_pts = np.asarray(outside.vertex_coords, dtype=float).reshape(-1, 3)
            outside_pv = pv.UnstructuredGrid(outside.cells, outside.types, outside_pts)
            outside_file = f"{args.write_prefix}_outside.vtu"
            outside_pv.save(outside_file)
            print(f"Wrote outside VTK -> {outside_file}")
        except Exception as e:  # pragma: no cover
            print("Failed to write outside mesh:", e)

    if args.screenshot:
        pl.show(screenshot=args.screenshot)
    else:
        pl.show()


if __name__ == "__main__":
    main()
