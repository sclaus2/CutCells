import argparse

import cutcells
import numpy as np
import pyvista as pv


def level_set(x):
    r0 = 0.6
    r = np.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])

    s = r0 - r
    for k in range(0, 5):
        xk = r0 / np.sqrt(5.0) * 2 * np.cos(2 * k * np.pi / 5)
        yk = r0 / np.sqrt(5.0) * 2 * np.sin(2 * k * np.pi / 5)
        zk = r0 / np.sqrt(5.0)
        tmp = (
            (x[0] - xk) * (x[0] - xk)
            + (x[1] - yk) * (x[1] - yk)
            + (x[2] - zk) * (x[2] - zk)
        )
        s += 2 * np.exp(-tmp / (0.04))

    for k in range(5, 10):
        xk = r0 / np.sqrt(5.0) * 2 * np.cos((2 * (k - 5) - 1) * np.pi / 5)
        yk = r0 / np.sqrt(5.0) * 2 * np.sin((2 * (k - 5) - 1) * np.pi / 5)
        zk = -r0 / np.sqrt(5.0)
        tmp = (
            (x[0] - xk) * (x[0] - xk)
            + (x[1] - yk) * (x[1] - yk)
            + (x[2] - zk) * (x[2] - zk)
        )
        s += 2 * np.exp(-tmp / (0.04))

    for k in range(10, 11):
        xk = 0
        yk = 0
        zk = r0
        tmp = (
            (x[0] - xk) * (x[0] - xk)
            + (x[1] - yk) * (x[1] - yk)
            + (x[2] - zk) * (x[2] - zk)
        )
        s += 2 * np.exp(-tmp / (0.04))

    for k in range(11, 12):
        xk = 0
        yk = 0
        zk = -r0
        tmp = (
            (x[0] - xk) * (x[0] - xk)
            + (x[1] - yk) * (x[1] - yk)
            + (x[2] - zk) * (x[2] - zk)
        )
        s += 2 * np.exp(-tmp / (0.04))

    return -1.0 * s


def create_hex_box_mesh(x0, y0, z0, x1, y1, z1, nx, ny, nz) -> pv.UnstructuredGrid:
    xs = np.linspace(x0, x1, num=nx)
    ys = np.linspace(y0, y1, num=ny)
    zs = np.linspace(z0, z1, num=nz)

    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing="ij")
    sg = pv.StructuredGrid(xx, yy, zz)
    return sg.cast_to_unstructured_grid()


def create_cut_mesh(
    grid: pv.UnstructuredGrid, triangulate: bool
) -> pv.UnstructuredGrid:
    points = grid.points

    ls_values = np.zeros(len(points), dtype=float)
    for j, point in enumerate(points):
        ls_values[j] = level_set(point)

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

    pv_cut = pv.UnstructuredGrid(
        cut_mesh.cells,
        cut_mesh.types,
        cut_mesh.vertex_coords,
    )
    extract = grid.extract_cells(inside_cells)

    return extract.merge(pv_cut)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cut a popcorn level set from a hexahedral (VTK_HEXAHEDRON) mesh."
    )
    parser.add_argument(
        "--n", type=int, default=21, help="Number of points per axis (>= 2)."
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
        default="popcorn_hex_mesh3D.png",
        help="Output screenshot filename (set to empty string to disable).",
    )
    parser.add_argument(
        "--show",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Show an interactive window (disable for headless runs).",
    )
    args = parser.parse_args()

    grid = create_hex_box_mesh(-1, -1, -1, 1, 1, 1, args.n, args.n, args.n)
    mesh = create_cut_mesh(grid, triangulate=args.triangulate)

    # Quick diagnostics (helps distinguish "too many cells" from bad geometry).
    print(
        f"input: n_points={grid.n_points}, n_cells={grid.n_cells}, unique_celltypes={sorted(set(map(int, grid.celltypes)))}"
    )
    print(f"output: n_points={mesh.n_points}, n_cells={mesh.n_cells}")

    pl = pv.Plotter(off_screen=not args.show)
    if args.show_grid:
        pl.add_mesh(grid, show_edges=True, style="wireframe")
    pl.add_mesh(mesh, show_edges=args.show_edges)
    pl.camera_position = "xz"

    if args.screenshot:
        pl.show(screenshot=args.screenshot)
    else:
        pl.show()


if __name__ == "__main__":
    main()
