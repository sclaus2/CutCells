import argparse

import cutcells
import numpy as np
import pyvista as pv


VTK_HEXAHEDRON = 12


def level_set(x: np.ndarray) -> float:
    # A smooth but nontrivial implicit surface: sphere with bumps.
    r0 = 0.6
    r = np.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])

    s = r0 - r

    for k in range(5):
        xk = r0 / np.sqrt(5.0) * 2.0 * np.cos(2.0 * k * np.pi / 5.0)
        yk = r0 / np.sqrt(5.0) * 2.0 * np.sin(2.0 * k * np.pi / 5.0)
        zk = r0 / np.sqrt(5.0)
        d2 = (x[0] - xk) ** 2 + (x[1] - yk) ** 2 + (x[2] - zk) ** 2
        s += 2.0 * np.exp(-d2 / 0.04)

    for k in range(5, 10):
        kk = k - 5
        xk = r0 / np.sqrt(5.0) * 2.0 * np.cos((2.0 * kk - 1.0) * np.pi / 5.0)
        yk = r0 / np.sqrt(5.0) * 2.0 * np.sin((2.0 * kk - 1.0) * np.pi / 5.0)
        zk = -r0 / np.sqrt(5.0)
        d2 = (x[0] - xk) ** 2 + (x[1] - yk) ** 2 + (x[2] - zk) ** 2
        s += 2.0 * np.exp(-d2 / 0.04)

    # Flip sign so phi<0 picks the "inside" of the bumpy shape.
    return -s


def _vertex_id(i: int, j: int, k: int, nx: int, ny: int) -> int:
    return i + (nx + 1) * (j + (ny + 1) * k)


def create_structured_hex_grid(x0, y0, z0, x1, y1, z1, nx, ny, nz):
    xs = np.linspace(x0, x1, nx + 1)
    ys = np.linspace(y0, y1, ny + 1)
    zs = np.linspace(z0, z1, nz + 1)

    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing="ij")
    points = np.c_[xx.ravel(), yy.ravel(), zz.ravel()].astype(float)

    num_cells = nx * ny * nz
    connectivity = np.empty(num_cells * 8, dtype=np.int32)
    offset = np.empty(num_cells, dtype=np.int32)
    celltypes = np.full(num_cells, VTK_HEXAHEDRON, dtype=np.int32)

    c = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                v0 = _vertex_id(i, j, k, nx, ny)
                v1 = _vertex_id(i + 1, j, k, nx, ny)
                v2 = _vertex_id(i + 1, j + 1, k, nx, ny)
                v3 = _vertex_id(i, j + 1, k, nx, ny)
                v4 = _vertex_id(i, j, k + 1, nx, ny)
                v5 = _vertex_id(i + 1, j, k + 1, nx, ny)
                v6 = _vertex_id(i + 1, j + 1, k + 1, nx, ny)
                v7 = _vertex_id(i, j + 1, k + 1, nx, ny)

                offset[c] = 8 * c
                connectivity[8 * c : 8 * c + 8] = [v0, v1, v2, v3, v4, v5, v6, v7]
                c += 1

    return points, connectivity, offset, celltypes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=14, help="Cells per axis")
    parser.add_argument(
        "--triangulate",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Triangulate interface quads (default: True)",
    )
    args = parser.parse_args()

    n = args.n
    points, connectivity, offset, celltypes = create_structured_hex_grid(
        -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, n, n, n
    )

    ls_values = np.array([level_set(p) for p in points], dtype=float)

    cut_mesh = cutcells.cut_vtk_mesh(
        ls_values,
        points,
        connectivity,
        offset,
        celltypes,
        "phi<0",
        triangulate=args.triangulate,
    )

    # PyVista expects VTK "cells" layout: [n0, ids..., n1, ids..., ...]
    cells = np.empty(celltypes.size * 9, dtype=np.int32)
    cells[0::9] = 8
    cells.reshape(-1, 9)[:, 1:] = connectivity.reshape(-1, 8)

    pv_in = pv.UnstructuredGrid(cells, celltypes, points)
    pv_cut = pv.UnstructuredGrid(
        cut_mesh.cells,
        cut_mesh.types,
        cut_mesh.vertex_coords,
    )

    pl = pv.Plotter()
    pl.add_mesh(pv_in, show_edges=True, style="wireframe", opacity=0.15)
    pl.add_mesh(pv_cut, show_edges=True)
    pl.camera_position = "xy"
    pl.show()


if __name__ == "__main__":
    main()
