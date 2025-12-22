import argparse

import cutcells
import numpy as np
import pyvista as pv


def unit_square_pyramid_vertices() -> np.ndarray:
    # VTK_PYRAMID / CutCells pyramid vertex order:
    # base quad: 0:(0,0,0) 1:(1,0,0) 2:(1,1,0) 3:(0,1,0)
    # apex:      4:(0.5,0.5,1)
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ],
        dtype=float,
    )


def case_id_from_ls_values(ls_values: np.ndarray) -> int:
    # CutCells convention: a vertex contributes bit i if phi_i < 0
    case_id = 0
    for i, v in enumerate(ls_values.tolist()):
        if v < 0:
            case_id |= 1 << i
    return case_id


def is_strictly_inside_unit_pyramid(p: np.ndarray, tol: float = 1e-12) -> bool:
    # Strict interior of the unit square pyramid with apex at (0.5,0.5,1).
    # For a horizontal slice at height z, the cross-section is a square centered at (0.5,0.5)
    # with half-width 0.5*(1-z).
    x, y, z = float(p[0]), float(p[1]), float(p[2])
    if not (z > tol and z < 1.0 - tol):
        return False
    half = 0.5 * (1.0 - z)
    return (abs(x - 0.5) < (half - tol)) and (abs(y - 0.5) < (half - tol))


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Demonstrate a pyramid (VTK_PYRAMID) volume cut that introduces an interior special point (VTK's N0)."
        )
    )
    parser.add_argument(
        "--screenshot",
        default="pyramid_n0_demo.png",
        help="Output screenshot filename (set to empty string to disable).",
    )
    parser.add_argument(
        "--show",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Show an interactive window (disable for headless runs).",
    )
    args = parser.parse_args()

    vertex_coordinates = unit_square_pyramid_vertices()

    # Pick a concrete case that (per generated pyramid volume tables) uses token 200 (N0).
    # Pyramid has 5 vertices => 32 cases; case id 19 corresponds to mask bits set at vertices {0,1,4}.
    ls_values = np.array([-1.0, -1.0, 1.0, 1.0, -1.0], dtype=float)
    case_id = case_id_from_ls_values(ls_values)

    cut_cell = cutcells.cut(
        cutcells.CellType.pyramid,
        vertex_coordinates,
        3,
        ls_values,
        "phi<0",
        False,
    )

    points = np.asarray(cut_cell.vertex_coords)
    grid = pv.UnstructuredGrid(cut_cell.connectivity, cut_cell.types, points)

    interior_mask = np.array(
        [is_strictly_inside_unit_pyramid(p) for p in points], dtype=bool
    )

    print(f"pyramid case_id (phi<0 bits): {case_id}")
    print(f"cut_cell: n_vertices={points.shape[0]}, n_cells={grid.n_cells}")
    print(f"strictly interior points: {int(interior_mask.sum())}")
    if interior_mask.any():
        print("example interior point(s) (typically includes N0):")
        for p in points[interior_mask][:8]:
            print("  ", p)

    pl = pv.Plotter(off_screen=not args.show)

    # Reference pyramid wireframe
    faces = np.hstack(
        [
            [4, 0, 1, 2, 3],
            [3, 0, 1, 4],
            [3, 1, 2, 4],
            [3, 2, 3, 4],
            [3, 3, 0, 4],
        ]
    ).astype(np.int64)
    surf = pv.PolyData(vertex_coordinates, faces)
    pl.add_mesh(surf, style="wireframe", color="black", line_width=1)

    pl.add_mesh(grid, show_edges=True, opacity=0.75)

    pc = pv.PolyData(points)
    pc["interior"] = interior_mask.astype(np.int32)
    pl.add_mesh(
        pc,
        scalars="interior",
        render_points_as_spheres=True,
        point_size=14,
        cmap=["dodgerblue", "crimson"],
    )

    pl.camera_position = "iso"

    screenshot = args.screenshot if args.screenshot else None
    pl.show(screenshot=screenshot)


if __name__ == "__main__":
    main()
