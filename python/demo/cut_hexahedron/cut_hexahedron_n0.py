import argparse

import cutcells
import numpy as np
import pyvista as pv


def unit_cube_vertices() -> np.ndarray:
    # VTK_HEXAHEDRON / CutCells hexahedron vertex order:
    # 0:(0,0,0) 1:(1,0,0) 2:(1,1,0) 3:(0,1,0)
    # 4:(0,0,1) 5:(1,0,1) 6:(1,1,1) 7:(0,1,1)
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
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


def on_unit_cube_edge(p: np.ndarray, tol: float = 1e-12) -> bool:
    # A point is on an edge if at least two coordinates are on the boundary {0,1}
    on_boundary = np.isclose(p, 0.0, atol=tol) | np.isclose(p, 1.0, atol=tol)
    return int(on_boundary.sum()) >= 2


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Demonstrate a hexahedron interface cut that introduces an interior special point (VTK's N0)."
    )
    parser.add_argument(
        "--screenshot",
        default="n0_demo.png",
        help="Output screenshot filename (set to empty string to disable).",
    )
    parser.add_argument(
        "--show",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Show an interactive window (disable for headless runs).",
    )
    args = parser.parse_args()

    vertex_coordinates = unit_cube_vertices()

    # This sign pattern corresponds to a VTK hexahedron case known to use N0 in the table.
    # (Used in the regression test: vertices 0,1,6 are inside (phi<0), others outside.)
    ls_values = np.array([-1.0, -1.1, 1.2, 1.3, 1.4, 1.5, -1.6, 1.7], dtype=float)
    case_id = case_id_from_ls_values(ls_values)

    cut_cell = cutcells.cut(
        cutcells.CellType.hexahedron,
        vertex_coordinates,
        3,
        ls_values,
        "phi=0",
        False,
    )

    points = np.asarray(cut_cell.vertex_coords, dtype=float).reshape(-1, 3)
    grid = pv.UnstructuredGrid(cut_cell.connectivity, cut_cell.types, points)

    edge_mask = np.array([on_unit_cube_edge(p) for p in points], dtype=bool)
    interior_points = points[~edge_mask]

    print(f"hexahedron case_id (phi<0 bits): {case_id}")
    print(f"cut_cell: n_vertices={points.shape[0]}, n_cells={grid.n_cells}")
    print(f"points on unit-cube edges: {int(edge_mask.sum())}")
    print(f"points not on edges (interior/face): {int((~edge_mask).sum())}")
    if interior_points.size:
        print("example non-edge point(s) (typically includes N0):")
        for p in interior_points[:8]:
            print("  ", p)

    # Visualization
    pl = pv.Plotter(off_screen=not args.show)

    # Show the reference cube as wireframe
    cube = pv.Cube(center=(0.5, 0.5, 0.5), x_length=1.0, y_length=1.0, z_length=1.0)
    pl.add_mesh(cube, style="wireframe", color="black", line_width=1)

    # Show the interface mesh
    pl.add_mesh(grid, show_edges=True, opacity=0.75)

    # Highlight interior/non-edge points (where N0 typically appears)
    point_cloud = pv.PolyData(points)
    point_cloud["non_edge"] = (~edge_mask).astype(np.int32)
    pl.add_mesh(
        point_cloud,
        scalars="non_edge",
        render_points_as_spheres=True,
        point_size=14,
        cmap=["dodgerblue", "crimson"],
    )

    pl.camera_position = "iso"

    screenshot = args.screenshot if args.screenshot else None
    pl.show(screenshot=screenshot)


if __name__ == "__main__":
    main()
