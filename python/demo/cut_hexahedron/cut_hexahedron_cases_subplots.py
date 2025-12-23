import argparse
import math

import cutcells
import numpy as np
import pyvista as pv


# VTK cell type ids (used by PyVista UnstructuredGrid)
VTK_TETRA = 10
VTK_HEXAHEDRON = 12
VTK_WEDGE = 13
VTK_PYRAMID = 14


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


def ls_values_from_case_id(case_id: int) -> np.ndarray:
    # CutCells convention: bit i is set when phi_i < 0.
    # We use a simple +/-1 pattern.
    ls = np.empty(8, dtype=float)
    for i in range(8):
        ls[i] = -1.0 if ((case_id >> i) & 1) else 1.0
    return ls


def on_unit_cube_edge(p: np.ndarray, tol: float = 1e-12) -> bool:
    # A point is on an edge if at least two coordinates are on the boundary {0,1}
    on_boundary = np.isclose(p, 0.0, atol=tol) | np.isclose(p, 1.0, atol=tol)
    return int(on_boundary.sum()) >= 2


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Visualize VTK_HEXAHEDRON interface-cut cases in 3x5 interactive subplots (15 at a time). "
            "Use 'n'/'p' keys to page through batches."
        )
    )
    parser.add_argument(
        "--batch",
        type=int,
        default=0,
        help="0-based batch index (each batch shows 15 intersected cases).",
    )
    parser.add_argument(
        "--screenshot",
        default="",
        help="Optional screenshot filename (empty disables).",
    )
    parser.add_argument(
        "--show",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Show an interactive window (disable for headless runs).",
    )
    args = parser.parse_args()

    vertex_coordinates = unit_cube_vertices()

    # Exclude non-intersected masks.
    case_ids = [i for i in range(256) if i not in (0, 255)]
    batch_size = 15
    num_batches = int(math.ceil(len(case_ids) / batch_size))

    batch = int(args.batch)
    if batch < 0:
        batch = 0
    if batch >= num_batches:
        batch = num_batches - 1

    pl = pv.Plotter(shape=(3, 5), off_screen=not args.show)

    cube = pv.Cube(center=(0.5, 0.5, 0.5), x_length=1.0, y_length=1.0, z_length=1.0)

    def render_batch(batch_index: int) -> None:
        pl.clear()

        start = batch_index * batch_size
        for i in range(batch_size):
            global_idx = start + i
            row = i // 5
            col = i % 5
            pl.subplot(row, col)

            if global_idx >= len(case_ids):
                pl.add_text("(no case)", font_size=10)
                pl.add_mesh(cube, style="wireframe", color="black", line_width=1)
                pl.camera_position = "iso"
                continue

            case_id = case_ids[global_idx]
            ls_values = ls_values_from_case_id(case_id)

            try:
                # Inside volume (phi<0)
                inside_cell = cutcells.cut(
                    cutcells.CellType.hexahedron,
                    vertex_coordinates,
                    3,
                    ls_values,
                    "phi<0",
                    False,
                )
                inside_points = np.asarray(inside_cell.vertex_coords)
                inside_grid = pv.UnstructuredGrid(
                    inside_cell.connectivity, inside_cell.types, inside_points
                )

                # Interface surface (phi=0)
                cut_cell = cutcells.cut(
                    cutcells.CellType.hexahedron,
                    vertex_coordinates,
                    3,
                    ls_values,
                    "phi=0",
                    False,
                )

                points = np.asarray(cut_cell.vertex_coords)
                grid = pv.UnstructuredGrid(
                    cut_cell.connectivity, cut_cell.types, points
                )

                edge_mask = np.array([on_unit_cube_edge(p) for p in points], dtype=bool)

                pl.add_text(
                    f"case {case_id}",
                    font_size=10,
                )
                pl.add_mesh(cube, style="wireframe", color="black", line_width=1)

                # Render inside composition behind the interface
                inside_types = np.asarray(inside_grid.celltypes, dtype=int)
                type_specs = [
                    (VTK_HEXAHEDRON, "hex", "seagreen"),
                    (VTK_WEDGE, "prism", "gold"),
                    (VTK_PYRAMID, "pyr", "orchid"),
                    (VTK_TETRA, "tet", "lightskyblue"),
                ]

                counts = []
                for vtk_id, label, color in type_specs:
                    cell_ids = np.where(inside_types == vtk_id)[0]
                    if cell_ids.size == 0:
                        continue
                    subgrid = inside_grid.extract_cells(cell_ids)
                    pl.add_mesh(
                        subgrid,
                        color=color,
                        opacity=0.25,
                        show_edges=True,
                        line_width=1,
                    )
                    counts.append(f"{label}:{int(cell_ids.size)}")

                if counts:
                    pl.add_text("\n".join(counts), position="lower_left", font_size=9)

                pl.add_mesh(grid, show_edges=True, opacity=0.75)

                pc = pv.PolyData(points)
                pc["non_edge"] = (~edge_mask).astype(np.int32)
                pl.add_mesh(
                    pc,
                    scalars="non_edge",
                    render_points_as_spheres=True,
                    point_size=10,
                    cmap=["dodgerblue", "crimson"],
                )
            except Exception as e:  # pragma: no cover
                pl.add_text(f"case {case_id}\nERROR\n{type(e).__name__}", font_size=10)
                pl.add_mesh(cube, style="wireframe", color="black", line_width=1)

            pl.camera_position = "iso"

        pl.render()
        print(
            f"batch {batch_index + 1}/{num_batches} (cases {start}..{min(start + batch_size - 1, len(case_ids) - 1)})"
        )

    def next_batch() -> None:
        nonlocal batch
        batch = (batch + 1) % num_batches
        render_batch(batch)

    def prev_batch() -> None:
        nonlocal batch
        batch = (batch - 1) % num_batches
        render_batch(batch)

    pl.add_key_event("n", next_batch)
    pl.add_key_event("p", prev_batch)

    render_batch(batch)

    screenshot = args.screenshot if args.screenshot else None
    pl.show(screenshot=screenshot)


if __name__ == "__main__":
    main()
