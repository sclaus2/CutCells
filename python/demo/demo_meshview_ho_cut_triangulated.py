#!/usr/bin/env python3
# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier: MIT
"""
Demo: MeshView + create_level_set + cut(mesh, ls) triangulated straight output.

This is the triangulated counterpart of demo_meshview_ho_cut.py:
1. build a pyvista tetrahedral mesh,
2. convert it to a MeshView,
3. build one polynomial level set with create_level_set(...),
4. call cut(mesh, ls),
5. select inside/interface/outside HOMeshPart objects,
6. triangulate the straight interface and straight cut-volume output,
7. test whether the triangulation-created edges carry a root of the level set,
8. write VTU files and plot the triangulated result in pyvista.

The edge check is important because a triangulation diagonal that crosses the
zero level set can later matter for quadrature mappings and Jacobian signs.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import numpy as np

import cutcells


CENTER = np.array([0.1, -0.05, 0.0], dtype=np.float64)
RADIUS = 0.55
VTK_TETRA = 10

EDGE_PATTERNS = {
    5: ((0, 1), (1, 2), (2, 0)),  # triangle
    9: ((0, 1), (1, 2), (2, 3), (3, 0)),  # quadrilateral
    10: ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)),  # tetrahedron
    13: (
        (0, 1),
        (1, 2),
        (2, 0),
        (3, 4),
        (4, 5),
        (5, 3),
        (0, 3),
        (1, 4),
        (2, 5),
    ),  # prism / wedge
}

TAG_LABELS = {
    int(cutcells.EdgeRootTag.no_root.value): "no_root",
    int(cutcells.EdgeRootTag.one_root.value): "one_root",
    int(cutcells.EdgeRootTag.multiple_roots.value): "multiple_roots",
    int(cutcells.EdgeRootTag.zero.value): "zero",
}


def phi_batch(X: np.ndarray) -> np.ndarray:
    return np.sqrt(
        (X[0] - CENTER[0]) ** 2
        + (X[1] - CENTER[1]) ** 2
        + (X[2] - CENTER[2]) ** 2
    ) - RADIUS


def structured_tetra_mesh(x0, y0, z0, x1, y1, z1, nx, ny, nz):
    xs = np.linspace(x0, x1, num=nx)
    ys = np.linspace(y0, y1, num=ny)
    zs = np.linspace(z0, z1, num=nz)
    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing="ij")
    points = np.c_[xx.ravel(), yy.ravel(), zz.ravel()].astype(np.float64, copy=False)

    def vid(i: int, j: int, k: int) -> int:
        return i + nx * (j + ny * k)

    tet_list = []
    for k in range(nz - 1):
        for j in range(ny - 1):
            for i in range(nx - 1):
                v000 = vid(i, j, k)
                v100 = vid(i + 1, j, k)
                v010 = vid(i, j + 1, k)
                v110 = vid(i + 1, j + 1, k)
                v001 = vid(i, j, k + 1)
                v101 = vid(i + 1, j, k + 1)
                v011 = vid(i, j + 1, k + 1)
                v111 = vid(i + 1, j + 1, k + 1)

                tet_list.extend(
                    (
                        [v000, v100, v110, v111],
                        [v000, v110, v010, v111],
                        [v000, v010, v011, v111],
                        [v000, v011, v001, v111],
                        [v000, v001, v101, v111],
                        [v000, v101, v100, v111],
                    )
                )

    connectivity = np.asarray(tet_list, dtype=np.int32).reshape(-1)
    offsets = np.arange(0, connectivity.size + 1, 4, dtype=np.int32)
    cell_types = np.full(offsets.size - 1, VTK_TETRA, dtype=np.int32)

    # Ensure positive tetra orientation for the affine maps used later.
    for cell_id in range(cell_types.size):
        start = offsets[cell_id]
        tet = connectivity[start : start + 4].copy()
        verts = points[tet]
        J = np.column_stack((verts[1] - verts[0], verts[2] - verts[0], verts[3] - verts[0]))
        if np.linalg.det(J) < 0.0:
            connectivity[start + 1], connectivity[start + 2] = (
                connectivity[start + 2],
                connectivity[start + 1],
            )

    mesh = cutcells.MeshView(points, connectivity, offsets, cell_types, tdim=3)
    return mesh, points, connectivity, offsets, cell_types


def mesh_to_pyvista(points: np.ndarray, connectivity: np.ndarray, offsets: np.ndarray, cell_types: np.ndarray, pv):
    cells = np.empty(cell_types.size * 5, dtype=np.int64)
    cells[0::5] = 4
    cells.reshape(-1, 5)[:, 1:] = connectivity.reshape(-1, 4)
    return pv.UnstructuredGrid(cells, np.asarray(cell_types, dtype=np.uint8), points)


def cutmesh_to_pyvista(cut_mesh, pv):
    return pv.UnstructuredGrid(
        np.array(cut_mesh.cells, dtype=np.int64, copy=True),
        np.array(cut_mesh.vtk_types, dtype=np.uint8, copy=True),
        np.array(cut_mesh.vertex_coords, dtype=np.float64, copy=True),
    )


def tetra_parent_vertices(mesh, cell_id: int) -> np.ndarray:
    connectivity = np.asarray(mesh.connectivity, dtype=np.int32)
    offsets = np.asarray(mesh.offsets, dtype=np.int32)
    coordinates = np.asarray(mesh.coordinates, dtype=np.float64)
    start = int(offsets[cell_id])
    end = int(offsets[cell_id + 1])
    vertex_ids = connectivity[start:end]
    return np.array(coordinates[vertex_ids], dtype=np.float64, copy=True)


def physical_to_parent_reference_tetra(parent_vertices: np.ndarray, x_phys: np.ndarray) -> np.ndarray:
    x0 = parent_vertices[0]
    J = np.column_stack(
        (
            parent_vertices[1] - x0,
            parent_vertices[2] - x0,
            parent_vertices[3] - x0,
        )
    )
    xi = np.linalg.solve(J, np.asarray(x_phys, dtype=np.float64) - x0)
    return np.array(xi, dtype=np.float64, copy=True)


def edge_key(parent_cell_id: int, xa: np.ndarray, xb: np.ndarray, digits: int = 12):
    pa = tuple(np.round(np.asarray(xa, dtype=np.float64), digits))
    pb = tuple(np.round(np.asarray(xb, dtype=np.float64), digits))
    return (int(parent_cell_id), tuple(sorted((pa, pb))))


def collect_edges_by_parent(cut_mesh) -> dict:
    points = np.asarray(cut_mesh.vertex_coords, dtype=np.float64)
    cells = np.asarray(cut_mesh.cells, dtype=np.int64)
    vtk_types = np.asarray(cut_mesh.vtk_types, dtype=np.int64)
    parent_map = np.asarray(cut_mesh.parent_map, dtype=np.int32)

    if parent_map.size != vtk_types.size:
        raise RuntimeError("cut_mesh.parent_map must contain one entry per output cell.")

    edges = {}
    cursor = 0
    for cell_id, vtk_type in enumerate(vtk_types):
        num_verts = int(cells[cursor])
        local_ids = cells[cursor + 1 : cursor + 1 + num_verts]
        cursor += num_verts + 1
        parent_cell_id = int(parent_map[cell_id])

        for local_a, local_b in EDGE_PATTERNS[int(vtk_type)]:
            xa = np.array(points[int(local_ids[local_a])], dtype=np.float64, copy=True)
            xb = np.array(points[int(local_ids[local_b])], dtype=np.float64, copy=True)
            edges[edge_key(parent_cell_id, xa, xb)] = {
                "parent_cell_id": parent_cell_id,
                "xa_phys": xa,
                "xb_phys": xb,
            }

    return edges


def classify_new_triangulation_edges(mesh, ls, base_mesh, triangulated_mesh) -> list[dict]:
    base_edges = collect_edges_by_parent(base_mesh)
    tri_edges = collect_edges_by_parent(triangulated_mesh)
    ls_cell_cache: dict[int, object] = {}
    parent_vertex_cache: dict[int, np.ndarray] = {}
    records = []

    for key in sorted(set(tri_edges) - set(base_edges)):
        record = dict(tri_edges[key])
        parent_cell_id = record["parent_cell_id"]

        if parent_cell_id not in ls_cell_cache:
            ls_cell_cache[parent_cell_id] = cutcells.make_cell_level_set(ls, parent_cell_id)
            parent_vertex_cache[parent_cell_id] = tetra_parent_vertices(mesh, parent_cell_id)

        ls_cell = ls_cell_cache[parent_cell_id]
        parent_vertices = parent_vertex_cache[parent_cell_id]
        xa_ref = physical_to_parent_reference_tetra(parent_vertices, record["xa_phys"])
        xb_ref = physical_to_parent_reference_tetra(parent_vertices, record["xb_phys"])

        edge_coeffs = np.asarray(
            cutcells.restrict_edge_bernstein_exact(
                ls_cell.cell_type,
                ls_cell.bernstein_order,
                np.asarray(ls_cell.bernstein_coeffs),
                xa_ref,
                xb_ref,
            ),
            dtype=np.float64,
        )
        tag, split_t = cutcells.classify_edge_roots(
            edge_coeffs,
            zero_tol=1.0e-12,
            sign_tol=1.0e-12,
            max_depth=20,
        )
        tag_value = int(tag.value)
        phi_a = float(
            cutcells.evaluate_bernstein(
                ls_cell.cell_type,
                ls_cell.bernstein_order,
                np.asarray(ls_cell.bernstein_coeffs),
                xa_ref,
            )
        )
        phi_b = float(
            cutcells.evaluate_bernstein(
                ls_cell.cell_type,
                ls_cell.bernstein_order,
                np.asarray(ls_cell.bernstein_coeffs),
                xb_ref,
            )
        )
        has_zero_endpoint = abs(phi_a) <= 1.0e-12 or abs(phi_b) <= 1.0e-12
        has_root = tag_value in {
            int(cutcells.EdgeRootTag.one_root.value),
            int(cutcells.EdgeRootTag.multiple_roots.value),
            int(cutcells.EdgeRootTag.zero.value),
        }
        has_interior_root = (
            tag_value in {
                int(cutcells.EdgeRootTag.one_root.value),
                int(cutcells.EdgeRootTag.multiple_roots.value),
            }
            and not has_zero_endpoint
        )

        record["xa_ref"] = xa_ref
        record["xb_ref"] = xb_ref
        record["phi_a"] = phi_a
        record["phi_b"] = phi_b
        record["edge_root_tag"] = tag_value
        record["edge_has_root"] = int(has_root)
        record["edge_has_zero_endpoint"] = int(has_zero_endpoint)
        record["edge_has_interior_root"] = int(has_interior_root)
        record["green_split_t"] = float(split_t) if split_t is not None else np.nan
        records.append(record)

    return records


def write_line_vtu(path: Path, records: list[dict]) -> None:
    if not records:
        return

    points = []
    connectivity = []
    offsets = []
    root_tags = []
    has_root = []
    parent_ids = []
    split_t = []

    for edge_id, record in enumerate(records):
        points.extend(record["xa_phys"].tolist())
        points.extend(record["xb_phys"].tolist())
        connectivity.extend((2 * edge_id, 2 * edge_id + 1))
        offsets.append(2 * edge_id + 2)
        root_tags.append(record["edge_root_tag"])
        has_root.append(record["edge_has_root"])
        parent_ids.append(record["parent_cell_id"])
        split_t.append(record["green_split_t"])

    point_text = " ".join(f"{value:.16g}" for value in points)
    connectivity_text = " ".join(str(value) for value in connectivity)
    offsets_text = " ".join(str(value) for value in offsets)
    types_text = " ".join("3" for _ in records)  # VTK_LINE
    root_tags_text = " ".join(str(value) for value in root_tags)
    has_root_text = " ".join(str(value) for value in has_root)
    parent_ids_text = " ".join(str(value) for value in parent_ids)
    split_t_text = " ".join(
        "nan" if np.isnan(value) else f"{value:.16g}" for value in split_t
    )

    xml = f"""<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="{2 * len(records)}" NumberOfCells="{len(records)}">
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          {point_text}
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          {connectivity_text}
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
          {offsets_text}
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">
          {types_text}
        </DataArray>
      </Cells>
      <CellData Scalars="edge_root_tag">
        <DataArray type="Int32" Name="edge_root_tag" format="ascii">
          {root_tags_text}
        </DataArray>
        <DataArray type="Int32" Name="edge_has_root" format="ascii">
          {has_root_text}
        </DataArray>
        <DataArray type="Int32" Name="parent_cell_id" format="ascii">
          {parent_ids_text}
        </DataArray>
        <DataArray type="Float64" Name="green_split_t" format="ascii">
          {split_t_text}
        </DataArray>
      </CellData>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
"""
    path.write_text(xml)


def summarize_records(records: list[dict]) -> dict[str, int]:
    counts = Counter(TAG_LABELS[record["edge_root_tag"]] for record in records)
    return dict(sorted(counts.items()))


def count_interior_root_edges(records: list[dict]) -> int:
    return int(sum(record["edge_has_interior_root"] for record in records))


def count_zero_endpoint_edges(records: list[dict]) -> int:
    return int(sum(record["edge_has_zero_endpoint"] for record in records))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Triangulated HO straight-output demo with triangulation-edge root checks."
    )
    parser.add_argument(
        "--n",
        type=int,
        default=8,
        help="Structured tetrahedral mesh resolution parameter.",
    )
    parser.add_argument(
        "--degree",
        type=int,
        default=2,
        help="Polynomial degree for create_level_set(...).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("demo_meshview_ho_cut_triangulated_output"),
        help="Directory for VTU files.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip interactive pyvista plotting.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Creating tetrahedral pyvista mesh with n={args.n} ...")
    mesh, points, connectivity, offsets, cell_types = structured_tetra_mesh(
        -1.0, -1.0, -1.0,
         1.0,  1.0,  1.0,
        args.n, args.n, args.n,
    )
    print(f"  tetra cells={mesh.num_cells()}, points={mesh.num_nodes()}")
    print(f"  MeshView gdim={mesh.gdim}, tdim={mesh.tdim}")

    print(f"Building polynomial level set with degree={args.degree} ...")
    ls = cutcells.create_level_set(mesh, phi_batch, degree=args.degree, name="phi")
    print(f"  level-set dofs={ls.mesh_data.num_dofs()}")

    print("Running cut(mesh, ls) ...")
    result = cutcells.cut(mesh, ls)
    negative = result["phi < 0"]
    interface = result["phi = 0"]
    positive = result["phi > 0"]

    print(f"  num_cut_cells={result.num_cut_cells}")
    print(f"  phi < 0  : cut={negative.num_cut_cells}, uncut={negative.num_uncut_cells}")
    print(f"  phi = 0  : cut={interface.num_cut_cells}, uncut={interface.num_uncut_cells}")
    print(f"  phi > 0  : cut={positive.num_cut_cells}, uncut={positive.num_uncut_cells}")

    print("Building straight visualization meshes from HOMeshPart ...")
    inside_full = negative.visualization_mesh(mode="full")
    inside_cut = negative.visualization_mesh(mode="cut_only")
    outside_full = positive.visualization_mesh(mode="full")
    outside_cut = positive.visualization_mesh(mode="cut_only")
    interface_mesh = interface.visualization_mesh(mode="cut_only")

    print(f"  inside full cells={len(np.asarray(inside_full.types))}")
    print(f"  inside cut-only cells={len(np.asarray(inside_cut.types))}")
    print(f"  interface cells={len(np.asarray(interface_mesh.types))}")
    print(f"  outside full cells={len(np.asarray(outside_full.types))}")
    print(f"  outside cut-only cells={len(np.asarray(outside_cut.types))}")

    print("Building reference meshes for the new-edge check ...")
    inside_full_base = negative.visualization_mesh(mode="full")
    inside_cut_base = negative.visualization_mesh(mode="cut_only")
    outside_full_base = positive.visualization_mesh(mode="full")
    outside_cut_base = positive.visualization_mesh(mode="cut_only")
    interface_base = interface.visualization_mesh(mode="cut_only")

    print("Checking whether triangulation-created edges carry level-set roots ...")
    edge_checks = {
        "phi_negative_full": classify_new_triangulation_edges(
            mesh, ls, inside_full_base, inside_full
        ),
        "phi_negative_cut_only": classify_new_triangulation_edges(
            mesh, ls, inside_cut_base, inside_cut
        ),
        "phi_interface": classify_new_triangulation_edges(
            mesh, ls, interface_base, interface_mesh
        ),
        "phi_positive_full": classify_new_triangulation_edges(
            mesh, ls, outside_full_base, outside_full
        ),
        "phi_positive_cut_only": classify_new_triangulation_edges(
            mesh, ls, outside_cut_base, outside_cut
        ),
    }

    for name, records in edge_checks.items():
        print(
            f"  {name}: new_edges={len(records)}, "
            f"interior_root_edges={count_interior_root_edges(records)}, "
            f"zero_endpoint_edges={count_zero_endpoint_edges(records)}, "
            f"tags={summarize_records(records)}"
        )

    for stem, records in edge_checks.items():
        if not records:
            continue
        path = args.output_dir / f"{stem}_triangulation_new_edges.vtu"
        write_line_vtu(path, records)
        print(f"  wrote {path}")

    mesh_outputs = [
        (negative, args.output_dir / "phi_negative_full_triangulated.vtu", "full"),
        (negative, args.output_dir / "phi_negative_cut_only_triangulated.vtu", "cut_only"),
        (interface, args.output_dir / "phi_interface_triangulated.vtu", "cut_only"),
        (positive, args.output_dir / "phi_positive_full_triangulated.vtu", "full"),
        (positive, args.output_dir / "phi_positive_cut_only_triangulated.vtu", "cut_only"),
    ]

    for part, path, mode in mesh_outputs:
        part.write_vtu(str(path), mode=mode)
        print(f"  wrote {path}")

    if args.no_plot:
        return

    try:
        import pyvista as pv
    except Exception as exc:
        raise SystemExit(
            "pyvista is required for plotting this demo. "
            "The VTU files were still written successfully.\n"
            f"Import error: {exc}"
        )

    pv_inside_full = cutmesh_to_pyvista(inside_full, pv)
    pv_interface = cutmesh_to_pyvista(interface_mesh, pv)
    pv_outside_full = cutmesh_to_pyvista(outside_full, pv)
    grid = mesh_to_pyvista(points, connectivity, offsets, cell_types, pv)

    plotter = pv.Plotter(shape=(1, 3), title="CutCells HO triangulated straight-output demo")
    plotter.subplot(0, 0)
    plotter.add_title("phi < 0, mode='full', triangulated")
    plotter.add_mesh(grid, style="wireframe", color="lightgrey", opacity=0.18)
    plotter.add_mesh(pv_inside_full, color="steelblue", opacity=0.78, show_edges=True)

    plotter.subplot(0, 1)
    plotter.add_title("phi = 0, mode='cut_only', triangulated")
    plotter.add_mesh(grid, style="wireframe", color="lightgrey", opacity=0.12)
    plotter.add_mesh(pv_interface, color="crimson", opacity=0.95, show_edges=True)

    plotter.subplot(0, 2)
    plotter.add_title("phi > 0, mode='full', triangulated")
    plotter.add_mesh(grid, style="wireframe", color="lightgrey", opacity=0.18)
    plotter.add_mesh(pv_outside_full, color="burlywood", opacity=0.78, show_edges=True)

    plotter.link_views()
    plotter.show()


if __name__ == "__main__":
    main()
