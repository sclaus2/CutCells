# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier:    MIT
"""
Demo: MeshView + LevelSetFunction Python wrappers
==================================================

Demonstrates:
  * Building a ``cutcells.MeshView`` from a pyvista ``UnstructuredGrid``.
  * Building a ``cutcells.LevelSetFunction`` from a plain Python callable.
  * Calling ``cutcells.cut_mesh_view`` to obtain cut meshes for
      – the interface  (phi = 0)
      – the inside     (phi < 0)
  * Visualising the results with pyvista.

The background mesh is a 2-D Delaunay triangulation of a square [-1, 1]^2.
The level-set is a circle of radius 0.7 centred at the origin:
    phi(x) = sqrt(x0^2 + x1^2) - 0.7
"""

import argparse

import numpy as np

try:
    import pyvista as pv
except ImportError as exc:
    raise SystemExit(
        "pyvista is required for this demo.  "
        "Install with:  python -m pip install pyvista\n"
        f"Import error: {exc}"
    )

import cutcells
from cutcells import (
    mesh_from_pyvista,
    rectangle_triangle_mesh,
)


# ---------------------------------------------------------------------------
# Level-set definition
# ---------------------------------------------------------------------------

RADIUS = 0.7
CENTER = np.array([0.0, 0.0, 0.0])


def circle_phi(x: np.ndarray) -> float:
    """Signed-distance to a circle of radius RADIUS centred at CENTER.

    ``x`` is a 1-D NumPy array of length 3 (pyvista always stores 3-D coords).
    Returns a float: negative inside the circle, positive outside.
    """
    return float(np.sqrt(np.sum((x - CENTER) ** 2)) - RADIUS)


# ---------------------------------------------------------------------------
# Helper: build MeshView from a pyvista UnstructuredGrid
# ---------------------------------------------------------------------------
#
# NOTE: Using mesh_from_pyvista from cutcells.mesh_utils module


# ---------------------------------------------------------------------------
# Mesh generation
# ---------------------------------------------------------------------------
#
# NOTE: Using rectangle_triangle_mesh from cutcells.mesh_utils module


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Demo: cut_mesh_view with MeshView and LevelSetFunction"
    )
    parser.add_argument(
        "--n",
        type=int,
        default=30,
        help="Grid points per axis (default: 30)",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip interactive visualisation (useful for CI / off-screen runs)",
    )
    args = parser.parse_args()

    N = int(args.n)

    # ------------------------------------------------------------------
    # 1.  Build the background pyvista grid and wrap it in a MeshView
    # ------------------------------------------------------------------
    print(f"Creating {N}×{N} Delaunay triangle mesh …")
    grid = rectangle_triangle_mesh(-1.0, -1.0, 1.0, 1.0, N, N)
    print(f"  nodes: {grid.n_points},  cells: {grid.n_cells}")

    mesh_view = mesh_from_pyvista(grid, tdim=2)
    print(
        f"  MeshView gdim={mesh_view.gdim}, tdim={mesh_view.tdim}, "
        f"num_nodes={mesh_view.num_nodes()}, num_cells={mesh_view.num_cells()}"
    )

    # ------------------------------------------------------------------
    # 2.  Build a LevelSetFunction from a Python callable
    # ------------------------------------------------------------------
    # The callable receives a length-3 NumPy array (x) and may optionally
    # accept a second argument cell_id (ignored here).
    level_set = cutcells.LevelSetFunction(
        value=circle_phi,
        gdim=mesh_view.gdim,  # 3 (pyvista always stores 3-D coords)
    )
    print(
        f"\nLevelSetFunction  has_value={level_set.has_value()}, "
        f"has_nodal_values={level_set.has_nodal_values()}"
    )

    # ------------------------------------------------------------------
    # 3.  Cut the mesh  –  interface (phi = 0)
    # ------------------------------------------------------------------
    print("\nCutting mesh for phi=0 (interface) …")
    cut_interface = cutcells.cut_mesh_view(
        mesh_view, level_set, "phi=0", triangulate=True
    )
    print(
        f"  Interface cut mesh: {len(cut_interface.vtk_types)} cells, "
        f"{len(cut_interface.vertex_coords)} vertices"
    )

    # ------------------------------------------------------------------
    # 4.  Cut the mesh  –  interior (phi < 0)
    # ------------------------------------------------------------------
    print("Cutting mesh for phi<0 (inside) …")
    cut_inside = cutcells.cut_mesh_view(mesh_view, level_set, "phi<0", triangulate=True)
    print(
        f"  Inside cut mesh:    {len(cut_inside.vtk_types)} cells, "
        f"{len(cut_inside.vertex_coords)} vertices"
    )

    # ------------------------------------------------------------------
    # 5.  (Optional) also build nodal values and use LevelSetFunction
    #     with nodal_values instead of value callable  – same result
    # ------------------------------------------------------------------
    print("\nBuilding nodal level-set values manually …")
    ls_nodal = np.array([circle_phi(p) for p in grid.points], dtype=np.float64)
    level_set_nodal = cutcells.LevelSetFunction(
        nodal_values=ls_nodal,
        gdim=mesh_view.gdim,
    )
    print(
        f"  LevelSetFunction (nodal)  has_nodal_values={level_set_nodal.has_nodal_values()}"
    )

    cut_inside_nodal = cutcells.cut_mesh_view(
        mesh_view, level_set_nodal, "phi<0", triangulate=True
    )
    assert len(cut_inside_nodal.vtk_types) == len(cut_inside.vtk_types), (
        "Nodal and callable level-set should give identical cut meshes"
    )
    print("  Nodal and callable results agree ✓")

    # ------------------------------------------------------------------
    # 6.  Runtime Quadrature and Physical Mapping
    # ------------------------------------------------------------------
    print("\nComputing runtime quadrature for phi<0 …")
    order = 2
    # Use the flat arrays from the grid
    points_flat = np.asarray(grid.points, dtype=np.float64).ravel()
    connectivity = np.asarray(grid.cell_connectivity, dtype=np.int32)
    offsets = np.asarray(grid.offset, dtype=np.int32)
    celltypes = np.asarray(grid.celltypes, dtype=np.int32)

    # Calculate quadrature rules in reference space for phi<0
    q_rules_inside = cutcells.runtime_quadrature(
        ls_nodal,
        points_flat,
        connectivity,
        offsets,
        celltypes,
        "phi<0",
        triangulate=True,
        order=order,
    )
    print(f"  Inside: Generated {len(q_rules_inside.weights)} total quadrature points.")

    # Calculate quadrature rules in reference space for phi=0
    print("Computing runtime quadrature for phi=0 …")
    q_rules_interface = cutcells.runtime_quadrature(
        ls_nodal,
        points_flat,
        connectivity,
        offsets,
        celltypes,
        "phi=0",
        triangulate=True,
        order=order,
    )
    print(
        f"  Interface: Generated {len(q_rules_interface.weights)} total quadrature points."
    )

    # Map to physical space
    print("  Mapping quadrature points to physical space …")
    q_points_inside_phys = cutcells.physical_points(
        q_rules_inside, points_flat, connectivity, offsets, celltypes
    ).reshape(-1, 3)

    q_points_interface_phys = cutcells.physical_points(
        q_rules_interface, points_flat, connectivity, offsets, celltypes
    ).reshape(-1, 3)

    # ------------------------------------------------------------------
    # 7.  Visualise
    # ------------------------------------------------------------------
    if args.no_plot:
        print("\n--no-plot specified, skipping visualisation.")
        return

    # -- wrap cut results as pyvista grids --
    pv_interface = pv.UnstructuredGrid(
        cut_interface.cells,
        cut_interface.vtk_types,
        np.asarray(cut_interface.vertex_coords, dtype=np.float64),
    )
    pv_inside = pv.UnstructuredGrid(
        cut_inside.cells,
        cut_inside.vtk_types,
        np.asarray(cut_inside.vertex_coords, dtype=np.float64),
    )

    # -- inside background cells (entirely inside the circle) --
    points_flat = np.asarray(grid.points, dtype=np.float64).ravel()
    connectivity = np.asarray(grid.cell_connectivity, dtype=np.int32)
    offsets = np.asarray(grid.offset, dtype=np.int32)
    celltypes = np.asarray(grid.celltypes, dtype=np.int32)
    inside_ids = cutcells.locate_cells(
        ls_nodal, points_flat, connectivity, offsets, celltypes, "phi<0"
    )
    pv_inside_bg = grid.extract_cells(inside_ids)

    # -- combined inside view: background inside cells + cut cells --
    pv_combined = pv_inside_bg.merge(pv_inside)

    # -- wrap quadrature points as pyvista points --
    pv_q_points_inside = pv.PolyData(q_points_inside_phys)
    pv_q_points_interface = pv.PolyData(q_points_interface_phys)

    # Plot
    pl = pv.Plotter(shape=(1, 2), title="CutCells – MeshView + LevelSetFunction demo")

    pl.subplot(0, 0)
    pl.add_title("Interface (phi = 0) + Quad Points")
    pl.add_mesh(grid, style="wireframe", color="lightgrey", opacity=0.4)
    pl.add_mesh(pv_interface, show_edges=True, color="red")
    pl.add_mesh(
        pv_q_points_interface,
        color="yellow",
        point_size=10,
        render_points_as_spheres=True,
        label=f"Interface Points ({len(q_rules_interface.weights)})",
    )
    pl.add_legend()
    pl.view_xy()

    pl.subplot(0, 1)
    pl.add_title("Inside domain (phi < 0) + Quad Points")
    pl.add_mesh(grid, style="wireframe", color="lightgrey", opacity=0.4)
    pl.add_mesh(pv_combined, show_edges=True, color="steelblue", opacity=0.8)
    pl.add_mesh(
        pv_q_points_inside,
        color="orange",
        point_size=10,
        render_points_as_spheres=True,
        label=f"Inside Points ({len(q_rules_inside.weights)})",
    )
    pl.add_legend()
    pl.view_xy()

    pl.show()


if __name__ == "__main__":
    main()
