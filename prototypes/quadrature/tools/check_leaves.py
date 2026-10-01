# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Checks and renders the prototype's leaf-cell .vtu files with VTK itself.

Needs vtk (and pyvista for `render`), e.g. the vtk-env environment:

    python tools/check_leaves.py ordering build/ordering
    python tools/check_leaves.py measure build/leaves/<file>.vtu <exact total>
    python tools/check_leaves.py render <interface.vtu> <volume.vtu> <out.png> [<cut_cells.vtu>]

ordering: every file written by test_leaf_ordering holds one Lagrange cell whose
nodes sit at their own parametric coordinates, so a correct node order makes VTK's
interpolation the identity: x(pcoords) == pcoords at random points.

measure: total size of all cells as VTK computes it (cells are triangulated
through their nodes, so curved cells come out slightly small), compared with an
exact total.
"""

import glob
import sys

import numpy as np
import vtk


def read(path):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(path)
    reader.Update()
    return reader.GetOutput()


def ordering(directory):
    rng = np.random.default_rng(0)
    failures = 0
    for path in sorted(glob.glob(f"{directory}/lagrange_*_p*.vtu")):
        grid = read(path)
        cell = grid.GetCell(0)
        dim = cell.GetCellDimension()
        worst = 0.0
        for _ in range(20):
            pc = [*rng.random(dim), *([0.0] * (3 - dim))]
            x = [0.0, 0.0, 0.0]
            weights = [0.0] * cell.GetNumberOfPoints()
            cell.EvaluateLocation(vtk.reference(0), pc, x, weights)
            worst = max(worst, max(abs(a - b) for a, b in zip(x, pc)))
        ok = worst < 1e-12
        failures += not ok
        print(f"{path.split('/')[-1]:24s} {cell.GetClassName():28s} max |x - pcoords| = {worst:.1e} {'ok' if ok else 'WRONG ORDER'}")
    return failures


def measure(path, exact):
    grid = read(path)
    sizes = vtk.vtkCellSizeFilter()
    sizes.SetInputData(grid)
    sizes.Update()
    data = sizes.GetOutput().GetCellData()
    total = 0.0
    negative = 0
    for name in ("Area", "Volume"):
        values = data.GetArray(name)
        if values is None:
            continue
        v = np.array([values.GetValue(i) for i in range(values.GetNumberOfTuples())])
        negative += int(np.sum(v < 0))
        total += float(np.sum(np.abs(v)))
    types = {grid.GetCellType(i) for i in range(grid.GetNumberOfCells())}
    print(f"{path.split('/')[-1]}: {grid.GetNumberOfCells()} cells, types {sorted(types)}, "
          f"total {total:.10f}, exact {exact:.10f}, rel. difference {abs(total - exact) / exact:.1e}, "
          f"negative-size cells {negative}")
    return 0


def render(interface_path, volume_path, png, cells_path=None):
    """Left: the whole sphere, every interface leaf in its own colour. Middle and
    right: one cut tet with its interface leaves and its phi < 0 leaves (shrunk)."""
    import pyvista as pv

    pv.OFF_SCREEN = True
    surface = pv.read(interface_path)
    volume = pv.read(volume_path)
    cells = pv.read(cells_path) if cells_path else None

    def coloured(mesh):
        mesh.cell_data["leaf"] = (np.arange(mesh.n_cells) * 2654435761 % 4096) / 4096.0
        return mesh

    def curved(mesh):
        return coloured(mesh).extract_surface(nonlinear_subdivision=4)

    def leaf_edges(subdivided):
        return subdivided.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False,
                                                non_manifold_edges=False)

    # one cut cell near a chosen point of the sphere, with a rich decomposition
    centre = np.array([0.0123, -0.0371, 0.0217])
    normal = np.array([1.0, -1.0, 0.6]) / np.linalg.norm([1.0, -1.0, 0.6])
    point = centre + 0.7 * normal
    parents = surface.cell_data["parent"]
    near = np.linalg.norm(surface.cell_centers().points - point, axis=1) < 0.25
    candidates, counts = np.unique(parents[near], return_counts=True)
    parent = candidates[np.argmax(counts)]

    def of_parent(mesh):
        return mesh.extract_cells(np.where(mesh.cell_data["parent"] == parent)[0])

    side = np.cross(normal, [0.0, 0.0, 1.0])
    view = normal + 0.9 * side / np.linalg.norm(side) + np.array([0.0, 0.0, 0.5])

    plotter = pv.Plotter(shape=(1, 3), window_size=(2100, 700), border=False)
    plotter.subplot(0, 0)
    plotter.add_text("interface leaves of all cut tets (Lagrange quadrilaterals, degree 3)", font_size=9)
    whole = curved(surface)
    plotter.add_mesh(whole, scalars="leaf", cmap="tab20", show_scalar_bar=False, smooth_shading=True)
    plotter.add_mesh(leaf_edges(whole), color="black", line_width=1)
    plotter.view_vector(normal)

    plotter.subplot(0, 1)
    plotter.add_text(f"cut tet {parent}: interface leaves", font_size=9)
    patch = curved(of_parent(surface))
    plotter.add_mesh(patch, scalars="leaf", cmap="tab20", show_scalar_bar=False, smooth_shading=True)
    plotter.add_mesh(leaf_edges(patch), color="black", line_width=3)
    if cells is not None:
        plotter.add_mesh(of_parent(cells).extract_all_edges(), color="black", line_width=2)
    plotter.view_vector(view)

    plotter.subplot(0, 2)
    plotter.add_text(f"cut tet {parent}: phi < 0 leaves (Lagrange hexahedra, shrunk 0.8)", font_size=9)
    leaves = coloured(of_parent(volume)).shrink(0.8)
    plotter.add_mesh(leaves.extract_surface(nonlinear_subdivision=4), scalars="leaf", cmap="tab20",
                     show_scalar_bar=False)
    if cells is not None:
        plotter.add_mesh(of_parent(cells).extract_all_edges(), color="black", line_width=2)
    plotter.view_vector(view)
    plotter.screenshot(png)
    print(f"wrote {png}: cut tet {parent}, {patch.n_cells and of_parent(surface).n_cells} interface leaves, "
          f"{of_parent(volume).n_cells} volume leaves")
    return 0


if __name__ == "__main__":
    command = sys.argv[1]
    if command == "ordering":
        sys.exit(ordering(sys.argv[2]))
    if command == "measure":
        sys.exit(measure(sys.argv[2], float(sys.argv[3])))
    if command == "render":
        sys.exit(render(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5] if len(sys.argv) > 5 else None))
    raise SystemExit(__doc__)
