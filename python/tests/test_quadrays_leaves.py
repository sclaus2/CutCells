# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Leaf cells of the quadrays decomposition: geometry, and VTK's reading of them.

The VTK checks need the vtk package and are skipped without it:
- node order: a leaf of a box cut by a plane is an affine image of the
  reference cell, so VTK's own interpolation of the Lagrange cell must be
  affine in the parametric coordinates; a wrong node order breaks that;
- coverage: VTK's cell sizes of the leaves converge to the exact volume.
"""

import math

import numpy as np
import pytest

import cutcells

from test_quadrays import CENTRE, RADIUS, box_mesh, sphere

UNIT_HEX = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]],
                    dtype=np.float64)


def box_leaf(selection: str, degree: int):
    """The leaf of the unit hexahedron cut by the plane x = 0.6: a box or a square."""
    # degree-1 tensor Bernstein coefficients are vertex values, in CutCells' order (i * 2 + j) * 2 + k
    values = np.array([i - 0.6 for i in (0, 1) for j in (0, 1) for k in (0, 1)])
    return cutcells.quadrays_cell_leaves(cutcells.CellType.hexahedron, UNIT_HEX, 1, values, selection,
                                         leaf_degree=degree)


def test_leaves_of_a_sphere():
    mesh = box_mesh("tet", 4)
    result = cutcells.cut(mesh, cutcells.create_level_set(mesh, sphere, degree=2, name="phi"))
    surface = cutcells.quadrays_leaves(result["phi = 0"], degree=3)
    volume = cutcells.quadrays_leaves(result["phi < 0"], degree=3, mode="full")
    assert surface.n_cells() > 0 and set(np.unique(surface.vtk_types)) == {70}
    assert set(np.unique(volume.vtk_types)) == {10, 72}  # uncut cells as linear tetrahedra
    assert volume.offsets[-1] == len(volume.connectivity)
    distance = np.linalg.norm(surface.points - CENTRE, axis=1)
    np.testing.assert_allclose(distance, RADIUS, atol=1e-14)
    assert np.all(np.linalg.norm(volume.points - CENTRE, axis=1) <= RADIUS + 1e-12)


def test_write_leaves(tmp_path):
    leaves = box_leaf("phi < 0", 3)
    path = tmp_path / "leaves.vtu"
    cutcells.write_quadrays_leaves(str(path), leaves)
    text = path.read_text()
    assert 'version="2.2"' in text and "HigherOrderDegrees" in text


def evaluate(cell, pcoords, vtk):
    x = [0.0, 0.0, 0.0]
    weights = [0.0] * cell.GetNumberOfPoints()
    cell.EvaluateLocation(vtk.reference(0), list(pcoords), x, weights)
    return np.array(x)


def read(path, vtk):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    return reader.GetOutput()


@pytest.mark.parametrize("selection", ["phi < 0", "phi = 0"])
@pytest.mark.parametrize("degree", [1, 2, 3, 4])
def test_vtk_node_order(tmp_path, selection, degree):
    vtk = pytest.importorskip("vtk")
    leaves = box_leaf(selection, degree)
    assert leaves.n_cells() == 1
    path = tmp_path / "leaf.vtu"
    cutcells.write_quadrays_leaves(str(path), leaves)
    cell = read(path, vtk).GetCell(0)
    dim = cell.GetCellDimension()
    origin = evaluate(cell, [0.0, 0.0, 0.0], vtk)
    axes = [evaluate(cell, np.eye(3)[d], vtk) - origin for d in range(dim)]
    rng = np.random.default_rng(0)
    for _ in range(20):
        pc = np.zeros(3)
        pc[:dim] = rng.random(dim)
        expected = origin + sum(pc[d] * axes[d] for d in range(dim))
        np.testing.assert_allclose(evaluate(cell, pc, vtk), expected, atol=1e-12)


def test_vtk_coverage(tmp_path):
    vtk = pytest.importorskip("vtk")
    mesh = box_mesh("tet", 8)
    result = cutcells.cut(mesh, cutcells.create_level_set(mesh, sphere, degree=2, name="phi"))
    exact = 4.0 / 3.0 * math.pi * RADIUS**3
    differences = []
    for degree in (2, 4):
        path = tmp_path / f"volume_{degree}.vtu"
        cutcells.write_quadrays_leaves(str(path), cutcells.quadrays_leaves(result["phi < 0"], degree=degree))
        sizes = vtk.vtkCellSizeFilter()
        sizes.SetInputData(read(path, vtk))
        sizes.Update()
        data = sizes.GetOutput().GetCellData().GetArray("Volume")
        values = np.array([data.GetValue(i) for i in range(data.GetNumberOfTuples())])
        assert np.all(values >= 0)
        differences.append(abs(values.sum() - exact) / exact)
    # VTK measures curved cells by linear subdivision: the difference falls like 1/p^2
    assert differences[1] < 0.5 * differences[0] and differences[1] < 2e-3
