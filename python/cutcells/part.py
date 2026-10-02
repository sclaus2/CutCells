# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""The front end without AdaptCell.

``cut(mesh, level_sets)`` classifies every cell by every level set, Pk
(LevelSetFunction with dof values) or analytic (AnalyticLevelSet), by the level
sets' own bounds. ``result["phi1 < 0 and phi2 = 0"]`` selects a MeshPart, whose
``quadrature``, ``visualization_mesh`` and ``write_vtu`` come from a backend::

    import cutcells

    result = cutcells.part.cut(mesh, cutcells.analytic_sphere([0, 0, 0], 0.7))
    rules = result["phi < 0"].quadrature(order=5, backend="quadrays")
    result["phi = 0"].write_vtu("sphere.vtu", mode="cut_only")

It sits beside HOCutResult and HOMeshPart until it replaces them (phase 5).
"""

from ._cutcellscpp import part as _part

CutResult = _part.CutResult
CutResult_float32 = _part.CutResult_float32
CutResult_float64 = _part.CutResult_float64
MeshPart = _part.MeshPart
MeshPart_float32 = _part.MeshPart_float32
MeshPart_float64 = _part.MeshPart_float64
cut = _part.cut
cut_float32 = _part.cut_float32
cut_float64 = _part.cut_float64
