# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""The types of the front end.

``cutcells.cut(mesh, level_sets)`` classifies every cell by every level set, Pk
(LevelSetFunction with dof values) or analytic (AnalyticLevelSet), by the level
sets' own bounds, and returns a CutResult. ``result["phi1 < 0 and phi2 = 0"]``
selects a MeshPart, whose ``quadrature``, ``visualization_mesh`` and
``write_vtu`` come from a backend: quadrays by default (options
QuadraysOptions), or the lookup tables on Pk-iso-P1 templates,
``backend="lut"`` (options LutOptions), for the whole result or per call::

    import cutcells

    result = cutcells.cut(mesh, cutcells.analytic_sphere([0, 0, 0], 0.7))
    rules = result["phi < 0"].quadrature(order=5)
    result["phi = 0"].write_vtu("sphere.vtu", mode="cut_only")
    straight = result["phi < 0"].quadrature(order=2, backend="lut",
                                            options=cutcells.LutOptions(template_order=3))

``cut``, ``cut_float32`` and ``cut_float64`` here are deprecated aliases of
``cutcells.cut``.
"""

import warnings as _warnings

from ._cutcellscpp import part as _part

CutResult = _part.CutResult
CutResult_float32 = _part.CutResult_float32
CutResult_float64 = _part.CutResult_float64
MeshPart = _part.MeshPart
MeshPart_float32 = _part.MeshPart_float32
MeshPart_float64 = _part.MeshPart_float64


def _retired(name):
    target = getattr(_part, name)

    def alias(*args, **kwargs):
        _warnings.warn(f"cutcells.part.{name} is deprecated; use cutcells.cut", DeprecationWarning, stacklevel=2)
        return target(*args, **kwargs)

    alias.__name__ = alias.__qualname__ = name
    alias.__doc__ = f"Deprecated: use cutcells.cut.\n\n{target.__doc__}"
    return alias


cut = _retired("cut")
cut_float32 = _retired("cut_float32")
cut_float64 = _retired("cut_float64")
