# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""ShapeForest shapes as analytic level sets for the quadrays backend.

ShapeForest is optional. Its shapes and expressions are lowered to tapes
(register code) by ShapeForest itself, and the tapes reach CutCells as arrays
(:func:`cutcells.analytic_level_set_from_tape`), so CutCells does not depend on
ShapeForest. A tape object needs nothing from ShapeForest at all::

    import shapeforest as sf
    import cutcells
    from cutcells import shapeforest as csf

    phi = csf.analytic_level_set(sf.sphere(0.7).translate(0.1, 0.0, 0.0))
    result = cutcells.cut(mesh, phi)
    rules = result["phi = 0"].quadrature(order=5)
"""

import numpy as np

from ._cutcellscpp import analytic_level_set_from_tape


def to_tape(shape):
    """The ShapeForest tape of a shape (``shape.to_tape()``), an expression
    (``shapeforest.expr_to_tape``) or a tape."""
    if hasattr(shape, "instrs") and hasattr(shape, "n_regs"):
        return shape
    if hasattr(shape, "to_tape"):
        return shape.to_tape()
    from shapeforest.ir.expr_to_tape import expr_to_tape

    return expr_to_tape(shape)


def _parameter_values(tape, parameters):
    parameters = {} if parameters is None else parameters
    missing = [name for name, _ in tape.extra_inputs if name not in parameters]
    if missing:
        raise ValueError(f"the tape has parameters without values: {missing}")
    return [float(parameters[name]) for name, _ in tape.extra_inputs]


def tape_arrays(shape, parameters=None):
    """The arguments of :func:`cutcells.analytic_level_set_from_tape` for a
    shape, an expression or a tape; ``parameters`` gives values to the tape's
    named parameters."""
    tape = to_tape(shape)
    values = _parameter_values(tape, parameters)
    instrs = tape.instrs

    def column(name):
        return np.array([getattr(i, name) for i in instrs], dtype=np.int32)

    return dict(
        op=np.array([int(i.op) for i in instrs], dtype=np.uint8),
        a=column("a"),
        b=column("b"),
        c=column("c"),
        out=column("out"),
        imm=np.array([0.0 if i.imm is None else float(i.imm) for i in instrs], dtype=np.float64),
        n_registers=int(tape.n_regs),
        output=int(tape.out_reg),
        inputs=np.array([tape.x_reg, tape.y_reg, tape.z_reg], dtype=np.int32),
        extra_registers=np.array([reg for _, reg in tape.extra_inputs], dtype=np.int32),
        extra_values=np.array(values, dtype=np.float64),
    )


def analytic_level_set(shape, parameters=None):
    """A ShapeForest shape, expression or tape as a
    :class:`cutcells.AnalyticLevelSet`, evaluated by CutCells' interpreter of
    ShapeForest tapes (values, gradients and Taylor-model bounds)."""
    return analytic_level_set_from_tape(**tape_arrays(shape, parameters))


def write_tape(shape, path, parameters=None):
    """Write the tape of a shape in the text format that the benchmark
    ``quadrays_study --tape`` reads: "shapeforest-tape 1", then
    "n_regs out x y z n", then one line "op a b c out imm" per instruction.
    Parameters become constant instructions at the start."""
    arrays = tape_arrays(shape, parameters)
    lines = [
        (0, -1, -1, -1, int(reg), repr(float(value)))
        for reg, value in zip(arrays["extra_registers"], arrays["extra_values"])
    ]
    for k in range(len(arrays["op"])):
        lines.append(
            (
                int(arrays["op"][k]),
                int(arrays["a"][k]),
                int(arrays["b"][k]),
                int(arrays["c"][k]),
                int(arrays["out"][k]),
                repr(float(arrays["imm"][k])),
            )
        )
    x, y, z = (int(r) for r in arrays["inputs"])
    with open(path, "w") as out:
        out.write("shapeforest-tape 1\n")
        out.write(f"{arrays['n_registers']} {arrays['output']} {x} {y} {z} {len(lines)}\n")
        for line in lines:
            out.write(" ".join(str(v) for v in line) + "\n")
