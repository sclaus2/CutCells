# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Writes a ShapeForest level set as a text tape for the prototype (quadrature_study --tape).

    python tools/export_tape.py sphere-sdf 0.0123,-0.0371,0.0217 0.7 build/sphere_sdf.tape
    python tools/export_tape.py sphere-quadratic 0.0123,-0.0371,0.0217 0.7 build/sphere_quadratic.tape

Shapes: sphere-sdf (ShapeForest's sphere, |x - c| - r) and sphere-quadratic
(|x - c|^2 - r^2, the polynomial the other generators use).

Format: a header line "shapeforest-tape 1", then "n_regs out_reg x_reg y_reg z_reg
n_instr", then one line per instruction "op a b c out imm" with ShapeForest's
opcodes (shapeforest/_cpp/opcodes.hpp).
"""

import sys

import shapeforest as sf
from shapeforest.ir.expr_to_tape import expr_to_tape


def shape_tape(name, centre, radius):
    cx, cy, cz = centre
    if name == "sphere-sdf":
        return sf.sphere(radius).translate(cx, cy, cz).to_tape()
    if name == "sphere-quadratic":
        X, Y, Z = sf.X, sf.Y, sf.Z
        expr = (X - cx) * (X - cx) + (Y - cy) * (Y - cy) + (Z - cz) * (Z - cz) - radius * radius
        return expr_to_tape(expr)
    raise SystemExit(f"unknown shape {name}")


def write_tape(tape, path):
    if tape.extra_inputs:
        raise SystemExit("tapes with parameters are not supported")
    with open(path, "w") as out:
        out.write("shapeforest-tape 1\n")
        out.write(f"{tape.n_regs} {tape.out_reg} {tape.x_reg} {tape.y_reg} {tape.z_reg} {len(tape.instrs)}\n")
        for ins in tape.instrs:
            imm = 0.0 if ins.imm is None else ins.imm
            out.write(f"{int(ins.op)} {ins.a} {ins.b} {ins.c} {ins.out} {imm!r}\n")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    centre = [float(v) for v in sys.argv[2].split(",")]
    write_tape(shape_tape(sys.argv[1], centre, float(sys.argv[3])), sys.argv[4])
