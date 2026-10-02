# Quadrature prototype

Compares ways of generating quadrature rules on cut cells before the new
engine is designed. It is standalone: it builds against the vendored algoim
and CutCells' selection-expression parser, nothing else.

Background: the hand-off note "CutCells quadrature reimplementation" and its
tab "3D surface stall: findings", and the per-cell study in
`CutFEMx/docs/plans/clipped_box_quadrature/surface_stall/`.

## Keeping today's query interface

Users select parts of the domain with expressions such as
`result["phi1 < 0 and phi2 = 0"]` and ask the part for quadrature. The
prototype keeps that shape so a later engine can become another `backend` of
`HOMeshPart.quadrature(...)`:

- Parts are parsed and compiled with CutCells' own `parse_selection_expr` and
  `compile_selection_expr`. Generators take the compiled `SelectionTerm`, with
  its zero, negative and positive bitmasks, not a hard-coded "inside".
- Rules come back as points in the parent cell's reference coordinates with
  physical weights, as in `quadrature::QuadratureRules`.
- For now a term may only constrain the first level set (`part_kind`).
  Several level sets are the next step.

## Layout

| File | Contents |
| --- | --- |
| `src/clipped_box.h/.cpp` | a cell as the unit box clipped by half-spaces; maps to physical and parent-reference coordinates; sub-boxes; the exact clipped polytope |
| `src/generators.h/.cpp` | the algoim-based generators; also compiles `certified.inl`, because algoim's headers can only be included in one translation unit |
| `src/certified.h/.inl` | the certify-and-bisect engine |
| `src/vtk_output.h` | point-cloud `.vtu` writer |
| `src/leaf_mesh.h` | leaf meshes (CSR), VTK Lagrange node order, `.vtu` writer |
| `src/test_leaf_ordering.cpp` | writes Lagrange cells for checking the node order against VTK |
| `tools/check_leaves.py` | checks leaf files with VTK (node order, sizes) and renders them; needs vtk and pyvista, e.g. the `vtk-env` environment |
| `src/exact_reference.h` | exact sphere area and ball volume inside a convex polytope (about 1e-14); convex polytope cut by a plane |
| `src/study.cpp` | comparison driver |
| `src/robustness.cpp` | robustness driver: degenerate placements, scaling, several components, singular points; every cell, checked against exact values (RESULTS.md, v1.2) |
| `src/test_mesh.h` | the test meshes (hexes and Kuhn tets of [-1, 1]^3) |
| `src/tape.h` | ShapeForest tapes: reader and `evaluate<T>`, one templated interpreter for values, gradients (`Dual`) and bounds |
| `src/tape_algoim.h` | branch helpers and a rigorous `tsqrt` for algoim's Taylor-model intervals |
| `tools/export_tape.py` | writes a ShapeForest shape as a tape for `quadrature_study --tape` (run with an environment that has ShapeForest, e.g. `fenicsx0.11`) |
| `src/test_exact_reference.cpp` | checks the exact reference |
| `CMakeLists.txt` | also writes a patched copy of algoim's `quadrature_multipoly.hpp` to the build tree; `third_party/` is untouched |

The algoim patch has three parts: a fix for the two-polynomial constructor
with user masks (it does not compile upstream: `maskEmpty` needs `detail::`),
a constructor taking any number of polynomials, and `force_k[N]`, which
overrides the height direction chosen at dimension N.

## Generators

| Name | What it does |
| --- | --- |
| `algoim-auto` | algoim's 2022 multi-polynomial engine with AutoMixed, as CutCells' `"algoim"` backend does today; tets are boxes clipped by their facet plane |
| `algoim-gl` | the same with Gauss-Legendre everywhere |
| `gl-cellmask` | plus algoim's masks restricted to sub-cells that meet the clipped cell |
| `alpha` | plus height axes chosen by Cui et al.'s angle indicator alpha on the clipped region |
| `split` | `gl-cellmask` plus splitting the box along its longest axis while the best alpha >= 0.99, up to 3 levels; algoim still chooses the axes |
| `alpha-split` | `alpha` plus the same splitting |
| `quadgen` | algoim's 2015 engine (interval arithmetic and bisection), CutCells' `"algoim_general"`; hexahedra and spheres only |
| `certify`, `certify:<margin>` | the new certify-and-bisect engine (`src/certified.inl`), margin 0.25 by default; clip planes handled directly, no resultants; bounds that are never active are pruned; at level 2 a diagonal frame is tried before bisecting |

## Metrics

For each cut cell the rule is compared with the exact value of the selected
part. Reported per run:

- per-cell L1: sum over cut cells of |cell error|, relative to the ball volume
  (volume parts) or the sphere area (interface);
- worst cell: largest |cell error| relative to that cell's exact value;
- global: error of the total, which hides per-cell errors that cancel;
- points and microseconds per cut cell. Timings come from a shared machine, so
  compare only runs made together.

## Build and run

Use the `fenicsx0.11` environment, as for the rest of the FEniCSx-pr workspace.

```bash
E=/Users/sclaus/miniforge3/envs/fenicsx0.11
env CONDA_PREFIX=$E CMAKE_PREFIX_PATH=$E CC=$E/bin/clang CXX=$E/bin/clang++ $E/bin/cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
$E/bin/cmake --build build -j 4
```

```bash
./build/test_exact_reference
./build/quadrature_study --mesh tet --n 16,32 --q 3,5 --gen algoim-auto,alpha-split --part "phi < 0" --part "phi = 0" --csv tet.csv
```

Options: `--mesh hex|tet`, `--n`, `--q`, `--centre x,y,z` (default off the
grid's symmetry, 0.0123,-0.0371,0.0217), `--radius`, `--gen`, `--part`
(repeatable), `--csv`, `--vtk <prefix>` (point clouds), `--leaves <prefix>` and
`--leaf-degree p` (leaf cells of the certify engine as Lagrange cells; see
VISUALIZATION.md), `--plane` (planar level set; every rule must be exact).

`--tape <file>` reads the level set as a ShapeForest tape describing the same
sphere (`python tools/export_tape.py sphere-sdf 0.0123,-0.0371,0.0217 0.7
build/sphere_sdf.tape`); certify and quadgen then evaluate it directly.

For the certify engine: `--diagnose` prints why certification failed (by level,
blocking function and cause) and the worst cell; `--only <cell>` restricts a run
to one cell index; `--masks M` and `--no-diagonal` set
`CertifyOptions::mask_subdivisions` and `diagonal_frames`. With the environment
variable `CERTIFY_TRACE` set, every box is printed to stderr with its margins.

Configure with the environment variables above: without `CONDA_PREFIX`, CMake
picks Apple's Accelerate, which lacks the LAPACKE symbols algoim needs.

## Results and visualisation

See `RESULTS.md` and `VISUALIZATION.md`.
