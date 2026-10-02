# quadrays in the library

The prototype's certify engine (`README.md`, `RESULTS.md` here) is now the
`quadrays` backend: `cpp/src/quadrays/`, namespace `cutcells::quadrays`. It
reproduces the prototype's numbers to every printed digit and runs 2.4 to 3.3
times faster per cut cell.

## Files

| File | Contents |
| --- | --- |
| `box_bernstein.h/.cpp` | tensor Bernstein forms on boxes: a cell's level set converted exactly (a rational matrix per degree for simplices, computed once), restriction to affine images, derivatives, margins, roots on lines |
| `clipped_box.h/.cpp` | a cell as the unit box clipped by half-spaces, from its type and vertices; maps to physical and reference coordinates |
| `engine.h/.cpp` | dimension reduction with certified height directions, bisection, the diagonal frame; `Options`, `Stats`, `part_of` |
| `rules.h/.cpp` | rules for one cell and one selection term into `quadrature::QuadratureRules` |
| `leaves.h/.cpp` | the decomposition as VTK Lagrange cells; `io::write_lagrange_vtk` writes them (VTK 9.1 node order) |

From Python, `part.quadrature(order, mode, backend="quadrays")` uses it, with
`order` Gauss-Legendre points per segment as for `backend="algoim"`; uncut cells
of volume parts get the straight rules of degree `order`. `quadrays_quadrature`
and `quadrays_leaves` take a `QuadraysOptions`; `quadrays_cell_rules` and
`quadrays_cell_leaves` work on one cell.

## Differences from the prototype

- **No interpolation, no LAPACK.** Restrictions of the level set to faces, clip
  planes, frames and lines are computed by de Casteljau's algorithm with an
  affine argument; along axes this is one blossom per fiber.
- **Sign test.** A Bernstein form counts as possibly vanishing only if its
  coefficients take both signs. With exact restrictions, a level set touching a
  face at a box corner has an exact zero coefficient there; the prototype's
  interpolation turned it into noise of either sign. The engine then no longer
  bisects towards tangency points: 0 instead of 199 bisections for the sphere
  tangent to grid planes on hexahedra, at the same accuracy.
- **Roots on lines.** A root exactly on a split point of the root isolation is
  found once, and a piece whose end is a root is subdivided further instead of
  dropping the root inside it.
- **Diagnostics.** `Options::diagnose` and `Stats::causes` are kept; the
  `CERTIFY_TRACE` printout is not.

## Results (2026-10-02)

Sphere of radius 0.7 off the grid, n = 32, q = 5, margin 0.25; per-cell L1 / worst
cell; points per cut cell (volume / interface):

| Mesh | Interface | Volume | Points | Bisections |
| --- | --- | --- | --- | --- |
| tets | 2.2e-9 / 7.7e-6 | 1.2e-11 / 5.3e-7 | 366.3 / 38.9 | 604 |
| hexes | 1.3e-9 / 4.2e-6 | 1.1e-12 / 2.9e-8 | 225.8 / 32.3 | 25 |

Microseconds per cut cell, prototype and library run one after another three
times on a loaded machine (load average 38 to 105), smallest of the three:

| | hex volume | hex interface | tet volume | tet interface |
| --- | --- | --- | --- | --- |
| prototype | 490 | 349 | 814 | 582 |
| quadrays | 149 | 111 | 269 | 243 |

Planes are integrated exactly to 3e-15 per cell. In the robustness report
(batch 1 on n = 16, batch 2 on n = 8, q = 3) no run has a failure, a non-finite
value, a negative weight, or a point outside its cell or on the wrong side; the
prototype had wrong-side points for the cone with its apex on a vertex and for
the double root. Of 96 rows, 78 are identical to the prototype's; the others are
the tangent placements, the cones (fewer points), and interfaces nobody owns
(in mesh faces, the double root). VTK 9.6 reads the leaves in the intended
order (Lagrange hexahedra and quadrilaterals of degrees 1 to 4).

## Build and run

```bash
E=$HOME/miniforge3/envs/fenicsx0.11
env CONDA_PREFIX=$E CMAKE_PREFIX_PATH=$E CC=$E/bin/clang CXX=$E/bin/clang++ $E/bin/cmake -S cpp -B build-quadrays -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DCMAKE_OSX_DEPLOYMENT_TARGET=13.4
$E/bin/cmake --build build-quadrays -j 4 && $E/bin/ctest --test-dir build-quadrays
$E/bin/cmake --install build-quadrays --prefix build-quadrays/install
env CONDA_PREFIX=$E CC=$E/bin/clang CXX=$E/bin/clang++ $E/bin/cmake -S benchmarks -B build-quadrays/benchmarks -DCMAKE_PREFIX_PATH="$PWD/build-quadrays/install;$E" -DCMAKE_BUILD_TYPE=Release -DCUTCELLS_WITH_ALGOIM=ON -DCMAKE_OSX_DEPLOYMENT_TARGET=13.4
$E/bin/cmake --build build-quadrays/benchmarks
build-quadrays/benchmarks/quadrays/quadrays_study --mesh tet --n 32 --q 5 --gen quadrays,alpha-split
build-quadrays/benchmarks/quadrays/quadrays_robustness_report --list
```
