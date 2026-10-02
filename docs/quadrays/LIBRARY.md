# quadrays in the library

The prototype's certify engine (`README.md`, `RESULTS.md` here) is now the
`quadrays` backend: `cpp/src/quadrays/`, namespace `cutcells::quadrays`. It
reproduces the prototype's numbers to every printed digit and runs 2.4 to 3.3
times faster per cut cell. Phase 2 added analytic level sets from any geometry
library (below).

## Files

| File | Contents |
| --- | --- |
| `box_bernstein.h/.cpp` | tensor Bernstein forms on boxes: a cell's level set converted exactly (a rational matrix per degree for simplices, computed once), restriction to affine images, derivatives, margins, roots on lines |
| `clipped_box.h/.cpp` | a cell as the unit box clipped by half-spaces, from its type and vertices; maps to physical and reference coordinates |
| `engine.h/.cpp` | dimension reduction with certified height directions, bisection, the diagonal frame; `Options`, `Stats`, `part_of` |
| `rules.h/.cpp` | rules for one cell and one selection term into `quadrature::QuadratureRules` |
| `leaves.h/.cpp` | the decomposition as VTK Lagrange cells; `io::write_lagrange_vtk` writes them (VTK 9.1 node order) |
| `analytic.h/.cpp` | the interface of analytic level sets: a C struct of function pointers (value, gradient, bounds over a box, optional Taylor models over a parallelepiped) and a context; `analytic_level_set(functor)` fills it from an algoim-style functor |
| `taylor.h` | first-order Taylor models with rigorous remainders (sqrt, log, exp, sin, cos, division), `Dual<V, N>` for derivatives, branch helpers for abs, min and max |
| `source.h/.cpp` | the level set of a cell as the engine reads it: a Bernstein form or an analytic level set; values, gradients, Taylor bounds on affine images, roots on lines |
| `adapters/shapeforest_tape.h` | ShapeForest tapes (register code) run on every scalar type of `taylor.h`; no build dependency |

The front end of phase 3, `cutcells.part`, uses quadrays without AdaptCell
(`FRONTEND.md`). From Python, `part.quadrature(order, mode, backend="quadrays")` uses it, with
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

## Analytic level sets (phase 2)

An analytic level set reaches quadrays through `AnalyticLevelSet` (`analytic.h`),
a plain struct: `value(x)`, `gradient(x)`, `box_bounds(lo, hi)` (intervals of the
value and of the gradient over a box) and, optionally, `taylor_bounds` (first-order
Taylor models of the value and of the derivatives along the axes of a
parallelepiped), each with a `void* context`. Any library can fill it, from C++ or
through a Python capsule, without CutCells' templates. `analytic_level_set(f)`
fills it from a functor `template <typename T> T operator()(const std::array<T, 3>&)`
that it evaluates with `double`, `Dual<double, 3>` and
`Dual<Taylor<double, m>, m>`. The engine takes its margins from these models instead
of Bernstein coefficients; the maps from a box at any level to physical space are
affine, so the models represent them exactly. Roots on lines are isolated with
one-dimensional models.

From Python:

- `cutcells.AnalyticLevelSet(value, gradient, box_bounds, taylor_bounds=None)` from
  Python callables (slow, for experiments);
  `AnalyticLevelSet.from_capsule(capsule)` from a capsule named
  `cutcells.AnalyticLevelSet` that points to the C struct (fast, for compiled
  libraries); `level_set.capsule` gives one;
- `cutcells.analytic_sphere(centre, radius, signed_distance=True)`, a compiled
  functor;
- `cutcells.shapeforest.analytic_level_set(shape)` for ShapeForest shapes,
  expressions and tapes (`analytic_level_set_from_tape` takes the arrays);
  `cutcells.shapeforest.write_tape` writes the text format the benchmark reads;
- `cutcells.cut(mesh, level_set, degree=2)` and `create_level_set(mesh, level_set,
  degree)`: tetrahedra and hexahedra are classified as inside, outside or cut by
  the level set's own bounds (`cell_sign` in `source.h`: Taylor models over the
  cell, its corners, bisection), and `backend="quadrays"` integrates it. The
  interpolant of degree `degree` only feeds today's AdaptCells, which the straight
  and algoim backends and the straight visualisation use;
- `quadrays_cell_rules` and `quadrays_cell_leaves` take an `AnalyticLevelSet` in
  place of the degree and the Bernstein coefficients.

### Results (2026-10-02)

The sphere of radius 0.7 off the grid as the signed distance |x - c| - r, a C++
functor through the interface; q = 5, margin 0.25; per-cell L1 / worst cell;
points per cut cell (volume / interface):

| Mesh | Level set | Interface | Volume | Points | Bisections |
| --- | --- | --- | --- | --- | --- |
| hexes n = 32 | distance | 1.3e-9 / 4.2e-6 | 1.1e-12 / 2.9e-8 | 226.1 / 32.3 | 28 |
| hexes n = 32 | quadratic, Bernstein | 1.3e-9 / 4.2e-6 | 1.1e-12 / 2.9e-8 | 225.8 / 32.3 | 25 |
| hexes n = 32 | distance, algoim 2015 | 6.9e-8 / 2.4e-4 | 8.9e-12 / 2.9e-8 | 224.4 / 32.0 | |
| tets n = 32 | distance | 1.5e-9 / 7.7e-6 | 8.5e-12 / 5.3e-7 | 376.3 / 40.2 | 962 |
| tets n = 32 | quadratic, Bernstein | 2.2e-9 / 7.7e-6 | 1.2e-11 / 5.3e-7 | 366.3 / 38.9 | 604 |
| tets n = 16 | distance | 6.7e-9 / 2.4e-6 | 9.8e-11 / 1.1e-6 | 1,259.3 / 127.2 | 7,178 |
| tets n = 16 | quadratic, Bernstein | 1.4e-8 / 2.7e-6 | 2.2e-10 / 2.8e-6 | 599.1 / 62.4 | 1,974 |

The distance reproduces the polynomial sphere's errors on hexes and is better on
tets, where its looser bounds bisect more. The prototype's tape numbers (v1.3) are
reproduced to every printed digit. ShapeForest's sphere as a tape gives exactly the
functor's numbers; |x - c|^2 - r^2 as a functor gives exactly the Bernstein path's.

Microseconds per cut cell at n = 32, q = 5, three interleaved runs, smallest kept
(load average 3 to 4); algoim's 2015 engine evaluates the same functor or tape,
with its bounds through `taylor.h` (algoim's sqrt leaves out part of the
remainder):

| | hex volume | hex interface | tet volume | tet interface |
| --- | --- | --- | --- | --- |
| quadrays, distance functor | 36.6 | 23.0 | 61.2 | 38.0 |
| algoim 2015, same functor | 25.5 | 17.7 | | |
| quadrays, ShapeForest tape | 63.8 | 38.5 | 105.7 | 60.8 |
| algoim 2015, same tape | 52.0 | 44.2 | | |
| quadrays, Bernstein | 60.5 | 44.0 | 111.0 | 83.5 |

On hexes quadrays takes 1.3 to 1.4 times algoim's time on the functor and 0.9 to
1.2 times on the tape. The robustness cases also run through the interface
(`test_robustness`): no failure, non-finite value, negative weight, point outside
or wrong side; planes as exact as with Bernstein coefficients (1e-15 per cell, 8e-12
for the plane 1e-12 from faces on tets); scaled level sets give the unscaled totals.

A ball of radius 0.26 about the centre of a cell of side 0.5 reaches 0.01 into
its six neighbours: every vertex lies outside it, so a P1 interpolant sees no ball
at all, and the caps lie away from all vertices. `cut()` finds the 7 cut hexahedra
(24 tetrahedra) for any degree, and quadrays integrates the ball to 5e-11 (volume)
and 6e-10 (area) on hexahedra, 1e-10 and 2e-9 on tetrahedra; the straight backend
on the same cut sees no ball with degree 1 and 37% (hexahedra) or 54%
(tetrahedra) of it with degree 2. On
n = 32 hexahedra the classification finds exactly the 2,354 cells the sphere of
the study cuts (the degree-2 interpolant flags 2,360), in the same time.

Two things made that cell work:

- **Bounds of the value alone.** Around the centre of a distance function the
  argument of the square root reaches 0. There is no first-order model there, but
  the values still lie in [0, sqrt(max)]; only the derivative (a quotient by the
  square root) loses its bound. `taylor_bounds` and `box_bounds` return 2 in that
  case, and the engine drops boxes whose value has one sign instead of bisecting
  them.
- **Deeper limits.** Taylor models of the distance's derivative are quotients by
  |x - c|, loose until a box is small against the radius; here that takes about
  10 bisections per level. The defaults are now `max_depth` 12 and
  `max_bisections` 1024 (the prototype's 8 and 256, tuned on Bernstein bounds,
  left the area 1.4e-2 off). Runs that never reached the old limits are
  unchanged, the n = 32 studies above included. Hard cells improve: touching
  spheres on n = 4 tetrahedra go from 7.1e-2 to 6.8e-4 per-cell L1 on the area
  (Bernstein; 2.4e-2 to 9.6e-4 analytic), the torus area from 6.3e-4 to 1.8e-5
  (analytic 7.6e-3 to 1.6e-5); `test_robustness` takes 96 s instead of 62.

### Found on the way

- `append_rules` reserved exactly the room it needed on every call, which made
  `HOMeshPart.quadrature(backend="quadrays")` quadratic in the number of cut cells
  (23 s instead of 0.5 s for 640 tetrahedra). The C++ tests and benchmarks use one
  rule set per cell and were not affected.
- A line lying in the zero set (a plane in mesh faces) has no certain sign and no
  certain slope anywhere, so isolating its roots by bisection took 2^40 steps.
  Intervals where the models bound |phi| below the engine's zero tolerance now end
  the bisection, and a line visits at most 4,096 intervals.
- The vertex a tetrahedron's box starts from matters for Taylor bounds: with a
  middle vertex of a Kuhn tetrahedron's path as the origin, the distance sphere at
  n = 8 needs 3,300 to 3,800 bisections, with an end vertex 7,300 to 9,300 (the box
  is less skewed). Choosing the vertex is one way to tighter bounds on tets
  (phase 6).

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
build-quadrays/benchmarks/quadrays/quadrays_study --mesh hex --n 32 --q 5 --analytic distance --gen quadrays,quadgen
PYTHONPATH=build-quadrays/python $E/bin/python -c "import shapeforest as sf; from cutcells import shapeforest as csf; csf.write_tape(sf.sphere(0.7).translate(0.0123, -0.0371, 0.0217), 'sphere.tape')"
build-quadrays/benchmarks/quadrays/quadrays_study --mesh hex --n 32 --q 5 --tape sphere.tape --gen quadrays,quadgen
```
