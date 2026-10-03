# quadrays in the library

The prototype's certify engine (`README.md`, `RESULTS.md` here) is now the
`quadrays` backend: `cpp/src/quadrays/`, namespace `cutcells::quadrays`. It
reproduces the prototype's numbers to every printed digit and runs 2.4 to 3.3
times faster per cut cell. Phase 2 added analytic level sets from any geometry
library, phase 6 2D cells, prisms, pyramids, several level sets per cell and
two roots per line (below).

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
| `../compression/compress.h/.cpp` | positive rules reduced to at most as many points as a polynomial space has moments (below) |

The front end (`cutcells.cut`, `cutcells.part`; `FRONTEND.md`) integrates with
quadrays when a part's quadrature asks for `backend="quadrays"`, with `order`
Gauss-Legendre points per segment of each height line; whole cells and faces in
a zero set get the reference rules exact for degree 2 order - 1.
`quadrays_quadrature` and `quadrays_leaves` take a part and a `QuadraysOptions`;
`quadrays_cell_rules` and `quadrays_cell_leaves` work on one cell.

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
  degree)`: cells (tetrahedra and hexahedra, since phase 6 all quadrays takes)
  are classified as inside, outside or cut by
  the level set's own bounds (`cell_sign` in `source.h`: Taylor models over the
  cell, its corners, bisection), and `backend="quadrays"` integrates it. Until
  phase 5 the interpolant of degree `degree` fed the AdaptCells of the straight
  and algoim backends; since then the lookup tables also read the analytic level
  set itself, at their template's vertices;
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
(24 tetrahedra; 18 since phase 3 bounds tetrahedra themselves) for any degree, and quadrays integrates the ball to 5e-11 (volume)
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

## Compressed rules (2026-10-02)

quadrays certifies boxes and puts q^D points in each, so a cut cell ends with
hundreds or thousands of points, all of which the assembler visits.
`cpp/src/compression/` (namespace `cutcells::compression`) shrinks any positive
rule to at most as many points as a polynomial space has moments, with
positive weights that integrate that space exactly as the original rule did
(Caratheodory-Tchakaloff). For Q_k elements on affine hexahedra the stiffness
and mass integrands lie in Q_2k, so `compress_rules(rules, 2k)` changes
nothing in the assembled matrices up to rounding; on affine tetrahedra with
P_k elements use `space="total"`. Outside the space, a compressed rule is less
accurate than the original rule.

The points are recombined in groups (Tchernychova and Lyons): 2M groups of
consecutive points, each replaced by its weighted mean moment vector, are
reduced to M by a Caratheodory step, which halves the points; the last at most
2M points are reduced directly. A Caratheodory step takes the null space from
a pivoted Householder QR, P [-R11^-1 R12; I], and moves the weights along one
null vector at a time until a weight reaches 0. The basis is Legendre on the
points' bounding box; directions the points do not resolve (a piece flat in
one coordinate) drop out of the rank, and such pieces keep fewer points. A
least-squares solve on the chosen points polishes the weights if they stay
positive. Rules with a negative weight are copied; rules run in parallel with
OpenMP.

Turbine blade, h = 3 cm hexahedra, the DOLFINx Q4 interpolant of the
ShapeForest level set, 1,774 cut cells, Q_4 (125 moments). Errors are per cell
against quadrays at q = 8, relative to the cell volume; times on an i7-7920HQ
under load:

| Rule | Points per cut cell | Volume error | Time |
| --- | --- | --- | --- |
| Algoim, q = 4 (CutFEMx) | 323 | 4.4e-3 | 314 s |
| quadrays, q = 4 | 2,149 | 1.6e-5 | 3.3 s |
| quadrays, q = 4, compressed | 100 | 1.6e-5 | + 5.1 s (14.8 s on one thread) |
| quadrays, q = 4, margin 0.02 | 827 | 3.3e-4 | 1.0 s |
| quadrays, q = 4, margin 0.02, compressed | 97 | 3.3e-4 | + 3.0 s (13.1 s on one thread) |

A random Q_4 polynomial integrates to the same value before and after
compression (largest moment residual 1e-11 of a rule's weight).

On the cases of phase 6 (q = 4; n = 8, quadrilaterals and triangles n = 16),
points per cut cell before and after compression to degree 2 (P1 and Q1
elements) and degree 4 (P2 and Q2), total degree on simplices, tensor on the
other cells:

| Case | Volume | Degree 2 | Degree 4 | Interface | Degree 2 | Degree 4 |
| --- | --- | --- | --- | --- | --- | --- |
| hexahedra, sphere | 485 | 27 | 118 | 71 | 25 | 63 |
| tetrahedra, sphere | 1,102 | 10 | 35 | 141 | 9 | 24 |
| tetrahedra, P2 sphere | 1,008 | 10 | 35 | 128 | 9 | 24 |
| prisms, sphere | 390 | 27 | 112 | 55 | 23 | 47 |
| pyramids, sphere | 423 | 27 | 115 | 63 | 24 | 48 |
| hexahedra, lens of two balls | 825 | 27 | 105 | 137 | 24 | 61 |
| tetrahedra, lens | 4,772 | 10 | 35 | 383 | 9 | 24 |
| quadrilaterals / triangles, circle | 20 / 24 | 9 / 6 | 18 / 15 | 4 | 4 | 4 |

Moments stay exact to 6e-11 of a rule's weight. A function outside the space,
exp(x + y + z) cos 2x in reference coordinates, moves by up to 7e-4 of a rule's
weight at degree 4 and 5e-2 at degree 2. Most interface rules on hexahedra have
fewer points than Q_4 has moments and are kept. Pyramids' elements are
rational, so no polynomial space makes their compression exact for them.
Compressing to degree 4 took 0.01 to 0.7 s per mesh, about as long as quadrays
took to make the rules.

## Phase 6: 2D cells, prisms and pyramids, several level sets, two roots

### 2D cells

Triangles and quadrilaterals are boxes in the plane u2 = 0 (`ClippedBox::tdim`
= 2): the engine starts at level 2; a triangle is the unit square clipped by
u0 + u1 <= 1. Lines are exact on every cell (largest error 7.8e-16 per cell,
relative to h^2 or h). The circle of radius 0.7 off the grid's symmetry, n = 8,
per-cell L1 of the disk and the circle, relative to their totals:

| | q = 3 | q = 5 | q = 8 |
| --- | --- | --- | --- |
| quadrilaterals | 9.0e-7 / 5.1e-6 | 4.9e-10 / 5.2e-9 | 1.5e-14 / 2.5e-13 |
| triangles | 5.0e-7 / 1.1e-5 | 7.7e-10 / 3.2e-8 | 9.6e-14 / 6.7e-12 |

(P2 coefficients; the analytic circle agrees to two digits.) The robustness
cases in 2D (circles through vertices, tangent to edges, caps of 1e-6 and
1e-12, scaled by 1e+-150 and 1e+-200; lines in edges and through vertices; two
circles touching; an annulus 1e-3 wide; two lines crossing in a vertex and in a
cell) pass every check on n = 16 with q = 3, with a per-cell L1 of at most 8e-6.

### Prisms and pyramids

A prism is the unit box clipped by u0 + u1 <= 1, a pyramid the box clipped by
u0 + u2 <= 1 and u1 + u2 <= 1. Pk level sets need Bernstein forms on them,
which `bernstein.h` now has: on a prism the triangle's basis times the
interval's (the space of Basix's prism elements); on a pyramid the functions
B^(n-k)_i(s) B^(n-k)_j(t) B^n_k(z) of s = x / (1 - z), t = y / (1 - z) (Chan
and Warburton), which span the rational space of Basix's pyramid elements (P1
holds xy / (1 - z), not xy). A prism's form converts to the box layer by layer
like a triangle's; for a pyramid the engine reads (1 - z)^n phi, a polynomial of
degree (n, n, 2n) with phi's sign and zero set inside the pyramid and, on the
zero set, its normal: the apex, where the factor vanishes, is the only zero it
adds. That keeps the map affine and planes exact; reading the form on the cube
of (s, t, z) instead, with a collapsed map, made planes curved surfaces (2.8e-5
per cell at q = 3).

Planes are exact as P1, P2 and analytic level sets (largest error 3.2e-15). The
sphere of radius 0.7, n = 8, per-cell L1 of the ball and the sphere:

| | q = 3 | q = 5 | q = 8 |
| --- | --- | --- | --- |
| prisms, P2 and analytic | 5.9e-6 / 4.4e-5 | 1.6e-8 / 3.4e-7 | 1.3e-11 / 7.5e-10 |
| pyramids, analytic | 8.9e-7 / 7.9e-6 | 1.3e-9 / 3.7e-8 | 6.8e-13 / 4.2e-11 |
| pyramids, P2 | 2.2e-9 / 6.9e-8 | 5.2e-13 / 1.8e-10 | 1.1e-14 / 1.1e-13 |

A pyramid's P2 form costs: where a zero set passes near the apex, it nears the
zero face z = 1 of (1 - z)^n phi there, and the engine bisects towards it (346
bisections in the worst cell; 3.5 s for the 3,072 pyramids at n = 8, 0.6 s with
the analytic sphere). The robustness cases of the 3D tests pass on prisms and
pyramids, with Bernstein forms and analytic level sets.

### Several level sets in one cell

The terms of a selection are requirement masks over the cell's level sets;
points are filtered by them, and an interface part integrates the zero set of
one level set. Where two level sets cut a certified box, their zero sets meet
on a ridge: the base splits there, where one level set vanishes on the other's
zero set. That function, s(y) = on(y, r(y)) with r the root of under along the
height line, is a nested height function (Beck and Kummer): its root is
followed half the height range beyond the box, so that s stays smooth where the
root leaves it; its bounds come from the slab the root sweeps over the box,
about the root where under's derivative changes sign over the range. The
non-interface level sets of an interface part need no margins (passive), but
their restrictions to the bounds stay, so that ridges through faces are found.
Where several level sets meet, the base splits into regions with one active
lower and one upper bound; where no axis suits both zero sets, a rotated frame
between their normals does (always in 2D, from depth 6 in 3D); where three zero
sets meet, the crossings of two surface functions give the corners. An
analytic level set that is linear on the cell (its Hessian bounds vanish) is
read as its form of degree 1.

Ball and half-space per cell, n = 8, q = 5, per-cell L1 relative to the totals:

| | cap | its sphere | disk | union |
| --- | --- | --- | --- | --- |
| hexahedra | 2.3e-9 | 4.0e-8 | 5.6e-12 | 1.5e-9 |
| tetrahedra, P2 and P1 | 7.1e-10 | 3.4e-8 | 4.4e-14 | 1.7e-10 |
| tetrahedra, analytic | 9.1e-10 | 4.3e-8 | 3.9e-14 | 2.1e-10 |

A disk and a half-plane in 2D reach 1e-15 at q = 8. Exact totals, n = 8, q = 5:
the lens of two balls within 1.0e-7, the Steinmetz bicylinder within 1.5e-6
(its cylinders touch where their intersection curves cross). The robustness
cases (planes through vertices and 1e-3 off faces, a sphere with its tangent
plane, two balls touching, a plane in grid faces with a sphere through
vertices) pass every check.

### Two roots per line

Bisection rarely separates two sheets of one level set, so the roots decide:
a direction qualifies where the derivative along it is monotone
(second-derivative margin) and the extreme points on a grid of lines, with the
planes through them, show that the sheets never merge in the box; each line
then splits at its extreme point. The shell 1e-3 wide on the octant of n = 8,
q = 3, per-cell L1 of volume and area:

| | two roots | one root per line |
| --- | --- | --- |
| hexahedra, P4 | 1.0e-5 / 7.4e-4 | 1.7e-2 / 3.4e-2 |
| tetrahedra, P4 | 4.4e-6 / 2.7e-4 | 1.1e-5 / 1.3e-4 |
| hexahedra, analytic | 1.0e-2 / 2.5e-2 | 1.6e-2 / 3.0e-2 |
| tetrahedra, analytic | 4.6e-6 / 2.7e-4 | 1.1e-5 / 1.7e-4 |

The shell 1e-6 wide passes the checks with 5.5e-2; analytic shells stay weak
(Taylor models of the product separate the sheets late).

### Taylor sub-boxes

The margins of an analytic level set come from Taylor models over the 2^D
sub-boxes on which it may vanish (`taylor_subdivisions`); sign and size still
come from the whole box. On the tetrahedra at n = 16, q = 5: 1275 instead of
7144 bisections, 123 instead of 266 microseconds per cut cell, L1 4.0e-10
instead of 9.8e-11: the certified boxes are larger, their accuracy then set by
q and the margin as for Bernstein forms.

### Found on the way

- Analytic roots exactly on a split point of the root isolation were lost: two
  lines crossing in a vertex lost a box's share of their length.
- A bisection exactly through an interface leaves it on the face the halves
  share, where neither has a root: the engine now splits at 5/8 there.
- Surface functions with Taylor sub-boxes needed the root's slab over the whole
  height range: a ball and a plane on tetrahedra took 53 s instead of 0.3 s, a
  distance ball with a plane 7 minutes instead of 4 s on n = 4. A slab about the
  root now holds it where the range does not.
- The fold surfaces of the first two-root attempt are gone.
- Level-set mesh data: pyramids take degree 2, prisms get their interior nodes
  from degree 3 in the reference points too, and two of the pyramid's triangle
  faces were wrong (VTK's base order with Basix numbering; no dofs on them below
  degree 3).

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
