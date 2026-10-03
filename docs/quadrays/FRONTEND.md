# The front end (phases 3 to 6)

`cpp/src/part/`, namespace `cutcells::part`, Python `cutcells.cut` and
`cutcells.part`. It classifies cells by their level sets, selects mesh parts by
expressions, and hands their quadrature and visualisation to a backend:
quadrays, or the lookup tables on Pk-iso-P1 templates (`cpp/src/lut/`, phase 4).
It needs no interpolant of an analytic level set. Since phase 5 it is the only
front end: AdaptCell, its certification and `HOCutResult`/`HOMeshPart` of old
are gone, and those names now denote this front end's classes.

```python
import cutcells

result = cutcells.cut(mesh, cutcells.analytic_sphere([0, 0, 0], 0.7))
rules = result["phi < 0"].quadrature(order=2)                     # lookup tables
curved = result["phi < 0"].quadrature(order=5, backend="quadrays")
result["phi = 0"].write_vtu("sphere.vtu", mode="cut_only")
```

`cut(mesh, level_sets)` takes LevelSetFunctions with dof values (Pk) and
AnalyticLevelSets, alone or in a list; analytic ones are named `phi`, or `phi1`,
`phi2`, ... in a list, unless `names` are given.

## Files

| File | Contents |
| --- | --- |
| `cell_source.h/.cpp` | one cell and one level set as quadrays reads them: the cell's clipped box and a Bernstein form on it, or the analytic level set |
| `classify.h/.cpp` | the sign of a level set on a cell by bounds; `bisect_tetrahedron` splits simplex Bernstein forms |
| `cut_result.h/.cpp` | `CutResult` (a domain per cell and level set, the cut cells, the faces lying in zero sets with their owners) and `cut()` |
| `mesh_part.h/.cpp` | `MeshPart` (whole cells, cut cells, zero faces) and `select()`; `term_on_cell` |
| `output.h/.cpp` | `quadrature_rules`, `visualization_mesh` and `write_vtu`; the type of the options (`quadrays::Options`, `lut::Options`) picks the backend |

## Classification

Every cell is inside (phi < 0), outside (phi > 0) or cut, for every level set,
by bounds of the level set itself:

- **Hexahedra and quadrilaterals:** Bernstein coefficients of the box form of
  a Pk level set, or Taylor models of an analytic one, over the box and its
  halves (`quadrays::cell_sign`).
- **Tetrahedra and triangles:** bounds over the simplex itself and its halves,
  split at the midpoint of the longest edge: the simplex Bernstein coefficients
  of a Pk level set (de Casteljau along the edge), or Taylor models over the
  parallelepiped of a piece's edges whose linear part is bounded at the
  piece's vertices (intervals over the piece's bounding box when the level
  set has no Taylor models). Bounds over the tetrahedron's box would cover six
  times its volume and flag cells the level set only crosses outside it.
- **Prisms and pyramids:** as hexahedra, over their clipped boxes and halves
  (sub-boxes that miss the cell drop out). A level set vanishing on an edge
  where the clip meets the box leaves the sign unproven: such cells count as
  cut (a plane in mesh faces flags the pyramids touching it along an edge).

A sign the bounds cannot prove after `max_depth` (12) splits counts as cut; the
backend then finds what the cell holds. On the sphere of radius 0.7 at n = 8,
the P2 level set flags exactly the cells the sphere crosses (138 hexahedra, 627
tetrahedra, each with interface); the signed distance flags 4 tetrahedra more.
A ball of radius 0.26 inside a cell of side 0.5 flags exactly the 18
tetrahedra that dense sampling finds crossed.

Seconds for `cut()` on a level set |x - c|^2 - r^2 + 0.05 sin(4x) cos(3y) sin(5z),
interpolated to degree k, and on the analytic distance (load average 8):

| | k = 1 | k = 2 | k = 3 | k = 4 | distance |
| --- | --- | --- | --- | --- | --- |
| 32,768 hexahedra, `part.cut` | 0.05 | 0.11 | 0.33 | 1.34 | 0.11 |
| same, AdaptCell `cut` | 0.46 | 1.71 | 5.69 | 18.69 | |
| 24,576 tetrahedra, `part.cut` | 0.04 | 0.03 | 0.09 | 0.16 | 0.07 |
| same, AdaptCell `cut` | 0.42 | 0.93 | 2.39 | 8.18 | |

Classification costs 1 to 41 microseconds per cell; the subdivision depth of the
sign tests does not show at degree 4.

## Parts

`result[expr]` evaluates each term of the expression on each cell from its
domains: a term holds on the whole cell, on a piece of it bounded by the level
sets that cut it, or nowhere. A cell is in the part whole if some term holds on
all of it, otherwise it is a cut cell of the part if some term holds on a piece
(`term_on_cell`). Both backends take any number of level sets per cell (quadrays
since phase 6); quadrays integrates volumes and interfaces, the curves where two
level sets vanish come from the lookup tables.

## Faces lying in a zero set

Where a level set vanishes on a whole mesh facet (a plane through grid faces, a
line through grid edges in 2D), no cell is cut there, and no engine integrated
that interface before. `cut()` finds such facets and gives each to one cell:
the cell on its negative side (by its domain, or for a cut cell by the
derivative into the cell), else the lower cell index. On tetrahedra and
hexahedra the level set is 0 at the face's vertices and, by its bounds, on the
face, and so on the other cells quadrays takes (phase 6) by the bounds of their
box forms on the facet. An interface part then integrates the face once, with the reference
rule of the face merged into the owner's rule, and shows it as a linear face.
For the plane x = 0.25 at n = 8 that is 64 quadrilaterals or 128 triangles, and
the interface measures 4 to rounding, as do the volumes below and above it.

## Quadrature and visualisation

`part.quadrature(order, mode, backend="quadrays", options)` gives one rule per
cell, cells ascending; points in the cell's reference coordinates, physical
weights. quadrays' `order` counts Gauss points per segment; whole cells, zero
faces and the lookup tables' pieces get the reference rules exact for degree
2 order - 1 (at most 10). `options` is a `QuadraysOptions` for quadrays and a
`LutOptions` for `backend="lut"` (None: the defaults).
`part.visualization_mesh(mode, backend, degree)` gives quadrays' leaves as
Lagrange cells (a `QuadraysLeafMesh`) or the lookup tables' straight pieces (a
`CutMesh`), zero faces as linear faces and, with mode `full`, the whole cells as
linear cells; `write_vtu` writes them.

## The lookup-table backend (phase 4)

`cpp/src/lut/`, namespace `cutcells::lut`:

| File | Contents |
| --- | --- |
| `cell_pieces.h/.cpp` | `cut_cell`: the straight pieces of a cell on its Pk-iso-P1 template, each with its side of every level set; `template_vertices` |
| `piece_rules.h/.cpp` | `CellMap` (affine, or multilinear on quadrilaterals and hexahedra), `push_forward`, and `append_piece_rule`: rules on straight pieces through the cell's map |

For each cut cell of a part, `backend="lut"`

1. takes the level sets that cut the cell and that the expression names;
2. picks the Pk-iso-P1 template of order k: `LutOptions.template_order`, by
   default the highest degree among them (2 for analytic level sets), 1 to 4;
3. evaluates them at the template's vertices: the Pk polynomial, or the
   analytic level set itself at the mapped point;
4. cuts each sub-cell with the lookup tables (`cell::cut`) by the level sets'
   P1 interpolants, one level set after the other, both sides kept; the zero
   sets the expression asks for are cut by the other level sets;
5. integrates the pieces that some term selects: by their sides for the
   cutting level sets, by the cell's domains for the others. Every piece is
   counted once, so unions and complements add up exactly.

Values within 64 eps max|v| of 0 count as positive; they are moved to that
bound, so the tables' case masks (which treat values within 2 eps |v_0| as 0)
agree with the classification. A level set vanishing on a face of the sub-cells
is integrated once, by the sub-cell below it, as zero faces of mesh cells are.

What the tables need, found on the way:

- The tables of quadrilaterals, hexahedra, prisms and pyramids take and give
  Basix vertex order; those of triangles and tetrahedra give their
  quadrilaterals in cyclic (VTK) order. The backend turns them into Basix order.
- `cut_tetrahedron`'s table listed the prism of case 13 (vertex 1 alone on its
  side) with its top triangle rotated: a twisted wedge, which split into
  tetrahedra measured 0.094 instead of 0.133 in one case; fixed in
  `cpp/src/cut_tetrahedron.cpp`. Its triangulated output was right before (it
  rebuilds the prism from the intersected edges).
- The hexahedron's tables cut parallelepipeds exactly for affine values, but
  not the other hexahedra that a first cut leaves (errors up to 0.09 of the
  unit cube in 400 pairs of random planes), and for multilinear values their
  two sides do not fit together (up to 0.8% of a cell). The backend cuts
  leaves that are not parallelograms or parallelepipeds as simplices, and
  hexahedral sub-cells on which some level set is not affine as their Kuhn
  tetrahedra (for affine values the tables give the same measures as those).
  The quadrilateral's tables fit together, also at saddles.
- Rules: simplices take the simplex rules, parallelograms and parallelepipeds
  the tensor rules, other quadrilaterals and hexahedra, prisms and pyramids are
  split into simplices; so on pieces with planar faces in affine cells the
  rules are exact for polynomials of degree 2 order - 1. Areas come from cross
  products: the Gram determinant lost half the digits on the slivers next to
  values moved off 0 (6.6e-10 instead of 2e-15).

### Against the straight backend of AdaptCell

The straight backend that phase 5 removed refined cut cells with the same
iso-Pk templates and cut them with the same tables. Its numbers for these cases
are kept in `python/tests/test_part_lut.py`. On [-1, 1]^3 and [-1, 1]^2 at n = 6 and 8 the lookup tables give every
cell the same measure and first moments, to 1e-15, for P1 and P2 spheres
(circles) on hexahedra, tetrahedra, quadrilaterals and triangles, volumes and
interfaces, and for a P2 sphere and a plane crossing in cells on tetrahedra,
quadrilaterals and triangles (`python/tests/test_part_lut.py`). Two
differences: the straight backend refuses two level sets on hexahedra
(`refine_red_on_ambiguous_cells: unsupported cell type`), and where the
values on hexahedra are multilinear, its two sides do not fill the cells: for
the Q1 and Q2 interpolants of the sphere's distance at n = 8, 48 and 94 of
the 138 cut cells miss up to 1.3e-4 of their volume 0.0156 (0.8%), which the
lookup tables fill to rounding.

### Template order

The sphere of radius 0.7 by its analytic distance, its values at the template's
vertices, relative errors of the ball's volume and the sphere's area
(`python/tests/test_part_lut.py`, `cpp/tests/lut/test_lut.cpp`); hexahedra and
tetrahedra give the same numbers here, since the distance is multilinear on the
hexahedra's sub-cells, which are cut as Kuhn tetrahedra:

| n | k = 1 | k = 2 | k = 3 | k = 4 | rate in k |
| --- | --- | --- | --- | --- | --- |
| 4 | 2.5e-1 / 1.4e-1 | 6.3e-2 / 3.3e-2 | 2.8e-2 / 1.5e-2 | 1.6e-2 / 8.2e-3 | 1.98 / 2.05 |
| 8 | 6.3e-2 / 3.3e-2 | 1.6e-2 / 8.2e-3 | 7.1e-3 / 3.7e-3 | 4.0e-3 / 2.1e-3 | 2.00 / 2.01 |
| 16 | 1.6e-2 / 8.2e-3 | 4.0e-3 / 2.1e-3 | 1.8e-3 / 9.1e-4 | 1.0e-3 / 5.1e-4 | 2.00 / 2.00 |

The errors fall as (h / k)^2: a mesh of size h with template order k gives the
errors of size h / k with order 1.

### Time

A sphere interpolated to degree k, `phi < 0` (mode full) and `phi = 0`, seconds
for the classification and both rules (best of three, load average 6 to 8):

| | AdaptCell `cut` + straight | `part.cut` + lut |
| --- | --- | --- |
| 32,768 hexahedra, P1 | 0.57 + 0.08 | 0.06 + 0.25 |
| same, P2 | 2.16 + 0.24 | 0.13 + 0.72 |
| same, P3 | 6.34 + 0.63 | 0.34 + 1.37 |
| 24,576 tetrahedra, P1 | 0.12 + 0.04 | 0.02 + 0.07 |
| same, P2 | 0.67 + 0.07 | 0.03 + 0.20 |
| same, P3 | 2.20 + 0.17 | 0.06 + 0.52 |

The lookup tables cut the cells when the rules are asked for, on every call;
the straight backend read the pieces AdaptCell had made in `cut`.

### Curves and points

Where a part's terms name two zero sets ("a = 0 and b = 0": curves in 3D,
points in 2D), the zero pieces of the first are also cut by the zero set of the
second; every piece carries the mask of the zero sets it lies in. A curve on a
face between two sub-cells is counted once, by the sub-cell below the second
level set. Two planes give their line, in every cell type and template order, to
5e-15 (`cpp/tests/lut/test_lut.cpp`).

## Phase 5: AdaptCell retired

- **Removed from the library:** `adapt_cell`, `refine_cell`, `entity_numbering`,
  `cell_certification`, `edge_certification`, `ho_cut_mesh` (`HOCutCells`,
  `ParentCellClassification`, the AdaptCell `cut`) and `ho_mesh_part_output`.
  `iso_refine` keeps its templates without `apply_iso_refine`, and
  `level_set_cell` no longer includes AdaptCell.
- **Moved:** the lookup tables (`cut_<cell>`, their generated tables,
  `triangulation`, the midpoint splits, `cell_subdivision`, `cut_cell`,
  `cut_mesh`, `iso_refine`) into `cpp/src/lut/` (installed under
  `include/cutcells/lut/`); `algoim_quadrature` and `edge_root` into
  `benchmarks/quadrays/` as `mesh_part_algoim` (`benchmarks::algoim_rules` on
  the front end's parts, checked by `mesh_part_algoim_check` with
  `CUTCELLS_WITH_ALGOIM`). The library no longer includes algoim or links
  LAPACK.
- **Python:** `cutcells.cut(mesh, level_sets)` (and `ho_cut`) is `part.cut`
  with the lookup tables as the default backend ("straight" names them too) and
  the keywords of the former cut(): `triangulate`, `triangulation`
  ('classical', 'midpoint'), `cut_approximation` ('auto', 'linear', 'iso_p1')
  with `cut_approximation_order`, `degree` (the template order for analytic
  level sets) and `name`. `result.backend` and `result.options` set the default
  of the parts selected afterwards; a call can name another. `HOCutResult` and
  `HOMeshPart` are `part.CutResult` and `part.MeshPart`, with `parent_cell_ids`,
  `cell_domains`, `num_level_sets`, `cut_cell_ids` and `uncut_cell_ids` kept.
  Gone: `AdaptCell`, `adapt_cell()`, the certification and refinement
  functions, their tags, the AdaptCell tuning keywords
  (`max_refinement_iterations`, `edge_max_depth`, `linear_fast_path`) and the
  algoim backends.
- **Behaviour:** one rule per cell (the straight backend gave one per piece);
  the lookup tables' `order` counts Gauss points per direction (rules exact for
  degree 2 order - 1, the straight backend took the degree); hexahedra with
  multilinear values fill their cells exactly (above); 1D meshes and curves go
  through the lookup tables.
- **Fixed on the way:** `create_level_set_mesh_data(mesh, degree)` gave intervals
  their vertex dofs twice.

## Phase 6: 2D cells, prisms and pyramids, several level sets per cell

- **Cells:** quadrays takes triangles and quadrilaterals in 2D, tetrahedra,
  hexahedra, prisms and pyramids in 3D (`quadrays_takes`), with Pk and analytic
  level sets. Pk level sets on prisms and pyramids have Bernstein forms now
  (`bernstein.h`); `create_level_set` takes pyramids up to degree 2 (Basix's
  rational pyramid space). Faces lying in a zero set are found on edges (2D)
  and faces (3D).
- **Several level sets per cell:** a part's terms on a cell are grouped by the
  level set whose zero set they integrate (`cell_terms`): volumes in one engine
  run, interfaces one run per zero level set. Only the curves (points in 2D)
  where two level sets vanish are left to the lookup tables.
- **Python:** the quadrays options gain `split_bounds`, `rotation_depth`,
  `taylor_subdivisions` and `two_roots_depth`, the stats `two_roots` and
  `surfaces`.

## Limits

- quadrays: affine cells (parallelograms, parallelepipeds, prisms and pyramids
  with parallelogram faces); volumes and interfaces, with any number of level
  sets per cell. Pyramids with Pk level sets bisect towards their apex where a
  zero set passes near it.
- The lookup tables: intervals, triangles, quadrilaterals, tetrahedra and
  hexahedra (the cells of the iso-P1 templates), affine or multilinear cell
  maps; volumes, interfaces and the curves (points in 2D) where two level sets
  vanish. No prisms or pyramids yet.
- Pk level sets need dof values; level sets with nodal values only are refused.
  Pyramids take degree 2 at most.
