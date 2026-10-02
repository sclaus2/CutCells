# The front end without AdaptCell (phase 3)

`cpp/src/part/`, namespace `cutcells::part`, Python `cutcells.part`. It classifies
cells by their level sets, selects mesh parts by expressions, and hands their
quadrature and visualisation to a backend; quadrays today, the lookup tables in
phase 4. It needs no AdaptCell and no interpolant of an analytic level set. It
sits beside `HOCutResult` and `HOMeshPart` until it replaces them (phase 5).

```python
import cutcells

result = cutcells.part.cut(mesh, cutcells.analytic_sphere([0, 0, 0], 0.7))
rules = result["phi < 0"].quadrature(order=5, backend="quadrays")
result["phi = 0"].write_vtu("sphere.vtu", mode="cut_only", degree=3)
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
| `output.h/.cpp` | `quadrature_rules`, `visualization_mesh` and `write_vtu`, by backend |

## Classification

Every cell is inside (phi < 0), outside (phi > 0) or cut, for every level set,
by bounds of the level set itself:

- **Hexahedra:** Bernstein coefficients of the box form of a Pk level set, or
  Taylor models of an analytic one, over the box and its halves
  (`quadrays::cell_sign`).
- **Tetrahedra:** bounds over the tetrahedron itself and its halves, split at
  the midpoint of the longest edge: the simplex Bernstein coefficients of a Pk
  level set (de Casteljau along the edge), or Taylor models over the
  parallelepiped of a piece's edges whose linear part is bounded at the
  piece's four vertices (intervals over the piece's bounding box when the level
  set has no Taylor models). Bounds over the tetrahedron's box would cover six
  times its volume and flag cells the level set only crosses outside it.
- **Other cells:** the signs of their Bernstein coefficients, without
  subdivision (quadrays takes tetrahedra and hexahedra only; analytic level sets
  are refused there).

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
(`term_on_cell`). With quadrays, the pieces of a cut cell must be bounded by one
level set: a term whose two level sets cut the same cell, or terms bounded by
different level sets in one cell, raise an error until phase 6; terms on one
level set ("phi < 0 or phi > 0") are fine.

## Faces lying in a zero set

Where a level set vanishes on a whole mesh face (a plane through grid faces),
no cell is cut there, and no engine integrated that interface before. `cut()`
finds such faces (the level set is 0 at the face's vertices and, by its bounds,
on the face) and gives each to one cell: the cell on its negative side (by its
domain, or for a cut cell by the derivative into the cell), else the lower cell
index. An interface part then integrates the face once, with the reference
rule of the face merged into the owner's rule, and shows it as a linear face.
For the plane x = 0.25 at n = 8 that is 64 quadrilaterals or 128 triangles, and
the interface measures 4 to rounding, as do the volumes below and above it.

## Quadrature and visualisation

`part.quadrature(order, mode, backend="quadrays", options)` gives one rule per
cell, cells ascending; points in the cell's reference coordinates, physical
weights. quadrays' `order` counts Gauss points per segment; whole cells and
zero faces get the reference rules exact for degree 2 order - 1 (at most 10).
`part.visualization_mesh(mode, backend, degree)` gives quadrays' leaves as
Lagrange cells, zero faces as linear faces and, with mode `full`, the whole cells
as linear cells; `write_vtu` writes them.

## Limits

- Affine cells, as everywhere in quadrays; tetrahedra and hexahedra in 3D for
  analytic level sets and the quadrays backend (2D cells, prisms and pyramids
  come in phase 6).
- Pk level sets need dof values; level sets with nodal values only are refused.
- Several level sets in one cell: phase 6.
