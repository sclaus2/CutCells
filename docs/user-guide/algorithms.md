# Algorithms and scope

CutCells classifies every cell of a background mesh by every level set, selects
mesh parts by expressions, and integrates or visualises the selected parts
with one of two backends. It takes Pk level sets with degree-of-freedom values
and analytic level sets, and any number of level sets per cell.

## Supported parent cell types

| Dimension | Cells |
| ---: | --- |
| 1 | interval (lookup tables) |
| 2 | triangle, quadrilateral |
| 3 | tetrahedron, hexahedron, prism, pyramid |

Pk level sets of any degree are supported; pyramids take degree 2 at most.

## Classification

`cut(mesh, level_sets)` classifies each cell as inside (`phi < 0`), outside
(`phi > 0`) or cut, for every level set, by bounds of the level set itself:
Bernstein coefficients of a Pk level set, or Taylor models of an analytic one,
over the cell and its halves. A cell whose sign the bounds cannot prove after
`max_depth` bisections counts as cut, and the backend finds out what it holds.

Where a level set vanishes on a whole mesh face (a plane through grid faces,
for example), no cell is cut there. `cut` finds such faces and gives each one
to a single cell, so interface parts integrate it exactly once.

## Backends

The parts of a result integrate and visualise with a backend:

- **`"quadrays"`** (the default): a certify-and-bisect engine that reduces the
  dimension one height direction at a time. It integrates the level sets
  themselves, curved pieces included, to high order. It takes any number of
  level sets per cell, follows the edges where two zero sets meet and the
  corners where three meet, and works on affine cells: triangles,
  parallelograms, tetrahedra, parallelepipeds, and prisms and pyramids with
  parallelogram faces. It integrates volumes, interfaces, and the curves
  where two level sets vanish (points with unit weights in 2D).
- **`"lut"`** (also `"straight"`): the lookup tables. A cut cell is subdivided
  by a Pk-iso-P1 template, and each sub-cell is cut by the level sets' P1
  interpolants. The pieces are straight; they include the curves (points in
  2D) where two level sets vanish.

Choose the backend when cutting, or per part:

```python
result = cutcells.cut(mesh, phi)                    # quadrays
straight = cutcells.cut(mesh, phi, backend="lut")   # lookup tables
rules = result["phi < 0"].quadrature(order=4, backend="lut")
```

Options go through `QuadraysOptions` and `LutOptions` (see the
[API overview](../reference/python-api.md)).

## Triangulation (lookup tables)

The lookup tables can split their pieces into simplices:

```python
result = cutcells.cut(mesh, phi, backend="lut", triangulate=True)
# or: cutcells.cut(mesh, phi, backend="lut", options=cutcells.LutOptions(triangulate=True))
```

`triangulate=True` gives simplex pieces, often the most convenient choice for
visualisation. `triangulate=False` keeps the native cut topology where the
tables allow it. These keywords belong to the lookup tables; given without
`backend="lut"`, they raise.

## Generated tables and third-party data

Several tables are generated offline:

- VTK-derived clip and cut case tables for hexahedra, prisms and pyramids.
- Basix-generated quadrature tables.

A pinned Algoim snapshot under `third_party/algoim` is used only by the
comparison drivers in `benchmarks/`; the library does not depend on it. The
repository includes license and provenance information under
`third_party/THIRD_PARTY_NOTICES.md`.

## Relationship to CutFEMx

CutFEMx uses CutCells as its cut-geometry backend. CutCells is intentionally
independent of DOLFINx assembly concepts: it provides geometry, selected mesh
parts, parent mappings and quadrature data. CutFEMx combines those objects
with FEniCSx forms, runtime quadrature providers, ghost penalties and
degree-of-freedom constraints.
