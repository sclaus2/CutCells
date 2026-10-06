# Algorithms and scope

CutCells classifies every cell of a background mesh by every level set, selects
mesh parts by expressions, and integrates or visualises the selected parts
with one of two backends. It takes Pk level sets with degree-of-freedom values
and analytic level sets, and any number of level sets per cell.

## Supported parent cell types

| Dimension | Cells |
| ---: | --- |
| 1 | interval |
| 2 | triangle, quadrilateral |
| 3 | tetrahedron, hexahedron, prism, pyramid |

Pk level sets of any degree are supported; pyramids take degree 2 at most.

A mesh's cells may lie in a space of higher dimension than their own: the
facets of another mesh taken as cells (triangles and quadrilaterals in 3D,
intervals in 2D or 3D), its edges, or a surface mesh. Each such cell is cut and
integrated as the cell it is: quadrays maps it to an orthonormal frame of its
plane or line, so that areas and lengths are measured on it, points are its
own reference coordinates, and an analytic level set's derivative normal to it
does not count. Its own facets, edges and vertices (below) are those of the
cell: the edges of a triangle in 3D, the vertices of an interval.

## Classification

`cut(mesh, level_sets)` classifies each cell as inside (`phi < 0`), outside
(`phi > 0`) or cut, for every level set, by bounds of the level set itself:
Bernstein coefficients of a Pk level set, or Taylor models of an analytic one,
over the cell and its halves. A cell whose sign the bounds cannot prove after
`max_depth` bisections counts as cut, and the backend finds out what it holds.

## Facets, edges and vertices

Every point of the mesh lies in the interior of exactly one entity: a cell, a
facet, an edge or a vertex. A piece of a part is integrated on that entity,
once, and cells keep to their interiors.

Where a level set vanishes on a whole mesh facet (a plane through grid faces,
for example), no cell is cut there. `cut` finds the facets, edges and
vertices on which a level set vanishes, each with every cell that shares it.
It decides once per entity, from the level set's restriction to the entity:
the restriction of its Bernstein form, or Taylor bounds over the entity itself
for an analytic level set. One threshold, 64 eps relative to
the level set's largest magnitude over the cells at its zero set, decides
both that and whether the cells beside the entity touch zero there, so a
plane slightly off grid faces is either on the faces or cuts the cells next
to them, never both. It counts as surface only where it changes sign across
the entity (or on the mesh boundary).

Rules on those entities are made once, in the entity's reference cell, and
come into the cells as the selection's `entities` option says (by default
1/n to each of the n cells sharing an entity). The backends keep to the
cells' interiors: the lookup tables drop pieces lying in a facet or an edge on
which their zero sets vanish, quadrays leaves level sets vanishing on its
boxes' faces and curves lying in a cell's faces to the entities, and the
points where two curves cross on 2D cells go to the edge or vertex holding
them.

## Backends

The parts of a result integrate and visualise with a backend:

- **`"quadrays"`** (the default): a certify-and-bisect engine that reduces the
  dimension one height direction at a time. It integrates the level sets
  themselves, curved pieces included, to high order. It takes any number of
  level sets per cell, follows the edges where two zero sets meet and the
  corners where three meet, and works on affine cells: intervals, triangles,
  parallelograms, tetrahedra, parallelepipeds, and prisms and pyramids with
  parallelogram faces. It integrates volumes, interfaces, and the curves
  where two level sets vanish (points with unit weights in 2D). On intervals
  an interface is points with unit weights.
- **`"lut"`** (also `"straight"`): the lookup tables. A cut cell is subdivided
  by a Pk-iso-P1 template, and each sub-cell is cut by the level sets' P1
  interpolants. The pieces are straight; they include the curves (points in
  2D) where two level sets vanish.

Choose the backend when cutting, or per part:

```python
result = cutcells.cut(mesh, phi)                    # quadrays
straight = cutcells.cut(mesh, phi, backend="lut")   # lookup tables
rules = result[phi < 0].quadrature(order=4, backend="lut")
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
