# Quadrature and output

A `MeshPart` gives quadrature rules, a visualisation mesh, or a VTU file.

## Quadrature rules

```python
inside = result[phi < 0]
interface = result[phi == 0]

volume_rules = inside.quadrature(order=4, mode="full")
interface_rules = interface.quadrature(order=4, mode="cut_only")

print(volume_rules.points.reshape(-1, volume_rules.tdim).shape)
print(volume_rules.weights.sum())
```

`QuadratureRules` stores flat arrays for points, weights, offsets and parent
cells: one rule per cell, points in the cell's reference coordinates, weights
physical. On a mesh of facets in a space of higher dimension, points are the
facets' reference coordinates (`tdim` of them) and weights are measured on the
facets.

`order` is the polynomial degree integrated exactly on flat pieces, from 1 to
10, with both backends. On curved pieces quadrays converges with the order.
Whole cells, the facets, edges and vertices a part holds (below) and straight
pieces use reference rules of degree `order`. Where such a rule has a negative
weight (the tetrahedron's rules of degree 3, 7 and 8), the lowest higher
degree with positive weights is used instead, so all weights stay positive.
quadrays chooses its Gauss-Legendre points per segment from the order; setting
`points_per_segment` on a `QuadraysOptions` fixes them instead, for studies.

Modes:

- `mode="full"`: cut cells plus the whole cells of a volume part.
- `mode="cut_only"`: only cut cells; use it for interfaces and for inspecting
  local decompositions.

Both modes include the rules on facets, edges and vertices, as `entities`
says (below).

## Rules on facets, edges and vertices

A piece of a part can lie exactly on a facet, an edge or a vertex of the mesh:
a plane through grid faces, two grid planes meeting along a mesh edge, two
lines crossing at a mesh vertex. Such a piece belongs to every cell sharing
that entity. CutCells integrates it once, on the entity, and the cells keep to
their interiors:

```python
result = cutcells.cut(mesh, plane)          # x = 0.25, on grid faces
face = result[plane == 0]

on = face.entity_quadrature(order=4)        # rules on the entities, once each
f = on.facets                               # faces in 3D, edges in 2D
f.points        # in each facet's reference cell, shape (n, entity_dim)
f.weights       # physical
f.cells         # per rule, the cells sharing the facet (CSR, f.cell_offsets)
f.local_index   # the facet's number in each of them (Basix numbering)
f.negative      # per sharing cell, bit i: the cell lies where level set i < 0
f.cell_points   # the points in each sharing cell's reference coordinates
                # (CSR, f.cell_point_offsets)

rules = face.quadrature(order=4)            # entities="split", the default
```

`entity_quadrature` returns an `EntityQuadrature`: `facets`, `ridges`
(codimension 2: edges in 3D, vertices in 2D) and, in 3D, `vertices`, each an
`EntityRules` with one rule per entity. On a mesh of intervals the facets are
the vertices, and on a mesh of facets in 3D (a surface mesh) its edges. A rule's points lie in the reference
cell of the entity's type with its vertices in the entity's chart: its nodes
sorted, a quadrilateral's starting at its lowest node, then that node's lower
and higher neighbours. `cell_points` gives the same points in each sharing
cell, so no orientation tables are needed.

Which pieces are entity rules:

- interfaces lying on facets, where the level set vanishes on the facet and
  changes sign across it (or the facet lies on the mesh boundary);
- in 3D, curves `a = 0 and b = 0` along mesh edges where both vanish;
- in 2D, the points where two curves cross at a mesh vertex or on a mesh edge.

A level set that vanishes on a facet without changing sign across it, such as
`min(box1, box2)` on the face two boxes share, is no surface there. Where
another level set of the selection cuts the entity (a disk `plane == 0` with
`ball < 0`, its circle on grid faces), the piece is a part of the entity,
which is not integrated yet.

`quadrature(order, entities=...)` folds the entity rules into the cell rules:

| `entities` | Each sharing cell gets |
| --- | --- |
| `"split"` (default) | the rule with weight 1/n for n sharing cells: 1/2 per cell for a facet between two cells, 1 on the mesh boundary |
| `"negative"`, `"positive"` | the cells on that side of the zero level sets share it; on the mesh boundary, all |
| `"owner"` | the lowest sharing cell takes it whole |
| `"exclude"` | nothing: only the pieces inside cells |

Each piece counts once in total with every option. `split` suits integrands
that are continuous across the entity; for one-sided traces or jumps, take
the cells' sides from `entity_quadrature` instead.

## Backends

The default backend is `quadrays`; `backend="lut"` (also `"straight"`) gives
straight pieces:

```python
rules = inside.quadrature(order=4)                 # the part's backend, quadrays by default
straight = inside.quadrature(order=4, backend="lut")

options = cutcells.QuadraysOptions()
options.margin = 0.3
tuned = inside.quadrature(order=4, options=options)
```

The backend of a part is the backend given to `cut`, or set on the result
before selecting (`result.backend = "lut"`). Both backends integrate the curves
where two level sets vanish, with points of unit weight in 2D. Neither finds a
curve where only one of the two level sets lies in mesh faces (a plane through
grid faces cut by a sphere, say); where both do, the curve runs along mesh
edges and comes from the entity rules.

## Fewer points

`compress_rules` replaces every rule by at most as many of its points as the
polynomial space of a given degree has, keeping that degree exact:

```python
small, stats = cutcells.compress_rules(volume_rules, degree=4)
```

## Visualisation meshes

```python
inside_mesh = inside.visualization_mesh(mode="full", degree=2)
interface_mesh = interface.visualization_mesh(mode="cut_only")
```

With quadrays, the result is a `QuadraysLeafMesh` of Lagrange cells of the
given degree (arrays `points`, `connectivity`, `offsets`, `vtk_types`,
`parent`). With the lookup tables, it is a `CutMesh` of straight cells, which
`cutcells.cutmesh_to_pyvista` converts to PyVista. Both draw each piece on a
facet, an edge or a vertex once, as a linear cell whose parent is the lowest
cell sharing it.

## VTU output

Write selected parts directly:

```python
inside.write_vtu("inside_full.vtu", mode="full", degree=2)
inside.write_vtu("inside_cut_only.vtu", mode="cut_only")
interface.write_vtu("interface.vtu", mode="cut_only")
```

For level-set diagnostics on the background mesh, use:

```python
cutcells.write_level_set_vtu("level_set.vtu", phi, field_name="phi")
```

See `python/demo/demo_mesh_parts.py` and `python/demo/demo_two_level_sets.py`
for complete workflows that write several selected parts.

## Lower-level paths

`runtime_quadrature` works directly on VTK-style arrays and level-set values,
with straight cuts:

```python
rules = cutcells.runtime_quadrature(
    ls_vals,
    points_flat,
    connectivity,
    offsets,
    vtk_types,
    "phi<0",
    triangulate=True,
    order=4,
)
physical = cutcells.physical_points(
    rules,
    points_flat,
    connectivity,
    offsets,
    vtk_types,
).reshape(-1, 3)
```

It is useful when coupling to an assembly backend that already owns VTK-style
mesh arrays.
