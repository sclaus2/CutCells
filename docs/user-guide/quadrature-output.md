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
physical.

`order` is the polynomial degree integrated exactly on flat pieces, from 1 to
10, with both backends. On curved pieces quadrays converges with the order.
Whole cells, zero faces and straight pieces use reference rules of degree
`order`. Where such a rule has a negative weight (the tetrahedron's rules of
degree 3, 7 and 8), the lowest higher degree with positive weights is used
instead, so all weights stay positive.
quadrays chooses its Gauss-Legendre points per segment from the order; setting
`points_per_segment` on a `QuadraysOptions` fixes them instead, for studies.

Modes:

- `mode="full"`: cut cells plus the whole cells of a volume part.
- `mode="cut_only"`: only cut cells; use it for interfaces and for inspecting
  local decompositions.

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
curve where one of the two level sets lies along mesh faces (a plane through
grid faces, say).

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
`cutcells.cutmesh_to_pyvista` converts to PyVista.

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
