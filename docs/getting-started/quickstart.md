# Quickstart

This example cuts a structured triangular mesh by a circular level set and
integrates the selected region.

```python
import numpy as np
import cutcells

mesh = cutcells.structured_triangle_mesh_view(
    -1.0, -1.0,
     1.0,  1.0,
    41, 41,
)

def phi(X):
    # create_level_set passes X with shape (gdim, num_points).
    return np.sqrt(X[0] ** 2 + X[1] ** 2) - 0.55

level_set = cutcells.create_level_set(mesh, phi, degree=2, name="phi")
result = cutcells.cut(mesh, level_set)

inside = result["phi < 0"]
interface = result["phi = 0"]
outside = result["phi > 0"]

print("cut parent cells:", result.num_cut_cells)
print("inside:   ", inside.num_cut_cells, "cut cells")
print("interface:", interface.num_cut_cells, "cut cells")
print("outside:  ", outside.num_cut_cells, "cut cells")

inside_rules = inside.quadrature(order=4, mode="full")
interface_rules = interface.quadrature(order=4, mode="cut_only")

print("area approximation:     ", inside_rules.weights.sum())
print("perimeter approximation:", interface_rules.weights.sum())
```

The important pieces are:

- `MeshView`: a lightweight background mesh view.
- `create_level_set`: interpolation of an implicit function on the mesh.
- `cut`: classification of every cell by every level set.
- `result["phi < 0"]`: selection of a region, interface, or intersection.
- `MeshPart.quadrature`: quadrature rules for selected mesh parts; `order` is
  the polynomial degree integrated exactly on flat pieces.

By default the parts integrate with quadrays, which follows the curved level
set to high order. For straight pieces, cut with the lookup tables:
`cutcells.cut(mesh, level_set, backend="lut")`.

## Exact geometry

An analytic level set skips the interpolant altogether:

```python
circle = cutcells.analytic_sphere([0.0, 0.0, 0.0], 0.55)
result = cutcells.cut(mesh, circle)   # named "phi"
area = result["phi < 0"].quadrature(order=4).weights.sum()
```

## Write a visualisation mesh

The selected part can also be written as VTU:

```python
inside.write_vtu("inside.vtu", mode="full", degree=2)
interface.write_vtu("interface.vtu", mode="cut_only")
```

With quadrays the cut pieces are Lagrange cells of the given `degree`. Use
`mode="full"` for a volume region when you want both uncut cells and cut
pieces. Use `mode="cut_only"` when inspecting only the cells touched by the
level set or when writing interface parts.

## PyVista input

If your mesh already lives in PyVista, convert it to a `MeshView`:

```python
from cutcells import box_tetrahedron_mesh, mesh_from_pyvista

grid = box_tetrahedron_mesh(-1, -1, -1, 1, 1, 1, 8, 8, 8)
mesh = mesh_from_pyvista(grid, tdim=3)
```

See `python/demo/demo_mesh_parts.py` for a complete PyVista-based workflow.
