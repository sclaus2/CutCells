# Mesh views

`MeshView` is the Python representation of a background mesh. It is deliberately
small: CutCells needs coordinates, cell connectivity, offsets, cell types, and
topological dimension.

## Array layout

Create a mesh view directly from NumPy arrays:

```python
import numpy as np
import cutcells

coordinates = np.array(
    [
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
    ],
    dtype=np.float64,
)
connectivity = np.array([0, 1, 2], dtype=np.int32)
offsets = np.array([0, 3], dtype=np.int32)
cell_types = np.array([5], dtype=np.int32)  # VTK_TRIANGLE

mesh = cutcells.MeshView(
    coordinates=coordinates,
    connectivity=connectivity,
    offsets=offsets,
    cell_types=cell_types,
    tdim=2,
)
```

The arrays follow the VTK unstructured-grid convention:

- `coordinates`: shape `(num_nodes, gdim)`.
- `connectivity`: flattened cell-to-node connectivity.
- `offsets`: start offset for each cell plus a final sentinel offset.
- `cell_types`: one VTK cell type ID per cell.
- `tdim`: topological dimension of the background cells.

CutCells currently uses these VTK cell IDs:

| Cell | VTK ID | Dimension |
| --- | ---: | ---: |
| interval | 3 | 1 |
| triangle | 5 | 2 |
| quadrilateral | 9 | 2 |
| tetrahedron | 10 | 3 |
| hexahedron | 12 | 3 |
| prism/wedge | 13 | 3 |
| pyramid | 14 | 3 |

Use C-contiguous arrays with numeric dtypes (`float64` for coordinates and
`int32` for connectivity metadata) unless you have a reason to use the
`float32` bindings.

## Structured mesh helpers

CutCells provides lightweight helpers for common examples:

```python
mesh = cutcells.structured_triangle_mesh_view(-1, -1, 1, 1, 41, 41)
```

The following helpers return PyVista grids and require `pyvista`:

```python
quad_grid = cutcells.rectangle_quad_mesh(-1, -1, 1, 1, 20, 20)
tet_grid = cutcells.box_tetrahedron_mesh(-1, -1, -1, 1, 1, 1, 8, 8, 8)
hex_grid = cutcells.box_hex_mesh(-1, -1, -1, 1, 1, 1, 8, 8, 8)
```

Convert a PyVista `UnstructuredGrid` to a `MeshView` with:

```python
mesh = cutcells.mesh_from_pyvista(tet_grid, tdim=3)
```

If `tdim` is not supplied, `mesh_from_pyvista` tries to infer it from the first
cell type.

## Inspecting a mesh

Useful `MeshView` properties and methods:

```python
print(mesh.gdim, mesh.tdim)
print(mesh.num_nodes(), mesh.num_cells())
print(mesh.coordinates.shape)
print(mesh.connectivity.shape, mesh.offsets.shape)
print(mesh.cell_num_nodes(0))
print(mesh.cell_node(0, 0))
```

`MeshView` is intended as an adapter at the CutCells boundary. Keep ownership
and lifetime of source mesh data clear in your application, especially when
building mesh views from external libraries.
