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

## Where the data lives

A `MeshView` does not copy the mesh. It points into the NumPy arrays it was
given and keeps a reference to them, so they stay alive as long as the view
does. Its `coordinates`, `connectivity` and `offsets` properties are those
same arrays, and editing an array in place changes the mesh the view
describes. Arrays of another dtype or layout, such as `int64` connectivity,
are converted once when the view is made; the view then holds the converted
copy. The view converts the VTK cell types into its own small array.

This makes a `MeshView` cheap to put over a mesh that already lives
elsewhere, such as the geometry of a DOLFINx mesh. The arrays of a mesh built
with `rectangle_mesh` or `box_mesh` are held by the view alone, so there the
view is in effect the mesh. `mesh_from_pyvista` may share the grid's point
array.

## Facets and surface meshes

The coordinates may have more components than the cells' dimension: the
facets of a mesh, or a surface mesh, are triangles or quadrilaterals with
`tdim=2` in 3D, or intervals with `tdim=1` in 2D or 3D. Each cell is cut and
integrated as the cell it is, so weights are areas and lengths on the cells
and points are in their own reference coordinates:

```python
# two triangles of the plane z = x / 2 in 3D
coordinates = np.array(
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.5], [1.0, 1.0, 0.5], [0.0, 1.0, 0.0]]
)
facets = cutcells.MeshView(
    coordinates=coordinates,
    connectivity=np.array([0, 1, 2, 0, 2, 3], dtype=np.int32),
    offsets=np.array([0, 3, 6], dtype=np.int32),
    cell_types=np.array([5, 5], dtype=np.int32),
    tdim=2,
)
ball = cutcells.analytic_sphere([0.5, 0.5, 0.25], 0.3)
result = cutcells.cut(facets, ball)
rules = result[ball < 0].quadrature(order=4)   # rules.tdim == 2
area = rules.weights.sum()                      # the disk on the plane
```

## Structured mesh helpers

CutCells provides lightweight helpers for common examples. They build a
`MeshView` on a rectangle or a box with NumPy alone; the arguments are the
corners and the numbers of grid points per direction:

```python
tri = cutcells.rectangle_mesh(-1, -1, 1, 1, 41, 41, cell_type="triangle")
quad = cutcells.rectangle_mesh(-1, -1, 1, 1, 41, 41, cell_type="quadrilateral")
tet = cutcells.box_mesh(-1, -1, -1, 1, 1, 1, 13, 13, 13, cell_type="tetrahedron")
hexes = cutcells.box_mesh(-1, -1, -1, 1, 1, 1, 13, 13, 13, cell_type="hexahedron")
prisms = cutcells.box_mesh(-1, -1, -1, 1, 1, 1, 13, 13, 13, cell_type="prism")
pyramids = cutcells.box_mesh(-1, -1, -1, 1, 1, 1, 13, 13, 13, cell_type="pyramid")
```

The tetrahedral mesh splits every grid cube into six tetrahedra around its
diagonal, so it is conforming. The prismatic mesh splits every cube into two
prisms along z, the pyramidal mesh into six pyramids with their apex at its
centre, a node of its own.

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

The other way, `cutcells.to_pyvista` converts a `MeshView`, or the
visualisation mesh of a selected part, to a PyVista grid:

```python
grid = cutcells.to_pyvista(mesh)
pieces = cutcells.to_pyvista(part.visualization_mesh(mode="cut_only", degree=3))
surface = pieces.extract_surface(nonlinear_subdivision=3, algorithm="dataset_surface")
```

The pieces of quadrays are Lagrange cells; `nonlinear_subdivision` draws them
curved.

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
