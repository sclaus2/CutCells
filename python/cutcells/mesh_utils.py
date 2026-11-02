"""Structured mesh generation and pyvista ↔ MeshView conversion utilities."""

import numpy as np

try:
    import pyvista as pv
except ImportError:
    pv = None


# ---------------------------------------------------------------------------
# VTK cell-type → topological dimension lookup
# ---------------------------------------------------------------------------

_VTK_TDIM = {
    3: 1,  # VTK_LINE
    5: 2,  # VTK_TRIANGLE
    9: 2,  # VTK_QUAD
    10: 3,  # VTK_TETRA
    12: 3,  # VTK_HEXAHEDRON
    13: 3,  # VTK_WEDGE (prism)
    14: 3,  # VTK_PYRAMID
}


# ---------------------------------------------------------------------------
# pyvista ↔ MeshView
# ---------------------------------------------------------------------------


def mesh_from_pyvista(grid, *, tdim=None):
    """Create a cutcells.MeshView from a pyvista UnstructuredGrid.

    Parameters
    ----------
    grid : pyvista.UnstructuredGrid
        Source mesh.
    tdim : int, optional
        Topological dimension. If None, inferred from VTK cell types.

    Returns
    -------
    cutcells.MeshView
    """
    from . import MeshView

    if pv is None:
        raise ImportError("pyvista is required for mesh_from_pyvista")

    coordinates = np.asarray(grid.points, dtype=np.float64)
    connectivity = np.asarray(grid.cell_connectivity, dtype=np.int32)
    offsets = np.asarray(grid.offset, dtype=np.int32)
    cell_types_np = np.asarray(grid.celltypes, dtype=np.int32)

    if tdim is None:
        first_type = int(cell_types_np[0])
        tdim = _VTK_TDIM.get(first_type)
        if tdim is None:
            raise ValueError(
                f"Cannot infer tdim from VTK cell type {first_type}. "
                "Pass tdim explicitly."
            )

    return MeshView(
        coordinates=coordinates,
        connectivity=connectivity,
        offsets=offsets,
        cell_types=cell_types_np,
        tdim=tdim,
    )


def safe_part_name(expr):
    """Return a filesystem-safe compact name for a level-set selection."""
    return (
        str(expr)
        .replace(" ", "")
        .replace("<", "lt")
        .replace(">", "gt")
        .replace("=", "eq")
    )


def cutmesh_to_pyvista(cut_mesh):
    """Convert a cutcells CutMesh to a pyvista.UnstructuredGrid."""
    if pv is None:
        raise ImportError("pyvista is required for cutmesh_to_pyvista")

    points = np.asarray(cut_mesh.vertex_coords, dtype=np.float64)
    if points.ndim == 2 and points.shape[1] == 2:
        points = np.column_stack([points, np.zeros(points.shape[0])])

    return pv.UnstructuredGrid(
        np.asarray(cut_mesh.cells, dtype=np.int64),
        np.asarray(cut_mesh.vtk_types, dtype=np.uint8),
        points,
    )


# cell::type → VTK cell type: point, interval, triangle, tetrahedron,
# quadrilateral, hexahedron, prism, pyramid.
_CELL_TYPE_TO_VTK = np.array([1, 3, 5, 10, 9, 12, 13, 14], dtype=np.uint8)


def to_pyvista(mesh):
    """Convert a CutCells mesh to a pyvista.UnstructuredGrid.

    ``mesh`` is a ``MeshView`` (the background mesh), a ``CutMesh`` (the
    straight pieces of the lookup tables) or a ``QuadraysLeafMesh`` (the
    curved Lagrange pieces of quadrays, as from
    ``part.visualization_mesh(...)``). Points in 2D get ``z = 0``. A
    ``QuadraysLeafMesh`` keeps its parent cells as the cell array
    ``"parent"``.
    """
    if pv is None:
        raise ImportError("pyvista is required for to_pyvista")

    if hasattr(mesh, "vertex_coords"):
        return cutmesh_to_pyvista(mesh)

    if hasattr(mesh, "vtk_types"):
        points = np.asarray(mesh.points, dtype=np.float64)
        vtk_types = np.asarray(mesh.vtk_types, dtype=np.uint8)
    else:
        points = np.asarray(mesh.coordinates, dtype=np.float64)
        vtk_types = _CELL_TYPE_TO_VTK[np.asarray(mesh.cell_types, dtype=np.int64)]
    if points.shape[1] < 3:
        points = np.column_stack([points, np.zeros((points.shape[0], 3 - points.shape[1]))])

    connectivity = np.asarray(mesh.connectivity, dtype=np.int64)
    offsets = np.asarray(mesh.offsets, dtype=np.int64)
    sizes = np.diff(offsets)
    cells = np.insert(connectivity, offsets[:-1], sizes)

    grid = pv.UnstructuredGrid(cells, vtk_types, points)
    if hasattr(mesh, "parent"):
        grid.cell_data["parent"] = np.asarray(mesh.parent)
    return grid


def structured_triangle_mesh_view(x0, y0, x1, y1, nx, ny):
    """Pure NumPy structured triangular MeshView on [x0,x1] x [y0,y1].

    Each quadrilateral grid cell is split along the bottom-left to top-right
    diagonal. This helper does not require pyvista.
    """
    from . import MeshView

    x = np.linspace(x0, x1, num=nx)
    y = np.linspace(y0, y1, num=ny)
    xx, yy = np.meshgrid(x, y)
    coordinates = np.column_stack([xx.ravel(), yy.ravel()]).astype(np.float64)

    connectivity = []
    cell_types = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            v0 = j * nx + i
            v1 = v0 + 1
            v3 = v0 + nx
            v2 = v3 + 1
            connectivity.extend([v0, v1, v2])
            connectivity.extend([v0, v2, v3])
            cell_types.extend([5, 5])  # VTK_TRIANGLE

    connectivity = np.asarray(connectivity, dtype=np.int32)
    offsets = np.arange(0, connectivity.size + 1, 3, dtype=np.int32)
    cell_types = np.asarray(cell_types, dtype=np.int32)

    return MeshView(coordinates, connectivity, offsets, cell_types, tdim=2)


def _grid_points(lower, upper, n):
    axes = [np.linspace(a, b, num=k) for a, b, k in zip(lower, upper, n)]
    grids = np.meshgrid(*axes, indexing="ij")
    # Node (i, j[, k]) has index i + nx * (j + ny * k): x runs fastest.
    return np.column_stack([g.transpose().ravel() for g in grids]).astype(np.float64)


def _grid_corners(n):
    """Per grid cell, the node indices of its corners (x fastest)."""
    shape = [k - 1 for k in n]
    idx = np.indices(shape[::-1]).reshape(len(n), -1)[::-1]
    strides = np.cumprod([1] + list(n[:-1]))
    base = sum(i * s for i, s in zip(idx, strides))
    corners = {}
    for bits in np.ndindex(*(2,) * len(n)):
        corners[bits[::-1]] = base + sum(b * s for b, s in zip(bits[::-1], strides))
    return corners


def _mesh_view(coordinates, cells, vtk_type, tdim):
    from . import MeshView

    cells = np.asarray(cells, dtype=np.int32)
    connectivity = np.ascontiguousarray(cells.ravel())
    offsets = np.arange(0, connectivity.size + 1, cells.shape[1], dtype=np.int32)
    cell_types = np.full(cells.shape[0], vtk_type, dtype=np.int32)
    return MeshView(coordinates, connectivity, offsets, cell_types, tdim=tdim)


def structured_quadrilateral_mesh_view(x0, y0, x1, y1, nx, ny):
    """Pure NumPy structured quadrilateral MeshView on [x0,x1] x [y0,y1].

    ``nx`` and ``ny`` are numbers of grid points, so the mesh has
    ``(nx - 1) * (ny - 1)`` quadrilaterals. This helper does not require
    pyvista.
    """
    n = (nx, ny)
    c = _grid_corners(n)
    cells = np.column_stack([c[0, 0], c[1, 0], c[1, 1], c[0, 1]])
    return _mesh_view(_grid_points((x0, y0), (x1, y1), n), cells, 9, 2)  # VTK_QUAD


def structured_tetrahedron_mesh_view(x0, y0, z0, x1, y1, z1, nx, ny, nz):
    """Pure NumPy structured tetrahedral MeshView on a box.

    Each grid hexahedron is split into six tetrahedra around its diagonal
    from the lowest to the highest corner, so the mesh is conforming and has
    ``6 * (nx - 1) * (ny - 1) * (nz - 1)`` tetrahedra, all positively
    oriented. This helper does not require pyvista.
    """
    n = (nx, ny, nz)
    c = _grid_corners(n)
    cells = []
    for order in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)):
        corner = [0, 0, 0]
        tet = [c[tuple(corner)]]
        for axis in order:
            corner[axis] = 1
            tet.append(c[tuple(corner)])
        # Odd permutations of the axes give negatively oriented tetrahedra.
        if order in ((0, 2, 1), (1, 0, 2), (2, 1, 0)):
            tet[1], tet[2] = tet[2], tet[1]
        cells.append(np.column_stack(tet))
    cells = np.stack(cells, axis=1).reshape(-1, 4)
    return _mesh_view(_grid_points((x0, y0, z0), (x1, y1, z1), n), cells, 10, 3)  # VTK_TETRA


def structured_hexahedron_mesh_view(x0, y0, z0, x1, y1, z1, nx, ny, nz):
    """Pure NumPy structured hexahedral MeshView on a box.

    ``nx``, ``ny`` and ``nz`` are numbers of grid points, so the mesh has
    ``(nx - 1) * (ny - 1) * (nz - 1)`` hexahedra in VTK node order. This
    helper does not require pyvista.
    """
    n = (nx, ny, nz)
    c = _grid_corners(n)
    cells = np.column_stack([
        c[0, 0, 0], c[1, 0, 0], c[1, 1, 0], c[0, 1, 0],
        c[0, 0, 1], c[1, 0, 1], c[1, 1, 1], c[0, 1, 1],
    ])
    return _mesh_view(_grid_points((x0, y0, z0), (x1, y1, z1), n), cells, 12, 3)  # VTK_HEXAHEDRON


# ---------------------------------------------------------------------------
# 2D structured meshes
# ---------------------------------------------------------------------------


def rectangle_triangle_mesh(x0, y0, x1, y1, nx, ny):
    """Structured triangular mesh on [x0,x1] × [y0,y1].

    Creates a regular grid of nx × ny points and triangulates it with
    Delaunay. Returns a pyvista UnstructuredGrid.
    """
    if pv is None:
        raise ImportError("pyvista is required")
    x = np.linspace(x0, x1, num=nx)
    y = np.linspace(y0, y1, num=ny)
    xx, yy, zz = np.meshgrid(x, y, [0.0])
    points = np.c_[xx.ravel(), yy.ravel(), zz.ravel()]
    poly = pv.PolyData(points).delaunay_2d()
    return pv.UnstructuredGrid(poly)


def rectangle_quad_mesh(x0, y0, x1, y1, nx, ny):
    """Structured quadrilateral mesh on [x0,x1] × [y0,y1].

    Returns a pyvista UnstructuredGrid with (nx-1)*(ny-1) quads.
    """
    if pv is None:
        raise ImportError("pyvista is required")
    x = np.linspace(x0, x1, num=nx)
    y = np.linspace(y0, y1, num=ny)
    xx, yy = np.meshgrid(x, y)
    points = np.c_[xx.ravel(), yy.ravel(), np.zeros(nx * ny)]

    cells = []
    celltypes = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            v0 = j * nx + i
            v1 = v0 + 1
            v2 = v1 + nx
            v3 = v0 + nx
            cells.extend([4, v0, v1, v2, v3])
            celltypes.append(9)  # VTK_QUAD

    return pv.UnstructuredGrid(
        np.array(cells, dtype=np.int64),
        np.array(celltypes, dtype=np.uint8),
        points,
    )


# ---------------------------------------------------------------------------
# 3D structured meshes
# ---------------------------------------------------------------------------


def box_tetrahedron_mesh(x0, y0, z0, x1, y1, z1, nx, ny, nz):
    """Structured tetrahedral mesh on [x0,x1] × [y0,y1] × [z0,z1].

    Creates a regular grid and uses Delaunay tetrahedralization.
    Returns a pyvista UnstructuredGrid.
    """
    if pv is None:
        raise ImportError("pyvista is required")
    x = np.linspace(x0, x1, num=nx)
    y = np.linspace(y0, y1, num=ny)
    z = np.linspace(z0, z1, num=nz)
    xx, yy, zz = np.meshgrid(x, y, z)
    points = np.c_[xx.ravel(), yy.ravel(), zz.ravel()]
    poly = pv.PolyData(points).delaunay_3d()
    return pv.UnstructuredGrid(poly)


def box_hex_mesh(x0, y0, z0, x1, y1, z1, nx, ny, nz):
    """Structured hexahedral mesh on [x0,x1] × [y0,y1] × [z0,z1].

    Returns a pyvista UnstructuredGrid with (nx-1)*(ny-1)*(nz-1) hexahedra.
    """
    if pv is None:
        raise ImportError("pyvista is required")
    sg = pv.ImageData(
        dimensions=(nx, ny, nz),
        origin=(x0, y0, z0),
        spacing=(
            (x1 - x0) / (nx - 1),
            (y1 - y0) / (ny - 1),
            (z1 - z0) / (nz - 1),
        ),
    )
    return sg.cast_to_unstructured_grid()
