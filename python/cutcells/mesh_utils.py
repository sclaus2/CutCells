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
