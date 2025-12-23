import numpy as np
import pyvista as pv
import cutcells


import numpy as np


def level_set(xs, R0=1.0, a=0.25, k=6, center=(0.0, 0.0), sdf_like=True, eps=1e-12):
    """
    Flower-shaped level set.
    Negative inside, positive outside.

    Parameters
    ----------
    x, y : array_like
        Coordinates (scalars or numpy arrays).
    R0 : float
        Mean radius.
    a : float
        Petal amplitude (keep a < R0 to avoid self-intersections).
    k : int
        Number of petals.
    center : (float, float)
        Center of the flower.
    sdf_like : bool
        If True, applies a local scaling so that |grad(phi)| ~ 1 near the interface.
    eps : float
        Small number to avoid division by zero.

    Returns
    -------
    phi : ndarray
        Level set values (approx signed distance).
    """
    x = xs[0]
    y = xs[1]
    cx, cy = center
    dx = np.asarray(x) - cx
    dy = np.asarray(y) - cy

    r = np.hypot(dx, dy)
    theta = np.arctan2(dy, dx)

    R = R0 + a * np.cos(k * theta)
    phi = r - R  # simple (radial) signed distance approximation

    if not sdf_like:
        return phi

    # First-order correction: scale so that |∇phi| ≈ 1 near the boundary
    r_safe = np.maximum(r, eps)
    denom = np.maximum(dx * dx + dy * dy, eps)

    dtheta_dx = -dy / denom
    dtheta_dy = dx / denom

    dR_dtheta = -a * k * np.sin(k * theta)

    dphi_dx = dx / r_safe - dR_dtheta * dtheta_dx
    dphi_dy = dy / r_safe - dR_dtheta * dtheta_dy

    g = np.sqrt(dphi_dx * dphi_dx + dphi_dy * dphi_dy)
    g = np.maximum(g, eps)

    return phi / g  # more “signed-distance-like” near φ=0


def create_cut_mesh(grid):
    points = grid.points
    ls_values = np.zeros(len(points))
    j = 0
    for point in points:
        ls_values[j] = level_set(point)
        j = j + 1

    cut_mesh = cutcells.cut_vtk_mesh(
        ls_values, points, grid.cell_connectivity, grid.offset, grid.celltypes, "phi<0"
    )
    inside_cells = cutcells.locate_cells(
        ls_values, points, grid.cell_connectivity, grid.offset, grid.celltypes, "phi<0"
    )

    pv_cut = pv.UnstructuredGrid(
        cut_mesh.cells, cut_mesh.types, np.asarray(cut_mesh.vertex_coords)
    )
    extract = grid.extract_cells(inside_cells)

    return extract.merge(pv_cut)


# Mesh parameters
nx, ny = 21, 21
x = np.linspace(-1.5, 1.5, nx + 1)
y = np.linspace(-1.5, 1.5, ny + 1)
X, Y = np.meshgrid(x, y)
points = np.column_stack([X.ravel(), Y.ravel()])

cells = []
celltypes = []
vtk_tri = pv.CellType.TRIANGLE
vtk_quad = pv.CellType.QUAD

# Lower half: triangles, Upper half: quads
for j in range(ny):
    for i in range(nx):
        n0 = j * (nx + 1) + i
        n1 = n0 + 1
        n2 = n0 + (nx + 1)
        n3 = n2 + 1
        y_center = (Y[j, i] + Y[j + 1, i] + Y[j, i + 1] + Y[j + 1, i + 1]) / 4
        if y_center < 0:
            # Two triangles per cell
            cells.append([3, n0, n1, n2])
            celltypes.append(vtk_tri)
            cells.append([3, n1, n3, n2])
            celltypes.append(vtk_tri)
        else:
            # One quad per cell
            cells.append([4, n0, n1, n3, n2])
            celltypes.append(vtk_quad)

# Flatten cells for VTK
flat_cells = np.array([v for cell in cells for v in cell], dtype=np.int64)
celltypes = np.array(celltypes, dtype=np.uint8)

# Build a VTK mesh and cut it (ensures parent_map is populated).
points3 = np.c_[points, np.zeros((points.shape[0],), dtype=points.dtype)]
grid = pv.UnstructuredGrid(flat_cells, celltypes, points3)
mesh = create_cut_mesh(grid)

mesh.plot(cpos="xy", show_edges=True)

pl = pv.Plotter()
pl.add_mesh(grid, show_edges=True, color="white")
pl.add_mesh(mesh, show_edges=True)
pl.camera_position = "xy"
pl.show(screenshot="hybrid_flower.png")
