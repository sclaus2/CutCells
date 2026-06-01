"""
Quadrature point visualisation — 3D (tetrahedral mesh)
======================================================
Creates a background tetrahedral mesh inside the box [-1, 1]^3 and cuts it
with a spherical level-set phi(x) = ||x|| - r.

For every qualifying cell runtime_quadrature produces reference-space
quadrature rules; physical_points maps them to the physical domain.

The plot shows:
  - clipped background mesh (light grey wireframe)
  - cut cells (red, semi-transparent)  — one half of the sphere, clipped
  - quadrature points (glyphs sized and coloured by weight)
  - a translucent sphere surface at the interface
"""

import numpy as np
import pyvista as pv
import cutcells

RADIUS = 0.65
N = 9  # grid resolution
ORDER = 3  # quadrature order
CUT_TYPE = "phi<0"


# ---------------------------------------------------------------------------
# Build mesh
# ---------------------------------------------------------------------------


def sphere_ls(pts: np.ndarray) -> np.ndarray:
    """Signed distance to a sphere of radius RADIUS centred at the origin."""
    return np.sqrt((pts**2).sum(axis=1)) - RADIUS


def build_box_mesh(N: int) -> pv.UnstructuredGrid:
    x = np.linspace(-1.0, 1.0, N)
    y = np.linspace(-1.0, 1.0, N)
    z = np.linspace(-1.0, 1.0, N)
    xx, yy, zz = np.meshgrid(x, y, z)
    pts = np.c_[xx.ravel(), yy.ravel(), zz.ravel()]
    grid = pv.UnstructuredGrid(pv.PolyData(pts).delaunay_3d())
    return grid


grid = build_box_mesh(N)

points = np.asarray(grid.points, dtype=np.float64)
ls_values = sphere_ls(points)
points_flat = points.ravel()
connectivity = np.asarray(grid.cell_connectivity, dtype=np.int32)
offset = np.asarray(grid.offset, dtype=np.int32)
celltypes = np.asarray(grid.celltypes, dtype=np.int32)

# ---------------------------------------------------------------------------
# Quadrature via runtime_quadrature
# ---------------------------------------------------------------------------

rules = cutcells.runtime_quadrature(
    ls_values, points_flat, connectivity, offset, celltypes, CUT_TYPE, True, ORDER
)

phys_pts = cutcells.physical_points(
    rules, points_flat, connectivity, offset, celltypes
).reshape(-1, 3)

weights = np.asarray(rules.weights)
parent_map = np.asarray(rules.parent_map)

exact = (4.0 / 3.0) * np.pi * RADIUS**3
print(f"Mesh: {grid.n_cells} cells  (N = {N})")
print(f"Quadrature rules: {len(parent_map)} rules,  {len(weights)} total points")
print(f"∫1 dV ≈ {weights.sum():.6f}  (exact (4/3)π r³ = {exact:.6f})")

# ---------------------------------------------------------------------------
# Classify cells for colouring
# ---------------------------------------------------------------------------

inside_ids = cutcells.locate_cells(
    ls_values, points_flat, connectivity, offset, celltypes, CUT_TYPE
)
# Cut cells: at least one vertex inside, at least one outside
n_cells = grid.n_cells
cell_min_ls = np.zeros(n_cells)
cell_max_ls = np.zeros(n_cells)
for ci in range(n_cells):
    vids = connectivity[offset[ci] : offset[ci + 1]]
    cell_min_ls[ci] = ls_values[vids].min()
    cell_max_ls[ci] = ls_values[vids].max()
cut_parent_ids = np.where((cell_min_ls < 0) & (cell_max_ls > 0))[0].astype(np.int64)

# ---------------------------------------------------------------------------
# PyVista visualisation
# ---------------------------------------------------------------------------

# Clipping plane (y = 0: show x > 0 half for interior view)
clip_normal = (0.0, -1.0, 0.0)
clip_origin = (0.0, 0.0, 0.0)

# Analytic sphere surface
sphere_surf = pv.Sphere(radius=RADIUS, theta_resolution=60, phi_resolution=60)

# Quadrature points on the selected half
mask = phys_pts[:, 1] >= 0.0
phys_half = phys_pts[mask]
weights_half = weights[mask]

qpts_cloud = pv.PolyData(phys_half)
qpts_cloud["weight"] = weights_half

glyph_geom = pv.Sphere(radius=1.0, theta_resolution=8, phi_resolution=8)
w_max = weights_half.max() if len(weights_half) else 1.0
qpts_cloud["radius"] = 0.008 + 0.018 * (weights_half / w_max)

pl = pv.Plotter(window_size=(1300, 1000))
pl.set_background("white")

# Clipped background mesh (wireframe)
grid_clipped = grid.clip(normal=clip_normal, origin=clip_origin, invert=True)
pl.add_mesh(
    grid_clipped, style="wireframe", color="#cccccc", line_width=0.4, opacity=0.5
)

# Inside cells (clipped half)
if len(inside_ids) > 0:
    inside_mesh = grid.extract_cells(inside_ids)
    inside_clipped = inside_mesh.clip(
        normal=clip_normal, origin=clip_origin, invert=True
    )
    pl.add_mesh(
        inside_clipped,
        color="#80b3ff",
        opacity=0.25,
        show_edges=True,
        edge_color="steelblue",
        line_width=0.6,
        label="inside cells",
    )

# Cut cells (clipped half)
if len(cut_parent_ids) > 0:
    cut_mesh = grid.extract_cells(cut_parent_ids)
    cut_clipped = cut_mesh.clip(normal=clip_normal, origin=clip_origin, invert=True)
    pl.add_mesh(
        cut_clipped,
        color="#ffaaaa",
        opacity=0.40,
        show_edges=True,
        edge_color="#cc4444",
        line_width=0.8,
        label="cut cells",
    )

# Sphere surface (translucent)
sphere_clipped = sphere_surf.clip(normal=clip_normal, origin=clip_origin, invert=True)
pl.add_mesh(
    sphere_clipped,
    color="lightyellow",
    opacity=0.35,
    show_edges=False,
    label="level set",
)
# Outline ring at the clip plane
sphere_ring = sphere_surf.clip(normal=clip_normal, origin=clip_origin)
pl.add_mesh(sphere_ring, style="wireframe", color="gold", line_width=1.2, opacity=0.6)

# Quadrature point glyphs (visible half only)
if len(phys_half) > 0:
    glyphs = qpts_cloud.glyph(geom=glyph_geom, scale="radius", orient=False)
    pl.add_mesh(
        glyphs,
        scalars="weight",
        cmap="plasma",
        show_scalar_bar=True,
        scalar_bar_args={
            "title": "weight",
            "n_labels": 4,
            "position_x": 0.82,
            "position_y": 0.25,
        },
        label="quadrature pts",
    )

pl.add_legend(bcolor="white", border=True, size=(0.25, 0.20), loc="upper right")
pl.add_title(
    f"Quadrature points — 3D  (r = {RADIUS}, order = {ORDER},  "
    f"∫1 dV ≈ {weights.sum():.4f}  [exact {exact:.4f}])",
    font_size=10,
)
pl.camera_position = [
    (2.5, 2.0, 1.8),
    (0.0, 0.0, 0.0),
    (0.0, 0.0, 1.0),
]
pl.show()
