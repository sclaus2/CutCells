"""
Quadrature point visualisation — 2D (triangulated mesh, z = 0)
==============================================================
Creates a background triangular mesh in the square [-1, 1]^2 embedded in
the z = 0 plane and cuts it with a circular level-set phi(x) = ||x|| - r.

Two quadrature integrals are computed and visualised:

  Area integral  (phi < 0):  ∫_{Ω_h} 1 dA  ≈  π r²
    · fully-inside triangles  → standard reference Gauss rule
    · cut triangles           → cut-cell sub-triangulation rule
  Perimeter integral (phi = 0):  ∫_{Γ_h} 1 ds  ≈  2π r
    · cut triangles           → interface edge rule

Both sets of quadrature points are displayed. Physical points are obtained
by mapping 2D reference coordinates through the affine triangle map to the
3D z = 0 plane.
"""

import numpy as np
import pyvista as pv
import cutcells

RADIUS = 0.7
N = 22  # grid resolution (NxN before Delaunay)
ORDER = 2  # quadrature order


# ---------------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------------
def circle_ls(pts: np.ndarray) -> np.ndarray:
    return np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2) - RADIUS


def build_mesh(N: int) -> pv.UnstructuredGrid:
    x = np.linspace(-1.0, 1.0, N)
    y = np.linspace(-1.0, 1.0, N)
    xx, yy, zz = np.meshgrid(x, y, [0.0])
    pts = np.c_[xx.ravel(), yy.ravel(), zz.ravel()]
    return pv.UnstructuredGrid(pv.PolyData(pts).delaunay_2d())


grid = build_mesh(N)
points = np.asarray(grid.points, dtype=np.float64)
ls_values = circle_ls(points)
points_flat = points.ravel()
connectivity = np.asarray(grid.cell_connectivity, dtype=np.int32)
offset = np.asarray(grid.offset, dtype=np.int32)
celltypes = np.asarray(grid.celltypes, dtype=np.int32)

# ---------------------------------------------------------------------------
# Area quadrature  (phi < 0  — inside + cut cells)
# ---------------------------------------------------------------------------
rules_area = cutcells.runtime_quadrature(
    ls_values,
    points_flat,
    connectivity,
    offset,
    celltypes,
    "phi<0",
    True,
    ORDER,
)
phys_area = cutcells.physical_points(
    rules_area, points_flat, connectivity, offset, celltypes
).reshape(-1, 3)
w_area = np.asarray(rules_area.weights)
n_area = len(np.asarray(rules_area.parent_map))

# ---------------------------------------------------------------------------
# Perimeter quadrature  (phi = 0  — interface)
# ---------------------------------------------------------------------------
rules_peri = cutcells.runtime_quadrature(
    ls_values,
    points_flat,
    connectivity,
    offset,
    celltypes,
    "phi=0",
    True,
    ORDER,
)
phys_peri = cutcells.physical_points(
    rules_peri, points_flat, connectivity, offset, celltypes
).reshape(-1, 3)
w_peri = np.asarray(rules_peri.weights)
n_peri = len(np.asarray(rules_peri.parent_map))

exact_area = np.pi * RADIUS**2
exact_peri = 2.0 * np.pi * RADIUS
print(f"Mesh : {grid.n_cells} cells  ({N}×{N} Delaunay 2D)")
print(
    f"Area : {n_area:4d} rules, {len(w_area):5d} qpts  "
    f"∫1 dA = {w_area.sum():.6f}  (π r² = {exact_area:.6f}  "
    f"err = {abs(w_area.sum() - exact_area):.4f})"
)
print(
    f"Peri : {n_peri:4d} rules, {len(w_peri):5d} qpts  "
    f"∫1 ds = {w_peri.sum():.6f}  (2π r  = {exact_peri:.6f}  "
    f"err = {abs(w_peri.sum() - exact_peri):.4f})"
)

# ---------------------------------------------------------------------------
# Classify cells and split area q-pts by parent cell type
# ---------------------------------------------------------------------------
inside_ids = cutcells.locate_cells(
    ls_values, points_flat, connectivity, offset, celltypes, "phi<0"
)
inside_set = set(inside_ids.tolist())

n_cells = grid.n_cells
cell_min_ls = np.zeros(n_cells)
cell_max_ls = np.zeros(n_cells)
for ci in range(n_cells):
    vids = connectivity[offset[ci] : offset[ci + 1]]
    cell_min_ls[ci] = ls_values[vids].min()
    cell_max_ls[ci] = ls_values[vids].max()
cut_ids = np.where((cell_min_ls < 0) & (cell_max_ls > 0))[0]

# Split area quadrature points into inside-cell pts and cut-cell pts so the
# two populations can be shown in distinct colours (not via a weight colormap
# which confusingly mixes cell-type information with weight magnitude).
area_parent_map = np.asarray(rules_area.parent_map)
area_offset = np.asarray(rules_area.offset)
inside_mask = np.zeros(len(w_area), dtype=bool)
for rule_i, parent in enumerate(area_parent_map):
    q0, q1 = int(area_offset[rule_i]), int(area_offset[rule_i + 1])
    if parent in inside_set:
        inside_mask[q0:q1] = True

phys_area_inside = phys_area[inside_mask]  # full cells → reference Gauss rule
phys_area_cut = phys_area[~inside_mask]  # cut cells → sub-triangulation rule

# ---------------------------------------------------------------------------
# Build cut sub-meshes for visual verification
# ---------------------------------------------------------------------------
# Volume sub-triangulation of cut region (phi < 0 part of each cut triangle)
cut_vol = cutcells.cut_vtk_mesh(
    ls_values, points_flat, connectivity, offset, celltypes, "phi<0"
)
pv_cut_vol = pv.UnstructuredGrid(
    np.asarray(cut_vol.cells),
    np.asarray(cut_vol.vtk_types),
    np.asarray(cut_vol.vertex_coords),
)

# Interface edges (phi = 0 — one interval per cut triangle)
cut_iface = cutcells.cut_vtk_mesh(
    ls_values, points_flat, connectivity, offset, celltypes, "phi=0"
)
pv_cut_iface = pv.UnstructuredGrid(
    np.asarray(cut_iface.cells),
    np.asarray(cut_iface.vtk_types),
    np.asarray(cut_iface.vertex_coords),
)

# ---------------------------------------------------------------------------
# PyVista visualisation
# ---------------------------------------------------------------------------
POINT_SIZE = 8  # pixels, same for all point clouds → easy visual comparison

pl = pv.Plotter(window_size=(1200, 900))
pl.set_background("white")

# Background wireframe
pl.add_mesh(grid, style="wireframe", color="lightgrey", line_width=0.5, opacity=0.6)

# Inside cells  (blue fill)
if len(inside_ids):
    pl.add_mesh(
        grid.extract_cells(inside_ids),
        color="#80b3ff",
        opacity=0.25,
        show_edges=True,
        edge_color="steelblue",
        line_width=0.5,
        label="inside cells",
    )

# Cut sub-triangles of the volume part (phi<0)  — dark green outline, transparent fill
pl.add_mesh(
    pv_cut_vol,
    color="#00aa44",
    opacity=0.20,
    show_edges=True,
    edge_color="#007730",
    line_width=1.2,
    label="cut sub-tris (phi<0)",
)

# Interface edges (phi=0)  — thick magenta lines + their nodes
pl.add_mesh(
    pv_cut_iface,
    color="magenta",
    line_width=3.0,
    label="interface edges (phi=0)",
)
pl.add_points(
    pv_cut_iface.points,
    color="darkviolet",
    point_size=10,
    render_points_as_spheres=True,
    label=f"interface nodes ({pv_cut_iface.n_points})",
)

# Cut volume sub-triangle nodes
pl.add_points(
    pv_cut_vol.points,
    color="#007730",
    point_size=7,
    render_points_as_spheres=True,
    label=f"cut vol nodes ({pv_cut_vol.n_points})",
)

# Area q-pts in fully-inside cells  (green)
if len(phys_area_inside):
    pl.add_points(
        phys_area_inside,
        color="limegreen",
        point_size=POINT_SIZE,
        render_points_as_spheres=True,
        label=f"area q-pts inside ({len(phys_area_inside)})",
    )

# Area q-pts in cut cells  (orange-red)
if len(phys_area_cut):
    pl.add_points(
        phys_area_cut,
        color="orangered",
        point_size=POINT_SIZE,
        render_points_as_spheres=True,
        label=f"area q-pts cut ({len(phys_area_cut)})",
    )

# Perimeter q-pts  (dark blue)
if len(phys_peri):
    pl.add_points(
        phys_peri,
        color="navy",
        point_size=POINT_SIZE,
        render_points_as_spheres=True,
        label=f"perimeter q-pts ({len(phys_peri)})",
    )

# Circle outline
theta = np.linspace(0.0, 2.0 * np.pi, 512)
circle_pts = np.c_[RADIUS * np.cos(theta), RADIUS * np.sin(theta), np.zeros(512)]
pl.add_mesh(
    pv.Spline(circle_pts, 512), color="black", line_width=2.0, label="level set φ = 0"
)

pl.add_legend(bcolor="white", border=True, size=(0.42, 0.34))
pl.camera_position = "xy"
pl.add_title(
    f"Quadrature — 2D  r={RADIUS}  order={ORDER}  "
    f"∫dA≈{w_area.sum():.4f} [π r²={exact_area:.4f}]  "
    f"∫ds≈{w_peri.sum():.4f} [2πr={exact_peri:.4f}]",
    font_size=9,
)
pl.show()
