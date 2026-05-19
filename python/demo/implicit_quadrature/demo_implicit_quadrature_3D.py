"""
Implicit quadrature demo on the reference tetrahedron.

The level set is the analytic sphere

    phi(x, y, z) = x^2 + y^2 + z^2 - R^2

on K = {(x,y,z): x,y,z >= 0, x + y + z <= 1}. For R < 1/sqrt(3),
the negative domain is an octant ball and the interface is an octant sphere:

    |phi < 0| = pi R^3 / 6,
    |phi = 0| = pi R^2 / 2.

If PyVista is available, this writes

    results/implicit_quadrature_tetra_sphere_points.vtp

and opens an interactive scene unless --no-show is passed.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

import cutcells

try:
    import pyvista as pv
except ImportError:  # pragma: no cover - optional demo dependency
    pv = None


def make_reference_tetrahedron():
    coordinates = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )
    connectivity = np.array([0, 1, 2, 3], dtype=np.int32)
    offsets = np.array([0, 4], dtype=np.int32)
    cell_types = np.array([10], dtype=np.int32)
    return cutcells.MeshView(coordinates, connectivity, offsets, cell_types, tdim=3)


def physical_points(mesh, rules):
    vtk_type = np.full(mesh.num_cells(), 10, dtype=np.int32)
    return cutcells.physical_points(
        rules,
        np.asarray(mesh.coordinates, dtype=np.float64).ravel(),
        np.asarray(mesh.connectivity, dtype=np.int32),
        np.asarray(mesh.offsets, dtype=np.int32),
        vtk_type,
    ).reshape(-1, 3)


def octant_sphere_patch(radius: float, n_theta: int = 48, n_phi: int = 32):
    if pv is None:
        return None

    theta = np.linspace(0.0, math.pi / 2.0, n_theta)
    phi = np.linspace(0.0, math.pi / 2.0, n_phi)
    points = []
    for p in phi:
        sp = math.sin(p)
        cp = math.cos(p)
        for t in theta:
            points.append(
                [
                    radius * sp * math.cos(t),
                    radius * sp * math.sin(t),
                    radius * cp,
                ]
            )

    faces = []
    for i in range(n_phi - 1):
        for j in range(n_theta - 1):
            a = i * n_theta + j
            b = (i + 1) * n_theta + j
            c = (i + 1) * n_theta + j + 1
            d = i * n_theta + j + 1
            faces.extend([4, a, b, c, d])

    return pv.PolyData(np.asarray(points, dtype=np.float64), np.asarray(faces))


def make_tetra_grid(coordinates):
    cells = np.array([4, 0, 1, 2, 3], dtype=np.int64)
    cell_types = np.array([10], dtype=np.uint8)
    return pv.UnstructuredGrid(cells, cell_types, coordinates)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--radius", type=float, default=0.45)
    parser.add_argument("--order", type=int, default=8)
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument("--screenshot", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    args = parser.parse_args()

    mesh = make_reference_tetrahedron()
    radius = args.radius
    ls = cutcells.interpolate_level_set(
        mesh,
        lambda X: X[:, 0] ** 2 + X[:, 1] ** 2 + X[:, 2] ** 2 - radius**2,
        degree=2,
        name="phi",
    )
    result = cutcells.cut(mesh, ls)

    q_negative = result["phi < 0"].implicit_quadrature(
        order=args.order, mode="cut_only"
    )
    q_positive = result["phi > 0"].implicit_quadrature(
        order=args.order, mode="cut_only"
    )
    q_interface = result["phi = 0"].implicit_quadrature(
        order=args.order, mode="cut_only"
    )

    negative = float(np.asarray(q_negative.weights).sum())
    positive = float(np.asarray(q_positive.weights).sum())
    interface = float(np.asarray(q_interface.weights).sum())
    exact_negative = math.pi * radius**3 / 6.0
    exact_positive = 1.0 / 6.0 - exact_negative
    exact_interface = math.pi * radius**2 / 2.0

    pts_interface_ref = np.asarray(q_interface.points).reshape(-1, q_interface.tdim)
    residual = np.max(np.abs(np.sum(pts_interface_ref**2, axis=1) - radius**2))

    print(f"order                  : {args.order}")
    print(f"negative volume        : {negative:.12e}")
    print(f"negative exact/error   : {exact_negative:.12e} / {abs(negative - exact_negative):.3e}")
    print(f"positive volume        : {positive:.12e}")
    print(f"positive exact/error   : {exact_positive:.12e} / {abs(positive - exact_positive):.3e}")
    print(f"interface area         : {interface:.12e}")
    print(f"interface exact/error  : {exact_interface:.12e} / {abs(interface - exact_interface):.3e}")
    print(f"max interface residual : {residual:.3e}")

    if pv is None:
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pts_negative = physical_points(mesh, q_negative)
    pts_positive = physical_points(mesh, q_positive)
    pts_interface = physical_points(mesh, q_interface)

    clouds = []
    for part_id, points, weights in [
        (0, pts_negative, np.asarray(q_negative.weights).copy()),
        (1, pts_positive, np.asarray(q_positive.weights).copy()),
        (2, pts_interface, np.asarray(q_interface.weights).copy()),
    ]:
        if len(points) == 0:
            continue
        cloud = pv.PolyData(points)
        cloud["weight"] = weights
        cloud["part_id"] = np.full(points.shape[0], part_id, dtype=np.int32)
        clouds.append(cloud)

    if clouds:
        merged = clouds[0]
        for cloud in clouds[1:]:
            merged = merged.merge(cloud)
        output_path = args.output_dir / "implicit_quadrature_tetra_sphere_points.vtp"
        merged.save(output_path)
        print(f"wrote {output_path}")

    if args.no_show and args.screenshot is None:
        return

    tetra = make_tetra_grid(np.asarray(mesh.coordinates, dtype=np.float64))
    sphere_patch = octant_sphere_patch(radius)
    off_screen = args.no_show or args.screenshot is not None
    plotter = pv.Plotter(window_size=(1300, 1000), off_screen=off_screen)
    plotter.set_background("white")
    plotter.add_mesh(
        tetra,
        color="white",
        opacity=0.12,
        show_edges=True,
        edge_color="black",
        line_width=2,
        label="reference tetrahedron",
    )
    if sphere_patch is not None:
        plotter.add_mesh(
            sphere_patch,
            color="lightyellow",
            opacity=0.42,
            smooth_shading=True,
            label="exact octant sphere",
        )
        plotter.add_mesh(
            sphere_patch.extract_feature_edges(
                boundary_edges=True,
                feature_edges=False,
                manifold_edges=False,
                non_manifold_edges=False,
            ),
            color="gold",
            line_width=3,
        )

    if len(pts_negative):
        plotter.add_points(
            pts_negative,
            color="royalblue",
            point_size=8,
            render_points_as_spheres=True,
            label="negative volume points",
        )
    if len(pts_positive):
        plotter.add_points(
            pts_positive,
            color="lightgrey",
            point_size=6,
            render_points_as_spheres=True,
            label="positive volume points",
        )
    if len(pts_interface):
        plotter.add_points(
            pts_interface,
            color="crimson",
            point_size=12,
            render_points_as_spheres=True,
            label="interface points",
        )

    plotter.add_legend(bcolor="white", border=True, loc="upper right")
    plotter.add_title(
        f"Implicit quadrature on a spherical cut tetrahedron, order {args.order}",
        font_size=10,
    )
    plotter.camera_position = [
        (1.75, 1.65, 1.35),
        (0.22, 0.22, 0.22),
        (0.0, 0.0, 1.0),
    ]
    if args.screenshot is not None:
        args.screenshot.parent.mkdir(parents=True, exist_ok=True)
        plotter.show(screenshot=str(args.screenshot))
        print(f"wrote {args.screenshot}")
    else:
        plotter.show()


if __name__ == "__main__":
    main()
