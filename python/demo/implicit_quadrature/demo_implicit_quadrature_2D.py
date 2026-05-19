"""
Implicit quadrature demo on the reference triangle.

The default level set is the analytic circle

    phi(x, y) = x^2 + y^2 - R^2

on K = {(x, y): x >= 0, y >= 0, x + y <= 1}.  For R < 1/sqrt(2) the
negative domain is a quarter disk and the interface is a quarter circle, so
the exact references are pi*R^2/4 and pi*R/2.

The `--case height-recovery` mode uses phi(x, y) = x - a and deliberately
sets the sign-oriented chart transversality threshold too high, so the
projected height-direction recovery path is required.

If PyVista is available, this writes:

    results/implicit_quadrature_triangle_points.vtp

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


def make_reference_triangle():
    coordinates = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        dtype=np.float64,
    )
    connectivity = np.array([0, 1, 2], dtype=np.int32)
    offsets = np.array([0, 3], dtype=np.int32)
    cell_types = np.array([5], dtype=np.int32)
    return cutcells.MeshView(coordinates, connectivity, offsets, cell_types, tdim=2)


def physical_points(mesh, rules):
    vtk_type = np.full(
        mesh.num_cells(), 5 if mesh.tdim == 2 else 10, dtype=np.int32
    )
    return cutcells.physical_points(
        rules,
        np.asarray(mesh.coordinates, dtype=np.float64).ravel(),
        np.asarray(mesh.connectivity, dtype=np.int32),
        np.asarray(mesh.offsets, dtype=np.int32),
        vtk_type,
    ).reshape(-1, 3)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case",
        choices=("circle", "height-recovery"),
        default="circle",
    )
    parser.add_argument("--radius", type=float, default=0.45)
    parser.add_argument("--cut-at", type=float, default=0.2)
    parser.add_argument("--order", type=int, default=8)
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    args = parser.parse_args()

    mesh = make_reference_triangle()
    opts = cutcells.ImplicitQuadratureOptions()

    if args.case == "circle":
        radius = args.radius
        ls = cutcells.interpolate_level_set(
            mesh,
            lambda X: X[:, 0] ** 2 + X[:, 1] ** 2 - radius**2,
            degree=2,
            name="phi",
        )
        exact_negative = math.pi * radius**2 / 4.0
        exact_positive = 0.5 - exact_negative
        exact_interface = math.pi * radius / 2.0
        title = "Implicit quadrature on a circular cut triangle"
    else:
        cut_at = args.cut_at
        ls = cutcells.interpolate_level_set(
            mesh,
            lambda X: X[:, 0] - cut_at,
            degree=1,
            name="phi",
        )
        opts.min_transversality = 1.0
        opts.min_height_transversality = 0.05
        opts.enable_height_direction_recovery = True
        exact_negative = cut_at * (1.0 - cut_at) + 0.5 * cut_at**2
        exact_positive = 0.5 - exact_negative
        exact_interface = 1.0 - cut_at
        title = "Projected height-direction recovery on a linear cut triangle"

    result = cutcells.cut(mesh, ls)

    q_negative = result["phi < 0"].implicit_quadrature(
        order=args.order, mode="cut_only", options=opts
    )
    q_positive = result["phi > 0"].implicit_quadrature(
        order=args.order, mode="cut_only", options=opts
    )
    q_interface = result["phi = 0"].implicit_quadrature(
        order=args.order, mode="cut_only", options=opts
    )

    negative = float(np.asarray(q_negative.weights).sum())
    positive = float(np.asarray(q_positive.weights).sum())
    interface = float(np.asarray(q_interface.weights).sum())
    pts_interface_ref = np.asarray(q_interface.points).reshape(-1, q_interface.tdim)
    if args.case == "circle":
        radius = args.radius
        residual = np.max(
            np.abs(
                pts_interface_ref[:, 0] ** 2
                + pts_interface_ref[:, 1] ** 2
                - radius**2
            )
        )
    else:
        residual = np.max(np.abs(pts_interface_ref[:, 0] - args.cut_at))

    print(f"case                   : {args.case}")
    print(f"order                  : {args.order}")
    print(f"negative area          : {negative:.12e}")
    print(f"negative exact/error   : {exact_negative:.12e} / {abs(negative - exact_negative):.3e}")
    print(f"positive area          : {positive:.12e}")
    print(f"positive exact/error   : {exact_positive:.12e} / {abs(positive - exact_positive):.3e}")
    print(f"interface length       : {interface:.12e}")
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
        cloud = pv.PolyData(points)
        cloud["weight"] = weights
        cloud["part_id"] = np.full(points.shape[0], part_id, dtype=np.int32)
        clouds.append(cloud)

    merged = clouds[0]
    for cloud in clouds[1:]:
        merged = merged.merge(cloud)
    output_path = args.output_dir / f"implicit_quadrature_triangle_{args.case}_points.vtp"
    merged.save(output_path)
    print(f"wrote {output_path}")

    if args.no_show:
        return

    triangle = pv.PolyData(
        np.asarray(mesh.coordinates, dtype=np.float64),
        np.array([3, 0, 1, 2], dtype=np.int64),
    )
    if args.case == "circle":
        radius = args.radius
        theta = np.linspace(0.0, math.pi / 2.0, 160)
        interface_curve = np.column_stack(
            [radius * np.cos(theta), radius * np.sin(theta), np.zeros_like(theta)]
        )
    else:
        cut_at = args.cut_at
        interface_curve = np.array(
            [[cut_at, 0.0, 0.0], [cut_at, 1.0 - cut_at, 0.0]],
            dtype=np.float64,
        )

    plotter = pv.Plotter()
    plotter.set_background("white")
    plotter.add_mesh(triangle, style="wireframe", color="black", line_width=2)
    plotter.add_points(pts_negative, color="royalblue", point_size=10, render_points_as_spheres=True)
    plotter.add_points(pts_positive, color="lightgrey", point_size=8, render_points_as_spheres=True)
    plotter.add_points(pts_interface, color="crimson", point_size=14, render_points_as_spheres=True)
    plotter.add_lines(interface_curve, color="crimson", width=3)
    plotter.add_title(title, font_size=10)
    plotter.show()


if __name__ == "__main__":
    main()
