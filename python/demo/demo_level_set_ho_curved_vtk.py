# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier:    MIT

import argparse
from pathlib import Path

import numpy as np

import cutcells


def _parse_degrees(arg: str) -> list[int]:
    if arg == "all":
        return [2, 3, 4]
    return [int(arg)]


def _triangle_basix_nodes(degree: int) -> np.ndarray:
    points: list[tuple[float, float, float]] = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
    ]

    for k in range(1, degree):
        t = k / degree
        points.append((1.0 - t, t, 0.0))
    for k in range(1, degree):
        t = k / degree
        points.append((0.0, t, 0.0))
    for k in range(1, degree):
        t = k / degree
        points.append((t, 0.0, 0.0))

    for j in range(1, degree - 1):
        for i in range(1, degree - j):
            points.append((i / degree, j / degree, 0.0))

    return np.asarray(points, dtype=np.float64)


def _tetra_basix_nodes(degree: int) -> np.ndarray:
    points: list[tuple[float, float, float]] = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    ]

    for k in range(1, degree):
        t = k / degree
        points.append((0.0, 1.0 - t, t))
    for k in range(1, degree):
        t = k / degree
        points.append((1.0 - t, 0.0, t))
    for k in range(1, degree):
        t = k / degree
        points.append((1.0 - t, t, 0.0))
    for k in range(1, degree):
        t = k / degree
        points.append((0.0, 0.0, t))
    for k in range(1, degree):
        t = k / degree
        points.append((0.0, t, 0.0))
    for k in range(1, degree):
        t = k / degree
        points.append((t, 0.0, 0.0))

    faces = (
        ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
        ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
        ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
        ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
    )
    for x0, x1, x2 in faces:
        for j in range(1, degree - 1):
            for i in range(1, degree - j):
                u = i / degree
                v = j / degree
                w0 = 1.0 - u - v
                points.append(
                    (
                        w0 * x0[0] + u * x1[0] + v * x2[0],
                        w0 * x0[1] + u * x1[1] + v * x2[1],
                        w0 * x0[2] + u * x1[2] + v * x2[2],
                    )
                )

    for k in range(1, degree - 2):
        for j in range(1, degree - k - 1):
            for i in range(1, degree - j - k):
                points.append((i / degree, j / degree, k / degree))

    return np.asarray(points, dtype=np.float64)


def _curved_triangle_nodes(degree: int) -> np.ndarray:
    X = _triangle_basix_nodes(degree).copy()
    x = X[:, 0]
    y = X[:, 1]
    l0 = 1.0 - x - y
    l1 = x
    l2 = y
    X[:, 0] = x + 0.18 * l0 * l2 - 0.10 * l1 * l2
    X[:, 1] = y + 0.24 * l0 * l1 + 0.08 * l1 * l2
    return X


def _curved_tetra_nodes(degree: int) -> np.ndarray:
    X = _tetra_basix_nodes(degree).copy()
    x = X[:, 0]
    y = X[:, 1]
    z = X[:, 2]
    l0 = 1.0 - x - y - z
    l1 = x
    l2 = y
    l3 = z
    X[:, 0] = x + 0.10 * l0 * l2 - 0.05 * l1 * l3
    X[:, 1] = y + 0.10 * l0 * l1 + 0.05 * l2 * l3
    X[:, 2] = z + 0.18 * l0 * l1 + 0.10 * l0 * l2 - 0.08 * l1 * l2
    return X


def _phi(X: np.ndarray) -> np.ndarray:
    return np.sqrt((X[:, 0] - 0.35) ** 2 + (X[:, 1] - 0.30) ** 2 + (X[:, 2] - 0.20) ** 2) - 0.28


def _make_level_set(cell_family: str, degree: int) -> cutcells.LevelSetFunction:
    if cell_family == "triangle":
        dof_coordinates = _curved_triangle_nodes(degree)
        vtk_type = np.array([5], dtype=np.int32)
        tdim = 2
        expected_local_dofs = (degree + 1) * (degree + 2) // 2
    elif cell_family == "tetra":
        dof_coordinates = _curved_tetra_nodes(degree)
        vtk_type = np.array([10], dtype=np.int32)
        tdim = 3
        expected_local_dofs = (degree + 1) * (degree + 2) * (degree + 3) // 6
    else:
        raise ValueError(f"Unsupported cell_family={cell_family}")

    if dof_coordinates.shape[0] != expected_local_dofs:
        raise RuntimeError(
            f"{cell_family} p{degree}: generated {dof_coordinates.shape[0]} dofs, "
            f"expected {expected_local_dofs}"
        )

    cell_dofs = np.arange(expected_local_dofs, dtype=np.int32)
    cell_offsets = np.array([0, expected_local_dofs], dtype=np.int32)
    mesh_data = cutcells.create_level_set_mesh_data(
        dof_coordinates=dof_coordinates,
        cell_dofs=cell_dofs,
        cell_offsets=cell_offsets,
        degree=degree,
        tdim=tdim,
        cell_types=vtk_type,
    )
    dof_values = _phi(dof_coordinates)
    return cutcells.create_level_set_function(mesh_data, dof_values)


def _run_case(cell_family: str, degree: int, output_dir: Path) -> None:
    ls = _make_level_set(cell_family, degree)
    output_path = output_dir / f"curved_{cell_family}_p{degree}.vtu"
    cutcells.write_level_set_vtu(str(output_path), ls, field_name="phi")

    print(f"[{cell_family} p{degree}] cell_dofs={ls.mesh_data.cell_num_dofs(0)}")
    print(f"[{cell_family} p{degree}] num_dofs={ls.mesh_data.num_dofs()}")
    print(f"[{cell_family} p{degree}] wrote {output_path}")
    print(
        f"[{cell_family} p{degree}] ParaView: set 'Nonlinear Subdivision Level' > 0 "
        "to visualize the curved high-order cell."
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Write single curved high-order Lagrange cells for ParaView inspection."
    )
    parser.add_argument(
        "--cell-family",
        choices=["triangle", "tetra", "both"],
        default="both",
        help="Cell family to write.",
    )
    parser.add_argument(
        "--degree",
        choices=["2", "3", "4", "all"],
        default="all",
        help="Polynomial degree.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("demo_level_set_ho_curved_vtu"),
        help="Output directory for VTU files.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    degrees = _parse_degrees(args.degree)
    families = ["triangle", "tetra"] if args.cell_family == "both" else [args.cell_family]

    for family in families:
        for degree in degrees:
            _run_case(family, degree, args.output_dir)


if __name__ == "__main__":
    main()
