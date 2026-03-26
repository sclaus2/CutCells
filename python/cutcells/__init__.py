# Copyright (c) 2022 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier:    MIT
import numpy as np
import importlib.machinery
import importlib.util
import pathlib
import sys


def _ensure_binary_extension_loaded():
    try:
        import cutcells._cutcellscpp  # noqa: F401

        return
    except ModuleNotFoundError:
        pass

    pkg_dir = pathlib.Path(__file__).resolve().parent
    seen = set()
    for raw_entry in sys.path:
        entry = pathlib.Path(raw_entry or ".").resolve()
        if entry in seen:
            continue
        seen.add(entry)

        candidate_dir = entry / "cutcells"
        if candidate_dir == pkg_dir:
            continue
        if not candidate_dir.is_dir():
            continue

        for suffix in importlib.machinery.EXTENSION_SUFFIXES:
            candidate = candidate_dir / f"_cutcellscpp{suffix}"
            if not candidate.is_file():
                continue
            spec = importlib.util.spec_from_file_location(
                "cutcells._cutcellscpp", candidate
            )
            if spec is None or spec.loader is None:
                continue
            module = importlib.util.module_from_spec(spec)
            sys.modules["cutcells._cutcellscpp"] = module
            spec.loader.exec_module(module)
            return

    raise ModuleNotFoundError(
        "cutcells extension module not found. Install the C++ CutCells library in "
        "the active env, run `pip install .`, then rerun the demo."
    )


_ensure_binary_extension_loaded()


try:
    from ._cutcellscpp import (
        classify_cell_domain,
        CellType,
        EdgeRootMethod,
        LocalLevelSetBackend,
        CutCell_float32,
        CutCell_float64,
        CutCells_float32,
        CutCells_float64,
        CutMesh_float32,
        CutMesh_float64,
        QuadratureRules_float32,
        QuadratureRules_float64,
        MeshView,
        MeshView_float32,
        MeshView_float64,
        LevelSetFunction,
        LevelSetFunction_float32,
        LevelSetFunction_float64,
        cut,
        cut_from_cached_roots,
        higher_order_cut,
        create_cut_mesh,
        locate_cells,
        cut_vtk_mesh,
        cut_mesh_view,
        csr_to_vtk_cells,
        compute_physical_cut_vertices,
        complete_from_physical,
        make_quadrature,
        runtime_quadrature,
        physical_points,
        RefinementTemplate,
        iso_p1_template,
        iso_p1_ref_coords,
        LocalMesh,
        LocalMesh_float32,
        LocalMesh_float64,
        BernsteinCell,
        BernsteinCell_float32,
        BernsteinCell_float64,
        lagrange_to_bernstein_cell,
        lagrange_to_bernstein_cell_float32,
        lagrange_to_bernstein_cell_float64,
        init_local_mesh_from_template,
        init_local_mesh_from_template_float32,
        init_local_mesh_from_template_float64,
        init_local_mesh_from_cell,
        init_local_mesh_from_cell_float32,
        init_local_mesh_from_cell_float64,
        classify_edges_on_local_mesh,
        classify_edges_on_local_mesh_float32,
        classify_edges_on_local_mesh_float64,
        compute_roots_on_local_mesh,
        compute_roots_on_local_mesh_float32,
        compute_roots_on_local_mesh_float64,
        decompose_local_mesh,
        decompose_local_mesh_float32,
        decompose_local_mesh_float64,
        CurvedGlobalMesh,
        CurvedGlobalMesh_float32,
        CurvedGlobalMesh_float64,
        CurvedCutMeshResult,
        CurvedCutMeshResult_float32,
        CurvedCutMeshResult_float64,
        cut_vtk_mesh_curved,
        cut_vtk_mesh_curved_float32,
        cut_vtk_mesh_curved_float64,
        assemble_curved_interface_mesh,
        assemble_curved_interface_mesh_float32,
        assemble_curved_interface_mesh_float64,
        assemble_curved_interface_mesh_from_cache,
        assemble_curved_interface_mesh_from_cache_float32,
        assemble_curved_interface_mesh_from_cache_float64,
        assemble_curved_volume_mesh,
        assemble_curved_volume_mesh_float32,
        assemble_curved_volume_mesh_float64,
        find_root_on_segment,
        find_root_on_segment_float32,
        find_root_on_segment_float64,
        build_zero_entities,
        build_zero_entities_float32,
        build_zero_entities_float64,
        assign_zero_entity_ownership,
        assign_zero_entity_ownership_float32,
        assign_zero_entity_ownership_float64,
        build_interface_entities,
        build_interface_entities_float32,
        build_interface_entities_float64,
        build_zero_chains,
        build_zero_chains_float32,
        build_zero_chains_float64,
        build_zero_patches,
        build_zero_patches_float32,
        build_zero_patches_float64,
        curve_zero_entities,
        curve_zero_entities_float32,
        curve_zero_entities_float64,
        curve_interface_entities,
        curve_interface_entities_float32,
        curve_interface_entities_float64,
        append_interface_quadrature_curved,
        append_interface_quadrature_curved_float32,
        append_interface_quadrature_curved_float64,
        append_volume_quadrature_curved,
        append_volume_quadrature_curved_float32,
        append_volume_quadrature_curved_float64,
        make_quadrature_curved,
        make_quadrature_curved_physical_points,
        make_quadrature_curved_physical_points_float32,
        make_quadrature_curved_physical_points_float64,
        make_quadrature_curved_float32,
        make_quadrature_curved_float64,
    )
except ModuleNotFoundError as exc:
    raise ModuleNotFoundError(
        "cutcells extension module not found. Install the C++ CutCells library in "
        "the active env, run `pip install .`, then rerun the demo."
    ) from exc


_cut_vtk_mesh_curved_cpp = cut_vtk_mesh_curved


def _vtk_arrays_from_grid(grid):
    points = np.asarray(grid.points, dtype=np.float64).reshape(-1)
    connectivity = np.asarray(grid.cell_connectivity, dtype=np.int32).reshape(-1)
    offset = np.asarray(grid.offset, dtype=np.int32).reshape(-1)
    vtk_type = np.asarray(grid.celltypes, dtype=np.int32).reshape(-1)
    return points, connectivity, offset, vtk_type


def cut_vtk_mesh_curved(
    grid_or_points,
    f,
    grad=None,
    backend="bernstein",
    geom_order=4,
    vis_subdivision=3,
    tol=1e-12,
    connectivity=None,
    offset=None,
    vtk_type=None,
):
    if connectivity is None or offset is None or vtk_type is None:
        points, connectivity, offset, vtk_type = _vtk_arrays_from_grid(grid_or_points)
    else:
        points = np.asarray(grid_or_points, dtype=np.float64).reshape(-1)
        connectivity = np.asarray(connectivity, dtype=np.int32).reshape(-1)
        offset = np.asarray(offset, dtype=np.int32).reshape(-1)
        vtk_type = np.asarray(vtk_type, dtype=np.int32).reshape(-1)

    return _cut_vtk_mesh_curved_cpp(
        points,
        connectivity,
        offset,
        vtk_type,
        f,
        grad,
        backend=backend,
        geom_order=geom_order,
        vis_subdivision=vis_subdivision,
        tol=tol,
    )
