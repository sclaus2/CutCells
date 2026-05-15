# Copyright (c) 2022 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
#
# SPDX-License-Identifier:    MIT

import importlib as _importlib
from importlib import util as _importlib_util
from pathlib import Path as _Path
import ctypes as _ctypes
import sys as _sys
import sysconfig as _sysconfig


def _load_cpp_module():
    build_dir = _Path(__file__).resolve().parents[1] / "build"
    candidates = sorted(build_dir.glob("_cutcellscpp*.so"))
    if candidates:
        ext_suffix = _sysconfig.get_config_var("EXT_SUFFIX")
        if ext_suffix:
            matching = [candidate for candidate in candidates if candidate.name.endswith(ext_suffix)]
            if matching:
                candidates = matching

        # If we're running from a source tree, prefer an up-to-date in-tree build.
        # The repo keeps `cutcells/wrapper.cpp` next to this file; if the extension
        # is older than the wrapper source or the C++ library sources, it's stale.
        wrapper_cpp = _Path(__file__).with_name("wrapper.cpp")
        source_roots = [wrapper_cpp]
        cpp_src = _Path(__file__).resolve().parents[2] / "cpp" / "src"
        if cpp_src.exists():
            try:
                source_roots.extend(cpp_src.glob("*.cpp"))
                source_roots.extend(cpp_src.glob("*.h"))
            except OSError:
                pass
        source_roots = [p for p in source_roots if p.exists()]
        if source_roots:
            try:
                newest_source_ts = max(p.stat().st_mtime for p in source_roots)
                candidates = [
                    c for c in candidates
                    if c.exists() and c.stat().st_mtime >= newest_source_ts
                ]
            except OSError:
                pass

        # Choose the newest candidate (by mtime) to avoid picking an old leftover.
        if candidates:
            try:
                candidates = sorted(candidates, key=lambda p: p.stat().st_mtime, reverse=True)
            except OSError:
                pass

        lib_candidates = [
            build_dir / "cutcells_cpp" / "src" / "libcutcells.dylib",
            _Path(__file__).resolve().parents[2]
            / "cpp"
            / "build"
            / "src"
            / "libcutcells.dylib",
        ]
        for lib in lib_candidates:
            if lib.exists():
                _ctypes.CDLL(str(lib), mode=_ctypes.RTLD_GLOBAL)
                break

        if candidates:
            try:
                spec = _importlib_util.spec_from_file_location(
                    "cutcells._cutcellscpp", candidates[0]
                )
                if spec is None or spec.loader is None:
                    raise ImportError(f"Could not load extension from {candidates[0]}")

                module = _importlib_util.module_from_spec(spec)
                _sys.modules["cutcells._cutcellscpp"] = module
                spec.loader.exec_module(module)

                # Basic API sanity check for common dev/test usage. If this fails,
                # fall back to the installed extension module.
                if not (
                    hasattr(module, "HOMeshPart_float64")
                    and hasattr(module.HOMeshPart_float64, "visualization_mesh")
                ):
                    raise ImportError(f"Stale in-tree extension loaded from {candidates[0]}")

                return module
            except Exception:
                _sys.modules.pop("cutcells._cutcellscpp", None)

    # If the extension was installed into the package directory (normal wheel
    # install), prefer loading it from right next to this __init__.py.
    # The sys.path scan below intentionally skips this package dir to avoid
    # re-import loops in a source-tree checkout, so without this branch an
    # installed package cannot load its own extension.
    this_pkg_dir = _Path(__file__).resolve().parent
    ext_suffix = _sysconfig.get_config_var("EXT_SUFFIX")
    local_sos = []
    try:
        for so in this_pkg_dir.glob("_cutcellscpp*.so"):
            if ext_suffix and not so.name.endswith(ext_suffix):
                continue
            local_sos.append(so)
    except Exception:
        local_sos = []

    if local_sos:
        local_sos = sorted(local_sos, key=lambda p: p.stat().st_mtime, reverse=True)
        so_path = local_sos[0]

        lib_candidates = [
            this_pkg_dir / "libcutcells.dylib",
            this_pkg_dir.parent / "libcutcells.dylib",
            this_pkg_dir.parent / "lib" / "libcutcells.dylib",
        ]
        for lib in lib_candidates:
            try:
                if lib.exists():
                    _ctypes.CDLL(str(lib), mode=_ctypes.RTLD_GLOBAL)
                    break
            except Exception:
                pass

        return _importlib.import_module("cutcells._cutcellscpp")

    # Fall back to an installed extension (e.g. site-packages) when running from
    # a source tree without an up-to-date in-tree build.
    installed = []
    for entry in list(_sys.path):
        try:
            pkg_dir = _Path(entry) / "cutcells"
            if not pkg_dir.exists() or pkg_dir.resolve() == this_pkg_dir:
                continue
            for so in pkg_dir.glob("_cutcellscpp*.so"):
                if ext_suffix and not so.name.endswith(ext_suffix):
                    continue
                installed.append(so)
        except Exception:
            continue

    if installed:
        installed = sorted(installed, key=lambda p: p.stat().st_mtime, reverse=True)
        so_path = installed[0]

        # Try to preload the companion dylib if it was packaged separately.
        lib_candidates = [
            so_path.parent / "libcutcells.dylib",
            so_path.parent.parent / "libcutcells.dylib",
            so_path.parent.parent / "lib" / "libcutcells.dylib",
        ]
        for lib in lib_candidates:
            try:
                if lib.exists():
                    _ctypes.CDLL(str(lib), mode=_ctypes.RTLD_GLOBAL)
                    break
            except Exception:
                pass

        spec = _importlib_util.spec_from_file_location("cutcells._cutcellscpp", so_path)
        if spec is None or spec.loader is None:
            raise ImportError(f"Could not load installed extension from {so_path}")
        module = _importlib_util.module_from_spec(spec)
        _sys.modules["cutcells._cutcellscpp"] = module
        spec.loader.exec_module(module)
        return module

    raise ImportError(
        "Could not load cutcells extension: no in-tree build found and no installed "
        "cutcells/_cutcellscpp*.so found on sys.path."
    )


_cutcellscpp = _load_cpp_module()

classify_cell_domain = _cutcellscpp.classify_cell_domain
CellType = _cutcellscpp.CellType
CutCell_float32 = _cutcellscpp.CutCell_float32
CutCell_float64 = _cutcellscpp.CutCell_float64
CutCells_float32 = _cutcellscpp.CutCells_float32
CutCells_float64 = _cutcellscpp.CutCells_float64
CutMesh_float32 = _cutcellscpp.CutMesh_float32
CutMesh_float64 = _cutcellscpp.CutMesh_float64
QuadratureRules_float32 = _cutcellscpp.QuadratureRules_float32
QuadratureRules_float64 = _cutcellscpp.QuadratureRules_float64
EdgeRootTag = _cutcellscpp.EdgeRootTag
CellCertTag = _cutcellscpp.CellCertTag
FaceCertTag = _cutcellscpp.FaceCertTag
MeshView = _cutcellscpp.MeshView
MeshView_float32 = _cutcellscpp.MeshView_float32
MeshView_float64 = _cutcellscpp.MeshView_float64
AdaptCell = _cutcellscpp.AdaptCell
AdaptCell_float32 = _cutcellscpp.AdaptCell_float32
AdaptCell_float64 = _cutcellscpp.AdaptCell_float64
LevelSetCell = _cutcellscpp.LevelSetCell
LevelSetCell_float32 = _cutcellscpp.LevelSetCell_float32
LevelSetCell_float64 = _cutcellscpp.LevelSetCell_float64
LevelSetMeshData = _cutcellscpp.LevelSetMeshData
LevelSetMeshData_float32 = _cutcellscpp.LevelSetMeshData_float32
LevelSetMeshData_float64 = _cutcellscpp.LevelSetMeshData_float64
LevelSetFunction = _cutcellscpp.LevelSetFunction
LevelSetFunction_float32 = _cutcellscpp.LevelSetFunction_float32
LevelSetFunction_float64 = _cutcellscpp.LevelSetFunction_float64
create_level_set_mesh_data = _cutcellscpp.create_level_set_mesh_data
create_level_set_function = _cutcellscpp.create_level_set_function
create_level_set = _cutcellscpp.create_level_set
interpolate_level_set = _cutcellscpp.interpolate_level_set
make_adapt_cell = _cutcellscpp.make_adapt_cell
build_edges = _cutcellscpp.build_edges
make_cell_level_set = _cutcellscpp.make_cell_level_set
evaluate_bernstein = _cutcellscpp.evaluate_bernstein
extract_parent_edge_bernstein = _cutcellscpp.extract_parent_edge_bernstein
restrict_edge_bernstein_exact = _cutcellscpp.restrict_edge_bernstein_exact
classify_edge_roots = _cutcellscpp.classify_edge_roots
classify_new_edges = _cutcellscpp.classify_new_edges
fill_all_vertex_signs_from_level_set = _cutcellscpp.fill_all_vertex_signs_from_level_set
classify_leaf_cells = _cutcellscpp.classify_leaf_cells
process_ready_to_cut_cells = _cutcellscpp.process_ready_to_cut_cells
refine_ready_cell_on_largest_midpoint_value = (
    _cutcellscpp.refine_ready_cell_on_largest_midpoint_value
)
refine_green_on_multiple_root_edges = _cutcellscpp.refine_green_on_multiple_root_edges
refine_red_on_ambiguous_cells = _cutcellscpp.refine_red_on_ambiguous_cells
certify_and_refine = _cutcellscpp.certify_and_refine
certify_refine_and_process_ready_cells = (
    _cutcellscpp.certify_refine_and_process_ready_cells
)
cut = _cutcellscpp.cut
higher_order_cut = _cutcellscpp.higher_order_cut
create_cut_mesh = _cutcellscpp.create_cut_mesh
locate_cells = _cutcellscpp.locate_cells
cut_vtk_mesh = _cutcellscpp.cut_vtk_mesh
cut_mesh_view = _cutcellscpp.cut_mesh_view
csr_to_vtk_cells = _cutcellscpp.csr_to_vtk_cells
compute_physical_cut_vertices = _cutcellscpp.compute_physical_cut_vertices
complete_from_physical = _cutcellscpp.complete_from_physical
make_quadrature = _cutcellscpp.make_quadrature
runtime_quadrature = _cutcellscpp.runtime_quadrature
physical_points = _cutcellscpp.physical_points
write_vtk = _cutcellscpp.write_vtk
write_level_set_vtu = _cutcellscpp.write_level_set_vtu
ho_cut = _cutcellscpp.ho_cut
HOCutResult = _cutcellscpp.HOCutResult
HOMeshPart = _cutcellscpp.HOMeshPart

from .mesh_utils import (
    cutmesh_to_pyvista,
    mesh_from_pyvista,
    rectangle_triangle_mesh,
    rectangle_quad_mesh,
    safe_part_name,
    structured_triangle_mesh_view,
    box_tetrahedron_mesh,
    box_hex_mesh,
)

try:
    from . import triangulation_analysis as triangulation_analysis
    from .triangulation_analysis import (
        CellCase,
        analyze_all_cases,
        analyze_single_cell_case,
        classify_new_triangulation_edges,
        classify_new_triangulation_simplices,
        enumerate_cases,
        single_cell_level_set,
        single_cell_mesh,
        summarize_analysis,
    )
except (ModuleNotFoundError, ImportError) as exc:
    if isinstance(exc, ModuleNotFoundError) and exc.name != f"{__name__}.triangulation_analysis":
        raise
