# VTK-Based Iso-Refine Table Generator

This directory contains tooling for generating `iso_refine` refinement tables by:

1. tessellating higher-order VTK reference cells in a C++ helper,
2. mapping VTK local node ids to Basix interpolation-node ordering,
3. emitting C++ `RefinementTemplate` array blocks.

## Files

- `tablegen/tools/vtk_iso_refine_dump.cpp`
  C++ helper that dumps tessellated simplices from a VTK higher-order cell as JSON.
- `tablegen/scripts/gen_iso_refine_from_vtk.py`
  Python script that remaps VTK ids to Basix ordering and emits C++ blocks.

## Build helper

```bash
cmake -S tablegen/tools -B tablegen/tools/build
cmake --build tablegen/tools/build --target vtk_iso_refine_dump
```

If your VTK CMake package is not found or is broken, set `CMAKE_PREFIX_PATH` to a valid VTK installation.

## Generate blocks

Run in an environment that has `basix`:

```bash
python \
  tablegen/scripts/gen_iso_refine_from_vtk.py \
  --helper tablegen/tools/build/vtk_iso_refine_dump \
  --cell triangle quadrilateral tetrahedron hexahedron prism \
  --orders 2 3 4 
```

Default output path:
- `cpp/src/generated/iso_refine_vtk_blocks.h`

Use a custom output path if needed:
```bash
python tablegen/scripts/gen_iso_refine_from_vtk.py ... --out /tmp/iso_refine_vtk_blocks.h
```

For pyramid:

- VTK helper currently supports pyramid only at order `2` in this build.
- For `p=3,4`, use existing in-repo tables or a different VTK build exposing high-order pyramid tessellation.

## Notes

- The emitted blocks are intended for manual integration into `cpp/src/iso_refine.cpp`.
- Child connectivity from VTK tessellation is simplex-based (`triangle`/`tetrahedron`) by default.
- Always validate generated tables with project tests after integration.
