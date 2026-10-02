# Third-Party Notices

## VTK (Visualization Toolkit)

CutCells includes *generated* clipping/cutting case tables that are derived from VTK's `vtkTableBasedClipCases.h` (TableBasedClip).

VTK is licensed under the BSD 3-Clause License.

- License text: `third_party/VTK-Copyright.txt`
- Upstream: https://github.com/Kitware/VTK
- Current tablegen default ref for future regeneration: `v9.4.2`
- Existing generated headers that record `master` were produced before an exact
  VTK commit was recorded; regenerate with `tablegen/scripts/gen_tables.py` to
  record a pinned ref in the header comments.

## Basix

CutCells includes *generated* quadrature lookup tables produced offline with
Basix `make_quadrature`.

Basix is licensed under the MIT License.

- License text: `third_party/Basix-LICENSE.txt`
- Upstream: https://github.com/FEniCS/basix

## Algoim

CutCells vendors a pinned snapshot of Algoim under `third_party/algoim`
for the comparison drivers in `benchmarks/`; the library and the Python
package do not contain it.

Algoim is licensed under a BSD-style license.

- License text: `third_party/algoim/LICENSE`
- Upstream: https://github.com/algoim/algoim
- Vendored commit: `third_party/algoim/UPSTREAM_COMMIT`
- Local patches: `third_party/algoim/PATCHES.md`
