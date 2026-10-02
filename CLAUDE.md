# CutCells — Coding Style Guide for AI Agents

> This document describes the conventions of the CutCells library so that
> generated code blends in with the existing codebase.  Follow these rules
> strictly when adding new files, functions, or Python bindings.

Please always use the conda env fenicsx0.11, with the tool paths, compilers and
discovery variables that the workspace's `AGENTS.md` (one level above this
repository) gives.  Do not use fenicsxv10 or other legacy environments, and do
not install CutCells into base.  Section 17 gives the repository layout and the
build and test commands.

---

## 1  File-level boilerplate

Every file starts with the ONERA copyright header.  Use the year the file
is created and keep the MIT SPDX tag.

```cpp
// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT
```

Headers always use `#pragma once` (no include guards).

---

## 2  Language standard and compiler features

- C++20.  Use `<concepts>`, `<span>`, `std::floating_point`, `std::integral`.
- Templates are constrained with `template <std::floating_point T>` (not SFINAE).
- Prefer `std::span<const T>` for read-only contiguous input.
  Mutable output goes into a reference to `std::vector<T>` or a pre-sized
  `std::span<T>`.
- Use `int` / `int32_t` for indices and sizes in the API.
  Use `std::size_t` only for STL container arithmetic in implementation bodies.
- Avoid `auto` in function signatures.  `auto` is fine inside function bodies
  where the type is obvious.

---

## 3  Namespace layout

All code lives under the `cutcells` namespace.  Sub-namespaces are used
for logical grouping:

| Namespace                       | Purpose |
|---------------------------------|---------|
| `cutcells`                      | Top-level structs (`LocalMesh`, `MeshView`, `LevelSetFunction`, `BernsteinCell`, enums) |
| `cutcells::cell`                | Cell-level operations: `CutCell`, `cut`, `volume`, `type`, `domain`, cell topology helpers |
| `cutcells::cell::edge_root`     | Root-finding on edges (`method`, `RootSolveInfo`, solver functions) |
| `cutcells::mesh`                | Mesh-level structures and operations (`CutCells`, `CutMesh`, `cut_vtk_mesh`) |
| `cutcells::quadrature`          | Quadrature rules (`QuadratureRules`, `append_quadrature`, `runtime_quadrature`) |
| `cutcells::math`                | Small vector helpers (`Vec3`, `cross`, `dot`, `distance`) |

When adding a new logical area, prefer adding a new sub-namespace over
polluting existing ones.  The namespace declaration style uses either the
fully-qualified form `namespace cutcells::cell {` or the nested form:

```cpp
namespace cutcells
{
    namespace cell
    {
        // ...
    }
}
```

Both styles appear in the codebase; the fully-qualified form is preferred
for new code.

---

## 4  Architecture: flat structs + free functions

This is the **central design principle**.  Follow it strictly.

### 4.1  Data lives in plain structs

Structs are plain-old-data containers with public fields.  They may contain
convenience accessors (`n_vertices()`, `n_cells()`, `has_value()`) but
**no mutating methods, no constructors with logic, no virtual dispatch**.

Field naming uses a leading underscore for "legacy" structs (`CutCell`,
`CutMesh`, `QuadratureRules`) and plain names for newer structs
(`LocalMesh`, `MeshView`, `LevelSetFunction`).  **For new code, use plain
names without leading underscore.**

```cpp
// ✓  New struct style
template <std::floating_point T>
struct MyData
{
    int gdim = 0;
    int tdim = 0;
    std::vector<T> vertex_coords;
    std::vector<int32_t> connectivity;
    std::vector<int32_t> offsets;

    int n_vertices() const { return gdim > 0 ? static_cast<int>(vertex_coords.size()) / gdim : 0; }
    int n_cells()    const { return offsets.empty() ? 0 : static_cast<int>(offsets.size()) - 1; }
};
```

Use CSR layout for variable-length connectivity:

- `connectivity` — flattened vertex indices
- `offsets`      — size `n_cells + 1`, `offsets[0] = 0`

### 4.2  Logic lives in free functions

All algorithmic work is done in free functions that take the struct by
reference (const or mutable).  Functions that produce new data typically
write into an output reference parameter rather than returning by value
(but returning by value is also acceptable for lightweight results).

```cpp
// ✓  Preferred: output parameter
template <std::floating_point T>
void compute_foo(const MyData<T>& input, MyData<T>& output);

// ✓  Acceptable: return by value for small / self-contained results
template <std::floating_point T>
T compute_volume(const MyData<T>& data);
```

**Never add member functions that mutate state beyond trivial accessors.**

### 4.3  Why this matters

- Flat arrays are trivially wrapped in nanobind (zero-copy to NumPy).
- Free functions compose without class hierarchies.
- Easier to test: construct a struct, call a function, inspect the struct.

---

## 5  Enums

Use `enum class` with explicit underlying type when relevant.

```cpp
enum class EdgeState : uint8_t
{
    no_root = 0,
    one_root = 1,
    multiple_roots = 2,
};
```

Provide `_to_str` / `string_to_` conversion free functions when the enum
crosses the Python boundary.

---

## 6  Header / source split

- **Header (.h):** struct definitions, enum definitions, free-function
  declarations, small `inline` helpers, and template function declarations.
- **Source (.cpp):** template explicit instantiations and non-trivial
  function bodies.

Template functions that are only instantiated for `float` and `double`
are declared in the header and explicitly instantiated at the bottom of
the `.cpp`:

```cpp
// my_feature.cpp  — bottom of file
template void my_func<float>(...);
template void my_func<double>(...);
```

Purely header-only utilities (like `span_math.h`, `cell_types.h`) are fine
for small helpers.

---

## 7  Naming conventions

| Entity            | Convention          | Example |
|-------------------|---------------------|---------|
| Namespace         | `snake_case`        | `cutcells::quadrature` |
| Struct / Class    | `PascalCase`        | `CutCell`, `LocalMesh`, `QuadratureRules` |
| Enum class        | `PascalCase`        | `EdgeState`, `LocalLevelSetBackend` |
| Enum value        | `snake_case`        | `one_root`, `nodal_signs` |
| Free function     | `snake_case`        | `compute_edge_root`, `classify_local_edges` |
| Member accessor   | `snake_case`        | `n_vertices()`, `has_value()` |
| Template param    | single uppercase    | `T`, `I`, `N` |
| Struct field (new)| `snake_case`        | `vertex_x`, `cell_offsets`, `gdim` |
| Struct field (old)| `_snake_case`       | `_vertex_coords`, `_types` |
| Local variable    | `snake_case`        | `n_cells`, `edge_id` |
| Constant          | `snake_case`        | used inline or as `constexpr` |

---

## 8  Function signature patterns

### Input spans

Read-only contiguous arrays are passed as `std::span<const T>`:

```cpp
template <std::floating_point T>
void cut(const type cell_type,
         const std::span<const T> vertex_coordinates,
         const int gdim,
         const std::span<const T> ls_values,
         const std::string& cut_type_str,
         CutCell<T>& cut_cell,
         bool triangulate = false);
```

### Scalar parameters

Pass by value (`int gdim`, `bool triangulate`) or by const ref for strings
(`const std::string& cut_type_str`).  Do **not** pass scalars by const
reference (e.g. avoid `const int& gdim` in new code — some legacy code
does this but it is not preferred).

### Output parameters

Mutable struct references as the last non-default parameter:

```cpp
template <std::floating_point T>
void init_local_mesh_from_template(
    LocalMesh<T>&             mesh,       // output
    const RefinementTemplate& tpl,        // input
    std::span<const T>        parent_cell_coords,
    cell::type                parent_cell_type,
    int                       parent_cell_id,
    int                       n_level_sets = 1);
```

### Return values

Functions returning a single scalar or small struct may return by value:

```cpp
template <std::floating_point T>
T volume(const CutCell<T>& cut_cell);

template <std::floating_point T>
CutCell<T> merge(std::vector<CutCell<T>> cut_cell_vec);
```

---

## 9  Implementation patterns

### Anonymous namespace for file-local helpers

```cpp
namespace
{
  struct MergeVertexKey { ... };
  // file-local helpers that don't leak into the header
}
```

### Error handling

Use `throw std::invalid_argument(...)` or `throw std::runtime_error(...)`
for precondition violations.  No custom exception types.  Use `assert()`
only for internal invariants that should never fail.

### Tolerances

Pass as a parameter with a default:

```cpp
T tol = static_cast<T>(1e-14)
```

Use `static_cast<T>(...)` to convert literal doubles into the template
type.

---

## 10  Python bindings (nanobind) — **this is critical**

The Python wrapper lives in `wrapper.cpp` and uses nanobind.  The user
interacts **primarily through Python**.  Efficient, zero-copy bindings are
essential.

### 10.1  Module structure

All bindings are registered in a single `NB_MODULE(_cutcellscpp, m)` block.
Type-generic bindings are factored into template helper functions that are
called for each scalar type:

```cpp
template <typename T>
void declare_float(nb::module_& m, std::string type)
{
    // register all T-dependent free functions and struct bindings here
}

NB_MODULE(_cutcellscpp, m)
{
    // enums, non-templated types
    nb::enum_<cell::type>(m, "CellType")
      .value("triangle", cell::type::triangle)
      ...;

    declare_float<float>(m, "float32");
    declare_float<double>(m, "float64");

    // default aliases for double
    if constexpr (std::is_same_v<T, double>)
    {
      m.attr("MyStruct") = m.attr("MyStruct_float64");
      m.attr("my_function") = m.attr("my_function_float64");
    }
}
```

### 10.2  Array type aliases

```cpp
template <typename T>
using ndarray1 = nb::ndarray<const T, nb::numpy, nb::shape<-1>, nb::c_contig>;

template <typename T>
using ndarray2 = nb::ndarray<const T, nb::numpy, nb::shape<-1, -1>, nb::c_contig>;
```

### 10.3  Converting C++ vectors to NumPy (zero-copy)

Use the `as_nbarray` helper that moves the vector into a capsule-owned
heap allocation:

```cpp
template <typename V>
auto as_nbarray(V&& x, std::size_t ndim, const std::size_t* shape)
{
  using _V = std::decay_t<V>;
  _V* ptr = new _V(std::move(x));
  return nb::ndarray<typename _V::value_type, nb::numpy>(
      ptr->data(), ndim, shape,
      nb::capsule(ptr, [](void* p) noexcept { delete (_V*)p; }));
}
```

Overloads exist for 1-D (`as_nbarray(std::move(vec))`) and shaped
(`as_nbarray(std::move(vec), {rows, cols})`).

### 10.4  Converting NumPy to std::span (zero-copy input)

Inside lambda bindings, convert input arrays to spans:

```cpp
m.def("my_function",
    [](const ndarray1<T>& coords, const ndarray1<int>& connectivity) {
        std::span<const T> coords_span(coords.data(), coords.size());
        std::span<const int> conn_span(connectivity.data(), connectivity.size());
        // call C++ free function
    });
```

Or more concisely:

```cpp
std::span(coords.data(), coords.size())
```

### 10.5  GIL management

**Release the GIL** around pure C++ computation.  Re-acquire only when
calling back into Python:

```cpp
m.def("expensive_cut", [](const ndarray1<T>& ls_vals, ...) {
    // ... prepare spans while GIL is held ...
    nb::gil_scoped_release release;
    return mesh::cut_vtk_mesh<T>(...);
});
```

If a C++ function calls a Python callback (e.g. `LevelSetFunction`),
the callback wrapper must re-acquire the GIL:

```cpp
value_fn = [value_callable](const T* x, int cell_id) -> T {
    nb::gil_scoped_acquire gil;
    // ... call Python ...
};
```

### 10.6  Struct binding pattern

Bind structs as nanobind classes.  Expose flat arrays as read-only
properties returning `nb::ndarray` views (zero-copy, backed by
`reference_internal` or capsule ownership):

```cpp
nb::class_<MyStruct<T>>(m, name.c_str())
  .def(nb::init<>())
  .def_prop_ro("gdim", [](const MyStruct<T>& self) { return self.gdim; })
  .def_prop_ro("vertex_coords",
    [](const MyStruct<T>& self) {
      return nb::ndarray<const T, nb::numpy, nb::shape<-1>, nb::c_contig>(
        self.vertex_coords.data(), {self.vertex_coords.size()},
        nb::cast(self, nb::rv_policy::reference));
    },
    nb::rv_policy::reference_internal)
  .def("n_cells", &MyStruct<T>::n_cells);
```

### 10.7  Free function binding pattern

Bind C++ free functions as Python module-level functions using lambdas
that convert numpy arrays to spans:

```cpp
m.def("my_function",
    [](const ndarray1<T>& input_array, int param) {
        nb::gil_scoped_release release;
        return do_work<T>(
            std::span(input_array.data(), input_array.size()),
            param);
    },
    nb::arg("input_array"),
    nb::arg("param"),
    "Docstring for my_function.");
```

### 10.8  Naming in Python

- Struct: `MyStruct_float32`, `MyStruct_float64`, plus a default alias
  `MyStruct = MyStruct_float64`.
- Free function: `my_function_float32`, `my_function_float64`, plus a
  default alias `my_function = my_function_float64`.
- All Python-facing names use `snake_case`.
- Enums: `PascalCase` class name, `snake_case` values.

### 10.9  `__init__.py` re-exports

Every public symbol must be re-exported in `__init__.py`:

```python
from ._cutcellscpp import (
    MyStruct,
    MyStruct_float32,
    MyStruct_float64,
    my_function,
    my_function_float32,
    my_function_float64,
)
```

---

## 11  CMakeLists.txt integration

New source files are added to `target_sources(cutcells PRIVATE ...)` and
new headers to the `HEADERS` list:

```cmake
set(HEADERS
  ...
  ${CMAKE_CURRENT_SOURCE_DIR}/my_feature.h
  PARENT_SCOPE)

target_sources(cutcells PRIVATE
  ...
  ${CMAKE_CURRENT_SOURCE_DIR}/my_feature.cpp
)
```

---

## 12  Data layout conventions

- **Coordinates** are flat: `[x0_0, x0_1, ..., x0_{gdim-1}, x1_0, ...]`.
  Size = `n_vertices * gdim`.
- **Connectivity** is CSR: flat vertex indices with a separate offsets
  array.  `offsets[0] = 0`, `offsets.size() = n_cells + 1`.
- **Per-entity arrays** are parallel: `edge_state[i]` corresponds to
  edge `i` defined by `edge_vertices[2*i], edge_vertices[2*i+1]`.
- **Level-set values** per vertex are stored interleaved when there are
  multiple level sets: `[phi_0(v0), phi_1(v0), ..., phi_0(v1), ...]`.
- **Bit masks** (`uint64_t`) for up to 64 level sets: bit `i` set means
  the condition holds for level set `i`.

---

## 13  Documentation

- Use `/// @brief` Doxygen-style comments for public API functions.
- Multi-line doc comments use `///` on each line.
- Parameters documented with `/// @param name  description`.
- Keep comments concise; the code should be self-documenting through
  naming.
- Section separators use:

```cpp
// ============================================================================
// Section Name
// ============================================================================
```

---

## 14  Testing

- Test executables are added conditionally with `BUILD_TESTING`.
- Tests are standalone `.cpp` files that call the C++ API directly.
- Python-level testing uses pytest and is the primary test path since
  the user interacts through Python.

---

## 15  Complete example: adding a new feature

Suppose you are adding a `compute_normals` feature.

### `compute_normals.h`

```cpp
// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <span>
#include <vector>

#include "cut_cell.h"
#include "cell_types.h"

namespace cutcells::cell
{

/// Compute outward unit normals for each sub-cell facet of a cut cell.
///
/// @param cut_cell  cut cell with _vertex_coords_phys filled
/// @param normals   output flat array, size = n_facets * gdim
template <std::floating_point T>
void compute_normals(const CutCell<T>& cut_cell,
                     std::vector<T>& normals);

/// Convenience overload returning the result.
template <std::floating_point T>
std::vector<T> compute_normals(const CutCell<T>& cut_cell);

} // namespace cutcells::cell
```

### `compute_normals.cpp`

```cpp
// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "compute_normals.h"
#include "span_math.h"

#include <cassert>
#include <cmath>

namespace cutcells::cell
{

template <std::floating_point T>
void compute_normals(const CutCell<T>& cut_cell,
                     std::vector<T>& normals)
{
    const int gdim = cut_cell._gdim;
    // ... implementation ...
}

template <std::floating_point T>
std::vector<T> compute_normals(const CutCell<T>& cut_cell)
{
    std::vector<T> normals;
    compute_normals(cut_cell, normals);
    return normals;
}

// Explicit instantiations
template void compute_normals<float>(const CutCell<float>&, std::vector<float>&);
template void compute_normals<double>(const CutCell<double>&, std::vector<double>&);
template std::vector<float> compute_normals<float>(const CutCell<float>&);
template std::vector<double> compute_normals<double>(const CutCell<double>&);

} // namespace cutcells::cell
```

### In `wrapper.cpp` (inside `declare_float<T>`)

```cpp
m.def("compute_normals",
    [](const cell::CutCell<T>& cut_cell) {
        std::vector<T> normals;
        {
            nb::gil_scoped_release release;
            normals = cell::compute_normals(cut_cell);
        }
        const std::size_t n = normals.size() / static_cast<std::size_t>(cut_cell._gdim);
        return as_nbarray(std::move(normals), {n, static_cast<std::size_t>(cut_cell._gdim)});
    },
    nb::arg("cut_cell"),
    "Compute outward unit normals for each sub-cell facet.");
```

### In `__init__.py`

```python
from ._cutcellscpp import (
    ...
    compute_normals,
    compute_normals_float32,
    compute_normals_float64,
)
```

### In `CMakeLists.txt`

```cmake
set(HEADERS
  ...
  ${CMAKE_CURRENT_SOURCE_DIR}/compute_normals.h
  PARENT_SCOPE)

target_sources(cutcells PRIVATE
  ...
  ${CMAKE_CURRENT_SOURCE_DIR}/compute_normals.cpp
)
```

---

## 16  Common mistakes to avoid

1. **Do not add member functions with logic to structs.**  Use free functions.
2. **Do not use `std::unique_ptr` or inheritance hierarchies.**  Flat structs + free functions.
3. **Do not forget explicit template instantiations** in `.cpp` files.
4. **Do not copy NumPy data unnecessarily.**  Use `std::span` for input, `as_nbarray(std::move(...))` for output.
5. **Do not hold the GIL during expensive C++ computation.**  Release it.
6. **Do not use `size_t` in public API signatures.**  Use `int` or `int32_t`.
7. **Do not create Python classes with methods when a free function will do.**
   The Python API mirrors the C++ one: free functions operating on struct objects.
8. **Do not forget the double alias** (`_float32` / `_float64` + default).
9. **Do not add `www.` or external dependencies** without checking CMake integration.
10. **Do not use `const int&` for scalar parameters in new code.**  Pass by value.

---

## 17  Repository layout, build and tests

- **Library code only in `cpp/src`.**  New code goes into sub-folders with
  their own `CMakeLists.txt`, added from `cpp/src/CMakeLists.txt` (for example
  `cpp/src/quadrays/`, namespace `cutcells::quadrays`).  A header keeps its
  path below `cpp/src` when installed (`cpp/src/quadrays/rules.h` becomes
  `<cutcells/quadrays/rules.h>`), so headers in a sub-folder include library
  headers with `../`, e.g. `#include "../cell_types.h"`.
- **C++ tests in `cpp/tests/<folder>/`**: plain executables that return
  non-zero on failure, registered with `cutcells_add_test` and built with
  `-DBUILD_TESTING=ON`; no test framework.  Exact references and test meshes
  shared with the benchmarks live in `cpp/tests/<folder>/support/`.
- **Python tests** in `python/tests` (pytest), **user demos** in
  `python/demo`, **comparison drivers** in `benchmarks/` (its own CMake project
  that finds the installed CutCells), **notes and results** in `docs/`.
- **Algoim is optional, for tests and benchmarks only.**  New library code
  never includes it; drivers that compare against it build with
  `-DCUTCELLS_WITH_ALGOIM=ON` in `benchmarks/`.  Geometry libraries such as
  ShapeForest are optional adapters behind a generic interface, never build
  dependencies.

Build and run the C++ tests with the `fenicsx0.11` environment:

```bash
E=$HOME/miniforge3/envs/fenicsx0.11
env CONDA_PREFIX=$E CMAKE_PREFIX_PATH=$E CC=$E/bin/clang CXX=$E/bin/clang++ $E/bin/cmake -S cpp -B build-quadrays -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DCMAKE_OSX_DEPLOYMENT_TARGET=13.4
$E/bin/cmake --build build-quadrays -j 4 && $E/bin/ctest --test-dir build-quadrays
```

The benchmarks configure against an installed CutCells, found through
`CMAKE_PREFIX_PATH` (other CutCells installs may be on the search path, so put
the wanted prefix first):

```bash
$E/bin/cmake --install build-quadrays --prefix build-quadrays/install
env CONDA_PREFIX=$E CC=$E/bin/clang CXX=$E/bin/clang++ $E/bin/cmake -S benchmarks -B build-benchmarks -DCMAKE_PREFIX_PATH="$PWD/build-quadrays/install;$E" -DCMAKE_BUILD_TYPE=Release -DCUTCELLS_WITH_ALGOIM=ON -DCMAKE_OSX_DEPLOYMENT_TARGET=13.4
```

Without `CONDA_PREFIX`, CMake picks Apple's Accelerate, which lacks the
LAPACKE symbols algoim needs.
