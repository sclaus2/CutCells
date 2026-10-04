# Installation

CutCells can be used as a standalone C++ library, a Python package, or a
geometry backend for CutFEMx. Downstream packages expect the C++ CMake package
and the Python extension to come from the same environment.

## Create an environment

For a standalone CutCells development environment:

```bash
conda create -n cutcells-dev -c conda-forge \
  python=3.12 \
  cmake ninja pkg-config compilers \
  scikit-build-core nanobind \
  numpy pytest pyvista
conda activate cutcells-dev
```

Tie CMake, Python, and the compiler to the active environment before building:

```bash
export CONDA_PREFIX="${CONDA_PREFIX:?activate the conda environment first}"
export CMAKE_PREFIX_PATH="$CONDA_PREFIX"
export PYTHONNOUSERSITE=1
export CC="$CONDA_PREFIX/bin/clang"
export CXX="$CONDA_PREFIX/bin/clang++"
```

On macOS with the conda-forge clang/libc++ toolchain, also use:

```bash
export CMAKE_ARGS="-DCMAKE_OSX_DEPLOYMENT_TARGET=13.4"
export CXXFLAGS="-D_LIBCPP_ENABLE_EXPERIMENTAL -D_LIBCPP_HAS_NO_EXPERIMENTAL_TZDB -D_LIBCPP_HAS_NO_EXPERIMENTAL_SYNCSTREAM -D_LIBCPP_HAS_NO_INCOMPLETE_PSTL"
```

The same build shape also works inside a FEniCSx 0.11 environment used by
CutFEMx. Keep `CONDA_PREFIX`, `CMAKE_PREFIX_PATH`, `CC`, and `CXX` pointed at
that FEniCSx environment.

## Build the C++ core

From the CutCells repository root:

```bash
cmake -G Ninja -S cpp -B build/cutcells \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=13.4
cmake --build build/cutcells --target install
```

On Linux, omit `-DCMAKE_OSX_DEPLOYMENT_TARGET=13.4`.

Add `-DBUILD_TESTING=ON` to build the C++ tests, then run them with
`ctest --test-dir build/cutcells`.

The library itself needs nothing beyond a C++20 compiler. Algoim is used only
by the comparison drivers in `benchmarks/`, a separate CMake project that finds
the installed CutCells and takes `-DCUTCELLS_WITH_ALGOIM=ON`.

## Build the Python package

Build the Python extension against the installed C++ package:

```bash
python -m pip install --no-build-isolation --no-deps --force-reinstall \
  ./python
```

`--no-build-isolation` is important because the extension must use the CMake
package installed in the active environment. `--no-deps` prevents `pip` from
changing the tested dependency stack.

## Verify the install

```bash
python - <<'PY'
import importlib.metadata as metadata
import cutcells

print("cutcells", metadata.version("cutcells"))
print("module", cutcells.__file__)
PY
```

For a stronger smoke test, run the Python tests:

```bash
python -m pytest -q python/tests
```

## Build these docs

Install documentation dependencies once:

```bash
cd docs
python -m pip install -r requirements.txt
```

Build HTML:

```bash
bash build_docs.sh clean
bash build_docs.sh build
```

The generated site is written to `docs/_build/html/index.html`. To serve it:

```bash
bash build_docs.sh serve --port 8000
```
