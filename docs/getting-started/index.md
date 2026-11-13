# Getting started

```{toctree}
:maxdepth: 1

installation
quickstart
```

`pip install cutcells` installs the Python package with its C++ core. To
develop CutCells, or to build downstream packages such as CutFEMx against it,
install the C++ library into the active environment and build the Python
package against that installed CMake package. The quickstart uses only NumPy
and the structured mesh helpers, so it is useful as a smoke test before the
PyVista-based [tutorials](../tutorials/index.md).
