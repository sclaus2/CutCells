# CutCells

CutCells cuts the cells of a mesh by one or more level sets `phi(x)` and
returns the inside, outside and interface pieces with explicit connectivity and
quadrature rules. It is a geometry kernel for unfitted finite element methods
(CutFEM): the C++20 core has no dependencies beyond a compiler, and the Python
interface passes arrays to and from NumPy without copying where the layout
allows.

```console
pip install cutcells
```

Wheels are published for Linux and macOS. Other platforms build from the source
distribution, which needs a C++20 compiler and CMake.

## Example

Cut a triangle mesh by a circle and integrate the inside and the interface:

```python
import numpy as np
import cutcells

mesh = cutcells.structured_triangle_mesh_view(-1.0, -1.0, 1.0, 1.0, 41, 41)

def phi(X):
    # X has shape (gdim, num_points)
    return np.sqrt(X[0] ** 2 + X[1] ** 2) - 0.55

level_set = cutcells.create_level_set(mesh, phi, degree=2, name="phi")
result = cutcells.cut(mesh, level_set)

inside = result[level_set < 0]
interface = result[level_set == 0]

print("area     ", inside.quadrature(order=4, mode="full").weights.sum())
print("perimeter", interface.quadrature(order=4, mode="cut_only").weights.sum())
```

Parent cells can be intervals, triangles, quadrilaterals, tetrahedra,
hexahedra, prisms and pyramids, with any number of level sets per cell. Level
sets are interpolated (any degree) or analytic.

## Links

- Source, issues and the C++ library: <https://github.com/sclaus2/CutCells>
- Licence: MIT. Third-party notices for generated VTK-derived tables and
  Basix-generated quadrature tables ship with the package.
