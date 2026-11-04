# Tutorials

These tutorials are structured as a software-demo book. Each chapter starts
from a runnable script in `python/demo`, states the geometry and the exact
values it is checked against, then follows the CutCells workflow used in the
code: a background mesh, level sets, `cut`, selections, and quadrature rules,
visualisation meshes or VTU files of the selected parts.

The figures are the PyVista views of the demos themselves, rendered off
screen; the 3D chapters add views you can turn and zoom. The demos need NumPy
only; PyVista draws the views (matplotlib the convergence plot) and is skipped
with `--no-plot`:

```bash
pip install cutcells pyvista matplotlib
python python/demo/demo_disk.py
```

## One Level Set

```{toctree}
:maxdepth: 1

disk
convergence
sphere
flower
ellipsoid
torus
mesh-aligned
```

1. Start with [Disk in a square](disk.md) for the whole workflow in 2D: a
   level set, the three parts it makes, areas, lengths and an integral of a
   function.
2. Use [Convergence](convergence.md) to compare straight pieces, curved
   pieces of an interpolated level set, and exact geometry under refinement.
3. Move to 3D with [Ball in a tetrahedral mesh](sphere.md): exact and
   interpolated level sets, curved Lagrange output.
4. [Flower](flower.md) shows a non-convex curve on quadrilaterals, and what each
   degree of the level set gains.
5. [Ellipsoid](ellipsoid.md) integrates functions: physical points, moments,
   normals and the divergence theorem.
6. [Torus](torus.md) uses hexahedra and an analytic level set written in Python.
7. [Plane on mesh facets](mesh-aligned.md) covers surfaces that lie on the mesh,
   integrated once, and the cells that take them.

## Several Level Sets

Several level sets per cell are under active development: expect their rules
and costs to change more than those of a single level set.

```{toctree}
:maxdepth: 1

sphere-plane
napkin-ring
gyroid
```

1. [Ball and plane](sphere-plane.md) selects volumes, surfaces and curves with
   conditions on two level sets.
2. [Napkin ring](napkin-ring.md) gives a shape with sharp edges as a level set
   made of level sets.
3. [Gyroid](gyroid.md) cuts a triply periodic minimal surface on hexahedra, and
   selects a sheet between two level sets.
