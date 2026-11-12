# Tutorials

These tutorials are structured as a software-demo book. Each chapter starts
from a runnable script in `python/demo`, states the geometry and the exact
values it is checked against, then follows the CutCells workflow used in the
code: a background mesh, level sets, `cut`, selections, and quadrature rules,
visualisation meshes or VTU files of the selected parts.

The figures are the PyVista views of the demos themselves, rendered off
screen; the 3D chapters add views you can turn and zoom. The demos need NumPy
only; PyVista and matplotlib draw the views and plots, and `--no-plot` skips
them:

```bash
pip install cutcells pyvista matplotlib
python python/demo/demo_disk.py
```

## One Level Set

The first three chapters show the workflow. The next four are measurements:
what each choice costs and what it buys, for one level set. The last two cover
what comes after the rule.

```{toctree}
:maxdepth: 1

disk
convergence
sphere
degree-order
cell-shapes
points-accuracy
compression
ellipsoid
mesh-aligned
```

1. Start with [Disk in a square](disk.md) for the whole workflow in 2D: a
   level set, the three parts it makes, areas, lengths and an integral of a
   function.
2. Use [Convergence](convergence.md) to compare straight pieces, curved
   pieces of an interpolated level set, and exact geometry under refinement.
3. Move to 3D with [Ball in a tetrahedral mesh](sphere.md): exact and
   interpolated level sets, curved Lagrange output.
4. [Degree and order](degree-order.md) turns the two dials separately: which
   error each one controls, and what a higher degree costs.
5. [Cell shapes](cell-shapes.md) cuts one ball on triangles, quadrilaterals,
   tetrahedra, hexahedra, prisms and pyramids, with both backends.
6. [Points and accuracy](points-accuracy.md) compares lookup tables, quadrays
   and Algoim by the error they reach for a number of points.
7. [Compression](compression.md) reduces a rule to the moments you need, and
   measures when that pays and what it gives up.
8. [Integrating functions](ellipsoid.md) uses a rule on an ellipsoid: physical
   points, moments, normals and the divergence theorem.
9. [Plane on mesh facets](mesh-aligned.md) covers surfaces that lie on the mesh,
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
