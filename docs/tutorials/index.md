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

```{toctree}
:maxdepth: 1
:caption: Tutorial Book

disk
convergence
sphere
sphere-plane
napkin-ring
gyroid
```

## Suggested Session Flow

1. Start with [Disk in a square](disk.md) for the whole workflow in 2D: a
   level set, the three parts it makes, areas, lengths and an integral of a
   function.
2. Use [Convergence](convergence.md) to compare straight pieces, curved
   pieces of an interpolated level set, and exact geometry under refinement.
3. Move to 3D with [Ball in a tetrahedral mesh](sphere.md): exact and
   interpolated level sets, curved Lagrange output.
4. Use [Ball and plane](sphere-plane.md) to select volumes, surfaces and
   curves with conditions on two level sets.
5. Use [Napkin ring](napkin-ring.md) for a level set made of level sets, with
   sharp edges and the faces of its surface.
6. Finish with [Gyroid](gyroid.md): a triply periodic minimal surface on
   hexahedra, and a sheet selected between two level sets.
