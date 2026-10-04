# User guide

```{toctree}
:maxdepth: 1

mesh-views
level-sets
selections
quadrature-output
algorithms
```

The main workflow is:

1. Represent the background mesh as a `MeshView`.
2. Build one or more named level sets: Pk level sets on that mesh, analytic
   ones, or level sets made of them.
3. Call `cut(mesh, level_set)` or `cut(mesh, [level_sets...])`; parts integrate
   with quadrays unless you pass `backend="lut"`.
4. Select mesh parts with conditions built from the level sets, such as
   `phi < 0` or `(sphere == 0) & (plane < 0)`, or with their strings,
   `"phi < 0"` or `"sphere = 0 and plane < 0"`.
5. Request quadrature rules, visualisation meshes, or VTU output from the
   selected part.
