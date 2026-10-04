# Level sets

A level set defines implicit geometry through a scalar function `phi(x)`.
CutCells uses the sign convention:

- `phi < 0`: inside region.
- `phi > 0`: outside region.
- `phi = 0`: embedded interface.

Selections refer to level sets: a condition such as `phi < 0` by the object
itself, a string such as `"phi < 0"` by its name. `cut` takes two kinds:

- **Pk level sets** (`LevelSetFunction` with degree-of-freedom values),
  interpolated on the background mesh;
- **analytic level sets** (`AnalyticLevelSet`), evaluated exactly, with
  bounds that classify cells and drive the quadrays backend.

A level set can also be **made of level sets** of these kinds, for shapes
with sharp edges (see below).

## Batched interpolation

The primary helper is `create_level_set`:

```python
import numpy as np
import cutcells

mesh = cutcells.structured_triangle_mesh_view(-1, -1, 1, 1, 41, 41)

def circle(X):
    # X has shape (gdim, num_dofs).
    return np.sqrt(X[0] ** 2 + X[1] ** 2) - 0.5

phi = cutcells.create_level_set(mesh, circle, degree=1, name="phi")
```

`create_level_set` builds level-set mesh data, evaluates the callback once on
all interpolation points, and returns a `LevelSetFunction`. The callback must
return a one-dimensional array of length `num_dofs`.

For applications that prefer one point per row, use `interpolate_level_set`:

```python
def circle_rows(X):
    # X has shape (num_dofs, gdim).
    return np.sqrt(X[:, 0] ** 2 + X[:, 1] ** 2) - 0.5

phi = cutcells.interpolate_level_set(mesh, circle_rows, degree=1, name="phi")
```

Both helpers create Pk level sets on the background mesh, of any degree
(pyramids take degree 2 at most).

## Reusing level-set mesh data

When several level sets share the same interpolation space, build the mesh data
once and then create functions from arrays:

```python
mesh_data = cutcells.create_level_set_mesh_data(mesh, degree=1, merge_tol=-1.0)
values = np.asarray(circle_rows(mesh_data.dof_coordinates), dtype=np.float64)
phi = cutcells.create_level_set_function(mesh_data, values, name="phi")
```

This path is useful when level-set values come from another finite element code
or from a cached array.

## Analytic level sets

An analytic level set gives values, gradients, and bounds of both over boxes.
`cut` classifies cells by those bounds, and the quadrays backend integrates
the level set itself, with no interpolant:

```python
ball = cutcells.analytic_sphere([0.0, 0.0, 0.0], 0.5)  # |x - c| - r
result = cutcells.cut(mesh, ball)                       # named "phi"
named = cutcells.analytic_sphere([0.0, 0.0, 0.0], 0.5, name="ball")
```

Other sources:

- `cutcells.AnalyticLevelSet(value, gradient, box_bounds, taylor_bounds=None)`
  from Python callables (slow, for experiments);
- `cutcells.AnalyticLevelSet.from_capsule(capsule)` from a capsule holding the
  C struct of function pointers (fast, for compiled geometry libraries);
- `cutcells.shapeforest.analytic_level_set(shape)` for ShapeForest shapes,
  expressions and tapes (ShapeForest is optional).

Each takes `name=`; `cut` calls an analytic level set without one `phi`.

## Multiple level sets

Each level set needs a stable name, which strings and `str(condition)` use:

```python
sphere = cutcells.create_level_set(mesh, sphere_phi, degree=1, name="sphere")
plane = cutcells.create_level_set(mesh, plane_phi, degree=1, name="plane")

result = cutcells.cut(mesh, [sphere, plane])
disk = result[(plane == 0) & (sphere < 0)]
curve = result[(plane == 0) & (sphere == 0)]
```

Analytic level sets without a name are called `phi1`, `phi2`, ... in a list.
`names` overrides all names: `cutcells.cut(mesh, [ball, plane], names=["ball",
"plane"])`. The result stores `level_set_names` and, per level set and cell,
the cell's domain (`result.domains`).

## Level sets made of level sets

Shapes built with unions, intersections and differences have sharp edges and
corners. Given as one level set, they lose accuracy along the edges. Given as
the smooth level sets they are made of, the edges are exact. A
`LevelSetFunction` takes them as `level_sets`, and `inside`, the selection of
them where it is negative:

```python
ball = cutcells.analytic_sphere(c, 0.8, signed_distance=False, name="ball")
cyl = cutcells.create_level_set(mesh, cylinder, degree=2, name="cyl")

# a napkin ring: the ball outside the cylinder
ring = cutcells.LevelSetFunction(level_sets=[ball, cyl],
                                 inside=(ball < 0) & (cyl > 0), name="ring")
ring.inside           # "ball < 0 and cyl > 0"

result = cutcells.cut(mesh, [ring, fluid])
result.level_set_names   # ['ring.ball', 'ring.cyl', 'fluid']
result.insides           # {'ring': 'ball < 0 and cyl > 0'}

solid = result[ring < 0]
surface = result[ring == 0]
hole = surface.on("cyl")    # the face on the cylinder
```

- **Its level sets** are `LevelSetFunction`s of its precision, or
  `AnalyticLevelSet`s with a name. Their names must differ, and each is a
  word without spaces, parentheses, `<`, `>` or `=`. `ring.level_sets` returns
  them as given, so conditions on them select in a cut: `result[cyl < 0]`.
- **Its inside** is a condition built from them, or a string with their names,
  of clauses `< 0` and `> 0` joined with `and` and `or` (`&` and `|` in
  conditions), with parentheses as needed.
- **Its own value** is optional. `value=` takes a callable or an
  `AnalyticLevelSet`, such as the whole shape, for code that evaluates the
  level set; `cut` uses its level sets.
- **In a cut**, its level sets take its place, named after it (`ring.ball`).
  A cut takes up to 64 of them, counting all levels.
- **Nesting:** a level set made of level sets can be one of the `level_sets` of
  another; its own are then named `outer.ring.ball`.

Selections of it expand into selections of its level sets; see
[Cut results and selections](selections.md).
