# Python API overview

The top-level `cutcells` package re-exports the primary C++ bindings and Python
mesh utilities. `CutResult` and `MeshPart` live in `cutcells.part`, `Condition`
in `cutcells.conditions` (re-exported as `cutcells.Condition`).

## Mesh construction

| Name | Purpose |
| --- | --- |
| `MeshView` | Background mesh view from coordinates, connectivity, offsets, cell types, and `tdim`. |
| `mesh_from_pyvista(grid, tdim=None)` | Convert a PyVista `UnstructuredGrid` to `MeshView`. |
| `structured_triangle_mesh_view(...)` | Pure NumPy structured triangular mesh helper. |
| `rectangle_triangle_mesh(...)` | PyVista triangular rectangle mesh helper. |
| `rectangle_quad_mesh(...)` | PyVista quadrilateral rectangle mesh helper. |
| `box_tetrahedron_mesh(...)` | PyVista tetrahedral box mesh helper. |
| `box_hex_mesh(...)` | PyVista hexahedral box mesh helper. |

## Level sets

| Name | Purpose |
| --- | --- |
| `LevelSetFunction` | Scalar level-set function object; `cut` needs it with degree-of-freedom values, an analytic value, or level sets it is made of. `name` (read-only) gives its name. |
| `LevelSetFunction(value=None, grad=None, nodal_values=None, gdim=0, name=None, level_sets=None, inside=None)` | A level set from a callable or an `AnalyticLevelSet` (`value`), nodal values, or level sets: with `level_sets` (`LevelSetFunction`s of its precision or named `AnalyticLevelSet`s) and `inside` (a condition on them, or a string with their names), a level set made of level sets; its own value is then optional. |
| `level_set.level_sets`, `level_set.inside` | The level sets it is made of (copies), or `[]`; its inside, or `None`. |
| `LevelSetMeshData` | Interpolation mesh data for level-set degrees of freedom. |
| `create_level_set_mesh_data(mesh, degree, merge_tol=-1.0)` | Build level-set interpolation data. |
| `create_level_set_function(mesh_data, dof_values, name="phi")` | Build a named level set from mesh data and values. |
| `create_level_set(mesh, phi, degree, name="phi")` | Interpolate a batched callback where `X.shape == (gdim, num_dofs)`. |
| `interpolate_level_set(mesh, phi, degree, name="phi")` | Interpolate a batched callback where `X.shape == (num_dofs, gdim)`. |
| `AnalyticLevelSet(value, gradient, box_bounds, taylor_bounds=None, name=None)` | Analytic level set from Python callables (slow, for experiments). `name` (read-only) gives its name, or `None`. |
| `AnalyticLevelSet.from_capsule(capsule, name=None)` | Analytic level set from a capsule holding its C struct (fast, for compiled libraries); `level_set.capsule` gives one. |
| `analytic_sphere(centre, radius, signed_distance=True, name=None)` | Compiled sphere, `\|x - centre\| - radius` (a disk in 2D). |
| `analytic_level_set_from_tape(..., name=None)` | Analytic level set from the arrays of a ShapeForest tape. |
| `cutcells.shapeforest.analytic_level_set(shape, parameters=None, name=None)` | Analytic level set from a ShapeForest shape, expression or tape. |

## Conditions

| Name | Purpose |
| --- | --- |
| `level_set < 0`, `level_set > 0`, `level_set == 0` | A `Condition`, for `LevelSetFunction` and `AnalyticLevelSet`; a number other than 0 raises. For a level set made of level sets: its inside, the rest, and its faces. |
| `a & b`, `a \| b` | Conditions combined: and, or. Each comparison needs parentheses. |
| `Condition` | A selection built from level sets. `str(condition)` writes it as a selection string with the level sets' own names, `condition.text(names)` with `names[level_set]`. |
| `cutcells.conditions.selection_text(selection, level_sets, names, owner=None)` | The string of a selection for a cut of `level_sets` named `names`: a string as it is, a condition with its level sets matched by identity. For libraries that pass selections on as strings. With `owner`, the inside of the level set of that name, made of `level_sets`. |

## Cutting

| Name | Purpose |
| --- | --- |
| `cut(mesh, level_set, ...)` | Classify every cell of a `MeshView` by one level set; returns a `CutResult`. |
| `cut(mesh, [level_sets], ...)` | The same with several level sets. |
| `ho_cut(...)` | Deprecated alias of `cut`. |
| `cut(cell_type, vertex_coordinates, gdim, ls_values, cut_type_str, ...)` | Low-level single-cell cutting API (straight cuts). |
| `higher_order_cut(...)` | Low-level single-cell cutting of a second-order cell. |
| `cut_vtk_mesh(...)` | Low-level whole-mesh cutting from flat VTK arrays. |
| `locate_cells(...)` | Locate parent cells matching a sign or interface selection. |

`cut(mesh, level_sets, ...)` keyword arguments:

| Argument | Default | Meaning |
| --- | --- | --- |
| `names` | `None` | Names of the level sets, in order; by default their own, and `phi`, or `phi1`, `phi2`, ... in a list, for analytic ones without a name. |
| `name` | `None` | Name of a single level set (not with `names`). |
| `max_depth` | `12` | Bisections of a cell before an unproven sign counts as cut. |
| `backend` | `"quadrays"` | Default backend of the selected parts; `"lut"` (also `"straight"`) gives straight pieces. |
| `options` | `None` | `QuadraysOptions` or `LutOptions` for that backend. |
| `triangulate` | `False` | Lookup tables only: split pieces into simplices. |
| `triangulation` | `"classical"` | Lookup tables only: `"classical"` or `"midpoint"`. |
| `cut_approximation` | `"auto"` | Lookup tables only: template order from the level sets' degree (`"auto"`), `"linear"`, or `"iso_p1"`. |
| `cut_approximation_order` | `1` | Lookup tables only: template order for `"iso_p1"`. |
| `degree` | `None` | Lookup tables only: template order for analytic level sets with `"auto"`. |

The lookup-table keywords raise when given with another backend.

## Cut result objects

`CutResult`:

| Member | Meaning |
| --- | --- |
| `num_cells` | Number of background cells. |
| `num_cut_cells` | Number of cells that some level set cuts. |
| `num_level_sets` | Number of level sets, a level set made of level sets counted as its own. |
| `level_set_names` | Names used in selection strings; the level sets a level set made of level sets is made of are named after it (`ring.ball`). |
| `insides` | The level sets made of level sets, by name, with their insides (`{"ring": "ball < 0 and cyl > 0"}`). |
| `cut_cells` | Cells that some level set cuts, ascending. |
| `domains` | Per level set and cell: 0 inside, 1 cut, 2 outside; shape `(num_level_sets, num_cells)`. |
| `zero_faces` | Faces lying in a zero set, each with the cell that owns it. |
| `backend`, `options` | Backend and options of the parts selected from now on (settable). |
| `result[selection]` | Select a `MeshPart` with a condition (`phi < 0`, `(a > 0) & (b == 0)`) or a string (`"phi < 0"`, `"a > 0 and b = 0"`, `"(a < 0 or b < 0) and c > 0"`). |

`MeshPart`:

| Member | Meaning |
| --- | --- |
| `dim` | Selected entity dimension. |
| `num_cut_cells` | Number of selected cut cells. |
| `num_uncut_cells` | Number of cells wholly in a volume part. |
| `cut_cells` | Selected cut cells. |
| `uncut_cells` | Cells wholly in a volume part. |
| `zero_faces` | Faces lying in a zero set that belong to the part. |
| `on(*names)` | The faces of a surface on the level sets named, such as `result[ring == 0].on("cyl")`: full or short names, or a level set made of level sets for all its faces. |
| `backend`, `options` | The part's backend and options. |
| `quadrature(order=3, mode="full", backend=None, options=None)` | Quadrature rules with positive weights; `order` is the polynomial degree integrated exactly on flat pieces (1 to 10). |
| `visualization_mesh(mode="full", backend=None, degree=3, options=None)` | `QuadraysLeafMesh` of Lagrange cells (quadrays) or `CutMesh` of straight cells (lookup tables). |
| `write_vtu(filename, mode="full", backend=None, degree=3, options=None)` | Write the visualisation mesh to a `.vtu` file. |

`mode` is `"full"` (cut cells plus whole cells) or `"cut_only"`. `backend=None`
uses the part's backend.

## Backend options

`QuadraysOptions()` (set attributes after construction):

| Attribute | Meaning |
| --- | --- |
| `margin` | Required `\|d_k psi\| / \|grad psi\|` on a box for a height direction; also the accuracy control. |
| `max_depth`, `max_bisections` | Bisection limits per level and per cell. |
| `points_per_segment` | Gauss-Legendre points per segment; 0 chooses them from the order. |
| `diagnose` | Record why each bisection happened (`QuadraysStats.causes`). |
| `prune_bounds`, `split_bounds`, `diagonal_frames`, `rotation_depth`, `mask_subdivisions`, `taylor_subdivisions`, `two_roots_depth` | Engine controls; see the docstrings. |

`LutOptions(template_order=0, triangulate=False, triangulation="classical")`:

| Attribute | Meaning |
| --- | --- |
| `template_order` | Order k of the Pk-iso-P1 template that subdivides a cut cell, 1 to 4; 0: the highest degree of the level sets that cut the cell, 2 for analytic ones. |
| `triangulate` | Split the cut pieces into simplices. |

## Quadrature

| Name | Purpose |
| --- | --- |
| `QuadratureRules` | Flat quadrature storage for points, weights, offsets, and parent cells. |
| `compress_rules(rules, degree, space="tensor")` | Rules with at most as many points as the polynomial space has; returns `(QuadratureRules, CompressionStats)`. |
| `make_quadrature(cut_cells, order)` | Generate rules for low-level cut cells. |
| `runtime_quadrature(...)` | Generate rules directly from VTK-style mesh arrays and level-set values (straight cuts). |
| `physical_points(rules, points, connectivity, offsets, vtk_types)` | Map reference quadrature points to physical space. |
| `quadrays_cell_rules(...)`, `quadrays_cell_leaves(...)` | quadrays on a single cell. |

## Output and conversion

| Name | Purpose |
| --- | --- |
| `cutmesh_to_pyvista(cut_mesh)` | Convert a CutCells `CutMesh` to PyVista. |
| `QuadraysLeafMesh` | Lagrange cells of the quadrays decomposition: `points`, `connectivity`, `offsets`, `vtk_types`, `parent`, `degree`. |
| `write_vtk(filename, vertex_coordinates, connectivity, offsets, element_types, gdim)` | Low-level VTK/VTU writer for flat arrays. |
| `write_level_set_vtu(filename, level_set, field_name="phi")` | Write level-set diagnostics from a `LevelSetFunction`. |
| `safe_part_name(expr)` | Convert a selection expression to a compact file-safe token. |

## Precision-specific aliases

The binding exports `float32` and `float64` variants for many native objects,
for example `MeshView_float64`, `CutCell_float64`, `QuadratureRules_float64`
and `cutcells.part.MeshPart_float64`. The unsuffixed aliases such as
`MeshView`, `LevelSetFunction` and `cutcells.part.MeshPart` are the recommended
default entry points for ordinary Python workflows. The former names
`HOCutResult` and `HOMeshPart` still work but are deprecated.
