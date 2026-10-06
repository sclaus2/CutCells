# Cut results and selections

`cut(mesh, level_set)` returns a `CutResult`. Index the result with a
condition built from the level set to obtain a `MeshPart`:

```python
result = cutcells.cut(mesh, phi)

inside = result[phi < 0]
outside = result[phi > 0]
interface = result[phi == 0]

result["phi < 0"]   # the same as inside, as a string
```

`inside`, `outside` and `interface` are selected mesh parts. They expose
counts, cell IDs, quadrature rules, visualisation meshes and VTU output.

## Conditions

A level set compared with `0` is a condition, and `&` (and) and `|` (or)
combine conditions, as boolean masks do in pandas:

```python
result = cutcells.cut(mesh, [body, fluid], names=["body", "fluid"])

solid = result[body < 0]
around = result[(body > 0) & (fluid < 0)]
surface = result[body == 0]
edge = result[(body == 0) & (fluid == 0)]
```

- The clauses are `< 0`, `> 0` and `== 0`; a number other than 0 raises.
- `&` binds more tightly than `<`, so each comparison needs parentheses.
- There is no `~`: `body > 0` is the opposite of `body < 0`. A condition has no
  truth value, so `and`, `or` and `not` raise.
- `<=`, `>=` and `!=` give no conditions; `<=` and `>=` remain in strings.
- The result names the level sets of a condition as `cut` named them, matched
  by identity. A level set that was not given to `cut` raises.
- Level sets still compare and hash by identity with each other, so they work
  in lists and sets and as dictionary keys.

`str(condition)` writes the selection string with the level sets' own names,
for example `"body > 0 and fluid < 0"`. An analytic level set without a name
(`name=`) prints as a placeholder. The type is `cutcells.Condition`.

## Strings

Selection strings are the shorthand, as `df.query` is in pandas. They use the
names of the level sets and comparisons against zero:

```text
expression: term ('or' term)*
term:       factor ('and' factor)*
factor:     clause | '(' expression ')'
clause:     name ('<' | '<=' | '>' | '>=' | '=' | '==') '0'
```

Examples:

```python
result["phi < 0"]
result["phi = 0"]
result["phi == 0"]
result["sphere < 0 and plane > 0"]
result["sphere = 0 and plane < 0"]
result["sphere < 0 or plane < 0"]
```

The parser is case-insensitive for the `and` and `or` keywords; `and` binds
more tightly than `or`, and parentheses group, as in
`"(sphere < 0 or plane < 0) and box > 0"`. A name is a word without spaces,
parentheses, `<`, `>` or `=`. The right-hand side must be exactly `0`. The
`<=` and `>=` operators select the same sides as `<` and `>`, and `==` the same
as `=`.

## Level sets made of level sets

A level set made of level sets (see [Level sets](level-sets.md)) is selected
like any other. Its selections expand into selections of its level sets. For
`ring` with inside `F`:

- `ring < 0` is `F`;
- `ring > 0` is the opposite of `F`: `and` and `or` swap, and so do `< 0` and
  `> 0`;
- `ring == 0` is its faces: for each of its level sets `a`, the part of
  `a = 0` where `F` with `a` negative differs from `F` with `a` positive.

```python
ring = cutcells.LevelSetFunction(level_sets=[ball, cyl],
                                 inside=(ball < 0) & (cyl > 0), name="ring")
result = cutcells.cut(mesh, [ring, fluid])

solid = result[ring < 0]                    # ring.ball < 0 and ring.cyl > 0
around = result[(ring > 0) & (fluid < 0)]   # (ring.ball > 0 or ring.cyl < 0) and fluid < 0
surface = result[ring == 0]                 # the sphere zone and the inner cylinder
```

A level set used with both signs gives no face where the inside holds on both
of its sides. In two boxes side by side whose common face is one level set
`p`, the inside of their union holds with `p` negative and with `p` positive
alike, so that face is not part of the union's surface.

**Faces:** `part.on(*names)` keeps the faces of a surface on the level sets
named. Several names keep several faces, and the name of a level set made of
level sets keeps all its faces:

```python
hole = surface.on("cyl")             # the inner cylinder
zone = surface.on("ring.ball")       # the sphere zone, by its full name
both = surface.on("ball", "cyl")     # the whole surface
```

**Names:** in `result.level_set_names`, the level sets of `ring` have full
names (`ring.ball`, `ring.cyl`), and strings may use them. Short names
(`"cyl < 0"`) work wherever only one level set has that short name. In
`.on(...)`, a short name resolves among the part's faces first.
`result.insides` maps each level set made of level sets to its inside.

Conditions on its level sets select them by their full names:
`result[cyl < 0]` is `result["ring.cyl < 0"]`, at any depth. A level set given
to two level sets made of level sets of one cut raises; select it with a string
then.

## Dimension inference

The selected entity dimension is inferred from equality clauses:

- No equality clause selects volume cells of dimension `tdim`.
- One equality clause selects an interface of dimension `tdim - 1`.
- Two equality clauses select an intersection of dimension `tdim - 2`.

For example, in a tetrahedral mesh:

```python
volume = result[(sphere < 0) & (plane < 0)]     # dim 3
disk = result[(plane == 0) & (sphere < 0)]      # dim 2
curve = result[(plane == 0) & (sphere == 0)]    # dim 1
```

All `or` terms of one selection (joined by `|`, or by `or` in a string) must
infer the same dimension. Curves get line rules, and points in 2D get unit
weights: `curve.quadrature(order=4)`.

## Pieces on facets, edges and vertices

A selection of lower dimension can hold whole facets, edges or vertices of
the mesh: the interface `plane == 0` where the plane runs through grid faces,
the curve `(a == 0) & (b == 0)` where both vanish along mesh edges, the point
where two lines cross at a mesh vertex in 2D. A term holds such an entity
whole when

- every level set of its `== 0` clauses vanishes on the entity and changes
  sign across it: some cell sharing the entity lies on its negative side and
  some on its positive side (on the mesh boundary, one side is enough);
- none of its `< 0` or `> 0` clauses is on a level set vanishing there;
- each of its other clauses holds on a cell sharing the entity, whole.

A level set that vanishes on an entity without changing sign across it is no
surface there: `min(box1, box2)` vanishes on the face two boxes share and is
negative on both sides. `part.entities(dim)` lists the entities a part holds,
and `part.entity_quadrature(order)` gives their rules (see
[Quadrature and output](quadrature-output.md)).

## Result fields

Useful `CutResult` fields:

```python
print(result.num_cells)
print(result.num_cut_cells)
print(result.num_level_sets)
print(list(result.level_set_names))
print(result.cut_cells)       # cells that some level set cuts
print(result.domains.shape)   # (num_level_sets, num_cells): 0 inside, 1 cut, 2 outside
print(result.entities(2))     # the facets of a 3D mesh on which a level set vanishes
print(result.zero_faces)      # facets lying in a zero set, with one owning cell each
print(result.insides)         # level sets made of level sets, with their insides
print(result.backend)         # the default backend of parts selected from now on
```

`result.entities(dim)` is an `EntityTable` of the entities of dimension `dim`
(vertices, edges, facets) on which a level set vanishes. Each comes with
every cell that shares it, its number in each (Basix numbering), the level sets that
vanish on it (`vanishing`), those of them that count there (`zero`: they
change sign across it, or it lies on the mesh boundary), and each sharing
cell's side of them (`negative`, `positive`). Its arrays are CSR:

```python
facets = result.entities(2)
for e in range(facets.num_entities):
    nodes = facets.nodes[facets.node_offsets[e]:facets.node_offsets[e + 1]]
    cells = facets.cells[facets.cell_offsets[e]:facets.cell_offsets[e + 1]]
```

A level set vanishes on an entity if the bound of its restriction there (its
Bernstein coefficients restricted to the entity, Taylor bounds over the
entity itself for an analytic level set) is within 64 eps of the level set's
largest magnitude over the cells at its zero set. The same threshold decides
whether a cell beside the entity touches zero there, so that no zero set is
both on an entity and in a cell next to it. `result.zero_faces` lists the
facets in a zero set with one owning cell each, the cell on the negative side
(else the lowest).

Useful `MeshPart` fields:

```python
part = result[phi < 0]
print(part.dim)
print(part.num_cut_cells, part.num_uncut_cells)
print(part.cut_cells)
print(part.uncut_cells)
print(part.backend)
print(result[phi == 0].entities(2))   # facets it holds whole: indices into result.entities(2)
```

`num_uncut_cells` matters for volume selections: `mode="full"` includes those
uncut parent cells, while `mode="cut_only"` restricts output and quadrature to
cut cells.
