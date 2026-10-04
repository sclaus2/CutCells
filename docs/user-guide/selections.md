# Cut results and selections

`cut(mesh, level_set)` returns a `CutResult`. Index the result with a
selection expression to obtain a `MeshPart`.

```python
result = cutcells.cut(mesh, phi)

inside = result["phi < 0"]
outside = result["phi > 0"]
interface = result["phi = 0"]
```

`inside`, `outside` and `interface` are selected mesh parts. They expose
counts, cell IDs, quadrature rules, visualisation meshes and VTU output.

## Expression grammar

Selection expressions use named level sets and comparisons against zero:

```text
term:       clause ('and' clause)*
expression: term ('or' term)*
clause:     name ('<' | '<=' | '>' | '>=' | '=') '0'
```

Examples:

```python
result["phi < 0"]
result["phi = 0"]
result["sphere < 0 and plane > 0"]
result["sphere = 0 and plane < 0"]
result["sphere < 0 or plane < 0"]
```

The parser is case-insensitive for the `and` and `or` keywords. The right-hand
side must be exactly `0`. The `<=` and `>=` operators select the same sides as
`<` and `>`.

## Dimension inference

The selected entity dimension is inferred from equality clauses:

- No equality clause selects volume cells of dimension `tdim`.
- One equality clause selects an interface of dimension `tdim - 1`.
- Two equality clauses select an intersection of dimension `tdim - 2`.

For example, in a tetrahedral mesh:

```python
volume = result["sphere < 0 and plane < 0"]      # dim 3
disk = result["plane = 0 and sphere < 0"]        # dim 2
curve = result["plane = 0 and sphere = 0"]       # dim 1
```

All `or` terms in one expression must infer the same dimension. Curves get
line rules, and points in 2D get unit weights: `curve.quadrature(order=4)`.

## Result fields

Useful `CutResult` fields:

```python
print(result.num_cells)
print(result.num_cut_cells)
print(result.num_level_sets)
print(list(result.level_set_names))
print(result.cut_cells)       # cells that some level set cuts
print(result.domains.shape)   # (num_level_sets, num_cells): 0 inside, 1 cut, 2 outside
print(result.zero_faces)      # faces lying in a zero set, with their owning cells
print(result.backend)         # the default backend of parts selected from now on
```

Useful `MeshPart` fields:

```python
part = result["phi < 0"]
print(part.dim)
print(part.num_cut_cells, part.num_uncut_cells)
print(part.cut_cells)
print(part.uncut_cells)
print(part.backend)
```

`num_uncut_cells` matters for volume selections: `mode="full"` includes those
uncut parent cells, while `mode="cut_only"` restricts output and quadrature to
cut cells.
