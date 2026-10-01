# Visualising cuts from the quadrature engine

Today `HOMeshPart` gets its visualisation mesh from AdaptCell's sub-cells
(`part.write_vtu(...)`, straight or curved). The new engine produces quadrature
points only, so it needs its own visual output. The recommendation is to draw the
cut from the engine's own decomposition: what you see is then exactly what is
integrated.

## 1. Point clouds (implemented)

`quadrature_study --vtk <prefix>` writes every rule as a VTK point cloud (`.vtu`,
VTK_VERTEX cells) with point data `weight` and `cell`, one file per generator,
part, n and q. In ParaView, colour by `weight` or glyph by it.

Use it to check that points stay inside the cell and the selected part, to see
where bisection concentrates points, and to find cells with too few points. Files
are ASCII and grow with the point count (n = 6, q = 3, certify: 36 MB for
phi < 0); binary output is the next step for larger runs.

## 2. Leaf cells of the decomposition (proposed)

The engine already parametrises every piece it integrates. A volume leaf is a
column over a column over an interval:

- x1 in [a, b], an interval of the 1D level;
- x2 in [L2(x1), U2(x1)], one segment of the level-2 height line;
- x3 in [L3(x1, x2), U3(x1, x2)], one segment of the level-3 height line.

Each bound is a box face, a clip plane or a root of a level set, and the
decomposition guarantees that the same bounds apply across the whole leaf. The
map s in [0, 1]^3 -> x is the one whose image of the Gauss-Legendre grid is the
quadrature rule. Running the same recursion with p + 1 Gauss-Lobatto nodes per
direction instead gives the nodes of a Lagrange hexahedron of order p
(VTK_LAGRANGE_HEXAHEDRON). An interface leaf is the graph of the root over a
level-2 leaf: a Lagrange quadrilateral (VTK_LAGRANGE_QUADRILATERAL). ParaView,
VTK 9 and pyvista render both as curved cells.

Why this option:

- **Consistent by construction.** The visual mesh is the integration domain; a
  gap or an overlap in the picture is a quadrature bug.
- **Parts come for free.** Each volume leaf has one sign per level set, the sign
  of its segment. A part such as `result["phi1 < 0 and phi2 > 0"]` is the set
  of leaves whose signs satisfy the selection term, the same masks the
  generators already take.
- **Curved without a mesher.** No separate curved sub-cell construction; the
  order follows the requested p.
- **Debug data per leaf**: parent cell, level-set signs, bisection depth,
  certified or not, the bounds that define it.

Caveats:

- Leaves are columns aligned with the chosen height directions. They do not
  conform across cells; that is fine for display.
- Segments shrink to zero at tangencies and where three bounds meet, giving
  collapsed cells. Display is acceptable; leaves below a volume threshold can be
  dropped.
- Typical size: 5 to 20 leaves per cut cell, 64 nodes per leaf at p = 3.

Implementation sketch: a second emitter in the same recursion. Each level
records, instead of Gauss-Legendre points, the index of the interval or segment
and its two bound functions. At the top level the leaf is evaluated on the
(p + 1)^3 Gauss-Lobatto grid and its nodes reordered into VTK's Lagrange
ordering.

## 3. Whole meshes and the Python side

Uncut cells are exported as themselves, straight or curved by their geometry map,
and cut cells as their leaves per part, in one `.vtu` (or VTKHDF) per part, as
`HOMeshPart.write_vtu()` offers today. The binding can also return
`(points, connectivity, offsets, cell_types)` arrays for
`pyvista.UnstructuredGrid`, zero-copy like the other bindings.

This answers one of the hand-off note's open questions: if the engine draws its
own leaves, AdaptCell is not needed for visualisation, and the picture can never
disagree with the quadrature.
