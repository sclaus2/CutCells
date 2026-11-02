# Convergence

This tutorial follows `python/demo/demo_convergence.py`. It measures how the
area and the perimeter of a disk converge under mesh refinement, for three
ways of describing and integrating the disk.

```{raw} html
<figure class="tutorial-figure">
  <img class="tutorial-image" src="../_static/tutorials/convergence-plot.png" alt="Convergence of area and perimeter errors">
  <figcaption>Area and perimeter errors against the mesh size h, on triangular meshes of [-1, 1]^2.</figcaption>
</figure>
```

## Three Descriptions

The disk is the one of the previous chapter, $\{|x-c|<R\}$ with $R=0.6$, on
triangular meshes of size $h=2/(n-1)$.

`P1 + lut`
: A P1 level set, integrated on the straight pieces of the lookup tables. The
  pieces are polygons inscribed in the circle; the error decays as $h^2$.

`P3 + quadrays`
: A P3 level set, integrated by quadrays on the curved pieces. quadrays
  integrates the region of the interpolant to the quadrature order, so the
  error is the interpolation error of the level set, $h^{k+1}=h^4$.

`exact + quadrays`
: The analytic circle `cutcells.analytic_sphere`, integrated by quadrays. The
  geometry is exact, only the quadrature error remains, and it falls to
  round-off.

## Implementation

The three level sets differ only in how they are made and in the backend
passed to `cut`:

```python
def level_set(mesh, method):
    if method == "P1 + lut":
        return cutcells.create_level_set(mesh, circle, degree=1, name="phi"), "lut"
    if method == "P3 + quadrays":
        return cutcells.create_level_set(mesh, circle, degree=3, name="phi"), "quadrays"
    exact = cutcells.analytic_sphere([CENTER[0], CENTER[1], 0.0], RADIUS, name="phi")
    return exact, "quadrays"
```

and the rest of the study is the same for all of them:

```python
for n in sizes:
    mesh = cutcells.structured_triangle_mesh_view(-1.0, -1.0, 1.0, 1.0, n, n)
    for method in METHODS:
        phi, backend = level_set(mesh, method)
        result = cutcells.cut(mesh, phi, backend=backend)
        area = result[phi < 0].quadrature(order=order).weights.sum()
        perimeter = result[phi == 0].quadrature(order=order).weights.sum()
```

An analytic sphere in 2D is a disk; its centre takes three coordinates, of
which the third is ignored.

## Results

```text
area error (rate)
        h                  P1 + lut             P3 + quadrays          exact + quadrays
   0.2500         3.271e-02                 1.721e-05                 5.087e-08
   0.1250         8.127e-03 ( 2.01)         9.477e-07 ( 4.18)         3.029e-10 ( 7.39)
   0.0625         2.085e-03 ( 1.96)         2.737e-08 ( 5.11)         2.314e-12 ( 7.03)
   0.0312         5.136e-04 ( 2.02)         1.330e-09 ( 4.36)         6.439e-15 ( 8.49)

perimeter error (rate)
        h                  P1 + lut             P3 + quadrays          exact + quadrays
   0.2500         3.170e-02                 3.004e-05                 1.612e-06
   0.1250         7.739e-03 ( 2.03)         1.570e-06 ( 4.26)         1.026e-08 ( 7.30)
   0.0625         1.921e-03 ( 2.01)         4.526e-08 ( 5.12)         8.548e-11 ( 6.91)
   0.0312         4.796e-04 ( 2.00)         2.210e-09 ( 4.36)         1.696e-13 ( 8.98)
```

The straight pieces converge with rate 2, the interpolated level set with
rate 4, and exact geometry reaches round-off on the third mesh.

:::{admonition} Reading the rates
:class: practice

The rate of the interpolated level set fluctuates around $k+1$ on single
steps: its error comes from the interpolation error along the curve, which
partly cancels on coarse meshes. Over a longer range the rate settles: for
P3, from $n=17$ to $n=513$ the errors are 9.5e-7, 2.7e-8, 1.3e-9, 1.9e-10,
1.3e-11 and 7.8e-13, rate 4.0 on average and 4.04 on the last step. Changing
the quadrature order, the points per segment or the quadrays margin does not
change these errors: they belong to the interpolant, not to the quadrature.
:::

:::{admonition} Curved pieces in 3D
:class: practice

On tetrahedra and hexahedra, quadrays chooses its Gauss-Legendre points per
segment from `order`, which is exact on flat pieces. On curved pieces of
coarse meshes that can leave a quadrature error that does not shrink with
$h$ at first. Set the points per segment for accuracy studies:

```python
options = cutcells.QuadraysOptions()
options.points_per_segment = 12
rules = part.quadrature(order=6, options=options)
```
:::

## Run the Demo

```bash
python python/demo/demo_convergence.py
python python/demo/demo_convergence.py --sizes 17 33 65 129 257 --no-plot
```

Options: `--sizes` grid points per side of the meshes, `--order` of the
quadrature, `--no-plot` to skip the matplotlib window and
`--screenshot file.png` to save the plot.
