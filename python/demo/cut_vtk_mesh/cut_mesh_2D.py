import cutcells
import numpy as np
import pyvista as pv
import time

def level_set(x):
  r = 0.7
  c = np.array([0,0,0])
  value = 0
  for i in range(0,3):
    value = value  + (x[i]-c[i])**2
  value = np.sqrt(value)-r
  return value

def create_rectangle_mesh(x0,y0,x1,y1,Nx,Ny):
  x = np.linspace(x0, x1, num=Nx)
  y = np.linspace(y0,y1, num=Ny)
  #produce grid of points by tensor product
  xx, yy, zz = np.meshgrid(x, y, [0])
  points = np.c_[xx.reshape(-1), yy.reshape(-1), zz.reshape(-1)]
  poly_points = pv.PolyData(points)
  poly_mesh = poly_points.delaunay_2d()
  grid = pv.UnstructuredGrid(poly_mesh)

  return grid

def create_cut_mesh(grid):
  points = grid.points
  ls_values = np.zeros(len(points))
  j = 0
  for point in points:
    ls_values[j] = level_set(point)
    j = j+1

  cut_mesh = cutcells.cut_vtk_mesh(ls_values,points,grid.cell_connectivity,grid.offset,grid.celltypes,"phi<0")
  inside_cells = cutcells.locate_cells(ls_values,points,grid.cell_connectivity,grid.offset,grid.celltypes,"phi<0")

  pv_cut =  pv.UnstructuredGrid(cut_mesh.cells,cut_mesh.types,cut_mesh.vertex_coords)
  extract = grid.extract_cells(inside_cells)

  return extract.merge(pv_cut)


N = 22
grid = create_rectangle_mesh(-1,-1,1,1,N,N)

mesh = create_cut_mesh(grid)

mesh.plot(cpos="xy", show_edges=True)

pl = pv.Plotter()
pl.add_mesh(grid, show_edges=True, color = 'white')
pl.add_mesh(mesh, show_edges=True)
pl.camera_position = 'xy'
pl.show()


