import cutcells
import numpy as np
import pyvista as pv

def level_set(x):
  r0 = 0.6
  r  = np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])

  sum = r0 - r
  for k in range(0,5):
    xk =  r0/np.sqrt(5.0)*2*np.cos(2*k*np.pi/5)
    yk =  r0/np.sqrt(5.0)*2*np.sin(2*k*np.pi/5)
    zk =  r0/np.sqrt(5.0)
    tmp =  (x[0] - xk)*(x[0] - xk) +  (x[1] - yk)*(x[1] - yk) \
            +  (x[2] - zk)*(x[2] - zk)
    sum += 2*np.exp(-tmp/(0.04))

  for k in range(5,10):
    xk =  r0/np.sqrt(5.0)*2*np.cos((2*(k-5)-1)*np.pi/5)
    yk =  r0/np.sqrt(5.0)*2*np.sin((2*(k-5)-1)*np.pi/5)
    zk =  -r0/np.sqrt(5.0);
    tmp =  (x[0] - xk)*(x[0] - xk) +  (x[1] - yk)*(x[1] - yk) \
            +  (x[2] - zk)*(x[2] - zk)
    sum += 2*np.exp(-tmp/(0.04))

  for k in range(10,11):
    xk =  0
    yk =  0
    zk =  r0
    tmp =  (x[0] - xk)*(x[0] - xk) +  (x[1] - yk)*(x[1] - yk) \
            +  (x[2] - zk)*(x[2] - zk)
    sum += 2*np.exp(-tmp/(0.04))

  for k in range(11,12):
    xk =  0
    yk =  0
    zk =  -r0
    tmp =  (x[0] - xk)*(x[0] - xk) +  (x[1] - yk)*(x[1] - yk) \
            +  (x[2] - zk)*(x[2] - zk);
    sum += 2*np.exp(-tmp/(0.04));

  value = -1.0*sum;

  return value

def create_box_mesh(x0,y0,z0,x1,y1,z1,Nx,Ny,Nz):
  x = np.linspace(x0, x1, num=Nx)
  y = np.linspace(y0,y1, num=Ny)
  z = np.linspace(z0,z1, num=Nz)
  #produce grid of points by tensor product
  xx, yy, zz = np.meshgrid(x, y, z)
  points = np.c_[xx.reshape(-1), yy.reshape(-1), zz.reshape(-1)]
  poly_points = pv.PolyData(points)
  poly_mesh = poly_points.delaunay_3d()
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

N = 21
grid = create_box_mesh(-1,-1,-1, 1,1,1, N,N,N)
mesh = create_cut_mesh(grid)

mesh.plot(cpos="xy", show_edges=True)

pl = pv.Plotter()
pl.add_mesh(grid, show_edges=True, style = 'wireframe')
pl.add_mesh(mesh, show_edges=True)
pl.camera_position = 'xy'
pl.show()
pl.screenshot('mesh3D.png')



