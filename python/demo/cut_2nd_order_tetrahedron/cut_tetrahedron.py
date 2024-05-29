import cutcells
import numpy as np
import pyvista as pv

subdivision = np.array([[2,3], [1,3], [1,2], [0,3], [0,2], [0,1]])

#e01 : 4 , e12: 5, e02: 6, e03: 7, e13: 8, e23: 9
#subdivision_vtk = np.array([[0,1], [1,2], [0,2], [0,3], [1,3], [2,3]])

def level_set(x, c, r):
  value = 0
  for i in range(0,3):
    value = value  + (x[i]-c[i])**2
  value = np.sqrt(value)-r
  return value

c = np.array([0,0,0])
r = 0.6

tetra_vertices = np.array([[0.0, 0.0, 0.0],
  [1.0, 0.0, 0.0],
  [0.0, 1.0, 0.0],
  [0.0, 0.0, 1.0] ])

mid_points = np.empty([6, 3])

#obtain edge nodes with vtk ordering
for i in range(0,6):
  mid_points[i] = tetra_vertices[subdivision[i,0]] + (tetra_vertices[subdivision[i,1]]-tetra_vertices[subdivision[i,0]])/2.0

# add mid-edge nodes to tetrahedral nodes
tetra_nodes = np.concatenate((tetra_vertices, mid_points), axis=0)
print(tetra_nodes)

ls_values=np.empty(10)
idx=0
for x in tetra_nodes:
  ls_values[idx] = level_set(x, c, r)
  idx = idx+1

print(ls_values)

cell_type = cutcells.CellType.tetrahedron
triangulate = True
gdim = 3

cut_cell_int = cutcells.higher_order_cut(cell_type, tetra_nodes,  gdim, ls_values, "phi<0", triangulate)
print(cut_cell_int.str())
cut_cell_int.write_vtk("interior.vtu")

cut_cell_ext = cutcells.higher_order_cut(cell_type, tetra_nodes,  gdim, ls_values, "phi>0", triangulate)
print(cut_cell_ext.str())
cut_cell_ext.write_vtk("exterior.vtu")

#pv.OFF_SCREEN = True
#pv.start_xvfb()

grid_int = pv.UnstructuredGrid(cut_cell_int.connectivity, cut_cell_int.types, cut_cell_int.vertex_coords)
grid_ext = pv.UnstructuredGrid(cut_cell_ext.connectivity, cut_cell_ext.types, cut_cell_ext.vertex_coords)

plotter = pv.Plotter() #off_screen=True
plotter.set_background('white', top='white')
plotter.add_mesh(grid_int, color="blue",show_edges=True)
plotter.add_mesh(grid_ext, color="red",show_edges=True)
plotter.show() #screenshot='cut_triangle.png'
