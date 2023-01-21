import cutcells
import numpy as np
import pyvista as pv 
import vtk

def cut_cell_type_to_pyvista(cut_cell_type):
    if(cut_cell_type==cutcells.CellType.triangle):
        return pv.CellType.TRIANGLE;

ls_values = np.array([0.1,-0.1,0.2])
vertex_coordinates = np.array([0.,0.,1.,0.,1.,1.])

domain_id = cutcells.classify_cell_domain(ls_values)

print ("domain_id =", domain_id)

cell_type = cutcells.CellType.triangle
triangulate = True
gdim = 2

cut_cell_int = cutcells.cut(cell_type, vertex_coordinates,  gdim, ls_values, "phi<0", triangulate)
#print(cut_cell.str())
cut_cell_int.write_vtk("interior.vtu")

cut_cell_ext = cutcells.cut(cell_type, vertex_coordinates,  gdim, ls_values, "phi>0", triangulate)
#print(cut_cell.str())
cut_cell_ext.write_vtk("exterior.vtu")

#add offset
# connectivity_ext = np.zeros(cut_cell_ext.connectivity.size,dtype=np.uint32)
# connectivity_ext[0] = cut_cell_ext.connectivity[0]
# for i in range(1,connectivity_ext.size):
#     connectivity_ext[i] = cut_cell_ext.connectivity[i]+cut_cell_int.connectivity[0]

# cells = np.append(cut_cell_int.connectivity, connectivity_ext,axis=0)  
# points = np.append(cut_cell_int.vertex_coords, cut_cell_ext.vertex_coords,axis=0)
# celltypes = np.append(cut_cell_int.types,cut_cell_ext.types,axis=0)

# print(cells)
# print(points)
# print(celltypes)

#pv.OFF_SCREEN = True
pv.start_xvfb()

#grid = pv.UnstructuredGrid(cells, celltypes, points)

grid_int = pv.UnstructuredGrid(cut_cell_int.connectivity, cut_cell_int.types, cut_cell_int.vertex_coords)
grid_ext = pv.UnstructuredGrid(cut_cell_ext.connectivity, cut_cell_ext.types, cut_cell_ext.vertex_coords)

plotter = pv.Plotter(off_screen=True)
plotter.add_mesh(grid_int, color="blue",show_edges=True)
plotter.add_mesh(grid_ext, color="red",show_edges=True)
plotter.show(screenshot='cut_triangle.png')