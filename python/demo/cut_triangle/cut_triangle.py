import cutcells
import numpy as np
import pyvista as pv 

def cut_cell_type_to_pyvista(cut_cell_type):
    if(cut_cell_type==cutcells.CellType.triangle):
        return pv.CellType.TRIANGLE;

ls_values = np.array([0.1,-0.1,0.2])
vertex_coordinates = np.array([0.,0.,1.,0.,1.,1.])

cell_type = cutcells.CellType.triangle
triangulate = True
gdim = 2

cut_cell_int = cutcells.cut(cell_type, vertex_coordinates,  gdim, ls_values, "phi<0", triangulate)
print(cut_cell_int.str())
cut_cell_int.write_vtk("interior.vtu")
vol = cut_cell_int.volume()
print("Volume(int)=", vol)

cut_cell_ext = cutcells.cut(cell_type, vertex_coordinates,  gdim, ls_values, "phi>0", triangulate)
print(cut_cell.str())
cut_cell_ext.write_vtk("exterior.vtu")
vol2 = cut_cell_int.volume()
print("Volume(ext)=", vol2)

pv.start_xvfb()

grid_int = pv.UnstructuredGrid(cut_cell_int.connectivity, cut_cell_int.types, cut_cell_int.vertex_coords)
grid_ext = pv.UnstructuredGrid(cut_cell_ext.connectivity, cut_cell_ext.types, cut_cell_ext.vertex_coords)

plotter = pv.Plotter(off_screen=True)
plotter.set_background('white', top='white')
plotter.add_mesh(grid_int, color="blue",show_edges=True)
plotter.add_mesh(grid_ext, color="red",show_edges=True)
plotter.show(screenshot='cut_triangle.png')
