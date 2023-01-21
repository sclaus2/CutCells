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
triangulate = False
gdim = 2

cut_cell = cutcells.cut(cell_type, vertex_coordinates,  gdim, ls_values, "phi<0", triangulate)
#print(cut_cell.str())
cut_cell.write_vtk("interior.vtu")

print("vertex coordinates=",cut_cell.vertex_coords)
print("connectivity=",cut_cell.connectivity)
print("Cell type=", cut_cell.types)

cut_cell = cutcells.cut(cell_type, vertex_coordinates,  gdim, ls_values, "phi>0", triangulate)
#print(cut_cell.str())
cut_cell.write_vtk("exterior.vtu")

cells = cut_cell.connectivity
points = cut_cell.vertex_coords
celltypes = cut_cell.types

pv.OFF_SCREEN = True
grid = pv.UnstructuredGrid(cells, celltypes, points)

plotter = pv.Plotter(off_screen=True)
plotter.add_mesh(grid, color="orange")
plotter.show(screenshot='cut_triangle.png')