import cutcells
import numpy as np

ls_values = np.array([0.1,-0.1,0.2, 0.4])
vertex_coordinates = np.array([1.,1.,1., 1.,-1., -1., -1, 1., -1., -1., -1, 1])

domain_id = cutcells.classify_cell_domain(ls_values)

print ("domain_id =", domain_id)

cell_type = cutcells.CellType.tetrahedron
triangulate = False
gdim = 3

cut_cell = cutcells.cut(cell_type, vertex_coordinates,  gdim, ls_values, "phi<0", triangulate)
print(cut_cell.str())
cut_cell.write_vtk("interior.vtu")

cut_cell = cutcells.cut(cell_type, vertex_coordinates,  gdim, ls_values, "phi>0", triangulate)
print(cut_cell.str())
cut_cell.write_vtk("exterior.vtu")

cut_cell = cutcells.cut(cell_type, vertex_coordinates,  gdim, ls_values, "phi=0", triangulate)
print(cut_cell.str())
cut_cell.write_vtk("interface.vtu")