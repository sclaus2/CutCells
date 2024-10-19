#include <cutcells/cut_cell.h>
#include <cutcells/cell_types.h>
#include <cutcells/cell_flags.h>

#include <cutcells/write_tikz.h>
#include <cutcells/write_vtk.h>
#include <vector>
#include <span>
#include <iostream>
#include <string>

using namespace cutcells;

using T = double;

int main()
{
    int gdim = 2;
    std::vector<T> ls_values = {-0.1,0.2,-0.3};

    cell::type cell_type = cell::type::triangle;
    std::vector<T> vertex_coordinates = {0.,0.,1.,0.,1.,1.};
    std::vector<std::vector<int>> bg_elements(1);
    bg_elements[0] = {0,1,2};

    cell::domain cell_domain = cell::classify_cell_domain<T>(ls_values);

    std::cout << "Cell domain=" << domain_type_to_string(cell_domain) << std::endl;

    if(cell_domain == cell::domain::intersected)
    {
        cell::CutCell<T> cut_cell;
        cell::cut<T>(cell_type, vertex_coordinates, gdim, ls_values, "phi=0", cut_cell);
        std::string fname = "interface.tex";
        io::write_tikz(fname,cut_cell._vertex_coords,cut_cell._connectivity,vertex_coordinates,bg_elements,ls_values,gdim);

        cell::cut<T>(cell_type, vertex_coordinates, gdim, ls_values, "phi<0", cut_cell);
        fname = "interior.tex";
        io::write_tikz(fname,cut_cell._vertex_coords,cut_cell._connectivity,vertex_coordinates,bg_elements,ls_values,gdim);
        fname = "interior.vtu";
        io::write_vtk(fname,cut_cell);

        cell::cut<T>(cell_type, vertex_coordinates, gdim, ls_values, "phi>0", cut_cell, false);
        fname = "exterior.tex";
        io::write_tikz(fname,cut_cell._vertex_coords,cut_cell._connectivity,vertex_coordinates,bg_elements,ls_values,gdim);
        fname = "exterior.vtu";
        io::write_vtk(fname,cut_cell);
    }
}