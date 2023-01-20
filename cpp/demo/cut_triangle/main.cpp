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

int main()
{
    int gdim = 2;
    std::vector<double> ls_values = {-0.1,0.2,-0.3};

    cell::type cell_type = cell::type::triangle;
    std::vector<double> vertex_coordinates = {0.,0.,1.,0.,1.,1.};
    std::vector<std::vector<int>> bg_elements(1);
    bg_elements[0] = {0,1,2};

    cell::domain cell_domain = cell::classify_cell_domain(ls_values);

    std::cout << "Cell domain=" << domain_type_to_string(cell_domain) << std::endl;

    if(cell_domain == cell::domain::intersected)
    {
        cell::CutCell cut_cell;
        cell::cut(cell_type, vertex_coordinates, gdim, ls_values, "phi=0", cut_cell);
        std::string fname = "./results/interface.tex";
        io::write_tikz(fname,cut_cell._vertex_coords,cut_cell._connectivity,vertex_coordinates,bg_elements,ls_values,gdim);

        cell::cut(cell_type, vertex_coordinates, gdim, ls_values, "phi<0", cut_cell);
        fname = "./results/interior.tex";
        io::write_tikz(fname,cut_cell._vertex_coords,cut_cell._connectivity,vertex_coordinates,bg_elements,ls_values,gdim);
        fname = "./results/interior.vtu";
        io::write_vtk(fname,cut_cell);

        cell::cut(cell_type, vertex_coordinates, gdim, ls_values, "phi>0", cut_cell, false);
        fname = "./results/exterior.tex";
        io::write_tikz(fname,cut_cell._vertex_coords,cut_cell._connectivity,vertex_coordinates,bg_elements,ls_values,gdim);
        fname = "./results/exterior.vtu";
        io::write_vtk(fname,cut_cell);
    }
}