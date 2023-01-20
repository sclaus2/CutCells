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
    int gdim = 3;
    std::vector<double> ls_values = {-0.1,-0.2,0.3,0.4};

    cell::type cell_type = cell::type::tetrahedron;
    std::vector<double> vertex_coordinates = {1.,1.,1., 1.,-1., -1., -1, 1., -1., -1., -1, 1};
    //(1,1,1), (1,−1,−1), (−1,1,−1), (−1,−1,1)
    std::vector<std::vector<int>> bg_elements(1);
    bg_elements[0] = {0,1, 2, 3};

    bool triangulate = true;

    cell::domain cell_domain = cell::classify_cell_domain(ls_values);

    std::cout << "Cell domain=" << domain_type_to_string(cell_domain) << std::endl;
    std::cout << "Flag=" << cell::get_entity_flag(ls_values, false) << std::endl;

    if(cell_domain == cell::domain::intersected)
    {
        cell::CutCell cut_cell;
        cell::cut(cell_type, vertex_coordinates, gdim, ls_values, "phi=0", cut_cell, triangulate);
        std::string fname = "./results/interface.tex";
        io::write_tikz(fname,cut_cell._vertex_coords,cut_cell._connectivity,vertex_coordinates,bg_elements,ls_values,gdim);
        fname = "./results/interface.vtu";
        io::write_vtk(fname,cut_cell._vertex_coords,cut_cell._connectivity,cut_cell._types,ls_values,gdim);

        cell::cut(cell_type, vertex_coordinates, gdim, ls_values, "phi<0", cut_cell, triangulate);
        fname = "./results/interior.tex";
        io::write_tikz(fname,cut_cell._vertex_coords,cut_cell._connectivity,vertex_coordinates,bg_elements,ls_values,gdim);
        fname = "./results/interior.vtu";
        io::write_vtk(fname,cut_cell._vertex_coords,cut_cell._connectivity,cut_cell._types,ls_values,gdim);

        cell::cut(cell_type, vertex_coordinates, gdim, ls_values, "phi>0", cut_cell, triangulate);
        fname = "./results/exterior.tex";
        io::write_tikz(fname,cut_cell._vertex_coords,cut_cell._connectivity,vertex_coordinates,bg_elements,ls_values,gdim);
        fname = "./results/exterior.vtu";
        io::write_vtk(fname,cut_cell._vertex_coords,cut_cell._connectivity,cut_cell._types,ls_values,gdim);
    }
}