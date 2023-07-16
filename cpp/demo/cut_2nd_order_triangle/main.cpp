#include <cutcells/cut_cell.h>
#include <cutcells/cell_types.h>
#include <cutcells/cell_flags.h>
#include <cutcells/cell_subdivision.h>

#include <cutcells/write_tikz.h>
#include <cutcells/write_vtk.h>
#include <vector>
#include <span>
#include <iostream>
#include <string>

using namespace cutcells;

int main()
{
    std::size_t gdim = 2;
    bool triangulate = true;
    //second order triangle test case
    std::vector<double> ls_values = {-0.1,0.1,0.2,-0.3, 0.4, 0.2};

    cell::type cell_type = cell::type::triangle;
    std::vector<double> vertex_coordinates = {0.,0.,1.,0.,0.,1.,0.5,0.5,0.,0.5,0.5,0.0};

    cutcells::cell::CutCell cut_cell_interface =  cutcells::cell::higher_order_cut(cell_type,  vertex_coordinates,  gdim,
      ls_values, "phi=0",triangulate);

    std::string fname = "interface.vtu";
    cutcells::io::write_vtk(fname,cut_cell_interface);

    cutcells::cell::CutCell cut_cell_interior =  cutcells::cell::higher_order_cut(cell_type,  vertex_coordinates,  gdim,
      ls_values, "phi<0",triangulate);

    fname = "interior.vtu";
    cutcells::io::write_vtk(fname,cut_cell_interior);

    cutcells::cell::CutCell cut_cell_exterior =  cutcells::cell::higher_order_cut(cell_type,  vertex_coordinates,  gdim,
      ls_values, "phi>0",triangulate);

    fname = "exterior.vtu";
    cutcells::io::write_vtk(fname,cut_cell_exterior);
}