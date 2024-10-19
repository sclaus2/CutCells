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

using T=double;

int main()
{
    std::size_t gdim = 3;
    //second order triangle test case
    std::vector<T> ls_values = {-0.1,-0.1,-0.2,-0.3, 0.4, 0.2, 0.3, 0.4, 0.2, 0.3};

    cutcells::cell::type cell_type = cutcells::cell::type::tetrahedron;
    std::vector<T> vertex_coordinates = {0.,0.,0.,
    1.,0.,0.,
    0.,1.,0.,
    0.,0.,1.,
    0., 0.5,0.5,
    0.5,0.,0.5,
    0.5,0.5,0.,
    0.,0.,0.5,
    0.,0.5,0.,
    0.5,0.,0.};

    cutcells::cell::CutCell<T> cut_cell_interface =  cutcells::cell::higher_order_cut<T>(cell_type,  vertex_coordinates,  gdim,
      ls_values, "phi=0");

    std::string fname = "interface.vtu";
    cutcells::io::write_vtk(fname,cut_cell_interface);

    cutcells::cell::CutCell<T> cut_cell_interior =  cutcells::cell::higher_order_cut<T>(cell_type,  vertex_coordinates,  gdim,
      ls_values, "phi<0");

    fname = "interior.vtu";
    cutcells::io::write_vtk(fname,cut_cell_interior);

    cutcells::cell::CutCell<T> cut_cell_exterior =  cutcells::cell::higher_order_cut<T>(cell_type,  vertex_coordinates,  gdim,
      ls_values, "phi>0");

    fname = "exterior.vtu";
    cutcells::io::write_vtk(fname,cut_cell_exterior);


}