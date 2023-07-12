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
    //second order triangle test case
    std::vector<double> ls_values = {-0.1,0.1,0.2,-0.3, 0.4, 0.2};

    cell::type cell_type = cell::type::triangle;
    std::vector<double> vertex_coordinates = {0.,0.,1.,0.,0.,1.,0.5,0.5,0.,0.5,0.5,0.0};

    cell::domain cell_domain = cell::classify_cell_domain(ls_values);

    std::cout << "Cell domain=" << domain_type_to_string(cell_domain) << std::endl;

    if(cell_domain == cell::domain::intersected)
    {
      std::size_t num_sub_cells = cutcells::cell::triangle_subdivision_table.size();
      std::cout<< "Number of sub cells=" << num_sub_cells << std::endl;

      std::vector<cell::CutCell> sub_interface_cells(num_sub_cells);
      std::vector<cell::CutCell> sub_interior_cells(num_sub_cells);
      std::vector<cell::CutCell> sub_exterior_cells(num_sub_cells);

      //Iterate over sub cells
      for(std::size_t i=0;i<num_sub_cells;i++)
      {
        auto sub_triangle = cutcells::cell::triangle_subdivision_table[i];
        std::size_t num_vertices = sub_triangle.size();

        std::vector<double> sub_ls_values(num_vertices);
        std::vector<double> sub_vertex_coordinates(num_vertices*gdim);

        for(std::size_t j=0;j<num_vertices;j++)
        {
          std::size_t vertex_id = sub_triangle[j];
          sub_ls_values[j] = ls_values[vertex_id];

          for(std::size_t k=0;k<gdim;k++)
          {
            sub_vertex_coordinates[j*gdim+k] = vertex_coordinates[vertex_id*gdim+k];
          }
        }

        cell::cut(cell_type, sub_vertex_coordinates, gdim, sub_ls_values, "phi=0", sub_interface_cells[i]);
        sub_interface_cells[i]._parent_cell_index.push_back(0);

        cell::cut(cell_type, sub_vertex_coordinates, gdim, sub_ls_values, "phi<0", sub_interior_cells[i]);
        sub_interior_cells[i]._parent_cell_index.push_back(0);

        cell::cut(cell_type, sub_vertex_coordinates, gdim, sub_ls_values, "phi>0", sub_exterior_cells[i]);
        sub_exterior_cells[i]._parent_cell_index.push_back(0);
      }

      cell::CutCell merged_interface_cell = cell::merge(sub_interface_cells);
      cell::str(merged_interface_cell);

      cell::CutCell merged_interior_cell = cell::merge(sub_interior_cells);
      cell::str(merged_interior_cell);

      cell::CutCell merged_exterior_cell = cell::merge(sub_exterior_cells);
      cell::str(merged_exterior_cell);

      std::string fname = "interface.vtu";
      io::write_vtk(fname,merged_interface_cell);

      fname = "interior.vtu";
      io::write_vtk(fname,merged_interior_cell);

      fname = "exterior.vtu";
      io::write_vtk(fname,merged_exterior_cell);

    }
}