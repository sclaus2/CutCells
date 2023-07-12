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

int main()
{
    std::size_t gdim = 3;
    //second order triangle test case
    std::vector<double> ls_values = {-0.1,-0.1,-0.2,-0.3, 0.4, 0.2, 0.3, 0.4, 0.2, 0.3};

    cutcells::cell::type cell_type = cutcells::cell::type::tetrahedron;
    std::vector<double> vertex_coordinates = {0.,0.,0.,
    1.,0.,0.,
    0.,1.,0.,
    0.,0.,1.,
    0., 0.5,0.5,
    0.5,0.,0.5,
    0.5,0.5,0.,
    0.,0.,0.5,
    0.,0.5,0.,
    0.5,0.,0.};

    cutcells::cell::domain cell_domain = cutcells::cell::classify_cell_domain(ls_values);

    std::cout << "Cell domain=" << domain_type_to_string(cell_domain) << std::endl;

    if(cell_domain == cutcells::cell::domain::intersected)
    {
      std::size_t num_sub_cells = 0;
      switch(cell_type)
      {
        case cutcells::cell::type::triangle: {num_sub_cells = cutcells::cell::triangle_subdivision_table.size();
                                       break;}
        case cutcells::cell::type::tetrahedron: {num_sub_cells = cutcells::cell::tetrahedron_subdivision_table.size();
                                       break;}
      }

      std::cout<< "Number of sub cells=" << num_sub_cells << std::endl;

      std::vector<cutcells::cell::CutCell> sub_interface_cells;
      std::vector<cutcells::cell::CutCell> sub_interior_cells;
      std::vector<cutcells::cell::CutCell> sub_exterior_cells;

      int sub_interface_cell_id = 0;
      int sub_interior_cell_id = 0;
      int sub_exterior_cell_id = 0;

      //Iterate over sub cells
      for(std::size_t i=0;i<num_sub_cells;i++)
      {
        const std::size_t num_vertices = cutcells::cell::get_num_vertices(cell_type);
        std::span<int> sub_tet;

        switch(cell_type)
        {
          case cutcells::cell::type::triangle: {sub_tet = cutcells::cell::triangle_subdivision_table[i];
                                        break;}
          case cutcells::cell::type::tetrahedron: {sub_tet = cutcells::cell::tetrahedron_subdivision_table[i];
                                        break;}
        }

        std::vector<double> sub_ls_values(num_vertices);
        std::vector<double> sub_vertex_coordinates(num_vertices*gdim);

        for(std::size_t j=0;j<num_vertices;j++)
        {
          std::size_t vertex_id = sub_tet[j];
          sub_ls_values[j] = ls_values[vertex_id];

          for(std::size_t k=0;k<gdim;k++)
          {
            sub_vertex_coordinates[j*gdim+k] = vertex_coordinates[vertex_id*gdim+k];
          }
        }

        //Determine if sub_cell is intersected or inside or outside
        cutcells::cell::domain cell_domain = cutcells::cell::classify_cell_domain(sub_ls_values);

        switch(cell_domain)
        {
          case cutcells::cell::domain::intersected:
          {
              cutcells::cell::CutCell tmp_cut_cell;
              cutcells::cell::cut(cell_type, sub_vertex_coordinates, gdim, sub_ls_values, "phi=0", tmp_cut_cell);
              sub_interface_cells.push_back(tmp_cut_cell);
              sub_interface_cells[sub_interface_cell_id]._parent_cell_index.push_back(0);
              sub_interface_cell_id++;

              cutcells::cell::cut(cell_type, sub_vertex_coordinates, gdim, sub_ls_values, "phi<0", tmp_cut_cell);
              sub_interior_cells.push_back(tmp_cut_cell);
              sub_interior_cells[sub_interior_cell_id]._parent_cell_index.push_back(0);
              sub_interior_cell_id++;

              cutcells::cell::cut(cell_type, sub_vertex_coordinates, gdim, sub_ls_values, "phi>0", tmp_cut_cell);
              sub_exterior_cells.push_back(tmp_cut_cell);
              sub_exterior_cells[sub_exterior_cell_id]._parent_cell_index.push_back(0);
              sub_exterior_cell_id++;
              break;
          }
          case cutcells::cell::domain::inside:
          {
            cutcells::cell::CutCell tmp_cut_cell = cutcells::cell::create_cut_cell(cell_type, sub_vertex_coordinates, gdim);
            sub_interior_cells.push_back(tmp_cut_cell);
            sub_interior_cells[sub_interior_cell_id]._parent_cell_index.push_back(0);
            sub_interior_cell_id++;
            break;
          }
          case cutcells::cell::domain::outside:
          {
              cutcells::cell::CutCell tmp_cut_cell = cutcells::cell::create_cut_cell(cell_type, sub_vertex_coordinates, gdim);
              sub_exterior_cells.push_back(tmp_cut_cell);
              sub_exterior_cells[sub_exterior_cell_id]._parent_cell_index.push_back(0);
              sub_exterior_cell_id++;
              break;
          }
        }
      }

      cutcells::cell::CutCell merged_interface_cell = cutcells::cell::merge(sub_interface_cells);
      cutcells::cell::str(merged_interface_cell);

      cutcells::cell::CutCell merged_interior_cell = cutcells::cell::merge(sub_interior_cells);
      cutcells::cell::str(merged_interior_cell);

      cutcells::cell::CutCell merged_exterior_cell = cutcells::cell::merge(sub_exterior_cells);
      cutcells::cell::str(merged_exterior_cell);

      std::string fname = "interface.vtu";
      cutcells::io::write_vtk(fname,merged_interface_cell);

      fname = "interior.vtu";
      cutcells::io::write_vtk(fname,merged_interior_cell);

      fname = "exterior.vtu";
      cutcells::io::write_vtk(fname,merged_exterior_cell);

    }
}