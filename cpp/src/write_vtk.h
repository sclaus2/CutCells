// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    LGPL-3.0-or-later
#pragma once

 
#include <iostream>
#include <fstream>
#include <string>

#include "cell_types.h"

namespace cutcells::io
{
    inline write_vtk(std::string filename, const std::span<const double> element_vertex_coords,  
                     const std::vector<std::vector<int>> elements,
                     const std::span<cell::type> element_types, 
                     const std::span<const double> ls_values, const int gdim)
    {
        std::ofstream ofs;
        ofs.open(filename.c_str(), std::ios::out );  
 
        if(ofs)
        {
            int num_points = element_vertex_coords.size()/gdim;
            int num_cells = elements.size();

            ofs << "<?xml version=\"1.0\"?>\n" 
                << "<VTKFile type=\"UnstructuredGrid\" version=\"2.2\">\n" 
                << "\t<UnstructuredGrid>\n"
                << "\t\t<Piece NumberOfPoints=\"" << num_points 
                << "\" NumberOfCells=\"" << num_cells << "\">\n";

            ofs << "\t\t\t<Points>\n" 
                <<  "\t\t\t  <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"; 
            
            //coordinates
            for(int i=0;i<num_points;i++)
            {
                double x = element_vertex_coords[i*gdim];
                double y = element_vertex_coords[i*gdim+1];
                double z = 0.0;

                if(gdim==3)
                {
                   z = element_vertex_coords[i*gdim+2];
                }

                ofs << x << " " << y << " " << z << " ";
            }
    
            ofs <<  "</DataArray>\n"; 
            ofs << "\t\t\t</Points>\n"; 
            ofs << "\t\t\t<Cells>\n";
            ofs << "\t\t\t  <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
            for(auto& element: elements)
            {
                for(auto& vertex_id: element)
                {                    
                    ofs << vertex_id << " ";

                }
            }
            ofs <<  "</DataArray>\n"; 
            ofs << "\t\t\t  <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
            int offset = 0;
            for(auto& element: elements)
            {
                offset +=element.size();
                ofs << offset << " ";
            }
            ofs <<  "</DataArray>\n"; 
            ofs << "\t\t\t  <DataArray type=\"Int8\" Name=\"types\" format=\"ascii\">";
            for(auto& type: element_types)
            {
                ofs << static_cast<int>(cell::map_cell_type_to_vtk(type)) << " ";
            }
            ofs <<  "</DataArray>\n"; 
            ofs << "\t\t\t</Cells>\n";

            
            ofs << "\t\t</Piece>\n"
                << "\t</UnstructuredGrid>\n"
                << "</VTKFile>\n";

//               <PointData Scalars=”Temperature” Vectors=”Velocity”>
//     <DataArray Name=”Velocity” .../>
//     <DataArray Name=”Temperature” .../>
//     <DataArray Name=”Pressure” .../>
//   </PointData>
        
            ofs.close();
        }
        else std:: cout << "Unable to open file " << filename << " to write" << std::endl;
    };
}