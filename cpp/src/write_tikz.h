// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

 
#include <iostream>
#include <fstream>
#include <string>

namespace cutcells::io
{
    inline void write_tikz(std::string filename, const std::span<const double> element_vertex_coords,  const std::vector<std::vector<int>> elements,
    const std::span<const double> bg_vertex_coords,  const std::vector<std::vector<int>> bg_elements, const std::span<const double> ls_values, const int gdim)
    {
        std::ofstream ofs;
        ofs.open(filename.c_str(), std::ios::out );  

        double scale = 10;
 
        if(ofs)
        {
            ofs << "\\documentclass[tikz]{standalone}\n" 
                << "\\usepackage{pgflibraryshapes}\n" 
                << "\\usetikzlibrary{backgrounds}\n"
                << "\\usetikzlibrary{arrows}\n"
                << "\\begin{document}\n";
            ofs << "\\begin{tikzpicture}[scale =" << scale << "]\n" << std::endl;

                        //add coordinates for element
            for(auto& element: bg_elements)
            {
                for(auto& vertex_id: element)
                {
                    double x = bg_vertex_coords[vertex_id*gdim];
                    double y = bg_vertex_coords[vertex_id*gdim+1];
                    double z = 0.;

                    if(gdim==3)
                    {
                        z = bg_vertex_coords[vertex_id*gdim+2];
                    }

                    if(gdim==2)
                    {
                        ofs << "\\node[coordinate] (B" << vertex_id << ") at (" << x << "," << y << ") {};" << std::endl;
                    }
                    else
                    {
                        ofs << "\\node[coordinate] (B" << vertex_id << ") at (" << x << "," << y << ","<< z << ") {};" << std::endl;
                    }

                    
                }
            }

            for(auto& element: bg_elements)
            {
                ofs << "\\draw";
                for(auto& vertex_id: element)
                {                    
                    ofs << "(B" << vertex_id << ") -- ";

                }
                ofs << "(B" << element[0] << ");\n\n"; 
            }

            for(auto& element: bg_elements)
            {
                for(auto& vertex_id: element)
                {
                    double x = bg_vertex_coords[vertex_id*gdim];
                    double y = bg_vertex_coords[vertex_id*gdim+1];
                    double z = 0;

                    if(gdim==3)
                    {
                        z = bg_vertex_coords[vertex_id*gdim+2];
                    }

                    ofs << "\\node"; 
                    if(ls_values[vertex_id]<0)
                        ofs << "[circle, fill=black]"; 
                    else
                        ofs << "[circle, draw]";

                    if(gdim==2)
                    {
                        ofs << "(LS" << vertex_id << ") at (" << x << "," << y << ") {};" << std::endl;
                    }
                    else
                    {
                        ofs << "(LS" << vertex_id << ") at (" << x << "," << y <<  ","<< z << ") {};" << std::endl;
                    }

                }
            }

            ofs<<std::endl;
            int k=0;

            //add coordinates for sub element
            for(auto& element: elements)
            {
                for(auto& vertex_id: element)
                {
                    double x = element_vertex_coords[vertex_id*gdim];
                    double y = element_vertex_coords[vertex_id*gdim+1];
                    double z= 0.;

                    if(gdim==3)
                    {
                        z = element_vertex_coords[vertex_id*gdim+2];
                    }

                    if(gdim==2)
                    {
                        ofs << "\\coordinate[label=right:V" << k << "] (V" << k << ") at (" << x << "," << y << ");" << std::endl;
                    }
                    else
                    {
                        ofs << "\\coordinate[label=right:V" << k << "] (V" << k << ") at (" << x << "," << y << ","<< z << ");" << std::endl;
                    }

                    k++;
                }
            }

            k=0;
            for(auto& element: elements)
            {
                ofs << "\\draw[fill opacity=0.7,fill=green!80!blue,thick,draw=green!80!black]";
                for(auto& vertex_id: element)
                {                    
                    ofs << "(V" << k << ") -- ";
                    k++;
                }
                ofs << "(V" << 0 << ");\n\n"; 
            }

            ofs << "\\end{tikzpicture}\n\\end{document}" << std::endl;
        
            ofs.close();
        }
        else std:: cout << "Unable to open file " << filename << " to write" << std::endl;
    };
}