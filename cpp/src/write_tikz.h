// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

 
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <concepts>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <span>

#include "cut_cell.h"

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

    inline void write_tikz(std::string filename,
                           const std::span<const double> element_vertex_coords,
                           const std::span<const int> connectivity,
                           const std::span<const int> offsets,
                           const std::span<const double> bg_vertex_coords,
                           const std::vector<std::vector<int>> bg_elements,
                           const std::span<const double> ls_values,
                           const int gdim)
    {
        if (offsets.empty())
        {
            write_tikz(filename, element_vertex_coords, std::vector<std::vector<int>>{},
                       bg_vertex_coords, bg_elements, ls_values, gdim);
            return;
        }

        const bool has_terminal_offset = offsets.back() == static_cast<int>(connectivity.size());
        const int num_cells = has_terminal_offset
                                  ? static_cast<int>(offsets.size()) - 1
                                  : static_cast<int>(offsets.size());

        std::vector<std::vector<int>> elements;
        elements.reserve(std::max(num_cells, 0));

        for (int i = 0; i < num_cells; ++i)
        {
            const int begin = offsets[i];
            const int end = (i + 1 < static_cast<int>(offsets.size()))
                                ? offsets[i + 1]
                                : static_cast<int>(connectivity.size());

            if (begin < 0 || end < begin || end > static_cast<int>(connectivity.size()))
                throw std::runtime_error("Invalid CSR layout in write_tikz");

            elements.emplace_back(connectivity.begin() + begin, connectivity.begin() + end);
        }

        write_tikz(filename, element_vertex_coords, elements, bg_vertex_coords, bg_elements, ls_values, gdim);
    }

    inline void write_tikz(std::string filename,
                           const cell::CutCell<double>& cut_cell,
                           const std::span<const double> bg_vertex_coords,
                           const std::vector<std::vector<int>> bg_elements,
                           const std::span<const double> ls_values,
                           const int gdim)
    {
        write_tikz(filename,
                   std::span<const double>(cut_cell._vertex_coords.data(), cut_cell._vertex_coords.size()),
                   std::span<const int>(cut_cell._connectivity.data(), cut_cell._connectivity.size()),
                   std::span<const int>(cut_cell._offset.data(), cut_cell._offset.size()),
                   bg_vertex_coords,
                   bg_elements,
                   ls_values,
                   gdim);
    }
}