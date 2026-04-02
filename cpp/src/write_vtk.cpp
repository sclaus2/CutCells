// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include <iostream>
#include <fstream>
#include <cstdint>
#include <stdexcept>
#include <string>

#include "write_vtk.h"
#include "cell_types.h"


namespace cutcells::io
{
    namespace
    {
        constexpr int VTK_LAGRANGE_CURVE = 68;
        constexpr int VTK_LAGRANGE_TRIANGLE = 69;
        constexpr int VTK_LAGRANGE_QUADRILATERAL = 70;
        constexpr int VTK_LAGRANGE_TETRAHEDRON = 71;
        constexpr int VTK_LAGRANGE_HEXAHEDRON = 72;
        constexpr int VTK_LAGRANGE_WEDGE = 73;
        constexpr int VTK_LAGRANGE_PYRAMID = 74;

        int linear_node_count(cell::type ct)
        {
            return cell::get_num_vertices(ct);
        }

        int lagrange_node_count(cell::type ct, int order)
        {
            if (order < 1)
                throw std::invalid_argument("write_vtu: geometry order must be positive");

            switch (ct)
            {
                case cell::type::interval:
                    return order + 1;
                case cell::type::triangle:
                    return (order + 1) * (order + 2) / 2;
                case cell::type::quadrilateral:
                    return (order + 1) * (order + 1);
                case cell::type::tetrahedron:
                    return (order + 1) * (order + 2) * (order + 3) / 6;
                case cell::type::hexahedron:
                    return (order + 1) * (order + 1) * (order + 1);
                case cell::type::prism:
                    return (order + 1) * (order + 1) * (order + 2) / 2;
                case cell::type::pyramid:
                {
                    int count = 0;
                    for (int k = 0; k <= order; ++k)
                    {
                        const int layer = order - k + 1;
                        count += layer * layer;
                    }
                    return count;
                }
                default:
                    return linear_node_count(ct);
            }
        }

        int vtk_type_for_curved_cell(cell::type ct, int width, int geom_order)
        {
            const int linear_n = linear_node_count(ct);
            if (width == linear_n)
                return static_cast<int>(cell::map_cell_type_to_vtk(ct));

            const int lagrange_n = lagrange_node_count(ct, geom_order);
            if (width != lagrange_n)
            {
                throw std::invalid_argument(
                    "write_vtu: cell connectivity width does not match either linear or VTK Lagrange layout");
            }

            switch (ct)
            {
                case cell::type::interval:
                    return VTK_LAGRANGE_CURVE;
                case cell::type::triangle:
                    return VTK_LAGRANGE_TRIANGLE;
                case cell::type::quadrilateral:
                    return VTK_LAGRANGE_QUADRILATERAL;
                case cell::type::tetrahedron:
                    return VTK_LAGRANGE_TETRAHEDRON;
                case cell::type::hexahedron:
                    return VTK_LAGRANGE_HEXAHEDRON;
                case cell::type::prism:
                    return VTK_LAGRANGE_WEDGE;
                case cell::type::pyramid:
                    return VTK_LAGRANGE_PYRAMID;
                default:
                    throw std::invalid_argument("write_vtu: unsupported cell type for Lagrange export");
            }
        }
    }

    void write_vtk(std::string filename, const std::span<const double> element_vertex_coords,
                     const std::span<const int> connectivity,
                     const std::span<const int> offsets,
                     const std::span<cell::type> element_types,
                     const int gdim)
    {
        std::ofstream ofs;
        ofs.open(filename.c_str(), std::ios::out );

        if(ofs)
        {
            int num_points = element_vertex_coords.size()/gdim;
            int num_cells = static_cast<int>(offsets.size()) - 1;

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
            for(auto& vertex_id: connectivity)
                ofs << vertex_id << " ";
            ofs <<  "</DataArray>\n"; 
            ofs << "\t\t\t  <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
            for(std::size_t i=1;i<offsets.size();i++)
                ofs << offsets[i] << " ";
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

            ofs.close();
        }
        else std:: cout << "Unable to open file " << filename << " to write" << std::endl;
    }

    void write_vtk(std::string filename, cell::CutCell<double>& cut_cell)
    {
        //std::as_const()
        write_vtk(filename, cut_cell._vertex_coords, cut_cell._connectivity, cut_cell._offset,
                     cut_cell._types,
                     cut_cell._gdim);
    }

    template <std::floating_point T>
    void write_vtu(std::string filename, const mesh::CurvedGlobalMesh<T>& mesh)
    {
        if (mesh.gdim <= 0)
            throw std::invalid_argument("write_vtu: mesh.gdim must be positive");
        if (mesh.offsets.empty())
            throw std::invalid_argument("write_vtu: mesh.offsets must contain at least one entry");
        if (static_cast<int>(mesh.cell_types.size()) != mesh.n_cells())
            throw std::invalid_argument("write_vtu: cell_types size does not match offsets");
        if (mesh.offsets.front() != 0)
            throw std::invalid_argument("write_vtu: offsets must start at zero");
        if (mesh.offsets.back() != static_cast<int32_t>(mesh.connectivity.size()))
            throw std::invalid_argument("write_vtu: offsets do not match connectivity size");
        if (mesh.vertex_coords.size() % static_cast<std::size_t>(mesh.gdim) != 0)
            throw std::invalid_argument("write_vtu: vertex_coords size is not divisible by gdim");

        std::ofstream ofs(filename.c_str(), std::ios::out);
        if (!ofs)
        {
            throw std::runtime_error("write_vtu: unable to open file " + filename);
        }

        const int num_points = mesh.n_vertices();
        const int num_cells = mesh.n_cells();

        ofs << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"UnstructuredGrid\" version=\"2.2\">\n"
            << "\t<UnstructuredGrid>\n"
            << "\t\t<Piece NumberOfPoints=\"" << num_points
            << "\" NumberOfCells=\"" << num_cells << "\">\n";

        ofs << "\t\t\t<Points>\n"
            << "\t\t\t  <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
        for (int i = 0; i < num_points; ++i)
        {
            const std::size_t base = static_cast<std::size_t>(i * mesh.gdim);
            const double x = static_cast<double>(mesh.vertex_coords[base]);
            const double y = static_cast<double>(mesh.vertex_coords[base + 1]);
            const double z = (mesh.gdim == 3)
                ? static_cast<double>(mesh.vertex_coords[base + 2]) : 0.0;
            ofs << x << " " << y << " " << z << " ";
        }
        ofs << "</DataArray>\n";
        ofs << "\t\t\t</Points>\n";

        ofs << "\t\t\t<Cells>\n";
        ofs << "\t\t\t  <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
        for (const int32_t vid : mesh.connectivity)
            ofs << vid << " ";
        ofs << "</DataArray>\n";

        ofs << "\t\t\t  <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
        for (int c = 0; c < num_cells; ++c)
            ofs << mesh.offsets[static_cast<std::size_t>(c + 1)] << " ";
        ofs << "</DataArray>\n";

        ofs << "\t\t\t  <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
        for (int c = 0; c < num_cells; ++c)
        {
            const int start = mesh.offsets[static_cast<std::size_t>(c)];
            const int end = mesh.offsets[static_cast<std::size_t>(c + 1)];
            const int width = end - start;
            const int vtk_type = vtk_type_for_curved_cell(
                mesh.cell_types[static_cast<std::size_t>(c)],
                width,
                mesh.geom_order);
            ofs << vtk_type << " ";
        }
        ofs << "</DataArray>\n";
        ofs << "\t\t\t</Cells>\n";

        ofs << "\t\t</Piece>\n"
            << "\t</UnstructuredGrid>\n"
            << "</VTKFile>\n";
    }

    template void write_vtu<float>(std::string filename, const mesh::CurvedGlobalMesh<float>& mesh);
    template void write_vtu<double>(std::string filename, const mesh::CurvedGlobalMesh<double>& mesh);

}
