// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cutcells::proto
{

/// Quadrature points as a VTK point cloud (ASCII .vtu), for ParaView.
///
/// @param points   physical coordinates, flat, 3 per point
/// @param weights  physical weights, one per point
/// @param cells    background cell of each point
inline void write_point_cloud(const std::string& path, const std::vector<double>& points,
                              const std::vector<double>& weights, const std::vector<std::int32_t>& cells)
{
    const std::size_t n = weights.size();
    if (points.size() != 3 * n || cells.size() != n)
        throw std::runtime_error("write_point_cloud: inconsistent array sizes");
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("write_point_cloud: cannot open " + path);
    out.precision(17);
    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "<UnstructuredGrid>\n"
        << "<Piece NumberOfPoints=\"" << n << "\" NumberOfCells=\"" << n << "\">\n"
        << "<PointData Scalars=\"weight\">\n"
        << "<DataArray type=\"Float64\" Name=\"weight\" format=\"ascii\">\n";
    for (double w : weights)
        out << w << '\n';
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"cell\" format=\"ascii\">\n";
    for (std::int32_t c : cells)
        out << c << '\n';
    out << "</DataArray>\n</PointData>\n<Points>\n"
        << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i)
        out << points[3 * i] << ' ' << points[3 * i + 1] << ' ' << points[3 * i + 2] << '\n';
    out << "</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i)
        out << i << '\n';
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i)
        out << i + 1 << '\n';
    out << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i)
        out << "1\n"; // VTK_VERTEX
    out << "</DataArray>\n</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
}

} // namespace cutcells::proto
