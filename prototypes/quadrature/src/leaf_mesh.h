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

/// VTK cell types used for leaves and uncut cells.
constexpr std::uint8_t vtk_tetra = 10;
constexpr std::uint8_t vtk_hexahedron = 12;
constexpr std::uint8_t vtk_lagrange_quadrilateral = 70;
constexpr std::uint8_t vtk_lagrange_hexahedron = 72;

/// Cells for visualisation in CSR layout; Lagrange cells use VTK's node order.
struct LeafMesh
{
    std::vector<double> points;              ///< physical coordinates, 3 per node
    std::vector<std::int32_t> connectivity;  ///< node indices of all cells
    std::vector<std::int32_t> offsets = {0}; ///< size n_cells + 1, offsets[0] = 0
    std::vector<std::uint8_t> types;         ///< VTK cell type per cell
    std::vector<std::int32_t> parent;        ///< background cell per cell
    std::vector<std::int32_t> degree;        ///< polynomial degree per cell (1 for linear cells)

    int n_points() const { return static_cast<int>(points.size()) / 3; }
    int n_cells() const { return static_cast<int>(offsets.size()) - 1; }
};

/// Position of node (i, j, k), 0 <= i, j, k <= p, of a degree-p Lagrange hexahedron
/// in VTK's node order (vertices, edges, faces, interior), as in VTK 9.1 and later.
inline int vtk_lagrange_hex_index(int i, int j, int k, int p)
{
    const bool ib = i == 0 || i == p, jb = j == 0 || j == p, kb = k == 0 || k == p;
    const int nb = int(ib) + int(jb) + int(kb);
    if (nb == 3)
        return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
    int offset = 8;
    if (nb == 2)
    {
        if (!ib)
            return (i - 1) + (j ? 2 * (p - 1) : 0) + (k ? 4 * (p - 1) : 0) + offset;
        if (!jb)
            return (j - 1) + (i ? (p - 1) : 3 * (p - 1)) + (k ? 4 * (p - 1) : 0) + offset;
        offset += 8 * (p - 1);
        return (k - 1) + (p - 1) * (i ? (j ? 2 : 1) : (j ? 3 : 0)) + offset;
    }
    offset += 12 * (p - 1);
    const int f = (p - 1) * (p - 1);
    if (nb == 1)
    {
        if (ib)
            return (j - 1) + (p - 1) * (k - 1) + (i ? f : 0) + offset;
        offset += 2 * f;
        if (jb)
            return (i - 1) + (p - 1) * (k - 1) + (j ? f : 0) + offset;
        offset += 2 * f;
        return (i - 1) + (p - 1) * (j - 1) + (k ? f : 0) + offset;
    }
    offset += 6 * f;
    return offset + (i - 1) + (p - 1) * ((j - 1) + (p - 1) * (k - 1));
}

/// Position of node (i, j) of a degree-p Lagrange quadrilateral in VTK's node order.
inline int vtk_lagrange_quad_index(int i, int j, int p)
{
    const bool ib = i == 0 || i == p, jb = j == 0 || j == p;
    if (ib && jb)
        return i ? (j ? 2 : 1) : (j ? 3 : 0);
    const int offset = 4;
    if (!ib && jb)
        return (i - 1) + (j ? 2 * (p - 1) : 0) + offset;
    if (ib && !jb)
        return (j - 1) + (i ? (p - 1) : 3 * (p - 1)) + offset;
    return offset + 4 * (p - 1) + (i - 1) + (p - 1) * (j - 1);
}

/// Write as an ASCII .vtu. File version 2.2 tells VTK the higher-order node order is
/// the one of VTK 9.1 and later.
inline void write_leaf_mesh(const std::string& path, const LeafMesh& mesh)
{
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("write_leaf_mesh: cannot open " + path);
    out.precision(17);
    const int np = mesh.n_points(), nc = mesh.n_cells();
    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"2.2\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n"
        << "<Piece NumberOfPoints=\"" << np << "\" NumberOfCells=\"" << nc << "\">\n"
        << "<CellData HigherOrderDegrees=\"HigherOrderDegrees\">\n"
        << "<DataArray type=\"Int32\" Name=\"parent\" format=\"ascii\">\n";
    for (std::int32_t v : mesh.parent)
        out << v << '\n';
    out << "</DataArray>\n<DataArray type=\"Float64\" Name=\"HigherOrderDegrees\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int c = 0; c < nc; ++c)
    {
        const int d = mesh.degree[c];
        out << d << ' ' << d << ' ' << (mesh.types[c] == vtk_lagrange_quadrilateral ? 0 : d) << '\n';
    }
    out << "</DataArray>\n</CellData>\n<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < np; ++i)
        out << mesh.points[3 * i] << ' ' << mesh.points[3 * i + 1] << ' ' << mesh.points[3 * i + 2] << '\n';
    out << "</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (std::int32_t v : mesh.connectivity)
        out << v << '\n';
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int c = 1; c <= nc; ++c)
        out << mesh.offsets[c] << '\n';
    out << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::uint8_t t : mesh.types)
        out << int(t) << '\n';
    out << "</DataArray>\n</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
}

} // namespace cutcells::proto
