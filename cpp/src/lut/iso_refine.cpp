// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "iso_refine.h"

#include "../cell_topology.h"
#include "../reference_cell.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace cutcells
{
namespace
{

int key_ij(int i, int j)
{
    return i * 100 + j;
}

int key_ijk(int i, int j, int k)
{
    return i * 10000 + j * 100 + k;
}

double tetra_det(std::span<const double> coords,
                 std::array<int, 4> tet,
                 int tdim)
{
    auto x = [&](int v, int d) -> double
    {
        return coords[static_cast<std::size_t>(v * tdim + d)];
    };
    const double ax = x(tet[1], 0) - x(tet[0], 0);
    const double ay = x(tet[1], 1) - x(tet[0], 1);
    const double az = x(tet[1], 2) - x(tet[0], 2);
    const double bx = x(tet[2], 0) - x(tet[0], 0);
    const double by = x(tet[2], 1) - x(tet[0], 1);
    const double bz = x(tet[2], 2) - x(tet[0], 2);
    const double cx = x(tet[3], 0) - x(tet[0], 0);
    const double cy = x(tet[3], 1) - x(tet[0], 1);
    const double cz = x(tet[3], 2) - x(tet[0], 2);
    return ax * (by * cz - bz * cy)
         - ay * (bx * cz - bz * cx)
         + az * (bx * cy - by * cx);
}

void orient_tetrahedra(IsoRefineTemplate& tpl)
{
    if (tpl.tdim != 3)
        return;

    for (int c = 0; c < tpl.n_cells; ++c)
    {
        if (tpl.cell_types[static_cast<std::size_t>(c)] != cell::type::tetrahedron)
            continue;
        const auto offset = static_cast<std::size_t>(tpl.cell_offsets[static_cast<std::size_t>(c)]);
        std::array<int, 4> tet = {
            tpl.cell_connectivity[offset + 0],
            tpl.cell_connectivity[offset + 1],
            tpl.cell_connectivity[offset + 2],
            tpl.cell_connectivity[offset + 3]};
        if (tetra_det(std::span<const double>(tpl.ref_vertex_coords.data(),
                                              tpl.ref_vertex_coords.size()),
                      tet, 3) < 0.0)
            std::swap(tpl.cell_connectivity[offset + 2],
                      tpl.cell_connectivity[offset + 3]);
    }
}

IsoRefineTemplate make_template(cell::type parent_cell_type,
                                cell::type child_cell_type,
                                int tdim,
                                int vertices_per_cell,
                                std::vector<double> ref_coords,
                                std::vector<int> parent_dim,
                                std::vector<int> parent_id,
                                std::vector<int> cells)
{
    IsoRefineTemplate tpl;
    tpl.n_vertices = static_cast<int>(parent_dim.size());
    tpl.n_cells = static_cast<int>(cells.size()) / vertices_per_cell;
    tpl.tdim = tdim;
    tpl.vertices_per_cell = vertices_per_cell;
    tpl.parent_cell_type = parent_cell_type;
    tpl.child_cell_type = child_cell_type;
    tpl.ref_vertex_coords = std::move(ref_coords);
    tpl.vertex_parent_dim = std::move(parent_dim);
    tpl.vertex_parent_id = std::move(parent_id);
    tpl.cell_connectivity = std::move(cells);
    tpl.cell_types.assign(static_cast<std::size_t>(tpl.n_cells), child_cell_type);
    for (int c = 0; c <= tpl.n_cells; ++c)
        tpl.cell_offsets.push_back(c * vertices_per_cell);
    orient_tetrahedra(tpl);
    return tpl;
}

/// The parent entity of a point of a parent reference cell with planar faces:
/// (0, vertex), (1, edge), (2, face) or (3, 0), in CutCells' numbering of the
/// cell's edges and faces.
std::pair<int, int> parent_entity(cell::type parent, const double* x)
{
    const std::vector<double> ref = cell::reference_vertices<double>(parent);
    const int tdim = cell::get_tdim(parent), nv = cell::get_num_vertices(parent);
    const double tol = 1e-12;
    for (int v = 0; v < nv; ++v)
    {
        double d = 0;
        for (int i = 0; i < tdim; ++i)
            d = std::max(d, std::abs(x[i] - ref[static_cast<std::size_t>(v * tdim + i)]));
        if (d <= tol)
            return {0, v};
    }
    // the faces whose planes hold x
    std::vector<int> on;
    for (int f = 0; f < cell::num_faces(parent); ++f)
    {
        const std::span<const int> fv = cell::face_vertices(parent, f);
        std::array<std::array<double, 3>, 3> p{};
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i)
                p[j][i] = ref[static_cast<std::size_t>(fv[j] * tdim + i)];
        const std::array<double, 3> a = {p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]},
                                    b = {p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]};
        const std::array<double, 3> n = {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
        double dot = 0;
        for (int i = 0; i < 3; ++i)
            dot += n[i] * (x[i] - p[0][i]);
        if (std::abs(dot) <= tol)
            on.push_back(f);
    }
    if (on.empty())
        return {3, 0};
    if (on.size() == 1)
        return {2, on[0]};
    // an edge: the two vertices the first two faces share
    std::vector<int> shared;
    for (const int u : cell::face_vertices(parent, on[0]))
        for (const int w : cell::face_vertices(parent, on[1]))
            if (u == w)
                shared.push_back(u);
    const std::span<const std::array<int, 2>> edges = cell::edges(parent);
    for (std::size_t e = 0; e < edges.size(); ++e)
        if (shared.size() == 2
            && ((edges[e][0] == shared[0] && edges[e][1] == shared[1])
                || (edges[e][0] == shared[1] && edges[e][1] == shared[0])))
            return {1, static_cast<int>(e)};
    throw std::logic_error("iso_refine: a point on two faces but no edge");
}

/// A template from its vertices (each with its parent entity) and children of
/// mixed types.
IsoRefineTemplate make_mixed_template(cell::type parent_cell_type, std::vector<double> ref_coords,
                                      std::vector<cell::type> types, std::vector<int> cells)
{
    IsoRefineTemplate tpl;
    const int tdim = cell::get_tdim(parent_cell_type);
    tpl.tdim = tdim;
    tpl.n_vertices = static_cast<int>(ref_coords.size()) / tdim;
    tpl.parent_cell_type = parent_cell_type;
    tpl.ref_vertex_coords = std::move(ref_coords);
    for (int v = 0; v < tpl.n_vertices; ++v)
    {
        const auto [dim, id] = parent_entity(parent_cell_type, tpl.ref_vertex_coords.data() + v * tdim);
        tpl.vertex_parent_dim.push_back(dim);
        tpl.vertex_parent_id.push_back(id);
    }
    tpl.n_cells = static_cast<int>(types.size());
    tpl.cell_offsets.push_back(0);
    for (const cell::type t : types)
        tpl.cell_offsets.push_back(tpl.cell_offsets.back() + cell::get_num_vertices(t));
    tpl.cell_types = std::move(types);
    tpl.cell_connectivity = std::move(cells);
    if (std::all_of(tpl.cell_types.begin(), tpl.cell_types.end(),
                    [&](cell::type t) { return t == tpl.cell_types.front(); }))
    {
        tpl.child_cell_type = tpl.cell_types.front();
        tpl.vertices_per_cell = cell::get_num_vertices(tpl.child_cell_type);
    }
    orient_tetrahedra(tpl);
    return tpl;
}

IsoRefineTemplate make_p1_storage(cell::type cell_type)
{
    const int tdim = cell::get_tdim(cell_type);
    const int n_vertices = cell::get_num_vertices(cell_type);
    std::vector<double> coords;
    switch (cell_type)
    {
    case cell::type::interval:
        coords = {0.0, 1.0};
        break;
    case cell::type::triangle:
        coords = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0};
        break;
    case cell::type::quadrilateral:
        coords = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0};
        break;
    case cell::type::tetrahedron:
        coords = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                  0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
        break;
    case cell::type::hexahedron:
        coords = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                  0.0, 1.0, 0.0, 1.0, 1.0, 0.0,
                  0.0, 0.0, 1.0, 1.0, 0.0, 1.0,
                  0.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        break;
    case cell::type::prism:
        coords = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                  0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0};
        break;
    case cell::type::pyramid:
        coords = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                  1.0, 1.0, 0.0, 0.0, 0.0, 1.0};
        break;
    default:
        throw std::invalid_argument(
            "p1_template: unsupported cell type "
            + cell::cell_type_to_str(cell_type));
    }

    std::vector<int> parent_dim(static_cast<std::size_t>(n_vertices), 0);
    std::vector<int> parent_id(static_cast<std::size_t>(n_vertices));
    std::vector<int> cells(static_cast<std::size_t>(n_vertices));
    for (int i = 0; i < n_vertices; ++i)
    {
        parent_id[static_cast<std::size_t>(i)] = i;
        cells[static_cast<std::size_t>(i)] = i;
    }

    return make_template(cell_type, cell_type, tdim, n_vertices,
                         std::move(coords), std::move(parent_dim),
                         std::move(parent_id), std::move(cells));
}

IsoRefineTemplate make_interval_iso_p1_storage(int order)
{
    std::vector<double> coords;
    std::vector<int> parent_dim;
    std::vector<int> parent_id;
    std::vector<int> cells;
    coords.reserve(static_cast<std::size_t>(order + 1));
    parent_dim.reserve(static_cast<std::size_t>(order + 1));
    parent_id.reserve(static_cast<std::size_t>(order + 1));
    cells.reserve(static_cast<std::size_t>(2 * order));

    coords.push_back(0.0);
    coords.push_back(1.0);
    parent_dim.push_back(0);
    parent_dim.push_back(0);
    parent_id.push_back(0);
    parent_id.push_back(1);

    for (int i = 1; i < order; ++i)
    {
        coords.push_back(static_cast<double>(i) / static_cast<double>(order));
        parent_dim.push_back(1);
        parent_id.push_back(0);
    }

    auto vertex = [order](int i)
    {
        if (i == 0)
            return 0;
        if (i == order)
            return 1;
        return 1 + i;
    };

    for (int i = 0; i < order; ++i)
    {
        cells.push_back(vertex(i));
        cells.push_back(vertex(i + 1));
    }

    return make_template(cell::type::interval, cell::type::interval, 1, 2,
                         std::move(coords), std::move(parent_dim),
                         std::move(parent_id), std::move(cells));
}

IsoRefineTemplate make_triangle_iso_p1_storage(int order)
{
    const double h = 1.0 / static_cast<double>(order);
    std::vector<double> coords;
    std::vector<int> parent_dim;
    std::vector<int> parent_id;
    std::vector<int> cells;
    std::unordered_map<int, int> lattice_to_vertex;

    auto add_pt = [&](int i, int j, int pd, int pid)
    {
        const int id = static_cast<int>(parent_dim.size());
        coords.push_back(i * h);
        coords.push_back(j * h);
        parent_dim.push_back(pd);
        parent_id.push_back(pid);
        lattice_to_vertex.emplace(key_ij(i, j), id);
    };

    add_pt(0, 0, 0, 0);
    add_pt(order, 0, 0, 1);
    add_pt(0, order, 0, 2);

    // Basix triangle edge order: (1,2), (0,2), (0,1).
    for (int t = 1; t < order; ++t)
        add_pt(order - t, t, 1, 0);
    for (int t = 1; t < order; ++t)
        add_pt(0, t, 1, 1);
    for (int t = 1; t < order; ++t)
        add_pt(t, 0, 1, 2);

    for (int j = 1; j < order; ++j)
        for (int i = 1; i < order - j; ++i)
            add_pt(i, j, 2, 0);

    auto vertex = [&](int i, int j)
    {
        return lattice_to_vertex.at(key_ij(i, j));
    };

    for (int i = 0; i < order; ++i)
    {
        for (int j = 0; j < order - i; ++j)
        {
            cells.push_back(vertex(i, j));
            cells.push_back(vertex(i + 1, j));
            cells.push_back(vertex(i, j + 1));

            if (i + j <= order - 2)
            {
                cells.push_back(vertex(i + 1, j));
                cells.push_back(vertex(i + 1, j + 1));
                cells.push_back(vertex(i, j + 1));
            }
        }
    }

    return make_template(cell::type::triangle, cell::type::triangle, 2, 3,
                         std::move(coords), std::move(parent_dim),
                         std::move(parent_id), std::move(cells));
}

IsoRefineTemplate make_quadrilateral_iso_p1_storage(int order)
{
    const double h = 1.0 / static_cast<double>(order);
    std::vector<double> coords;
    std::vector<int> parent_dim;
    std::vector<int> parent_id;
    std::vector<int> cells;
    std::unordered_map<int, int> map_ij;

    auto add_pt = [&](int i, int j, int pd, int pid)
    {
        const int id = static_cast<int>(parent_dim.size());
        coords.push_back(i * h);
        coords.push_back(j * h);
        parent_dim.push_back(pd);
        parent_id.push_back(pid);
        map_ij.emplace(key_ij(i, j), id);
    };

    add_pt(0, 0, 0, 0);
    add_pt(order, 0, 0, 1);
    add_pt(0, order, 0, 2);
    add_pt(order, order, 0, 3);

    // Basix quadrilateral edge order: (0,1), (0,2), (1,3), (2,3).
    for (int t = 1; t < order; ++t)
        add_pt(t, 0, 1, 0);
    for (int t = 1; t < order; ++t)
        add_pt(0, t, 1, 1);
    for (int t = 1; t < order; ++t)
        add_pt(order, t, 1, 2);
    for (int t = 1; t < order; ++t)
        add_pt(t, order, 1, 3);

    for (int j = 1; j < order; ++j)
        for (int i = 1; i < order; ++i)
            add_pt(i, j, 2, 0);

    auto vertex = [&](int i, int j)
    {
        return map_ij.at(key_ij(i, j));
    };

    for (int j = 0; j < order; ++j)
    {
        for (int i = 0; i < order; ++i)
        {
            const int v00 = vertex(i, j);
            const int v10 = vertex(i + 1, j);
            const int v01 = vertex(i, j + 1);
            const int v11 = vertex(i + 1, j + 1);
            cells.insert(cells.end(), {v00, v10, v01, v11});
        }
    }

    return make_template(cell::type::quadrilateral, cell::type::quadrilateral, 2, 4,
                         std::move(coords), std::move(parent_dim),
                         std::move(parent_id), std::move(cells));
}

IsoRefineTemplate make_tetrahedron_iso_p1_storage(int order)
{
    if (order == 2)
    {
        return make_template(
            cell::type::tetrahedron, cell::type::tetrahedron, 3, 4,
            {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0.5, 0.5,
             0.5, 0, 0.5, 0.5, 0.5, 0, 0, 0, 0.5, 0, 0.5, 0,
             0.5, 0, 0},
            {0, 0, 0, 0, 1, 1, 1, 1, 1, 1},
            {0, 1, 2, 3, 0, 1, 2, 3, 4, 5},
            {0, 8, 7, 9, 9, 6, 5, 1, 8, 2, 4, 6, 7, 4, 3, 5,
             8, 9, 6, 5, 8, 6, 4, 5, 8, 4, 7, 5, 8, 7, 9, 5});
    }

    if (order == 3)
    {
        return make_template(
            cell::type::tetrahedron, cell::type::tetrahedron, 3, 4,
            {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0.66666666666666674,
             0.33333333333333331, 0, 0.33333333333333326, 0.66666666666666674,
             0.66666666666666674, 0, 0.33333333333333331, 0.33333333333333326,
             0, 0.66666666666666674, 0.66666666666666674, 0.33333333333333331,
             0, 0.33333333333333326, 0.66666666666666674, 0, 0, 0,
             0.33333333333333331, 0, 0, 0.66666666666666674, 0,
             0.33333333333333331, 0, 0, 0.66666666666666674, 0,
             0.33333333333333331, 0, 0, 0.66666666666666674, 0, 0,
             0.33333333333333343, 0.33333333333333331, 0.33333333333333331,
             0, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331,
             0, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331, 0},
            {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2},
            {0, 1, 2, 3, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0, 1, 2, 3},
            {0, 12, 10, 14, 15, 8, 6, 1, 13, 2, 4, 9, 11, 5, 3, 7,
             14, 19, 18, 15, 19, 9, 16, 8, 12, 13, 17, 19, 10, 17, 11, 18,
             18, 16, 7, 6, 17, 4, 5, 16, 12, 14, 19, 18, 12, 19, 17, 18,
             12, 17, 10, 18, 12, 10, 14, 18, 19, 15, 8, 6, 19, 8, 16, 6,
             19, 16, 18, 6, 19, 18, 15, 6, 13, 19, 9, 16, 13, 9, 4, 16,
             13, 4, 17, 16, 13, 17, 19, 16, 17, 18, 16, 7, 17, 16, 5, 7,
             17, 5, 11, 7, 17, 11, 18, 7, 16, 17, 19, 18});
    }

    return make_template(
        cell::type::tetrahedron, cell::type::tetrahedron, 3, 4,
        {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0.75, 0.25, 0, 0.5, 0.5,
         0, 0.25, 0.75, 0.75, 0, 0.25, 0.5, 0, 0.5, 0.25, 0, 0.75,
         0.75, 0.25, 0, 0.5, 0.5, 0, 0.25, 0.75, 0, 0, 0, 0.25,
         0, 0, 0.5, 0, 0, 0.75, 0, 0.25, 0, 0, 0.5, 0, 0,
         0.75, 0, 0.25, 0, 0, 0.5, 0, 0, 0.75, 0, 0, 0.5,
         0.25, 0.25, 0.25, 0.5, 0.25, 0.25, 0.25, 0.5, 0,
         0.25, 0.25, 0, 0.5, 0.25, 0, 0.25, 0.5, 0.25, 0,
         0.25, 0.5, 0, 0.25, 0.25, 0, 0.5, 0.25, 0.25, 0,
         0.5, 0.25, 0, 0.25, 0.5, 0, 0.25, 0.25, 0.25},
        {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3},
        {0, 1, 2, 3, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5,
         0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 0},
        {0, 16, 13, 19, 21, 10, 7, 1, 18, 2, 4, 12, 15, 6, 3, 9,
         19, 31, 28, 20, 20, 32, 29, 21, 32, 11, 22, 10, 33, 12, 23, 11,
         17, 18, 26, 33, 16, 17, 25, 31, 13, 25, 14, 28, 14, 27, 15, 30,
         29, 22, 8, 7, 30, 24, 9, 8, 26, 4, 5, 23, 27, 5, 6, 24,
         28, 34, 30, 29, 34, 23, 24, 22, 25, 26, 27, 34, 31, 33, 34, 32,
         16, 19, 31, 28, 16, 31, 25, 28, 16, 25, 13, 28, 16, 13, 19, 28,
         32, 21, 10, 7, 32, 10, 22, 7, 32, 22, 29, 7, 32, 29, 21, 7,
         18, 33, 12, 23, 18, 12, 4, 23, 18, 4, 26, 23, 18, 26, 33, 23,
         27, 30, 24, 9, 27, 24, 6, 9, 27, 6, 15, 9, 27, 15, 30, 9,
         31, 20, 32, 29, 31, 32, 34, 29, 31, 34, 28, 29, 31, 28, 20, 29,
         33, 32, 11, 22, 33, 11, 23, 22, 33, 23, 34, 22, 33, 34, 32, 22,
         17, 31, 33, 34, 17, 33, 26, 34, 17, 26, 25, 34, 17, 25, 31, 34,
         25, 28, 34, 30, 25, 34, 27, 30, 25, 27, 14, 30, 25, 14, 28, 30,
         34, 29, 22, 8, 34, 22, 24, 8, 34, 24, 30, 8, 34, 30, 29, 8,
         26, 34, 23, 24, 26, 23, 5, 24, 26, 5, 27, 24, 26, 27, 34, 24,
         34, 25, 31, 28, 22, 34, 32, 29, 23, 26, 33, 34, 24, 27, 34, 30});
}

IsoRefineTemplate make_hexahedron_iso_p1_storage(int order)
{
    const double h = 1.0 / static_cast<double>(order);
    std::vector<double> coords;
    std::vector<int> parent_dim;
    std::vector<int> parent_id;
    std::vector<int> cells;
    std::unordered_map<int, int> map_ijk;

    auto add_pt = [&](int i, int j, int k, int pd, int pid)
    {
        const int id = static_cast<int>(parent_dim.size());
        coords.insert(coords.end(), {i * h, j * h, k * h});
        parent_dim.push_back(pd);
        parent_id.push_back(pid);
        map_ijk.emplace(key_ijk(i, j, k), id);
    };

    add_pt(0, 0, 0, 0, 0); add_pt(order, 0, 0, 0, 1);
    add_pt(0, order, 0, 0, 2); add_pt(order, order, 0, 0, 3);
    add_pt(0, 0, order, 0, 4); add_pt(order, 0, order, 0, 5);
    add_pt(0, order, order, 0, 6); add_pt(order, order, order, 0, 7);

    for (int t = 1; t < order; ++t) add_pt(t, 0, 0, 1, 0);
    for (int t = 1; t < order; ++t) add_pt(0, t, 0, 1, 1);
    for (int t = 1; t < order; ++t) add_pt(0, 0, t, 1, 2);
    for (int t = 1; t < order; ++t) add_pt(order, t, 0, 1, 3);
    for (int t = 1; t < order; ++t) add_pt(order, 0, t, 1, 4);
    for (int t = 1; t < order; ++t) add_pt(t, order, 0, 1, 5);
    for (int t = 1; t < order; ++t) add_pt(0, order, t, 1, 6);
    for (int t = 1; t < order; ++t) add_pt(order, order, t, 1, 7);
    for (int t = 1; t < order; ++t) add_pt(t, 0, order, 1, 8);
    for (int t = 1; t < order; ++t) add_pt(0, t, order, 1, 9);
    for (int t = 1; t < order; ++t) add_pt(order, t, order, 1, 10);
    for (int t = 1; t < order; ++t) add_pt(t, order, order, 1, 11);

    for (int j = 1; j < order; ++j) for (int i = 1; i < order; ++i) add_pt(i, j, 0, 2, 0);
    for (int k = 1; k < order; ++k) for (int i = 1; i < order; ++i) add_pt(i, 0, k, 2, 1);
    for (int k = 1; k < order; ++k) for (int j = 1; j < order; ++j) add_pt(0, j, k, 2, 2);
    for (int k = 1; k < order; ++k) for (int j = 1; j < order; ++j) add_pt(order, j, k, 2, 3);
    for (int k = 1; k < order; ++k) for (int i = 1; i < order; ++i) add_pt(i, order, k, 2, 4);
    for (int j = 1; j < order; ++j) for (int i = 1; i < order; ++i) add_pt(i, j, order, 2, 5);

    for (int k = 1; k < order; ++k)
        for (int j = 1; j < order; ++j)
            for (int i = 1; i < order; ++i)
                add_pt(i, j, k, 3, 0);

    auto vertex = [&](int i, int j, int k)
    {
        return map_ijk.at(key_ijk(i, j, k));
    };

    for (int k = 0; k < order; ++k)
        for (int j = 0; j < order; ++j)
            for (int i = 0; i < order; ++i)
                cells.insert(cells.end(),
                             {vertex(i, j, k),
                              vertex(i + 1, j, k),
                              vertex(i, j + 1, k),
                              vertex(i + 1, j + 1, k),
                              vertex(i, j, k + 1),
                              vertex(i + 1, j, k + 1),
                              vertex(i, j + 1, k + 1),
                              vertex(i + 1, j + 1, k + 1)});

    return make_template(cell::type::hexahedron, cell::type::hexahedron, 3, 8,
                         std::move(coords), std::move(parent_dim),
                         std::move(parent_id), std::move(cells));
}

/// The prism as k^2 triangles of its base (as the triangle's template) in k
/// layers; the lattice points ordered by parent entity, vertices first.
IsoRefineTemplate make_prism_iso_p1_storage(int order)
{
    const double h = 1.0 / static_cast<double>(order);
    struct Point
    {
        int dim, id, i, j, l;
    };
    std::vector<Point> points;
    for (int l = 0; l <= order; ++l)
        for (int j = 0; j <= order; ++j)
            for (int i = 0; i <= order - j; ++i)
            {
                const double x[3] = {i * h, j * h, l * h};
                const auto [dim, id] = parent_entity(cell::type::prism, x);
                points.push_back({dim, id, i, j, l});
            }
    std::stable_sort(points.begin(), points.end(),
                     [](const Point& a, const Point& b) { return a.dim != b.dim ? a.dim < b.dim : a.id < b.id; });
    std::vector<double> coords;
    std::unordered_map<int, int> map_ijk;
    for (std::size_t v = 0; v < points.size(); ++v)
    {
        const Point& p = points[v];
        coords.insert(coords.end(), {p.i * h, p.j * h, p.l * h});
        map_ijk.emplace(key_ijk(p.i, p.j, p.l), static_cast<int>(v));
    }
    auto vertex = [&](int i, int j, int l) { return map_ijk.at(key_ijk(i, j, l)); };
    std::vector<cell::type> types;
    std::vector<int> cells;
    for (int l = 0; l < order; ++l)
        for (int i = 0; i < order; ++i)
            for (int j = 0; j < order - i; ++j)
            {
                const std::array<std::array<int, 2>, 3> up = {{{i, j}, {i + 1, j}, {i, j + 1}}};
                const std::array<std::array<int, 2>, 3> down = {{{i + 1, j}, {i + 1, j + 1}, {i, j + 1}}};
                for (int t = 0; t < (i + j <= order - 2 ? 2 : 1); ++t)
                {
                    const auto& tri = t == 0 ? up : down;
                    for (const int layer : {l, l + 1})
                        for (const auto& [a, b] : tri)
                            cells.push_back(vertex(a, b, layer));
                    types.push_back(cell::type::prism);
                }
            }
    return make_mixed_template(cell::type::prism, std::move(coords), std::move(types), std::move(cells));
}

/// The pyramid of order 2 on its 14 nodes (vertices, the midpoints of its
/// edges, the centre of its base): the pyramid on the square through the
/// midpoints of the slanted edges, the four pyramids on the quarters of the
/// base under that square's corners, the pyramid from that square down to the
/// base's centre, and four tetrahedra between them.
IsoRefineTemplate make_pyramid_iso_p1_storage()
{
    std::vector<double> coords = {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1};
    const std::vector<double> ref = cell::reference_vertices<double>(cell::type::pyramid);
    for (const std::array<int, 2>& e : cell::edges(cell::type::pyramid))
        for (int i = 0; i < 3; ++i)
            coords.push_back(0.5 * (ref[static_cast<std::size_t>(e[0] * 3 + i)]
                                    + ref[static_cast<std::size_t>(e[1] * 3 + i)]));
    coords.insert(coords.end(), {0.5, 0.5, 0});
    auto at = [&](double x, double y, double z)
    {
        for (std::size_t v = 0; v < coords.size() / 3; ++v)
            if (coords[3 * v] == x && coords[3 * v + 1] == y && coords[3 * v + 2] == z)
                return static_cast<int>(v);
        throw std::logic_error("iso_refine: no pyramid node there");
    };
    const int v0 = 0, v1 = 1, v2 = 2, v3 = 3, apex = 4;
    const int e01 = at(0.5, 0, 0), e02 = at(0, 0.5, 0), e13 = at(1, 0.5, 0), e23 = at(0.5, 1, 0);
    const int m0 = at(0, 0, 0.5), m1 = at(0.5, 0, 0.5), m2 = at(0, 0.5, 0.5), m3 = at(0.5, 0.5, 0.5);
    const int c = at(0.5, 0.5, 0);
    // pyramids: base in Basix order (v3 = v1 + v2 - v0), then the apex
    std::vector<int> cells = {m0, m1, m2, m3, apex, v0, e01, e02, c, m0, e01, v1, c, e13, m1,
                              e02, c, v2, e23, m2, c, e13, e23, v3, m3, m0, m1, m2, m3, c,
                              e01, c, m0, m1, e02, c, m0, m2, e13, c, m1, m3, e23, c, m2, m3};
    std::vector<cell::type> types(6, cell::type::pyramid);
    types.insert(types.end(), 4, cell::type::tetrahedron);
    return make_mixed_template(cell::type::pyramid, std::move(coords), std::move(types), std::move(cells));
}

} // namespace

const IsoRefineTemplate& p1_template(cell::type cell_type)
{
    switch (cell_type)
    {
    case cell::type::interval:
    {
        static const IsoRefineTemplate tpl = make_p1_storage(cell::type::interval);
        return tpl;
    }
    case cell::type::triangle:
    {
        static const IsoRefineTemplate tpl = make_p1_storage(cell::type::triangle);
        return tpl;
    }
    case cell::type::quadrilateral:
    {
        static const IsoRefineTemplate tpl = make_p1_storage(cell::type::quadrilateral);
        return tpl;
    }
    case cell::type::tetrahedron:
    {
        static const IsoRefineTemplate tpl = make_p1_storage(cell::type::tetrahedron);
        return tpl;
    }
    case cell::type::hexahedron:
    {
        static const IsoRefineTemplate tpl = make_p1_storage(cell::type::hexahedron);
        return tpl;
    }
    case cell::type::prism:
    {
        static const IsoRefineTemplate tpl = make_p1_storage(cell::type::prism);
        return tpl;
    }
    case cell::type::pyramid:
    {
        static const IsoRefineTemplate tpl = make_p1_storage(cell::type::pyramid);
        return tpl;
    }
    default:
        throw std::invalid_argument(
            "p1_template: unsupported cell type "
            + cell::cell_type_to_str(cell_type));
    }
}

const IsoRefineTemplate& iso_p1_template(cell::type cell_type, int order)
{
    if (order == 1)
        return p1_template(cell_type);
    if (order < 2 || order > 4)
        throw std::invalid_argument("iso_p1_template: order must be 1, 2, 3, or 4");

    switch (cell_type)
    {
    case cell::type::interval:
    {
        static const IsoRefineTemplate p2 = make_interval_iso_p1_storage(2);
        static const IsoRefineTemplate p3 = make_interval_iso_p1_storage(3);
        static const IsoRefineTemplate p4 = make_interval_iso_p1_storage(4);
        return order == 2 ? p2 : order == 3 ? p3 : p4;
    }
    case cell::type::triangle:
    {
        static const IsoRefineTemplate p2 = make_triangle_iso_p1_storage(2);
        static const IsoRefineTemplate p3 = make_triangle_iso_p1_storage(3);
        static const IsoRefineTemplate p4 = make_triangle_iso_p1_storage(4);
        return order == 2 ? p2 : order == 3 ? p3 : p4;
    }
    case cell::type::quadrilateral:
    {
        static const IsoRefineTemplate p2 = make_quadrilateral_iso_p1_storage(2);
        static const IsoRefineTemplate p3 = make_quadrilateral_iso_p1_storage(3);
        static const IsoRefineTemplate p4 = make_quadrilateral_iso_p1_storage(4);
        return order == 2 ? p2 : order == 3 ? p3 : p4;
    }
    case cell::type::tetrahedron:
    {
        static const IsoRefineTemplate p2 = make_tetrahedron_iso_p1_storage(2);
        static const IsoRefineTemplate p3 = make_tetrahedron_iso_p1_storage(3);
        static const IsoRefineTemplate p4 = make_tetrahedron_iso_p1_storage(4);
        return order == 2 ? p2 : order == 3 ? p3 : p4;
    }
    case cell::type::hexahedron:
    {
        static const IsoRefineTemplate p2 = make_hexahedron_iso_p1_storage(2);
        static const IsoRefineTemplate p3 = make_hexahedron_iso_p1_storage(3);
        static const IsoRefineTemplate p4 = make_hexahedron_iso_p1_storage(4);
        return order == 2 ? p2 : order == 3 ? p3 : p4;
    }
    case cell::type::prism:
    {
        static const IsoRefineTemplate p2 = make_prism_iso_p1_storage(2);
        static const IsoRefineTemplate p3 = make_prism_iso_p1_storage(3);
        static const IsoRefineTemplate p4 = make_prism_iso_p1_storage(4);
        return order == 2 ? p2 : order == 3 ? p3 : p4;
    }
    case cell::type::pyramid:
    {
        if (order > 2)
            throw std::invalid_argument("iso_p1_template: pyramids take orders 1 and 2");
        static const IsoRefineTemplate p2 = make_pyramid_iso_p1_storage();
        return p2;
    }
    default:
        throw std::invalid_argument(
            "iso_p1_template: unsupported cell type "
            + cell::cell_type_to_str(cell_type)
            + "; supported types are interval, triangle, tetrahedron, "
              "quadrilateral, hexahedron, prism and pyramid");
    }
}

std::span<const double> iso_p1_ref_coords(cell::type cell_type, int order)
{
    const auto& tpl = iso_p1_template(cell_type, order);
    return std::span<const double>(tpl.ref_vertex_coords.data(),
                                   tpl.ref_vertex_coords.size());
}

} // namespace cutcells
