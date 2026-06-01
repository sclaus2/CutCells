// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "iso_refine.h"

#include "cell_topology.h"
#include "reference_cell.h"

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
    if (tpl.child_cell_type != cell::type::tetrahedron || tpl.tdim != 3)
        return;

    for (int c = 0; c < tpl.n_cells; ++c)
    {
        const auto offset = static_cast<std::size_t>(4 * c);
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
    for (int t = 1; t < order; ++t) add_pt(order, t, 0, 1, 2);
    for (int t = 1; t < order; ++t) add_pt(t, order, 0, 1, 3);
    for (int t = 1; t < order; ++t) add_pt(0, 0, t, 1, 4);
    for (int t = 1; t < order; ++t) add_pt(order, 0, t, 1, 5);
    for (int t = 1; t < order; ++t) add_pt(0, order, t, 1, 6);
    for (int t = 1; t < order; ++t) add_pt(order, order, t, 1, 7);
    for (int t = 1; t < order; ++t) add_pt(t, 0, order, 1, 8);
    for (int t = 1; t < order; ++t) add_pt(0, t, order, 1, 9);
    for (int t = 1; t < order; ++t) add_pt(order, t, order, 1, 10);
    for (int t = 1; t < order; ++t) add_pt(t, order, order, 1, 11);

    for (int j = 1; j < order; ++j) for (int i = 1; i < order; ++i) add_pt(i, j, 0, 2, 4);
    for (int k = 1; k < order; ++k) for (int i = 1; i < order; ++i) add_pt(i, 0, k, 2, 2);
    for (int k = 1; k < order; ++k) for (int j = 1; j < order; ++j) add_pt(0, j, k, 2, 0);
    for (int k = 1; k < order; ++k) for (int j = 1; j < order; ++j) add_pt(order, j, k, 2, 1);
    for (int k = 1; k < order; ++k) for (int i = 1; i < order; ++i) add_pt(i, order, k, 2, 3);
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

void clear_entity_host_provenance(auto& ac, int dim)
{
    ac.entity_host_cell_id[dim].clear();
    ac.entity_host_cell_type[dim].clear();
    ac.entity_host_face_id[dim].clear();
    ac.entity_source_level_set[dim].clear();
    ac.entity_host_cell_vertices[dim].offsets.clear();
    ac.entity_host_cell_vertices[dim].indices.clear();
    ac.entity_host_cell_vertices[dim].offsets.push_back(std::int32_t(0));
}

template <std::floating_point T>
void append_entity_host_provenance(AdaptCell<T>& ac,
                                   int dim,
                                   int host_cell_id,
                                   cell::type host_cell_type,
                                   int host_face_id,
                                   int source_level_set,
                                   std::span<const std::int32_t> host_vertices)
{
    ac.entity_host_cell_id[dim].push_back(static_cast<std::int32_t>(host_cell_id));
    ac.entity_host_cell_type[dim].push_back(host_cell_type);
    ac.entity_host_face_id[dim].push_back(static_cast<std::int32_t>(host_face_id));
    ac.entity_source_level_set[dim].push_back(static_cast<std::int32_t>(source_level_set));
    for (const auto v : host_vertices)
        ac.entity_host_cell_vertices[dim].indices.push_back(v);
    ac.entity_host_cell_vertices[dim].offsets.push_back(
        static_cast<std::int32_t>(ac.entity_host_cell_vertices[dim].indices.size()));
}

template <std::floating_point T>
T edge_param_from_parent_reference(const AdaptCell<T>& ac,
                                   int parent_edge_id,
                                   std::span<const T> xi)
{
    const auto edges = cell::edges(ac.parent_cell_type);
    const auto edge = edges[static_cast<std::size_t>(parent_edge_id)];
    const auto ref_vertices = cell::reference_vertices<T>(ac.parent_cell_type);
    const T* a = ref_vertices.data() + static_cast<std::size_t>(edge[0] * ac.tdim);
    const T* b = ref_vertices.data() + static_cast<std::size_t>(edge[1] * ac.tdim);
    T ab2 = T(0);
    T ap_ab = T(0);
    for (int d = 0; d < ac.tdim; ++d)
    {
        const T ab = b[d] - a[d];
        ab2 += ab * ab;
        ap_ab += (xi[static_cast<std::size_t>(d)] - a[d]) * ab;
    }
    return ab2 > T(0) ? std::clamp(ap_ab / ab2, T(0), T(1)) : T(0);
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
    default:
        throw std::invalid_argument(
            "iso_p1_template: unsupported cell type "
            + cell::cell_type_to_str(cell_type)
            + "; supported v1 types are interval, triangle, tetrahedron, "
              "quadrilateral, and hexahedron");
    }
}

std::span<const double> iso_p1_ref_coords(cell::type cell_type, int order)
{
    const auto& tpl = iso_p1_template(cell_type, order);
    return std::span<const double>(tpl.ref_vertex_coords.data(),
                                   tpl.ref_vertex_coords.size());
}

template <std::floating_point T>
void apply_iso_refine(AdaptCell<T>& ac, const IsoRefineTemplate& tpl)
{
    if (tpl.parent_cell_type != ac.parent_cell_type)
    {
        throw std::invalid_argument(
            "apply_iso_refine: template cell type does not match AdaptCell parent");
    }
    if (tpl.tdim != ac.tdim)
        throw std::invalid_argument("apply_iso_refine: template tdim mismatch");

    const int tdim = ac.tdim;
    const int n_vertices = tpl.n_vertices;
    const int n_cells = tpl.n_cells;
    const int parent_corners = cell::get_num_vertices(ac.parent_cell_type);
    const std::vector<std::int32_t> old_parent_ids = ac.vertex_parent_id;

    ac.vertex_coords.clear();
    ac.vertex_coords.reserve(static_cast<std::size_t>(n_vertices * tdim));
    ac.vertex_parent_dim.clear();
    ac.vertex_parent_id.clear();
    ac.vertex_parent_param.clear();
    ac.vertex_parent_param_offset.clear();
    ac.vertex_parent_param_offset.push_back(0);
    ac.vertex_source_edge_id.assign(static_cast<std::size_t>(n_vertices), std::int32_t(-1));
    ac.zero_mask_per_vertex.assign(static_cast<std::size_t>(n_vertices), std::uint64_t(0));
    ac.negative_mask_per_vertex.assign(static_cast<std::size_t>(n_vertices), std::uint64_t(0));

    for (int v = 0; v < n_vertices; ++v)
    {
        std::vector<T> xi(static_cast<std::size_t>(tdim), T(0));
        for (int d = 0; d < tdim; ++d)
        {
            xi[static_cast<std::size_t>(d)] = static_cast<T>(
                tpl.ref_vertex_coords[static_cast<std::size_t>(v * tdim + d)]);
            ac.vertex_coords.push_back(xi[static_cast<std::size_t>(d)]);
        }

        const int dim = tpl.vertex_parent_dim[static_cast<std::size_t>(v)];
        int parent_id = tpl.vertex_parent_id[static_cast<std::size_t>(v)];
        std::vector<T> parent_param;

        if (dim == 0)
        {
            if (parent_id >= 0 && parent_id < parent_corners
                && parent_id < static_cast<int>(old_parent_ids.size()))
            {
                parent_id = old_parent_ids[static_cast<std::size_t>(parent_id)];
            }
        }
        else if (dim == 1)
        {
            parent_param.push_back(edge_param_from_parent_reference<T>(
                ac, parent_id, std::span<const T>(xi.data(), xi.size())));
        }
        else if (dim == tdim)
        {
            parent_id = ac.parent_cell_id;
            parent_param = xi;
        }
        else
        {
            parent_param.assign(xi.begin(), xi.begin() + dim);
        }

        ac.vertex_parent_dim.push_back(static_cast<std::int8_t>(dim));
        ac.vertex_parent_id.push_back(static_cast<std::int32_t>(parent_id));
        ac.vertex_parent_param.insert(ac.vertex_parent_param.end(),
                                      parent_param.begin(), parent_param.end());
        ac.vertex_parent_param_offset.push_back(
            static_cast<std::int32_t>(ac.vertex_parent_param.size()));
    }

    for (auto& types : ac.entity_types)
        types.clear();
    for (auto& adj : ac.entity_to_vertex)
    {
        adj.offsets.clear();
        adj.indices.clear();
        adj.offsets.push_back(0);
    }

    ac.entity_types[tdim].assign(static_cast<std::size_t>(n_cells),
                                 tpl.child_cell_type);
    ac.entity_to_vertex[tdim].indices.reserve(tpl.cell_connectivity.size());
    ac.entity_to_vertex[tdim].offsets.clear();
    ac.entity_to_vertex[tdim].offsets.push_back(0);
    for (int c = 0; c < n_cells; ++c)
    {
        for (int j = 0; j < tpl.vertices_per_cell; ++j)
        {
            ac.entity_to_vertex[tdim].indices.push_back(static_cast<std::int32_t>(
                tpl.cell_connectivity[static_cast<std::size_t>(
                    c * tpl.vertices_per_cell + j)]));
        }
        ac.entity_to_vertex[tdim].offsets.push_back(static_cast<std::int32_t>(
            ac.entity_to_vertex[tdim].indices.size()));
    }

    ac.cell_source_cell_id.assign(static_cast<std::size_t>(n_cells), std::int32_t(0));
    ac.cell_refinement_generation.assign(static_cast<std::size_t>(n_cells), std::int32_t(1));
    ac.cell_refinement_reason.assign(static_cast<std::size_t>(n_cells),
                                     CellRefinementReason::iso_refine);
    ac.cell_host_parent_cell_id.assign(static_cast<std::size_t>(n_cells),
                                       std::int32_t(ac.parent_cell_id));

    for (int dim = 0; dim < AdaptCell<T>::max_dim; ++dim)
        clear_entity_host_provenance(ac, dim);

    for (int c = 0; c < n_cells; ++c)
    {
        const auto verts = ac.entity_to_vertex[tdim][static_cast<std::int32_t>(c)];
        append_entity_host_provenance<T>(
            ac, tdim, c, tpl.child_cell_type, -1, -1, verts);
    }

    for (auto& row : ac.connectivity)
    {
        for (auto& conn : row)
        {
            conn.offsets.clear();
            conn.indices.clear();
        }
    }
    ac.has_connectivity = {};
    ac.zero_entity_dim.clear();
    ac.zero_entity_id.clear();
    ac.zero_entity_zero_mask.clear();
    ac.zero_entity_is_owned.clear();
    ac.zero_entity_parent_dim.clear();
    ac.zero_entity_parent_id.clear();
    ac.zero_entity_host_cell_id.clear();
    ac.zero_entity_host_cell_type.clear();
    ac.zero_entity_host_face_id.clear();
    ac.zero_entity_source_level_set.clear();
    ac.zero_entity_host_cell_vertices.offsets.clear();
    ac.zero_entity_host_cell_vertices.indices.clear();
    ac.zero_entity_host_cell_vertices.offsets.push_back(0);
    ++ac.zero_entity_version;

    ac.edge_root_tag.clear();
    ac.edge_root_tag_num_level_sets = 0;
    ac.edge_green_split_param.clear();
    ac.edge_green_split_has_value.clear();
    ac.edge_one_root_param.clear();
    ac.edge_one_root_vertex_id.clear();
    ac.edge_one_root_has_value.clear();
    ac.cell_cert_tag.clear();
    ac.cell_cert_tag_num_level_sets = 0;
    ac.face_cert_tag.clear();
    ac.face_cert_tag_num_level_sets = 0;

    build_edges(ac);
    if (tdim == 3)
        build_faces(ac);
    recompute_active_level_set_masks(ac, /*num_level_sets=*/0);
}

template void apply_iso_refine(AdaptCell<double>&, const IsoRefineTemplate&);
template void apply_iso_refine(AdaptCell<float>&, const IsoRefineTemplate&);

} // namespace cutcells
