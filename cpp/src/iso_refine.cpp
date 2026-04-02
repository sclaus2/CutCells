// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "iso_refine.h"
#include "generated/iso_refine_vtk_blocks.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace cutcells
{
namespace
{

// Basix triangle edge order is: (1,2), (0,2), (0,1).

inline int key_ij(int i, int j)
{
    return i * 100 + j;
}

inline int key_ijk(int i, int j, int k)
{
    return i * 10000 + j * 100 + k;
}

int vtk_interval_point_index(int i, int order)
{
    if (i == 0)
        return 0;
    if (i == order)
        return 1;
    return 1 + i;
}

int vtk_triangle_point_index(const std::array<int, 3>& bindex, int order)
{
    int index = 0;
    int max = order;
    int min = 0;
    const int bmin = std::min({bindex[0], bindex[1], bindex[2]});

    while (bmin > min)
    {
        index += 3 * order;
        max -= 2;
        min += 1;
        order -= 3;
    }

    for (int dim = 0; dim < 3; ++dim)
    {
        if (bindex[(dim + 2) % 3] == max)
            return index;
        ++index;
    }

    for (int dim = 0; dim < 3; ++dim)
    {
        if (bindex[(dim + 1) % 3] == min)
            return index + bindex[dim] - (min + 1);
        index += max - (min + 1);
    }

    return index;
}

int vtk_quadrilateral_point_index(int i, int j, int order)
{
    const bool ibdy = (i == 0 || i == order);
    const bool jbdy = (j == 0 || j == order);
    const int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0);

    if (nbdy == 2)
        return (i ? (j ? 2 : 1) : (j ? 3 : 0));

    int offset = 4;
    if (nbdy == 1)
    {
        if (!ibdy)
            return (i - 1) + (j ? order - 1 + order - 1 : 0) + offset;
        if (!jbdy)
            return (j - 1) + (i ? order - 1 : 2 * (order - 1) + order - 1) + offset;
    }

    offset += 4 * (order - 1);
    return offset + (i - 1) + (order - 1) * (j - 1);
}

std::vector<int> build_vtk_lagrange_to_basix_permutation(cell::type ct, int order)
{
    if (order < 1 || order > 4)
        throw std::invalid_argument("vtk_lagrange_to_basix_permutation supports orders 1..4");

    const std::span<const double> ref_coords = (order == 1)
        ? p1_ref_coords(ct)
        : iso_p1_ref_coords(ct, order);

    int tdim = 0;
    switch (ct)
    {
        case cell::type::interval:
            tdim = 1;
            break;
        case cell::type::triangle:
        case cell::type::quadrilateral:
            tdim = 2;
            break;
        default:
            throw std::invalid_argument(
                "vtk_lagrange_to_basix_permutation supports only interval/triangle/quadrilateral");
    }

    const int n_points = static_cast<int>(ref_coords.size()) / tdim;
    std::vector<int> perm(static_cast<std::size_t>(n_points), -1);

    for (int basix_id = 0; basix_id < n_points; ++basix_id)
    {
        const double* x = &ref_coords[static_cast<std::size_t>(basix_id * tdim)];
        int vtk_id = -1;

        switch (ct)
        {
            case cell::type::interval:
            {
                const int i = static_cast<int>(std::lround(x[0] * static_cast<double>(order)));
                vtk_id = vtk_interval_point_index(i, order);
                break;
            }
            case cell::type::triangle:
            {
                const int b0 = static_cast<int>(std::lround(x[0] * static_cast<double>(order)));
                const int b1 = static_cast<int>(std::lround(x[1] * static_cast<double>(order)));
                const int b2 = order - b0 - b1;
                vtk_id = vtk_triangle_point_index({b0, b1, b2}, order);
                break;
            }
            case cell::type::quadrilateral:
            {
                const int i = static_cast<int>(std::lround(x[0] * static_cast<double>(order)));
                const int j = static_cast<int>(std::lround(x[1] * static_cast<double>(order)));
                vtk_id = vtk_quadrilateral_point_index(i, j, order);
                break;
            }
            default:
                break;
        }

        if (vtk_id < 0 || vtk_id >= n_points)
            throw std::runtime_error("vtk_lagrange_to_basix_permutation produced invalid index");
        perm[static_cast<std::size_t>(vtk_id)] = basix_id;
    }

    for (int vtk_id = 0; vtk_id < n_points; ++vtk_id)
    {
        if (perm[static_cast<std::size_t>(vtk_id)] < 0)
            throw std::runtime_error("vtk_lagrange_to_basix_permutation is incomplete");
    }

    return perm;
}

std::vector<std::array<int, 3>> triangle_lattice_cells(int order)
{
    std::vector<std::array<int, 3>> tris;
    tris.reserve(order * order);
    for (int i = 0; i < order; ++i)
    {
        for (int j = 0; j < order - i; ++j)
        {
            tris.push_back({key_ij(i, j), key_ij(i + 1, j), key_ij(i, j + 1)});
            if (i + j <= order - 2)
                tris.push_back({key_ij(i + 1, j), key_ij(i + 1, j + 1), key_ij(i, j + 1)});
        }
    }
    return tris;
}

RefinementTemplate make_triangle_iso_p1_storage(
    int order,
    const std::vector<std::array<double, 2>>& basix_points)
{
    if (order < 2)
        throw std::invalid_argument("Triangle iso refine order must be >= 2");

    RefinementTemplate out;
    const int n_vertices = static_cast<int>(basix_points.size());
    out.ref_vertex_coords.reserve(n_vertices * 2);
    out.vertex_parent_dim.reserve(n_vertices);
    out.vertex_parent_id.reserve(n_vertices);

    // Build a lattice map (i,j) -> Basix vertex id where x=i/p, y=j/p.
    std::unordered_map<int, int> lattice_to_vertex;
    lattice_to_vertex.reserve(n_vertices);

    for (int vid = 0; vid < n_vertices; ++vid)
    {
        const double x = basix_points[vid][0];
        const double y = basix_points[vid][1];
        out.ref_vertex_coords.push_back(x);
        out.ref_vertex_coords.push_back(y);

        const int i = static_cast<int>(std::lround(x * static_cast<double>(order)));
        const int j = static_cast<int>(std::lround(y * static_cast<double>(order)));
        lattice_to_vertex.emplace(key_ij(i, j), vid);
    }

    // Parent entity tags.
    // Basix ordering:
    //   1. vertices
    //   2. edge dofs in edge order (1,2), (0,2), (0,1)
    //   3. interior dofs
    for (int vid = 0; vid < 3; ++vid)
    {
        out.vertex_parent_dim.push_back(0);
        out.vertex_parent_id.push_back(vid);
    }

    int cursor = 3;
    for (int e = 0; e < 3; ++e)
    {
        for (int t = 1; t < order; ++t)
        {
            out.vertex_parent_dim.push_back(1);
            out.vertex_parent_id.push_back(e);
            ++cursor;
        }
    }

    // Interior points.
    while (cursor < n_vertices)
    {
        out.vertex_parent_dim.push_back(2);
        out.vertex_parent_id.push_back(0);
        ++cursor;
    }

    // Uniform Pk-iso-P1 triangulation (p^2 triangles).
    // Lattice coordinates: (i,j), 0<=i, 0<=j, i+j<=p.
    // Cell triangles:
    //   (i,j) -> (i+1,j) -> (i,j+1)
    //   (i+1,j) -> (i+1,j+1) -> (i,j+1) when i+j <= p-2
    out.cell_connectivity.reserve(order * order * 3);
    auto vid = [&](int i, int j) -> int {
        const auto it = lattice_to_vertex.find(key_ij(i, j));
        if (it == lattice_to_vertex.end())
            throw std::runtime_error("Missing lattice vertex for triangle iso refinement");
        return it->second;
    };

    for (int i = 0; i < order; ++i)
    {
        for (int j = 0; j < order - i; ++j)
        {
            out.cell_connectivity.push_back(vid(i, j));
            out.cell_connectivity.push_back(vid(i + 1, j));
            out.cell_connectivity.push_back(vid(i, j + 1));

            if (i + j <= order - 2)
            {
                out.cell_connectivity.push_back(vid(i + 1, j));
                out.cell_connectivity.push_back(vid(i + 1, j + 1));
                out.cell_connectivity.push_back(vid(i, j + 1));
            }
        }
    }

    out.n_vertices = n_vertices;
    out.n_cells = static_cast<int>(out.cell_connectivity.size()) / 3;
    out.tdim = 2;
    out.vertices_per_cell = 3;
    out.bg_cell_type = cell::type::triangle;
    out.child_cell_type = cell::type::triangle;

    return out;
}

const std::vector<std::array<double, 2>>& p2_points()
{
    // Basix equispaced P2 triangle points, order verified against Basix 0.11.
    static const std::vector<std::array<double, 2>> pts = {
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.5, 0.5}, {0.0, 0.5}, {0.5, 0.0}
    };
    return pts;
}

const std::vector<std::array<double, 2>>& p3_points()
{
    // Basix equispaced P3 triangle points, order verified against Basix 0.11.
    static const std::vector<std::array<double, 2>> pts = {
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {2.0 / 3.0, 1.0 / 3.0}, {1.0 / 3.0, 2.0 / 3.0},
        {0.0, 1.0 / 3.0}, {0.0, 2.0 / 3.0},
        {1.0 / 3.0, 0.0}, {2.0 / 3.0, 0.0},
        {1.0 / 3.0, 1.0 / 3.0}
    };
    return pts;
}

const std::vector<std::array<double, 2>>& p4_points()
{
    // Basix equispaced P4 triangle points, order verified against Basix 0.11.
    static const std::vector<std::array<double, 2>> pts = {
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.75, 0.25}, {0.5, 0.5}, {0.25, 0.75},
        {0.0, 0.25}, {0.0, 0.5}, {0.0, 0.75},
        {0.25, 0.0}, {0.5, 0.0}, {0.75, 0.0},
        {0.25, 0.25}, {0.5, 0.25}, {0.25, 0.5}
    };
    return pts;
}

RefinementTemplate make_template(
    cell::type bg_cell_type,
    cell::type child_cell_type,
    int tdim,
    int vertices_per_cell,
    std::vector<double> ref_coords,
    std::vector<int> parent_dim,
    std::vector<int> parent_id,
    std::vector<int> cells)
{
    RefinementTemplate tpl{
        .n_vertices = static_cast<int>(parent_dim.size()),
        .n_cells = static_cast<int>(cells.size()) / vertices_per_cell,
        .tdim = tdim,
        .vertices_per_cell = vertices_per_cell,
        .bg_cell_type = bg_cell_type,
        .child_cell_type = child_cell_type,
        .ref_vertex_coords = std::move(ref_coords),
        .vertex_parent_dim = std::move(parent_dim),
        .vertex_parent_id = std::move(parent_id),
        .cell_connectivity = std::move(cells)
    };
    return tpl;
}

RefinementTemplate make_interval_iso_p1_storage(int order)
{
    std::vector<double> coords;
    std::vector<int> parent_dim;
    std::vector<int> parent_id;
    std::vector<int> cells;

    coords.reserve(order + 1);
    parent_dim.reserve(order + 1);
    parent_id.reserve(order + 1);
    cells.reserve(order * 2);

    // Vertices
    coords.push_back(0.0);
    coords.push_back(1.0);
    parent_dim.push_back(0);
    parent_dim.push_back(0);
    parent_id.push_back(0);
    parent_id.push_back(1);

    // Edge points
    for (int i = 1; i < order; ++i)
    {
        coords.push_back(static_cast<double>(i) / order);
        parent_dim.push_back(1);
        parent_id.push_back(0);
    }

    auto v = [&](int i) -> int
    {
        if (i == 0) return 0;
        if (i == order) return 1;
        return 1 + i;
    };

    for (int i = 0; i < order; ++i)
    {
        cells.push_back(v(i));
        cells.push_back(v(i + 1));
    }

    return make_template(
        cell::type::interval,
        cell::type::interval,
        1,
        2,
        std::move(coords),
        std::move(parent_dim),
        std::move(parent_id),
        std::move(cells));
}

RefinementTemplate make_quadrilateral_iso_p1_storage(int order)
{
    const double h = 1.0 / static_cast<double>(order);
    std::vector<double> coords;
    std::vector<int> parent_dim;
    std::vector<int> parent_id;
    std::vector<int> cells;
    std::unordered_map<int, int> map_ij;

    auto add_pt = [&](int i, int j, int pd, int pid)
    {
        const int idx = static_cast<int>(coords.size() / 2);
        coords.push_back(i * h);
        coords.push_back(j * h);
        parent_dim.push_back(pd);
        parent_id.push_back(pid);
        map_ij.emplace(key_ij(i, j), idx);
    };

    add_pt(0, 0, 0, 0);
    add_pt(order, 0, 0, 1);
    add_pt(0, order, 0, 2);
    add_pt(order, order, 0, 3);

    // Basix edge order in cell_topology.h
    for (int t = 1; t < order; ++t) add_pt(t, 0, 1, 0);          // (0,1)
    for (int t = 1; t < order; ++t) add_pt(0, t, 1, 1);          // (0,2)
    for (int t = 1; t < order; ++t) add_pt(order, t, 1, 2);      // (1,3)
    for (int t = 1; t < order; ++t) add_pt(t, order, 1, 3);      // (2,3)

    for (int j = 1; j < order; ++j)
        for (int i = 1; i < order; ++i)
            add_pt(i, j, 2, 0);

    auto vid = [&](int i, int j) -> int
    {
        return map_ij.at(key_ij(i, j));
    };

    for (int j = 0; j < order; ++j)
    {
        for (int i = 0; i < order; ++i)
        {
            cells.push_back(vid(i, j));
            cells.push_back(vid(i + 1, j));
            cells.push_back(vid(i, j + 1));
            cells.push_back(vid(i + 1, j + 1));
        }
    }

    return make_template(
        cell::type::quadrilateral,
        cell::type::quadrilateral,
        2,
        4,
        std::move(coords),
        std::move(parent_dim),
        std::move(parent_id),
        std::move(cells));
}

RefinementTemplate make_hexahedron_iso_p1_storage(int order)
{
    const double h = 1.0 / static_cast<double>(order);
    std::vector<double> coords;
    std::vector<int> parent_dim;
    std::vector<int> parent_id;
    std::vector<int> cells;
    std::unordered_map<int, int> map_ijk;

    auto add_pt = [&](int i, int j, int k, int pd, int pid)
    {
        const int idx = static_cast<int>(coords.size() / 3);
        coords.push_back(i * h);
        coords.push_back(j * h);
        coords.push_back(k * h);
        parent_dim.push_back(pd);
        parent_id.push_back(pid);
        map_ijk.emplace(key_ijk(i, j, k), idx);
    };

    // Vertices in Basix order
    add_pt(0, 0, 0, 0, 0); add_pt(order, 0, 0, 0, 1);
    add_pt(0, order, 0, 0, 2); add_pt(order, order, 0, 0, 3);
    add_pt(0, 0, order, 0, 4); add_pt(order, 0, order, 0, 5);
    add_pt(0, order, order, 0, 6); add_pt(order, order, order, 0, 7);

    // Edge points in Basix edge order
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

    // Face points in Basix face order
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

    auto vid = [&](int i, int j, int k) -> int
    {
        return map_ijk.at(key_ijk(i, j, k));
    };

    for (int k = 0; k < order; ++k)
    {
        for (int j = 0; j < order; ++j)
        {
            for (int i = 0; i < order; ++i)
            {
                cells.push_back(vid(i, j, k));
                cells.push_back(vid(i + 1, j, k));
                cells.push_back(vid(i, j + 1, k));
                cells.push_back(vid(i + 1, j + 1, k));
                cells.push_back(vid(i, j, k + 1));
                cells.push_back(vid(i + 1, j, k + 1));
                cells.push_back(vid(i, j + 1, k + 1));
                cells.push_back(vid(i + 1, j + 1, k + 1));
            }
        }
    }

    return make_template(
        cell::type::hexahedron,
        cell::type::hexahedron,
        3,
        8,
        std::move(coords),
        std::move(parent_dim),
        std::move(parent_id),
        std::move(cells));
}

RefinementTemplate make_prism_iso_p1_storage(int order)
{
    const double h = 1.0 / static_cast<double>(order);
    std::vector<double> coords;
    std::vector<int> parent_dim;
    std::vector<int> parent_id;
    std::vector<int> cells;
    std::unordered_map<int, int> map_ijk;

    auto add_pt = [&](int i, int j, int k, int pd, int pid)
    {
        const int idx = static_cast<int>(coords.size() / 3);
        coords.push_back(i * h);
        coords.push_back(j * h);
        coords.push_back(k * h);
        parent_dim.push_back(pd);
        parent_id.push_back(pid);
        map_ijk.emplace(key_ijk(i, j, k), idx);
    };

    // Vertices (Basix order)
    add_pt(0, 0, 0, 0, 0);          // 0
    add_pt(order, 0, 0, 0, 1);      // 1
    add_pt(0, order, 0, 0, 2);      // 2
    add_pt(0, 0, order, 0, 3);      // 3
    add_pt(order, 0, order, 0, 4);  // 4
    add_pt(0, order, order, 0, 5);  // 5

    // Edges in Basix order from cell_topology
    for (int t = 1; t < order; ++t) add_pt(t, 0, 0, 1, 0);              // (0,1)
    for (int t = 1; t < order; ++t) add_pt(0, t, 0, 1, 1);              // (0,2)
    for (int t = 1; t < order; ++t) add_pt(0, 0, t, 1, 2);              // (0,3)
    for (int t = 1; t < order; ++t) add_pt(order - t, t, 0, 1, 3);      // (1,2)
    for (int t = 1; t < order; ++t) add_pt(order, 0, t, 1, 4);          // (1,4)
    for (int t = 1; t < order; ++t) add_pt(0, order, t, 1, 5);          // (2,5)
    for (int t = 1; t < order; ++t) add_pt(t, 0, order, 1, 6);          // (3,4)
    for (int t = 1; t < order; ++t) add_pt(0, t, order, 1, 7);          // (3,5)
    for (int t = 1; t < order; ++t) add_pt(order - t, t, order, 1, 8);  // (4,5)

    // Face interiors in Basix order
    // Face 0: triangle [0,1,2] at z=0
    for (int j = 1; j < order; ++j)
        for (int i = 1; i + j < order; ++i)
            add_pt(i, j, 0, 2, 0);

    // Face 1: quad [0,1,3,4] (y=0)
    for (int k = 1; k < order; ++k)
        for (int i = 1; i < order; ++i)
            add_pt(i, 0, k, 2, 1);

    // Face 2: quad [0,2,3,5] (x=0)
    for (int k = 1; k < order; ++k)
        for (int j = 1; j < order; ++j)
            add_pt(0, j, k, 2, 2);

    // Face 3: quad [1,2,4,5] (x+y=order)
    for (int k = 1; k < order; ++k)
        for (int t = 1; t < order; ++t)
            add_pt(order - t, t, k, 2, 3);

    // Face 4: triangle [3,4,5] at z=1
    for (int j = 1; j < order; ++j)
        for (int i = 1; i + j < order; ++i)
            add_pt(i, j, order, 2, 4);

    // Cell interior
    for (int k = 1; k < order; ++k)
        for (int j = 1; j < order; ++j)
            for (int i = 1; i + j < order; ++i)
                add_pt(i, j, k, 3, 0);

    auto vid = [&](int i, int j, int k) -> int
    {
        return map_ijk.at(key_ijk(i, j, k));
    };

    const auto tri_cells = triangle_lattice_cells(order);
    cells.reserve(static_cast<std::size_t>(order) * tri_cells.size() * 6);
    for (int k = 0; k < order; ++k)
    {
        for (const auto& tri : tri_cells)
        {
            const int a_i = tri[0] / 100, a_j = tri[0] % 100;
            const int b_i = tri[1] / 100, b_j = tri[1] % 100;
            const int c_i = tri[2] / 100, c_j = tri[2] % 100;

            cells.push_back(vid(a_i, a_j, k));
            cells.push_back(vid(b_i, b_j, k));
            cells.push_back(vid(c_i, c_j, k));
            cells.push_back(vid(a_i, a_j, k + 1));
            cells.push_back(vid(b_i, b_j, k + 1));
            cells.push_back(vid(c_i, c_j, k + 1));
        }
    }

    return make_template(
        cell::type::prism,
        cell::type::prism,
        3,
        6,
        std::move(coords),
        std::move(parent_dim),
        std::move(parent_id),
        std::move(cells));
}

RefinementTemplate make_pyramid_iso_p1_storage(int order)
{
    const double h = 1.0 / static_cast<double>(order);
    std::vector<double> coords;
    std::vector<int> parent_dim;
    std::vector<int> parent_id;
    std::vector<int> cells;
    std::unordered_map<int, int> map_ijk;

    auto add_pt = [&](int i, int j, int k, int pd, int pid)
    {
        const int idx = static_cast<int>(coords.size() / 3);
        coords.push_back(i * h);
        coords.push_back(j * h);
        coords.push_back(k * h);
        parent_dim.push_back(pd);
        parent_id.push_back(pid);
        map_ijk.emplace(key_ijk(i, j, k), idx);
    };

    // Vertices
    add_pt(0, 0, 0, 0, 0);
    add_pt(order, 0, 0, 0, 1);
    add_pt(0, order, 0, 0, 2);
    add_pt(order, order, 0, 0, 3);
    add_pt(0, 0, order, 0, 4);

    // Edges in Basix order
    for (int t = 1; t < order; ++t) add_pt(t, 0, 0, 1, 0);                  // (0,1)
    for (int t = 1; t < order; ++t) add_pt(0, t, 0, 1, 1);                  // (0,2)
    for (int t = 1; t < order; ++t) add_pt(0, 0, t, 1, 2);                  // (0,4)
    for (int t = 1; t < order; ++t) add_pt(order, t, 0, 1, 3);              // (1,3)
    for (int t = 1; t < order; ++t) add_pt(order - t, 0, t, 1, 4);          // (1,4)
    for (int t = 1; t < order; ++t) add_pt(t, order, 0, 1, 5);              // (2,3)
    for (int t = 1; t < order; ++t) add_pt(0, order - t, t, 1, 6);          // (2,4)
    for (int t = 1; t < order; ++t) add_pt(order - t, order - t, t, 1, 7);  // (3,4)

    // Face interiors in Basix order
    // Face 0: base quad [0,1,2,3], z=0
    for (int j = 1; j < order; ++j)
        for (int i = 1; i < order; ++i)
            add_pt(i, j, 0, 2, 0);

    // Face 1: tri [0,1,4], y=0
    for (int k = 1; k < order; ++k)
        for (int i = 1; i + k < order; ++i)
            add_pt(i, 0, k, 2, 1);

    // Face 2: tri [0,2,4], x=0
    for (int k = 1; k < order; ++k)
        for (int j = 1; j + k < order; ++j)
            add_pt(0, j, k, 2, 2);

    // Face 3: tri [1,3,4]
    for (int k = 1; k < order; ++k)
        for (int t = 1; t + k < order; ++t)
            add_pt(order - k, t, k, 2, 3);

    // Face 4: tri [2,3,4]
    for (int k = 1; k < order; ++k)
        for (int t = 1; t + k < order; ++t)
            add_pt(t, order - k, k, 2, 4);

    // Cell interior
    for (int k = 1; k < order; ++k)
    {
        const int n = order - k;
        for (int j = 1; j < n; ++j)
            for (int i = 1; i < n; ++i)
                add_pt(i, j, k, 3, 0);
    }

    auto vid = [&](int i, int j, int k) -> int
    {
        return map_ijk.at(key_ijk(i, j, k));
    };

    // Layered pyramid decomposition.
    for (int k = 0; k < order; ++k)
    {
        const int n = order - k;
        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                cells.push_back(vid(i, j, k));
                cells.push_back(vid(i + 1, j, k));
                cells.push_back(vid(i, j + 1, k));
                cells.push_back(vid(i + 1, j + 1, k));
                cells.push_back(vid(i, j, k + 1));
            }
        }
    }

    return make_template(
        cell::type::pyramid,
        cell::type::pyramid,
        3,
        5,
        std::move(coords),
        std::move(parent_dim),
        std::move(parent_id),
        std::move(cells));
}

RefinementTemplate topology_only(const RefinementTemplate& in)
{
    RefinementTemplate out = in;
    out.ref_vertex_coords.clear();
    return out;
}

} // namespace

const RefinementTemplate& p2_iso_p1_triangle_template()
{
    static const RefinementTemplate s = topology_only(iso_refine_generated::triangle_p2_vtk_tess_basix_tpl);
    return s;
}

const RefinementTemplate& p3_iso_p1_triangle_template()
{
    static const RefinementTemplate s = topology_only(iso_refine_generated::triangle_p3_vtk_tess_basix_tpl);
    return s;
}

const RefinementTemplate& p4_iso_p1_triangle_template()
{
    static const RefinementTemplate s = topology_only(iso_refine_generated::triangle_p4_vtk_tess_basix_tpl);
    return s;
}

const RefinementTemplate& triangle_iso_p1_template(int order)
{
    switch (order)
    {
        case 2:
            return p2_iso_p1_triangle_template();
        case 3:
            return p3_iso_p1_triangle_template();
        case 4:
            return p4_iso_p1_triangle_template();
        default:
            throw std::invalid_argument("triangle_iso_p1_template supports only orders 2, 3, 4");
    }
}

const RefinementTemplate& p2_iso_p1_tetrahedron_template()
{
    static const RefinementTemplate s = topology_only(iso_refine_generated::tetrahedron_p2_vtk_tess_basix_tpl);
    return s;
}

const RefinementTemplate& p3_iso_p1_tetrahedron_template()
{
    static const RefinementTemplate s = topology_only(iso_refine_generated::tetrahedron_p3_vtk_tess_basix_tpl);
    return s;
}

const RefinementTemplate& p4_iso_p1_tetrahedron_template()
{
    static const RefinementTemplate s = topology_only(iso_refine_generated::tetrahedron_p4_vtk_tess_basix_tpl);
    return s;
}

const RefinementTemplate& tetrahedron_iso_p1_template(int order)
{
    switch (order)
    {
        case 2: return p2_iso_p1_tetrahedron_template();
        case 3: return p3_iso_p1_tetrahedron_template();
        case 4: return p4_iso_p1_tetrahedron_template();
        default:
            throw std::invalid_argument("tetrahedron_iso_p1_template supports orders 2, 3, 4");
    }
}

const RefinementTemplate& p2_iso_p1_quadrilateral_template()
{
    static const RefinementTemplate s = topology_only(iso_refine_generated::quadrilateral_p2_vtk_tess_basix_tpl);
    return s;
}

const RefinementTemplate& p3_iso_p1_quadrilateral_template()
{
    static const RefinementTemplate s = topology_only(iso_refine_generated::quadrilateral_p3_vtk_tess_basix_tpl);
    return s;
}

const RefinementTemplate& p4_iso_p1_quadrilateral_template()
{
    static const RefinementTemplate s = topology_only(iso_refine_generated::quadrilateral_p4_vtk_tess_basix_tpl);
    return s;
}

const RefinementTemplate& quadrilateral_iso_p1_template(int order)
{
    switch (order)
    {
        case 2: return p2_iso_p1_quadrilateral_template();
        case 3: return p3_iso_p1_quadrilateral_template();
        case 4: return p4_iso_p1_quadrilateral_template();
        default:
            throw std::invalid_argument("quadrilateral_iso_p1_template supports orders 2, 3, 4");
    }
}

const RefinementTemplate& p2_iso_p1_hexahedron_template()
{
    static const RefinementTemplate s = topology_only(iso_refine_generated::hexahedron_p2_vtk_tess_basix_tpl);
    return s;
}

const RefinementTemplate& p3_iso_p1_hexahedron_template()
{
    static const RefinementTemplate s = topology_only(iso_refine_generated::hexahedron_p3_vtk_tess_basix_tpl);
    return s;
}

const RefinementTemplate& p4_iso_p1_hexahedron_template()
{
    static const RefinementTemplate s = topology_only(iso_refine_generated::hexahedron_p4_vtk_tess_basix_tpl);
    return s;
}

const RefinementTemplate& hexahedron_iso_p1_template(int order)
{
    switch (order)
    {
        case 2: return p2_iso_p1_hexahedron_template();
        case 3: return p3_iso_p1_hexahedron_template();
        case 4: return p4_iso_p1_hexahedron_template();
        default:
            throw std::invalid_argument("hexahedron_iso_p1_template supports orders 2, 3, 4");
    }
}

const RefinementTemplate& p2_iso_p1_interval_template()
{
    static const RefinementTemplate full = make_interval_iso_p1_storage(2);
    static const RefinementTemplate s = topology_only(full);
    return s;
}

const RefinementTemplate& p3_iso_p1_interval_template()
{
    static const RefinementTemplate full = make_interval_iso_p1_storage(3);
    static const RefinementTemplate s = topology_only(full);
    return s;
}

const RefinementTemplate& p4_iso_p1_interval_template()
{
    static const RefinementTemplate full = make_interval_iso_p1_storage(4);
    static const RefinementTemplate s = topology_only(full);
    return s;
}

const RefinementTemplate& interval_iso_p1_template(int order)
{
    switch (order)
    {
        case 2: return p2_iso_p1_interval_template();
        case 3: return p3_iso_p1_interval_template();
        case 4: return p4_iso_p1_interval_template();
        default:
            throw std::invalid_argument("interval_iso_p1_template supports orders 2, 3, 4");
    }
}

const RefinementTemplate& p2_iso_p1_prism_template()
{
    static const RefinementTemplate s = topology_only(iso_refine_generated::prism_p2_vtk_tess_basix_tpl);
    return s;
}

const RefinementTemplate& p3_iso_p1_prism_template()
{
    static const RefinementTemplate s = topology_only(iso_refine_generated::prism_p3_vtk_tess_basix_tpl);
    return s;
}

const RefinementTemplate& p4_iso_p1_prism_template()
{
    static const RefinementTemplate s = topology_only(iso_refine_generated::prism_p4_vtk_tess_basix_tpl);
    return s;
}

const RefinementTemplate& prism_iso_p1_template(int order)
{
    switch (order)
    {
        case 2: return p2_iso_p1_prism_template();
        case 3: return p3_iso_p1_prism_template();
        case 4: return p4_iso_p1_prism_template();
        default:
            throw std::invalid_argument("prism_iso_p1_template supports orders 2, 3, 4");
    }
}

const RefinementTemplate& p2_iso_p1_pyramid_template()
{
    static const RefinementTemplate full = make_pyramid_iso_p1_storage(2);
    static const RefinementTemplate s = topology_only(full);
    return s;
}

const RefinementTemplate& p3_iso_p1_pyramid_template()
{
    static const RefinementTemplate full = make_pyramid_iso_p1_storage(3);
    static const RefinementTemplate s = topology_only(full);
    return s;
}

const RefinementTemplate& p4_iso_p1_pyramid_template()
{
    static const RefinementTemplate full = make_pyramid_iso_p1_storage(4);
    static const RefinementTemplate s = topology_only(full);
    return s;
}

const RefinementTemplate& pyramid_iso_p1_template(int order)
{
    switch (order)
    {
        case 2: return p2_iso_p1_pyramid_template();
        case 3: return p3_iso_p1_pyramid_template();
        case 4: return p4_iso_p1_pyramid_template();
        default:
            throw std::invalid_argument("pyramid_iso_p1_template supports orders 2, 3, 4");
    }
}

const RefinementTemplate& iso_p1_template(cell::type ct, int order)
{
    switch (ct)
    {
        case cell::type::interval:
            return interval_iso_p1_template(order);
        case cell::type::triangle:
            return triangle_iso_p1_template(order);
        case cell::type::tetrahedron:
            return tetrahedron_iso_p1_template(order);
        case cell::type::quadrilateral:
            return quadrilateral_iso_p1_template(order);
        case cell::type::hexahedron:
            return hexahedron_iso_p1_template(order);
        case cell::type::prism:
            return prism_iso_p1_template(order);
        case cell::type::pyramid:
            return pyramid_iso_p1_template(order);
        default:
            break;
    }
    throw std::invalid_argument(
        "iso_p1_template currently supports interval/triangle/tetrahedron/quadrilateral/hexahedron/prism/pyramid for orders 2,3,4");
}

const RefinementTemplate& p1_template(cell::type ct)
{
    switch (ct)
    {
        case cell::type::triangle:
        {
            static const RefinementTemplate full = make_template(
                cell::type::triangle, cell::type::triangle, 2, 3,
                {0.0, 0.0, 1.0, 0.0, 0.0, 1.0},
                {0, 0, 0},
                {0, 1, 2},
                {0, 1, 2});
            static const RefinementTemplate s = topology_only(full);
            return s;
        }
        case cell::type::quadrilateral:
        {
            static const RefinementTemplate full = make_template(
                cell::type::quadrilateral, cell::type::quadrilateral, 2, 4,
                {0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0},
                {0, 0, 0, 0},
                {0, 1, 2, 3},
                {0, 1, 2, 3});
            static const RefinementTemplate s = topology_only(full);
            return s;
        }
        case cell::type::tetrahedron:
        {
            static const RefinementTemplate full = make_template(
                cell::type::tetrahedron, cell::type::tetrahedron, 3, 4,
                {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},
                {0, 0, 0, 0},
                {0, 1, 2, 3},
                {0, 1, 2, 3});
            static const RefinementTemplate s = topology_only(full);
            return s;
        }
        case cell::type::hexahedron:
        {
            static const RefinementTemplate full = make_template(
                cell::type::hexahedron, cell::type::hexahedron, 3, 8,
                {
                    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0,
                    0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0
                },
                {0, 0, 0, 0, 0, 0, 0, 0},
                {0, 1, 2, 3, 4, 5, 6, 7},
                {0, 1, 2, 3, 4, 5, 6, 7});
            static const RefinementTemplate s = topology_only(full);
            return s;
        }
        default:
            throw std::invalid_argument(
                "p1_template currently supports triangle, tetrahedron, quadrilateral, hexahedron");
    }
}

std::span<const double> iso_p1_ref_coords(cell::type ct, int order)
{
    switch (ct)
    {
        case cell::type::triangle:
            switch (order)
            {
                case 2: return iso_refine_generated::triangle_p2_vtk_tess_basix_tpl.ref_vertex_coords;
                case 3: return iso_refine_generated::triangle_p3_vtk_tess_basix_tpl.ref_vertex_coords;
                case 4: return iso_refine_generated::triangle_p4_vtk_tess_basix_tpl.ref_vertex_coords;
                default: break;
            }
            break;
        case cell::type::quadrilateral:
            switch (order)
            {
                case 2: return iso_refine_generated::quadrilateral_p2_vtk_tess_basix_tpl.ref_vertex_coords;
                case 3: return iso_refine_generated::quadrilateral_p3_vtk_tess_basix_tpl.ref_vertex_coords;
                case 4: return iso_refine_generated::quadrilateral_p4_vtk_tess_basix_tpl.ref_vertex_coords;
                default: break;
            }
            break;
        case cell::type::tetrahedron:
            switch (order)
            {
                case 2: return iso_refine_generated::tetrahedron_p2_vtk_tess_basix_tpl.ref_vertex_coords;
                case 3: return iso_refine_generated::tetrahedron_p3_vtk_tess_basix_tpl.ref_vertex_coords;
                case 4: return iso_refine_generated::tetrahedron_p4_vtk_tess_basix_tpl.ref_vertex_coords;
                default: break;
            }
            break;
        case cell::type::hexahedron:
            switch (order)
            {
                case 2: return iso_refine_generated::hexahedron_p2_vtk_tess_basix_tpl.ref_vertex_coords;
                case 3: return iso_refine_generated::hexahedron_p3_vtk_tess_basix_tpl.ref_vertex_coords;
                case 4: return iso_refine_generated::hexahedron_p4_vtk_tess_basix_tpl.ref_vertex_coords;
                default: break;
            }
            break;
        case cell::type::prism:
            switch (order)
            {
                case 2: return iso_refine_generated::prism_p2_vtk_tess_basix_tpl.ref_vertex_coords;
                case 3: return iso_refine_generated::prism_p3_vtk_tess_basix_tpl.ref_vertex_coords;
                case 4: return iso_refine_generated::prism_p4_vtk_tess_basix_tpl.ref_vertex_coords;
                default: break;
            }
            break;
        case cell::type::interval:
            switch (order)
            {
                case 2:
                {
                    static const RefinementTemplate s = make_interval_iso_p1_storage(2);
                    return s.ref_vertex_coords;
                }
                case 3:
                {
                    static const RefinementTemplate s = make_interval_iso_p1_storage(3);
                    return s.ref_vertex_coords;
                }
                case 4:
                {
                    static const RefinementTemplate s = make_interval_iso_p1_storage(4);
                    return s.ref_vertex_coords;
                }
                default: break;
            }
            break;
        case cell::type::pyramid:
            switch (order)
            {
                case 2:
                {
                    static const RefinementTemplate s = make_pyramid_iso_p1_storage(2);
                    return s.ref_vertex_coords;
                }
                case 3:
                {
                    static const RefinementTemplate s = make_pyramid_iso_p1_storage(3);
                    return s.ref_vertex_coords;
                }
                case 4:
                {
                    static const RefinementTemplate s = make_pyramid_iso_p1_storage(4);
                    return s.ref_vertex_coords;
                }
                default: break;
            }
            break;
        default:
            break;
    }
    throw std::invalid_argument("iso_p1_ref_coords supports interval/triangle/tetrahedron/quadrilateral/hexahedron/prism/pyramid for orders 2,3,4");
}

std::span<const double> p1_ref_coords(cell::type ct)
{
    switch (ct)
    {
        case cell::type::interval:
        {
            static const std::vector<double> x = {0.0, 1.0};
            return x;
        }
        case cell::type::triangle:
        {
            static const std::vector<double> x = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0};
            return x;
        }
        case cell::type::quadrilateral:
        {
            static const std::vector<double> x = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0};
            return x;
        }
        case cell::type::tetrahedron:
        {
            static const std::vector<double> x = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
            return x;
        }
        case cell::type::hexahedron:
        {
            static const std::vector<double> x = {
                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0,
                0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0};
            return x;
        }
        default:
            break;
    }
    throw std::invalid_argument(
        "p1_ref_coords currently supports interval, triangle, tetrahedron, quadrilateral, hexahedron");
}

const std::vector<int>& vtk_lagrange_to_basix_permutation(cell::type ct, int order)
{
    struct PermutationCache
    {
        std::array<std::vector<int>, 5> interval;
        std::array<std::vector<int>, 5> triangle;
        std::array<std::vector<int>, 5> quadrilateral;
    };

    static const PermutationCache cache = []()
    {
        PermutationCache out;
        for (int p = 1; p <= 4; ++p)
        {
            out.interval[static_cast<std::size_t>(p)]
                = build_vtk_lagrange_to_basix_permutation(cell::type::interval, p);
            out.triangle[static_cast<std::size_t>(p)]
                = build_vtk_lagrange_to_basix_permutation(cell::type::triangle, p);
            out.quadrilateral[static_cast<std::size_t>(p)]
                = build_vtk_lagrange_to_basix_permutation(cell::type::quadrilateral, p);
        }
        return out;
    }();

    if (order < 1 || order > 4)
        throw std::invalid_argument("vtk_lagrange_to_basix_permutation supports orders 1..4");

    switch (ct)
    {
        case cell::type::interval:
            return cache.interval[static_cast<std::size_t>(order)];
        case cell::type::triangle:
            return cache.triangle[static_cast<std::size_t>(order)];
        case cell::type::quadrilateral:
            return cache.quadrilateral[static_cast<std::size_t>(order)];
        default:
            throw std::invalid_argument(
                "vtk_lagrange_to_basix_permutation supports only interval/triangle/quadrilateral");
    }
}

} // namespace cutcells
