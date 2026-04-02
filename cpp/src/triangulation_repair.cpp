// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "triangulation_repair.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <utility>

namespace cutcells::cell
{

// ============================================================================
// is_triangulation_diagonal
// ============================================================================

bool is_triangulation_diagonal(type parent_cell_type,
                               int32_t token_a,
                               int32_t token_b)
{
    const int n_edges = num_edges(parent_cell_type);
    const int n_verts = get_num_vertices(parent_cell_type);

    const bool a_root = (token_a >= 0 && token_a < n_edges);
    const bool a_orig = (token_a >= 100 && token_a < 100 + n_verts);
    const bool b_root = (token_b >= 0 && token_b < n_edges);
    const bool b_orig = (token_b >= 100 && token_b < 100 + n_verts);

    // Root-to-root: interface edge, not a diagonal
    if (a_root && b_root)
        return false;

    // Root-to-original: diagonal if the original vertex is NOT an endpoint
    // of the parent edge on which the root lives
    if (a_root && b_orig)
    {
        const auto& edge_def
            = edges(parent_cell_type)[static_cast<std::size_t>(token_a)];
        const int pv = token_b - 100;
        return !(pv == edge_def[0] || pv == edge_def[1]);
    }
    if (b_root && a_orig)
    {
        const auto& edge_def
            = edges(parent_cell_type)[static_cast<std::size_t>(token_b)];
        const int pv = token_a - 100;
        return !(pv == edge_def[0] || pv == edge_def[1]);
    }

    // Original-to-original: diagonal if they don't form a parent cell edge
    if (a_orig && b_orig)
    {
        const int pva = token_a - 100;
        const int pvb = token_b - 100;
        for (const auto& e : edges(parent_cell_type))
        {
            if ((pva == e[0] && pvb == e[1]) || (pva == e[1] && pvb == e[0]))
                return false;
        }
        return true;
    }

    // Unknown tokens: treat as not a diagonal
    return false;
}

// ============================================================================
// diagonal_crosses_interface
// ============================================================================

template <std::floating_point T>
bool diagonal_crosses_interface(
    const CutCell<T>& cut_cell,
    int va, int vb,
    std::function<T(const T*)> eval_phi,
    int expected_sign,
    T tol)
{
    const int gdim = cut_cell._gdim;

    // Sample at multiple points along the edge to detect crossings
    // that a single midpoint check would miss.
    static constexpr std::array<double, 5> samples = {0.2, 0.35, 0.5, 0.65, 0.8};

    for (const double t : samples)
    {
        std::array<T, 3> pt = {T(0), T(0), T(0)};
        for (int d = 0; d < gdim; ++d)
        {
            const T xa = cut_cell._vertex_coords[static_cast<std::size_t>(va * gdim + d)];
            const T xb = cut_cell._vertex_coords[static_cast<std::size_t>(vb * gdim + d)];
            pt[static_cast<std::size_t>(d)] = xa + static_cast<T>(t) * (xb - xa);
        }

        const T phi_val = eval_phi(pt.data());

        if (expected_sign < 0 && phi_val > tol)
            return true;
        if (expected_sign > 0 && phi_val < -tol)
            return true;
    }

    return false;
}

// ============================================================================
// swap_diagonal_2d
// ============================================================================

template <std::floating_point T>
bool swap_diagonal_2d(CutCell<T>& cut_cell, int cell_i, int cell_j)
{
    const int n_sub = num_cells(cut_cell);
    if (cell_i < 0 || cell_j < 0 || cell_i >= n_sub || cell_j >= n_sub)
        return false;

    // Both must be triangles
    if (cut_cell._types[static_cast<std::size_t>(cell_i)] != type::triangle
        || cut_cell._types[static_cast<std::size_t>(cell_j)] != type::triangle)
        return false;

    auto vi = cell_vertices(cut_cell, cell_i);
    auto vj = cell_vertices(cut_cell, cell_j);
    if (vi.size() != 3 || vj.size() != 3)
        return false;

    // Find shared (the old diagonal) and unshared vertices
    int shared_a = -1, shared_b = -1;
    int unshared_i = -1, unshared_j = -1;
    int n_shared = 0;

    for (int a = 0; a < 3; ++a)
    {
        bool found = false;
        for (int b = 0; b < 3; ++b)
        {
            if (vi[a] == vj[b])
            {
                found = true;
                if (n_shared == 0) shared_a = vi[a];
                else               shared_b = vi[a];
                n_shared++;
                break;
            }
        }
        if (!found) unshared_i = vi[a];
    }
    for (int b = 0; b < 3; ++b)
    {
        bool found = false;
        for (int a = 0; a < 3; ++a)
        {
            if (vj[b] == vi[a]) { found = true; break; }
        }
        if (!found) unshared_j = vj[b];
    }

    if (n_shared != 2 || shared_a < 0 || shared_b < 0
        || unshared_i < 0 || unshared_j < 0)
        return false;

    // Old diagonal: shared_a — shared_b
    // New diagonal: unshared_i — unshared_j
    //
    // New triangles (preserving winding):
    //   T_i': {shared_a, unshared_i, unshared_j}
    //   T_j': {shared_b, unshared_j, unshared_i}
    const int off_i = cut_cell._offset[static_cast<std::size_t>(cell_i)];
    const int off_j = cut_cell._offset[static_cast<std::size_t>(cell_j)];

    cut_cell._connectivity[static_cast<std::size_t>(off_i + 0)] = shared_a;
    cut_cell._connectivity[static_cast<std::size_t>(off_i + 1)] = unshared_i;
    cut_cell._connectivity[static_cast<std::size_t>(off_i + 2)] = unshared_j;

    cut_cell._connectivity[static_cast<std::size_t>(off_j + 0)] = shared_b;
    cut_cell._connectivity[static_cast<std::size_t>(off_j + 1)] = unshared_j;
    cut_cell._connectivity[static_cast<std::size_t>(off_j + 2)] = unshared_i;

    return true;
}

// ============================================================================
// swap_prism_triangulation_3d
// ============================================================================

template <std::floating_point T>
bool swap_prism_triangulation_3d(CutCell<T>& cut_cell,
                                 int cell_start,
                                 std::span<const int> prism_verts)
{
    if (prism_verts.size() != 6)
        return false;

    const int n_sub = num_cells(cut_cell);
    if (cell_start < 0 || cell_start + 2 >= n_sub)
        return false;

    for (int k = 0; k < 3; ++k)
    {
        if (cut_cell._types[static_cast<std::size_t>(cell_start + k)]
            != type::tetrahedron)
            return false;
    }

    // Prism vertices: v0,v1,v2 = bottom, v3,v4,v5 = top
    const int v0 = prism_verts[0], v1 = prism_verts[1], v2 = prism_verts[2];
    const int v3 = prism_verts[3], v4 = prism_verts[4], v5 = prism_verts[5];

    // Alternative decomposition B:
    //   T0': {v0, v1, v2, v4}
    //   T1': {v0, v4, v2, v5}
    //   T2': {v0, v4, v5, v3}
    const int off0 = cut_cell._offset[static_cast<std::size_t>(cell_start)];
    const int off1 = cut_cell._offset[static_cast<std::size_t>(cell_start + 1)];
    const int off2 = cut_cell._offset[static_cast<std::size_t>(cell_start + 2)];

    cut_cell._connectivity[static_cast<std::size_t>(off0 + 0)] = v0;
    cut_cell._connectivity[static_cast<std::size_t>(off0 + 1)] = v1;
    cut_cell._connectivity[static_cast<std::size_t>(off0 + 2)] = v2;
    cut_cell._connectivity[static_cast<std::size_t>(off0 + 3)] = v4;

    cut_cell._connectivity[static_cast<std::size_t>(off1 + 0)] = v0;
    cut_cell._connectivity[static_cast<std::size_t>(off1 + 1)] = v4;
    cut_cell._connectivity[static_cast<std::size_t>(off1 + 2)] = v2;
    cut_cell._connectivity[static_cast<std::size_t>(off1 + 3)] = v5;

    cut_cell._connectivity[static_cast<std::size_t>(off2 + 0)] = v0;
    cut_cell._connectivity[static_cast<std::size_t>(off2 + 1)] = v4;
    cut_cell._connectivity[static_cast<std::size_t>(off2 + 2)] = v5;
    cut_cell._connectivity[static_cast<std::size_t>(off2 + 3)] = v3;

    return true;
}

// ============================================================================
// repair_cut_cell_diagonals
// ============================================================================

template <std::floating_point T>
TriangulationRepairInfo repair_cut_cell_diagonals(
    CutCell<T>& cut_cell,
    type parent_cell_type,
    std::function<T(const T*)> eval_phi,
    int expected_sign,
    T tol,
    bool debug)
{
    TriangulationRepairInfo info;

    const int n_sub = num_cells(cut_cell);
    if (debug)
        std::cerr << "triangulation_repair: n_sub=" << n_sub << "\n";
    if (n_sub < 2) return info;

    const int n_verts
        = static_cast<int>(cut_cell._vertex_coords.size()) / cut_cell._gdim;
    if (debug)
    {
        std::cerr << "triangulation_repair: n_verts=" << n_verts
                  << " pe_size=" << cut_cell._vertex_parent_entity.size()
                  << " coords_size=" << cut_cell._vertex_coords.size()
                  << " gdim=" << cut_cell._gdim << "\n";
    }
    if (static_cast<int>(cut_cell._vertex_parent_entity.size()) != n_verts)
        return info;

    // Build edge → cell adjacency (edge key is ordered pair of vertex ids)
    using EdgeKey = std::pair<int, int>;
    std::map<EdgeKey, std::vector<int>> edge_to_cells;

    for (int c = 0; c < n_sub; ++c)
    {
        auto verts = cell_vertices(cut_cell, c);
        const auto ct = cut_cell._types[static_cast<std::size_t>(c)];
        const auto& cell_edge_defs = edges(ct);
        for (const auto& e : cell_edge_defs)
        {
            int v0 = verts[static_cast<std::size_t>(e[0])];
            int v1 = verts[static_cast<std::size_t>(e[1])];
            if (v0 > v1) std::swap(v0, v1);
            edge_to_cells[{v0, v1}].push_back(c);
        }
    }

    if (debug)
    {
        std::cerr << "triangulation_repair: built " << edge_to_cells.size()
                  << " unique edges\n";
        for (const auto& [ek, cells] : edge_to_cells)
        {
            auto [va, vb] = ek;
            std::cerr << "  edge (" << va << "," << vb << ") tok=("
                      << cut_cell._vertex_parent_entity[static_cast<std::size_t>(va)]
                      << ","
                      << cut_cell._vertex_parent_entity[static_cast<std::size_t>(vb)]
                      << ") cells=" << cells.size()
                      << "\n";
        }
    }

    // Find internal edges shared by exactly two sub-cells and check them.
    // An edge shared by 2 sub-cells is an interior diagonal of the triangulation
    // (boundary edges of the CutCell appear in only 1 sub-cell).
    // Skip interface edges (root-to-root) since they belong on the cut surface.
    const int n_pedges = num_edges(parent_cell_type);
    for (auto& [ek, cells] : edge_to_cells)
    {
        if (cells.size() != 2) continue;

        auto [va, vb] = ek;
        const int32_t tok_a
            = cut_cell._vertex_parent_entity[static_cast<std::size_t>(va)];
        const int32_t tok_b
            = cut_cell._vertex_parent_entity[static_cast<std::size_t>(vb)];

        // Skip interface edges (both endpoints are root vertices)
        const bool a_root = (tok_a >= 0 && tok_a < n_pedges);
        const bool b_root = (tok_b >= 0 && tok_b < n_pedges);
        if (a_root && b_root)
            continue;

        info.checked_diagonals++;

        if (!diagonal_crosses_interface(cut_cell, va, vb,
                                        eval_phi, expected_sign, tol))
            continue;

        info.invalid_diagonals++;
        if (debug)
        {
            std::cerr << "triangulation_repair: diagonal (" << va << ","
                      << vb << ") crosses interface, attempting swap\n";
        }

        // ---- 2D swap (two triangles → two triangles) -----------------
        if (cut_cell._tdim == 2
            && cut_cell._types[static_cast<std::size_t>(cells[0])] == type::triangle
            && cut_cell._types[static_cast<std::size_t>(cells[1])] == type::triangle)
        {
            // Swap the diagonal
            if (!swap_diagonal_2d(cut_cell, cells[0], cells[1]))
            {
                info.unresolved++;
                continue;
            }

            // Find the new shared edge (the new diagonal)
            auto vi_new = cell_vertices(cut_cell, cells[0]);
            auto vj_new = cell_vertices(cut_cell, cells[1]);
            int new_a = -1, new_b = -1;
            int ns = 0;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    if (vi_new[i] == vj_new[j])
                    {
                        if (ns == 0) new_a = vi_new[i];
                        else         new_b = vi_new[i];
                        ns++;
                    }

            if (ns == 2 && !diagonal_crosses_interface(
                    cut_cell, new_a, new_b, eval_phi, expected_sign, tol))
            {
                info.swapped_diagonals++;
                if (debug)
                {
                    std::cerr << "triangulation_repair: swap succeeded, "
                              << "new diagonal (" << new_a << "," << new_b
                              << ") is valid\n";
                }
            }
            else
            {
                // New diagonal also bad — revert
                swap_diagonal_2d(cut_cell, cells[0], cells[1]);
                info.unresolved++;
                if (debug)
                {
                    std::cerr << "triangulation_repair: swap failed "
                              << "(new diagonal also crosses), "
                              << "marking as unresolved\n";
                }
            }
        }
        else
        {
            // 3D or non-triangle: mark as unresolved for now
            info.unresolved++;
            if (debug)
            {
                std::cerr << "triangulation_repair: 3D diagonal repair "
                          << "not yet automatic, marking as unresolved\n";
            }
        }
    }

    return info;
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template bool diagonal_crosses_interface<float>(
    const CutCell<float>&, int, int,
    std::function<float(const float*)>, int, float);
template bool diagonal_crosses_interface<double>(
    const CutCell<double>&, int, int,
    std::function<double(const double*)>, int, double);

template bool swap_diagonal_2d<float>(CutCell<float>&, int, int);
template bool swap_diagonal_2d<double>(CutCell<double>&, int, int);

template bool swap_prism_triangulation_3d<float>(
    CutCell<float>&, int, std::span<const int>);
template bool swap_prism_triangulation_3d<double>(
    CutCell<double>&, int, std::span<const int>);

template TriangulationRepairInfo repair_cut_cell_diagonals<float>(
    CutCell<float>&, type, std::function<float(const float*)>, int, float, bool);
template TriangulationRepairInfo repair_cut_cell_diagonals<double>(
    CutCell<double>&, type, std::function<double(const double*)>, int, double, bool);

} // namespace cutcells::cell
