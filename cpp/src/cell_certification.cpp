// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "cell_certification.h"
#include "bernstein.h"
#include "cell_flags.h"
#include "cell_topology.h"
#include "cut_cell.h"
#include "cut_interval.h"
#include "cut_tetrahedron.h"
#include "cut_triangle.h"
#include "edge_certification.h"
#include "mapping.h"
#include "reference_cell.h"
#include "refine_cell.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace cutcells
{
namespace
{

template <std::floating_point T>
void set_vertex_sign_for_level_set(AdaptCell<T>& adapt_cell,
                                   int vertex_id,
                                   int level_set_id,
                                   T value,
                                   T zero_tol)
{
    const std::uint64_t bit = std::uint64_t(1) << level_set_id;
    auto& zero_mask = adapt_cell.zero_mask_per_vertex[static_cast<std::size_t>(vertex_id)];
    auto& negative_mask = adapt_cell.negative_mask_per_vertex[static_cast<std::size_t>(vertex_id)];

    zero_mask &= ~bit;
    negative_mask &= ~bit;

    if (std::fabs(value) <= zero_tol)
        zero_mask |= bit;
    else if (value < T(0))
        negative_mask |= bit;
}

template <std::floating_point T>
std::map<std::pair<int, int>, int> build_leaf_edge_lookup(const AdaptCell<T>& adapt_cell)
{
    std::map<std::pair<int, int>, int> edge_lookup;
    const int n_edges = adapt_cell.n_entities(1);
    for (int e = 0; e < n_edges; ++e)
    {
        auto ev = adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(e)];
        edge_lookup[{std::min(static_cast<int>(ev[0]), static_cast<int>(ev[1])),
                     std::max(static_cast<int>(ev[0]), static_cast<int>(ev[1]))}] = e;
    }
    return edge_lookup;
}

void append_top_cell_local(std::vector<cell::type>& types,
                           EntityAdjacency& adj,
                           cell::type ctype,
                           std::span<const int> verts)
{
    types.push_back(ctype);
    for (int v : verts)
        adj.indices.push_back(static_cast<std::int32_t>(v));
    adj.offsets.push_back(static_cast<std::int32_t>(adj.indices.size()));
}

int basix_edge_id_for_vertices(cell::type cell_type, int a, int b)
{
    const auto ledges = cell::edges(cell_type);
    for (std::size_t e = 0; e < ledges.size(); ++e)
    {
        const auto edge = ledges[e];
        if ((edge[0] == a && edge[1] == b) || (edge[0] == b && edge[1] == a))
            return static_cast<int>(e);
    }

    throw std::runtime_error("basix_edge_id_for_vertices: edge not found");
}

template <std::floating_point T>
std::vector<int> canonicalize_cut_part_cell_vertices(
    const cell::CutCell<T>& part,
    int part_cell_id,
    cell::type leaf_cell_type)
{
    auto part_vertices = cell::cell_vertices(part, part_cell_id);
    std::vector<int> canonical(part_vertices.begin(), part_vertices.end());

    const cell::type part_cell_type =
        part._types[static_cast<std::size_t>(part_cell_id)];

    if (leaf_cell_type == cell::type::triangle
        && part_cell_type == cell::type::quadrilateral)
    {
        std::map<int, int> token_to_part_vertex;
        std::array<bool, 3> present = {false, false, false};

        for (int part_lv : part_vertices)
        {
            const int token = cell::vtk_parent_entity_token_to_basix(
                leaf_cell_type,
                part._vertex_parent_entity[static_cast<std::size_t>(part_lv)]);
            token_to_part_vertex[token] = part_lv;
            if (token >= 100 && token < 103)
                present[static_cast<std::size_t>(token - 100)] = true;
        }

        int omitted = -1;
        for (int v = 0; v < 3; ++v)
        {
            if (!present[static_cast<std::size_t>(v)])
            {
                omitted = v;
                break;
            }
        }
        if (omitted < 0)
            throw std::runtime_error("canonicalize_cut_part_cell_vertices: invalid triangle quad");

        const auto kept_edge = cell::edges(cell::type::triangle)[static_cast<std::size_t>(omitted)];
        const int a = kept_edge[0];
        const int b = kept_edge[1];
        const int e_a = basix_edge_id_for_vertices(cell::type::triangle, a, omitted);
        const int e_b = basix_edge_id_for_vertices(cell::type::triangle, b, omitted);

        canonical = {
            token_to_part_vertex.contains(100 + a) ? token_to_part_vertex.at(100 + a)
                                                   : throw std::runtime_error("canonicalize_cut_part_cell_vertices: missing triangle vertex token a"),
            token_to_part_vertex.contains(100 + b) ? token_to_part_vertex.at(100 + b)
                                                   : throw std::runtime_error("canonicalize_cut_part_cell_vertices: missing triangle vertex token b"),
            token_to_part_vertex.contains(e_a) ? token_to_part_vertex.at(e_a)
                                               : throw std::runtime_error("canonicalize_cut_part_cell_vertices: missing triangle edge token ea"),
            token_to_part_vertex.contains(e_b) ? token_to_part_vertex.at(e_b)
                                               : throw std::runtime_error("canonicalize_cut_part_cell_vertices: missing triangle edge token eb")};
    }
    else if (leaf_cell_type == cell::type::tetrahedron
             && part_cell_type == cell::type::prism)
    {
        std::map<int, int> token_to_part_vertex;
        std::array<bool, 4> present = {false, false, false, false};

        for (int part_lv : part_vertices)
        {
            const int token = cell::vtk_parent_entity_token_to_basix(
                leaf_cell_type,
                part._vertex_parent_entity[static_cast<std::size_t>(part_lv)]);
            token_to_part_vertex[token] = part_lv;
            if (token >= 100 && token < 104)
                present[static_cast<std::size_t>(token - 100)] = true;
        }

        std::vector<int> present_vertices;
        std::vector<int> missing_vertices;
        for (int v = 0; v < 4; ++v)
        {
            if (present[static_cast<std::size_t>(v)])
                present_vertices.push_back(v);
            else
                missing_vertices.push_back(v);
        }

        if (present_vertices.size() == 3)
        {
            const int omitted = missing_vertices.front();
            const auto face = cell::tet_face(omitted);
            canonical.resize(6);
            for (int i = 0; i < 3; ++i)
            {
                const int fv = face[static_cast<std::size_t>(i)];
                const int edge_id = basix_edge_id_for_vertices(
                    cell::type::tetrahedron, fv, omitted);
                auto v_it = token_to_part_vertex.find(100 + fv);
                auto e_it = token_to_part_vertex.find(edge_id);
                if (v_it == token_to_part_vertex.end() || e_it == token_to_part_vertex.end())
                {
                    throw std::runtime_error(
                        "canonicalize_cut_part_cell_vertices: invalid 3-vertex prism token set");
                }
                canonical[static_cast<std::size_t>(i)] = v_it->second;
                canonical[static_cast<std::size_t>(3 + i)] = e_it->second;
            }
        }
        else if (present_vertices.size() == 2)
        {
            const int p = present_vertices[0];
            const int q = present_vertices[1];
            const int r = missing_vertices[0];
            const int s = missing_vertices[1];

            const int e_pr = basix_edge_id_for_vertices(cell::type::tetrahedron, p, r);
            const int e_ps = basix_edge_id_for_vertices(cell::type::tetrahedron, p, s);
            const int e_qr = basix_edge_id_for_vertices(cell::type::tetrahedron, q, r);
            const int e_qs = basix_edge_id_for_vertices(cell::type::tetrahedron, q, s);

            const std::array<int, 6> prism_tokens = {
                100 + p, e_pr, e_ps,
                100 + q, e_qr, e_qs};
            canonical.resize(6);
            for (int i = 0; i < 6; ++i)
            {
                auto it = token_to_part_vertex.find(prism_tokens[static_cast<std::size_t>(i)]);
                if (it == token_to_part_vertex.end())
                {
                    throw std::runtime_error(
                        "canonicalize_cut_part_cell_vertices: invalid 2-vertex prism token set");
                }
                canonical[static_cast<std::size_t>(i)] = it->second;
            }
        }
        else
        {
            throw std::runtime_error(
                "canonicalize_cut_part_cell_vertices: unsupported tetra prism token pattern");
        }
    }

    return canonical;
}

template <std::floating_point T, std::integral I>
std::vector<T> gather_leaf_cell_vertex_level_set_values(
    const AdaptCell<T>& adapt_cell,
    const LevelSetCell<T, I>& ls_cell,
    int cell_id)
{
    const int tdim = adapt_cell.tdim;
    auto cell_verts = adapt_cell.entity_to_vertex[tdim][static_cast<std::int32_t>(cell_id)];
    std::vector<T> values(cell_verts.size(), T(0));

    for (std::size_t j = 0; j < cell_verts.size(); ++j)
    {
        const int gv = static_cast<int>(cell_verts[j]);
        std::span<const T> xi(
            adapt_cell.vertex_coords.data()
                + static_cast<std::size_t>(gv) * static_cast<std::size_t>(tdim),
            static_cast<std::size_t>(tdim));
        values[j] = ls_cell.value(xi);
    }

    return values;
}

template <std::floating_point T>
int append_vertex_with_parent_info(AdaptCell<T>& adapt_cell,
                                   std::span<const T> coords,
                                   int parent_dim,
                                   int parent_id,
                                   std::span<const T> parent_param,
                                   int source_edge_id)
{
    const int new_v = adapt_cell.n_vertices();
    adapt_cell.vertex_coords.insert(adapt_cell.vertex_coords.end(),
                                    coords.begin(), coords.end());
    adapt_cell.vertex_parent_dim.push_back(static_cast<std::int8_t>(parent_dim));
    adapt_cell.vertex_parent_id.push_back(parent_id);
    const std::int32_t param_offset =
        static_cast<std::int32_t>(adapt_cell.vertex_parent_param.size());
    adapt_cell.vertex_parent_param.insert(adapt_cell.vertex_parent_param.end(),
                                          parent_param.begin(),
                                          parent_param.end());
    adapt_cell.vertex_parent_param_offset.push_back(
        param_offset + static_cast<std::int32_t>(parent_param.size()));
    adapt_cell.zero_mask_per_vertex.push_back(0);
    adapt_cell.negative_mask_per_vertex.push_back(0);
    adapt_cell.vertex_source_edge_id.push_back(source_edge_id);
    return new_v;
}

template <std::floating_point T>
void inherit_common_edge_sign_masks(AdaptCell<T>& adapt_cell,
                                    int vertex_id,
                                    int edge_id)
{
    if (edge_id < 0)
        return;

    auto edge_vertices = adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(edge_id)];
    if (edge_vertices.size() != 2)
        return;

    const auto v0 = static_cast<std::size_t>(edge_vertices[0]);
    const auto v1 = static_cast<std::size_t>(edge_vertices[1]);
    const std::uint64_t common_zero =
        adapt_cell.zero_mask_per_vertex[v0] & adapt_cell.zero_mask_per_vertex[v1];
    const std::uint64_t common_negative =
        (adapt_cell.negative_mask_per_vertex[v0]
         & adapt_cell.negative_mask_per_vertex[v1])
        & ~common_zero;

    auto& zero_mask = adapt_cell.zero_mask_per_vertex[static_cast<std::size_t>(vertex_id)];
    auto& negative_mask = adapt_cell.negative_mask_per_vertex[static_cast<std::size_t>(vertex_id)];
    zero_mask |= common_zero;
    negative_mask &= ~common_zero;
    negative_mask |= common_negative;
}

template <std::floating_point T, std::integral I>
void gather_adapt_edge_bernstein(const AdaptCell<T>& adapt_cell,
                                 const LevelSetCell<T, I>& ls_cell,
                                 int edge_id,
                                 std::vector<T>& edge_coeffs)
{
    auto verts = adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(edge_id)];
    const int tdim = adapt_cell.tdim;
    std::span<const T> xi_a(
        adapt_cell.vertex_coords.data()
            + static_cast<std::size_t>(verts[0]) * static_cast<std::size_t>(tdim),
        static_cast<std::size_t>(tdim));
    std::span<const T> xi_b(
        adapt_cell.vertex_coords.data()
            + static_cast<std::size_t>(verts[1]) * static_cast<std::size_t>(tdim),
        static_cast<std::size_t>(tdim));

    restrict_edge_bernstein_exact(
        ls_cell.cell_type, ls_cell.bernstein_order,
        std::span<const T>(ls_cell.bernstein_coeffs),
        xi_a, xi_b, edge_coeffs);
}

template <std::floating_point T>
bool linear_one_root_parameter_from_endpoint_values(T phi0,
                                                    T phi1,
                                                    T zero_tol,
                                                    T& root_t)
{
    const bool left_zero = std::fabs(phi0) <= zero_tol;
    const bool right_zero = std::fabs(phi1) <= zero_tol;

    root_t = T(0);
    if (left_zero && right_zero)
        return false;
    if (left_zero)
        return true;
    if (right_zero)
    {
        root_t = T(1);
        return true;
    }

    const bool brackets_zero = (phi0 < T(0) && phi1 > T(0))
                            || (phi0 > T(0) && phi1 < T(0));
    if (!brackets_zero)
        return false;

    const T denom = phi0 - phi1;
    const T denom_tol =
        std::max(zero_tol, T(64) * std::numeric_limits<T>::epsilon());
    if (std::fabs(denom) <= denom_tol)
        return false;

    root_t = std::clamp(phi0 / denom, T(0), T(1));
    return true;
}

/// Returns true if the Bernstein expansion has a sign-definite partial
/// derivative in at least one reference direction on the subcell.
/// A monotone function with no boundary root cannot have an interior zero.
template <std::floating_point T>
bool has_monotone_direction(cell::type subcell_type, int degree,
                            std::span<const T> subcell_coeffs, T sign_tol)
{
    const int tdim = cell::get_tdim(subcell_type);
    std::vector<T> deriv_coeffs;
    for (int dir = 0; dir < tdim; ++dir)
    {
        bernstein::derivative_coefficients(subcell_type, degree,
                                           subcell_coeffs, dir, deriv_coeffs);
        if (deriv_coeffs.empty())
            continue;
        std::span<const T> dc(deriv_coeffs);
        if (bernstein_all_positive(dc, sign_tol)
            || bernstein_all_negative(dc, sign_tol))
        {
            return true;
        }
    }
    return false;
}

/// Check whether any edge incident to cell_id has a root (one_root, multiple_roots, or zero).
template <std::floating_point T>
bool cell_has_any_edge_root(const AdaptCell<T>& adapt_cell,
                            int level_set_id, int cell_id)
{
    const int tdim = adapt_cell.tdim;
    auto cell_verts = adapt_cell.entity_to_vertex[tdim][static_cast<std::int32_t>(cell_id)];
    const int n_edges = adapt_cell.n_entities(1);

    for (int e = 0; e < n_edges; ++e)
    {
        auto ev = adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(e)];
        bool v0_in = false, v1_in = false;
        for (auto cv : cell_verts)
        {
            if (cv == ev[0]) v0_in = true;
            if (cv == ev[1]) v1_in = true;
        }
        if (!v0_in || !v1_in)
            continue;

        const EdgeRootTag etag = adapt_cell.get_edge_root_tag(level_set_id, e);
        if (etag == EdgeRootTag::one_root
            || etag == EdgeRootTag::multiple_roots
            || etag == EdgeRootTag::zero)
        {
            return true;
        }
    }
    return false;
}

template <std::floating_point T>
bool vertex_parameter_on_parent_edge(const AdaptCell<T>& adapt_cell,
                                     int vertex_id,
                                     int parent_edge_id,
                                     T& t)
{
    const int dim = adapt_cell.vertex_parent_dim[static_cast<std::size_t>(vertex_id)];
    if (dim == 0)
    {
        const auto parent_edge =
            cell::edges(adapt_cell.parent_cell_type)[static_cast<std::size_t>(parent_edge_id)];
        if (parent_edge[0] == vertex_id)
        {
            t = T(0);
            return true;
        }
        if (parent_edge[1] == vertex_id)
        {
            t = T(1);
            return true;
        }
        return false;
    }

    if (dim != 1)
        return false;
    if (adapt_cell.vertex_parent_id[static_cast<std::size_t>(vertex_id)] != parent_edge_id)
        return false;

    const auto begin =
        static_cast<std::size_t>(adapt_cell.vertex_parent_param_offset[static_cast<std::size_t>(vertex_id)]);
    const auto end =
        static_cast<std::size_t>(adapt_cell.vertex_parent_param_offset[static_cast<std::size_t>(vertex_id + 1)]);
    if (end - begin != 1)
        return false;

    t = adapt_cell.vertex_parent_param[begin];
    return true;
}

template <std::floating_point T, std::integral I>
int ensure_one_root_vertex_on_edge(AdaptCell<T>& adapt_cell,
                                   const LevelSetCell<T, I>& ls_cell,
                                   int level_set_id,
                                   int edge_id,
                                   T zero_tol,
                                   T sign_tol,
                                   int edge_max_depth)
{
    const int n_edges = adapt_cell.n_entities(1);
    if (adapt_cell.edge_one_root_has_value.size()
        < static_cast<std::size_t>((level_set_id + 1) * n_edges))
    {
        adapt_cell.resize_one_root_data(level_set_id + 1);
    }

    const auto idx = static_cast<std::size_t>(level_set_id * n_edges + edge_id);
    if (adapt_cell.edge_one_root_vertex_id[idx] >= 0)
        return adapt_cell.edge_one_root_vertex_id[idx];

    if (adapt_cell.get_edge_root_tag(level_set_id, edge_id) != EdgeRootTag::one_root)
        return -1;

    (void)sign_tol;
    (void)edge_max_depth;

    T root_t = T(0);
    bool has_linear_root = adapt_cell.edge_one_root_has_value[idx] != 0;
    if (has_linear_root)
    {
        root_t = adapt_cell.edge_one_root_param[idx];
    }
    else
    {
        std::vector<T> edge_coeffs;
        gather_adapt_edge_bernstein(adapt_cell, ls_cell, edge_id, edge_coeffs);
        has_linear_root =
            edge_coeffs.size() >= 2
            && linear_one_root_parameter_from_endpoint_values<T>(
                edge_coeffs.front(), edge_coeffs.back(), zero_tol, root_t);
    }

    if (!has_linear_root)
    {
        throw std::runtime_error(
            "ensure_one_root_vertex_on_edge: failed to compute linear edge intersection");
    }

    auto verts = adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(edge_id)];
    const int v0 = static_cast<int>(verts[0]);
    const int v1 = static_cast<int>(verts[1]);
    const std::uint64_t bit = std::uint64_t(1) << level_set_id;
    const T endpoint_tol = std::max(zero_tol, T(64) * std::numeric_limits<T>::epsilon());
    if (std::fabs(root_t) <= endpoint_tol
        && (adapt_cell.zero_mask_per_vertex[static_cast<std::size_t>(v0)] & bit) != 0)
    {
        adapt_cell.edge_one_root_param[idx] = T(0);
        adapt_cell.edge_one_root_vertex_id[idx] = v0;
        adapt_cell.edge_one_root_has_value[idx] = 1;
        return v0;
    }
    if (std::fabs(root_t - T(1)) <= endpoint_tol
        && (adapt_cell.zero_mask_per_vertex[static_cast<std::size_t>(v1)] & bit) != 0)
    {
        adapt_cell.edge_one_root_param[idx] = T(1);
        adapt_cell.edge_one_root_vertex_id[idx] = v1;
        adapt_cell.edge_one_root_has_value[idx] = 1;
        return v1;
    }

    const int tdim = adapt_cell.tdim;
    std::vector<T> coords(static_cast<std::size_t>(tdim), T(0));
    for (int d = 0; d < tdim; ++d)
    {
        const T x0 = adapt_cell.vertex_coords[static_cast<std::size_t>(v0 * tdim + d)];
        const T x1 = adapt_cell.vertex_coords[static_cast<std::size_t>(v1 * tdim + d)];
        coords[static_cast<std::size_t>(d)] = (T(1) - root_t) * x0 + root_t * x1;
    }

    int parent_dim = adapt_cell.tdim;
    int parent_id = adapt_cell.parent_cell_id;
    std::vector<T> parent_param(coords.begin(), coords.end());

    int parent_edge_id = -1;
    T parent_t0 = T(0);
    T parent_t1 = T(1);
    if (edge_is_on_single_parent_edge(adapt_cell, edge_id, parent_edge_id)
        && vertex_parameter_on_parent_edge(adapt_cell, v0, parent_edge_id, parent_t0)
        && vertex_parameter_on_parent_edge(adapt_cell, v1, parent_edge_id, parent_t1))
    {
        parent_dim = 1;
        parent_id = parent_edge_id;
        parent_param.assign(
            1, (T(1) - root_t) * parent_t0 + root_t * parent_t1);
    }

    const int vertex_id = append_vertex_with_parent_info(
        adapt_cell,
        std::span<const T>(coords),
        parent_dim,
        parent_id,
        std::span<const T>(parent_param),
        edge_id);

    inherit_common_edge_sign_masks(adapt_cell, vertex_id, edge_id);
    set_vertex_sign_for_level_set(adapt_cell, vertex_id, level_set_id, T(0), zero_tol);
    adapt_cell.edge_one_root_param[idx] = root_t;
    adapt_cell.edge_one_root_vertex_id[idx] = vertex_id;
    adapt_cell.edge_one_root_has_value[idx] = 1;
    return vertex_id;
}

template <std::floating_point T, std::integral I>
void append_ready_cut_part_cells(
    AdaptCell<T>& adapt_cell,
    const LevelSetCell<T, I>& ls_cell,
    int level_set_id,
    int old_cell_id,
    const cell::CutCell<T>& part,
    cell::type leaf_cell_type,
    std::span<const std::int32_t> old_cell_vertices,
    std::span<const int> old_edge_ids_by_local_edge,
    std::vector<cell::type>& new_types,
    EntityAdjacency& new_cells,
    std::vector<int>& old_cell_ids_for_new_cells,
    std::vector<int>& source_cell_ids_for_new_cells,
    std::vector<CellRefinementReason>& refinement_reasons_for_new_cells,
    std::vector<CellCertTag>& explicit_current_ls_tags,
    std::map<int, int>& token_to_vertex,
    T zero_tol,
    CellCertTag side_tag)
{
    const int parent_tdim = adapt_cell.tdim;

    for (std::size_t lv = 0; lv < part._vertex_parent_entity.size(); ++lv)
    {
        const int raw_token = part._vertex_parent_entity[lv];
        const int token = cell::vtk_parent_entity_token_to_basix(leaf_cell_type, raw_token);
        if (token_to_vertex.contains(token))
            continue;

        int vertex_id = -1;
        if (token >= 100 && token < 200)
        {
            const int local_vid = token - 100;
            vertex_id = old_cell_vertices[static_cast<std::size_t>(local_vid)];
        }
        else
        {
            std::span<const T> x(
                part._vertex_coords.data()
                    + lv * static_cast<std::size_t>(parent_tdim),
                static_cast<std::size_t>(parent_tdim));

            int parent_dim = parent_tdim;
            int parent_id = -1;
            std::vector<T> parent_param;

            int old_edge_id_for_token = -1;
            if (token >= 0 && token < 100)
            {
                const int local_edge_id = token;
                if (local_edge_id >= 0
                    && local_edge_id < static_cast<int>(old_edge_ids_by_local_edge.size()))
                {
                    const int old_edge_id = old_edge_ids_by_local_edge[static_cast<std::size_t>(local_edge_id)];
                    old_edge_id_for_token = old_edge_id;
                    int parent_edge_id = -1;
                    if (old_edge_id >= 0
                        && edge_is_on_single_parent_edge(adapt_cell, old_edge_id, parent_edge_id))
                    {
                        parent_dim = 1;
                        parent_id = parent_edge_id;

                        auto edge_vertices =
                            adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(old_edge_id)];
                        T parent_t0 = T(0);
                        T parent_t1 = T(1);
                        if (edge_vertices.size() == 2
                            && vertex_parameter_on_parent_edge(
                                adapt_cell, edge_vertices[0], parent_edge_id, parent_t0)
                            && vertex_parameter_on_parent_edge(
                                adapt_cell, edge_vertices[1], parent_edge_id, parent_t1))
                        {
                            parent_param.assign(1, T(0.5) * (parent_t0 + parent_t1));
                        }
                    }
                }
            }
            if (parent_dim == parent_tdim && parent_id < 0)
            {
                parent_id = adapt_cell.parent_cell_id;
                parent_param.assign(x.begin(), x.end());
            }

            vertex_id = append_vertex_with_parent_info(
                adapt_cell, x, parent_dim, parent_id,
                std::span<const T>(parent_param),
                old_edge_id_for_token);

            if (token >= 0 && token < 100)
            {
                const int local_edge_id = token;
                if (old_edge_id_for_token >= 0)
                    inherit_common_edge_sign_masks(adapt_cell, vertex_id, old_edge_id_for_token);

                const EdgeRootTag edge_tag =
                    old_edge_id_for_token >= 0
                        ? adapt_cell.get_edge_root_tag(level_set_id, old_edge_id_for_token)
                        : EdgeRootTag::not_classified;
                const bool is_level_set_root =
                    local_edge_id >= 0
                    && local_edge_id < static_cast<int>(old_edge_ids_by_local_edge.size())
                    && edge_tag == EdgeRootTag::one_root;
                const bool is_level_set_zero_edge =
                    local_edge_id >= 0
                    && local_edge_id < static_cast<int>(old_edge_ids_by_local_edge.size())
                    && edge_tag == EdgeRootTag::zero;

                if (is_level_set_root || is_level_set_zero_edge)
                {
                    // Straight cut vertices lie on the straight interface by
                    // construction. They should participate in phi=0
                    // extraction even when the polynomial level set evaluated at
                    // this affine point is not exactly zero.
                    set_vertex_sign_for_level_set(
                        adapt_cell, vertex_id, level_set_id, T(0), zero_tol);
                }
                else
                {
                    // LUT triangulation may insert midpoints on uncut parent
                    // edges. They inherit parent-edge provenance, but are not
                    // zero vertices for this level set.
                    const T value = ls_cell.value(x);
                    set_vertex_sign_for_level_set(
                        adapt_cell, vertex_id, level_set_id, value, zero_tol);
                }
            }
            else
            {
                const T value = ls_cell.value(x);
                set_vertex_sign_for_level_set(
                    adapt_cell, vertex_id, level_set_id, value, zero_tol);
            }
        }

        token_to_vertex[token] = vertex_id;
    }

    const int n_part_cells = cell::num_cells(part);
    for (int pc = 0; pc < n_part_cells; ++pc)
    {
        const std::vector<int> part_vertices =
            canonicalize_cut_part_cell_vertices(part, pc, leaf_cell_type);
        std::vector<int> mapped(part_vertices.size(), -1);
        for (std::size_t j = 0; j < part_vertices.size(); ++j)
        {
            const int part_lv = part_vertices[j];
            const int raw_token = part._vertex_parent_entity[static_cast<std::size_t>(part_lv)];
            const int token = cell::vtk_parent_entity_token_to_basix(leaf_cell_type, raw_token);
            auto token_it = token_to_vertex.find(token);
            if (token_it == token_to_vertex.end())
            {
                throw std::runtime_error(
                    "append_ready_cut_part_cells: missing mapped token "
                    + std::to_string(token) + " raw_token=" + std::to_string(raw_token)
                    + " leaf_type=" + cell_type_to_str(leaf_cell_type)
                    + " part_cell_type="
                    + cell_type_to_str(part._types[static_cast<std::size_t>(pc)]));
            }
            mapped[j] = token_it->second;
        }

        append_top_cell_local(new_types, new_cells,
                              part._types[static_cast<std::size_t>(pc)],
                              std::span<const int>(mapped));
        old_cell_ids_for_new_cells.push_back(-1);
        source_cell_ids_for_new_cells.push_back(old_cell_id);
        refinement_reasons_for_new_cells.push_back(CellRefinementReason::cut_level_set);
        explicit_current_ls_tags.push_back(side_tag);
    }
}

template <std::floating_point T>
CellCertTag classify_ready_to_cut_topology(const AdaptCell<T>& adapt_cell,
                                           int level_set_id,
                                           int cell_id)
{
    const int tdim = adapt_cell.tdim;
    auto cell_verts = adapt_cell.entity_to_vertex[tdim][static_cast<std::int32_t>(cell_id)];
    const cell::type subcell_type = adapt_cell.entity_types[tdim][static_cast<std::size_t>(cell_id)];

    std::set<int> cut_point_tokens;
    bool has_multiple_roots = false;
    bool has_zero_edge = false;
    const std::uint64_t bit = std::uint64_t(1) << level_set_id;

    std::map<int, int> local_vertex_by_global;
    for (std::size_t lv = 0; lv < cell_verts.size(); ++lv)
    {
        const int gv = static_cast<int>(cell_verts[lv]);
        local_vertex_by_global[gv] = static_cast<int>(lv);
        if ((adapt_cell.zero_mask_per_vertex[static_cast<std::size_t>(gv)] & bit) != 0)
            cut_point_tokens.insert(100 + static_cast<int>(lv));
    }

    const int n_edges = adapt_cell.n_entities(1);
    for (int e = 0; e < n_edges; ++e)
    {
        auto edge_verts = adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(e)];
        auto it0 = local_vertex_by_global.find(static_cast<int>(edge_verts[0]));
        auto it1 = local_vertex_by_global.find(static_cast<int>(edge_verts[1]));
        if (it0 == local_vertex_by_global.end() || it1 == local_vertex_by_global.end())
            continue;

        const EdgeRootTag etag = adapt_cell.get_edge_root_tag(level_set_id, e);
        if (etag == EdgeRootTag::multiple_roots)
        {
            has_multiple_roots = true;
        }
        else if (etag == EdgeRootTag::zero)
        {
            has_zero_edge = true;
        }
        else if (etag == EdgeRootTag::one_root)
        {
            const bool v0_zero =
                (adapt_cell.zero_mask_per_vertex[static_cast<std::size_t>(edge_verts[0])] & bit) != 0;
            const bool v1_zero =
                (adapt_cell.zero_mask_per_vertex[static_cast<std::size_t>(edge_verts[1])] & bit) != 0;

            if (v0_zero && !v1_zero)
                cut_point_tokens.insert(100 + it0->second);
            else if (v1_zero && !v0_zero)
                cut_point_tokens.insert(100 + it1->second);
            else if (v0_zero && v1_zero)
                has_zero_edge = true;
            else
            {
                const int local_edge_id = basix_edge_id_for_vertices(
                    subcell_type, it0->second, it1->second);
                cut_point_tokens.insert(local_edge_id);
            }
        }
    }

    if (has_multiple_roots)
        return CellCertTag::cut;

    if (subcell_type == cell::type::interval && cut_point_tokens.size() == 1)
        return CellCertTag::ready_to_cut;

    if (subcell_type == cell::type::triangle && cut_point_tokens.size() == 2)
        return CellCertTag::ready_to_cut;

    if (subcell_type == cell::type::tetrahedron
        && (cut_point_tokens.size() == 3 || cut_point_tokens.size() == 4))
    {
        return CellCertTag::ready_to_cut;
    }

    if (!cut_point_tokens.empty())
        return CellCertTag::cut;

    if (has_zero_edge)
        return CellCertTag::not_classified;

    return CellCertTag::not_classified;
}

/// Check whether any edge incident to face_id has a root (one_root, multiple_roots, or zero).
template <std::floating_point T>
bool face_has_any_edge_root(const AdaptCell<T>& adapt_cell,
                            int level_set_id, int face_id)
{
    auto face_verts = adapt_cell.entity_to_vertex[2][static_cast<std::int32_t>(face_id)];
    const int n_edges = adapt_cell.n_entities(1);

    for (int e = 0; e < n_edges; ++e)
    {
        auto ev = adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(e)];
        bool v0_in = false, v1_in = false;
        for (auto fv : face_verts)
        {
            if (fv == ev[0]) v0_in = true;
            if (fv == ev[1]) v1_in = true;
        }
        if (!v0_in || !v1_in)
            continue;

        const EdgeRootTag etag = adapt_cell.get_edge_root_tag(level_set_id, e);
        if (etag == EdgeRootTag::one_root
            || etag == EdgeRootTag::multiple_roots
            || etag == EdgeRootTag::zero)
        {
            return true;
        }
    }
    return false;
}

} // anonymous namespace

// =====================================================================
// restrict_subcell_bernstein_exact
// =====================================================================

template <std::floating_point T>
void restrict_subcell_bernstein_exact(cell::type parent_cell_type,
                                      int degree,
                                      std::span<const T> parent_coeffs,
                                      cell::type subcell_type,
                                      std::span<const T> subcell_vertex_coords,
                                      std::vector<T>& subcell_coeffs)
{
    // Strategy: generate reference points in the subcell, map them to the
    // parent reference frame using the subcell vertex coordinates as an
    // affine map, evaluate the parent Bernstein at each mapped point,
    // then convert from nodal values to Bernstein on the subcell.

    const int subcell_tdim = cell::get_tdim(subcell_type);
    const int parent_tdim = cell::get_tdim(parent_cell_type);
    const int n_dofs = bernstein::num_polynomials(subcell_type, degree);

    // Number of subcell vertices.
    int n_subcell_verts = 0;
    switch (subcell_type)
    {
    case cell::type::interval:      n_subcell_verts = 2; break;
    case cell::type::triangle:      n_subcell_verts = 3; break;
    case cell::type::tetrahedron:   n_subcell_verts = 4; break;
    case cell::type::quadrilateral: n_subcell_verts = 4; break;
    case cell::type::hexahedron:    n_subcell_verts = 8; break;
    default:
        throw std::runtime_error(
            "restrict_subcell_bernstein_exact: unsupported subcell type");
    }

    // Generate equispaced reference points on the subcell.
    // For simplices: use lattice of order = degree.
    // For tensor-product: use tensor-product grid of order = degree.
    // These serve as Lagrange interpolation nodes.
    std::vector<T> ref_points_flat(static_cast<std::size_t>(n_dofs * subcell_tdim));

    if (bernstein::is_simplex(subcell_type))
    {
        // Equispaced lattice on the reference simplex.
        int idx = 0;
        if (subcell_tdim == 1)
        {
            for (int i = 0; i <= degree; ++i, ++idx)
                ref_points_flat[static_cast<std::size_t>(idx)] =
                    (degree > 0) ? T(i) / T(degree) : T(0);
        }
        else if (subcell_tdim == 2)
        {
            for (int j = 0; j <= degree; ++j)
                for (int i = 0; i <= degree - j; ++i, ++idx)
                {
                    ref_points_flat[static_cast<std::size_t>(idx * 2)] =
                        (degree > 0) ? T(i) / T(degree) : T(0);
                    ref_points_flat[static_cast<std::size_t>(idx * 2 + 1)] =
                        (degree > 0) ? T(j) / T(degree) : T(0);
                }
        }
        else if (subcell_tdim == 3)
        {
            for (int k = 0; k <= degree; ++k)
                for (int j = 0; j <= degree - k; ++j)
                    for (int i = 0; i <= degree - k - j; ++i, ++idx)
                    {
                        ref_points_flat[static_cast<std::size_t>(idx * 3)] =
                            (degree > 0) ? T(i) / T(degree) : T(0);
                        ref_points_flat[static_cast<std::size_t>(idx * 3 + 1)] =
                            (degree > 0) ? T(j) / T(degree) : T(0);
                        ref_points_flat[static_cast<std::size_t>(idx * 3 + 2)] =
                            (degree > 0) ? T(k) / T(degree) : T(0);
                    }
        }
    }
    else // tensor-product
    {
        int idx = 0;
        if (subcell_tdim == 2)
        {
            for (int iy = 0; iy <= degree; ++iy)
                for (int ix = 0; ix <= degree; ++ix, ++idx)
                {
                    ref_points_flat[static_cast<std::size_t>(idx * 2)] =
                        (degree > 0) ? T(ix) / T(degree) : T(0);
                    ref_points_flat[static_cast<std::size_t>(idx * 2 + 1)] =
                        (degree > 0) ? T(iy) / T(degree) : T(0);
                }
        }
        else if (subcell_tdim == 3)
        {
            for (int iz = 0; iz <= degree; ++iz)
                for (int iy = 0; iy <= degree; ++iy)
                    for (int ix = 0; ix <= degree; ++ix, ++idx)
                    {
                        ref_points_flat[static_cast<std::size_t>(idx * 3)] =
                            (degree > 0) ? T(ix) / T(degree) : T(0);
                        ref_points_flat[static_cast<std::size_t>(idx * 3 + 1)] =
                            (degree > 0) ? T(iy) / T(degree) : T(0);
                        ref_points_flat[static_cast<std::size_t>(idx * 3 + 2)] =
                            (degree > 0) ? T(iz) / T(degree) : T(0);
                    }
        }
    }

    // Map each subcell reference point to parent reference coordinates,
    // then evaluate the parent polynomial.
    std::vector<T> nodal_values(static_cast<std::size_t>(n_dofs));
    std::vector<T> xi_parent(static_cast<std::size_t>(parent_tdim));

    for (int d = 0; d < n_dofs; ++d)
    {
        // Map subcell reference point → parent reference.
        // For simplex subcells: xi_parent = (1 - sum(xi_sub)) * v0 + xi_sub[0]*v1 + xi_sub[1]*v2 + ...
        // For tensor-product subcells: bilinear/trilinear map.
        std::fill(xi_parent.begin(), xi_parent.end(), T(0));

        if (bernstein::is_simplex(subcell_type))
        {
            // Barycentric coordinates: lambda_0 = 1 - sum(xi), lambda_i = xi[i-1]
            T lambda0 = T(1);
            for (int s = 0; s < subcell_tdim; ++s)
                lambda0 -= ref_points_flat[static_cast<std::size_t>(d * subcell_tdim + s)];

            for (int c = 0; c < parent_tdim; ++c)
                xi_parent[static_cast<std::size_t>(c)] =
                    lambda0 * subcell_vertex_coords[static_cast<std::size_t>(0 * parent_tdim + c)];

            for (int v = 1; v < n_subcell_verts; ++v)
            {
                T lambda = ref_points_flat[static_cast<std::size_t>(d * subcell_tdim + (v - 1))];
                for (int c = 0; c < parent_tdim; ++c)
                    xi_parent[static_cast<std::size_t>(c)] +=
                        lambda * subcell_vertex_coords[static_cast<std::size_t>(v * parent_tdim + c)];
            }
        }
        else
        {
            // Tensor-product: bilinear/trilinear map.
            if (subcell_tdim == 2)
            {
                T u = ref_points_flat[static_cast<std::size_t>(d * 2)];
                T v = ref_points_flat[static_cast<std::size_t>(d * 2 + 1)];
                // Bilinear: (1-u)(1-v)*v0 + u(1-v)*v1 + uv*v2 + (1-u)v*v3
                T w00 = (T(1) - u) * (T(1) - v);
                T w10 = u * (T(1) - v);
                T w11 = u * v;
                T w01 = (T(1) - u) * v;
                for (int c = 0; c < parent_tdim; ++c)
                    xi_parent[static_cast<std::size_t>(c)] =
                        w00 * subcell_vertex_coords[static_cast<std::size_t>(0 * parent_tdim + c)]
                      + w10 * subcell_vertex_coords[static_cast<std::size_t>(1 * parent_tdim + c)]
                      + w11 * subcell_vertex_coords[static_cast<std::size_t>(2 * parent_tdim + c)]
                      + w01 * subcell_vertex_coords[static_cast<std::size_t>(3 * parent_tdim + c)];
            }
            else if (subcell_tdim == 3)
            {
                T u = ref_points_flat[static_cast<std::size_t>(d * 3)];
                T v = ref_points_flat[static_cast<std::size_t>(d * 3 + 1)];
                T w = ref_points_flat[static_cast<std::size_t>(d * 3 + 2)];
                T weights[8] = {
                    (T(1)-u)*(T(1)-v)*(T(1)-w),
                     u      *(T(1)-v)*(T(1)-w),
                     u      * v      *(T(1)-w),
                    (T(1)-u)* v      *(T(1)-w),
                    (T(1)-u)*(T(1)-v)* w,
                     u      *(T(1)-v)* w,
                     u      * v      * w,
                    (T(1)-u)* v      * w
                };
                for (int c = 0; c < parent_tdim; ++c)
                {
                    xi_parent[static_cast<std::size_t>(c)] = T(0);
                    for (int vv = 0; vv < 8; ++vv)
                        xi_parent[static_cast<std::size_t>(c)] +=
                            weights[vv] * subcell_vertex_coords[
                                static_cast<std::size_t>(vv * parent_tdim + c)];
                }
            }
        }

        nodal_values[static_cast<std::size_t>(d)] =
            bernstein::evaluate(parent_cell_type, degree, parent_coeffs,
                                std::span<const T>(xi_parent));
    }

    // Convert sampled nodal values to Bernstein on the subcell.
    bernstein::lagrange_to_bernstein(
        subcell_type, degree,
        std::span<const T>(ref_points_flat),
        std::span<const T>(nodal_values),
        subcell_coeffs);
}

// =====================================================================
// classify_leaf_cell
// =====================================================================

template <std::floating_point T, std::integral I>
CellCertTag classify_leaf_cell(const AdaptCell<T>& adapt_cell,
                               const LevelSetCell<T, I>& ls_cell,
                               int level_set_id,
                               int cell_id,
                               T zero_tol, T sign_tol)
{
    // A. Check the incident-edge topology first.
    {
        const CellCertTag edge_topology_tag =
            classify_ready_to_cut_topology(adapt_cell, level_set_id, cell_id);
        if (edge_topology_tag != CellCertTag::not_classified)
            return edge_topology_tag;
    }

    // B. No edge roots: use exact subcell Bernstein restriction.
    const int tdim = adapt_cell.tdim;
    auto cell_verts = adapt_cell.entity_to_vertex[tdim][static_cast<std::int32_t>(cell_id)];
    // Get the subcell vertex coordinates in parent reference frame.
    const int parent_tdim = adapt_cell.tdim;
    const int n_cell_verts = static_cast<int>(cell_verts.size());
    std::vector<T> subcell_vertex_coords(
        static_cast<std::size_t>(n_cell_verts * parent_tdim));

    for (int v = 0; v < n_cell_verts; ++v)
    {
        int gv = cell_verts[static_cast<std::size_t>(v)];
        for (int d = 0; d < parent_tdim; ++d)
            subcell_vertex_coords[static_cast<std::size_t>(v * parent_tdim + d)] =
                adapt_cell.vertex_coords[static_cast<std::size_t>(gv * parent_tdim + d)];
    }

    cell::type subcell_type = adapt_cell.entity_types[tdim][static_cast<std::size_t>(cell_id)];

    std::vector<T> subcell_coeffs;
    restrict_subcell_bernstein_exact(
        ls_cell.cell_type, ls_cell.bernstein_order,
        std::span<const T>(ls_cell.bernstein_coeffs),
        subcell_type,
        std::span<const T>(subcell_vertex_coords),
        subcell_coeffs);

    // Sign-hull classification.
    std::span<const T> sc(subcell_coeffs);

    if (bernstein_all_zero(sc, zero_tol))
    {
        // All Bernstein coefficients are zero within tolerance.
        // Evaluate the level set at the cell centroid to disambiguate
        // between a truly zero cell and a cell that is inside/outside
        // but whose vertices all happen to lie on the zero set.
        std::vector<T> xi_centroid(static_cast<std::size_t>(parent_tdim), T(0));
        for (int v = 0; v < n_cell_verts; ++v)
        {
            int gv = cell_verts[static_cast<std::size_t>(v)];
            for (int d = 0; d < parent_tdim; ++d)
                xi_centroid[static_cast<std::size_t>(d)] +=
                    adapt_cell.vertex_coords[static_cast<std::size_t>(gv * parent_tdim + d)];
        }
        for (int d = 0; d < parent_tdim; ++d)
            xi_centroid[static_cast<std::size_t>(d)] /= T(n_cell_verts);

        const T centroid_val = ls_cell.value(
            std::span<const T>(xi_centroid.data(),
                               static_cast<std::size_t>(parent_tdim)));

        if (std::fabs(centroid_val) <= zero_tol)
            return CellCertTag::zero;
        return (centroid_val > T(0)) ? CellCertTag::positive : CellCertTag::negative;
    }
    if (bernstein_all_positive(sc, sign_tol))
        return CellCertTag::positive;
    if (bernstein_all_negative(sc, sign_tol))
        return CellCertTag::negative;

    // Mixed Bernstein signs — apply additional filters before declaring ambiguous.

    // B4a. MONOTONICITY FILTER (cheap, all dimensions):
    //   If phi has a sign-definite partial derivative in some direction on
    //   the subcell AND no incident edge has a root, then phi cannot have
    //   an interior zero-crossing. Use vertex evaluation to determine sign.
    const bool no_edge_root = !cell_has_any_edge_root(adapt_cell, level_set_id, cell_id);

    if (no_edge_root
        && has_monotone_direction(subcell_type, ls_cell.bernstein_order, sc, sign_tol))
    {
        // Evaluate phi at the first vertex of the subcell.
        std::span<const T> xi_v0(
            adapt_cell.vertex_coords.data()
                + static_cast<std::size_t>(cell_verts[0]) * static_cast<std::size_t>(parent_tdim),
            static_cast<std::size_t>(parent_tdim));
        const T val = bernstein::evaluate(
            ls_cell.cell_type, ls_cell.bernstein_order,
            std::span<const T>(ls_cell.bernstein_coeffs), xi_v0);
        return (val >= T(0)) ? CellCertTag::positive : CellCertTag::negative;
    }

    // B4b. FACE CERTIFICATION FILTER (3D only):
    //   For each face of the subcell, restrict the subcell Bernstein to
    //   the face and check if the face can have a root. If ALL faces are
    //   sign-definite and no edge has a root, the cell is sign-definite.
    if (tdim == 3 && no_edge_root && subcell_type == cell::type::tetrahedron)
    {
        const int n_faces = cell::num_faces(subcell_type);
        bool all_faces_sign_definite = true;

        for (int fi = 0; fi < n_faces; ++fi)
        {
            auto face_v = cell::tet_face(fi);
            // Build face vertex coordinates in parent reference frame
            // (using the subcell's local vertex ordering).
            std::vector<T> face_vertex_coords(static_cast<std::size_t>(3 * parent_tdim));
            for (int fv = 0; fv < 3; ++fv)
            {
                const int local_v = face_v[static_cast<std::size_t>(fv)];
                for (int d = 0; d < parent_tdim; ++d)
                    face_vertex_coords[static_cast<std::size_t>(fv * parent_tdim + d)] =
                        subcell_vertex_coords[static_cast<std::size_t>(local_v * parent_tdim + d)];
            }

            // Restrict the parent polynomial to the face.
            std::vector<T> face_coeffs;
            restrict_subcell_bernstein_exact(
                ls_cell.cell_type, ls_cell.bernstein_order,
                std::span<const T>(ls_cell.bernstein_coeffs),
                cell::type::triangle,
                std::span<const T>(face_vertex_coords),
                face_coeffs);

            std::span<const T> fc(face_coeffs);
            if (bernstein_all_positive(fc, sign_tol)
                || bernstein_all_negative(fc, sign_tol)
                || bernstein_all_zero(fc, zero_tol))
            {
                continue; // This face is sign-definite or zero — no root here.
            }

            // Face has mixed Bernstein signs — possible root.
            all_faces_sign_definite = false;
            break;
        }

        if (all_faces_sign_definite)
        {
            // No root can enter the cell through any face; evaluate at a vertex.
            std::span<const T> xi_v0(
                adapt_cell.vertex_coords.data()
                    + static_cast<std::size_t>(cell_verts[0]) * static_cast<std::size_t>(parent_tdim),
                static_cast<std::size_t>(parent_tdim));
            const T val = bernstein::evaluate(
                ls_cell.cell_type, ls_cell.bernstein_order,
                std::span<const T>(ls_cell.bernstein_coeffs), xi_v0);
            return (val >= T(0)) ? CellCertTag::positive : CellCertTag::negative;
        }
    }

    return CellCertTag::ambiguous;
}

// =====================================================================
// classify_leaf_cells
// =====================================================================

template <std::floating_point T, std::integral I>
void classify_leaf_cells(AdaptCell<T>& adapt_cell,
                         const LevelSetCell<T, I>& ls_cell,
                         int level_set_id,
                         T zero_tol, T sign_tol)
{
    const int tdim = adapt_cell.tdim;
    const int n_cells = adapt_cell.n_entities(tdim);

    // Ensure tag storage.
    if (adapt_cell.cell_cert_tag_num_level_sets <= level_set_id)
        adapt_cell.resize_cell_cert_tags(level_set_id + 1);

    for (int c = 0; c < n_cells; ++c)
    {
        if (adapt_cell.get_cell_cert_tag(level_set_id, c)
            != CellCertTag::not_classified)
            continue;

        CellCertTag tag = classify_leaf_cell(
            adapt_cell, ls_cell, level_set_id, c, zero_tol, sign_tol);
        adapt_cell.set_cell_cert_tag(level_set_id, c, tag);
    }
}

// =====================================================================
// classify_leaf_face
// =====================================================================

template <std::floating_point T, std::integral I>
FaceCertTag classify_leaf_face(const AdaptCell<T>& adapt_cell,
                               const LevelSetCell<T, I>& ls_cell,
                               int level_set_id,
                               int face_id,
                               T zero_tol, T sign_tol)
{
    // Only meaningful for 3D cells.
    if (adapt_cell.tdim != 3)
        return FaceCertTag::not_classified;

    // A. Check incident-edge root topology.
    const bool has_root = face_has_any_edge_root(adapt_cell, level_set_id, face_id);
    if (has_root)
        return FaceCertTag::cut;

    // B. Restrict the parent Bernstein to the face.
    auto face_verts = adapt_cell.entity_to_vertex[2][static_cast<std::int32_t>(face_id)];
    const int parent_tdim = adapt_cell.tdim;
    const int n_face_verts = static_cast<int>(face_verts.size());
    std::vector<T> face_vertex_coords(
        static_cast<std::size_t>(n_face_verts * parent_tdim));

    for (int v = 0; v < n_face_verts; ++v)
    {
        int gv = face_verts[static_cast<std::size_t>(v)];
        for (int d = 0; d < parent_tdim; ++d)
            face_vertex_coords[static_cast<std::size_t>(v * parent_tdim + d)] =
                adapt_cell.vertex_coords[static_cast<std::size_t>(gv * parent_tdim + d)];
    }

    cell::type face_type = adapt_cell.entity_types[2][static_cast<std::size_t>(face_id)];

    std::vector<T> face_coeffs;
    restrict_subcell_bernstein_exact(
        ls_cell.cell_type, ls_cell.bernstein_order,
        std::span<const T>(ls_cell.bernstein_coeffs),
        face_type,
        std::span<const T>(face_vertex_coords),
        face_coeffs);

    // Sign-hull classification.
    std::span<const T> fc(face_coeffs);

    if (bernstein_all_zero(fc, zero_tol))
    {
        // All Bernstein coefficients are zero — evaluate centroid.
        std::vector<T> xi_centroid(static_cast<std::size_t>(parent_tdim), T(0));
        for (int v = 0; v < n_face_verts; ++v)
        {
            int gv = face_verts[static_cast<std::size_t>(v)];
            for (int d = 0; d < parent_tdim; ++d)
                xi_centroid[static_cast<std::size_t>(d)] +=
                    adapt_cell.vertex_coords[static_cast<std::size_t>(gv * parent_tdim + d)];
        }
        for (int d = 0; d < parent_tdim; ++d)
            xi_centroid[static_cast<std::size_t>(d)] /= T(n_face_verts);

        const T centroid_val = ls_cell.value(
            std::span<const T>(xi_centroid.data(),
                               static_cast<std::size_t>(parent_tdim)));

        if (std::fabs(centroid_val) <= zero_tol)
            return FaceCertTag::zero;
        return (centroid_val > T(0)) ? FaceCertTag::positive : FaceCertTag::negative;
    }
    if (bernstein_all_positive(fc, sign_tol))
        return FaceCertTag::positive;
    if (bernstein_all_negative(fc, sign_tol))
        return FaceCertTag::negative;

    // Mixed Bernstein signs — apply monotonicity filter.
    if (has_monotone_direction(face_type, ls_cell.bernstein_order, fc, sign_tol))
    {
        // Evaluate phi at the first vertex of the face.
        std::span<const T> xi_v0(
            adapt_cell.vertex_coords.data()
                + static_cast<std::size_t>(face_verts[0]) * static_cast<std::size_t>(parent_tdim),
            static_cast<std::size_t>(parent_tdim));
        const T val = bernstein::evaluate(
            ls_cell.cell_type, ls_cell.bernstein_order,
            std::span<const T>(ls_cell.bernstein_coeffs), xi_v0);
        return (val >= T(0)) ? FaceCertTag::positive : FaceCertTag::negative;
    }

    return FaceCertTag::ambiguous;
}

// =====================================================================
// classify_leaf_faces
// =====================================================================

template <std::floating_point T, std::integral I>
void classify_leaf_faces(AdaptCell<T>& adapt_cell,
                         const LevelSetCell<T, I>& ls_cell,
                         int level_set_id,
                         T zero_tol, T sign_tol)
{
    if (adapt_cell.tdim != 3)
        return;

    const int n_faces = adapt_cell.n_entities(2);
    if (n_faces == 0)
        return;

    // Ensure tag storage.
    if (adapt_cell.face_cert_tag_num_level_sets <= level_set_id)
        adapt_cell.resize_face_cert_tags(level_set_id + 1);

    for (int f = 0; f < n_faces; ++f)
    {
        if (adapt_cell.get_face_cert_tag(level_set_id, f)
            != FaceCertTag::not_classified)
            continue;

        FaceCertTag tag = classify_leaf_face(
            adapt_cell, ls_cell, level_set_id, f, zero_tol, sign_tol);
        adapt_cell.set_face_cert_tag(level_set_id, f, tag);
    }
}

// =====================================================================
// fill_all_vertex_signs_from_level_set
// =====================================================================

template <std::floating_point T, std::integral I>
void fill_all_vertex_signs_from_level_set(AdaptCell<T>& adapt_cell,
                                          const LevelSetCell<T, I>& ls_cell,
                                          int level_set_id,
                                          T zero_tol)
{
    const int n_vertices = adapt_cell.n_vertices();
    const int tdim = adapt_cell.tdim;
    for (int v = 0; v < n_vertices; ++v)
    {
        const std::uint64_t bit = std::uint64_t(1) << level_set_id;
        if ((adapt_cell.zero_mask_per_vertex[static_cast<std::size_t>(v)] & bit) != 0)
            continue;

        std::span<const T> xi(
            adapt_cell.vertex_coords.data()
                + static_cast<std::size_t>(v) * static_cast<std::size_t>(tdim),
            static_cast<std::size_t>(tdim));
        set_vertex_sign_for_level_set(
            adapt_cell, v, level_set_id, ls_cell.value(xi), zero_tol);
    }
}

// =====================================================================
// process_ready_to_cut_cells
// =====================================================================

template <std::floating_point T, std::integral I>
void process_ready_to_cut_cells(AdaptCell<T>& adapt_cell,
                                const LevelSetCell<T, I>& ls_cell,
                                int level_set_id,
                                T zero_tol,
                                T sign_tol,
                                int edge_max_depth,
                                bool triangulate_cut_parts)
{
    process_ready_to_cut_cells(
        adapt_cell, ls_cell, level_set_id, zero_tol, sign_tol, edge_max_depth,
        cell::triangulation_strategy_from_bool(triangulate_cut_parts));
}

template <std::floating_point T, std::integral I>
void process_ready_to_cut_cells(AdaptCell<T>& adapt_cell,
                                const LevelSetCell<T, I>& ls_cell,
                                int level_set_id,
                                T zero_tol,
                                T sign_tol,
                                int edge_max_depth,
                                cell::TriangulationStrategy triangulation_strategy)
{
    const int tdim = adapt_cell.tdim;
    const int n_cells = adapt_cell.n_entities(tdim);

    bool has_ready = false;
    for (int c = 0; c < n_cells; ++c)
    {
        if (adapt_cell.get_cell_cert_tag(level_set_id, c) == CellCertTag::ready_to_cut)
        {
            has_ready = true;
            break;
        }
    }
    if (!has_ready)
        return;

    const std::vector<cell::type> old_types = adapt_cell.entity_types[tdim];
    const EntityAdjacency old_cells = adapt_cell.entity_to_vertex[tdim];
    const auto edge_lookup = build_leaf_edge_lookup(adapt_cell);

    EntityAdjacency new_cells;
    new_cells.offsets.push_back(0);
    std::vector<cell::type> new_types;
    std::vector<int> old_cell_ids_for_new_cells;
    std::vector<int> source_cell_ids_for_new_cells;
    std::vector<CellRefinementReason> refinement_reasons_for_new_cells;
    std::vector<CellCertTag> explicit_current_ls_tags;

    for (int c = 0; c < n_cells; ++c)
    {
        auto old_cell_vertices = old_cells[static_cast<std::int32_t>(c)];
        const cell::type leaf_cell_type = old_types[static_cast<std::size_t>(c)];

        if (adapt_cell.get_cell_cert_tag(level_set_id, c) != CellCertTag::ready_to_cut)
        {
            std::vector<int> copy(old_cell_vertices.begin(), old_cell_vertices.end());
            append_top_cell_local(new_types, new_cells, leaf_cell_type, std::span<const int>(copy));
            old_cell_ids_for_new_cells.push_back(c);
            source_cell_ids_for_new_cells.push_back(c);
            refinement_reasons_for_new_cells.push_back(CellRefinementReason::none);
            explicit_current_ls_tags.push_back(CellCertTag::not_classified);
            continue;
        }

        if (leaf_cell_type != cell::type::interval
            && leaf_cell_type != cell::type::triangle
            && leaf_cell_type != cell::type::tetrahedron)
        {
            throw std::runtime_error(
                "process_ready_to_cut_cells: ready_to_cut only implemented for interval, triangle, and tetrahedron leaves");
        }

        std::vector<T> vertex_coords(
            static_cast<std::size_t>(old_cell_vertices.size() * tdim), T(0));
        for (std::size_t j = 0; j < old_cell_vertices.size(); ++j)
        {
            const int gv = old_cell_vertices[j];
            for (int d = 0; d < tdim; ++d)
            {
                vertex_coords[static_cast<std::size_t>(j * tdim + d)] =
                    adapt_cell.vertex_coords[static_cast<std::size_t>(gv * tdim + d)];
            }
        }
        const std::vector<T> ls_values =
            gather_leaf_cell_vertex_level_set_values(adapt_cell, ls_cell, c);

        const std::uint64_t current_level_set_bit = std::uint64_t(1) << level_set_id;
        std::vector<bool> current_level_set_zero(old_cell_vertices.size(), false);
        bool has_strict_negative = false;
        bool has_strict_positive = false;
        bool has_zero = false;
        for (std::size_t j = 0; j < ls_values.size(); ++j)
        {
            const T value = ls_values[j];
            const int gv = old_cell_vertices[j];
            current_level_set_zero[j] =
                (adapt_cell.zero_mask_per_vertex[static_cast<std::size_t>(gv)]
                 & current_level_set_bit) != 0;
            has_strict_negative = has_strict_negative || value < -zero_tol;
            has_strict_positive = has_strict_positive || value > zero_tol;
            has_zero = has_zero || current_level_set_zero[j];
        }
        if (!(has_strict_negative && has_strict_positive))
        {
            std::vector<int> copy(old_cell_vertices.begin(), old_cell_vertices.end());
            append_top_cell_local(new_types, new_cells, leaf_cell_type, std::span<const int>(copy));
            old_cell_ids_for_new_cells.push_back(c);
            CellCertTag copied_tag = CellCertTag::zero;
            if (has_strict_negative)
                copied_tag = CellCertTag::negative;
            else if (has_strict_positive)
                copied_tag = CellCertTag::positive;
            else if (!has_zero)
                copied_tag = CellCertTag::not_classified;
            explicit_current_ls_tags.push_back(copied_tag);
            continue;
        }

        std::array<int, 6> old_edge_ids_by_local_edge;
        old_edge_ids_by_local_edge.fill(-1);
        const auto ledges = cell::edges(leaf_cell_type);
        for (std::size_t le = 0; le < ledges.size(); ++le)
        {
            const int a = old_cell_vertices[static_cast<std::size_t>(ledges[le][0])];
            const int b = old_cell_vertices[static_cast<std::size_t>(ledges[le][1])];
            const std::pair<int, int> key = {std::min(a, b), std::max(a, b)};
            auto edge_it = edge_lookup.find(key);
            if (edge_it == edge_lookup.end())
            {
                throw std::runtime_error(
                    "process_ready_to_cut_cells: missing leaf edge for ready-to-cut cell "
                    + std::to_string(c) + " type="
                    + cell_type_to_str(leaf_cell_type) + " local_edge="
                    + std::to_string(le) + " key=("
                    + std::to_string(key.first) + ","
                    + std::to_string(key.second) + ")");
            }
            old_edge_ids_by_local_edge[le] = edge_it->second;
        }

        std::map<int, int> token_to_vertex;
        for (std::size_t j = 0; j < old_cell_vertices.size(); ++j)
            token_to_vertex[100 + static_cast<int>(j)] = old_cell_vertices[j];
        for (std::size_t le = 0; le < ledges.size(); ++le)
        {
            const int old_edge_id = old_edge_ids_by_local_edge[le];
            if (old_edge_id < 0)
                continue;
            if (adapt_cell.get_edge_root_tag(level_set_id, old_edge_id)
                != EdgeRootTag::one_root)
            {
                continue;
            }

            const int root_vertex_id = ensure_one_root_vertex_on_edge(
                adapt_cell, ls_cell, level_set_id, old_edge_id,
                zero_tol, sign_tol, edge_max_depth);
            if (root_vertex_id >= 0)
                token_to_vertex[static_cast<int>(le)] = root_vertex_id;
        }

        if (leaf_cell_type == cell::type::interval)
        {
            cell::CutCell<T> negative_part;
            cell::CutCell<T> positive_part;
            cell::interval::cut(
                std::span<const T>(vertex_coords),
                tdim,
                std::span<const T>(ls_values),
                "phi<0",
                negative_part);
            cell::interval::cut(
                std::span<const T>(vertex_coords),
                tdim,
                std::span<const T>(ls_values),
                "phi>0",
                positive_part);

            append_ready_cut_part_cells(
                adapt_cell, ls_cell, level_set_id, c, negative_part, leaf_cell_type,
                old_cell_vertices, std::span<const int>(old_edge_ids_by_local_edge.data(), 1),
                new_types, new_cells, old_cell_ids_for_new_cells,
                source_cell_ids_for_new_cells, refinement_reasons_for_new_cells,
                explicit_current_ls_tags,
                token_to_vertex, zero_tol, CellCertTag::negative);
            append_ready_cut_part_cells(
                adapt_cell, ls_cell, level_set_id, c, positive_part, leaf_cell_type,
                old_cell_vertices, std::span<const int>(old_edge_ids_by_local_edge.data(), 1),
                new_types, new_cells, old_cell_ids_for_new_cells,
                source_cell_ids_for_new_cells, refinement_reasons_for_new_cells,
                explicit_current_ls_tags,
                token_to_vertex, zero_tol, CellCertTag::positive);
        }
        else if (leaf_cell_type == cell::type::triangle)
        {
            if (has_zero)
            {
                auto append_zero_vertex_triangle_side =
                    [&](bool negative_side, CellCertTag side_tag)
                {
                    std::vector<int> tokens;
                    auto append_token = [&](int token)
                    {
                        if (tokens.empty() || tokens.back() != token)
                            tokens.push_back(token);
                    };

                    for (int lv = 0; lv < 3; ++lv)
                    {
                        const int next = (lv + 1) % 3;
                        const T vi = ls_values[static_cast<std::size_t>(lv)];
                        const T vj = ls_values[static_cast<std::size_t>(next)];
                        const bool zero_i =
                            current_level_set_zero[static_cast<std::size_t>(lv)];
                        const bool inside_i = negative_side
                                                  ? (zero_i || vi < -zero_tol)
                                                  : (zero_i || vi > zero_tol);
                        if (inside_i)
                            append_token(100 + lv);

                        const bool strict_cross =
                            (vi < -zero_tol && vj > zero_tol)
                            || (vi > zero_tol && vj < -zero_tol);
                        if (strict_cross)
                        {
                            const int local_edge =
                                basix_edge_id_for_vertices(leaf_cell_type, lv, next);
                            append_token(local_edge);
                        }
                    }

                    if (tokens.size() > 1 && tokens.front() == tokens.back())
                        tokens.pop_back();
                    if (tokens.size() < 3)
                        return;
                    if (tokens.size() != 3)
                    {
                        throw std::runtime_error(
                            "process_ready_to_cut_cells: zero-vertex triangle clipping "
                            "expected a triangle");
                    }

                    std::vector<int> mapped(tokens.size(), -1);
                    for (std::size_t j = 0; j < tokens.size(); ++j)
                    {
                        auto token_it = token_to_vertex.find(tokens[j]);
                        if (token_it == token_to_vertex.end())
                        {
                            throw std::runtime_error(
                                "process_ready_to_cut_cells: missing zero-vertex "
                                "triangle token " + std::to_string(tokens[j]));
                        }
                        mapped[j] = token_it->second;
                    }

                    append_top_cell_local(
                        new_types, new_cells, cell::type::triangle,
                        std::span<const int>(mapped));
                    old_cell_ids_for_new_cells.push_back(-1);
                    source_cell_ids_for_new_cells.push_back(c);
                    refinement_reasons_for_new_cells.push_back(
                        CellRefinementReason::cut_level_set);
                    explicit_current_ls_tags.push_back(side_tag);
                };

                append_zero_vertex_triangle_side(true, CellCertTag::negative);
                append_zero_vertex_triangle_side(false, CellCertTag::positive);
                continue;
            }

            cell::CutCell<T> negative_part;
            cell::CutCell<T> positive_part;
            cell::triangle::cut(
                std::span<const T>(vertex_coords),
                tdim,
                std::span<const T>(ls_values),
                "phi<0",
                negative_part,
                triangulation_strategy);
            cell::triangle::cut(
                std::span<const T>(vertex_coords),
                tdim,
                std::span<const T>(ls_values),
                "phi>0",
                positive_part,
                triangulation_strategy);

            append_ready_cut_part_cells(
                adapt_cell, ls_cell, level_set_id, c, negative_part, leaf_cell_type,
                old_cell_vertices, std::span<const int>(old_edge_ids_by_local_edge.data(), 3),
                new_types, new_cells, old_cell_ids_for_new_cells,
                source_cell_ids_for_new_cells, refinement_reasons_for_new_cells,
                explicit_current_ls_tags,
                token_to_vertex, zero_tol, CellCertTag::negative);
            append_ready_cut_part_cells(
                adapt_cell, ls_cell, level_set_id, c, positive_part, leaf_cell_type,
                old_cell_vertices, std::span<const int>(old_edge_ids_by_local_edge.data(), 3),
                new_types, new_cells, old_cell_ids_for_new_cells,
                source_cell_ids_for_new_cells, refinement_reasons_for_new_cells,
                explicit_current_ls_tags,
                token_to_vertex, zero_tol, CellCertTag::positive);
        }
        else
        {
            if (has_zero)
            {
                auto append_zero_vertex_tetra_side =
                    [&](bool negative_side, CellCertTag side_tag)
                {
                    std::vector<int> zero_vertices;
                    std::vector<int> inside_vertices;
                    std::vector<int> outside_vertices;

                    for (int lv = 0; lv < 4; ++lv)
                    {
                        const T value = ls_values[static_cast<std::size_t>(lv)];
                        if (current_level_set_zero[static_cast<std::size_t>(lv)])
                        {
                            zero_vertices.push_back(lv);
                            continue;
                        }

                        const bool inside =
                            negative_side ? (value < -zero_tol) : (value > zero_tol);
                        if (inside)
                            inside_vertices.push_back(lv);
                        else
                            outside_vertices.push_back(lv);
                    }

                    auto edge_token = [&](int a, int b)
                    {
                        const int token =
                            basix_edge_id_for_vertices(leaf_cell_type, a, b);
                        if (token < 0)
                        {
                            throw std::runtime_error(
                                "process_ready_to_cut_cells: missing tetra edge token");
                        }
                        return token;
                    };

                    auto append_tet_tokens = [&](const std::array<int, 4>& tokens)
                    {
                        std::array<int, 4> mapped;
                        for (std::size_t j = 0; j < tokens.size(); ++j)
                        {
                            auto token_it = token_to_vertex.find(tokens[j]);
                            if (token_it == token_to_vertex.end())
                            {
                                throw std::runtime_error(
                                    "process_ready_to_cut_cells: missing zero-vertex "
                                    "tetra token " + std::to_string(tokens[j]));
                            }
                            mapped[j] = token_it->second;
                        }

                        std::set<int> unique(mapped.begin(), mapped.end());
                        if (unique.size() != mapped.size())
                            return;

                        append_top_cell_local(
                            new_types, new_cells, cell::type::tetrahedron,
                            std::span<const int>(mapped));
                        old_cell_ids_for_new_cells.push_back(-1);
                        source_cell_ids_for_new_cells.push_back(c);
                        refinement_reasons_for_new_cells.push_back(
                            CellRefinementReason::cut_level_set);
                        explicit_current_ls_tags.push_back(side_tag);
                    };

                    if (zero_vertices.size() == 1)
                    {
                        const int z = zero_vertices[0];
                        if (inside_vertices.size() == 1
                            && outside_vertices.size() == 2)
                        {
                            const int a = inside_vertices[0];
                            append_tet_tokens({
                                100 + z,
                                100 + a,
                                edge_token(a, outside_vertices[0]),
                                edge_token(a, outside_vertices[1]),
                            });
                            return;
                        }

                        if (inside_vertices.size() == 2
                            && outside_vertices.size() == 1)
                        {
                            const int a = inside_vertices[0];
                            const int b = inside_vertices[1];
                            const int o = outside_vertices[0];
                            const int ra = edge_token(a, o);
                            const int rb = edge_token(b, o);
                            append_tet_tokens({100 + z, 100 + a, 100 + b, rb});
                            append_tet_tokens({100 + z, 100 + a, rb, ra});
                            return;
                        }
                    }

                    if (zero_vertices.size() == 2
                        && inside_vertices.size() == 1
                        && outside_vertices.size() == 1)
                    {
                        append_tet_tokens({
                            100 + zero_vertices[0],
                            100 + zero_vertices[1],
                            100 + inside_vertices[0],
                            edge_token(inside_vertices[0], outside_vertices[0]),
                        });
                        return;
                    }

                    throw std::runtime_error(
                        "process_ready_to_cut_cells: unsupported zero-vertex "
                        "tetrahedron clipping case");
                };

                append_zero_vertex_tetra_side(true, CellCertTag::negative);
                append_zero_vertex_tetra_side(false, CellCertTag::positive);
                continue;
            }

            cell::CutCell<T> negative_part;
            cell::CutCell<T> positive_part;
            cell::tetrahedron::cut(
                std::span<const T>(vertex_coords), tdim,
                std::span<const T>(ls_values), "phi<0",
                negative_part, triangulation_strategy);
            cell::tetrahedron::cut(
                std::span<const T>(vertex_coords), tdim,
                std::span<const T>(ls_values), "phi>0",
                positive_part, triangulation_strategy);

            append_ready_cut_part_cells(
                adapt_cell, ls_cell, level_set_id, c, negative_part, leaf_cell_type,
                old_cell_vertices, std::span<const int>(old_edge_ids_by_local_edge.data(), 6),
                new_types, new_cells, old_cell_ids_for_new_cells,
                source_cell_ids_for_new_cells, refinement_reasons_for_new_cells,
                explicit_current_ls_tags,
                token_to_vertex, zero_tol, CellCertTag::negative);
            append_ready_cut_part_cells(
                adapt_cell, ls_cell, level_set_id, c, positive_part, leaf_cell_type,
                old_cell_vertices, std::span<const int>(old_edge_ids_by_local_edge.data(), 6),
                new_types, new_cells, old_cell_ids_for_new_cells,
                source_cell_ids_for_new_cells, refinement_reasons_for_new_cells,
                explicit_current_ls_tags,
                token_to_vertex, zero_tol, CellCertTag::positive);
        }
    }

    apply_topology_update_preserve_certification(
        adapt_cell, std::move(new_types), std::move(new_cells),
        std::span<const int>(old_cell_ids_for_new_cells),
        std::span<const int>(source_cell_ids_for_new_cells),
        std::span<const CellRefinementReason>(refinement_reasons_for_new_cells));

    for (int c = 0; c < adapt_cell.n_entities(tdim); ++c)
    {
        if (c < static_cast<int>(refinement_reasons_for_new_cells.size())
            && refinement_reasons_for_new_cells[static_cast<std::size_t>(c)]
                   == CellRefinementReason::cut_level_set
            && c < static_cast<int>(adapt_cell.entity_source_level_set[tdim].size()))
        {
            adapt_cell.entity_source_level_set[tdim][static_cast<std::size_t>(c)] =
                static_cast<std::int32_t>(level_set_id);
        }
    }

    const int new_num_cells = adapt_cell.n_entities(tdim);
    for (int c = 0; c < new_num_cells; ++c)
    {
        if (explicit_current_ls_tags[static_cast<std::size_t>(c)] == CellCertTag::not_classified)
            continue;
        adapt_cell.set_cell_cert_tag(
            level_set_id, c, explicit_current_ls_tags[static_cast<std::size_t>(c)]);
    }

}

template <std::floating_point T, std::integral I>
bool refine_ready_cell_on_largest_midpoint_value(
    AdaptCell<T>& adapt_cell,
    const LevelSetCell<T, I>& ls_cell,
    int level_set_id,
    int cell_id)
{
    const int tdim = adapt_cell.tdim;
    if (cell_id < 0 || cell_id >= adapt_cell.n_entities(tdim))
        return false;

    const auto edge_lookup = build_leaf_edge_lookup(adapt_cell);
    auto verts = adapt_cell.entity_to_vertex[tdim][static_cast<std::int32_t>(cell_id)];
    const cell::type ctype = adapt_cell.entity_types[tdim][static_cast<std::size_t>(cell_id)];
    const auto ledges = cell::edges(ctype);

    int selected_edge = -1;
    T selected_value = -std::numeric_limits<T>::infinity();
    for (const auto& le : ledges)
    {
        const int a = verts[static_cast<std::size_t>(le[0])];
        const int b = verts[static_cast<std::size_t>(le[1])];
        const std::pair<int, int> key = {std::min(a, b), std::max(a, b)};
        auto it = edge_lookup.find(key);
        if (it == edge_lookup.end())
            continue;

        std::vector<T> midpoint(static_cast<std::size_t>(tdim), T(0));
        for (int d = 0; d < tdim; ++d)
        {
            midpoint[static_cast<std::size_t>(d)] =
                T(0.5)
                * (adapt_cell.vertex_coords[static_cast<std::size_t>(a * tdim + d)]
                   + adapt_cell.vertex_coords[static_cast<std::size_t>(b * tdim + d)]);
        }
        const T value = ls_cell.value(
            std::span<const T>(midpoint.data(), midpoint.size()));
        if (selected_edge < 0 || value > selected_value)
        {
            selected_edge = it->second;
            selected_value = value;
        }
    }

    if (selected_edge < 0)
        return false;

    const int n_edges = adapt_cell.n_entities(1);
    if (adapt_cell.edge_root_tag_num_level_sets <= level_set_id)
        adapt_cell.resize_edge_root_tags(level_set_id + 1);
    if (adapt_cell.edge_green_split_has_value.size()
        < static_cast<std::size_t>((level_set_id + 1) * n_edges))
    {
        adapt_cell.resize_green_split_data(level_set_id + 1);
    }

    const auto idx =
        static_cast<std::size_t>(level_set_id * n_edges + selected_edge);
    adapt_cell.set_edge_root_tag(
        level_set_id, selected_edge, EdgeRootTag::multiple_roots);
    adapt_cell.edge_green_split_param[idx] = T(0.5);
    adapt_cell.edge_green_split_has_value[idx] = 1;

    return refine_green_on_multiple_root_edges(adapt_cell, level_set_id);
}

// =====================================================================
// certify_and_refine
// =====================================================================

template <std::floating_point T, std::integral I>
void certify_and_refine(AdaptCell<T>& adapt_cell,
                        const LevelSetCell<T, I>& ls_cell,
                        int level_set_id,
                        int max_iterations,
                        T zero_tol, T sign_tol,
                        int edge_max_depth)
{
    for (int iter = 0; iter < max_iterations; ++iter)
    {
        // 1. Classify edges.
        classify_new_edges(adapt_cell, ls_cell, level_set_id,
                           zero_tol, sign_tol, edge_max_depth);

        // 1b. Build faces and classify them (3D only).
        if (adapt_cell.tdim == 3)
        {
            build_faces(adapt_cell);
            classify_leaf_faces(adapt_cell, ls_cell, level_set_id,
                                zero_tol, sign_tol);
        }

        // 2. Classify cells.
        classify_leaf_cells(adapt_cell, ls_cell, level_set_id,
                            zero_tol, sign_tol);

        // 3. Green refinement: multiple_roots edges.
        bool did_green = false;
        const int n_edges = adapt_cell.n_entities(1);
        for (int e = 0; e < n_edges; ++e)
        {
            if (adapt_cell.get_edge_root_tag(level_set_id, e)
                == EdgeRootTag::multiple_roots)
            {
                did_green = true;
                break;
            }
        }

        if (did_green)
        {
            if (refine_green_on_multiple_root_edges(adapt_cell, level_set_id))
            {
                fill_all_vertex_signs_from_level_set(
                    adapt_cell, ls_cell, level_set_id, zero_tol);
                continue;
            }
        }

        // 4. Red refinement: ambiguous cells.
        bool did_red = false;
        const int tdim = adapt_cell.tdim;
        const int n_cells = adapt_cell.n_entities(tdim);
        for (int c = 0; c < n_cells; ++c)
        {
            if (adapt_cell.get_cell_cert_tag(level_set_id, c)
                == CellCertTag::ambiguous)
            {
                did_red = true;
                break;
            }
        }

        if (did_red)
        {
            if (refine_red_on_ambiguous_cells(adapt_cell, level_set_id))
            {
                fill_all_vertex_signs_from_level_set(
                    adapt_cell, ls_cell, level_set_id, zero_tol);
                continue;
            }
        }

        // 5. No refinement needed — stop.
        break;
    }
}

template <std::floating_point T, std::integral I>
void certify_refine_and_process_ready_cells(AdaptCell<T>& adapt_cell,
                                            const LevelSetCell<T, I>& ls_cell,
                                            int level_set_id,
                                            int max_iterations,
                                            T zero_tol, T sign_tol,
                                            int edge_max_depth,
                                            bool triangulate_cut_parts)
{
    certify_refine_and_process_ready_cells(
        adapt_cell, ls_cell, level_set_id, max_iterations, zero_tol, sign_tol,
        edge_max_depth,
        cell::triangulation_strategy_from_bool(triangulate_cut_parts));
}

template <std::floating_point T, std::integral I>
void certify_refine_and_process_ready_cells(AdaptCell<T>& adapt_cell,
                                            const LevelSetCell<T, I>& ls_cell,
                                            int level_set_id,
                                            int max_iterations,
                                            T zero_tol, T sign_tol,
                                            int edge_max_depth,
                                            cell::TriangulationStrategy triangulation_strategy)
{
    fill_all_vertex_signs_from_level_set(adapt_cell, ls_cell, level_set_id, zero_tol);
    certify_and_refine(adapt_cell, ls_cell, level_set_id,
                       max_iterations, zero_tol, sign_tol, edge_max_depth);
    fill_all_vertex_signs_from_level_set(adapt_cell, ls_cell, level_set_id, zero_tol);
    process_ready_to_cut_cells(adapt_cell, ls_cell, level_set_id,
                               zero_tol, sign_tol, edge_max_depth,
                               triangulation_strategy);
    fill_all_vertex_signs_from_level_set(adapt_cell, ls_cell, level_set_id, zero_tol);
    build_edges(adapt_cell);
    if (adapt_cell.tdim == 3)
        build_faces(adapt_cell);
    recompute_active_level_set_masks(adapt_cell, level_set_id + 1);
    rebuild_zero_entity_inventory(adapt_cell);
}

// =====================================================================
// Explicit template instantiations
// =====================================================================

template void restrict_subcell_bernstein_exact(cell::type, int,
                                               std::span<const double>,
                                               cell::type,
                                               std::span<const double>,
                                               std::vector<double>&);
template void restrict_subcell_bernstein_exact(cell::type, int,
                                               std::span<const float>,
                                               cell::type,
                                               std::span<const float>,
                                               std::vector<float>&);

template CellCertTag classify_leaf_cell(const AdaptCell<double>&,
                                        const LevelSetCell<double, int>&,
                                        int, int, double, double);
template CellCertTag classify_leaf_cell(const AdaptCell<float>&,
                                        const LevelSetCell<float, int>&,
                                        int, int, float, float);
template CellCertTag classify_leaf_cell(const AdaptCell<double>&,
                                        const LevelSetCell<double, long>&,
                                        int, int, double, double);

template void classify_leaf_cells(AdaptCell<double>&,
                                  const LevelSetCell<double, int>&,
                                  int, double, double);
template void classify_leaf_cells(AdaptCell<float>&,
                                  const LevelSetCell<float, int>&,
                                  int, float, float);
template void classify_leaf_cells(AdaptCell<double>&,
                                  const LevelSetCell<double, long>&,
                                  int, double, double);

template FaceCertTag classify_leaf_face(const AdaptCell<double>&,
                                        const LevelSetCell<double, int>&,
                                        int, int, double, double);
template FaceCertTag classify_leaf_face(const AdaptCell<float>&,
                                        const LevelSetCell<float, int>&,
                                        int, int, float, float);
template FaceCertTag classify_leaf_face(const AdaptCell<double>&,
                                        const LevelSetCell<double, long>&,
                                        int, int, double, double);

template void classify_leaf_faces(AdaptCell<double>&,
                                  const LevelSetCell<double, int>&,
                                  int, double, double);
template void classify_leaf_faces(AdaptCell<float>&,
                                  const LevelSetCell<float, int>&,
                                  int, float, float);
template void classify_leaf_faces(AdaptCell<double>&,
                                  const LevelSetCell<double, long>&,
                                  int, double, double);

template void fill_all_vertex_signs_from_level_set(AdaptCell<double>&,
                                                   const LevelSetCell<double, int>&,
                                                   int, double);
template void fill_all_vertex_signs_from_level_set(AdaptCell<float>&,
                                                   const LevelSetCell<float, int>&,
                                                   int, float);
template void fill_all_vertex_signs_from_level_set(AdaptCell<double>&,
                                                   const LevelSetCell<double, long>&,
                                                   int, double);
template void fill_all_vertex_signs_from_level_set(AdaptCell<float>&,
                                                   const LevelSetCell<float, long>&,
                                                   int, float);

template void process_ready_to_cut_cells(AdaptCell<double>&,
                                         const LevelSetCell<double, int>&,
                                         int, double, double, int,
                                         cell::TriangulationStrategy);
template void process_ready_to_cut_cells(AdaptCell<float>&,
                                         const LevelSetCell<float, int>&,
                                         int, float, float, int,
                                         cell::TriangulationStrategy);
template void process_ready_to_cut_cells(AdaptCell<double>&,
                                         const LevelSetCell<double, long>&,
                                         int, double, double, int,
                                         cell::TriangulationStrategy);
template void process_ready_to_cut_cells(AdaptCell<float>&,
                                         const LevelSetCell<float, long>&,
                                         int, float, float, int,
                                         cell::TriangulationStrategy);

template void process_ready_to_cut_cells(AdaptCell<double>&,
                                         const LevelSetCell<double, int>&,
                                         int, double, double, int, bool);
template void process_ready_to_cut_cells(AdaptCell<float>&,
                                         const LevelSetCell<float, int>&,
                                         int, float, float, int, bool);
template void process_ready_to_cut_cells(AdaptCell<double>&,
                                         const LevelSetCell<double, long>&,
                                         int, double, double, int, bool);
template void process_ready_to_cut_cells(AdaptCell<float>&,
                                         const LevelSetCell<float, long>&,
                                         int, float, float, int, bool);

template bool refine_ready_cell_on_largest_midpoint_value(
    AdaptCell<double>&,
    const LevelSetCell<double, int>&,
    int,
    int);
template bool refine_ready_cell_on_largest_midpoint_value(
    AdaptCell<float>&,
    const LevelSetCell<float, int>&,
    int,
    int);
template bool refine_ready_cell_on_largest_midpoint_value(
    AdaptCell<double>&,
    const LevelSetCell<double, long>&,
    int,
    int);
template bool refine_ready_cell_on_largest_midpoint_value(
    AdaptCell<float>&,
    const LevelSetCell<float, long>&,
    int,
    int);

template void certify_and_refine(AdaptCell<double>&,
                                 const LevelSetCell<double, int>&,
                                 int, int, double, double, int);
template void certify_and_refine(AdaptCell<float>&,
                                 const LevelSetCell<float, int>&,
                                 int, int, float, float, int);
template void certify_and_refine(AdaptCell<double>&,
                                 const LevelSetCell<double, long>&,
                                 int, int, double, double, int);

template void certify_refine_and_process_ready_cells(AdaptCell<double>&,
                                                     const LevelSetCell<double, int>&,
                                                     int, int, double, double, int,
                                                     cell::TriangulationStrategy);
template void certify_refine_and_process_ready_cells(AdaptCell<float>&,
                                                     const LevelSetCell<float, int>&,
                                                     int, int, float, float, int,
                                                     cell::TriangulationStrategy);
template void certify_refine_and_process_ready_cells(AdaptCell<double>&,
                                                     const LevelSetCell<double, long>&,
                                                     int, int, double, double, int,
                                                     cell::TriangulationStrategy);
template void certify_refine_and_process_ready_cells(AdaptCell<float>&,
                                                     const LevelSetCell<float, long>&,
                                                     int, int, float, float, int,
                                                     cell::TriangulationStrategy);

template void certify_refine_and_process_ready_cells(AdaptCell<double>&,
                                                     const LevelSetCell<double, int>&,
                                                     int, int, double, double, int, bool);
template void certify_refine_and_process_ready_cells(AdaptCell<float>&,
                                                     const LevelSetCell<float, int>&,
                                                     int, int, float, float, int, bool);
template void certify_refine_and_process_ready_cells(AdaptCell<double>&,
                                                     const LevelSetCell<double, long>&,
                                                     int, int, double, double, int, bool);
template void certify_refine_and_process_ready_cells(AdaptCell<float>&,
                                                     const LevelSetCell<float, long>&,
                                                     int, int, float, float, int, bool);

} // namespace cutcells
