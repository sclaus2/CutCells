// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "cell_certification.h"
#include "bernstein.h"
#include "cell_flags.h"
#include "cell_topology.h"
#include "curving.h"
#include "cut_cell.h"
#include "cut_tetrahedron.h"
#include "cut_triangle.h"
#include "edge_certification.h"
#include "geometric_quantity.h"
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

    std::vector<T> edge_coeffs;
    gather_adapt_edge_bernstein(adapt_cell, ls_cell, edge_id, edge_coeffs);

    T root_t = T(0);
    if (!locate_one_root_parameter(std::span<const T>(edge_coeffs),
                                   zero_tol, sign_tol, edge_max_depth, root_t))
    {
        throw std::runtime_error(
            "ensure_one_root_vertex_on_edge: failed to localize unique edge root");
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
                    // extraction even when the curved level set evaluated at
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

template <std::floating_point T>
std::span<const T> point_span(const std::vector<T>& points,
                              int point_id,
                              int dim)
{
    return std::span<const T>(
        points.data()
            + static_cast<std::size_t>(point_id)
                  * static_cast<std::size_t>(dim),
        static_cast<std::size_t>(dim));
}

template <std::floating_point T>
geom::ParentEntity common_parent_entity_for_points(
    cell::type parent_cell_type,
    const std::vector<T>& points,
    int dim,
    T tol)
{
    const int npoints = static_cast<int>(
        points.size() / static_cast<std::size_t>(dim));

    if (npoints == 0)
        return {-1, -1};

    const auto edges = cell::edges(parent_cell_type);
    for (int e = 0; e < static_cast<int>(edges.size()); ++e)
    {
        bool all_on_edge = true;
        for (int i = 0; i < npoints; ++i)
            all_on_edge = all_on_edge
                       && geom::point_on_parent_edge<T>(
                              parent_cell_type, e, point_span<T>(points, i, dim), tol);
        if (all_on_edge)
            return {1, e};
    }

    if (cell::get_tdim(parent_cell_type) == 3)
    {
        for (int f = 0; f < cell::num_faces(parent_cell_type); ++f)
        {
            bool all_on_face = true;
            for (int i = 0; i < npoints; ++i)
                all_on_face = all_on_face
                           && geom::point_on_parent_face<T>(
                                  parent_cell_type, f, point_span<T>(points, i, dim), tol);
            if (all_on_face)
                return {2, f};
        }
    }

    return {cell::get_tdim(parent_cell_type), -1};
}

template <std::floating_point T>
geom::ParentEntity common_parent_face_or_cell_for_points(
    cell::type parent_cell_type,
    const std::vector<T>& points,
    int dim,
    T tol)
{
    const int npoints = static_cast<int>(
        points.size() / static_cast<std::size_t>(dim));
    if (cell::get_tdim(parent_cell_type) == 3)
    {
        for (int f = 0; f < cell::num_faces(parent_cell_type); ++f)
        {
            bool all_on_face = true;
            for (int i = 0; i < npoints; ++i)
                all_on_face = all_on_face
                           && geom::point_on_parent_face<T>(
                                  parent_cell_type, f, point_span<T>(points, i, dim), tol);
            if (all_on_face)
                return {2, f};
        }
    }
    return {cell::get_tdim(parent_cell_type), -1};
}

template <std::floating_point T>
void mark_graph_failure(ReadyCellGraphDiagnostics<T>& diagnostics,
                        int cell_id,
                        graph_criteria::FailureReason reason)
{
    diagnostics.accepted = false;
    ++diagnostics.failed_checks;
    if (diagnostics.first_failed_cell < 0)
    {
        diagnostics.first_failed_cell = cell_id;
        diagnostics.first_failure_reason = reason;
    }
}

template <std::floating_point T>
void mark_graph_refinement_request(ReadyCellGraphDiagnostics<T>& diagnostics,
                                   int entity_dim,
                                   int entity_id)
{
    if (entity_dim < 0 || entity_id < 0)
        return;
    if (diagnostics.first_requested_refinement_entity_dim < 0)
    {
        diagnostics.first_requested_refinement_entity_dim = entity_dim;
        diagnostics.first_requested_refinement_entity_id = entity_id;
    }
}

template <std::floating_point T>
void mark_graph_refinement_point(ReadyCellGraphDiagnostics<T>& diagnostics,
                                 std::span<const T> point)
{
    if (!diagnostics.first_requested_refinement_point.empty()
        || point.empty())
    {
        return;
    }
    diagnostics.first_requested_refinement_point.assign(
        point.begin(), point.end());
}

template <std::floating_point T>
void observe_graph_direction(ReadyCellGraphDiagnostics<T>& diagnostics,
                             const graph_criteria::DirectionReport<T>& report)
{
    diagnostics.min_true_transversality =
        std::min(diagnostics.min_true_transversality,
                 report.metrics.true_transversality);
    diagnostics.min_host_normal_alignment =
        std::min(diagnostics.min_host_normal_alignment,
                 report.metrics.host_normal_alignment);
    diagnostics.max_drift_amplification =
        std::max(diagnostics.max_drift_amplification,
                 report.metrics.drift_amplification);
    diagnostics.max_relative_correction_distance =
        std::max(diagnostics.max_relative_correction_distance,
                 report.metrics.relative_correction_distance);
    diagnostics.max_relative_tangential_shift =
        std::max(diagnostics.max_relative_tangential_shift,
                 report.metrics.relative_tangential_shift);
}

template <std::floating_point T>
void observe_graph_node(ReadyCellGraphDiagnostics<T>& diagnostics,
                        const GraphNodeDiagnostics<T>& node)
{
    diagnostics.nodes.push_back(node);
    if (std::isfinite(node.level_set_gradient_host_alignment))
    {
        diagnostics.min_level_set_gradient_host_alignment =
            std::min(diagnostics.min_level_set_gradient_host_alignment,
                     node.level_set_gradient_host_alignment);
    }
    if (!node.accepted)
    {
        mark_graph_refinement_request<T>(
            diagnostics,
            node.requested_refinement_entity_dim,
            node.requested_refinement_entity_id);
    }
}

template <std::floating_point T>
void observe_failed_graph_projection(
    ReadyCellGraphDiagnostics<T>& diagnostics,
    const curving::ProjectionDiagnostic<T>& projection)
{
    if (!diagnostics.first_failed_projection_seed.empty())
        return;
    diagnostics.first_failed_projection_seed = projection.seed;
    diagnostics.first_failed_projection_direction = projection.direction;
    diagnostics.first_failed_projection_clip_lo = projection.clip_lo;
    diagnostics.first_failed_projection_clip_hi = projection.clip_hi;
    diagnostics.first_failed_projection_root_t = projection.root_t;
}

template <std::floating_point T>
void observe_failed_face_orientation(
    ReadyCellGraphDiagnostics<T>& diagnostics,
    const graph_criteria::FaceQualityReport<T>& quality)
{
    if (diagnostics.first_failed_face_triangle_index < 0)
    {
        diagnostics.first_failed_face_triangle_index =
            quality.failed_triangle_index;
        diagnostics.first_failed_face_area_ratio =
            quality.failed_surface_jacobian_ratio;
    }
}

template <std::floating_point T>
int find_adapt_edge_by_vertices(const AdaptCell<T>& adapt_cell, int a, int b)
{
    if (a < 0 || b < 0)
        return -1;
    const int n_edges = adapt_cell.n_entities(1);
    for (int e = 0; e < n_edges; ++e)
    {
        auto ev = adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(e)];
        if (ev.size() != 2)
            continue;
        const int e0 = static_cast<int>(ev[0]);
        const int e1 = static_cast<int>(ev[1]);
        if ((e0 == a && e1 == b) || (e0 == b && e1 == a))
            return e;
    }
    return -1;
}

template <std::floating_point T>
int graph_refinement_edge_from_failed_face_segment(
    const AdaptCell<T>& adapt_cell,
    int a,
    int b)
{
    int edge_id = find_adapt_edge_by_vertices<T>(adapt_cell, a, b);
    if (edge_id >= 0)
        return edge_id;

    auto source_edge = [&](int v) -> int
    {
        if (v < 0
            || v >= static_cast<int>(adapt_cell.vertex_source_edge_id.size()))
        {
            return -1;
        }
        const int source = adapt_cell.vertex_source_edge_id[static_cast<std::size_t>(v)];
        return (source >= 0 && source < adapt_cell.n_entities(1)) ? source : -1;
    };

    edge_id = source_edge(a);
    if (edge_id >= 0)
        return edge_id;
    return source_edge(b);
}

template <std::floating_point T>
T graph_alignment_to_host_normal(std::span<const T> direction,
                                 std::span<const T> host_normal,
                                 T tol)
{
    if (direction.size() != host_normal.size())
        return std::numeric_limits<T>::quiet_NaN();
    const auto alignment = geom::alignment<T>(direction, host_normal, tol);
    if (alignment.degenerate)
        return T(0);
    return std::fabs(alignment.cosine);
}

template <std::floating_point T>
T graph_angle_to_tangent_degrees(T host_alignment)
{
    if (!std::isfinite(host_alignment))
        return std::numeric_limits<T>::quiet_NaN();
    const T clipped = std::clamp(host_alignment, T(0), T(1));
    const T pi = std::acos(T(-1));
    return std::asin(clipped) * T(180) / pi;
}

template <std::floating_point T>
bool solve_graph_dense_small(std::vector<T> A,
                             std::vector<T> b,
                             int n,
                             std::vector<T>& x,
                             T tol)
{
    if (n <= 0 || static_cast<int>(A.size()) != n * n
        || static_cast<int>(b.size()) != n)
    {
        return false;
    }

    for (int k = 0; k < n; ++k)
    {
        int pivot = k;
        T pivot_abs = std::fabs(A[static_cast<std::size_t>(k * n + k)]);
        for (int i = k + 1; i < n; ++i)
        {
            const T candidate =
                std::fabs(A[static_cast<std::size_t>(i * n + k)]);
            if (candidate > pivot_abs)
            {
                pivot = i;
                pivot_abs = candidate;
            }
        }
        if (pivot_abs <= tol)
            return false;

        if (pivot != k)
        {
            for (int j = k; j < n; ++j)
            {
                std::swap(A[static_cast<std::size_t>(k * n + j)],
                          A[static_cast<std::size_t>(pivot * n + j)]);
            }
            std::swap(b[static_cast<std::size_t>(k)],
                      b[static_cast<std::size_t>(pivot)]);
        }

        const T diag = A[static_cast<std::size_t>(k * n + k)];
        for (int i = k + 1; i < n; ++i)
        {
            const T factor = A[static_cast<std::size_t>(i * n + k)] / diag;
            A[static_cast<std::size_t>(i * n + k)] = T(0);
            for (int j = k + 1; j < n; ++j)
            {
                A[static_cast<std::size_t>(i * n + j)] -=
                    factor * A[static_cast<std::size_t>(k * n + j)];
            }
            b[static_cast<std::size_t>(i)] -=
                factor * b[static_cast<std::size_t>(k)];
        }
    }

    x.assign(static_cast<std::size_t>(n), T(0));
    for (int i = n - 1; i >= 0; --i)
    {
        T value = b[static_cast<std::size_t>(i)];
        for (int j = i + 1; j < n; ++j)
            value -= A[static_cast<std::size_t>(i * n + j)]
                   * x[static_cast<std::size_t>(j)];
        const T diag = A[static_cast<std::size_t>(i * n + i)];
        if (std::fabs(diag) <= tol)
            return false;
        x[static_cast<std::size_t>(i)] = value / diag;
    }
    return true;
}

template <std::floating_point T, std::integral I>
bool graph_affine_jacobian(const LevelSetCell<T, I>& ls_cell,
                           std::vector<T>& jacobian)
{
    const int tdim = ls_cell.tdim;
    const int gdim = ls_cell.gdim;
    if (tdim <= 0 || tdim > 3 || gdim <= 0
        || ls_cell.parent_vertex_coords.empty())
    {
        return false;
    }

    const auto cols = cell::jacobian_col_indices(ls_cell.cell_type);
    jacobian.assign(static_cast<std::size_t>(gdim * tdim), T(0));
    for (int a = 0; a < tdim; ++a)
    {
        const int va = cols[static_cast<std::size_t>(a)];
        if (va < 0)
            return false;
        for (int r = 0; r < gdim; ++r)
        {
            jacobian[static_cast<std::size_t>(r * tdim + a)] =
                ls_cell.parent_vertex_coords[
                    static_cast<std::size_t>(va * gdim + r)]
              - ls_cell.parent_vertex_coords[static_cast<std::size_t>(r)];
        }
    }
    return true;
}

template <std::floating_point T, std::integral I>
T graph_metric_alignment_to_host_normal(
    const LevelSetCell<T, I>& ls_cell,
    std::span<const T> direction_ref,
    const graph_criteria::HostFrame<T>& host,
    T tol)
{
    if (static_cast<int>(direction_ref.size()) != ls_cell.tdim)
    {
        return std::numeric_limits<T>::quiet_NaN();
    }

    std::vector<T> jacobian;
    if (!graph_affine_jacobian<T, I>(ls_cell, jacobian))
        return graph_alignment_to_host_normal<T>(
            direction_ref,
            std::span<const T>(host.normal.data(), host.normal.size()),
            tol);

    std::vector<T> direction_phys(static_cast<std::size_t>(ls_cell.gdim), T(0));
    for (int r = 0; r < ls_cell.gdim; ++r)
    {
        for (int a = 0; a < ls_cell.tdim; ++a)
        {
            const T j_ra = jacobian[static_cast<std::size_t>(
                r * ls_cell.tdim + a)];
            direction_phys[static_cast<std::size_t>(r)] +=
                j_ra * direction_ref[static_cast<std::size_t>(a)];
        }
    }

    const T direction_norm = geom::norm<T>(
        std::span<const T>(direction_phys.data(), direction_phys.size()));
    if (direction_norm <= tol)
        return T(0);

    if (host.dimension == graph_criteria::HostDimension::edge)
    {
        if (static_cast<int>(host.tangent.size()) != ls_cell.tdim)
            return std::numeric_limits<T>::quiet_NaN();

        std::vector<T> tangent_phys(static_cast<std::size_t>(ls_cell.gdim), T(0));
        for (int r = 0; r < ls_cell.gdim; ++r)
        {
            for (int a = 0; a < ls_cell.tdim; ++a)
            {
                tangent_phys[static_cast<std::size_t>(r)] +=
                    jacobian[static_cast<std::size_t>(r * ls_cell.tdim + a)]
                  * host.tangent[static_cast<std::size_t>(a)];
            }
        }

        const T tangent_norm = geom::norm<T>(
            std::span<const T>(tangent_phys.data(), tangent_phys.size()));
        if (tangent_norm <= tol)
            return std::numeric_limits<T>::quiet_NaN();

        const T cosine = geom::dot<T>(
            std::span<const T>(direction_phys.data(), direction_phys.size()),
            std::span<const T>(tangent_phys.data(), tangent_phys.size()))
          / (direction_norm * tangent_norm);
        const T clipped = std::clamp(cosine, T(-1), T(1));
        return std::sqrt(std::max(T(0), T(1) - clipped * clipped));
    }

    if (static_cast<int>(host.normal.size()) != ls_cell.tdim
        || ls_cell.gdim != ls_cell.tdim)
    {
        return graph_alignment_to_host_normal<T>(
            direction_ref,
            std::span<const T>(host.normal.data(), host.normal.size()),
            tol);
    }

    std::vector<T> jt(static_cast<std::size_t>(ls_cell.tdim * ls_cell.tdim), T(0));
    for (int a = 0; a < ls_cell.tdim; ++a)
    {
        for (int r = 0; r < ls_cell.gdim; ++r)
        {
            jt[static_cast<std::size_t>(a * ls_cell.tdim + r)] =
                jacobian[static_cast<std::size_t>(r * ls_cell.tdim + a)];
        }
    }
    std::vector<T> normal_phys;
    std::vector<T> rhs(host.normal.begin(), host.normal.end());
    if (!solve_graph_dense_small<T>(
            std::move(jt), std::move(rhs), ls_cell.tdim, normal_phys, tol))
    {
        return graph_alignment_to_host_normal<T>(
            direction_ref,
            std::span<const T>(host.normal.data(), host.normal.size()),
            tol);
    }

    return graph_alignment_to_host_normal<T>(
        std::span<const T>(direction_phys.data(), direction_phys.size()),
        std::span<const T>(normal_phys.data(), normal_phys.size()),
        tol);
}

template <std::floating_point T>
std::vector<T> graph_pull_back_physical_direction(
    std::span<const T> jacobian,
    int gdim,
    int tdim,
    std::span<const T> direction_phys,
    std::span<const T> fallback_ref,
    T tol)
{
    std::vector<T> gram(static_cast<std::size_t>(tdim * tdim), T(0));
    std::vector<T> rhs(static_cast<std::size_t>(tdim), T(0));
    for (int a = 0; a < tdim; ++a)
    {
        for (int r = 0; r < gdim; ++r)
        {
            rhs[static_cast<std::size_t>(a)] +=
                jacobian[static_cast<std::size_t>(r * tdim + a)]
              * direction_phys[static_cast<std::size_t>(r)];
        }
        for (int b = 0; b < tdim; ++b)
        {
            T value = T(0);
            for (int r = 0; r < gdim; ++r)
            {
                value += jacobian[static_cast<std::size_t>(r * tdim + a)]
                       * jacobian[static_cast<std::size_t>(r * tdim + b)];
            }
            gram[static_cast<std::size_t>(a * tdim + b)] = value;
        }
    }

    std::vector<T> out;
    if (solve_graph_dense_small<T>(
            std::move(gram), std::move(rhs), tdim, out, tol))
    {
        return out;
    }
    return std::vector<T>(fallback_ref.begin(), fallback_ref.end());
}

template <std::floating_point T, std::integral I>
std::vector<T> graph_parent_entity_metric_gradient_direction(
    const LevelSetCell<T, I>& ls_cell,
    geom::ParentEntity parent_entity,
    std::span<const T> grad_ref,
    T tol)
{
    const auto fallback = geom::restricted_level_set_gradient_in_parent_frame<T>(
        ls_cell.cell_type, parent_entity, grad_ref, tol);
    std::vector<T> fallback_ref =
        fallback.value.empty()
            ? std::vector<T>(grad_ref.begin(), grad_ref.end())
            : fallback.value;

    std::vector<T> jacobian;
    if (!graph_affine_jacobian<T, I>(ls_cell, jacobian))
        return fallback_ref;

    const int tdim = ls_cell.tdim;
    const int gdim = ls_cell.gdim;

    std::vector<T> grad_phys;
    if (gdim == tdim)
    {
        std::vector<T> jt(static_cast<std::size_t>(tdim * tdim), T(0));
        for (int a = 0; a < tdim; ++a)
        {
            for (int r = 0; r < gdim; ++r)
            {
                jt[static_cast<std::size_t>(a * tdim + r)] =
                    jacobian[static_cast<std::size_t>(r * tdim + a)];
            }
        }
        std::vector<T> rhs(grad_ref.begin(), grad_ref.end());
        if (!solve_graph_dense_small<T>(
                std::move(jt), std::move(rhs), tdim, grad_phys, tol))
        {
            return fallback_ref;
        }
    }
    else
    {
        grad_phys.assign(static_cast<std::size_t>(gdim), T(0));
        for (int r = 0; r < gdim; ++r)
        {
            for (int a = 0; a < tdim; ++a)
            {
                grad_phys[static_cast<std::size_t>(r)] +=
                    jacobian[static_cast<std::size_t>(r * tdim + a)]
                  * fallback_ref[static_cast<std::size_t>(a)];
            }
        }
    }

    if (parent_entity.dim == 2 && gdim == 3)
    {
        const auto face = cell::face_vertices(ls_cell.cell_type, parent_entity.id);
        if (face.size() >= 3)
        {
            std::array<T, 3> e0 = {T(0), T(0), T(0)};
            std::array<T, 3> e1 = {T(0), T(0), T(0)};
            for (int d = 0; d < 3; ++d)
            {
                e0[static_cast<std::size_t>(d)] =
                    ls_cell.parent_vertex_coords[
                        static_cast<std::size_t>(face[1] * gdim + d)]
                  - ls_cell.parent_vertex_coords[
                        static_cast<std::size_t>(face[0] * gdim + d)];
                e1[static_cast<std::size_t>(d)] =
                    ls_cell.parent_vertex_coords[
                        static_cast<std::size_t>(face[2] * gdim + d)]
                  - ls_cell.parent_vertex_coords[
                        static_cast<std::size_t>(face[0] * gdim + d)];
            }
            std::array<T, 3> normal = {
                e0[1] * e1[2] - e0[2] * e1[1],
                e0[2] * e1[0] - e0[0] * e1[2],
                e0[0] * e1[1] - e0[1] * e1[0]};
            const T nn = normal[0] * normal[0]
                       + normal[1] * normal[1]
                       + normal[2] * normal[2];
            if (nn > tol * tol)
            {
                const T c = (grad_phys[0] * normal[0]
                           + grad_phys[1] * normal[1]
                           + grad_phys[2] * normal[2]) / nn;
                for (int d = 0; d < 3; ++d)
                    grad_phys[static_cast<std::size_t>(d)] -= c * normal[d];
            }
        }
    }
    else if (parent_entity.dim == 1)
    {
        const auto edges = cell::edges(ls_cell.cell_type);
        if (parent_entity.id >= 0
            && parent_entity.id < static_cast<int>(edges.size()))
        {
            const auto edge = edges[static_cast<std::size_t>(parent_entity.id)];
            std::vector<T> tangent(static_cast<std::size_t>(gdim), T(0));
            for (int d = 0; d < gdim; ++d)
            {
                tangent[static_cast<std::size_t>(d)] =
                    ls_cell.parent_vertex_coords[
                        static_cast<std::size_t>(edge[1] * gdim + d)]
                  - ls_cell.parent_vertex_coords[
                        static_cast<std::size_t>(edge[0] * gdim + d)];
            }
            const T tt = geom::dot<T>(
                std::span<const T>(tangent.data(), tangent.size()),
                std::span<const T>(tangent.data(), tangent.size()));
            if (tt > tol * tol)
            {
                const T c = geom::dot<T>(
                    std::span<const T>(grad_phys.data(), grad_phys.size()),
                    std::span<const T>(tangent.data(), tangent.size())) / tt;
                for (int d = 0; d < gdim; ++d)
                    grad_phys[static_cast<std::size_t>(d)] =
                        c * tangent[static_cast<std::size_t>(d)];
            }
        }
    }
    else if (parent_entity.dim == 0)
    {
        return std::vector<T>(static_cast<std::size_t>(tdim), T(0));
    }

    auto pulled = graph_pull_back_physical_direction<T>(
        std::span<const T>(jacobian.data(), jacobian.size()),
        gdim,
        tdim,
        std::span<const T>(grad_phys.data(), grad_phys.size()),
        std::span<const T>(fallback_ref.data(), fallback_ref.size()),
        tol);
    const auto admissible = geom::admissible_direction_in_parent_frame<T>(
        ls_cell.cell_type,
        parent_entity,
        std::span<const T>(pulled.data(), pulled.size()),
        tol);
    return admissible.value.empty() ? pulled : admissible.value;
}

template <std::floating_point T>
curving::CurvingOptions<T> graph_curving_options(
    const ReadyCellGraphOptions<T>& graph_options)
{
    curving::CurvingOptions<T> options;
    options.geometry_order = std::max(graph_options.geometry_order, 1);
    options.direction_mode =
        (graph_options.projection_mode
         == GraphProjectionMode::straight_zero_entity_normal)
            ? curving::CurvingDirectionMode::straight_zero_entity_normal
            : curving::CurvingDirectionMode::level_set_gradient;
    options.domain_tol = graph_options.criteria.tolerance;
    options.active_face_tol = graph_options.criteria.tolerance;
    return options;
}

inline graph_criteria::FailureReason graph_failure_from_curving(
    curving::CurvingFailureCode code)
{
    switch (code)
    {
    case curving::CurvingFailureCode::none:
    case curving::CurvingFailureCode::exact_vertex:
    case curving::CurvingFailureCode::boundary_from_edge:
    case curving::CurvingFailureCode::small_entity_kept_straight:
        return graph_criteria::FailureReason::none;
    case curving::CurvingFailureCode::invalid_constraint_count:
    case curving::CurvingFailureCode::missing_level_set_cell:
    case curving::CurvingFailureCode::empty_zero_mask:
    case curving::CurvingFailureCode::unsupported_entity:
    case curving::CurvingFailureCode::missing_boundary_edge:
    case curving::CurvingFailureCode::boundary_edge_failed:
        return graph_criteria::FailureReason::invalid_input;
    case curving::CurvingFailureCode::no_host_interval:
    case curving::CurvingFailureCode::outside_host_domain:
        return graph_criteria::FailureReason::root_segment_leaves_parent_entity;
    case curving::CurvingFailureCode::singular_gradient_system:
        return graph_criteria::FailureReason::degenerate_direction;
    case curving::CurvingFailureCode::no_sign_changing_bracket:
    case curving::CurvingFailureCode::brent_failed:
    case curving::CurvingFailureCode::line_search_failed:
    case curving::CurvingFailureCode::max_iterations:
    case curving::CurvingFailureCode::projection_failed:
    case curving::CurvingFailureCode::closest_face_retry_failed:
    case curving::CurvingFailureCode::constrained_newton_failed:
        return graph_criteria::FailureReason::root_not_on_search_line;
    }
    return graph_criteria::FailureReason::root_not_on_search_line;
}

template <std::floating_point T, std::integral I>
bool project_graph_point(
    const AdaptCell<T>& adapt_cell,
    int local_zero_entity_id,
    const LevelSetCell<T, I>& ls_cell,
    geom::ParentEntity admissible_parent_entity,
    const graph_criteria::HostFrame<T>& host,
    std::span<const T> host_point,
    GraphNodeKind node_kind,
    int node_index,
    int requested_refinement_edge_id,
    const ReadyCellGraphOptions<T>& graph_options,
    ReadyCellGraphDiagnostics<T>& diagnostics,
    int source_cell_id,
    std::vector<T>& corrected_point)
{
    GraphNodeDiagnostics<T> node;
    node.node_index = node_index;
    node.node_kind = node_kind;
    node.parent_entity_dim = admissible_parent_entity.dim;
    node.parent_entity_id = admissible_parent_entity.id;
    node.seed.assign(host_point.begin(), host_point.end());
    node.straight_helper_normal.assign(host.normal.begin(), host.normal.end());
    if (requested_refinement_edge_id >= 0)
    {
        node.requested_refinement_entity_dim = 1;
        node.requested_refinement_entity_id = requested_refinement_edge_id;
    }
    node.selected_direction_kind =
        (graph_options.projection_mode
         == GraphProjectionMode::straight_zero_entity_normal)
            ? graph_criteria::DirectionKind::projected_straight_host_normal
            : graph_criteria::DirectionKind::projected_level_set_gradient;

    std::vector<T> gradient(static_cast<std::size_t>(ls_cell.tdim), T(0));
    curving::reference_level_set_gradient<T, I>(
        ls_cell,
        host_point,
        std::span<T>(gradient.data(), gradient.size()));
    auto graph_gradient_direction =
        graph_parent_entity_metric_gradient_direction<T, I>(
            ls_cell,
            admissible_parent_entity,
            std::span<const T>(gradient.data(), gradient.size()),
            graph_options.criteria.tolerance);
    node.level_set_gradient_direction = graph_gradient_direction;
    node.level_set_gradient_host_alignment =
        graph_metric_alignment_to_host_normal<T, I>(
        ls_cell,
        std::span<const T>(
            graph_gradient_direction.data(), graph_gradient_direction.size()),
        host,
        graph_options.criteria.tolerance);
    node.level_set_gradient_angle_to_tangent_deg =
        graph_angle_to_tangent_degrees<T>(
            node.level_set_gradient_host_alignment);
    if (std::isfinite(node.level_set_gradient_host_alignment)
        && node.level_set_gradient_host_alignment
               < graph_options.min_level_set_gradient_host_alignment)
    {
        node.accepted = false;
        node.failure_reason =
            graph_criteria::FailureReason::direction_too_tangential_to_host;
        observe_graph_node<T>(diagnostics, node);
        mark_graph_failure<T>(
            diagnostics, source_cell_id, node.failure_reason);
        return false;
    }

    std::vector<T> straight_normal_direction;
    const auto admissible_normal =
        geom::admissible_direction_in_parent_frame<T>(
            ls_cell.cell_type,
            admissible_parent_entity,
            std::span<const T>(host.normal.data(), host.normal.size()),
            graph_options.criteria.tolerance);
    if (!admissible_normal.degenerate())
        straight_normal_direction = admissible_normal.value;

    const T grad_norm = geom::norm<T>(
        std::span<const T>(
            graph_gradient_direction.data(), graph_gradient_direction.size()));
    const T normal_norm = geom::norm<T>(
        std::span<const T>(
            straight_normal_direction.data(), straight_normal_direction.size()));
    std::vector<T> selected_policy_direction;
    if (graph_options.projection_mode == GraphProjectionMode::level_set_gradient)
    {
        if (grad_norm > graph_options.criteria.tolerance)
        {
            selected_policy_direction = graph_gradient_direction;
            node.selected_direction_kind =
                graph_criteria::DirectionKind::projected_level_set_gradient;
        }
        else
        {
            selected_policy_direction = straight_normal_direction;
            node.selected_direction_kind =
                graph_criteria::DirectionKind::projected_straight_host_normal;
            node.fallback_used = true;
        }
    }
    else
    {
        if (normal_norm > graph_options.criteria.tolerance)
        {
            selected_policy_direction = straight_normal_direction;
            node.selected_direction_kind =
                graph_criteria::DirectionKind::projected_straight_host_normal;
        }
        else
        {
            selected_policy_direction = graph_gradient_direction;
            node.selected_direction_kind =
                graph_criteria::DirectionKind::projected_level_set_gradient;
            node.fallback_used = true;
        }
    }

    const auto curving_options = graph_curving_options<T>(graph_options);
    auto projection = curving::project_seed_to_zero_entity_diagnostic<T, I>(
        adapt_cell,
        local_zero_entity_id,
        ls_cell,
        host_point,
        curving_options);

    if (!projection.accepted)
    {
        node.accepted = false;
        node.failure_reason = graph_failure_from_curving(projection.failure_code);
        node.selected_direction = selected_policy_direction;
        node.corrected = projection.projected;
        observe_graph_node<T>(diagnostics, node);
        observe_failed_graph_projection<T>(diagnostics, projection);
        mark_graph_failure<T>(diagnostics, source_cell_id, node.failure_reason);
        return false;
    }

    corrected_point = std::move(projection.projected);
    node.corrected = corrected_point;
    if (corrected_point.size() != host_point.size())
    {
        node.accepted = false;
        node.failure_reason = graph_criteria::FailureReason::invalid_input;
        observe_graph_node<T>(diagnostics, node);
        mark_graph_failure<T>(
            diagnostics,
            source_cell_id,
            graph_criteria::FailureReason::invalid_input);
        return false;
    }

    std::vector<T> direction = std::move(selected_policy_direction);
    if (geom::norm<T>(std::span<const T>(direction.data(), direction.size()))
        <= graph_options.criteria.tolerance)
    {
        direction.resize(host_point.size(), T(0));
        for (std::size_t d = 0; d < host_point.size(); ++d)
            direction[d] = corrected_point[d] - host_point[d];
    }

    if (geom::norm<T>(std::span<const T>(direction.data(), direction.size()))
        <= graph_options.criteria.tolerance)
    {
        node.accepted = true;
        node.failure_reason = graph_criteria::FailureReason::none;
        node.selected_direction = direction;
        observe_graph_node<T>(diagnostics, node);
        return true;
    }
    node.selected_direction = direction;

    auto report = graph_criteria::evaluate_direction<T>(
        ls_cell.cell_type,
        admissible_parent_entity,
        host,
        host_point,
        std::span<const T>(
            graph_gradient_direction.data(), graph_gradient_direction.size()),
        std::span<const T>(direction.data(), direction.size()),
        std::span<const T>(corrected_point.data(), corrected_point.size()),
        node.selected_direction_kind,
        graph_options.criteria);
    const T metric_selected_alignment =
        graph_metric_alignment_to_host_normal<T, I>(
            ls_cell,
            std::span<const T>(direction.data(), direction.size()),
            host,
            graph_options.criteria.tolerance);
    if (std::isfinite(metric_selected_alignment))
    {
        report.metrics.host_normal_alignment = metric_selected_alignment;
        report.metrics.drift_amplification =
            graph_criteria::drift_from_alignment<T>(
                metric_selected_alignment,
                graph_options.criteria.tolerance);
        if ((report.failure_reason
             == graph_criteria::FailureReason::excessive_drift_amplification
             || report.failure_reason
                    == graph_criteria::FailureReason::direction_too_tangential_to_host)
            && report.metrics.drift_amplification
                   <= graph_options.criteria.max_drift_amplification
            && report.metrics.host_normal_alignment
                   >= graph_options.criteria.min_host_normal_alignment)
        {
            report.accepted = true;
            report.failure_reason = graph_criteria::FailureReason::none;
        }
    }
    observe_graph_direction<T>(diagnostics, report);
    node.true_transversality = report.metrics.true_transversality;
    node.selected_host_alignment = report.metrics.host_normal_alignment;
    node.drift_amplification = report.metrics.drift_amplification;
    node.relative_correction_distance =
        report.metrics.relative_correction_distance;
    node.relative_tangential_shift =
        report.metrics.relative_tangential_shift;
    if (!report.accepted)
    {
        if (report.failure_reason
            == graph_criteria::FailureReason::root_not_on_search_line)
        {
            node.accepted = true;
            node.failure_reason = graph_criteria::FailureReason::none;
            observe_graph_node<T>(diagnostics, node);
            return true;
        }
        node.accepted = false;
        node.failure_reason = report.failure_reason;
        observe_graph_node<T>(diagnostics, node);
        observe_failed_graph_projection<T>(diagnostics, projection);
        mark_graph_failure<T>(
            diagnostics, source_cell_id, report.failure_reason);
        return false;
    }
    node.accepted = true;
    node.failure_reason = graph_criteria::FailureReason::none;
    observe_graph_node<T>(diagnostics, node);
    return true;
}

template <std::floating_point T>
bool build_edge_host_frame(std::span<const T> a,
                           std::span<const T> b,
                           std::span<const T> zero_face_normal,
                           graph_criteria::HostFrame<T>& host,
                           T tol)
{
    const auto tangent = geom::segment_tangent<T>(a, b, true, tol);
    if (tangent.degenerate())
        return false;

    host.dimension = graph_criteria::HostDimension::edge;
    host.tangent = tangent.value;
    host.h = tangent.norm;

    if (a.size() == 2)
    {
        const auto normal = geom::segment_normal<T>(a, b, true, tol);
        if (normal.degenerate())
            return false;
        host.normal = normal.value;
        return true;
    }

    if (zero_face_normal.size() == 3)
    {
        const auto normal = geom::in_face_segment_normal<T>(
            a, b, zero_face_normal, true, tol);
        if (normal.degenerate())
            return false;
        host.normal = normal.value;
        return true;
    }

    return false;
}

template <std::floating_point T>
bool build_face_host_frame(const std::vector<T>& face_vertices,
                           graph_criteria::HostFrame<T>& host,
                           T tol)
{
    const int dim = 3;
    const int nverts = static_cast<int>(
        face_vertices.size() / static_cast<std::size_t>(dim));
    if (nverts < 3)
        return false;

    const auto normal = geom::face_normal<T>(
        point_span<T>(face_vertices, 0, dim),
        point_span<T>(face_vertices, 1, dim),
        point_span<T>(face_vertices, 2, dim),
        true,
        tol);
    if (normal.degenerate())
        return false;

    T h = T(0);
    for (int i = 0; i < nverts; ++i)
    {
        const auto a = point_span<T>(face_vertices, i, dim);
        const auto b = point_span<T>(face_vertices, (i + 1) % nverts, dim);
        h = std::max(h, geom::norm<T>(
            std::span<const T>(geom::subtract<T>(b, a).data(), dim)));
    }

    host.dimension = graph_criteria::HostDimension::face;
    host.normal = normal.value;
    host.tangent.clear();
    host.h = h;
    return h > tol;
}

template <std::floating_point T>
std::vector<int> zero_entity_host_cell_vertex_ids(const AdaptCell<T>& ac,
                                                  int local_zero_entity_id)
{
    if (local_zero_entity_id < 0
        || local_zero_entity_id >= ac.zero_entity_host_cell_vertices.size())
    {
        return {};
    }
    const auto vertices = ac.zero_entity_host_cell_vertices[
        static_cast<std::int32_t>(local_zero_entity_id)];
    std::vector<int> out;
    out.reserve(vertices.size());
    for (const auto v : vertices)
        out.push_back(static_cast<int>(v));
    return out;
}

template <std::floating_point T>
bool point_in_host_triangle(const AdaptCell<T>& ac,
                            std::span<const int> triangle_vertices,
                            std::span<const T> point,
                            T tol)
{
    if (ac.tdim != 3 || triangle_vertices.size() != 3 || point.size() != 3)
        return false;

    const auto a = point_span<T>(ac.vertex_coords, triangle_vertices[0], 3);
    const auto b = point_span<T>(ac.vertex_coords, triangle_vertices[1], 3);
    const auto c = point_span<T>(ac.vertex_coords, triangle_vertices[2], 3);
    std::array<T, 3> v0 = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    std::array<T, 3> v1 = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    std::array<T, 3> v2 = {point[0] - a[0], point[1] - a[1], point[2] - a[2]};
    const auto normal = geom::cross<T>(
        std::span<const T>(v0.data(), 3),
        std::span<const T>(v1.data(), 3));
    const T area2 = geom::norm<T>(
        std::span<const T>(normal.data(), normal.size()));
    if (area2 <= tol)
        return false;

    const T plane =
        std::fabs(geom::dot<T>(
            std::span<const T>(v2.data(), 3),
            std::span<const T>(normal.data(), normal.size()))) / area2;
    if (plane > tol)
        return false;

    const T d00 = geom::dot<T>(
        std::span<const T>(v0.data(), 3),
        std::span<const T>(v0.data(), 3));
    const T d01 = geom::dot<T>(
        std::span<const T>(v0.data(), 3),
        std::span<const T>(v1.data(), 3));
    const T d11 = geom::dot<T>(
        std::span<const T>(v1.data(), 3),
        std::span<const T>(v1.data(), 3));
    const T d20 = geom::dot<T>(
        std::span<const T>(v2.data(), 3),
        std::span<const T>(v0.data(), 3));
    const T d21 = geom::dot<T>(
        std::span<const T>(v2.data(), 3),
        std::span<const T>(v1.data(), 3));
    const T denom = d00 * d11 - d01 * d01;
    if (std::fabs(denom) <= tol * tol)
        return false;

    const T v = (d11 * d20 - d01 * d21) / denom;
    const T w = (d00 * d21 - d01 * d20) / denom;
    const T u = T(1) - v - w;
    return u >= -tol && v >= -tol && w >= -tol
        && u <= T(1) + tol && v <= T(1) + tol && w <= T(1) + tol;
}

template <std::floating_point T>
std::vector<T> host_boundary_face_normal_for_zero_face_edge(
    const AdaptCell<T>& ac,
    int zero_face_id,
    std::span<const T> edge_a,
    std::span<const T> edge_b,
    T tol,
    int* host_face_id = nullptr)
{
    if (zero_face_id < 0 || zero_face_id >= ac.n_zero_entities()
        || ac.tdim != 3 || edge_a.size() != 3 || edge_b.size() != 3)
    {
        return {};
    }
    if (zero_face_id >= static_cast<int>(ac.zero_entity_host_cell_type.size()))
        return {};

    const auto host_vertices =
        zero_entity_host_cell_vertex_ids<T>(ac, zero_face_id);
    if (host_vertices.empty())
        return {};

    const cell::type host_type =
        ac.zero_entity_host_cell_type[static_cast<std::size_t>(zero_face_id)];
    for (int f = 0; f < cell::num_faces(host_type); ++f)
    {
        const auto local_face = cell::face_vertices(host_type, f);
        if (local_face.size() != 3)
            continue;
        std::array<int, 3> face_vertices = {
            host_vertices[static_cast<std::size_t>(local_face[0])],
            host_vertices[static_cast<std::size_t>(local_face[1])],
            host_vertices[static_cast<std::size_t>(local_face[2])]
        };
        if (!point_in_host_triangle<T>(
                ac,
                std::span<const int>(face_vertices.data(), face_vertices.size()),
                edge_a,
                tol)
            || !point_in_host_triangle<T>(
                ac,
                std::span<const int>(face_vertices.data(), face_vertices.size()),
                edge_b,
                tol))
        {
            continue;
        }

        const auto normal = geom::face_normal<T>(
            point_span<T>(ac.vertex_coords, face_vertices[0], 3),
            point_span<T>(ac.vertex_coords, face_vertices[1], 3),
            point_span<T>(ac.vertex_coords, face_vertices[2], 3),
            true,
            tol);
        if (!normal.degenerate())
        {
            if (host_face_id != nullptr)
                *host_face_id = f;
            return normal.value;
        }
    }
    return {};
}

template <std::floating_point T>
int find_zero_edge_entity_by_vertices(const AdaptCell<T>& ac,
                                      int vertex_a,
                                      int vertex_b,
                                      std::uint64_t zero_mask)
{
    for (int z = 0; z < ac.n_zero_entities(); ++z)
    {
        if (ac.zero_entity_dim[static_cast<std::size_t>(z)] != 1)
            continue;
        if ((ac.zero_entity_zero_mask[static_cast<std::size_t>(z)] & zero_mask)
            != zero_mask)
        {
            continue;
        }

        const int edge_id = ac.zero_entity_id[static_cast<std::size_t>(z)];
        const auto edge_vertices =
            ac.entity_to_vertex[1][static_cast<std::int32_t>(edge_id)];
        if (edge_vertices.size() != 2)
            continue;
        const int a = static_cast<int>(edge_vertices[0]);
        const int b = static_cast<int>(edge_vertices[1]);
        if ((a == vertex_a && b == vertex_b)
            || (a == vertex_b && b == vertex_a))
        {
            return z;
        }
    }
    return -1;
}

template <std::floating_point T>
void override_zero_edge_host_from_zero_face(AdaptCell<T>& ac,
                                            int zero_edge_id,
                                            int zero_face_id,
                                            int host_face_id)
{
    if (zero_edge_id < 0 || zero_face_id < 0
        || zero_edge_id >= ac.n_zero_entities()
        || zero_face_id >= ac.n_zero_entities())
    {
        return;
    }
    if (zero_edge_id >= static_cast<int>(ac.zero_entity_host_cell_id.size())
        || zero_face_id >= static_cast<int>(ac.zero_entity_host_cell_id.size()))
    {
        return;
    }

    const auto face_host_vertices =
        zero_entity_host_cell_vertex_ids<T>(ac, zero_face_id);
    if (face_host_vertices.empty())
        return;

    ac.zero_entity_host_cell_id[static_cast<std::size_t>(zero_edge_id)] =
        ac.zero_entity_host_cell_id[static_cast<std::size_t>(zero_face_id)];
    ac.zero_entity_host_cell_type[static_cast<std::size_t>(zero_edge_id)] =
        ac.zero_entity_host_cell_type[static_cast<std::size_t>(zero_face_id)];
    ac.zero_entity_host_face_id[static_cast<std::size_t>(zero_edge_id)] =
        static_cast<std::int32_t>(host_face_id);
    ac.zero_entity_source_level_set[static_cast<std::size_t>(zero_edge_id)] =
        ac.zero_entity_source_level_set[static_cast<std::size_t>(zero_face_id)];

    EntityAdjacency updated;
    updated.offsets.push_back(std::int32_t(0));
    for (int z = 0; z < ac.n_zero_entities(); ++z)
    {
        if (z == zero_edge_id)
        {
            for (const int v : face_host_vertices)
                updated.indices.push_back(static_cast<std::int32_t>(v));
        }
        else if (z < ac.zero_entity_host_cell_vertices.size())
        {
            const auto old = ac.zero_entity_host_cell_vertices[
                static_cast<std::int32_t>(z)];
            for (const auto v : old)
                updated.indices.push_back(v);
        }
        updated.offsets.push_back(
            static_cast<std::int32_t>(updated.indices.size()));
    }
    ac.zero_entity_host_cell_vertices = std::move(updated);
}

template <std::floating_point T, std::integral I>
bool check_zero_edge_graph(const AdaptCell<T>& adapt_cell,
                           int local_zero_entity_id,
                           const LevelSetCell<T, I>& ls_cell,
                           std::span<const T> a,
                           std::span<const T> b,
                           std::span<const T> zero_face_normal,
                           bool face_boundary_edge_check,
                           int requested_refinement_edge_id,
                           const ReadyCellGraphOptions<T>& graph_options,
                           ReadyCellGraphDiagnostics<T>& diagnostics,
                           int source_cell_id)
{
    ++diagnostics.checked_edges;

    std::vector<T> host_points;
    const int order = std::max(graph_options.geometry_order, 1);
    for (int k = 0; k <= order; ++k)
    {
        const T s = T(k) / T(order);
        for (std::size_t d = 0; d < a.size(); ++d)
            host_points.push_back((T(1) - s) * a[d] + s * b[d]);
    }

    const auto parent_entity =
        (ls_cell.tdim == 3 && zero_face_normal.size() == 3)
            ? common_parent_face_or_cell_for_points<T>(
                  ls_cell.cell_type,
                  host_points,
                  ls_cell.tdim,
                  graph_options.criteria.tolerance)
            : common_parent_entity_for_points<T>(
                  ls_cell.cell_type,
                  host_points,
                  ls_cell.tdim,
                  graph_options.criteria.tolerance);

    std::vector<T> frame_normal(zero_face_normal.begin(), zero_face_normal.end());
    if (ls_cell.tdim == 3 && parent_entity.dim == 2)
    {
        const auto parent_normal = geom::parent_face_normal<T>(
            ls_cell.cell_type,
            parent_entity.id,
            true,
            graph_options.criteria.tolerance);
        if (!parent_normal.degenerate())
            frame_normal = parent_normal.value;
    }

    graph_criteria::HostFrame<T> host;
    if (!build_edge_host_frame<T>(
            a,
            b,
            std::span<const T>(frame_normal.data(), frame_normal.size()),
            host,
            graph_options.criteria.tolerance))
    {
        mark_graph_failure<T>(
            diagnostics,
            source_cell_id,
            graph_criteria::FailureReason::invalid_host_frame);
        return false;
    }

    std::vector<T> corrected_points;
    corrected_points.reserve(host_points.size());
    for (int i = 0; i <= order; ++i)
    {
        if (i == 0 || i == order)
        {
            auto endpoint = point_span<T>(host_points, i, ls_cell.tdim);
            corrected_points.insert(
                corrected_points.end(), endpoint.begin(), endpoint.end());
            continue;
        }

        const GraphNodeKind node_kind =
            face_boundary_edge_check
                ? GraphNodeKind::face_boundary_edge
                : GraphNodeKind::edge_interior;
        std::vector<T> corrected;
        if (!project_graph_point<T, I>(
                adapt_cell,
                local_zero_entity_id,
                ls_cell,
                parent_entity,
                host,
                point_span<T>(host_points, i, ls_cell.tdim),
                node_kind,
                i,
                requested_refinement_edge_id,
                graph_options,
                diagnostics,
                source_cell_id,
                corrected))
        {
            return false;
        }
        corrected_points.insert(
            corrected_points.end(), corrected.begin(), corrected.end());
    }

    const auto ordering = graph_criteria::evaluate_projected_edge_ordering<T>(
        std::span<const T>(host_points.data(), host_points.size()),
        std::span<const T>(corrected_points.data(), corrected_points.size()),
        ls_cell.tdim,
        std::span<const T>(host.tangent.data(), host.tangent.size()),
        host.h,
        graph_options.criteria);
    diagnostics.min_edge_gap_ratio =
        std::min(diagnostics.min_edge_gap_ratio, ordering.minimum_gap_ratio);
    if (!ordering.accepted)
    {
        mark_graph_failure<T>(
            diagnostics, source_cell_id, ordering.failure_reason);
        return false;
    }

    return true;
}

template <std::floating_point T>
std::vector<T> graph_push_forward_vector(const std::vector<T>& jacobian,
                                         int gdim,
                                         int tdim,
                                         std::span<const T> vector_ref);

template <std::floating_point T, std::integral I>
bool graph_physical_level_set_gradient(
    const LevelSetCell<T, I>& ls_cell,
    std::span<const T> point_ref,
    std::vector<T>& gradient_phys,
    T tol)
{
    if (static_cast<int>(point_ref.size()) != ls_cell.tdim)
        return false;

    std::vector<T> gradient_ref(static_cast<std::size_t>(ls_cell.tdim), T(0));
    curving::reference_level_set_gradient<T, I>(
        ls_cell,
        point_ref,
        std::span<T>(gradient_ref.data(), gradient_ref.size()));

    std::vector<T> jacobian;
    if (!graph_affine_jacobian<T, I>(ls_cell, jacobian))
        return false;

    if (ls_cell.gdim == ls_cell.tdim)
    {
        std::vector<T> jt(
            static_cast<std::size_t>(ls_cell.tdim * ls_cell.tdim), T(0));
        for (int a = 0; a < ls_cell.tdim; ++a)
        {
            for (int r = 0; r < ls_cell.gdim; ++r)
            {
                jt[static_cast<std::size_t>(a * ls_cell.tdim + r)] =
                    jacobian[static_cast<std::size_t>(r * ls_cell.tdim + a)];
            }
        }
        std::vector<T> rhs(gradient_ref.begin(), gradient_ref.end());
        return solve_graph_dense_small<T>(
            std::move(jt), std::move(rhs), ls_cell.tdim, gradient_phys, tol);
    }

    gradient_phys.assign(static_cast<std::size_t>(ls_cell.gdim), T(0));
    for (int r = 0; r < ls_cell.gdim; ++r)
    {
        for (int a = 0; a < ls_cell.tdim; ++a)
        {
            gradient_phys[static_cast<std::size_t>(r)] +=
                jacobian[static_cast<std::size_t>(r * ls_cell.tdim + a)]
              * gradient_ref[static_cast<std::size_t>(a)];
        }
    }
    return geom::norm<T>(
               std::span<const T>(gradient_phys.data(), gradient_phys.size()))
           > tol;
}

template <std::floating_point T, std::integral I>
int best_orthogonal_surface_edge_refinement_id(
    const AdaptCell<T>& adapt_cell,
    int local_zero_entity_id,
    const LevelSetCell<T, I>& ls_cell,
    std::span<const int> ordered_vertices,
    const ReadyCellGraphOptions<T>& graph_options,
    std::vector<T>* selected_midpoint_ref = nullptr)
{
    const int zdim =
        adapt_cell.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    const int zid =
        adapt_cell.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim != 2 || zid < 0 || zid >= adapt_cell.n_entities(2))
        return -1;

    std::vector<T> jacobian;
    if (!graph_affine_jacobian<T, I>(ls_cell, jacobian))
        return -1;

    const cell::type zero_face_type =
        adapt_cell.entity_types[2][static_cast<std::size_t>(zid)];
    const auto zero_face_edges = cell::edges(zero_face_type);
    const T tol = graph_options.criteria.tolerance;

    int selected_edge_id = -1;
    T selected_abs_cosine = std::numeric_limits<T>::infinity();
    T selected_length = T(0);
    for (const auto& zero_face_edge : zero_face_edges)
    {
        const int ia = zero_face_edge[0];
        const int ib = zero_face_edge[1];
        if (ia < 0 || ib < 0
            || ia >= static_cast<int>(ordered_vertices.size())
            || ib >= static_cast<int>(ordered_vertices.size()))
        {
            continue;
        }

        const int a = ordered_vertices[static_cast<std::size_t>(ia)];
        const int b = ordered_vertices[static_cast<std::size_t>(ib)];
        const auto pa = point_span<T>(adapt_cell.vertex_coords, a, ls_cell.tdim);
        const auto pb = point_span<T>(adapt_cell.vertex_coords, b, ls_cell.tdim);

        std::vector<T> midpoint(static_cast<std::size_t>(ls_cell.tdim), T(0));
        std::vector<T> tangent_ref(static_cast<std::size_t>(ls_cell.tdim), T(0));
        for (int d = 0; d < ls_cell.tdim; ++d)
        {
            midpoint[static_cast<std::size_t>(d)] = T(0.5) * (
                pa[static_cast<std::size_t>(d)]
              + pb[static_cast<std::size_t>(d)]);
            tangent_ref[static_cast<std::size_t>(d)] =
                pb[static_cast<std::size_t>(d)]
              - pa[static_cast<std::size_t>(d)];
        }

        std::vector<T> tangent_phys(static_cast<std::size_t>(ls_cell.gdim), T(0));
        for (int r = 0; r < ls_cell.gdim; ++r)
        {
            for (int d = 0; d < ls_cell.tdim; ++d)
            {
                tangent_phys[static_cast<std::size_t>(r)] +=
                    jacobian[static_cast<std::size_t>(r * ls_cell.tdim + d)]
                  * tangent_ref[static_cast<std::size_t>(d)];
            }
        }

        std::vector<T> gradient_phys;
        if (!graph_physical_level_set_gradient<T, I>(
                ls_cell,
                std::span<const T>(midpoint.data(), midpoint.size()),
                gradient_phys,
                tol))
        {
            continue;
        }

        const T tangent_norm = geom::norm<T>(
            std::span<const T>(tangent_phys.data(), tangent_phys.size()));
        const T gradient_norm = geom::norm<T>(
            std::span<const T>(gradient_phys.data(), gradient_phys.size()));
        if (tangent_norm <= tol || gradient_norm <= tol)
            continue;

        const T abs_cosine = std::fabs(geom::dot<T>(
            std::span<const T>(gradient_phys.data(), gradient_phys.size()),
            std::span<const T>(tangent_phys.data(), tangent_phys.size()))
            / (gradient_norm * tangent_norm));
        const int edge_id =
            graph_refinement_edge_from_failed_face_segment<T>(adapt_cell, a, b);
        if (edge_id < 0)
            continue;

        if (selected_edge_id < 0
            || abs_cosine < selected_abs_cosine
            || (std::fabs(abs_cosine - selected_abs_cosine) <= tol
                && tangent_norm > selected_length))
        {
            selected_edge_id = edge_id;
            selected_abs_cosine = abs_cosine;
            selected_length = tangent_norm;
            if (selected_midpoint_ref != nullptr)
                *selected_midpoint_ref = midpoint;
        }
    }

    return selected_edge_id;
}

template <std::floating_point T, std::integral I>
int best_midpoint_residual_surface_edge_refinement_id(
    const AdaptCell<T>& adapt_cell,
    int local_zero_entity_id,
    const LevelSetCell<T, I>& ls_cell,
    std::span<const int> ordered_vertices,
    const ReadyCellGraphOptions<T>& graph_options,
    std::vector<T>* selected_midpoint_ref = nullptr)
{
    const int zdim =
        adapt_cell.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    const int zid =
        adapt_cell.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim != 2 || zid < 0 || zid >= adapt_cell.n_entities(2))
        return -1;

    std::vector<T> jacobian;
    if (!graph_affine_jacobian<T, I>(ls_cell, jacobian))
        return -1;

    const cell::type zero_face_type =
        adapt_cell.entity_types[2][static_cast<std::size_t>(zid)];
    const auto zero_face_edges = cell::edges(zero_face_type);
    const T tol = graph_options.criteria.tolerance;

    int selected_edge_id = -1;
    T selected_score = -std::numeric_limits<T>::infinity();
    T selected_length = T(0);
    for (const auto& zero_face_edge : zero_face_edges)
    {
        const int ia = zero_face_edge[0];
        const int ib = zero_face_edge[1];
        if (ia < 0 || ib < 0
            || ia >= static_cast<int>(ordered_vertices.size())
            || ib >= static_cast<int>(ordered_vertices.size()))
        {
            continue;
        }

        const int a = ordered_vertices[static_cast<std::size_t>(ia)];
        const int b = ordered_vertices[static_cast<std::size_t>(ib)];
        const auto pa = point_span<T>(adapt_cell.vertex_coords, a, ls_cell.tdim);
        const auto pb = point_span<T>(adapt_cell.vertex_coords, b, ls_cell.tdim);

        std::vector<T> midpoint(static_cast<std::size_t>(ls_cell.tdim), T(0));
        std::vector<T> tangent_ref(static_cast<std::size_t>(ls_cell.tdim), T(0));
        for (int d = 0; d < ls_cell.tdim; ++d)
        {
            midpoint[static_cast<std::size_t>(d)] = T(0.5) * (
                pa[static_cast<std::size_t>(d)]
              + pb[static_cast<std::size_t>(d)]);
            tangent_ref[static_cast<std::size_t>(d)] =
                pb[static_cast<std::size_t>(d)]
              - pa[static_cast<std::size_t>(d)];
        }

        const auto tangent_phys = graph_push_forward_vector<T>(
            jacobian, ls_cell.gdim, ls_cell.tdim,
            std::span<const T>(tangent_ref.data(), tangent_ref.size()));
        const T tangent_norm = geom::norm<T>(
            std::span<const T>(tangent_phys.data(), tangent_phys.size()));
        if (tangent_norm <= tol)
            continue;

        std::vector<T> gradient_phys;
        T gradient_norm = T(1);
        if (graph_physical_level_set_gradient<T, I>(
                ls_cell,
                std::span<const T>(midpoint.data(), midpoint.size()),
                gradient_phys,
                tol))
        {
            gradient_norm = std::max(
                geom::norm<T>(
                    std::span<const T>(
                        gradient_phys.data(), gradient_phys.size())),
                tol);
        }

        const T score = std::fabs(ls_cell.value(
                            std::span<const T>(
                                midpoint.data(), midpoint.size())))
                      / gradient_norm;
        const int edge_id =
            graph_refinement_edge_from_failed_face_segment<T>(adapt_cell, a, b);
        if (edge_id < 0)
            continue;

        if (selected_edge_id < 0
            || score > selected_score
            || (std::fabs(score - selected_score) <= tol
                && tangent_norm > selected_length))
        {
            selected_edge_id = edge_id;
            selected_score = score;
            selected_length = tangent_norm;
            if (selected_midpoint_ref != nullptr)
                *selected_midpoint_ref = midpoint;
        }
    }

    return selected_edge_id;
}

template <std::floating_point T, std::integral I>
int best_normal_variation_surface_edge_refinement_id(
    const AdaptCell<T>& adapt_cell,
    int local_zero_entity_id,
    const LevelSetCell<T, I>& ls_cell,
    std::span<const int> ordered_vertices,
    const ReadyCellGraphOptions<T>& graph_options,
    std::vector<T>* selected_midpoint_ref = nullptr)
{
    const int zdim =
        adapt_cell.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    const int zid =
        adapt_cell.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim != 2 || zid < 0 || zid >= adapt_cell.n_entities(2))
        return -1;

    std::vector<T> jacobian;
    if (!graph_affine_jacobian<T, I>(ls_cell, jacobian))
        return -1;

    const cell::type zero_face_type =
        adapt_cell.entity_types[2][static_cast<std::size_t>(zid)];
    const auto zero_face_edges = cell::edges(zero_face_type);
    const T tol = graph_options.criteria.tolerance;

    int selected_edge_id = -1;
    T selected_score = -std::numeric_limits<T>::infinity();
    T selected_length = T(0);
    for (const auto& zero_face_edge : zero_face_edges)
    {
        const int ia = zero_face_edge[0];
        const int ib = zero_face_edge[1];
        if (ia < 0 || ib < 0
            || ia >= static_cast<int>(ordered_vertices.size())
            || ib >= static_cast<int>(ordered_vertices.size()))
        {
            continue;
        }

        const int a = ordered_vertices[static_cast<std::size_t>(ia)];
        const int b = ordered_vertices[static_cast<std::size_t>(ib)];
        const auto pa = point_span<T>(adapt_cell.vertex_coords, a, ls_cell.tdim);
        const auto pb = point_span<T>(adapt_cell.vertex_coords, b, ls_cell.tdim);

        std::vector<T> midpoint(static_cast<std::size_t>(ls_cell.tdim), T(0));
        std::vector<T> tangent_ref(static_cast<std::size_t>(ls_cell.tdim), T(0));
        for (int d = 0; d < ls_cell.tdim; ++d)
        {
            midpoint[static_cast<std::size_t>(d)] = T(0.5) * (
                pa[static_cast<std::size_t>(d)]
              + pb[static_cast<std::size_t>(d)]);
            tangent_ref[static_cast<std::size_t>(d)] =
                pb[static_cast<std::size_t>(d)]
              - pa[static_cast<std::size_t>(d)];
        }

        const auto tangent_phys = graph_push_forward_vector<T>(
            jacobian, ls_cell.gdim, ls_cell.tdim,
            std::span<const T>(tangent_ref.data(), tangent_ref.size()));
        const T tangent_norm = geom::norm<T>(
            std::span<const T>(tangent_phys.data(), tangent_phys.size()));
        if (tangent_norm <= tol)
            continue;

        std::vector<T> ga;
        std::vector<T> gb;
        if (!graph_physical_level_set_gradient<T, I>(
                ls_cell, pa, ga, tol)
            || !graph_physical_level_set_gradient<T, I>(
                ls_cell, pb, gb, tol))
        {
            continue;
        }
        const T ga_norm = geom::norm<T>(
            std::span<const T>(ga.data(), ga.size()));
        const T gb_norm = geom::norm<T>(
            std::span<const T>(gb.data(), gb.size()));
        if (ga_norm <= tol || gb_norm <= tol)
            continue;

        const T abs_cosine = std::fabs(geom::dot<T>(
            std::span<const T>(ga.data(), ga.size()),
            std::span<const T>(gb.data(), gb.size()))
            / (ga_norm * gb_norm));
        const T score = (T(1) - std::clamp(abs_cosine, T(0), T(1)))
                      * tangent_norm;
        const int edge_id =
            graph_refinement_edge_from_failed_face_segment<T>(adapt_cell, a, b);
        if (edge_id < 0)
            continue;

        if (selected_edge_id < 0
            || score > selected_score
            || (std::fabs(score - selected_score) <= tol
                && tangent_norm > selected_length))
        {
            selected_edge_id = edge_id;
            selected_score = score;
            selected_length = tangent_norm;
            if (selected_midpoint_ref != nullptr)
                *selected_midpoint_ref = midpoint;
        }
    }

    return selected_edge_id;
}

template <std::floating_point T>
std::vector<T> graph_push_forward_vector(const std::vector<T>& jacobian,
                                         int gdim,
                                         int tdim,
                                         std::span<const T> vector_ref)
{
    std::vector<T> out(static_cast<std::size_t>(gdim), T(0));
    if (static_cast<int>(vector_ref.size()) != tdim
        || static_cast<int>(jacobian.size()) != gdim * tdim)
    {
        return out;
    }
    for (int r = 0; r < gdim; ++r)
    {
        for (int a = 0; a < tdim; ++a)
        {
            out[static_cast<std::size_t>(r)] +=
                jacobian[static_cast<std::size_t>(r * tdim + a)]
              * vector_ref[static_cast<std::size_t>(a)];
        }
    }
    return out;
}

template <std::floating_point T>
std::pair<T, T> graph_legendre_value_and_derivative(int order, T x)
{
    if (order == 0)
        return {T(1), T(0)};
    if (order == 1)
        return {x, T(1)};

    T pm2 = T(1);
    T pm1 = x;
    for (int n = 2; n <= order; ++n)
    {
        const T p = ((T(2 * n - 1) * x * pm1) - T(n - 1) * pm2) / T(n);
        pm2 = pm1;
        pm1 = p;
    }

    const T denom = T(1) - x * x;
    if (std::abs(denom) <= T(64) * std::numeric_limits<T>::epsilon())
        return {pm1, T(0)};
    return {pm1, T(order) * (pm2 - x * pm1) / denom};
}

template <std::floating_point T>
std::vector<T> graph_gll_parameters(int order)
{
    order = std::max(order, 1);
    std::vector<T> params(static_cast<std::size_t>(order + 1), T(0));
    params.front() = T(0);
    params.back() = T(1);
    if (order == 1)
        return params;

    const T pi = std::acos(T(-1));
    const T eps = T(128) * std::numeric_limits<T>::epsilon();
    for (int i = 1; i < order; ++i)
    {
        T x = -std::cos(pi * T(i) / T(order));
        for (int iter = 0; iter < 32; ++iter)
        {
            const auto [p, dp] = graph_legendre_value_and_derivative<T>(order, x);
            const T denom = T(1) - x * x;
            if (std::abs(denom) <= eps)
                break;
            const T d2p = (T(2) * x * dp - T(order * (order + 1)) * p) / denom;
            if (std::abs(d2p) <= eps)
                break;
            const T step = dp / d2p;
            x -= step;
            x = std::clamp(x, -T(1) + eps, T(1) - eps);
            if (std::abs(step) <= eps)
                break;
        }
        params[static_cast<std::size_t>(i)] = T(0.5) * (x + T(1));
    }
    return params;
}

template <std::floating_point T>
std::vector<T> graph_interpolation_parameters(int order,
                                              curving::NodeFamily family)
{
    order = std::max(order, 1);
    if (family == curving::NodeFamily::gll)
        return graph_gll_parameters<T>(order);

    std::vector<T> params(static_cast<std::size_t>(order + 1), T(0));
    for (int i = 0; i <= order; ++i)
        params[static_cast<std::size_t>(i)] = T(i) / T(order);
    return params;
}

template <std::floating_point T>
T graph_lagrange_basis_1d(int i, std::span<const T> params, T x)
{
    T value = T(1);
    const T xi = params[static_cast<std::size_t>(i)];
    for (int j = 0; j < static_cast<int>(params.size()); ++j)
    {
        if (j == i)
            continue;
        value *= (x - params[static_cast<std::size_t>(j)])
               / (xi - params[static_cast<std::size_t>(j)]);
    }
    return value;
}

template <std::floating_point T>
T graph_warp_factor(int order, T r)
{
    if (order <= 1)
        return T(0);

    std::vector<T> equispaced(static_cast<std::size_t>(order + 1), T(0));
    auto gll = graph_gll_parameters<T>(order);
    for (int i = 0; i <= order; ++i)
    {
        equispaced[static_cast<std::size_t>(i)] =
            -T(1) + T(2 * i) / T(order);
        gll[static_cast<std::size_t>(i)] =
            T(2) * gll[static_cast<std::size_t>(i)] - T(1);
    }

    T warp = T(0);
    for (int i = 0; i <= order; ++i)
    {
        const T Li = graph_lagrange_basis_1d<T>(
            i, std::span<const T>(equispaced.data(), equispaced.size()), r);
        warp += Li * (gll[static_cast<std::size_t>(i)]
                    - equispaced[static_cast<std::size_t>(i)]);
    }

    const T edge_factor = T(1) - r * r;
    if (std::abs(edge_factor) > T(64) * std::numeric_limits<T>::epsilon())
        warp /= edge_factor;
    return warp;
}

template <std::floating_point T>
std::array<T, 3> graph_equilateral_to_reference_barycentric(T x, T y)
{
    const T sqrt3 = std::sqrt(T(3));
    std::array<T, 3> w = {};
    w[2] = (sqrt3 * y + T(1)) / T(3);
    w[1] = (x + T(1) - w[2]) / T(2);
    w[0] = T(1) - w[1] - w[2];

    T sum = T(0);
    for (T& wi : w)
    {
        if (std::abs(wi) < T(256) * std::numeric_limits<T>::epsilon())
            wi = T(0);
        if (std::abs(wi - T(1)) < T(256) * std::numeric_limits<T>::epsilon())
            wi = T(1);
        wi = std::clamp(wi, T(0), T(1));
        sum += wi;
    }
    if (sum > T(0))
        for (T& wi : w)
            wi /= sum;
    return w;
}

template <std::floating_point T>
std::vector<std::array<T, 3>> graph_triangle_barycentric_nodes(
    int order,
    curving::NodeFamily family)
{
    order = std::max(order, 1);
    std::vector<std::array<T, 3>> nodes;
    nodes.reserve(static_cast<std::size_t>((order + 1) * (order + 2) / 2));

    if (family != curving::NodeFamily::gll)
    {
        for (int j = 0; j <= order; ++j)
        {
            for (int i = 0; i <= order - j; ++i)
            {
                const T u = T(i) / T(order);
                const T v = T(j) / T(order);
                nodes.push_back({T(1) - u - v, u, v});
            }
        }
        return nodes;
    }

    constexpr std::array<T, 16> alpha_opt = {
        T(0.0), T(0.0), T(1.4152), T(0.1001),
        T(0.2751), T(0.9800), T(1.0999), T(1.2832),
        T(1.3648), T(1.4773), T(1.4959), T(1.5743),
        T(1.5770), T(1.6223), T(1.6258), T(1.6530)};
    const T alpha = (order < static_cast<int>(alpha_opt.size()))
                      ? alpha_opt[static_cast<std::size_t>(order)]
                      : T(5) / T(3);
    const T sqrt3 = std::sqrt(T(3));
    const T cos120 = -T(0.5);
    const T sin120 = sqrt3 / T(2);
    const T cos240 = -T(0.5);
    const T sin240 = -sqrt3 / T(2);

    for (int j = 0; j <= order; ++j)
    {
        for (int i = 0; i <= order - j; ++i)
        {
            const T u = T(i) / T(order);
            const T v = T(j) / T(order);
            const T lambda0 = T(1) - u - v;
            const T lambda1 = u;
            const T lambda2 = v;

            T x = -lambda0 + lambda1;
            T y = (-lambda0 - lambda1 + T(2) * lambda2) / sqrt3;

            const T L1 = lambda2;
            const T L2 = lambda0;
            const T L3 = lambda1;
            const T warp1 = T(4) * L2 * L3 * graph_warp_factor<T>(order, L3 - L2)
                          * (T(1) + (alpha * L1) * (alpha * L1));
            const T warp2 = T(4) * L1 * L3 * graph_warp_factor<T>(order, L1 - L3)
                          * (T(1) + (alpha * L2) * (alpha * L2));
            const T warp3 = T(4) * L1 * L2 * graph_warp_factor<T>(order, L2 - L1)
                          * (T(1) + (alpha * L3) * (alpha * L3));

            x += warp1 + cos120 * warp2 + cos240 * warp3;
            y += sin120 * warp2 + sin240 * warp3;
            nodes.push_back(graph_equilateral_to_reference_barycentric<T>(x, y));
        }
    }
    return nodes;
}

template <std::floating_point T>
std::vector<T> graph_triangle_monomials(int order, std::span<const T> bary)
{
    const T u = bary[1];
    const T v = bary[2];
    std::vector<T> values;
    values.reserve(static_cast<std::size_t>((order + 1) * (order + 2) / 2));
    for (int total = 0; total <= order; ++total)
    {
        for (int j = 0; j <= total; ++j)
        {
            const int i = total - j;
            values.push_back(std::pow(u, i) * std::pow(v, j));
        }
    }
    return values;
}

template <std::floating_point T>
std::vector<T> graph_triangle_lagrange_basis(
    int order,
    curving::NodeFamily family,
    std::span<const T> bary)
{
    const auto nodes = graph_triangle_barycentric_nodes<T>(order, family);
    const int n = static_cast<int>(nodes.size());
    std::vector<T> matrix(static_cast<std::size_t>(n * n), T(0));
    for (int row = 0; row < n; ++row)
    {
        const auto mono = graph_triangle_monomials<T>(
            order,
            std::span<const T>(nodes[static_cast<std::size_t>(row)].data(), 3));
        for (int col = 0; col < n; ++col)
            matrix[static_cast<std::size_t>(col * n + row)] =
                mono[static_cast<std::size_t>(col)];
    }

    const auto rhs = graph_triangle_monomials<T>(order, bary);
    std::vector<T> basis;
    if (!solve_graph_dense_small<T>(std::move(matrix), rhs, n, basis,
                                    T(256) * std::numeric_limits<T>::epsilon()))
    {
        return {};
    }
    return basis;
}

template <std::floating_point T>
std::vector<T> graph_eval_curved_zero_face_ref(
    const curving::CurvedZeroEntityState<T>& state,
    cell::type entity_type,
    int order,
    curving::NodeFamily family,
    std::span<const T> xi)
{
    const int nodes_per_face =
        (entity_type == cell::type::quadrilateral)
            ? (order + 1) * (order + 1)
            : (order + 1) * (order + 2) / 2;
    if (nodes_per_face <= 0
        || state.ref_nodes.size() % static_cast<std::size_t>(nodes_per_face) != 0)
    {
        return {};
    }
    const int tdim = static_cast<int>(
        state.ref_nodes.size() / static_cast<std::size_t>(nodes_per_face));
    std::vector<T> out(static_cast<std::size_t>(tdim), T(0));

    if (entity_type == cell::type::quadrilateral)
    {
        const auto params = graph_interpolation_parameters<T>(order, family);
        const T u = xi[0];
        const T v = xi[1];
        int node = 0;
        for (int j = 0; j <= order; ++j)
        {
            const T Lj = graph_lagrange_basis_1d<T>(
                j, std::span<const T>(params.data(), params.size()), v);
            for (int i = 0; i <= order; ++i)
            {
                const T Li = graph_lagrange_basis_1d<T>(
                    i, std::span<const T>(params.data(), params.size()), u);
                const T L = Li * Lj;
                for (int d = 0; d < tdim; ++d)
                    out[static_cast<std::size_t>(d)] += L * state.ref_nodes[
                        static_cast<std::size_t>(node * tdim + d)];
                ++node;
            }
        }
        return out;
    }

    if (entity_type == cell::type::triangle)
    {
        const std::array<T, 3> bary = {T(1) - xi[0] - xi[1], xi[0], xi[1]};
        const auto basis = graph_triangle_lagrange_basis<T>(
            order, family, std::span<const T>(bary.data(), bary.size()));
        if (basis.empty())
            return {};
        for (int node = 0; node < static_cast<int>(basis.size()); ++node)
        {
            const T L = basis[static_cast<std::size_t>(node)];
            for (int d = 0; d < tdim; ++d)
                out[static_cast<std::size_t>(d)] += L * state.ref_nodes[
                    static_cast<std::size_t>(node * tdim + d)];
        }
        return out;
    }

    return {};
}

template <std::floating_point T>
std::vector<T> graph_eval_straight_zero_face_ref(
    const AdaptCell<T>& adapt_cell,
    std::span<const int> vertices,
    cell::type entity_type,
    std::span<const T> xi)
{
    std::vector<T> out(static_cast<std::size_t>(adapt_cell.tdim), T(0));
    std::array<T, 4> weights = {};
    int nweights = 0;
    if (entity_type == cell::type::triangle)
    {
        weights = {T(1) - xi[0] - xi[1], xi[0], xi[1], T(0)};
        nweights = 3;
    }
    else if (entity_type == cell::type::quadrilateral)
    {
        const T u = xi[0];
        const T v = xi[1];
        weights = {(T(1) - u) * (T(1) - v),
                   u * (T(1) - v),
                   (T(1) - u) * v,
                   u * v};
        nweights = 4;
    }
    else
    {
        return {};
    }

    if (static_cast<int>(vertices.size()) < nweights)
        return {};
    for (int i = 0; i < nweights; ++i)
    {
        const int vertex = vertices[static_cast<std::size_t>(i)];
        for (int d = 0; d < adapt_cell.tdim; ++d)
        {
            out[static_cast<std::size_t>(d)] +=
                weights[static_cast<std::size_t>(i)]
              * adapt_cell.vertex_coords[
                    static_cast<std::size_t>(vertex * adapt_cell.tdim + d)];
        }
    }
    return out;
}

template <std::floating_point T>
T graph_signed_surface_jacobian_ratio(std::span<const T> host_du_phys,
                                      std::span<const T> host_dv_phys,
                                      std::span<const T> corr_du_phys,
                                      std::span<const T> corr_dv_phys,
                                      T tol)
{
    if (host_du_phys.size() != 3 || host_dv_phys.size() != 3
        || corr_du_phys.size() != 3
        || corr_dv_phys.size() != 3)
    {
        return std::numeric_limits<T>::quiet_NaN();
    }

    const auto host_cross = geom::cross<T>(
        host_du_phys, host_dv_phys);
    const T denom =
        host_cross[0] * host_cross[0]
      + host_cross[1] * host_cross[1]
      + host_cross[2] * host_cross[2];
    if (denom <= tol * tol)
        return std::numeric_limits<T>::quiet_NaN();

    const auto corr_cross = geom::cross<T>(corr_du_phys, corr_dv_phys);
    const T numer =
        corr_cross[0] * host_cross[0]
      + corr_cross[1] * host_cross[1]
      + corr_cross[2] * host_cross[2];
    return numer / denom;
}

template <std::floating_point T, std::integral I>
bool check_zero_face_surface_jacobian(
    const AdaptCell<T>& adapt_cell,
    int local_zero_entity_id,
    const LevelSetCell<T, I>& ls_cell,
    std::span<const int> ordered_vertices,
    const ReadyCellGraphOptions<T>& graph_options,
    ReadyCellGraphDiagnostics<T>& diagnostics,
    int source_cell_id)
{
    graph_criteria::FaceQualityReport<T> report;
    report.minimum_surface_jacobian_ratio =
        std::numeric_limits<T>::infinity();

    const int dim = ls_cell.tdim;
    if (dim != 3 || ls_cell.gdim != 3)
    {
        mark_graph_failure<T>(
            diagnostics,
            source_cell_id,
            graph_criteria::FailureReason::invalid_input);
        return false;
    }

    const int zdim =
        adapt_cell.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    const int zid =
        adapt_cell.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim != 2 || zid < 0 || zid >= adapt_cell.n_entities(2))
    {
        mark_graph_failure<T>(
            diagnostics,
            source_cell_id,
            graph_criteria::FailureReason::invalid_input);
        return false;
    }

    const cell::type entity_type =
        adapt_cell.entity_types[2][static_cast<std::size_t>(zid)];
    if (entity_type != cell::type::triangle
        && entity_type != cell::type::quadrilateral)
    {
        mark_graph_failure<T>(
            diagnostics,
            source_cell_id,
            graph_criteria::FailureReason::invalid_input);
        return false;
    }

    const auto curving_options = graph_curving_options<T>(graph_options);
    curving::CurvingData<T, I> curving_data;
    const std::array<I, 1> parent_cell_ids = {static_cast<I>(ls_cell.cell_id)};
    const std::array<AdaptCell<T>, 1> adapt_cells = {adapt_cell};
    const std::array<LevelSetCell<T, I>, 1> level_set_cells = {ls_cell};
    const std::array<int, 2> ls_offsets = {0, 1};
    const auto& state = curving::ensure_curved<T, I>(
        curving_data,
        std::span<const I>(parent_cell_ids.data(), parent_cell_ids.size()),
        std::span<const AdaptCell<T>>(adapt_cells.data(), adapt_cells.size()),
        std::span<const LevelSetCell<T, I>>(
            level_set_cells.data(), level_set_cells.size()),
        std::span<const int>(ls_offsets.data(), ls_offsets.size()),
        0,
        local_zero_entity_id,
        curving_options);
    if (state.status != curving::CurvingStatus::curved)
    {
        mark_graph_failure<T>(
            diagnostics,
            source_cell_id,
            graph_criteria::FailureReason::invalid_input);
        return false;
    }

    std::vector<T> jacobian;
    if (!graph_affine_jacobian<T, I>(ls_cell, jacobian))
    {
        mark_graph_failure<T>(
            diagnostics,
            source_cell_id,
            graph_criteria::FailureReason::invalid_input);
        return false;
    }

    std::vector<std::array<T, 2>> samples;
    if (entity_type == cell::type::triangle)
    {
        samples = {{{T(1) / T(3), T(1) / T(3)},
                    {T(0.2), T(0.2)},
                    {T(0.6), T(0.2)},
                    {T(0.2), T(0.6)}}};
    }
    else
    {
        samples = {{{T(0.5), T(0.5)},
                    {T(0.25), T(0.25)},
                    {T(0.75), T(0.25)},
                    {T(0.25), T(0.75)},
                    {T(0.75), T(0.75)}}};
    }

    const T h = T(1.0e-5);
    auto eval_curved = [&](std::span<const T> xi)
    {
        return graph_eval_curved_zero_face_ref<T>(
            state, entity_type, curving_options.geometry_order,
            curving_options.node_family, xi);
    };
    auto eval_straight = [&](std::span<const T> xi)
    {
        return graph_eval_straight_zero_face_ref<T>(
            adapt_cell, ordered_vertices, entity_type, xi);
    };

    for (int sample_id = 0; sample_id < static_cast<int>(samples.size()); ++sample_id)
    {
        const auto xi = samples[static_cast<std::size_t>(sample_id)];
        const std::array<T, 2> xi_u_plus = {xi[0] + h, xi[1]};
        const std::array<T, 2> xi_u_minus = {xi[0] - h, xi[1]};
        const std::array<T, 2> xi_v_plus = {xi[0], xi[1] + h};
        const std::array<T, 2> xi_v_minus = {xi[0], xi[1] - h};

        const auto curved_u_plus = eval_curved(
            std::span<const T>(xi_u_plus.data(), xi_u_plus.size()));
        const auto curved_u_minus = eval_curved(
            std::span<const T>(xi_u_minus.data(), xi_u_minus.size()));
        const auto curved_v_plus = eval_curved(
            std::span<const T>(xi_v_plus.data(), xi_v_plus.size()));
        const auto curved_v_minus = eval_curved(
            std::span<const T>(xi_v_minus.data(), xi_v_minus.size()));
        const auto straight_u_plus = eval_straight(
            std::span<const T>(xi_u_plus.data(), xi_u_plus.size()));
        const auto straight_u_minus = eval_straight(
            std::span<const T>(xi_u_minus.data(), xi_u_minus.size()));
        const auto straight_v_plus = eval_straight(
            std::span<const T>(xi_v_plus.data(), xi_v_plus.size()));
        const auto straight_v_minus = eval_straight(
            std::span<const T>(xi_v_minus.data(), xi_v_minus.size()));

        if (curved_u_plus.empty() || curved_u_minus.empty()
            || curved_v_plus.empty() || curved_v_minus.empty()
            || straight_u_plus.empty() || straight_u_minus.empty()
            || straight_v_plus.empty() || straight_v_minus.empty())
        {
            mark_graph_failure<T>(
                diagnostics,
                source_cell_id,
                graph_criteria::FailureReason::invalid_input);
            return false;
        }

        std::vector<T> curved_du(static_cast<std::size_t>(dim), T(0));
        std::vector<T> curved_dv(static_cast<std::size_t>(dim), T(0));
        std::vector<T> straight_du(static_cast<std::size_t>(dim), T(0));
        std::vector<T> straight_dv(static_cast<std::size_t>(dim), T(0));
        for (int d = 0; d < dim; ++d)
        {
            curved_du[static_cast<std::size_t>(d)] =
                (curved_u_plus[static_cast<std::size_t>(d)]
               - curved_u_minus[static_cast<std::size_t>(d)]) / (T(2) * h);
            curved_dv[static_cast<std::size_t>(d)] =
                (curved_v_plus[static_cast<std::size_t>(d)]
               - curved_v_minus[static_cast<std::size_t>(d)]) / (T(2) * h);
            straight_du[static_cast<std::size_t>(d)] =
                (straight_u_plus[static_cast<std::size_t>(d)]
               - straight_u_minus[static_cast<std::size_t>(d)]) / (T(2) * h);
            straight_dv[static_cast<std::size_t>(d)] =
                (straight_v_plus[static_cast<std::size_t>(d)]
               - straight_v_minus[static_cast<std::size_t>(d)]) / (T(2) * h);
        }

        const auto curved_du_phys = graph_push_forward_vector<T>(
            jacobian, ls_cell.gdim, ls_cell.tdim,
            std::span<const T>(curved_du.data(), curved_du.size()));
        const auto curved_dv_phys = graph_push_forward_vector<T>(
            jacobian, ls_cell.gdim, ls_cell.tdim,
            std::span<const T>(curved_dv.data(), curved_dv.size()));
        const auto straight_du_phys = graph_push_forward_vector<T>(
            jacobian, ls_cell.gdim, ls_cell.tdim,
            std::span<const T>(straight_du.data(), straight_du.size()));
        const auto straight_dv_phys = graph_push_forward_vector<T>(
            jacobian, ls_cell.gdim, ls_cell.tdim,
            std::span<const T>(straight_dv.data(), straight_dv.size()));

        const T ratio = graph_signed_surface_jacobian_ratio<T>(
            std::span<const T>(straight_du_phys.data(), straight_du_phys.size()),
            std::span<const T>(straight_dv_phys.data(), straight_dv_phys.size()),
            std::span<const T>(curved_du_phys.data(), curved_du_phys.size()),
            std::span<const T>(curved_dv_phys.data(), curved_dv_phys.size()),
            graph_options.criteria.tolerance);
        if (!std::isfinite(ratio))
        {
            mark_graph_failure<T>(
                diagnostics,
                source_cell_id,
                graph_criteria::FailureReason::face_degenerate);
            return false;
        }

        report.minimum_surface_jacobian_ratio =
            std::min(report.minimum_surface_jacobian_ratio, ratio);
        if (ratio <= graph_options.criteria.min_surface_jacobian_ratio)
        {
            report.accepted = false;
            report.failure_reason =
                graph_criteria::FailureReason::surface_jacobian_not_positive;
            report.failed_triangle_index = sample_id;
            report.failed_surface_jacobian_ratio = ratio;
            diagnostics.min_face_area_ratio =
                std::min(diagnostics.min_face_area_ratio,
                         report.minimum_surface_jacobian_ratio);
            observe_failed_face_orientation<T>(diagnostics, report);
            int requested_edge_id = -1;
            if (entity_type == cell::type::triangle
                && ordered_vertices.size() >= 3)
            {
                const std::array<T, 3> bary = {T(1) - xi[0] - xi[1], xi[0], xi[1]};
                int closest = 0;
                if (bary[1] < bary[closest])
                    closest = 1;
                if (bary[2] < bary[closest])
                    closest = 2;
                const int a = ordered_vertices[
                    static_cast<std::size_t>((closest + 1) % 3)];
                const int b = ordered_vertices[
                    static_cast<std::size_t>((closest + 2) % 3)];
                requested_edge_id =
                    graph_refinement_edge_from_failed_face_segment<T>(
                        adapt_cell, a, b);
            }
            else if (entity_type == cell::type::quadrilateral
                     && ordered_vertices.size() >= 4)
            {
                std::array<T, 4> dist = {xi[1], T(1) - xi[0], T(1) - xi[1], xi[0]};
                int side = 0;
                for (int i = 1; i < 4; ++i)
                    if (dist[static_cast<std::size_t>(i)]
                        < dist[static_cast<std::size_t>(side)])
                        side = i;
                const std::array<std::array<int, 2>, 4> side_vertices = {{
                    {{0, 1}}, {{1, 3}}, {{2, 3}}, {{0, 2}}}};
                const int a = ordered_vertices[static_cast<std::size_t>(
                    side_vertices[static_cast<std::size_t>(side)][0])];
                const int b = ordered_vertices[static_cast<std::size_t>(
                    side_vertices[static_cast<std::size_t>(side)][1])];
                requested_edge_id =
                    graph_refinement_edge_from_failed_face_segment<T>(
                        adapt_cell, a, b);
            }
            mark_graph_refinement_request<T>(
                diagnostics, 1, requested_edge_id);
            mark_graph_failure<T>(
                diagnostics, source_cell_id, report.failure_reason);
            return false;
        }
    }

    diagnostics.min_face_area_ratio =
        std::min(diagnostics.min_face_area_ratio,
                 report.minimum_surface_jacobian_ratio);
    return true;
}

template <std::floating_point T, std::integral I>
bool check_zero_face_graph(const AdaptCell<T>& adapt_cell,
                           int local_zero_entity_id,
                           const LevelSetCell<T, I>& ls_cell,
                           std::span<const int> ordered_vertices,
                           const std::vector<T>& vertex_coords,
                           const ReadyCellGraphOptions<T>& graph_options,
                           ReadyCellGraphDiagnostics<T>& diagnostics,
                           int source_cell_id)
{
    ++diagnostics.checked_faces;

    std::vector<T> face_vertices;
    face_vertices.reserve(ordered_vertices.size() * static_cast<std::size_t>(ls_cell.tdim));
    for (const int v : ordered_vertices)
    {
        auto p = point_span<T>(vertex_coords, v, ls_cell.tdim);
        face_vertices.insert(face_vertices.end(), p.begin(), p.end());
    }

    graph_criteria::HostFrame<T> host;
    if (!build_face_host_frame<T>(
            face_vertices, host, graph_options.criteria.tolerance))
    {
        mark_graph_failure<T>(
            diagnostics,
            source_cell_id,
            graph_criteria::FailureReason::invalid_host_frame);
        return false;
    }

    const int nverts = static_cast<int>(ordered_vertices.size());
    const int zdim =
        adapt_cell.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    const int zid =
        adapt_cell.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim != 2 || zid < 0 || zid >= adapt_cell.n_entities(2))
    {
        mark_graph_failure<T>(
            diagnostics,
            source_cell_id,
            graph_criteria::FailureReason::invalid_input);
        return false;
    }

    const cell::type zero_face_type =
        adapt_cell.entity_types[2][static_cast<std::size_t>(zid)];
    const auto zero_face_edges = cell::edges(zero_face_type);
    const geom::ParentEntity cell_interior{ls_cell.tdim, -1};
    int longest_face_edge_id = -1;
    T longest_face_edge_length = T(0);
    int face_interior_node_index = 1;
    for (const auto& zero_face_edge : zero_face_edges)
    {
        const int ia = zero_face_edge[0];
        const int ib = zero_face_edge[1];
        if (ia < 0 || ib < 0 || ia >= nverts || ib >= nverts)
        {
            mark_graph_failure<T>(
                diagnostics,
                source_cell_id,
                graph_criteria::FailureReason::invalid_input);
            return false;
        }

        const int a = ordered_vertices[static_cast<std::size_t>(ia)];
        const int b = ordered_vertices[static_cast<std::size_t>(ib)];
        const auto pa = point_span<T>(vertex_coords, a, ls_cell.tdim);
        const auto pb = point_span<T>(vertex_coords, b, ls_cell.tdim);
        const auto edge_delta = geom::subtract<T>(pb, pa);
        const T edge_length = geom::norm<T>(
            std::span<const T>(edge_delta.data(), edge_delta.size()));
        const int adapt_edge_id = find_adapt_edge_by_vertices<T>(
            adapt_cell, a, b);
        if (adapt_edge_id >= 0 && edge_length > longest_face_edge_length)
        {
            longest_face_edge_length = edge_length;
            longest_face_edge_id = adapt_edge_id;
        }

        const std::uint64_t face_zero_mask =
            adapt_cell.zero_entity_zero_mask[
                static_cast<std::size_t>(local_zero_entity_id)];
        const int zero_edge_id =
            find_zero_edge_entity_by_vertices<T>(
                adapt_cell, a, b, face_zero_mask);
        if (zero_edge_id < 0)
        {
            mark_graph_failure<T>(
                diagnostics,
                source_cell_id,
                graph_criteria::FailureReason::invalid_host_frame);
            return false;
        }

        int host_face_id = -1;
        auto boundary_host_normal =
            host_boundary_face_normal_for_zero_face_edge<T>(
                adapt_cell,
                local_zero_entity_id,
                pa,
                pb,
                graph_options.criteria.tolerance,
                &host_face_id);
        if (boundary_host_normal.empty() || host_face_id < 0)
        {
            // This is an interior zero edge, for example the artificial
            // diagonal introduced by triangulating a zero quadrilateral. It
            // still curves as a zero edge, but its admissible host is the
            // uncut subcell volume and graph checks use face-interior rules.
            AdaptCell<T> edge_context = adapt_cell;
            override_zero_edge_host_from_zero_face<T>(
                edge_context,
                zero_edge_id,
                local_zero_entity_id,
                -1);

            const int order = std::max(graph_options.geometry_order, 1);
            for (int k = 1; k < order; ++k)
            {
                const T s = T(k) / T(order);
                std::vector<T> seed(static_cast<std::size_t>(ls_cell.tdim), T(0));
                for (int d = 0; d < ls_cell.tdim; ++d)
                {
                    seed[static_cast<std::size_t>(d)] =
                        (T(1) - s) * pa[static_cast<std::size_t>(d)]
                      + s * pb[static_cast<std::size_t>(d)];
                }

                std::vector<T> corrected;
                if (!project_graph_point<T, I>(
                        edge_context,
                        zero_edge_id,
                        ls_cell,
                        cell_interior,
                        host,
                        std::span<const T>(seed.data(), seed.size()),
                        GraphNodeKind::face_interior,
                        face_interior_node_index++,
                        adapt_edge_id,
                        graph_options,
                        diagnostics,
                        source_cell_id,
                        corrected))
                {
                    return false;
                }
            }
            continue;
        }

        AdaptCell<T> edge_context = adapt_cell;
        override_zero_edge_host_from_zero_face<T>(
            edge_context,
            zero_edge_id,
            local_zero_entity_id,
            host_face_id);
        if (!check_zero_edge_graph<T, I>(
                edge_context,
                zero_edge_id,
                ls_cell,
                pa,
                pb,
                std::span<const T>(
                    boundary_host_normal.data(), boundary_host_normal.size()),
                true,
                adapt_edge_id,
                graph_options,
                diagnostics,
                source_cell_id))
        {
            return false;
        }
    }

    std::vector<T> face_center(static_cast<std::size_t>(ls_cell.tdim), T(0));
    for (const int v : ordered_vertices)
    {
        const auto p = point_span<T>(vertex_coords, v, ls_cell.tdim);
        for (int d = 0; d < ls_cell.tdim; ++d)
            face_center[static_cast<std::size_t>(d)] += p[static_cast<std::size_t>(d)];
    }
    for (T& value : face_center)
        value /= static_cast<T>(nverts);

    std::vector<T> corrected_face_center;
    int face_center_refinement_edge_id = longest_face_edge_id;
    std::vector<T> face_center_refinement_point;
    if (graph_options.refinement_mode
        == GraphRefinementMode::green_orthogonal_surface_edge)
    {
        const int orthogonal_edge_id =
            best_orthogonal_surface_edge_refinement_id<T, I>(
                adapt_cell,
                local_zero_entity_id,
                ls_cell,
                ordered_vertices,
                graph_options,
                &face_center_refinement_point);
        if (orthogonal_edge_id >= 0)
            face_center_refinement_edge_id = orthogonal_edge_id;
    }
    else if (graph_options.refinement_mode
             == GraphRefinementMode::green_midpoint_residual)
    {
        const int residual_edge_id =
            best_midpoint_residual_surface_edge_refinement_id<T, I>(
                adapt_cell,
                local_zero_entity_id,
                ls_cell,
                ordered_vertices,
                graph_options,
                &face_center_refinement_point);
        if (residual_edge_id >= 0)
            face_center_refinement_edge_id = residual_edge_id;
    }
    else if (graph_options.refinement_mode
             == GraphRefinementMode::green_normal_variation)
    {
        const int normal_variation_edge_id =
            best_normal_variation_surface_edge_refinement_id<T, I>(
                adapt_cell,
                local_zero_entity_id,
                ls_cell,
                ordered_vertices,
                graph_options,
                &face_center_refinement_point);
        if (normal_variation_edge_id >= 0)
            face_center_refinement_edge_id = normal_variation_edge_id;
    }
    if (!project_graph_point<T, I>(
            adapt_cell,
            local_zero_entity_id,
            ls_cell,
            cell_interior,
            host,
            std::span<const T>(face_center.data(), face_center.size()),
            GraphNodeKind::face_interior,
            0,
            face_center_refinement_edge_id,
            graph_options,
            diagnostics,
            source_cell_id,
            corrected_face_center))
    {
        mark_graph_refinement_point<T>(
            diagnostics,
            std::span<const T>(
                face_center_refinement_point.data(),
                face_center_refinement_point.size()));
        return false;
    }

    return check_zero_face_surface_jacobian<T, I>(
        adapt_cell,
        local_zero_entity_id,
        ls_cell,
        ordered_vertices,
        graph_options,
        diagnostics,
        source_cell_id);
}

template <std::floating_point T>
std::vector<int> adapt_zero_entity_vertex_ids(const AdaptCell<T>& ac,
                                              int local_zero_entity_id)
{
    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    const int zid = ac.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim == 0)
        return {zid};

    auto verts = ac.entity_to_vertex[zdim][static_cast<std::int32_t>(zid)];
    std::vector<int> out;
    out.reserve(verts.size());
    for (const auto v : verts)
        out.push_back(static_cast<int>(v));
    return out;
}

inline bool contains_vertex_id(std::span<const int> vertices, int vertex_id)
{
    for (const int v : vertices)
    {
        if (v == vertex_id)
            return true;
    }
    return false;
}

template <std::floating_point T>
std::vector<T> face_normal_for_zero_entity(const AdaptCell<T>& ac,
                                           int local_zero_entity_id,
                                           T tol)
{
    std::vector<T> normal(static_cast<std::size_t>(ac.tdim), T(0));
    if (ac.tdim != 3
        || ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)] != 2)
    {
        return normal;
    }

    const auto face_vertices =
        adapt_zero_entity_vertex_ids<T>(ac, local_zero_entity_id);
    std::vector<T> face_coords;
    face_coords.reserve(face_vertices.size() * 3);
    for (const int v : face_vertices)
    {
        auto p = point_span<T>(ac.vertex_coords, v, 3);
        face_coords.insert(face_coords.end(), p.begin(), p.end());
    }

    graph_criteria::HostFrame<T> host;
    if (!build_face_host_frame<T>(face_coords, host, tol))
        return normal;
    return host.normal;
}

template <std::floating_point T>
std::vector<T> face_normal_for_zero_edge(const AdaptCell<T>& ac,
                                         int local_zero_edge_id,
                                         int level_set_id,
                                         T tol)
{
    std::vector<T> normal(static_cast<std::size_t>(ac.tdim), T(0));
    if (ac.tdim != 3
        || ac.zero_entity_dim[static_cast<std::size_t>(local_zero_edge_id)] != 1)
    {
        return normal;
    }

    const auto edge_vertices =
        adapt_zero_entity_vertex_ids<T>(ac, local_zero_edge_id);
    if (edge_vertices.size() != 2)
        return normal;

    const std::uint64_t bit = std::uint64_t(1) << level_set_id;
    T best_norm = T(0);
    for (int z = 0; z < ac.n_zero_entities(); ++z)
    {
        if (ac.zero_entity_dim[static_cast<std::size_t>(z)] != 2)
            continue;
        if ((ac.zero_entity_zero_mask[static_cast<std::size_t>(z)] & bit) == 0)
            continue;

        const auto face_vertices = adapt_zero_entity_vertex_ids<T>(ac, z);
        if (!contains_vertex_id(std::span<const int>(face_vertices.data(), face_vertices.size()),
                                edge_vertices[0])
            || !contains_vertex_id(
                std::span<const int>(face_vertices.data(), face_vertices.size()),
                edge_vertices[1]))
        {
            continue;
        }

        auto candidate = face_normal_for_zero_entity<T>(ac, z, tol);
        const T n = geom::norm<T>(
            std::span<const T>(candidate.data(), candidate.size()));
        if (n > best_norm)
        {
            best_norm = n;
            normal = std::move(candidate);
        }
    }

    if (best_norm <= tol)
    std::fill(normal.begin(), normal.end(), T(0));
    return normal;
}

template <std::floating_point T>
ZeroEntityGraphDiagnostics<T> make_zero_entity_graph_record(
    int local_zero_entity_id,
    int level_set_id,
    int dimension,
    std::uint64_t zero_mask,
    int host_cell_id,
    cell::type host_cell_type,
    int host_face_id,
    int source_level_set,
    const ReadyCellGraphDiagnostics<T>& entity_diag)
{
    ZeroEntityGraphDiagnostics<T> record;
    record.local_zero_entity_id = local_zero_entity_id;
    record.level_set_id = level_set_id;
    record.dimension = dimension;
    record.zero_mask = zero_mask;
    record.host_cell_id = host_cell_id;
    record.host_cell_type = host_cell_type;
    record.host_face_id = host_face_id;
    record.source_level_set = source_level_set;
    record.accepted = entity_diag.accepted;
    record.checked_edges = entity_diag.checked_edges;
    record.checked_faces = entity_diag.checked_faces;
    record.failed_checks = entity_diag.failed_checks;
    record.failure_reason = entity_diag.accepted
                                ? graph_criteria::FailureReason::none
                                : entity_diag.first_failure_reason;
    record.min_true_transversality = entity_diag.min_true_transversality;
    record.min_host_normal_alignment = entity_diag.min_host_normal_alignment;
    record.max_drift_amplification = entity_diag.max_drift_amplification;
    record.max_relative_correction_distance =
        entity_diag.max_relative_correction_distance;
    record.max_relative_tangential_shift =
        entity_diag.max_relative_tangential_shift;
    record.min_edge_gap_ratio = entity_diag.min_edge_gap_ratio;
    record.min_face_area_ratio = entity_diag.min_face_area_ratio;
    record.min_level_set_gradient_host_alignment =
        entity_diag.min_level_set_gradient_host_alignment;
    record.failed_face_triangle_index =
        entity_diag.first_failed_face_triangle_index;
    record.failed_face_area_ratio =
        entity_diag.first_failed_face_area_ratio;
    record.failed_projection_seed = entity_diag.first_failed_projection_seed;
    record.failed_projection_direction =
        entity_diag.first_failed_projection_direction;
    record.failed_projection_clip_lo =
        entity_diag.first_failed_projection_clip_lo;
    record.failed_projection_clip_hi =
        entity_diag.first_failed_projection_clip_hi;
    record.failed_projection_root_t =
        entity_diag.first_failed_projection_root_t;
    record.requested_refinement_entity_dim =
        entity_diag.first_requested_refinement_entity_dim;
    record.requested_refinement_entity_id =
        entity_diag.first_requested_refinement_entity_id;
    record.requested_refinement_point =
        entity_diag.first_requested_refinement_point;
    record.nodes = entity_diag.nodes;
    return record;
}

template <std::floating_point T>
void accumulate_zero_entity_graph_record(ReadyCellGraphDiagnostics<T>& diagnostics,
                                         const ZeroEntityGraphDiagnostics<T>& record,
                                         int source_cell_id)
{
    diagnostics.checked_edges += record.checked_edges;
    diagnostics.checked_faces += record.checked_faces;
    diagnostics.min_true_transversality =
        std::min(diagnostics.min_true_transversality,
                 record.min_true_transversality);
    diagnostics.min_host_normal_alignment =
        std::min(diagnostics.min_host_normal_alignment,
                 record.min_host_normal_alignment);
    diagnostics.max_drift_amplification =
        std::max(diagnostics.max_drift_amplification,
                 record.max_drift_amplification);
    diagnostics.max_relative_correction_distance =
        std::max(diagnostics.max_relative_correction_distance,
                 record.max_relative_correction_distance);
    diagnostics.max_relative_tangential_shift =
        std::max(diagnostics.max_relative_tangential_shift,
                 record.max_relative_tangential_shift);
    diagnostics.min_edge_gap_ratio =
        std::min(diagnostics.min_edge_gap_ratio, record.min_edge_gap_ratio);
    diagnostics.min_face_area_ratio =
        std::min(diagnostics.min_face_area_ratio, record.min_face_area_ratio);
    diagnostics.min_level_set_gradient_host_alignment =
        std::min(diagnostics.min_level_set_gradient_host_alignment,
                 record.min_level_set_gradient_host_alignment);
    if (diagnostics.first_failed_face_triangle_index < 0
        && record.failed_face_triangle_index >= 0)
    {
        diagnostics.first_failed_face_triangle_index =
            record.failed_face_triangle_index;
        diagnostics.first_failed_face_area_ratio =
            record.failed_face_area_ratio;
    }
    for (const auto& node : record.nodes)
        diagnostics.nodes.push_back(node);
    if (!record.accepted)
    {
        mark_graph_refinement_request<T>(
            diagnostics,
            record.requested_refinement_entity_dim,
            record.requested_refinement_entity_id);
        mark_graph_refinement_point<T>(
            diagnostics,
            std::span<const T>(
                record.requested_refinement_point.data(),
                record.requested_refinement_point.size()));
        const int failed_cell_id =
            source_cell_id >= 0 ? source_cell_id : record.host_cell_id;
        mark_graph_failure<T>(diagnostics, failed_cell_id, record.failure_reason);
    }
}

template <std::floating_point T, std::integral I>
void populate_committed_zero_entity_graph_diagnostics(
    const AdaptCell<T>& adapt_cell,
    const LevelSetCell<T, I>& ls_cell,
    int level_set_id,
    const ReadyCellGraphOptions<T>& graph_options,
    ReadyCellGraphDiagnostics<T>& diagnostics)
{
    const int graph_refinements = diagnostics.graph_refinements;
    diagnostics = ReadyCellGraphDiagnostics<T>{};
    diagnostics.graph_refinements = graph_refinements;
    if (!graph_options.enabled)
        return;

    const std::uint64_t bit = std::uint64_t(1) << level_set_id;
    for (int z = 0; z < adapt_cell.n_zero_entities(); ++z)
    {
        if ((adapt_cell.zero_entity_zero_mask[static_cast<std::size_t>(z)] & bit) == 0)
            continue;

        const int zdim = adapt_cell.zero_entity_dim[static_cast<std::size_t>(z)];
        if (zdim != 1 && zdim != 2)
            continue;

        ReadyCellGraphDiagnostics<T> entity_diag;
        const auto vertices = adapt_zero_entity_vertex_ids<T>(adapt_cell, z);
        if (zdim == 1)
        {
            if (vertices.size() != 2)
            {
                mark_graph_failure<T>(
                    entity_diag, -1, graph_criteria::FailureReason::invalid_input);
            }
            else
            {
                const auto zero_face_normal =
                    face_normal_for_zero_edge<T>(
                        adapt_cell, z, level_set_id,
                        graph_options.criteria.tolerance);
                (void)check_zero_edge_graph<T, I>(
                    adapt_cell,
                    z,
                    ls_cell,
                    point_span<T>(adapt_cell.vertex_coords, vertices[0], adapt_cell.tdim),
                    point_span<T>(adapt_cell.vertex_coords, vertices[1], adapt_cell.tdim),
                    std::span<const T>(zero_face_normal.data(), zero_face_normal.size()),
                    false,
                    find_adapt_edge_by_vertices<T>(
                        adapt_cell, vertices[0], vertices[1]),
                    graph_options,
                    entity_diag,
                    -1);
            }
        }
        else
        {
            if (adapt_cell.tdim != 3 || vertices.size() < 3)
            {
                mark_graph_failure<T>(
                    entity_diag, -1, graph_criteria::FailureReason::invalid_input);
            }
            else
            {
                const auto face_normal =
                    face_normal_for_zero_entity<T>(
                        adapt_cell, z, graph_options.criteria.tolerance);
                const T normal_norm = geom::norm<T>(
                    std::span<const T>(face_normal.data(), face_normal.size()));
                if (normal_norm <= graph_options.criteria.tolerance)
                {
                    mark_graph_failure<T>(
                        entity_diag, -1,
                        graph_criteria::FailureReason::invalid_host_frame);
                }
                else
                {
                    (void)check_zero_face_graph<T, I>(
                        adapt_cell,
                        z,
                        ls_cell,
                        std::span<const int>(vertices.data(), vertices.size()),
                        adapt_cell.vertex_coords,
                        graph_options,
                        entity_diag,
                        -1);
                }
            }
        }

        auto record =
            make_zero_entity_graph_record<T>(
                z,
                level_set_id,
                zdim,
                adapt_cell.zero_entity_zero_mask[static_cast<std::size_t>(z)],
                z < static_cast<int>(adapt_cell.zero_entity_host_cell_id.size())
                    ? adapt_cell.zero_entity_host_cell_id[static_cast<std::size_t>(z)]
                    : -1,
                z < static_cast<int>(adapt_cell.zero_entity_host_cell_type.size())
                    ? adapt_cell.zero_entity_host_cell_type[static_cast<std::size_t>(z)]
                    : cell::type::point,
                z < static_cast<int>(adapt_cell.zero_entity_host_face_id.size())
                    ? adapt_cell.zero_entity_host_face_id[static_cast<std::size_t>(z)]
                    : -1,
                z < static_cast<int>(adapt_cell.zero_entity_source_level_set.size())
                    ? adapt_cell.zero_entity_source_level_set[static_cast<std::size_t>(z)]
                    : -1,
                entity_diag);
        accumulate_zero_entity_graph_record<T>(diagnostics, record, -1);
        diagnostics.zero_entities.push_back(std::move(record));
    }
}

template <std::floating_point T>
int find_top_cell_incident_to_edge(const AdaptCell<T>& adapt_cell, int edge_id)
{
    if (edge_id < 0 || edge_id >= adapt_cell.n_entities(1))
        return -1;

    const auto edge_vertices =
        adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(edge_id)];
    if (edge_vertices.size() != 2)
        return -1;

    const int a = static_cast<int>(edge_vertices[0]);
    const int b = static_cast<int>(edge_vertices[1]);
    const int tdim = adapt_cell.tdim;
    for (int c = 0; c < adapt_cell.n_entities(tdim); ++c)
    {
        const auto verts =
            adapt_cell.entity_to_vertex[tdim][static_cast<std::int32_t>(c)];
        bool has_a = false;
        bool has_b = false;
        for (const auto v : verts)
        {
            has_a = has_a || static_cast<int>(v) == a;
            has_b = has_b || static_cast<int>(v) == b;
        }
        if (has_a && has_b)
            return c;
    }
    return -1;
}

template <std::floating_point T, std::integral I>
ReadyCellGraphDiagnostics<T> check_processed_ready_to_cut_cell_graphs(
    const AdaptCell<T>& adapt_cell,
    const LevelSetCell<T, I>& ls_cell,
    int level_set_id,
    T zero_tol,
    T sign_tol,
    int edge_max_depth,
    bool triangulate_cut_parts,
    const ReadyCellGraphOptions<T>& graph_options)
{
    ReadyCellGraphDiagnostics<T> diagnostics;
    if (!graph_options.enabled)
        return diagnostics;

    bool has_ready = false;
    const int tdim = adapt_cell.tdim;
    for (int c = 0; c < adapt_cell.n_entities(tdim); ++c)
    {
        if (adapt_cell.get_cell_cert_tag(level_set_id, c)
            == CellCertTag::ready_to_cut)
        {
            has_ready = true;
            break;
        }
    }
    if (!has_ready)
        return diagnostics;

    AdaptCell<T> temporary = adapt_cell;
    process_ready_to_cut_cells(
        temporary,
        ls_cell,
        level_set_id,
        zero_tol,
        sign_tol,
        edge_max_depth,
        triangulate_cut_parts);
    fill_all_vertex_signs_from_level_set(
        temporary, ls_cell, level_set_id, zero_tol);
    build_edges(temporary);
    if (temporary.tdim == 3)
        build_faces(temporary);
    recompute_active_level_set_masks(temporary, level_set_id + 1);
    rebuild_zero_entity_inventory(temporary);

    populate_committed_zero_entity_graph_diagnostics(
        temporary, ls_cell, level_set_id, graph_options, diagnostics);

    if (diagnostics.zero_entities.empty())
    {
        mark_graph_failure<T>(
            diagnostics, -1, graph_criteria::FailureReason::invalid_input);
    }

    if (!diagnostics.accepted && diagnostics.first_failed_cell < 0
        && diagnostics.first_requested_refinement_entity_dim == 1)
    {
        diagnostics.first_failed_cell = find_top_cell_incident_to_edge<T>(
            adapt_cell, diagnostics.first_requested_refinement_entity_id);
    }

    return diagnostics;
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

        if (leaf_cell_type != cell::type::triangle
            && leaf_cell_type != cell::type::tetrahedron)
        {
            throw std::runtime_error(
                "process_ready_to_cut_cells: ready_to_cut only implemented for triangle and tetrahedron leaves");
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

        if (leaf_cell_type == cell::type::triangle)
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
                triangulate_cut_parts);
            cell::triangle::cut(
                std::span<const T>(vertex_coords),
                tdim,
                std::span<const T>(ls_values),
                "phi>0",
                positive_part,
                triangulate_cut_parts);

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
            cell::CutCell<T> negative_part;
            cell::CutCell<T> positive_part;
            cell::tetrahedron::cut(
                std::span<const T>(vertex_coords), tdim,
                std::span<const T>(ls_values), "phi<0",
                negative_part, triangulate_cut_parts);
            cell::tetrahedron::cut(
                std::span<const T>(vertex_coords), tdim,
                std::span<const T>(ls_values), "phi>0",
                positive_part, triangulate_cut_parts);

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
ReadyCellGraphDiagnostics<T> check_ready_to_cut_cell_graphs(
    const AdaptCell<T>& adapt_cell,
    const LevelSetCell<T, I>& ls_cell,
    int level_set_id,
    const ReadyCellGraphOptions<T>& graph_options)
{
    ReadyCellGraphDiagnostics<T> diagnostics;

    if (!graph_options.enabled)
        return diagnostics;

    const int tdim = adapt_cell.tdim;
    const int n_cells = adapt_cell.n_entities(tdim);

    for (int c = 0; c < n_cells; ++c)
    {
        if (adapt_cell.get_cell_cert_tag(level_set_id, c)
            != CellCertTag::ready_to_cut)
        {
            continue;
        }

        ++diagnostics.checked_cells;

        const cell::type leaf_cell_type =
            adapt_cell.entity_types[tdim][static_cast<std::size_t>(c)];
        if (leaf_cell_type != cell::type::triangle
            && leaf_cell_type != cell::type::tetrahedron)
        {
            mark_graph_failure<T>(
                diagnostics,
                c,
                graph_criteria::FailureReason::invalid_input);
            return diagnostics;
        }

        AdaptCell<T> temporary = adapt_cell;
        for (int other = 0; other < n_cells; ++other)
        {
            if (other != c
                && temporary.get_cell_cert_tag(level_set_id, other)
                    == CellCertTag::ready_to_cut)
            {
                temporary.set_cell_cert_tag(
                    level_set_id, other, CellCertTag::not_classified);
            }
        }

        process_ready_to_cut_cells(
            temporary,
            ls_cell,
            level_set_id,
            graph_options.criteria.tolerance,
            graph_options.criteria.tolerance,
            /*edge_max_depth=*/20,
            /*triangulate_cut_parts=*/false);

        ReadyCellGraphDiagnostics<T> temporary_zero_diag;
        populate_committed_zero_entity_graph_diagnostics(
            temporary,
            ls_cell,
            level_set_id,
            graph_options,
            temporary_zero_diag);

        bool checked_zero_entity = false;
        for (const auto& record : temporary_zero_diag.zero_entities)
        {
            checked_zero_entity = true;
            accumulate_zero_entity_graph_record<T>(diagnostics, record, c);
            if (!diagnostics.accepted)
                return diagnostics;
        }

        if (!checked_zero_entity)
        {
            mark_graph_failure<T>(
                diagnostics,
                c,
                graph_criteria::FailureReason::invalid_input);
            return diagnostics;
        }
    }

    return diagnostics;
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

template <std::floating_point T, std::integral I>
T graph_requested_edge_split_parameter(
    const AdaptCell<T>& adapt_cell,
    const LevelSetCell<T, I>& ls_cell,
    int level_set_id,
    int edge_id,
    T zero_tol,
    T sign_tol,
    int edge_max_depth)
{
    const T fallback = T(0.5);
    const int n_edges = adapt_cell.n_entities(1);
    if (edge_id < 0 || edge_id >= n_edges)
        return fallback;

    const T endpoint_tol =
        std::max(zero_tol, T(64) * std::numeric_limits<T>::epsilon());
    auto usable = [&](T t)
    {
        return std::isfinite(t) && t > endpoint_tol
            && t < T(1) - endpoint_tol;
    };

    const auto idx =
        static_cast<std::size_t>(level_set_id * n_edges + edge_id);
    if (idx < adapt_cell.edge_one_root_has_value.size()
        && adapt_cell.edge_one_root_has_value[idx]
        && idx < adapt_cell.edge_one_root_param.size()
        && usable(adapt_cell.edge_one_root_param[idx]))
    {
        return adapt_cell.edge_one_root_param[idx];
    }

    if (adapt_cell.get_edge_root_tag(level_set_id, edge_id)
        != EdgeRootTag::one_root)
    {
        return fallback;
    }

    std::vector<T> edge_coeffs;
    gather_adapt_edge_bernstein(adapt_cell, ls_cell, edge_id, edge_coeffs);

    T root_t = fallback;
    if (locate_one_root_parameter(
            std::span<const T>(edge_coeffs), zero_tol, sign_tol,
            edge_max_depth, root_t)
        && usable(root_t))
    {
        return root_t;
    }

    return fallback;
}

template <std::floating_point T, std::integral I>
bool refine_green_on_requested_graph_edge(AdaptCell<T>& adapt_cell,
                                          const LevelSetCell<T, I>& ls_cell,
                                          int level_set_id,
                                          int edge_id,
                                          T zero_tol,
                                          T sign_tol,
                                          int edge_max_depth)
{
    if (edge_id < 0 || edge_id >= adapt_cell.n_entities(1))
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
        static_cast<std::size_t>(level_set_id * n_edges + edge_id);
    adapt_cell.set_edge_root_tag(
        level_set_id, edge_id, EdgeRootTag::multiple_roots);
    adapt_cell.edge_green_split_param[idx] =
        graph_requested_edge_split_parameter<T, I>(
            adapt_cell, ls_cell, level_set_id, edge_id,
            zero_tol, sign_tol, edge_max_depth);
    adapt_cell.edge_green_split_has_value[idx] = 1;
    return refine_green_on_multiple_root_edges(adapt_cell, level_set_id);
}

template <std::floating_point T>
bool refine_green_on_requested_graph_edge_midpoint(AdaptCell<T>& adapt_cell,
                                                   int level_set_id,
                                                   int edge_id)
{
    if (edge_id < 0 || edge_id >= adapt_cell.n_entities(1))
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
        static_cast<std::size_t>(level_set_id * n_edges + edge_id);
    adapt_cell.set_edge_root_tag(
        level_set_id, edge_id, EdgeRootTag::multiple_roots);
    adapt_cell.edge_green_split_param[idx] = T(0.5);
    adapt_cell.edge_green_split_has_value[idx] = 1;
    return refine_green_on_multiple_root_edges(adapt_cell, level_set_id);
}

template <std::floating_point T>
T graph_tetra_signed_volume6(const AdaptCell<T>& adapt_cell,
                             std::span<const int> tet)
{
    if (adapt_cell.tdim != 3 || tet.size() != 4)
        return T(0);
    const auto p0 = point_span<T>(adapt_cell.vertex_coords, tet[0], 3);
    const auto p1 = point_span<T>(adapt_cell.vertex_coords, tet[1], 3);
    const auto p2 = point_span<T>(adapt_cell.vertex_coords, tet[2], 3);
    const auto p3 = point_span<T>(adapt_cell.vertex_coords, tet[3], 3);
    std::array<T, 3> a = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    std::array<T, 3> b = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
    std::array<T, 3> c = {p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]};
    const auto axb = geom::cross<T>(
        std::span<const T>(a.data(), a.size()),
        std::span<const T>(b.data(), b.size()));
    return geom::dot<T>(
        std::span<const T>(axb.data(), axb.size()),
        std::span<const T>(c.data(), c.size()));
}

template <std::floating_point T>
bool graph_point_barycentric_in_tetra(
    const AdaptCell<T>& adapt_cell,
    std::span<const std::int32_t> tet_vertices,
    std::span<const T> point,
    std::array<T, 4>& bary,
    T tol)
{
    if (adapt_cell.tdim != 3 || tet_vertices.size() != 4 || point.size() != 3)
        return false;

    const auto p0 = point_span<T>(adapt_cell.vertex_coords, tet_vertices[0], 3);
    std::vector<T> A(9, T(0));
    std::vector<T> rhs(3, T(0));
    for (int col = 0; col < 3; ++col)
    {
        const auto pc = point_span<T>(
            adapt_cell.vertex_coords, tet_vertices[static_cast<std::size_t>(col + 1)], 3);
        for (int row = 0; row < 3; ++row)
        {
            A[static_cast<std::size_t>(row * 3 + col)] =
                pc[static_cast<std::size_t>(row)] - p0[static_cast<std::size_t>(row)];
        }
    }
    for (int row = 0; row < 3; ++row)
        rhs[static_cast<std::size_t>(row)] =
            point[static_cast<std::size_t>(row)] - p0[static_cast<std::size_t>(row)];

    std::vector<T> x;
    if (!solve_graph_dense_small<T>(std::move(A), std::move(rhs), 3, x, tol))
        return false;

    bary[1] = x[0];
    bary[2] = x[1];
    bary[3] = x[2];
    bary[0] = T(1) - bary[1] - bary[2] - bary[3];

    T sum = T(0);
    for (const T value : bary)
    {
        if (value < -tol || value > T(1) + tol)
            return false;
        sum += value;
    }
    return std::fabs(sum - T(1)) <= T(16) * tol;
}

template <std::floating_point T, std::integral I>
bool refine_green_on_requested_graph_point(AdaptCell<T>& adapt_cell,
                                           const LevelSetCell<T, I>& ls_cell,
                                           int level_set_id,
                                           int failed_cell_id,
                                           std::span<const T> point,
                                           T zero_tol)
{
    if (adapt_cell.tdim != 3 || point.size() != 3)
        return false;

    const T tol = std::max(zero_tol, T(1024) * std::numeric_limits<T>::epsilon());
    struct SplitCell
    {
        int cell_id = -1;
        std::array<T, 4> bary = {};
    };
    std::vector<SplitCell> split_cells;

    auto consider_cell = [&](int c)
    {
        if (c < 0 || c >= adapt_cell.n_entities(3))
            return;
        if (adapt_cell.entity_types[3][static_cast<std::size_t>(c)]
            != cell::type::tetrahedron)
        {
            return;
        }
        const auto verts =
            adapt_cell.entity_to_vertex[3][static_cast<std::int32_t>(c)];
        std::array<T, 4> bary = {};
        if (graph_point_barycentric_in_tetra<T>(
                adapt_cell, verts, point, bary, tol))
        {
            split_cells.push_back({c, bary});
        }
    };

    consider_cell(failed_cell_id);
    if (split_cells.empty())
    {
        for (int c = 0; c < adapt_cell.n_entities(3); ++c)
            consider_cell(c);
    }
    if (split_cells.empty())
        return false;

    for (int v = 0; v < adapt_cell.n_vertices(); ++v)
    {
        const auto existing = point_span<T>(adapt_cell.vertex_coords, v, 3);
        const auto delta = geom::subtract<T>(existing, point);
        if (geom::norm<T>(std::span<const T>(delta.data(), delta.size())) <= tol)
            return false;
    }

    std::vector<T> parent_param(point.begin(), point.end());
    const int new_v = append_vertex_with_parent_info(
        adapt_cell,
        point,
        adapt_cell.tdim,
        adapt_cell.parent_cell_id,
        std::span<const T>(parent_param.data(), parent_param.size()),
        -1);
    set_vertex_sign_for_level_set(
        adapt_cell,
        new_v,
        level_set_id,
        ls_cell.value(point),
        zero_tol);

    std::set<int> split_cell_ids;
    for (const auto& split : split_cells)
        split_cell_ids.insert(split.cell_id);

    EntityAdjacency new_cells;
    new_cells.offsets.push_back(0);
    std::vector<cell::type> new_types;
    std::vector<int> old_cell_ids_for_new_cells;
    std::vector<int> source_cell_ids_for_new_cells;
    std::vector<CellRefinementReason> refinement_reasons_for_new_cells;

    const std::array<std::array<int, 3>, 4> opposite_faces = {{
        {{1, 2, 3}},
        {{0, 2, 3}},
        {{0, 1, 3}},
        {{0, 1, 2}},
    }};

    for (int c = 0; c < adapt_cell.n_entities(3); ++c)
    {
        const auto verts =
            adapt_cell.entity_to_vertex[3][static_cast<std::int32_t>(c)];
        if (split_cell_ids.find(c) == split_cell_ids.end())
        {
            std::vector<int> copy(verts.begin(), verts.end());
            append_top_cell_local(
                new_types, new_cells, adapt_cell.entity_types[3][static_cast<std::size_t>(c)],
                std::span<const int>(copy.data(), copy.size()));
            old_cell_ids_for_new_cells.push_back(c);
            source_cell_ids_for_new_cells.push_back(c);
            refinement_reasons_for_new_cells.push_back(CellRefinementReason::none);
            continue;
        }

        std::array<int, 4> old_tet = {
            static_cast<int>(verts[0]),
            static_cast<int>(verts[1]),
            static_cast<int>(verts[2]),
            static_cast<int>(verts[3])};
        const T old_volume = graph_tetra_signed_volume6<T>(
            adapt_cell,
            std::span<const int>(old_tet.data(), old_tet.size()));

        for (const auto& face : opposite_faces)
        {
            std::array<int, 4> child = {
                new_v,
                old_tet[static_cast<std::size_t>(face[0])],
                old_tet[static_cast<std::size_t>(face[1])],
                old_tet[static_cast<std::size_t>(face[2])]};
            T child_volume = graph_tetra_signed_volume6<T>(
                adapt_cell,
                std::span<const int>(child.data(), child.size()));
            if (std::fabs(child_volume) <= tol * std::max(T(1), std::fabs(old_volume)))
                continue;
            if (old_volume * child_volume < T(0))
            {
                std::swap(child[1], child[2]);
                child_volume = -child_volume;
            }
            (void)child_volume;
            append_top_cell_local(
                new_types, new_cells, cell::type::tetrahedron,
                std::span<const int>(child.data(), child.size()));
            old_cell_ids_for_new_cells.push_back(-1);
            source_cell_ids_for_new_cells.push_back(c);
            refinement_reasons_for_new_cells.push_back(
                CellRefinementReason::graph_green_edge);
        }
    }

    apply_topology_update_preserve_certification<T>(
        adapt_cell,
        std::move(new_types),
        std::move(new_cells),
        std::span<const int>(
            old_cell_ids_for_new_cells.data(), old_cell_ids_for_new_cells.size()),
        std::span<const int>(
            source_cell_ids_for_new_cells.data(), source_cell_ids_for_new_cells.size()),
        std::span<const CellRefinementReason>(
            refinement_reasons_for_new_cells.data(),
            refinement_reasons_for_new_cells.size()));
    return true;
}

template <std::floating_point T>
bool refine_red_on_graph_failed_cell(AdaptCell<T>& adapt_cell,
                                     int level_set_id,
                                     int cell_id)
{
    if (cell_id < 0 || cell_id >= adapt_cell.n_entities(adapt_cell.tdim))
        return false;

    adapt_cell.set_cell_cert_tag(
        level_set_id, cell_id, CellCertTag::ambiguous);
    return refine_red_on_ambiguous_cells(adapt_cell, level_set_id);
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
ReadyCellGraphDiagnostics<T> certify_refine_graph_check_and_process_ready_cells(
    AdaptCell<T>& adapt_cell,
    const LevelSetCell<T, I>& ls_cell,
    int level_set_id,
    int max_iterations,
    T zero_tol, T sign_tol,
    int edge_max_depth,
    bool triangulate_cut_parts,
    const ReadyCellGraphOptions<T>& graph_options)
{
    fill_all_vertex_signs_from_level_set(adapt_cell, ls_cell, level_set_id, zero_tol);

    ReadyCellGraphDiagnostics<T> diagnostics;
    for (int graph_iter = 0; graph_iter <= graph_options.max_refinements; ++graph_iter)
    {
        certify_and_refine(adapt_cell, ls_cell, level_set_id,
                           max_iterations, zero_tol, sign_tol, edge_max_depth);
        fill_all_vertex_signs_from_level_set(adapt_cell, ls_cell, level_set_id, zero_tol);

        diagnostics = check_processed_ready_to_cut_cell_graphs(
            adapt_cell, ls_cell, level_set_id,
            zero_tol, sign_tol, edge_max_depth,
            triangulate_cut_parts, graph_options);
        if (!diagnostics.accepted
            && diagnostics.first_requested_refinement_entity_dim != 1)
        {
            diagnostics = check_ready_to_cut_cell_graphs(
                adapt_cell, ls_cell, level_set_id, graph_options);
        }
        diagnostics.graph_refinements = graph_iter;

        if (diagnostics.accepted)
            break;

        if (graph_iter >= graph_options.max_refinements)
        {
            break;
        }

        bool refined = false;
        if (graph_options.refinement_mode
            == GraphRefinementMode::red_failed_cell)
        {
            refined = refine_red_on_graph_failed_cell<T>(
                adapt_cell, level_set_id, diagnostics.first_failed_cell);
        }
        else if (graph_options.refinement_mode
                 == GraphRefinementMode::green_orthogonal_surface_edge
                 && diagnostics.first_requested_refinement_entity_dim == 1)
        {
            if (!diagnostics.first_requested_refinement_point.empty())
            {
                refined = refine_green_on_requested_graph_point<T, I>(
                    adapt_cell,
                    ls_cell,
                    level_set_id,
                    diagnostics.first_failed_cell,
                    std::span<const T>(
                        diagnostics.first_requested_refinement_point.data(),
                        diagnostics.first_requested_refinement_point.size()),
                    zero_tol);
            }
            if (!refined)
            {
                refined = refine_green_on_requested_graph_edge_midpoint<T>(
                    adapt_cell,
                    level_set_id,
                    diagnostics.first_requested_refinement_entity_id);
            }
        }
        else if ((graph_options.refinement_mode
                  == GraphRefinementMode::green_midpoint_residual
                  || graph_options.refinement_mode
                         == GraphRefinementMode::green_normal_variation)
                 && !diagnostics.first_requested_refinement_point.empty())
        {
            refined = refine_green_on_requested_graph_point<T, I>(
                adapt_cell,
                ls_cell,
                level_set_id,
                diagnostics.first_failed_cell,
                std::span<const T>(
                    diagnostics.first_requested_refinement_point.data(),
                    diagnostics.first_requested_refinement_point.size()),
                zero_tol);
        }
        else if (diagnostics.first_requested_refinement_entity_dim == 1)
        {
            refined = refine_green_on_requested_graph_edge<T>(
                adapt_cell,
                ls_cell,
                level_set_id,
                diagnostics.first_requested_refinement_entity_id,
                zero_tol,
                sign_tol,
                edge_max_depth);
        }

        if (!refined && diagnostics.first_failed_cell >= 0)
        {
            refined = refine_red_on_graph_failed_cell<T>(
                adapt_cell, level_set_id, diagnostics.first_failed_cell);
        }

        if (!refined)
        {
            refined = refine_ready_cell_on_largest_midpoint_value(
                adapt_cell, ls_cell, level_set_id, diagnostics.first_failed_cell);
        }

        if (!refined)
        {
            break;
        }

        fill_all_vertex_signs_from_level_set(
            adapt_cell, ls_cell, level_set_id, zero_tol);
    }

    process_ready_to_cut_cells(adapt_cell, ls_cell, level_set_id,
                               zero_tol, sign_tol, edge_max_depth,
                               triangulate_cut_parts);
    fill_all_vertex_signs_from_level_set(
        adapt_cell, ls_cell, level_set_id, zero_tol);
    build_edges(adapt_cell);
    if (adapt_cell.tdim == 3)
        build_faces(adapt_cell);
    recompute_active_level_set_masks(adapt_cell, level_set_id + 1);
    rebuild_zero_entity_inventory(adapt_cell);
    populate_committed_zero_entity_graph_diagnostics(
        adapt_cell, ls_cell, level_set_id, graph_options, diagnostics);
    return diagnostics;
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
    const ReadyCellGraphOptions<T> graph_options;
    (void)certify_refine_graph_check_and_process_ready_cells(
        adapt_cell,
        ls_cell,
        level_set_id,
        max_iterations,
        zero_tol,
        sign_tol,
        edge_max_depth,
        triangulate_cut_parts,
        graph_options);
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

template ReadyCellGraphDiagnostics<double> check_ready_to_cut_cell_graphs(
    const AdaptCell<double>&,
    const LevelSetCell<double, int>&,
    int,
    const ReadyCellGraphOptions<double>&);
template ReadyCellGraphDiagnostics<float> check_ready_to_cut_cell_graphs(
    const AdaptCell<float>&,
    const LevelSetCell<float, int>&,
    int,
    const ReadyCellGraphOptions<float>&);
template ReadyCellGraphDiagnostics<double> check_ready_to_cut_cell_graphs(
    const AdaptCell<double>&,
    const LevelSetCell<double, long>&,
    int,
    const ReadyCellGraphOptions<double>&);
template ReadyCellGraphDiagnostics<float> check_ready_to_cut_cell_graphs(
    const AdaptCell<float>&,
    const LevelSetCell<float, long>&,
    int,
    const ReadyCellGraphOptions<float>&);

template void populate_committed_zero_entity_graph_diagnostics(
    const AdaptCell<double>&,
    const LevelSetCell<double, int>&,
    int,
    const ReadyCellGraphOptions<double>&,
    ReadyCellGraphDiagnostics<double>&);
template void populate_committed_zero_entity_graph_diagnostics(
    const AdaptCell<float>&,
    const LevelSetCell<float, int>&,
    int,
    const ReadyCellGraphOptions<float>&,
    ReadyCellGraphDiagnostics<float>&);
template void populate_committed_zero_entity_graph_diagnostics(
    const AdaptCell<double>&,
    const LevelSetCell<double, long>&,
    int,
    const ReadyCellGraphOptions<double>&,
    ReadyCellGraphDiagnostics<double>&);
template void populate_committed_zero_entity_graph_diagnostics(
    const AdaptCell<float>&,
    const LevelSetCell<float, long>&,
    int,
    const ReadyCellGraphOptions<float>&,
    ReadyCellGraphDiagnostics<float>&);

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

template ReadyCellGraphDiagnostics<double>
certify_refine_graph_check_and_process_ready_cells(
    AdaptCell<double>&,
    const LevelSetCell<double, int>&,
    int,
    int,
    double,
    double,
    int,
    bool,
    const ReadyCellGraphOptions<double>&);
template ReadyCellGraphDiagnostics<float>
certify_refine_graph_check_and_process_ready_cells(
    AdaptCell<float>&,
    const LevelSetCell<float, int>&,
    int,
    int,
    float,
    float,
    int,
    bool,
    const ReadyCellGraphOptions<float>&);
template ReadyCellGraphDiagnostics<double>
certify_refine_graph_check_and_process_ready_cells(
    AdaptCell<double>&,
    const LevelSetCell<double, long>&,
    int,
    int,
    double,
    double,
    int,
    bool,
    const ReadyCellGraphOptions<double>&);
template ReadyCellGraphDiagnostics<float>
certify_refine_graph_check_and_process_ready_cells(
    AdaptCell<float>&,
    const LevelSetCell<float, long>&,
    int,
    int,
    float,
    float,
    int,
    bool,
    const ReadyCellGraphOptions<float>&);

} // namespace cutcells
