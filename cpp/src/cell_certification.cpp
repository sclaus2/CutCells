// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "cell_certification.h"
#include "bernstein.h"
#include "cell_topology.h"
#include "cut_cell.h"
#include "cut_tetrahedron.h"
#include "cut_triangle.h"
#include "edge_certification.h"
#include "reference_cell.h"
#include "refine_cell.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <map>
#include <stdexcept>

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

template <std::floating_point T, std::integral I>
void gather_adapt_edge_bernstein(const AdaptCell<T>& adapt_cell,
                                 const LevelSetCell<T, I>& ls_cell,
                                 int edge_id,
                                 std::vector<T>& edge_coeffs)
{
    int parent_edge_id = -1;
    if (edge_is_on_single_parent_edge(adapt_cell, edge_id, parent_edge_id))
    {
        extract_parent_edge_bernstein(
            ls_cell.cell_type, ls_cell.bernstein_order,
            std::span<const T>(ls_cell.bernstein_coeffs),
            parent_edge_id, edge_coeffs);
        return;
    }

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

            if (token >= 0 && token < 100)
            {
                const int local_edge_id = token;
                if (local_edge_id >= 0
                    && local_edge_id < static_cast<int>(old_edge_ids_by_local_edge.size()))
                {
                    const int old_edge_id = old_edge_ids_by_local_edge[static_cast<std::size_t>(local_edge_id)];
                    int parent_edge_id = -1;
                    if (old_edge_id >= 0
                        && edge_is_on_single_parent_edge(adapt_cell, old_edge_id, parent_edge_id))
                    {
                        parent_dim = 1;
                        parent_id = parent_edge_id;
                    }
                }
            }

            vertex_id = append_vertex_with_parent_info(
                adapt_cell, x, parent_dim, parent_id,
                std::span<const T>(), /*source_edge_id=*/-1);

            const T value = ls_cell.value(x);
            set_vertex_sign_for_level_set(
                adapt_cell, vertex_id, level_set_id, value, zero_tol);
        }

        token_to_vertex[token] = vertex_id;
    }

    const int n_part_cells = cell::num_cells(part);
    for (int pc = 0; pc < n_part_cells; ++pc)
    {
        auto part_vertices = cell::cell_vertices(part, pc);
        std::vector<int> mapped(part_vertices.size(), -1);
        for (std::size_t j = 0; j < part_vertices.size(); ++j)
        {
            const int part_lv = part_vertices[j];
            const int raw_token = part._vertex_parent_entity[static_cast<std::size_t>(part_lv)];
            const int token = cell::vtk_parent_entity_token_to_basix(leaf_cell_type, raw_token);
            mapped[j] = token_to_vertex.at(token);
        }

        append_top_cell_local(new_types, new_cells,
                              part._types[static_cast<std::size_t>(pc)],
                              std::span<const int>(mapped));
        old_cell_ids_for_new_cells.push_back(-1);
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
    const int n_edges = adapt_cell.n_entities(1);

    int num_one_root_edges = 0;
    bool has_multiple_roots = false;
    bool has_zero_edge = false;

    for (int e = 0; e < n_edges; ++e)
    {
        auto edge_verts = adapt_cell.entity_to_vertex[1][static_cast<std::int32_t>(e)];
        bool v0_in = false;
        bool v1_in = false;
        for (auto cv : cell_verts)
        {
            if (cv == edge_verts[0]) v0_in = true;
            if (cv == edge_verts[1]) v1_in = true;
        }
        if (!v0_in || !v1_in)
            continue;

        const EdgeRootTag etag = adapt_cell.get_edge_root_tag(level_set_id, e);
        if (etag == EdgeRootTag::one_root)
            ++num_one_root_edges;
        else if (etag == EdgeRootTag::multiple_roots)
            has_multiple_roots = true;
        else if (etag == EdgeRootTag::zero)
            has_zero_edge = true;
    }

    if (has_multiple_roots || has_zero_edge)
        return CellCertTag::cut;

    if (subcell_type == cell::type::triangle && num_one_root_edges == 2)
        return CellCertTag::ready_to_cut;

    if (subcell_type == cell::type::tetrahedron
        && (num_one_root_edges == 3 || num_one_root_edges == 4))
    {
        return CellCertTag::ready_to_cut;
    }

    if (num_one_root_edges > 0)
        return CellCertTag::cut;

    return CellCertTag::not_classified;
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
        return CellCertTag::zero;
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
                                int edge_max_depth)
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

        std::array<int, 6> old_edge_ids_by_local_edge;
        old_edge_ids_by_local_edge.fill(-1);
        const auto ledges = cell::edges(leaf_cell_type);
        for (std::size_t le = 0; le < ledges.size(); ++le)
        {
            const int a = old_cell_vertices[static_cast<std::size_t>(ledges[le][0])];
            const int b = old_cell_vertices[static_cast<std::size_t>(ledges[le][1])];
            old_edge_ids_by_local_edge[le] =
                edge_lookup.at({std::min(a, b), std::max(a, b)});
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
            cell::CutCell<T> negative_part;
            cell::CutCell<T> positive_part;
            cell::triangle::cut(
                std::span<const T>(vertex_coords),
                tdim,
                std::span<const T>(ls_values),
                "phi<0",
                negative_part,
                /*triangulate=*/false);
            cell::triangle::cut(
                std::span<const T>(vertex_coords),
                tdim,
                std::span<const T>(ls_values),
                "phi>0",
                positive_part,
                /*triangulate=*/false);

            append_ready_cut_part_cells(
                adapt_cell, ls_cell, level_set_id, c, negative_part, leaf_cell_type,
                old_cell_vertices, std::span<const int>(old_edge_ids_by_local_edge.data(), 3),
                new_types, new_cells, old_cell_ids_for_new_cells, explicit_current_ls_tags,
                token_to_vertex, zero_tol, CellCertTag::negative);
            append_ready_cut_part_cells(
                adapt_cell, ls_cell, level_set_id, c, positive_part, leaf_cell_type,
                old_cell_vertices, std::span<const int>(old_edge_ids_by_local_edge.data(), 3),
                new_types, new_cells, old_cell_ids_for_new_cells, explicit_current_ls_tags,
                token_to_vertex, zero_tol, CellCertTag::positive);
        }
        else
        {
            cell::CutCell<T> negative_part;
            cell::CutCell<T> positive_part;
            cell::tetrahedron::cut(
                std::span<const T>(vertex_coords), tdim,
                std::span<const T>(ls_values), "phi<0",
                negative_part, /*triangulate=*/false);
            cell::tetrahedron::cut(
                std::span<const T>(vertex_coords), tdim,
                std::span<const T>(ls_values), "phi>0",
                positive_part, /*triangulate=*/false);

            append_ready_cut_part_cells(
                adapt_cell, ls_cell, level_set_id, c, negative_part, leaf_cell_type,
                old_cell_vertices, std::span<const int>(old_edge_ids_by_local_edge.data(), 6),
                new_types, new_cells, old_cell_ids_for_new_cells, explicit_current_ls_tags,
                token_to_vertex, zero_tol, CellCertTag::negative);
            append_ready_cut_part_cells(
                adapt_cell, ls_cell, level_set_id, c, positive_part, leaf_cell_type,
                old_cell_vertices, std::span<const int>(old_edge_ids_by_local_edge.data(), 6),
                new_types, new_cells, old_cell_ids_for_new_cells, explicit_current_ls_tags,
                token_to_vertex, zero_tol, CellCertTag::positive);
        }
    }

    apply_topology_update_preserve_certification(
        adapt_cell, std::move(new_types), std::move(new_cells),
        std::span<const int>(old_cell_ids_for_new_cells));

    const int new_num_cells = adapt_cell.n_entities(tdim);
    for (int c = 0; c < new_num_cells; ++c)
    {
        if (explicit_current_ls_tags[static_cast<std::size_t>(c)] == CellCertTag::not_classified)
            continue;
        adapt_cell.set_cell_cert_tag(
            level_set_id, c, explicit_current_ls_tags[static_cast<std::size_t>(c)]);
    }

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
                continue;
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
                continue;
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
                                            int edge_max_depth)
{
    fill_all_vertex_signs_from_level_set(adapt_cell, ls_cell, level_set_id, zero_tol);
    certify_and_refine(adapt_cell, ls_cell, level_set_id,
                       max_iterations, zero_tol, sign_tol, edge_max_depth);
    fill_all_vertex_signs_from_level_set(adapt_cell, ls_cell, level_set_id, zero_tol);
    process_ready_to_cut_cells(adapt_cell, ls_cell, level_set_id,
                               zero_tol, sign_tol, edge_max_depth);
    fill_all_vertex_signs_from_level_set(adapt_cell, ls_cell, level_set_id, zero_tol);
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
                                         int, double, double, int);
template void process_ready_to_cut_cells(AdaptCell<float>&,
                                         const LevelSetCell<float, int>&,
                                         int, float, float, int);
template void process_ready_to_cut_cells(AdaptCell<double>&,
                                         const LevelSetCell<double, long>&,
                                         int, double, double, int);
template void process_ready_to_cut_cells(AdaptCell<float>&,
                                         const LevelSetCell<float, long>&,
                                         int, float, float, int);

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
                                                     int, int, double, double, int);
template void certify_refine_and_process_ready_cells(AdaptCell<float>&,
                                                     const LevelSetCell<float, int>&,
                                                     int, int, float, float, int);
template void certify_refine_and_process_ready_cells(AdaptCell<double>&,
                                                     const LevelSetCell<double, long>&,
                                                     int, int, double, double, int);
template void certify_refine_and_process_ready_cells(AdaptCell<float>&,
                                                     const LevelSetCell<float, long>&,
                                                     int, int, float, float, int);

} // namespace cutcells
