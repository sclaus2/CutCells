// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "ho_mesh_part_output.h"

#include "mapping.h"
#include "quadrature_tables.h"
#include "reference_cell.h"
#include "cell_topology.h"
#include "quad_midpoint_split.h"
#include "prism_midpoint_split.h"
#include "triangulation.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace cutcells::output
{
namespace
{

struct SelectedEntity
{
    cell::type type = cell::type::point;
    std::vector<int> vertices;
};

inline bool is_simplex(cell::type cell_type)
{
    return cell_type == cell::type::interval
        || cell_type == cell::type::triangle
        || cell_type == cell::type::tetrahedron;
}

inline cell::type simplex_type_for_dim(int dim)
{
    switch (dim)
    {
        case 1:
            return cell::type::interval;
        case 2:
            return cell::type::triangle;
        case 3:
            return cell::type::tetrahedron;
        default:
            throw std::runtime_error("Unsupported simplex dimension");
    }
}

template <std::floating_point T, std::integral I>
std::vector<T> parent_cell_vertex_coords_vtk(const MeshView<T, I>& mesh, I cell_id)
{
    const auto ctype = mesh.cell_type(cell_id);
    const int nv = cell::get_num_vertices(ctype);
    std::vector<T> coords(static_cast<std::size_t>(nv * mesh.gdim), T(0));

    for (int vtk_v = 0; vtk_v < nv; ++vtk_v)
    {
        const int local_v = mesh.vtk_vertex_order
                                ? vtk_v
                                : cell::vtk_to_basix_vertex(ctype, vtk_v);
        const I node_id = mesh.cell_node(cell_id, static_cast<I>(local_v));
        const T* x = mesh.node(node_id);
        for (int d = 0; d < mesh.gdim; ++d)
            coords[static_cast<std::size_t>(vtk_v * mesh.gdim + d)] = x[d];
    }

    return coords;
}

template <std::floating_point T>
bool vertex_is_zero_for_level_set(const AdaptCell<T>& adapt_cell,
                                  int vertex_id,
                                  int level_set_id)
{
    const std::uint64_t bit = std::uint64_t(1) << level_set_id;
    return (adapt_cell.zero_mask_per_vertex[static_cast<std::size_t>(vertex_id)] & bit) != 0;
}

template <std::floating_point T>
bool cell_contains_all_vertices(std::span<const std::int32_t> cell_verts,
                                std::span<const int> entity_verts)
{
    for (const int v : entity_verts)
    {
        bool found = false;
        for (const auto cv : cell_verts)
        {
            if (cv == v)
            {
                found = true;
                break;
            }
        }
        if (!found)
            return false;
    }
    return true;
}

template <std::floating_point T>
void vertex_state_for_level_set(const AdaptCell<T>& ac,
                                int vertex_id,
                                int level_set_id,
                                bool& is_negative,
                                bool& is_positive,
                                bool& is_zero)
{
    const std::uint64_t bit = std::uint64_t(1) << level_set_id;
    const auto zm = ac.zero_mask_per_vertex[static_cast<std::size_t>(vertex_id)];
    const auto nm = ac.negative_mask_per_vertex[static_cast<std::size_t>(vertex_id)];
    is_zero = (zm & bit) != 0;
    is_negative = !is_zero && ((nm & bit) != 0);
    is_positive = !is_zero && !is_negative;
}

template <std::floating_point T, std::integral I>
bool leaf_cell_matches_sign_requirements(
    const AdaptCell<T>& ac,
    std::span<const std::int32_t> cell_verts,
    const SelectionExpr& expr,
    std::uint64_t cut_cell_active_mask,
    const BackgroundMeshData<T, I>& bg,
    I parent_cell_id)
{
    const int nls = std::min(bg.num_level_sets, 64);
    for (int li = 0; li < nls; ++li)
    {
        const std::uint64_t bit = std::uint64_t(1) << li;
        const bool require_neg = (expr.negative_required & bit) != 0;
        const bool require_pos = (expr.positive_required & bit) != 0;
        if (!require_neg && !require_pos)
            continue;

        const bool ls_is_active = (cut_cell_active_mask & bit) != 0;
        if (!ls_is_active)
        {
            const auto dom = bg.domain(li, parent_cell_id);
            if (require_neg && dom != cell::domain::inside)
                return false;
            if (require_pos && dom != cell::domain::outside)
                return false;
            continue;
        }

        bool has_neg = false;
        bool has_pos = false;
        for (const auto cv : cell_verts)
        {
            bool is_neg = false;
            bool is_pos = false;
            bool is_zero = false;
            vertex_state_for_level_set(ac, static_cast<int>(cv), li, is_neg, is_pos, is_zero);
            has_neg = has_neg || is_neg;
            has_pos = has_pos || is_pos;
        }

        if (require_neg)
        {
            if (has_pos || !has_neg)
                return false;
        }
        if (require_pos)
        {
            if (has_neg || !has_pos)
                return false;
        }
    }

    return true;
}

template <std::floating_point T, std::integral I>
bool zero_entity_matches(const HOMeshPart<T, I>& part,
                         const AdaptCell<T>& ac,
                         int zero_entity_index,
                         std::uint64_t cut_cell_active_mask,
                         I parent_cell_id)
{
    if (ac.zero_entity_dim[static_cast<std::size_t>(zero_entity_index)] != part.dim)
        return false;

    const auto zero_mask = ac.zero_entity_zero_mask[static_cast<std::size_t>(zero_entity_index)];
    if ((zero_mask & part.expr.zero_required) != part.expr.zero_required)
        return false;

    if (part.expr.negative_required == 0 && part.expr.positive_required == 0)
        return true;

    std::vector<int> zero_verts;
    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(zero_entity_index)];
    const int zid = ac.zero_entity_id[static_cast<std::size_t>(zero_entity_index)];
    if (zdim == 0)
    {
        zero_verts.push_back(zid);
    }
    else
    {
        auto verts = ac.entity_to_vertex[zdim][static_cast<std::int32_t>(zid)];
        zero_verts.reserve(verts.size());
        for (const auto v : verts)
            zero_verts.push_back(static_cast<int>(v));
    }

    const int tdim = ac.tdim;
    const int n_cells = ac.n_entities(tdim);
    for (int c = 0; c < n_cells; ++c)
    {
        auto cell_verts = ac.entity_to_vertex[tdim][static_cast<std::int32_t>(c)];
        if (!cell_contains_all_vertices<T>(cell_verts, std::span<const int>(zero_verts)))
            continue;

        if (leaf_cell_matches_sign_requirements(
                ac, cell_verts, part.expr, cut_cell_active_mask, *part.bg, parent_cell_id))
        {
            return true;
        }
    }

    return false;
}

template <std::floating_point T, std::integral I>
std::vector<SelectedEntity> selected_entities(const HOMeshPart<T, I>& part,
                                              const AdaptCell<T>& adapt_cell,
                                              int cut_cell_id)
{
    std::vector<SelectedEntity> entities;
    const I parent_cell_id =
        part.cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_cell_id)];
    const std::uint64_t cut_active_mask =
        part.cut_cells->active_level_set_mask[static_cast<std::size_t>(cut_cell_id)];

    if (part.dim == adapt_cell.tdim)
    {
        const int n_cells = adapt_cell.n_entities(adapt_cell.tdim);
        entities.reserve(static_cast<std::size_t>(n_cells));
        for (int c = 0; c < n_cells; ++c)
        {
            auto verts = adapt_cell.entity_to_vertex[adapt_cell.tdim][static_cast<std::int32_t>(c)];
            if (!leaf_cell_matches_sign_requirements(
                    adapt_cell, verts, part.expr, cut_active_mask, *part.bg, parent_cell_id))
            {
                continue;
            }

            SelectedEntity entity;
            entity.type = adapt_cell.entity_types[adapt_cell.tdim][static_cast<std::size_t>(c)];
            entity.vertices.assign(verts.begin(), verts.end());
            entities.push_back(std::move(entity));
        }
        return entities;
    }

    if ((cut_active_mask & part.expr.zero_required) != part.expr.zero_required)
        return entities;

    const int n_zero = adapt_cell.n_zero_entities();
    entities.reserve(static_cast<std::size_t>(n_zero));
    for (int z = 0; z < n_zero; ++z)
    {
        if (!zero_entity_matches(part, adapt_cell, z, cut_active_mask, parent_cell_id))
            continue;

        const int zdim = adapt_cell.zero_entity_dim[static_cast<std::size_t>(z)];
        const int zid = adapt_cell.zero_entity_id[static_cast<std::size_t>(z)];
        SelectedEntity entity;
        entity.type = (zdim == 0)
                          ? cell::type::point
                          : adapt_cell.entity_types[zdim][static_cast<std::size_t>(zid)];
        if (zdim == 0)
        {
            entity.vertices.push_back(zid);
        }
        else
        {
            auto verts = adapt_cell.entity_to_vertex[zdim][static_cast<std::int32_t>(zid)];
            entity.vertices.assign(verts.begin(), verts.end());
        }

        entities.push_back(std::move(entity));
    }

    return entities;
}

template <std::floating_point T>
std::vector<T> entity_reference_coords(const AdaptCell<T>& adapt_cell,
                                       std::span<const int> entity_vertices)
{
    std::vector<T> coords(
        static_cast<std::size_t>(entity_vertices.size() * adapt_cell.tdim), T(0));

    for (std::size_t j = 0; j < entity_vertices.size(); ++j)
    {
        const int gv = entity_vertices[j];
        for (int d = 0; d < adapt_cell.tdim; ++d)
        {
            coords[static_cast<std::size_t>(j * adapt_cell.tdim + d)] =
                adapt_cell.vertex_coords[static_cast<std::size_t>(gv * adapt_cell.tdim + d)];
        }
    }

    return coords;
}

template <std::floating_point T>
void gather_subcell_vertices(std::span<const T> coords,
                             int coord_dim,
                             std::span<const int> vertex_ids,
                             std::vector<T>& out)
{
    out.resize(static_cast<std::size_t>(vertex_ids.size() * coord_dim));
    for (std::size_t j = 0; j < vertex_ids.size(); ++j)
    {
        const int local_v = vertex_ids[j];
        for (int d = 0; d < coord_dim; ++d)
        {
            out[static_cast<std::size_t>(j * coord_dim + d)] =
                coords[static_cast<std::size_t>(local_v * coord_dim + d)];
        }
    }
}

template <std::floating_point T>
void map_canonical_to_subcell_points(const T* canonical_points,
                                     int num_points,
                                     int simplex_dim,
                                     const T* subcell_vertices,
                                     int parent_tdim,
                                     T* out_points)
{
    const T* v0 = subcell_vertices;

    for (int q = 0; q < num_points; ++q)
    {
        const T* X = canonical_points + q * simplex_dim;
        T* x = out_points + q * parent_tdim;

        for (int d = 0; d < parent_tdim; ++d)
            x[d] = v0[d];

        for (int i = 1; i <= simplex_dim; ++i)
        {
            const T* vi = subcell_vertices + i * parent_tdim;
            for (int d = 0; d < parent_tdim; ++d)
                x[d] += X[i - 1] * (vi[d] - v0[d]);
        }
    }
}

template <std::floating_point T>
T simplex_physical_measure(const T* vertices,
                           int simplex_dim,
                           int gdim)
{
    T J[9] = {};
    const T* v0 = vertices;

    for (int col = 0; col < simplex_dim; ++col)
    {
        const T* vi = vertices + (col + 1) * gdim;
        for (int row = 0; row < gdim; ++row)
            J[col * gdim + row] = vi[row] - v0[row];
    }

    if (simplex_dim == gdim)
    {
        if (simplex_dim == 1)
            return std::abs(J[0]);
        if (simplex_dim == 2)
            return std::abs(J[0] * J[3] - J[2] * J[1]);

        const T det =
            J[0] * (J[4] * J[8] - J[7] * J[5])
          - J[3] * (J[1] * J[8] - J[7] * J[2])
          + J[6] * (J[1] * J[5] - J[4] * J[2]);
        return std::abs(det);
    }

    T G[9] = {};
    for (int i = 0; i < simplex_dim; ++i)
    {
        for (int j = 0; j < simplex_dim; ++j)
        {
            T sum = 0;
            for (int k = 0; k < gdim; ++k)
                sum += J[i * gdim + k] * J[j * gdim + k];
            G[i * simplex_dim + j] = sum;
        }
    }

    if (simplex_dim == 1)
        return std::sqrt(G[0]);
    if (simplex_dim == 2)
        return std::sqrt(G[0] * G[3] - G[1] * G[2]);

    const T det =
        G[0] * (G[4] * G[8] - G[7] * G[5])
      - G[3] * (G[1] * G[8] - G[7] * G[2])
      + G[6] * (G[1] * G[5] - G[4] * G[2]);
    return std::sqrt(det);
}

template <std::floating_point T>
std::vector<T> reorder_vertex_coords_to_vtk(cell::type cell_type,
                                            std::span<const T> coords,
                                            int coord_dim)
{
    const auto perm = cell::basix_to_vtk_vertex_permutation(cell_type);
    return cell::permute_vertex_data(coords, coord_dim, perm);
}

inline std::vector<int> vtk_local_ids_from_basix(cell::type cell_type)
{
    const auto perm = cell::vtk_to_basix_vertex_permutation(cell_type);
    return std::vector<int>(perm.begin(), perm.end());
}

template <std::floating_point T>
void append_mesh_entity(mesh::CutMesh<T>& out,
                        std::span<const T> physical_coords,
                        int gdim,
                        cell::type cell_type,
                        int parent_cell_id,
                        bool triangulate,
                        bool input_is_basix,
                        cell::type source_parent_type,
                        std::span<const int> root_vertex_flags)
{
    if (out._gdim == 0)
        out._gdim = gdim;
    if (out._tdim == 0)
        out._tdim = cell::get_tdim(cell_type);

    // CutMesh stores vertex coordinates in basix ordering internally.
    // Basix-ordered inputs are stored as-is; VTK-ordered inputs are
    // permuted to basix ordering before storage.
    std::vector<T> reordered_coords;
    std::span<const T> output_coords = physical_coords;
    if (!input_is_basix && !is_simplex(cell_type))
    {
        const auto perm = cell::vtk_to_basix_vertex_permutation(cell_type);
        reordered_coords = cell::permute_vertex_data(physical_coords, gdim, perm);
        output_coords = std::span<const T>(reordered_coords.data(), reordered_coords.size());
    }

    const int nv = static_cast<int>(output_coords.size()) / gdim;
    const int vertex_base = out._num_vertices;

    if (triangulate
        && cell_type == cell::type::quadrilateral
        && source_parent_type == cell::type::triangle
        && root_vertex_flags.size() == static_cast<std::size_t>(nv))
    {
        std::array<int, 4> quad_tokens = {};
        for (int i = 0; i < nv; ++i)
        {
            quad_tokens[static_cast<std::size_t>(i)]
                = root_vertex_flags[static_cast<std::size_t>(i)] ? i : 100 + i;
        }

        int next_token_base = 200;
        auto split = cell::quad_midpoint::split_triangle_derived_quadrilateral<T>(
            output_coords,
            gdim,
            std::span<const int>(quad_tokens.data(), quad_tokens.size()),
            next_token_base);

        out._vertex_coords.insert(
            out._vertex_coords.end(), output_coords.begin(), output_coords.end());
        out._num_vertices += nv;

        const int midpoint_vertex_base = out._num_vertices;
        out._vertex_coords.insert(
            out._vertex_coords.end(),
            split.added_vertex_coords.begin(),
            split.added_vertex_coords.end());
        out._num_vertices += static_cast<int>(split.added_vertex_tokens.size());

        std::array<int, 256> token_to_vertex;
        token_to_vertex.fill(-1);
        for (int i = 0; i < nv; ++i)
            token_to_vertex[quad_tokens[static_cast<std::size_t>(i)]] = vertex_base + i;
        for (std::size_t i = 0; i < split.added_vertex_tokens.size(); ++i)
        {
            token_to_vertex[split.added_vertex_tokens[i]]
                = midpoint_vertex_base + static_cast<int>(i);
        }

        for (const auto& tri_tokens : split.triangles)
        {
            for (int k = 0; k < 3; ++k)
                out._connectivity.push_back(token_to_vertex[tri_tokens[static_cast<std::size_t>(k)]]);
            out._offset.push_back(static_cast<int>(out._connectivity.size()));
            out._types.push_back(cell::type::triangle);
            out._parent_map.push_back(parent_cell_id);
            out._num_cells += 1;
        }
        return;
    }

    if (triangulate
        && cell_type == cell::type::prism
        && source_parent_type == cell::type::tetrahedron
        && root_vertex_flags.size() == static_cast<std::size_t>(nv))
    {
        std::array<int, 6> prism_tokens = {};
        for (int i = 0; i < nv; ++i)
            prism_tokens[static_cast<std::size_t>(i)] = root_vertex_flags[static_cast<std::size_t>(i)] ? i : 100 + i;

        int next_token_base = 200;
        auto split = cell::prism_midpoint::split_tetra_derived_prism<T>(
            output_coords,
            gdim,
            std::span<const int>(prism_tokens.data(), prism_tokens.size()),
            next_token_base);

        out._vertex_coords.insert(
            out._vertex_coords.end(), output_coords.begin(), output_coords.end());
        out._num_vertices += nv;

        const int midpoint_vertex_base = out._num_vertices;
        out._vertex_coords.insert(
            out._vertex_coords.end(),
            split.added_vertex_coords.begin(),
            split.added_vertex_coords.end());
        out._num_vertices += static_cast<int>(split.added_vertex_tokens.size());

        std::array<int, 256> token_to_vertex;
        token_to_vertex.fill(-1);
        for (int i = 0; i < nv; ++i)
            token_to_vertex[prism_tokens[static_cast<std::size_t>(i)]] = vertex_base + i;
        for (std::size_t i = 0; i < split.added_vertex_tokens.size(); ++i)
            token_to_vertex[split.added_vertex_tokens[i]] = midpoint_vertex_base + static_cast<int>(i);

        for (const auto& tet_tokens : split.tets)
        {
            for (int k = 0; k < 4; ++k)
                out._connectivity.push_back(token_to_vertex[tet_tokens[static_cast<std::size_t>(k)]]);
            out._offset.push_back(static_cast<int>(out._connectivity.size()));
            out._types.push_back(cell::type::tetrahedron);
            out._parent_map.push_back(parent_cell_id);
            out._num_cells += 1;
        }
        return;
    }

    out._vertex_coords.insert(
        out._vertex_coords.end(), output_coords.begin(), output_coords.end());
    out._num_vertices += nv;

    if (triangulate && !is_simplex(cell_type) && cell::get_tdim(cell_type) >= 2)
    {
        std::vector<int> local_ids(static_cast<std::size_t>(nv));
        std::iota(local_ids.begin(), local_ids.end(), 0);

        std::vector<std::vector<int>> simplices;
        cell::triangulation(cell_type, local_ids.data(), simplices);
        const auto simplex_type = simplex_type_for_dim(cell::get_tdim(cell_type));

        for (const auto& simplex : simplices)
        {
            for (int lv : simplex)
                out._connectivity.push_back(vertex_base + lv);
            out._offset.push_back(static_cast<int>(out._connectivity.size()));
            out._types.push_back(simplex_type);
            out._parent_map.push_back(parent_cell_id);
            out._num_cells += 1;
        }
        return;
    }

    for (int lv = 0; lv < nv; ++lv)
        out._connectivity.push_back(vertex_base + lv);
    out._offset.push_back(static_cast<int>(out._connectivity.size()));
    out._types.push_back(cell_type);
    out._parent_map.push_back(parent_cell_id);
    out._num_cells += 1;
}

template <std::floating_point T>
void append_simplex_quadrature(quadrature::QuadratureRules<T>& rules,
                               cell::type simplex_type,
                               std::span<const T> ref_vertices,
                               std::span<const T> physical_vertices,
                               int parent_tdim,
                               int gdim,
                               int order)
{
    const auto ref_rule = quadrature::get_reference_rule<T>(simplex_type, order);
    const int num_points = ref_rule._num_points;
    const int simplex_dim = ref_rule._tdim;

    std::vector<T> mapped_ref_points(
        static_cast<std::size_t>(num_points * parent_tdim), T(0));
    map_canonical_to_subcell_points(
        ref_rule._points.data(),
        num_points,
        simplex_dim,
        ref_vertices.data(),
        parent_tdim,
        mapped_ref_points.data());

    rules._points.insert(
        rules._points.end(), mapped_ref_points.begin(), mapped_ref_points.end());

    const T measure = simplex_physical_measure(
        physical_vertices.data(), simplex_dim, gdim);
    for (int q = 0; q < num_points; ++q)
        rules._weights.push_back(ref_rule._weights[q] * measure);
}

template <std::floating_point T>
void append_entity_quadrature(quadrature::QuadratureRules<T>& rules,
                              cell::type cell_type,
                              std::span<const T> ref_vertices,
                              std::span<const T> physical_vertices,
                              int parent_tdim,
                              int gdim,
                              int parent_cell_id,
                              int order,
                              bool triangulate,
                              bool input_is_basix,
                              cell::type source_parent_type,
                              std::span<const int> root_vertex_flags)
{
    if (rules._tdim == 0)
        rules._tdim = parent_tdim;
    if (rules._offset.empty())
        rules._offset.push_back(0);

    std::span<const T> ref_use = ref_vertices;
    std::span<const T> phys_use = physical_vertices;
    std::vector<T> ref_vertices_basix;
    std::vector<T> phys_vertices_basix;
    if (triangulate && !input_is_basix && !is_simplex(cell_type))
    {
        const auto perm = cell::vtk_to_basix_vertex_permutation(cell_type);
        ref_vertices_basix = cell::permute_vertex_data(ref_vertices, parent_tdim, perm);
        phys_vertices_basix = cell::permute_vertex_data(physical_vertices, gdim, perm);
        ref_use = std::span<const T>(ref_vertices_basix.data(), ref_vertices_basix.size());
        phys_use = std::span<const T>(phys_vertices_basix.data(), phys_vertices_basix.size());
    }

    const int entity_dim = cell::get_tdim(cell_type);
    const int nv = static_cast<int>(ref_use.size()) / parent_tdim;

    if (triangulate
        && cell_type == cell::type::quadrilateral
        && source_parent_type == cell::type::triangle
        && root_vertex_flags.size() == static_cast<std::size_t>(nv))
    {
        std::array<int, 4> quad_tokens = {};
        for (int i = 0; i < nv; ++i)
        {
            quad_tokens[static_cast<std::size_t>(i)]
                = root_vertex_flags[static_cast<std::size_t>(i)] ? i : 100 + i;
        }

        int next_ref_token = 200;
        auto ref_split = cell::quad_midpoint::split_triangle_derived_quadrilateral<T>(
            ref_use,
            parent_tdim,
            std::span<const int>(quad_tokens.data(), quad_tokens.size()),
            next_ref_token);

        int next_phys_token = 200;
        auto phys_split = cell::quad_midpoint::split_triangle_derived_quadrilateral<T>(
            phys_use,
            gdim,
            std::span<const int>(quad_tokens.data(), quad_tokens.size()),
            next_phys_token);

        std::array<int, 256> token_to_local;
        token_to_local.fill(-1);
        for (int i = 0; i < nv; ++i)
            token_to_local[quad_tokens[static_cast<std::size_t>(i)]] = i;
        for (std::size_t i = 0; i < ref_split.added_vertex_tokens.size(); ++i)
            token_to_local[ref_split.added_vertex_tokens[i]] = nv + static_cast<int>(i);

        std::vector<T> ref_all(ref_use.begin(), ref_use.end());
        ref_all.insert(
            ref_all.end(),
            ref_split.added_vertex_coords.begin(),
            ref_split.added_vertex_coords.end());

        std::vector<T> phys_all(phys_use.begin(), phys_use.end());
        phys_all.insert(
            phys_all.end(),
            phys_split.added_vertex_coords.begin(),
            phys_split.added_vertex_coords.end());

        std::vector<T> ref_simplex;
        std::vector<T> phys_simplex;
        for (const auto& tri_tokens : ref_split.triangles)
        {
            const std::array<int, 3> local_ids = {
                token_to_local[tri_tokens[0]],
                token_to_local[tri_tokens[1]],
                token_to_local[tri_tokens[2]],
            };
            gather_subcell_vertices(
                std::span<const T>(ref_all.data(), ref_all.size()),
                parent_tdim,
                std::span<const int>(local_ids.data(), local_ids.size()),
                ref_simplex);
            gather_subcell_vertices(
                std::span<const T>(phys_all.data(), phys_all.size()),
                gdim,
                std::span<const int>(local_ids.data(), local_ids.size()),
                phys_simplex);
            append_simplex_quadrature(
                rules,
                cell::type::triangle,
                std::span<const T>(ref_simplex.data(), ref_simplex.size()),
                std::span<const T>(phys_simplex.data(), phys_simplex.size()),
                parent_tdim,
                gdim,
                order);
        }

        rules._parent_map.push_back(parent_cell_id);
        rules._offset.push_back(static_cast<int32_t>(rules._weights.size()));
        return;
    }

    if (triangulate
        && cell_type == cell::type::prism
        && source_parent_type == cell::type::tetrahedron
        && root_vertex_flags.size() == static_cast<std::size_t>(nv))
    {
        std::array<int, 6> prism_tokens = {};
        for (int i = 0; i < nv; ++i)
        {
            prism_tokens[static_cast<std::size_t>(i)]
                = root_vertex_flags[static_cast<std::size_t>(i)] ? i : 100 + i;
        }

        int next_ref_token = 200;
        auto ref_split = cell::prism_midpoint::split_tetra_derived_prism<T>(
            ref_use,
            parent_tdim,
            std::span<const int>(prism_tokens.data(), prism_tokens.size()),
            next_ref_token);

        int next_phys_token = 200;
        auto phys_split = cell::prism_midpoint::split_tetra_derived_prism<T>(
            phys_use,
            gdim,
            std::span<const int>(prism_tokens.data(), prism_tokens.size()),
            next_phys_token);

        std::array<int, 256> token_to_local;
        token_to_local.fill(-1);
        for (int i = 0; i < nv; ++i)
            token_to_local[prism_tokens[static_cast<std::size_t>(i)]] = i;
        for (std::size_t i = 0; i < ref_split.added_vertex_tokens.size(); ++i)
            token_to_local[ref_split.added_vertex_tokens[i]] = nv + static_cast<int>(i);

        std::vector<T> ref_all(ref_use.begin(), ref_use.end());
        ref_all.insert(
            ref_all.end(),
            ref_split.added_vertex_coords.begin(),
            ref_split.added_vertex_coords.end());

        std::vector<T> phys_all(phys_use.begin(), phys_use.end());
        phys_all.insert(
            phys_all.end(),
            phys_split.added_vertex_coords.begin(),
            phys_split.added_vertex_coords.end());

        std::vector<T> ref_simplex;
        std::vector<T> phys_simplex;
        for (const auto& tet_tokens : ref_split.tets)
        {
            const std::array<int, 4> local_ids = {
                token_to_local[tet_tokens[0]],
                token_to_local[tet_tokens[1]],
                token_to_local[tet_tokens[2]],
                token_to_local[tet_tokens[3]],
            };
            gather_subcell_vertices(
                std::span<const T>(ref_all.data(), ref_all.size()),
                parent_tdim,
                std::span<const int>(local_ids.data(), local_ids.size()),
                ref_simplex);
            gather_subcell_vertices(
                std::span<const T>(phys_all.data(), phys_all.size()),
                gdim,
                std::span<const int>(local_ids.data(), local_ids.size()),
                phys_simplex);
            append_simplex_quadrature(
                rules,
                cell::type::tetrahedron,
                std::span<const T>(ref_simplex.data(), ref_simplex.size()),
                std::span<const T>(phys_simplex.data(), phys_simplex.size()),
                parent_tdim,
                gdim,
                order);
        }

        rules._parent_map.push_back(parent_cell_id);
        rules._offset.push_back(static_cast<int32_t>(rules._weights.size()));
        return;
    }

    if (triangulate && !is_simplex(cell_type) && entity_dim >= 2)
    {
        std::vector<int> local_ids(static_cast<std::size_t>(nv));
        std::iota(local_ids.begin(), local_ids.end(), 0);

        std::vector<std::vector<int>> simplices;
        cell::triangulation(cell_type, local_ids.data(), simplices);
        const auto simplex_type = simplex_type_for_dim(entity_dim);

        std::vector<T> ref_simplex;
        std::vector<T> phys_simplex;
        for (const auto& simplex : simplices)
        {
            gather_subcell_vertices(
                ref_use, parent_tdim,
                std::span<const int>(simplex.data(), simplex.size()),
                ref_simplex);
            gather_subcell_vertices(
                phys_use, gdim,
                std::span<const int>(simplex.data(), simplex.size()),
                phys_simplex);
            append_simplex_quadrature(
                rules,
                simplex_type,
                std::span<const T>(ref_simplex.data(), ref_simplex.size()),
                std::span<const T>(phys_simplex.data(), phys_simplex.size()),
                parent_tdim,
                gdim,
                order);
        }
    }
    else if (is_simplex(cell_type))
    {
        append_simplex_quadrature(
            rules,
            cell_type,
            ref_vertices,
            physical_vertices,
            parent_tdim,
            gdim,
            order);
    }
    else
    {
        const auto ref_rule = quadrature::get_reference_rule<T>(cell_type, order);
        const T measure = cell::affine_volume_factor<T>(
            cell_type, phys_use.data(), gdim);
        rules._points.insert(
            rules._points.end(), ref_rule._points.begin(), ref_rule._points.end());
        for (int q = 0; q < ref_rule._num_points; ++q)
            rules._weights.push_back(ref_rule._weights[q] * measure);
    }

    rules._parent_map.push_back(parent_cell_id);
    rules._offset.push_back(static_cast<int32_t>(rules._weights.size()));
}

template <std::floating_point T, std::integral I>
void append_cut_entities(mesh::CutMesh<T>& out,
                         quadrature::QuadratureRules<T>* rules,
                         const HOMeshPart<T, I>& part,
                         bool triangulate,
                         int quadrature_order)
{
    const auto& mesh = *part.bg->mesh;
    for (std::int32_t cut_id : part.cut_cell_ids)
    {
        const auto& adapt_cell = part.cut_cells->adapt_cells[static_cast<std::size_t>(cut_id)];
        const auto entities = selected_entities(part, adapt_cell, cut_id);
        if (entities.empty())
            continue;

        const I parent_cell_id = part.cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_id)];
        const auto parent_vertex_coords = parent_cell_vertex_coords_vtk(mesh, parent_cell_id);

        for (const auto& entity : entities)
        {
            std::vector<int> root_vertex_flags;
            if ((entity.type == cell::type::prism
                 && adapt_cell.parent_cell_type == cell::type::tetrahedron)
                || (entity.type == cell::type::quadrilateral
                    && adapt_cell.parent_cell_type == cell::type::triangle))
            {
                root_vertex_flags.resize(entity.vertices.size(), 0);
                for (std::size_t j = 0; j < entity.vertices.size(); ++j)
                {
                    const auto zm = adapt_cell.zero_mask_per_vertex[
                        static_cast<std::size_t>(entity.vertices[j])];
                    bool is_root = false;
                    if (part.expr.zero_required != 0)
                    {
                        is_root = (zm & part.expr.zero_required) == part.expr.zero_required;
                    }
                    else
                    {
                        is_root = zm != 0;
                    }
                    root_vertex_flags[j] = is_root ? 1 : 0;
                }
            }

            const auto ref_coords = entity_reference_coords(
                adapt_cell, std::span<const int>(entity.vertices.data(), entity.vertices.size()));
            const auto phys_coords = cell::push_forward_affine_map<T>(
                adapt_cell.parent_cell_type,
                parent_vertex_coords,
                mesh.gdim,
                std::span<const T>(ref_coords.data(), ref_coords.size()));

            append_mesh_entity(
                out,
                std::span<const T>(phys_coords.data(), phys_coords.size()),
                mesh.gdim,
                entity.type,
                static_cast<int>(parent_cell_id),
                triangulate,
                /*input_is_basix=*/true,
                adapt_cell.parent_cell_type,
                std::span<const int>(root_vertex_flags.data(), root_vertex_flags.size()));

            if (rules != nullptr)
            {
                append_entity_quadrature(
                    *rules,
                    entity.type,
                    std::span<const T>(ref_coords.data(), ref_coords.size()),
                    std::span<const T>(phys_coords.data(), phys_coords.size()),
                    mesh.tdim,
                    mesh.gdim,
                    static_cast<int>(parent_cell_id),
                    quadrature_order,
                    triangulate || !is_simplex(entity.type),
                    /*input_is_basix=*/true,
                    adapt_cell.parent_cell_type,
                    std::span<const int>(root_vertex_flags.data(), root_vertex_flags.size()));
            }
        }
    }
}

template <std::floating_point T, std::integral I>
void append_uncut_volume_cells(mesh::CutMesh<T>& out,
                               quadrature::QuadratureRules<T>* rules,
                               const HOMeshPart<T, I>& part,
                               bool triangulate,
                               int quadrature_order)
{
    const auto& mesh = *part.bg->mesh;
    if (part.dim != mesh.tdim)
        return;

    for (I cell_id : part.uncut_cell_ids)
    {
        const auto ctype = mesh.cell_type(cell_id);
        const auto phys_coords = parent_cell_vertex_coords_vtk(mesh, cell_id);
        append_mesh_entity(
            out,
            std::span<const T>(phys_coords.data(), phys_coords.size()),
            mesh.gdim,
            ctype,
            static_cast<int>(cell_id),
            triangulate,
            /*input_is_basix=*/false,
            cell::type::point,
            std::span<const int>());

        if (rules != nullptr)
        {
            const auto ref_coords = cell::canonical_vertices<T>(ctype);
            append_entity_quadrature(
                *rules,
                ctype,
                std::span<const T>(ref_coords.data(), ref_coords.size()),
                std::span<const T>(phys_coords.data(), phys_coords.size()),
                mesh.tdim,
                mesh.gdim,
                static_cast<int>(cell_id),
                quadrature_order,
                triangulate,
                /*input_is_basix=*/false,
                cell::type::point,
                std::span<const int>());
        }
    }
}

} // namespace

template <std::floating_point T, std::integral I>
mesh::CutMesh<T> visualization_mesh(const HOMeshPart<T, I>& part,
                                    bool include_uncut_cells,
                                    bool triangulate)
{
    if (!part.cut_cells || !part.bg || !part.bg->mesh)
        throw std::runtime_error("HOMeshPart is not attached to cut-cell storage");

    mesh::CutMesh<T> out;
    out._gdim = part.bg->mesh->gdim;
    out._tdim = part.dim;
    out._offset.push_back(0);

    append_cut_entities(
        out,
        static_cast<quadrature::QuadratureRules<T>*>(nullptr),
        part,
        triangulate,
        /*quadrature_order=*/0);
    if (include_uncut_cells)
    {
        append_uncut_volume_cells(
            out,
            static_cast<quadrature::QuadratureRules<T>*>(nullptr),
            part,
            triangulate,
            /*quadrature_order=*/0);
    }

    return out;
}

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> quadrature_rules(const HOMeshPart<T, I>& part,
                                                int order,
                                                bool include_uncut_cells,
                                                bool triangulate)
{
    if (!part.cut_cells || !part.bg || !part.bg->mesh)
        throw std::runtime_error("HOMeshPart is not attached to cut-cell storage");
    if (part.bg->num_level_sets != 1)
    {
        throw std::runtime_error(
            "HOMeshPart output currently supports exactly one level set");
    }

    mesh::CutMesh<T> unused_mesh;
    quadrature::QuadratureRules<T> rules;
    rules._offset.push_back(0);

    append_cut_entities(unused_mesh, &rules, part, triangulate, order);
    if (include_uncut_cells)
        append_uncut_volume_cells(unused_mesh, &rules, part, triangulate, order);

    return rules;
}

template mesh::CutMesh<double> visualization_mesh(
    const HOMeshPart<double, int>&, bool, bool);
template mesh::CutMesh<float> visualization_mesh(
    const HOMeshPart<float, int>&, bool, bool);
template mesh::CutMesh<double> visualization_mesh(
    const HOMeshPart<double, long>&, bool, bool);
template mesh::CutMesh<float> visualization_mesh(
    const HOMeshPart<float, long>&, bool, bool);

template quadrature::QuadratureRules<double> quadrature_rules(
    const HOMeshPart<double, int>&, int, bool, bool);
template quadrature::QuadratureRules<float> quadrature_rules(
    const HOMeshPart<float, int>&, int, bool, bool);
template quadrature::QuadratureRules<double> quadrature_rules(
    const HOMeshPart<double, long>&, int, bool, bool);
template quadrature::QuadratureRules<float> quadrature_rules(
    const HOMeshPart<float, long>&, int, bool, bool);

} // namespace cutcells::output
