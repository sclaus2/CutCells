// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "ho_mesh_part_output.h"

#include "algoim_quadrature.h"
#include "cell_topology.h"
#include "mapping.h"
#include "quadrature_tables.h"
#include "reference_cell.h"
#include "triangulation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace cutcells::output
{
namespace
{

template <std::floating_point T, std::integral I>
std::vector<T> parent_cell_vertex_coords_vtk(const MeshView<T, I>& mesh,
                                             I cell_id)
{
    const auto ctype = mesh.cell_type(cell_id);
    const int nv = cell::get_num_vertices(ctype);
    std::vector<T> coords(static_cast<std::size_t>(nv * mesh.gdim), T(0));
    std::vector<I> cell_node_scratch;
    const auto parent_nodes = mesh.cell_nodes(cell_id, cell_node_scratch);

    for (int vtk_v = 0; vtk_v < nv; ++vtk_v)
    {
        const int local_v = mesh.vtk_vertex_order
                                ? vtk_v
                                : cell::vtk_to_basix_vertex(ctype, vtk_v);
        const I node_id = parent_nodes[static_cast<std::size_t>(local_v)];
        const T* x = mesh.node(node_id);
        for (int d = 0; d < mesh.gdim; ++d)
            coords[static_cast<std::size_t>(vtk_v * mesh.gdim + d)] = x[d];
    }
    return coords;
}

template <std::floating_point T, std::integral I>
std::vector<T> parent_cell_vertex_coords_basix(const MeshView<T, I>& mesh,
                                               I cell_id)
{
    const auto ctype = mesh.cell_type(cell_id);
    const int nv = cell::get_num_vertices(ctype);
    std::vector<T> coords(static_cast<std::size_t>(nv * mesh.gdim), T(0));
    std::vector<I> cell_node_scratch;
    const auto parent_nodes = mesh.cell_nodes(cell_id, cell_node_scratch);

    for (int basix_v = 0; basix_v < nv; ++basix_v)
    {
        const int local_v = mesh.vtk_vertex_order
                                ? cell::basix_to_vtk_vertex(ctype, basix_v)
                                : basix_v;
        const I node_id = parent_nodes[static_cast<std::size_t>(local_v)];
        const T* x = mesh.node(node_id);
        for (int d = 0; d < mesh.gdim; ++d)
            coords[static_cast<std::size_t>(basix_v * mesh.gdim + d)] = x[d];
    }
    return coords;
}

inline bool is_simplex(cell::type cell_type)
{
    using cell::type;
    return cell_type == type::point
        || cell_type == type::interval
        || cell_type == type::triangle
        || cell_type == type::tetrahedron;
}

inline cell::type simplex_type_for_dim(int dim)
{
    using cell::type;
    switch (dim)
    {
    case 0:
        return type::point;
    case 1:
        return type::interval;
    case 2:
        return type::triangle;
    case 3:
        return type::tetrahedron;
    default:
        throw std::runtime_error("Unsupported simplex dimension");
    }
}

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
bool vertex_state_for_level_set(const AdaptCell<T>& ac,
                                int vertex_id,
                                int level_set_id,
                                bool& is_negative,
                                bool& is_positive,
                                bool& is_zero)
{
    if (level_set_id < 0 || level_set_id >= 64)
        return false;

    const std::uint64_t bit = std::uint64_t(1) << level_set_id;
    const auto zm = ac.zero_mask_per_vertex[static_cast<std::size_t>(vertex_id)];
    const auto nm = ac.negative_mask_per_vertex[static_cast<std::size_t>(vertex_id)];
    is_zero = (zm & bit) != 0;
    is_negative = !is_zero && ((nm & bit) != 0);
    is_positive = !is_zero && !is_negative;
    return true;
}

template <std::floating_point T, std::integral I>
bool leaf_cell_matches_sign_requirements(
    const AdaptCell<T>& ac,
    std::span<const std::int32_t> cell_verts,
    const SelectionTerm& term,
    std::uint64_t cut_cell_active_mask,
    const ParentCellClassification<T, I>& parent_cells,
    I parent_cell_id)
{
    const int nls = std::min(parent_cells.num_level_sets, 64);
    for (int li = 0; li < nls; ++li)
    {
        const std::uint64_t bit = std::uint64_t(1) << li;
        const bool require_neg = (term.negative_required & bit) != 0;
        const bool require_pos = (term.positive_required & bit) != 0;
        if (!require_neg && !require_pos)
            continue;

        const bool ls_is_active = (cut_cell_active_mask & bit) != 0;
        if (!ls_is_active)
        {
            const auto dom = parent_cells.domain(li, parent_cell_id);
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
            vertex_state_for_level_set(
                ac, static_cast<int>(cv), li, is_neg, is_pos, is_zero);
            has_neg = has_neg || is_neg;
            has_pos = has_pos || is_pos;
        }

        if (require_neg && (has_pos || !has_neg))
            return false;
        if (require_pos && (has_neg || !has_pos))
            return false;
    }

    return true;
}

template <std::floating_point T, std::integral I>
bool leaf_cell_matches_selection_expr(
    const AdaptCell<T>& ac,
    std::span<const std::int32_t> cell_verts,
    const SelectionExpr& expr,
    std::uint64_t cut_cell_active_mask,
    const ParentCellClassification<T, I>& parent_cells,
    I parent_cell_id)
{
    for (const auto& term : expr.terms)
    {
        if (leaf_cell_matches_sign_requirements(
                ac, cell_verts, term, cut_cell_active_mask, parent_cells, parent_cell_id))
        {
            return true;
        }
    }
    return false;
}

template <std::floating_point T, std::integral I>
bool zero_entity_matches(
    const AdaptCell<T>& ac,
    int zero_entity_index,
    int target_dim,
    const SelectionTerm& term,
    std::uint64_t cut_cell_active_mask,
    const ParentCellClassification<T, I>& parent_cells,
    I parent_cell_id)
{
    if (ac.zero_entity_dim[static_cast<std::size_t>(zero_entity_index)]
        != target_dim)
    {
        return false;
    }

    const auto zero_mask =
        ac.zero_entity_zero_mask[static_cast<std::size_t>(zero_entity_index)];
    if ((cut_cell_active_mask & term.zero_required) != term.zero_required)
        return false;

    if ((zero_mask & term.zero_required) != term.zero_required)
        return false;

    if (term.negative_required == 0 && term.positive_required == 0)
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
        if (!cell_contains_all_vertices(
                cell_verts,
                std::span<const int>(zero_verts.data(), zero_verts.size())))
        {
            continue;
        }

        if (leaf_cell_matches_sign_requirements(
                ac, cell_verts, term, cut_cell_active_mask, parent_cells, parent_cell_id))
        {
            return true;
        }
    }

    return false;
}

template <std::floating_point T, std::integral I>
bool zero_entity_matches_selection_expr(
    const AdaptCell<T>& ac,
    int zero_entity_index,
    int target_dim,
    const SelectionExpr& expr,
    std::uint64_t cut_cell_active_mask,
    const ParentCellClassification<T, I>& parent_cells,
    I parent_cell_id)
{
    for (const auto& term : expr.terms)
    {
        if (zero_entity_matches(
                ac, zero_entity_index, target_dim, term,
                cut_cell_active_mask, parent_cells, parent_cell_id))
        {
            return true;
        }
    }
    return false;
}

struct SelectedEntity
{
    int local_zero_entity_id = -1;
    cell::type type = cell::type::point;
    std::vector<int> vertices;
};

template <std::floating_point T, std::integral I>
std::vector<SelectedEntity> selected_entities(const HOMeshPart<T, I>& part,
                                              const AdaptCell<T>& adapt_cell,
                                              int cut_cell_id)
{
    const I parent_cell_id =
        part.cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_cell_id)];
    const std::uint64_t cut_active_mask =
        part.cut_cells->active_level_set_mask[static_cast<std::size_t>(cut_cell_id)];

    std::vector<SelectedEntity> entities;
    if (part.dim == adapt_cell.tdim)
    {
        const int n_cells = adapt_cell.n_entities(adapt_cell.tdim);
        entities.reserve(static_cast<std::size_t>(n_cells));
        for (int c = 0; c < n_cells; ++c)
        {
            auto verts = adapt_cell.entity_to_vertex[adapt_cell.tdim][
                static_cast<std::int32_t>(c)];
            if (!leaf_cell_matches_selection_expr(
                    adapt_cell, verts, part.expr, cut_active_mask,
                    *part.parent_cells, parent_cell_id))
            {
                continue;
            }

            SelectedEntity entity;
            entity.type = adapt_cell.entity_types[adapt_cell.tdim][
                static_cast<std::size_t>(c)];
            entity.vertices.assign(verts.begin(), verts.end());
            entities.push_back(std::move(entity));
        }
        return entities;
    }

    const int n_zero = adapt_cell.n_zero_entities();
    entities.reserve(static_cast<std::size_t>(n_zero));
    for (int z = 0; z < n_zero; ++z)
    {
        if (!zero_entity_matches_selection_expr(
                adapt_cell, z, part.dim, part.expr, cut_active_mask,
                *part.parent_cells, parent_cell_id))
        {
            continue;
        }

        const int zdim = adapt_cell.zero_entity_dim[static_cast<std::size_t>(z)];
        const int zid = adapt_cell.zero_entity_id[static_cast<std::size_t>(z)];

        SelectedEntity entity;
        entity.local_zero_entity_id = z;
        if (zdim == 0)
        {
            entity.type = cell::type::point;
            entity.vertices.push_back(zid);
        }
        else
        {
            entity.type = adapt_cell.entity_types[zdim][static_cast<std::size_t>(zid)];
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
                adapt_cell.vertex_coords[static_cast<std::size_t>(
                    gv * adapt_cell.tdim + d)];
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
void tabulate_tensor_product_p1(cell::type cell_type,
                                const T* X,
                                std::vector<T>& shape,
                                std::vector<T>& dshape)
{
    if (cell_type == cell::type::quadrilateral)
    {
        const T u = X[0];
        const T v = X[1];
        shape = {
            (T(1) - u) * (T(1) - v),
            u * (T(1) - v),
            (T(1) - u) * v,
            u * v};
        dshape = {
            -(T(1) - v), -(T(1) - u),
             (T(1) - v), -u,
            -v,            T(1) - u,
             v,            u};
        return;
    }

    if (cell_type == cell::type::hexahedron)
    {
        const T u = X[0];
        const T v = X[1];
        const T w = X[2];
        const std::array<T, 2> nx = {T(1) - u, u};
        const std::array<T, 2> ny = {T(1) - v, v};
        const std::array<T, 2> nz = {T(1) - w, w};
        const std::array<T, 2> dnx = {-T(1), T(1)};
        const std::array<T, 2> dny = {-T(1), T(1)};
        const std::array<T, 2> dnz = {-T(1), T(1)};
        const std::array<std::array<int, 3>, 8> bits = {{
            {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
            {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}}};

        shape.assign(8, T(0));
        dshape.assign(8 * 3, T(0));
        for (std::size_t i = 0; i < bits.size(); ++i)
        {
            const int bx = bits[i][0];
            const int by = bits[i][1];
            const int bz = bits[i][2];
            shape[i] = nx[static_cast<std::size_t>(bx)]
                     * ny[static_cast<std::size_t>(by)]
                     * nz[static_cast<std::size_t>(bz)];
            dshape[i * 3 + 0] = dnx[static_cast<std::size_t>(bx)]
                              * ny[static_cast<std::size_t>(by)]
                              * nz[static_cast<std::size_t>(bz)];
            dshape[i * 3 + 1] = nx[static_cast<std::size_t>(bx)]
                              * dny[static_cast<std::size_t>(by)]
                              * nz[static_cast<std::size_t>(bz)];
            dshape[i * 3 + 2] = nx[static_cast<std::size_t>(bx)]
                              * ny[static_cast<std::size_t>(by)]
                              * dnz[static_cast<std::size_t>(bz)];
        }
        return;
    }

    throw std::invalid_argument(
        "tabulate_tensor_product_p1: unsupported cell type");
}

template <std::floating_point T>
T gram_measure(std::span<const T> jacobian,
               int entity_dim,
               int gdim)
{
    std::array<T, 9> gram = {};
    for (int i = 0; i < entity_dim; ++i)
    {
        for (int j = 0; j < entity_dim; ++j)
        {
            T sum = T(0);
            for (int d = 0; d < gdim; ++d)
                sum += jacobian[static_cast<std::size_t>(i * gdim + d)]
                     * jacobian[static_cast<std::size_t>(j * gdim + d)];
            gram[static_cast<std::size_t>(i * entity_dim + j)] = sum;
        }
    }

    T det = T(0);
    if (entity_dim == 1)
        det = gram[0];
    else if (entity_dim == 2)
        det = gram[0] * gram[3] - gram[1] * gram[2];
    else if (entity_dim == 3)
    {
        det =
            gram[0] * (gram[4] * gram[8] - gram[7] * gram[5])
          - gram[3] * (gram[1] * gram[8] - gram[7] * gram[2])
          + gram[6] * (gram[1] * gram[5] - gram[4] * gram[2]);
    }
    else
    {
        throw std::invalid_argument("gram_measure: unsupported entity dimension");
    }

    if (det < T(0) && std::abs(det) < T(100) * std::numeric_limits<T>::epsilon())
        det = T(0);
    if (det < T(0))
        throw std::runtime_error("gram_measure: negative Gram determinant");
    return std::sqrt(det);
}

template <std::floating_point T>
void append_tensor_product_quadrature(quadrature::QuadratureRules<T>& rules,
                                      cell::type cell_type,
                                      std::span<const T> ref_vertices,
                                      std::span<const T> physical_vertices,
                                      int parent_tdim,
                                      int gdim,
                                      int order)
{
    const auto ref_rule = quadrature::get_reference_rule<T>(cell_type, order);
    const int num_points = ref_rule._num_points;
    const int entity_dim = ref_rule._tdim;
    const int num_vertices = cell::get_num_vertices(cell_type);

    std::vector<T> shape;
    std::vector<T> dshape;
    std::vector<T> mapped_ref_point(static_cast<std::size_t>(parent_tdim), T(0));
    std::vector<T> jacobian(static_cast<std::size_t>(entity_dim * gdim), T(0));

    for (int q = 0; q < num_points; ++q)
    {
        const T* X = ref_rule._points.data()
            + static_cast<std::size_t>(q * entity_dim);
        tabulate_tensor_product_p1(cell_type, X, shape, dshape);

        std::fill(mapped_ref_point.begin(), mapped_ref_point.end(), T(0));
        std::fill(jacobian.begin(), jacobian.end(), T(0));

        for (int v = 0; v < num_vertices; ++v)
        {
            for (int d = 0; d < parent_tdim; ++d)
            {
                mapped_ref_point[static_cast<std::size_t>(d)] +=
                    shape[static_cast<std::size_t>(v)]
                    * ref_vertices[static_cast<std::size_t>(v * parent_tdim + d)];
            }

            for (int k = 0; k < entity_dim; ++k)
            {
                const T dN =
                    dshape[static_cast<std::size_t>(v * entity_dim + k)];
                for (int d = 0; d < gdim; ++d)
                {
                    jacobian[static_cast<std::size_t>(k * gdim + d)] +=
                        dN * physical_vertices[
                            static_cast<std::size_t>(v * gdim + d)];
                }
            }
        }

        rules._points.insert(
            rules._points.end(), mapped_ref_point.begin(), mapped_ref_point.end());
        rules._weights.push_back(
            ref_rule._weights[static_cast<std::size_t>(q)]
            * gram_measure(std::span<const T>(jacobian.data(), jacobian.size()),
                           entity_dim, gdim));
    }
}

template <std::floating_point T>
void append_mesh_entity(mesh::CutMesh<T>& out,
                        std::span<const T> physical_coords,
                        int gdim,
                        cell::type cell_type,
                        int parent_cell_id,
                        bool input_is_basix)
{
    if (out._gdim == 0)
        out._gdim = gdim;
    if (out._tdim == 0)
        out._tdim = cell::get_tdim(cell_type);

    std::vector<T> reordered_coords;
    std::span<const T> output_coords = physical_coords;
    if (!input_is_basix && !is_simplex(cell_type))
    {
        const auto perm = cell::vtk_to_basix_vertex_permutation(cell_type);
        reordered_coords = cell::permute_vertex_data(physical_coords, gdim, perm);
        output_coords = std::span<const T>(
            reordered_coords.data(), reordered_coords.size());
    }

    const int nv = static_cast<int>(output_coords.size()) / gdim;
    const int vertex_base = out._num_vertices;
    out._vertex_coords.insert(
        out._vertex_coords.end(), output_coords.begin(), output_coords.end());
    out._num_vertices += nv;

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
                              int order)
{
    if (rules._tdim == 0)
        rules._tdim = parent_tdim;
    if (rules._offset.empty())
        rules._offset.push_back(0);

    if (cell_type == cell::type::point)
    {
        rules._points.insert(
            rules._points.end(),
            ref_vertices.begin(),
            ref_vertices.begin() + parent_tdim);
        rules._weights.push_back(T(1));
        rules._parent_map.push_back(parent_cell_id);
        rules._offset.push_back(static_cast<std::int32_t>(rules._weights.size()));
        return;
    }

    const int entity_dim = cell::get_tdim(cell_type);
    if (is_simplex(cell_type))
    {
        append_simplex_quadrature(
            rules, cell_type, ref_vertices, physical_vertices,
            parent_tdim, gdim, order);
    }
    else if ((cell_type == cell::type::quadrilateral
              || cell_type == cell::type::hexahedron)
             && entity_dim == parent_tdim)
    {
        append_tensor_product_quadrature(
            rules, cell_type, ref_vertices, physical_vertices,
            parent_tdim, gdim, order);
    }
    else
    {
        std::vector<int> local_ids(
            static_cast<std::size_t>(ref_vertices.size() / parent_tdim));
        std::iota(local_ids.begin(), local_ids.end(), 0);
        std::vector<std::vector<int>> simplices;
        cell::triangulation(cell_type, local_ids.data(), simplices);
        const auto simplex_type = simplex_type_for_dim(entity_dim);

        std::vector<T> ref_simplex;
        std::vector<T> phys_simplex;
        for (const auto& simplex : simplices)
        {
            gather_subcell_vertices(
                ref_vertices, parent_tdim,
                std::span<const int>(simplex.data(), simplex.size()),
                ref_simplex);
            gather_subcell_vertices(
                physical_vertices, gdim,
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

    rules._parent_map.push_back(parent_cell_id);
    rules._offset.push_back(static_cast<std::int32_t>(rules._weights.size()));
}

template <std::floating_point T, std::integral I>
void append_cut_entities(mesh::CutMesh<T>& out,
                         quadrature::QuadratureRules<T>* rules,
                         const HOMeshPart<T, I>& part,
                         int quadrature_order)
{
    const auto& mesh = *part.mesh;
    for (std::int32_t cut_id : part.cut_cell_ids)
    {
        const auto& adapt_cell =
            part.cut_cells->adapt_cells[static_cast<std::size_t>(cut_id)];
        const auto entities = selected_entities(part, adapt_cell, cut_id);
        if (entities.empty())
            continue;

        const I parent_cell_id =
            part.cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_id)];
        const auto parent_vertex_coords =
            parent_cell_vertex_coords_basix(mesh, parent_cell_id);

        for (const auto& entity : entities)
        {
            const auto ref_coords = entity_reference_coords(
                adapt_cell,
                std::span<const int>(entity.vertices.data(), entity.vertices.size()));
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
                /*input_is_basix=*/true);

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
                    quadrature_order);
            }
        }
    }
}

template <std::floating_point T, std::integral I>
void append_uncut_volume_cells(mesh::CutMesh<T>& out,
                               quadrature::QuadratureRules<T>* rules,
                               const HOMeshPart<T, I>& part,
                               int quadrature_order)
{
    const auto& mesh = *part.mesh;
    if (part.dim != mesh.tdim)
        return;

    for (I cell_id : part.uncut_cell_ids)
    {
        const auto ctype = mesh.cell_type(cell_id);
        const auto vtk_phys_coords = parent_cell_vertex_coords_vtk(mesh, cell_id);
        append_mesh_entity(
            out,
            std::span<const T>(vtk_phys_coords.data(), vtk_phys_coords.size()),
            mesh.gdim,
            ctype,
            static_cast<int>(cell_id),
            /*input_is_basix=*/false);

        if (rules != nullptr)
        {
            const auto basix_phys_coords =
                parent_cell_vertex_coords_basix(mesh, cell_id);
            const auto ref_coords = cell::canonical_vertices<T>(ctype);
            append_entity_quadrature(
                *rules,
                ctype,
                std::span<const T>(ref_coords.data(), ref_coords.size()),
                std::span<const T>(basix_phys_coords.data(), basix_phys_coords.size()),
                mesh.tdim,
                mesh.gdim,
                static_cast<int>(cell_id),
                quadrature_order);
        }
    }
}

template <std::floating_point T>
T reference_cell_measure(cell::type cell_type)
{
    const auto rule = quadrature::get_reference_rule<T>(cell_type, 1);
    return std::accumulate(rule._weights.begin(), rule._weights.end(), T(0));
}

template <std::floating_point T>
T entity_reference_measure(cell::type cell_type,
                           std::span<const T> reference_vertices,
                           int parent_tdim)
{
    if (cell_type == cell::type::point)
        return T(1);

    return reference_cell_measure<T>(cell_type)
         * cell::affine_volume_factor<T>(
             cell_type, reference_vertices.data(), parent_tdim);
}

} // namespace

template <std::floating_point T, std::integral I>
std::vector<SelectedZeroEntityInfo> selected_zero_entity_infos(
    const HOMeshPart<T, I>& part)
{
    if (!part.cut_cells || !part.parent_cells)
        throw std::runtime_error("HOMeshPart is not attached to cut-cell storage");

    std::vector<SelectedZeroEntityInfo> infos;
    if (part.dim == part.cut_cells->tdim)
        return infos;

    for (std::int32_t cut_id : part.cut_cell_ids)
    {
        const auto& ac =
            part.cut_cells->adapt_cells[static_cast<std::size_t>(cut_id)];
        const I parent_cell_id =
            part.cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_id)];
        const std::uint64_t cut_active_mask =
            part.cut_cells->active_level_set_mask[static_cast<std::size_t>(cut_id)];

        const int n_zero = ac.n_zero_entities();
        for (int z = 0; z < n_zero; ++z)
        {
            if (!zero_entity_matches_selection_expr(
                    ac, z, part.dim, part.expr, cut_active_mask,
                    *part.parent_cells, parent_cell_id))
            {
                continue;
            }

            infos.push_back({
                cut_id,
                static_cast<std::int32_t>(parent_cell_id),
                static_cast<std::int32_t>(z),
                static_cast<std::int32_t>(part.dim)});
        }
    }

    return infos;
}

template <std::floating_point T, std::integral I>
mesh::CutMesh<T> visualization_mesh(const HOMeshPart<T, I>& part,
                                    bool include_uncut_cells)
{
    if (!part.cut_cells || !part.parent_cells || !part.mesh)
        throw std::runtime_error("HOMeshPart is not attached to cut-cell storage");

    mesh::CutMesh<T> out;
    out._gdim = part.mesh->gdim;
    out._tdim = part.dim;
    out._offset.push_back(0);

    append_cut_entities(
        out,
        static_cast<quadrature::QuadratureRules<T>*>(nullptr),
        part,
        /*quadrature_order=*/0);
    if (include_uncut_cells)
    {
        append_uncut_volume_cells(
            out,
            static_cast<quadrature::QuadratureRules<T>*>(nullptr),
            part,
            /*quadrature_order=*/0);
    }

    return out;
}

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> quadrature_rules(const HOMeshPart<T, I>& part,
                                                int order,
                                                bool include_uncut_cells)
{
    if (!part.mesh)
        throw std::runtime_error("HOMeshPart is not attached to a mesh");
    if ((!part.cut_cells || !part.parent_cells) && !part.cut_cell_ids.empty())
        throw std::runtime_error("HOMeshPart cut entities require cut-cell storage");

    mesh::CutMesh<T> unused_mesh;
    quadrature::QuadratureRules<T> rules;
    rules._offset.push_back(0);

    append_cut_entities(unused_mesh, &rules, part, order);
    if (include_uncut_cells)
        append_uncut_volume_cells(unused_mesh, &rules, part, order);

    return rules;
}

QuadratureBackend quadrature_backend_from_string(std::string_view backend)
{
    if (backend == "straight" || backend == "cutcells")
        return QuadratureBackend::Straight;
    if (backend == "algoim" || backend == "algoim_bernstein")
        return QuadratureBackend::AlgoimBernstein;
    if (backend == "algoim_general")
        return QuadratureBackend::AlgoimGeneral;

    throw std::runtime_error(
        "Unknown quadrature backend '" + std::string(backend)
        + "'. Expected 'straight', 'algoim', or 'algoim_general'.");
}

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> quadrature_rules(const HOMeshPart<T, I>& part,
                                                int order,
                                                bool include_uncut_cells,
                                                QuadratureBackend backend)
{
    switch (backend)
    {
    case QuadratureBackend::Straight:
        return quadrature_rules(part, order, include_uncut_cells);
    case QuadratureBackend::AlgoimBernstein:
        return algoim_quadrature_rules(part, order, include_uncut_cells);
    case QuadratureBackend::AlgoimGeneral:
        return algoim_general_quadrature_rules(part, order, include_uncut_cells);
    }

    throw std::runtime_error("Unhandled quadrature backend");
}

template <std::floating_point T, std::integral I>
std::vector<std::pair<std::string, quadrature::QuadratureRules<T>>>
paired_quadrature_rules(
    const std::vector<std::pair<std::string, HOMeshPart<T, I>>>& parts,
    int order,
    bool include_uncut_cells,
    QuadratureBackend backend)
{
    switch (backend)
    {
    case QuadratureBackend::Straight:
    case QuadratureBackend::AlgoimGeneral:
    {
        std::vector<std::pair<std::string, quadrature::QuadratureRules<T>>> out;
        out.reserve(parts.size());
        for (const auto& [name, part] : parts)
        {
            out.emplace_back(
                name,
                quadrature_rules(part, order, include_uncut_cells, backend));
        }
        return out;
    }
    case QuadratureBackend::AlgoimBernstein:
        return algoim_paired_quadrature_rules(
            parts, order, include_uncut_cells);
    }

    throw std::runtime_error("Unhandled quadrature backend");
}

template <std::floating_point T, std::integral I>
std::pair<std::vector<I>, std::vector<T>>
volume_fractions(const HOMeshPart<T, I>& part)
{
    if (!part.mesh)
        throw std::runtime_error("HOMeshPart is not attached to a mesh");
    if ((!part.cut_cells || !part.parent_cells) && !part.cut_cell_ids.empty())
        throw std::runtime_error("HOMeshPart cut entities require cut-cell storage");
    if (part.dim != part.mesh->tdim)
        throw std::runtime_error("volume_fractions expects a volume HOMeshPart");

    std::unordered_map<I, T> fraction_by_parent;

    if (part.cut_cells)
    {
        for (std::int32_t cut_id : part.cut_cell_ids)
        {
            const auto& adapt_cell =
                part.cut_cells->adapt_cells[static_cast<std::size_t>(cut_id)];
            const I parent_cell_id =
                part.cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_id)];
            const auto entities = selected_entities(part, adapt_cell, cut_id);
            if (entities.empty())
                continue;

            const T parent_measure =
                reference_cell_measure<T>(adapt_cell.parent_cell_type);
            if (parent_measure <= T(0))
                throw std::runtime_error(
                    "volume_fractions encountered a parent cell with zero "
                    "reference measure");

            T selected_fraction = T(0);
            for (const auto& entity : entities)
            {
                if (cell::get_tdim(entity.type) != part.mesh->tdim)
                {
                    throw std::runtime_error(
                        "volume_fractions encountered a non-volume selected "
                        "entity in a volume HOMeshPart");
                }

                const auto ref_coords = entity_reference_coords(
                    adapt_cell,
                    std::span<const int>(entity.vertices.data(),
                                         entity.vertices.size()));
                selected_fraction += entity_reference_measure<T>(
                    entity.type,
                    std::span<const T>(ref_coords.data(), ref_coords.size()),
                    adapt_cell.tdim)
                                     / parent_measure;
            }

            fraction_by_parent[parent_cell_id] += selected_fraction;
        }
    }

    for (I cell_id : part.uncut_cell_ids)
        fraction_by_parent[cell_id] += T(1);

    std::vector<I> parents;
    parents.reserve(fraction_by_parent.size());
    for (const auto& [parent, fraction] : fraction_by_parent)
        parents.push_back(parent);
    std::ranges::sort(parents);

    std::vector<T> fractions;
    fractions.reserve(parents.size());
    for (I parent : parents)
    {
        const T fraction = fraction_by_parent.at(parent);
        fractions.push_back(std::clamp(fraction, T(0), T(1)));
    }

    return {std::move(parents), std::move(fractions)};
}

template mesh::CutMesh<double> visualization_mesh(
    const HOMeshPart<double, int>&, bool);
template mesh::CutMesh<float> visualization_mesh(
    const HOMeshPart<float, int>&, bool);
template mesh::CutMesh<double> visualization_mesh(
    const HOMeshPart<double, long>&, bool);
template mesh::CutMesh<float> visualization_mesh(
    const HOMeshPart<float, long>&, bool);

template std::vector<SelectedZeroEntityInfo> selected_zero_entity_infos(
    const HOMeshPart<double, int>&);
template std::vector<SelectedZeroEntityInfo> selected_zero_entity_infos(
    const HOMeshPart<float, int>&);
template std::vector<SelectedZeroEntityInfo> selected_zero_entity_infos(
    const HOMeshPart<double, long>&);
template std::vector<SelectedZeroEntityInfo> selected_zero_entity_infos(
    const HOMeshPart<float, long>&);

template quadrature::QuadratureRules<double> quadrature_rules(
    const HOMeshPart<double, int>&, int, bool);
template quadrature::QuadratureRules<float> quadrature_rules(
    const HOMeshPart<float, int>&, int, bool);
template quadrature::QuadratureRules<double> quadrature_rules(
    const HOMeshPart<double, long>&, int, bool);
template quadrature::QuadratureRules<float> quadrature_rules(
    const HOMeshPart<float, long>&, int, bool);

template quadrature::QuadratureRules<double> quadrature_rules(
    const HOMeshPart<double, int>&, int, bool, QuadratureBackend);
template quadrature::QuadratureRules<float> quadrature_rules(
    const HOMeshPart<float, int>&, int, bool, QuadratureBackend);
template quadrature::QuadratureRules<double> quadrature_rules(
    const HOMeshPart<double, long>&, int, bool, QuadratureBackend);
template quadrature::QuadratureRules<float> quadrature_rules(
    const HOMeshPart<float, long>&, int, bool, QuadratureBackend);

template std::vector<std::pair<std::string, quadrature::QuadratureRules<double>>>
paired_quadrature_rules(
    const std::vector<std::pair<std::string, HOMeshPart<double, int>>>&,
    int, bool, QuadratureBackend);
template std::vector<std::pair<std::string, quadrature::QuadratureRules<float>>>
paired_quadrature_rules(
    const std::vector<std::pair<std::string, HOMeshPart<float, int>>>&,
    int, bool, QuadratureBackend);
template std::vector<std::pair<std::string, quadrature::QuadratureRules<double>>>
paired_quadrature_rules(
    const std::vector<std::pair<std::string, HOMeshPart<double, long>>>&,
    int, bool, QuadratureBackend);
template std::vector<std::pair<std::string, quadrature::QuadratureRules<float>>>
paired_quadrature_rules(
    const std::vector<std::pair<std::string, HOMeshPart<float, long>>>&,
    int, bool, QuadratureBackend);

template std::pair<std::vector<int>, std::vector<double>> volume_fractions(
    const HOMeshPart<double, int>&);
template std::pair<std::vector<int>, std::vector<float>> volume_fractions(
    const HOMeshPart<float, int>&);
template std::pair<std::vector<long>, std::vector<double>> volume_fractions(
    const HOMeshPart<double, long>&);
template std::pair<std::vector<long>, std::vector<float>> volume_fractions(
    const HOMeshPart<float, long>&);

} // namespace cutcells::output
