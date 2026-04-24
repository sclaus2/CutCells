// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "adapt_cell.h"
#include "cell_topology.h"
#include "reference_cell.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <map>
#include <stdexcept>

namespace cutcells
{

    template <std::floating_point T>
void build_faces(AdaptCell<T>& ac)
{
    if (ac.tdim < 3)
        return;

    // Clear existing 2D entity pool.
    ac.entity_types[2].clear();
    ac.entity_to_vertex[2].offsets.clear();
    ac.entity_to_vertex[2].indices.clear();
    ac.entity_to_vertex[2].offsets.push_back(std::int32_t(0));

    // Invalidate face cert tag storage: face count is about to change and
    // faces may be reordered, so all old tags are stale.
    ac.face_cert_tag.clear();
    ac.face_cert_tag_num_level_sets = 0;

    // Track unique faces by sorted vertex set → face_id.
    std::map<std::vector<std::int32_t>, int> face_map;

    const int n_cells = ac.n_entities(ac.tdim);
    for (int c = 0; c < n_cells; ++c)
    {
        const cell::type ctype = ac.entity_types[ac.tdim][static_cast<std::size_t>(c)];
        auto cell_verts = ac.entity_to_vertex[ac.tdim][static_cast<std::int32_t>(c)];
        const int nf = cell::num_faces(ctype);

        for (int fi = 0; fi < nf; ++fi)
        {
            auto local_fv = cell::face_vertices(ctype, fi);
            const int fsize = static_cast<int>(local_fv.size());

            // Build global-vertex face and sorted key
            std::vector<std::int32_t> global_fv(static_cast<std::size_t>(fsize));
            std::vector<std::int32_t> sorted_fv(static_cast<std::size_t>(fsize));
            for (int j = 0; j < fsize; ++j)
            {
                global_fv[static_cast<std::size_t>(j)] =
                    cell_verts[static_cast<std::size_t>(local_fv[static_cast<std::size_t>(j)])];
                sorted_fv[static_cast<std::size_t>(j)] = global_fv[static_cast<std::size_t>(j)];
            }
            std::sort(sorted_fv.begin(), sorted_fv.end());

            if (face_map.find(sorted_fv) == face_map.end())
            {
                face_map[sorted_fv] = ac.n_entities(2);
                ac.entity_types[2].push_back(cell::face_type(ctype, fi));
                for (int j = 0; j < fsize; ++j)
                    ac.entity_to_vertex[2].indices.push_back(global_fv[static_cast<std::size_t>(j)]);
                ac.entity_to_vertex[2].offsets.push_back(
                    static_cast<std::int32_t>(ac.entity_to_vertex[2].indices.size()));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// build_edges
// ---------------------------------------------------------------------------

template <std::floating_point T>
void build_edges(AdaptCell<T>& ac)
{
    const int tdim = ac.tdim;

    // Clear existing 1D entity pool.
    ac.entity_types[1].clear();
    ac.entity_to_vertex[1].offsets.clear();
    ac.entity_to_vertex[1].indices.clear();
    ac.entity_to_vertex[1].offsets.push_back(std::int32_t(0));

    // Track unique edges as (min_v, max_v) → edge_id.
    std::map<std::pair<std::int32_t, std::int32_t>, int> edge_map;

    const int n_cells = ac.n_entities(tdim);
    for (int c = 0; c < n_cells; ++c)
    {
        const cell::type ctype = ac.entity_types[tdim][static_cast<std::size_t>(c)];
        auto cell_verts = ac.entity_to_vertex[tdim][static_cast<std::int32_t>(c)];
        auto cell_edges = cell::edges(ctype);

        for (const auto& ce : cell_edges)
        {
            const std::int32_t lv0 = cell_verts[static_cast<std::size_t>(ce[0])];
            const std::int32_t lv1 = cell_verts[static_cast<std::size_t>(ce[1])];
            const auto key = std::make_pair(std::min(lv0, lv1), std::max(lv0, lv1));

            if (edge_map.find(key) == edge_map.end())
            {
                edge_map[key] = ac.n_entities(1);
                ac.entity_types[1].push_back(cell::type::interval);
                ac.entity_to_vertex[1].indices.push_back(lv0);
                ac.entity_to_vertex[1].indices.push_back(lv1);
                ac.entity_to_vertex[1].offsets.push_back(
                    static_cast<std::int32_t>(ac.entity_to_vertex[1].indices.size()));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// make_adapt_cell
// ---------------------------------------------------------------------------

template <std::floating_point T, std::integral I>
AdaptCell<T> make_adapt_cell(const MeshView<T, I>& mesh, I cell_id)
{
    if (!mesh.has_cell_types())
        throw std::runtime_error(
            "make_adapt_cell: MeshView must have cell types");

    const cell::type ctype = mesh.cell_type(cell_id);

    const int tdim = cell::get_tdim(ctype);
    const int gdim = mesh.gdim;
    const int nv   = cell::get_num_vertices(ctype);

    AdaptCell<T> ac;
    ac.gdim             = gdim;
    ac.tdim             = tdim;
    ac.parent_cell_type = ctype;
    ac.parent_cell_id   = static_cast<int>(cell_id);

    // Reference vertex coordinates: canonical corners of the reference element.
    // Stored as nv * tdim values.
    ac.vertex_coords = reference_vertices<T>(ctype);

    // Vertex provenance: all vertices coincide with parent-mesh vertices (dim 0).
    ac.vertex_parent_dim.assign(static_cast<std::size_t>(nv), std::int8_t(0));
    ac.vertex_parent_id.resize(static_cast<std::size_t>(nv));
    for (int v = 0; v < nv; ++v)
        ac.vertex_parent_id[static_cast<std::size_t>(v)] =
            static_cast<std::int32_t>(mesh.cell_node(cell_id, static_cast<I>(v)));

    // No parametric coordinates for parent vertices (empty param per vertex).
    ac.vertex_parent_param_offset.assign(static_cast<std::size_t>(nv + 1),
                                          std::int32_t(0));
    // vertex_parent_param stays empty.
    ac.vertex_source_edge_id.assign(static_cast<std::size_t>(nv), std::int32_t(-1));

    // Level-set sign masks start unset (bit not set = value is positive).
    ac.zero_mask_per_vertex.assign(static_cast<std::size_t>(nv),
                                    std::uint64_t(0));
    ac.negative_mask_per_vertex.assign(static_cast<std::size_t>(nv),
                                        std::uint64_t(0));

    // Entity pool at topological dimension: one entity = the full background cell.
    ac.entity_types[tdim].push_back(ctype);
    ac.entity_to_vertex[tdim].offsets = {std::int32_t(0), std::int32_t(nv)};
    ac.entity_to_vertex[tdim].indices.resize(static_cast<std::size_t>(nv));
    for (int v = 0; v < nv; ++v)
        ac.entity_to_vertex[tdim].indices[static_cast<std::size_t>(v)] =
            std::int32_t(v);
            
    // Create edges and faces (if tdim=3) in the entity pools 
    build_edges(ac);

    if (tdim == 3)
        build_faces(ac);

    recompute_active_level_set_masks(ac, /*num_level_sets=*/0);

    return ac;
}

// ---------------------------------------------------------------------------
// fill_vertex_signs
// ---------------------------------------------------------------------------

template <std::floating_point T>
void fill_vertex_signs(AdaptCell<T>& ac,
                       std::span<const T> vertex_ls_values,
                       int ls_index,
                       T tol){
    const int nv = ac.n_vertices();
    assert(static_cast<int>(vertex_ls_values.size()) >= nv
           && "fill_vertex_signs: vertex_ls_values too short");
    assert(ls_index >= 0 && ls_index < 64
           && "fill_vertex_signs: ls_index out of [0,63]");

    const std::uint64_t bit = std::uint64_t(1) << ls_index;

    for (int v = 0; v < nv; ++v)
    {
        const T val = vertex_ls_values[static_cast<std::size_t>(v)];
        if (std::fabs(val) < tol)
            ac.zero_mask_per_vertex[static_cast<std::size_t>(v)]     |= bit;
        else if (val < T(0))
            ac.negative_mask_per_vertex[static_cast<std::size_t>(v)] |= bit;
    }
}

// ---------------------------------------------------------------------------
// recompute_active_level_set_masks
// ---------------------------------------------------------------------------

template <std::floating_point T>
void recompute_active_level_set_masks(AdaptCell<T>& ac, int num_level_sets)
{
    const int tdim = ac.tdim;
    const int n_cells = ac.n_entities(tdim);

    ac.cell_active_level_set_mask.assign(static_cast<std::size_t>(n_cells),
                                         std::uint64_t(0));
    ac.active_level_set_mask = 0;

    if (num_level_sets <= 0 || n_cells <= 0)
        return;

    const int nls = std::min(num_level_sets, 64);

    for (int c = 0; c < n_cells; ++c)
    {
        auto verts = ac.entity_to_vertex[tdim][static_cast<std::int32_t>(c)];
        std::uint64_t leaf_mask = 0;

        for (int ls = 0; ls < nls; ++ls)
        {
            const std::uint64_t bit = std::uint64_t(1) << ls;
            bool has_neg = false;
            bool has_pos = false;
            bool has_zero = false;

            for (const auto v : verts)
            {
                const auto zv = ac.zero_mask_per_vertex[static_cast<std::size_t>(v)];
                const auto nv = ac.negative_mask_per_vertex[static_cast<std::size_t>(v)];
                if ((zv & bit) != 0)
                    has_zero = true;
                else if ((nv & bit) != 0)
                    has_neg = true;
                else
                    has_pos = true;
            }

            if (has_zero || (has_neg && has_pos))
                leaf_mask |= bit;
        }

        ac.cell_active_level_set_mask[static_cast<std::size_t>(c)] = leaf_mask;
        ac.active_level_set_mask |= leaf_mask;
    }
}

// ---------------------------------------------------------------------------
// rebuild_zero_entity_inventory
// ---------------------------------------------------------------------------

template <std::floating_point T>
void rebuild_zero_entity_inventory(AdaptCell<T>& ac)
{
    ac.zero_entity_dim.clear();
    ac.zero_entity_id.clear();
    ac.zero_entity_zero_mask.clear();
    ac.zero_entity_is_owned.clear();
    ac.zero_entity_parent_dim.clear();
    ac.zero_entity_parent_id.clear();

    const int n_vertices = ac.n_vertices();
    for (int v = 0; v < n_vertices; ++v)
    {
        const auto mask = ac.zero_mask_per_vertex[static_cast<std::size_t>(v)];
        if (mask == 0)
            continue;

        ac.zero_entity_dim.push_back(std::uint8_t(0));
        ac.zero_entity_id.push_back(static_cast<std::int32_t>(v));
        ac.zero_entity_zero_mask.push_back(mask);
        ac.zero_entity_is_owned.push_back(std::uint8_t(1));
        ac.zero_entity_parent_dim.push_back(
            ac.vertex_parent_dim.empty() ? std::int8_t(-1)
                                         : ac.vertex_parent_dim[static_cast<std::size_t>(v)]);
        ac.zero_entity_parent_id.push_back(
            ac.vertex_parent_id.empty() ? std::int32_t(-1)
                                        : ac.vertex_parent_id[static_cast<std::size_t>(v)]);
    }

    for (int dim = 1; dim < ac.tdim; ++dim)
    {
        const int n_entities = ac.n_entities(dim);
        for (int e = 0; e < n_entities; ++e)
        {
            auto verts = ac.entity_to_vertex[dim][static_cast<std::int32_t>(e)];
            if (verts.empty())
                continue;

            std::uint64_t mask = ~std::uint64_t(0);
            for (const auto v : verts)
            {
                mask &= ac.zero_mask_per_vertex[static_cast<std::size_t>(v)];
                if (mask == 0)
                    break;
            }
            if (mask == 0)
                continue;

            ac.zero_entity_dim.push_back(static_cast<std::uint8_t>(dim));
            ac.zero_entity_id.push_back(static_cast<std::int32_t>(e));
            ac.zero_entity_zero_mask.push_back(mask);
            ac.zero_entity_is_owned.push_back(std::uint8_t(1));
            ac.zero_entity_parent_dim.push_back(std::int8_t(-1));
            ac.zero_entity_parent_id.push_back(std::int32_t(-1));
        }
    }

    ++ac.zero_entity_version;
}

// ---------------------------------------------------------------------------
// Explicit template instantiations
// ---------------------------------------------------------------------------

template void build_edges(AdaptCell<double>&);
template void build_edges(AdaptCell<float>&);

template void build_faces(AdaptCell<double>&);
template void build_faces(AdaptCell<float>&);

template AdaptCell<double> make_adapt_cell(const MeshView<double, int>&,  int);
template AdaptCell<float>  make_adapt_cell(const MeshView<float,  int>&,  int);
template AdaptCell<double> make_adapt_cell(const MeshView<double, long>&, long);
template AdaptCell<float>  make_adapt_cell(const MeshView<float,  long>&, long);

template void fill_vertex_signs(AdaptCell<double>&, std::span<const double>, int, double);
template void fill_vertex_signs(AdaptCell<float>&,  std::span<const float>,  int, float);
template void recompute_active_level_set_masks(AdaptCell<double>&, int);
template void recompute_active_level_set_masks(AdaptCell<float>&, int);
template void rebuild_zero_entity_inventory(AdaptCell<double>&);
template void rebuild_zero_entity_inventory(AdaptCell<float>&);

} // namespace cutcells
