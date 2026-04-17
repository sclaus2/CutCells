// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "adapt_cell.h"
#include "cell_topology.h"
#include "reference_cell.h"

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

} // namespace cutcells
