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
#include <vector>

namespace cutcells
{
namespace
{

template <std::floating_point T>
void clear_entity_host_provenance(AdaptCell<T>& ac, int dim)
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
std::span<const std::int32_t> host_vertices_for_entity(const AdaptCell<T>& ac,
                                                       int dim,
                                                       int entity_id)
{
    if (dim >= 0 && dim < AdaptCell<T>::max_dim
        && entity_id >= 0
        && entity_id < ac.entity_host_cell_vertices[dim].size())
    {
        return ac.entity_host_cell_vertices[dim][static_cast<std::int32_t>(entity_id)];
    }
    return {};
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
T squared_distance_to_segment(std::span<const T> p,
                              std::span<const T> a,
                              std::span<const T> b,
                              int dim)
{
    T ab2 = T(0);
    T ap_ab = T(0);
    for (int d = 0; d < dim; ++d)
    {
        const T ab = b[static_cast<std::size_t>(d)] - a[static_cast<std::size_t>(d)];
        const T ap = p[static_cast<std::size_t>(d)] - a[static_cast<std::size_t>(d)];
        ab2 += ab * ab;
        ap_ab += ap * ab;
    }
    const T t = (ab2 > T(0)) ? std::clamp(ap_ab / ab2, T(0), T(1)) : T(0);
    T dist2 = T(0);
    for (int d = 0; d < dim; ++d)
    {
        const T x = a[static_cast<std::size_t>(d)]
                  + t * (b[static_cast<std::size_t>(d)] - a[static_cast<std::size_t>(d)]);
        const T r = p[static_cast<std::size_t>(d)] - x;
        dist2 += r * r;
    }
    return dist2;
}

template <std::floating_point T>
bool point_in_face_span(std::span<const T> p,
                        const std::vector<T>& ref_vertices,
                        int tdim,
                        std::span<const int> face_vertices,
                        T tol)
{
    if (tdim != 3 || face_vertices.size() < 3)
        return false;

    const T* a = ref_vertices.data() + static_cast<std::size_t>(face_vertices[0] * tdim);
    const T* b = ref_vertices.data() + static_cast<std::size_t>(face_vertices[1] * tdim);
    const T* c = ref_vertices.data() + static_cast<std::size_t>(face_vertices[2] * tdim);
    std::array<T, 3> u = {};
    std::array<T, 3> v = {};
    std::array<T, 3> w = {};
    for (int d = 0; d < 3; ++d)
    {
        u[static_cast<std::size_t>(d)] = b[d] - a[d];
        v[static_cast<std::size_t>(d)] = c[d] - a[d];
        w[static_cast<std::size_t>(d)] = p[static_cast<std::size_t>(d)] - a[d];
    }
    const std::array<T, 3> n = {
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0]};
    const T n2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
    if (n2 <= tol * tol)
        return false;
    const T signed_dist =
        (n[0] * w[0] + n[1] * w[1] + n[2] * w[2]) / std::sqrt(n2);
    return std::fabs(signed_dist) <= tol;
}

template <std::floating_point T>
std::pair<std::int8_t, std::int32_t>
infer_zero_entity_parent_host(const AdaptCell<T>& ac,
                              std::span<const std::int32_t> entity_vertices)
{
    const T tol = T(256) * std::numeric_limits<T>::epsilon();
    const int tdim = ac.tdim;
    const auto ref_vertices = cell::reference_vertices<T>(ac.parent_cell_type);

    if (tdim >= 1)
    {
        auto parent_edges = cell::edges(ac.parent_cell_type);
        for (int e = 0; e < static_cast<int>(parent_edges.size()); ++e)
        {
            const auto edge = parent_edges[static_cast<std::size_t>(e)];
            std::span<const T> a(
                ref_vertices.data() + static_cast<std::size_t>(edge[0] * tdim),
                static_cast<std::size_t>(tdim));
            std::span<const T> b(
                ref_vertices.data() + static_cast<std::size_t>(edge[1] * tdim),
                static_cast<std::size_t>(tdim));

            bool all_on_edge = true;
            for (const auto v : entity_vertices)
            {
                std::span<const T> p(
                    ac.vertex_coords.data() + static_cast<std::size_t>(v * tdim),
                    static_cast<std::size_t>(tdim));
                if (squared_distance_to_segment<T>(p, a, b, tdim) > tol * tol)
                {
                    all_on_edge = false;
                    break;
                }
            }
            if (all_on_edge)
                return {std::int8_t(1), std::int32_t(e)};
        }
    }

    if (tdim == 3)
    {
        const int nfaces = cell::num_faces(ac.parent_cell_type);
        for (int f = 0; f < nfaces; ++f)
        {
            auto face_vertices = cell::face_vertices(ac.parent_cell_type, f);
            bool all_on_face = true;
            for (const auto v : entity_vertices)
            {
                std::span<const T> p(
                    ac.vertex_coords.data() + static_cast<std::size_t>(v * tdim),
                    static_cast<std::size_t>(tdim));
                if (!point_in_face_span<T>(p, ref_vertices, tdim, face_vertices, tol))
                {
                    all_on_face = false;
                    break;
                }
            }
            if (all_on_face)
                return {std::int8_t(2), std::int32_t(f)};
        }
    }

    return {static_cast<std::int8_t>(tdim), static_cast<std::int32_t>(ac.parent_cell_id)};
}

} // namespace

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
    clear_entity_host_provenance<T>(ac, 2);

    // Invalidate face cert tag storage: face count is about to change and
    // faces may be reordered, so all old tags are stale.
    ac.face_cert_tag.clear();
    ac.face_cert_tag_num_level_sets = 0;

    // Track unique faces by sorted vertex set → face_id.
    std::map<std::vector<std::int32_t>, int> face_map;

    const int n_cells = ac.n_entities(ac.tdim);
    std::vector<std::vector<std::int32_t>> face_to_cells;
    std::vector<std::vector<std::int32_t>> cell_to_faces(static_cast<std::size_t>(n_cells));
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

            auto face_it = face_map.find(sorted_fv);
            int face_id = -1;
            if (face_it == face_map.end())
            {
                face_id = ac.n_entities(2);
                face_map[sorted_fv] = face_id;
                ac.entity_types[2].push_back(cell::face_type(ctype, fi));
                for (int j = 0; j < fsize; ++j)
                    ac.entity_to_vertex[2].indices.push_back(global_fv[static_cast<std::size_t>(j)]);
                ac.entity_to_vertex[2].offsets.push_back(
                    static_cast<std::int32_t>(ac.entity_to_vertex[2].indices.size()));
                face_to_cells.emplace_back();

                const int host_cell_id =
                    (ac.tdim < AdaptCell<T>::max_dim
                     && c < static_cast<int>(ac.entity_host_cell_id[ac.tdim].size()))
                        ? ac.entity_host_cell_id[ac.tdim][static_cast<std::size_t>(c)]
                    : (c < static_cast<int>(ac.cell_source_cell_id.size()))
                        ? ac.cell_source_cell_id[static_cast<std::size_t>(c)]
                        : c;
                const cell::type host_cell_type =
                    (ac.tdim < AdaptCell<T>::max_dim
                     && c < static_cast<int>(ac.entity_host_cell_type[ac.tdim].size()))
                        ? ac.entity_host_cell_type[ac.tdim][static_cast<std::size_t>(c)]
                        : ctype;
                const int source_level_set =
                    (ac.tdim < AdaptCell<T>::max_dim
                     && c < static_cast<int>(ac.entity_source_level_set[ac.tdim].size()))
                        ? ac.entity_source_level_set[ac.tdim][static_cast<std::size_t>(c)]
                        : -1;
                const auto host_vertices =
                    host_vertices_for_entity<T>(ac, ac.tdim, c);
                append_entity_host_provenance<T>(
                    ac, 2, host_cell_id, host_cell_type, -1, source_level_set,
                    host_vertices.empty() ? std::span<const std::int32_t>(cell_verts)
                                          : host_vertices);
            }
            else
                face_id = face_it->second;

            face_to_cells[static_cast<std::size_t>(face_id)].push_back(
                static_cast<std::int32_t>(c));
            cell_to_faces[static_cast<std::size_t>(c)].push_back(
                static_cast<std::int32_t>(face_id));
        }
    }

    auto fill_connectivity = [](EntityAdjacency& adjacency,
                                const std::vector<std::vector<std::int32_t>>& rows)
    {
        adjacency.indices.clear();
        adjacency.offsets.clear();
        adjacency.offsets.push_back(std::int32_t(0));
        for (const auto& row : rows)
        {
            for (const auto value : row)
                adjacency.indices.push_back(value);
            adjacency.offsets.push_back(
                static_cast<std::int32_t>(adjacency.indices.size()));
        }
    };
    fill_connectivity(ac.connectivity[2][ac.tdim], face_to_cells);
    ac.has_connectivity[2][ac.tdim] = 1;
    fill_connectivity(ac.connectivity[ac.tdim][2], cell_to_faces);
    ac.has_connectivity[ac.tdim][2] = 1;
}

// ---------------------------------------------------------------------------
// build_edges
// ---------------------------------------------------------------------------

template <std::floating_point T>
void build_edges(AdaptCell<T>& ac)
{
    const int tdim = ac.tdim;
    if (tdim == 1)
        return;

    // Clear existing 1D entity pool.
    ac.entity_types[1].clear();
    ac.entity_to_vertex[1].offsets.clear();
    ac.entity_to_vertex[1].indices.clear();
    ac.entity_to_vertex[1].offsets.push_back(std::int32_t(0));
    clear_entity_host_provenance<T>(ac, 1);

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

                const int host_cell_id =
                    (tdim < AdaptCell<T>::max_dim
                     && c < static_cast<int>(ac.entity_host_cell_id[tdim].size()))
                        ? ac.entity_host_cell_id[tdim][static_cast<std::size_t>(c)]
                    : (c < static_cast<int>(ac.cell_source_cell_id.size()))
                        ? ac.cell_source_cell_id[static_cast<std::size_t>(c)]
                        : c;
                const cell::type host_cell_type =
                    (tdim < AdaptCell<T>::max_dim
                     && c < static_cast<int>(ac.entity_host_cell_type[tdim].size()))
                        ? ac.entity_host_cell_type[tdim][static_cast<std::size_t>(c)]
                        : ctype;
                const int source_level_set =
                    (tdim < AdaptCell<T>::max_dim
                     && c < static_cast<int>(ac.entity_source_level_set[tdim].size()))
                        ? ac.entity_source_level_set[tdim][static_cast<std::size_t>(c)]
                        : -1;
                const auto host_vertices =
                    host_vertices_for_entity<T>(ac, tdim, c);
                append_entity_host_provenance<T>(
                    ac, 1, host_cell_id, host_cell_type, -1, source_level_set,
                    host_vertices.empty() ? std::span<const std::int32_t>(cell_verts)
                                          : host_vertices);
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
    std::vector<I> cell_node_scratch;
    const auto parent_nodes = mesh.cell_nodes(cell_id, cell_node_scratch);
    for (int v = 0; v < nv; ++v)
        ac.vertex_parent_id[static_cast<std::size_t>(v)] =
            static_cast<std::int32_t>(parent_nodes[static_cast<std::size_t>(v)]);

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

    ac.cell_source_cell_id.assign(1, std::int32_t(0));
    ac.cell_refinement_generation.assign(1, std::int32_t(0));
    ac.cell_refinement_reason.assign(1, CellRefinementReason::none);
    ac.cell_host_parent_cell_id.assign(1, ac.parent_cell_id);
    clear_entity_host_provenance<T>(ac, tdim);
    append_entity_host_provenance<T>(
        ac,
        tdim,
        0,
        ctype,
        -1,
        -1,
        std::span<const std::int32_t>(
            ac.entity_to_vertex[tdim].indices.data(),
            ac.entity_to_vertex[tdim].indices.size()));
            
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
    ac.zero_entity_host_cell_id.clear();
    ac.zero_entity_host_cell_type.clear();
    ac.zero_entity_host_face_id.clear();
    ac.zero_entity_source_level_set.clear();
    ac.zero_entity_host_cell_vertices.offsets.clear();
    ac.zero_entity_host_cell_vertices.indices.clear();
    ac.zero_entity_host_cell_vertices.offsets.push_back(std::int32_t(0));

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
        ac.zero_entity_host_cell_id.push_back(std::int32_t(-1));
        ac.zero_entity_host_cell_type.push_back(cell::type::point);
        ac.zero_entity_host_face_id.push_back(std::int32_t(-1));
        ac.zero_entity_source_level_set.push_back(std::int32_t(-1));
        ac.zero_entity_host_cell_vertices.offsets.push_back(
            static_cast<std::int32_t>(ac.zero_entity_host_cell_vertices.indices.size()));
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
            const auto host = infer_zero_entity_parent_host<T>(ac, verts);
            ac.zero_entity_parent_dim.push_back(host.first);
            ac.zero_entity_parent_id.push_back(host.second);
            if (e < static_cast<int>(ac.entity_host_cell_id[dim].size()))
            {
                ac.zero_entity_host_cell_id.push_back(
                    ac.entity_host_cell_id[dim][static_cast<std::size_t>(e)]);
                ac.zero_entity_host_cell_type.push_back(
                    ac.entity_host_cell_type[dim][static_cast<std::size_t>(e)]);
                ac.zero_entity_host_face_id.push_back(
                    ac.entity_host_face_id[dim][static_cast<std::size_t>(e)]);
                ac.zero_entity_source_level_set.push_back(
                    ac.entity_source_level_set[dim][static_cast<std::size_t>(e)]);
                const auto host_vertices =
                    ac.entity_host_cell_vertices[dim][static_cast<std::int32_t>(e)];
                for (const auto hv : host_vertices)
                    ac.zero_entity_host_cell_vertices.indices.push_back(hv);
            }
            else
            {
                ac.zero_entity_host_cell_id.push_back(std::int32_t(-1));
                ac.zero_entity_host_cell_type.push_back(cell::type::point);
                ac.zero_entity_host_face_id.push_back(std::int32_t(-1));
                ac.zero_entity_source_level_set.push_back(std::int32_t(-1));
            }
            ac.zero_entity_host_cell_vertices.offsets.push_back(
                static_cast<std::int32_t>(ac.zero_entity_host_cell_vertices.indices.size()));
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
