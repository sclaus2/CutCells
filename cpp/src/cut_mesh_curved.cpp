// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "cut_mesh_curved.h"

#include "bernstein_backend.h"
#include "cell_flags.h"
#include "local_mesh.h"
#include "mapping.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

template <std::floating_point T>
void init_curved_mesh(cutcells::mesh::CurvedGlobalMesh<T>& mesh,
                      int gdim,
                      int geom_order,
                      int level_set_id)
{
    mesh.vertex_coords.clear();
    mesh.connectivity.clear();
    mesh.offsets.clear();
    mesh.cell_types.clear();
    mesh.gdim = gdim;
    mesh.geom_order = geom_order;
    mesh.level_set_id = level_set_id;
    mesh.n_fallback_nodes = 0;
    mesh.offsets.push_back(0);
}

inline std::string to_lower_ascii(std::string_view text)
{
    std::string out(text);
    std::transform(out.begin(), out.end(), out.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return out;
}

inline bool bernstein_supported_cell(cutcells::cell::type ct)
{
    using cutcells::cell::type;
    return ct == type::interval || ct == type::triangle
           || ct == type::tetrahedron || ct == type::hexahedron;
}

template <std::floating_point T>
void ref_to_phys_affine_embedded(cutcells::cell::type ct,
                                 std::span<const T> parent_coords,
                                 int gdim,
                                 std::span<const T> x_ref,
                                 T* x_phys)
{
    if (gdim <= 0 || gdim > 3)
        throw std::invalid_argument("ref_to_phys_affine_embedded: gdim must be in [1, 3]");

    const int tdim = cutcells::cell::get_tdim(ct);
    if (static_cast<int>(x_ref.size()) != tdim)
        throw std::invalid_argument("ref_to_phys_affine_embedded: reference point has invalid dimension");

    const int p1_vertices = cutcells::cell::get_num_vertices(ct);
    if (static_cast<int>(parent_coords.size()) != p1_vertices * gdim)
        throw std::invalid_argument("ref_to_phys_affine_embedded: invalid parent_coords size");

    const T* x0 = parent_coords.data();
    for (int d = 0; d < gdim; ++d)
        x_phys[d] = x0[d];

    if (tdim == 0)
        return;

    const auto cols = cutcells::cell::jacobian_col_indices(ct);
    for (int k = 0; k < tdim; ++k)
    {
        const int vk = cols[static_cast<std::size_t>(k)];
        if (vk < 0 || vk >= p1_vertices)
            throw std::runtime_error("ref_to_phys_affine_embedded: invalid Jacobian column index");

        const T* xk = parent_coords.data() + static_cast<std::size_t>(vk * gdim);
        const T xi = x_ref[static_cast<std::size_t>(k)];
        for (int d = 0; d < gdim; ++d)
            x_phys[d] += xi * (xk[d] - x0[d]);
    }
}

enum class ParentCellDomain : uint8_t
{
    inside = 0,
    intersected = 1,
    outside = 2
};

template <std::floating_point T>
ParentCellDomain classify_parent_cell_bernstein(cutcells::cell::type ct,
                                                int degree,
                                                std::span<const T> nodal_values,
                                                T tol)
{
    const auto poly = cutcells::make_bernstein_cell<T>(ct, degree, nodal_values);

    T scale = T(0);
    for (const T c : poly.coeffs)
        scale = std::max(scale, std::abs(c));
    const T tol_eff = std::max(
        tol,
        std::numeric_limits<T>::epsilon() * std::max(T(1), scale) * T(32));

    bool all_neg = true;
    bool all_pos = true;
    for (const T c : poly.coeffs)
    {
        all_neg = all_neg && (c < -tol_eff);
        all_pos = all_pos && (c > tol_eff);
    }

    if (all_neg)
        return ParentCellDomain::inside;
    if (all_pos)
        return ParentCellDomain::outside;
    return ParentCellDomain::intersected;
}

template <std::floating_point T>
void sample_callable_on_parent_nodes(
    cutcells::cell::type                 ct,
    std::span<const T>                   parent_coords_p1,
    int                                  gdim,
    int                                  degree,
    const cutcells::LevelSetFunction<T>& level_set,
    int                                  cell_id,
    std::vector<T>&                      sampled_node_coords,
    std::vector<T>&                      sampled_nodal_values)
{
    const int tdim = cutcells::cell::get_tdim(ct);
    const int n_nodes = cutcells::lagrange_node_count(ct, degree);
    const auto ref_nodes = cutcells::bernstein_reference_nodes(ct, degree);
    if (static_cast<int>(ref_nodes.size()) != n_nodes * tdim)
        throw std::runtime_error("sample_callable_on_parent_nodes: unexpected reference-node count");

    sampled_node_coords.resize(static_cast<std::size_t>(n_nodes * gdim));
    sampled_nodal_values.resize(static_cast<std::size_t>(n_nodes));

    std::array<T, 3> x_phys = {T(0), T(0), T(0)};
    std::array<T, 3> x_query = {T(0), T(0), T(0)};

    for (int i = 0; i < n_nodes; ++i)
    {
        std::array<T, 3> x_ref = {T(0), T(0), T(0)};
        for (int d = 0; d < tdim; ++d)
        {
            x_ref[static_cast<std::size_t>(d)] = static_cast<T>(
                ref_nodes[static_cast<std::size_t>(i * tdim + d)]);
        }

        ref_to_phys_affine_embedded<T>(
            ct,
            parent_coords_p1,
            gdim,
            std::span<const T>(x_ref.data(), static_cast<std::size_t>(tdim)),
            x_phys.data());

        for (int d = 0; d < gdim; ++d)
        {
            const T value = x_phys[static_cast<std::size_t>(d)];
            sampled_node_coords[static_cast<std::size_t>(i * gdim + d)] = value;
            x_query[static_cast<std::size_t>(d)] = value;
        }
        for (int d = gdim; d < 3; ++d)
            x_query[static_cast<std::size_t>(d)] = T(0);

        sampled_nodal_values[static_cast<std::size_t>(i)] = level_set.value(
            x_query.data(), static_cast<int>(cell_id));
    }
}

template <std::floating_point T>
void append_cell_to_curved_mesh(cutcells::mesh::CurvedGlobalMesh<T>& mesh,
                                cutcells::cell::type                 ct,
                                std::span<const T>                   coords,
                                int                                  gdim)
{
    if (mesh.offsets.empty())
        mesh.offsets.push_back(0);

    const int n_vertices = static_cast<int>(coords.size()) / gdim;
    const int32_t v0 = static_cast<int32_t>(mesh.vertex_coords.size() / static_cast<std::size_t>(gdim));

    mesh.vertex_coords.insert(mesh.vertex_coords.end(), coords.begin(), coords.end());
    for (int i = 0; i < n_vertices; ++i)
        mesh.connectivity.push_back(v0 + static_cast<int32_t>(i));

    mesh.offsets.push_back(static_cast<int32_t>(mesh.connectivity.size()));
    mesh.cell_types.push_back(ct);
}

template <std::floating_point T>
void append_local_cells_by_domain(const cutcells::LocalMesh<T>&      local_mesh,
                                  cutcells::cell::domain             domain,
                                  cutcells::mesh::CurvedGlobalMesh<T>& out)
{
    const int gdim = local_mesh.gdim;
    const int nc = local_mesh.n_cells();

    for (int c = 0; c < nc; ++c)
    {
        if (local_mesh.cell_domain[static_cast<std::size_t>(c)] != static_cast<uint8_t>(domain))
            continue;

        const int c0 = local_mesh.cell_offsets[static_cast<std::size_t>(c)];
        const int c1 = local_mesh.cell_offsets[static_cast<std::size_t>(c + 1)];
        const int n_vertices = c1 - c0;
        const int32_t v0 = static_cast<int32_t>(out.vertex_coords.size() / static_cast<std::size_t>(gdim));

        for (int i = c0; i < c1; ++i)
        {
            const int lv = local_mesh.cell_vertices[static_cast<std::size_t>(i)];
            if (lv < 0 || lv >= local_mesh.n_vertices())
                throw std::runtime_error("append_local_cells_by_domain: local cell references invalid vertex id");
            const std::size_t base = static_cast<std::size_t>(lv * gdim);
            out.vertex_coords.insert(
                out.vertex_coords.end(),
                local_mesh.vertex_x.begin() + static_cast<std::ptrdiff_t>(base),
                local_mesh.vertex_x.begin() + static_cast<std::ptrdiff_t>(base + gdim));
        }

        for (int i = 0; i < n_vertices; ++i)
            out.connectivity.push_back(v0 + static_cast<int32_t>(i));
        out.offsets.push_back(static_cast<int32_t>(out.connectivity.size()));
        out.cell_types.push_back(local_mesh.cell_types[static_cast<std::size_t>(c)]);
    }
}

template <std::floating_point T>
cutcells::mesh::CutMesh<T> curved_to_cut_mesh_direct(
    const cutcells::mesh::CurvedGlobalMesh<T>& curved)
{
    cutcells::mesh::CutMesh<T> out;
    out._gdim = curved.gdim;
    out._tdim = curved.cell_types.empty() ? 0 : cutcells::cell::get_tdim(curved.cell_types.front());
    out._num_cells = curved.n_cells();
    out._num_vertices = curved.n_vertices();
    out._vertex_coords = curved.vertex_coords;

    out._connectivity.resize(curved.connectivity.size());
    std::transform(curved.connectivity.begin(), curved.connectivity.end(),
                   out._connectivity.begin(),
                   [](int32_t v) { return static_cast<int>(v); });

    out._offset.resize(curved.offsets.size());
    std::transform(curved.offsets.begin(), curved.offsets.end(),
                   out._offset.begin(),
                   [](int32_t v) { return static_cast<int>(v); });

    out._types = curved.cell_types;
    out._parent_map.assign(static_cast<std::size_t>(out._num_cells), -1);
    return out;
}

template <std::floating_point T>
std::vector<T> gauss_lobatto_nodes_01(int order)
{
    std::vector<T> nodes;
    if (order <= 1)
    {
        nodes = {T(0), T(1)};
        return nodes;
    }

    nodes.reserve(static_cast<std::size_t>(order + 1));
    nodes.push_back(T(0));
    const auto interior = cutcells::detail::gauss_lobatto_interior_points_1d<T>(order);
    nodes.insert(nodes.end(), interior.begin(), interior.end());
    nodes.push_back(T(1));
    return nodes;
}

template <std::floating_point T>
void lagrange_eval_nd(std::span<const T> t_nodes,
                      std::span<const T> values,
                      int                gdim,
                      T                  t,
                      T*                 out)
{
    const int n = static_cast<int>(t_nodes.size());
    for (int d = 0; d < gdim; ++d)
        out[d] = T(0);

    for (int i = 0; i < n; ++i)
    {
        T li = T(1);
        const T ti = t_nodes[static_cast<std::size_t>(i)];
        for (int j = 0; j < n; ++j)
        {
            if (i == j)
                continue;
            const T tj = t_nodes[static_cast<std::size_t>(j)];
            li *= (t - tj) / (ti - tj);
        }

        for (int d = 0; d < gdim; ++d)
            out[d] += li * values[static_cast<std::size_t>(i * gdim + d)];
    }
}

template <std::floating_point T>
void append_line_segments_from_polyline(std::span<const T> polyline,
                                        int                gdim,
                                        cutcells::mesh::CutMesh<T>& out)
{
    const int n = static_cast<int>(polyline.size()) / gdim;
    if (n < 2)
        return;

    if (out._offset.empty())
        out._offset.push_back(0);

    for (int i = 0; i < n - 1; ++i)
    {
        const int v0 = static_cast<int>(out._vertex_coords.size() / static_cast<std::size_t>(gdim));
        const std::size_t base0 = static_cast<std::size_t>(i * gdim);
        const std::size_t base1 = static_cast<std::size_t>((i + 1) * gdim);

        out._vertex_coords.insert(
            out._vertex_coords.end(),
            polyline.begin() + static_cast<std::ptrdiff_t>(base0),
            polyline.begin() + static_cast<std::ptrdiff_t>(base0 + gdim));
        out._vertex_coords.insert(
            out._vertex_coords.end(),
            polyline.begin() + static_cast<std::ptrdiff_t>(base1),
            polyline.begin() + static_cast<std::ptrdiff_t>(base1 + gdim));

        out._connectivity.push_back(v0);
        out._connectivity.push_back(v0 + 1);
        out._offset.push_back(static_cast<int>(out._connectivity.size()));
        out._types.push_back(cutcells::cell::type::interval);
        out._parent_map.push_back(-1);
    }
}

template <std::floating_point T>
cutcells::mesh::CutMesh<T> linearize_curved_interface_mesh(
    const cutcells::mesh::CurvedGlobalMesh<T>& curved,
    int                                        subdivision)
{
    cutcells::mesh::CutMesh<T> out;
    out._gdim = curved.gdim;
    out._tdim = curved.cell_types.empty() ? 0 : cutcells::cell::get_tdim(curved.cell_types.front());
    out._offset.push_back(0);

    const int gdim = curved.gdim;

    for (int c = 0; c < curved.n_cells(); ++c)
    {
        const int32_t c0 = curved.offsets[static_cast<std::size_t>(c)];
        const int32_t c1 = curved.offsets[static_cast<std::size_t>(c + 1)];
        const int n_nodes = static_cast<int>(c1 - c0);
        const auto ct = curved.cell_types[static_cast<std::size_t>(c)];

        if (ct == cutcells::cell::type::interval)
        {
            if (n_nodes < 2)
                continue;

            std::vector<T> cell_points(static_cast<std::size_t>(n_nodes * gdim), T(0));
            for (int i = 0; i < n_nodes; ++i)
            {
                const int32_t gv = curved.connectivity[static_cast<std::size_t>(c0 + i)];
                for (int d = 0; d < gdim; ++d)
                {
                    cell_points[static_cast<std::size_t>(i * gdim + d)] =
                        curved.vertex_coords[static_cast<std::size_t>(gv * gdim + d)];
                }
            }

            if (n_nodes == 2)
            {
                append_line_segments_from_polyline<T>(
                    std::span<const T>(cell_points.data(), cell_points.size()), gdim, out);
                continue;
            }

            const int order = n_nodes - 1;
            const int n_segments = std::max(1, subdivision * order);
            std::vector<T> sampled(static_cast<std::size_t>((n_segments + 1) * gdim), T(0));
            auto t_nodes = gauss_lobatto_nodes_01<T>(order);
            if (static_cast<int>(t_nodes.size()) != n_nodes)
            {
                t_nodes.resize(static_cast<std::size_t>(n_nodes));
                for (int i = 0; i < n_nodes; ++i)
                    t_nodes[static_cast<std::size_t>(i)] = static_cast<T>(i) / static_cast<T>(n_nodes - 1);
            }

            for (int i = 0; i <= n_segments; ++i)
            {
                const T t = static_cast<T>(i) / static_cast<T>(n_segments);
                lagrange_eval_nd<T>(
                    std::span<const T>(t_nodes.data(), t_nodes.size()),
                    std::span<const T>(cell_points.data(), cell_points.size()),
                    gdim,
                    t,
                    sampled.data() + static_cast<std::size_t>(i * gdim));
            }

            append_line_segments_from_polyline<T>(
                std::span<const T>(sampled.data(), sampled.size()), gdim, out);
            continue;
        }

        const int v0 = static_cast<int>(out._vertex_coords.size() / static_cast<std::size_t>(gdim));
        for (int i = 0; i < n_nodes; ++i)
        {
            const int32_t gv = curved.connectivity[static_cast<std::size_t>(c0 + i)];
            const std::size_t base = static_cast<std::size_t>(gv * gdim);
            out._vertex_coords.insert(
                out._vertex_coords.end(),
                curved.vertex_coords.begin() + static_cast<std::ptrdiff_t>(base),
                curved.vertex_coords.begin() + static_cast<std::ptrdiff_t>(base + gdim));
            out._connectivity.push_back(v0 + i);
        }
        out._offset.push_back(static_cast<int>(out._connectivity.size()));
        out._types.push_back(ct);
        out._parent_map.push_back(-1);
    }

    out._num_cells = static_cast<int>(out._types.size());
    out._num_vertices = static_cast<int>(out._vertex_coords.size() / static_cast<std::size_t>(out._gdim));
    if (out._types.empty())
        out._tdim = 0;
    return out;
}

template <std::floating_point T>
cutcells::mesh::CutMesh<T> linearize_local_interfaces(
    std::span<const cutcells::LocalMesh<T>*> local_meshes,
    int                                      level_set_id,
    int                                      subdivision)
{
    cutcells::mesh::CutMesh<T> out;
    out._offset.push_back(0);
    (void)subdivision;

    if (local_meshes.empty())
    {
        out._gdim = 0;
        out._tdim = 0;
        out._num_cells = 0;
        out._num_vertices = 0;
        return out;
    }

    out._gdim = local_meshes.front()->gdim;
    out._tdim = 1;

    const uint64_t zero_mask = uint64_t(1) << level_set_id;

    for (const cutcells::LocalMesh<T>* mesh_ptr : local_meshes)
    {
        if (mesh_ptr == nullptr)
            continue;

        const auto& mesh = *mesh_ptr;
        if (mesh.tdim != 2)
            continue;

        if (mesh.n_zero_entities() == 0)
            continue;

        if (mesh.n_zero_chains() == 0)
        {
            auto& nonconst_mesh = const_cast<cutcells::LocalMesh<T>&>(mesh);
            build_zero_chains(nonconst_mesh, zero_mask);
        }

        for (int chain = 0; chain < mesh.n_zero_chains(); ++chain)
        {
            if (mesh.zero_chain_zero_mask[static_cast<std::size_t>(chain)] != zero_mask)
                continue;

            const int c0 = mesh.zero_chain_offsets[static_cast<std::size_t>(chain)];
            const int c1 = mesh.zero_chain_offsets[static_cast<std::size_t>(chain + 1)];

            for (int ci = c0; ci < c1; ++ci)
            {
                const int ze = mesh.zero_chain_entity_ids[static_cast<std::size_t>(ci)];
                const bool reversed = mesh.zero_chain_entity_reversed[static_cast<std::size_t>(ci)] != 0;

                const int32_t lv0 = reversed
                    ? mesh.zero_entity_endpoint_v1[static_cast<std::size_t>(ze)]
                    : mesh.zero_entity_endpoint_v0[static_cast<std::size_t>(ze)];
                const int32_t lv1 = reversed
                    ? mesh.zero_entity_endpoint_v0[static_cast<std::size_t>(ze)]
                    : mesh.zero_entity_endpoint_v1[static_cast<std::size_t>(ze)];
                if (lv0 < 0 || lv1 < 0)
                    continue;

                std::vector<T> polyline(static_cast<std::size_t>(2 * mesh.gdim), T(0));
                for (int d = 0; d < mesh.gdim; ++d)
                {
                    polyline[static_cast<std::size_t>(d)] =
                        mesh.vertex_x[static_cast<std::size_t>(lv0 * mesh.gdim + d)];
                    polyline[static_cast<std::size_t>(mesh.gdim + d)] =
                        mesh.vertex_x[static_cast<std::size_t>(lv1 * mesh.gdim + d)];
                }

                append_line_segments_from_polyline<T>(
                    std::span<const T>(polyline.data(), polyline.size()), mesh.gdim, out);
            }
        }
    }

    out._num_cells = static_cast<int>(out._types.size());
    out._num_vertices = static_cast<int>(out._vertex_coords.size() / static_cast<std::size_t>(out._gdim));
    return out;
}

} // namespace

namespace cutcells::mesh
{

template <std::floating_point T>
CurvedCutMeshResult<T> cut_mesh_view_curved(
    const MeshView<T, int>&         mesh,
    const LevelSetFunction<T, int>& level_set,
    int                             geom_order,
    std::string_view                backend,
    int                             vis_subdivision,
    T                               tol)
{
    if (!mesh.has_cell_types())
        throw std::invalid_argument("cut_mesh_view_curved: mesh cell_types are required");
    if (!level_set.has_value())
        throw std::invalid_argument("cut_mesh_view_curved: level_set callable value(x) is required");
    if (geom_order < 1)
        throw std::invalid_argument("cut_mesh_view_curved: geom_order must be >= 1");
    if (vis_subdivision < 1)
        throw std::invalid_argument("cut_mesh_view_curved: vis_subdivision must be >= 1");

    const std::string backend_name = to_lower_ascii(backend);
    if (backend_name != "bernstein")
        throw std::invalid_argument("cut_mesh_view_curved: only backend='bernstein' is supported");

    if (mesh.gdim <= 0 || mesh.gdim > 3)
        throw std::invalid_argument("cut_mesh_view_curved: mesh.gdim must be in [1, 3]");

    CurvedCutMeshResult<T> result;
    init_curved_mesh(result.inside, mesh.gdim, 1, 0);
    init_curved_mesh(result.interface, mesh.gdim, geom_order, 0);
    init_curved_mesh(result.outside, mesh.gdim, 1, 0);

    const int n_cells = mesh.num_cells();
    result.local_meshes.reserve(static_cast<std::size_t>(n_cells));

    const int poly_order = std::max(geom_order, 2);

    for (int cell_id = 0; cell_id < n_cells; ++cell_id)
    {
        const auto ct = cell::map_vtk_type_to_cell_type(
            static_cast<cell::vtk_types>(mesh.cell_type(cell_id)));
        if (!bernstein_supported_cell(ct))
        {
            throw std::invalid_argument(
                "cut_mesh_view_curved: bernstein backend currently supports "
                "interval/triangle/tetrahedron/hexahedron");
        }

        const int p1_vertices = cell::get_num_vertices(ct);
        if (mesh.cell_num_nodes(cell_id) < p1_vertices)
            throw std::runtime_error("cut_mesh_view_curved: parent cell has too few nodes");

        std::vector<T> parent_coords_p1(static_cast<std::size_t>(p1_vertices * mesh.gdim), T(0));
        for (int i = 0; i < p1_vertices; ++i)
        {
            const int node_id = mesh.cell_node(cell_id, i);
            if (node_id < 0 || node_id >= mesh.num_nodes())
                throw std::runtime_error("cut_mesh_view_curved: invalid node id in connectivity");

            const T* x = mesh.node(node_id);
            for (int d = 0; d < mesh.gdim; ++d)
                parent_coords_p1[static_cast<std::size_t>(i * mesh.gdim + d)] = x[d];
        }

        std::vector<T> sampled_node_coords;
        std::vector<T> sampled_nodal_values;
        sample_callable_on_parent_nodes<T>(
            ct,
            std::span<const T>(parent_coords_p1.data(), parent_coords_p1.size()),
            mesh.gdim,
            poly_order,
            level_set,
            cell_id,
            sampled_node_coords,
            sampled_nodal_values);

        const auto domain = classify_parent_cell_bernstein<T>(
            ct,
            poly_order,
            std::span<const T>(sampled_nodal_values.data(), sampled_nodal_values.size()),
            tol);

        if (domain == ParentCellDomain::inside)
        {
            append_cell_to_curved_mesh<T>(
                result.inside,
                ct,
                std::span<const T>(parent_coords_p1.data(), parent_coords_p1.size()),
                mesh.gdim);
            result.n_parent_inside += 1;
            continue;
        }

        if (domain == ParentCellDomain::outside)
        {
            append_cell_to_curved_mesh<T>(
                result.outside,
                ct,
                std::span<const T>(parent_coords_p1.data(), parent_coords_p1.size()),
                mesh.gdim);
            result.n_parent_outside += 1;
            continue;
        }

        result.n_parent_intersected += 1;

        LevelSetFunction<T, int> local_level_set;
        local_level_set.nodal_values = std::span<const T>(sampled_nodal_values.data(), sampled_nodal_values.size());
        local_level_set.gdim = mesh.gdim;
        local_level_set.degree = poly_order;
        local_level_set.kind = LevelSetFunction<T, int>::Kind::fem_nodal;

        LocalMesh<T> local_mesh;
        init_local_mesh_from_cell<T>(
            local_mesh,
            std::span<const T>(parent_coords_p1.data(), parent_coords_p1.size()),
            ct,
            cell_id,
            1);

        decompose_local_mesh_with_backend<T, int>(
            local_mesh,
            local_level_set,
            LocalLevelSetBackend::bernstein,
            0,
            cell::edge_root::method::itp,
            true,
            6,
            tol,
            false);

        append_local_cells_by_domain<T>(local_mesh, cell::domain::inside, result.inside);
        append_local_cells_by_domain<T>(local_mesh, cell::domain::outside, result.outside);

        build_zero_entities<T>(local_mesh, 0);
        if (level_set.has_value() && level_set.has_gradient())
        {
            curve_zero_entities_with_backend<T, int>(
                local_mesh,
                level_set,
                LocalLevelSetBackend::analytical_callbacks,
                0,
                geom_order,
                tol);
        }
        else
        {
            curve_zero_entities_with_backend<T, int>(
                local_mesh,
                local_level_set,
                LocalLevelSetBackend::bernstein,
                0,
                geom_order,
                tol);
        }

        {
            const uint64_t zero_mask = uint64_t(1) << 0;
            if (local_mesh.tdim == 2)
                build_zero_chains(local_mesh, zero_mask);
            else if (local_mesh.tdim == 3)
                build_zero_patches(local_mesh, zero_mask);
        }

        result.local_meshes.push_back(std::move(local_mesh));
    }

    if (!result.local_meshes.empty())
    {
        std::vector<const LocalMesh<T>*> local_ptrs;
        local_ptrs.reserve(result.local_meshes.size());
        for (const LocalMesh<T>& local_mesh : result.local_meshes)
            local_ptrs.push_back(&local_mesh);

        result.interface = assemble_curved_interface_mesh<T>(
            std::span<const LocalMesh<T>*>(local_ptrs.data(), local_ptrs.size()),
            0,
            geom_order);
        result.interface_vis = linearize_local_interfaces<T>(
            std::span<const LocalMesh<T>*>(local_ptrs.data(), local_ptrs.size()),
            0,
            vis_subdivision);
    }

    result.inside_vis = curved_to_cut_mesh_direct<T>(result.inside);
    result.outside_vis = curved_to_cut_mesh_direct<T>(result.outside);
    if (result.interface_vis._offset.empty())
        result.interface_vis = linearize_curved_interface_mesh<T>(result.interface, vis_subdivision);

    return result;
}

template <std::floating_point T>
CurvedCutMeshResult<T> cut_vtk_mesh_curved(
    std::span<const T>              points,
    std::span<const int>            connectivity,
    std::span<const int>            offset,
    std::span<const int>            vtk_type,
    const LevelSetFunction<T, int>& level_set,
    int                             geom_order,
    std::string_view                backend,
    int                             vis_subdivision,
    T                               tol)
{
    if (offset.empty())
        throw std::invalid_argument("cut_vtk_mesh_curved: offset must not be empty");
    if (vtk_type.size() + 1 != offset.size())
        throw std::invalid_argument("cut_vtk_mesh_curved: offset size must be vtk_type size + 1");

    int max_node = -1;
    for (const int node : connectivity)
    {
        if (node < 0)
            throw std::invalid_argument("cut_vtk_mesh_curved: connectivity contains negative node index");
        max_node = std::max(max_node, node);
    }

    const int n_nodes = max_node + 1;
    if (n_nodes <= 0)
        throw std::invalid_argument("cut_vtk_mesh_curved: could not infer node count from connectivity");
    if (points.size() % static_cast<std::size_t>(n_nodes) != 0)
        throw std::invalid_argument("cut_vtk_mesh_curved: points size is not divisible by inferred node count");

    MeshView<T, int> mesh;
    mesh.gdim = static_cast<int>(points.size() / static_cast<std::size_t>(n_nodes));
    mesh.tdim = 0;
    mesh.coordinates = points;
    mesh.connectivity = connectivity;
    mesh.offsets = offset;
    mesh.cell_types = vtk_type;

    return cut_mesh_view_curved<T>(
        mesh,
        level_set,
        geom_order,
        backend,
        vis_subdivision,
        tol);
}

// explicit instantiations

template CurvedCutMeshResult<float> cut_mesh_view_curved<float>(
    const MeshView<float, int>&,
    const LevelSetFunction<float, int>&,
    int,
    std::string_view,
    int,
    float);

template CurvedCutMeshResult<double> cut_mesh_view_curved<double>(
    const MeshView<double, int>&,
    const LevelSetFunction<double, int>&,
    int,
    std::string_view,
    int,
    double);

template CurvedCutMeshResult<float> cut_vtk_mesh_curved<float>(
    std::span<const float>,
    std::span<const int>,
    std::span<const int>,
    std::span<const int>,
    const LevelSetFunction<float, int>&,
    int,
    std::string_view,
    int,
    float);

template CurvedCutMeshResult<double> cut_vtk_mesh_curved<double>(
    std::span<const double>,
    std::span<const int>,
    std::span<const int>,
    std::span<const int>,
    const LevelSetFunction<double, int>&,
    int,
    std::string_view,
    int,
    double);

} // namespace cutcells::mesh
