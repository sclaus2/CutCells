// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "curving.h"

#include "cell_topology.h"
#include "edge_root.h"
#include "geometric_quantity.h"
#include "mapping.h"
#include "reference_cell.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cutcells::curving
{
namespace
{

template <std::floating_point T>
bool solve_dense_small(std::vector<T> A, std::vector<T> b, int n, std::vector<T>& x);

template <std::floating_point T>
T local_dot(std::span<const T> a, std::span<const T> b)
{
    T value = T(0);
    for (std::size_t i = 0; i < a.size(); ++i)
        value += a[i] * b[i];
    return value;
}

template <std::floating_point T>
T local_norm(std::span<const T> a)
{
    return std::sqrt(local_dot<T>(a, a));
}

template <std::floating_point T>
std::array<T, 3> local_cross(std::span<const T> a, std::span<const T> b)
{
    return {a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
}

template <std::floating_point T>
std::span<const T> vertex_ref_point(const AdaptCell<T>& ac, int vertex_id)
{
    return std::span<const T>(
        ac.vertex_coords.data()
            + static_cast<std::size_t>(vertex_id)
                  * static_cast<std::size_t>(ac.tdim),
        static_cast<std::size_t>(ac.tdim));
}

template <class Int>
bool vertex_list_contains(std::span<const Int> vertices, int vertex)
{
    for (const auto v : vertices)
    {
        if (static_cast<int>(v) == vertex)
            return true;
    }
    return false;
}

template <std::floating_point T>
bool cell_contains_vertices(std::span<const std::int32_t> cell_vertices,
                            std::span<const int> query_vertices)
{
    for (const int v : query_vertices)
    {
        if (!vertex_list_contains<std::int32_t>(cell_vertices, v))
            return false;
    }
    return true;
}

template <std::floating_point T>
std::vector<int> local_zero_entity_vertex_ids(const AdaptCell<T>& ac,
                                              int local_zero_entity_id)
{
    if (local_zero_entity_id < 0
        || local_zero_entity_id >= ac.n_zero_entities())
    {
        return {};
    }

    const int zdim =
        ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    const int zid =
        ac.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim == 0)
        return {zid};

    auto verts = ac.entity_to_vertex[zdim][static_cast<std::int32_t>(zid)];
    std::vector<int> out;
    out.reserve(verts.size());
    for (const auto v : verts)
        out.push_back(static_cast<int>(v));
    return out;
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
    auto verts = ac.zero_entity_host_cell_vertices[
        static_cast<std::int32_t>(local_zero_entity_id)];
    std::vector<int> out;
    out.reserve(verts.size());
    for (const auto v : verts)
        out.push_back(static_cast<int>(v));
    return out;
}

template <std::floating_point T>
bool zero_entity_has_explicit_host_cell(const AdaptCell<T>& ac,
                                        int local_zero_entity_id)
{
    return !zero_entity_host_cell_vertex_ids<T>(ac, local_zero_entity_id).empty();
}

template <std::floating_point T>
void append_unique_host_vertex(const AdaptCell<T>& ac,
                               int vertex_id,
                               std::vector<int>& vertex_ids,
                               std::vector<T>& coords)
{
    if (vertex_list_contains<int>(
            std::span<const int>(vertex_ids.data(), vertex_ids.size()),
            vertex_id))
    {
        return;
    }
    vertex_ids.push_back(vertex_id);
    const auto p = vertex_ref_point<T>(ac, vertex_id);
    coords.insert(coords.end(), p.begin(), p.end());
}

template <std::floating_point T>
struct LocalHostDomain
{
    int dimension = -1;
    bool ordered_boundary = false;
    std::vector<T> vertices;

    bool valid(int tdim) const
    {
        if (dimension <= 0 || dimension > tdim)
            return false;
        const int nverts =
            static_cast<int>(vertices.size() / static_cast<std::size_t>(tdim));
        return nverts >= dimension + 1;
    }
};

template <std::floating_point T>
bool face_is_zero_for_mask(const AdaptCell<T>& ac,
                           std::span<const int> face_vertices,
                           std::uint64_t zero_mask)
{
    if (zero_mask == 0)
        return false;
    for (const int v : face_vertices)
    {
        if ((ac.zero_mask_per_vertex[static_cast<std::size_t>(v)] & zero_mask)
            != zero_mask)
        {
            return false;
        }
    }
    return true;
}

template <std::floating_point T>
LocalHostDomain<T> local_host_domain_for_zero_entity(
    const AdaptCell<T>& ac,
    int local_zero_entity_id)
{
    LocalHostDomain<T> host;
    if ((ac.tdim != 2 && ac.tdim != 3) || local_zero_entity_id < 0
        || local_zero_entity_id >= ac.n_zero_entities())
    {
        return host;
    }

    const int zdim =
        ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    const auto zverts =
        local_zero_entity_vertex_ids<T>(ac, local_zero_entity_id);
    if (zverts.empty())
        return host;

    const auto explicit_host_vertices =
        zero_entity_host_cell_vertex_ids<T>(ac, local_zero_entity_id);

    if (ac.tdim == 2 && zdim == 1 && !explicit_host_vertices.empty())
    {
        host.dimension = 2;
        host.ordered_boundary = true;
        for (const int vertex_id : explicit_host_vertices)
        {
            const auto p = vertex_ref_point<T>(ac, vertex_id);
            host.vertices.insert(host.vertices.end(), p.begin(), p.end());
        }
        return host;
    }

    if (zdim == 1
        && local_zero_entity_id
               < static_cast<int>(ac.zero_entity_host_face_id.size())
        && local_zero_entity_id
               < static_cast<int>(ac.zero_entity_host_cell_type.size())
        && !explicit_host_vertices.empty())
    {
        const int host_face_id =
            ac.zero_entity_host_face_id[static_cast<std::size_t>(local_zero_entity_id)];
        const cell::type host_type =
            ac.zero_entity_host_cell_type[static_cast<std::size_t>(local_zero_entity_id)];
        if (host_face_id >= 0
            && host_face_id < cell::num_faces(host_type))
        {
            auto face_vertices = cell::face_vertices(host_type, host_face_id);
            host.dimension = 2;
            host.ordered_boundary = true;
            for (const auto local_v : face_vertices)
            {
                if (local_v < 0
                    || local_v >= static_cast<int>(explicit_host_vertices.size()))
                {
                    host.vertices.clear();
                    host.dimension = -1;
                    return host;
                }
                const int vertex_id =
                    explicit_host_vertices[static_cast<std::size_t>(local_v)];
                const auto p = vertex_ref_point<T>(ac, vertex_id);
                host.vertices.insert(host.vertices.end(), p.begin(), p.end());
            }
            return host;
        }
    }

    if (zdim == 1 && !explicit_host_vertices.empty())
    {
        std::vector<int> ids;
        for (const int v : explicit_host_vertices)
            append_unique_host_vertex<T>(ac, v, ids, host.vertices);
        if (static_cast<int>(ids.size()) >= 4)
            host.dimension = 3;
        return host;
    }

    if (zdim == 2)
    {
        if (!explicit_host_vertices.empty())
        {
            std::vector<int> ids;
            for (const int v : explicit_host_vertices)
                append_unique_host_vertex<T>(ac, v, ids, host.vertices);
            if (static_cast<int>(ids.size()) >= 4)
                host.dimension = 3;
            return host;
        }

        const int n_cells = ac.n_entities(ac.tdim);
        std::vector<int> ids;
        for (int c = 0; c < n_cells; ++c)
        {
            const auto cell_vertices =
                ac.entity_to_vertex[ac.tdim][static_cast<std::int32_t>(c)];
            bool touches_zero_face = false;
            for (const int v : zverts)
                touches_zero_face = touches_zero_face
                                  || vertex_list_contains<std::int32_t>(
                                      cell_vertices, v);
            if (!touches_zero_face)
            {
                continue;
            }

            for (const auto v : cell_vertices)
                append_unique_host_vertex<T>(ac, static_cast<int>(v), ids, host.vertices);
        }

        if (static_cast<int>(ids.size()) >= 4)
            host.dimension = 3;
    }

    return host;
}

template <std::floating_point T>
bool clip_interval_by_halfspace(std::span<const T> x0,
                                std::span<const T> direction,
                                std::span<const T> normal,
                                T offset,
                                T tol,
                                T& lo,
                                T& hi)
{
    const T a = local_dot<T>(direction, normal);
    const T b = local_dot<T>(x0, normal) + offset;
    if (std::fabs(a) <= tol)
        return b >= -tol;

    const T t = (-tol - b) / a;
    if (a > T(0))
        lo = std::max(lo, t);
    else
        hi = std::min(hi, t);
    return lo <= hi;
}

template <std::floating_point T>
bool clip_line_interval_in_ordered_face(const LocalHostDomain<T>& host,
                                        int tdim,
                                        std::span<const T> x0,
                                        std::span<const T> direction,
                                        T tol,
                                        T& lo,
                                        T& hi)
{
    if (tdim != 3 || !host.ordered_boundary)
        return false;
    const int nverts =
        static_cast<int>(host.vertices.size() / static_cast<std::size_t>(tdim));
    if (nverts < 3)
        return false;

    const auto p0 = std::span<const T>(host.vertices.data(), 3);
    const auto p1 = std::span<const T>(host.vertices.data() + 3, 3);
    const auto p2 = std::span<const T>(host.vertices.data() + 6, 3);
    std::array<T, 3> e0 = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    std::array<T, 3> e1 = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
    auto normal = local_cross<T>(
        std::span<const T>(e0.data(), 3),
        std::span<const T>(e1.data(), 3));
    const T normal_norm =
        local_norm<T>(std::span<const T>(normal.data(), 3));
    const T direction_norm = local_norm<T>(direction);
    if (normal_norm <= tol || direction_norm <= tol)
        return false;
    const T plane_seed = (x0[0] - p0[0]) * normal[0]
                       + (x0[1] - p0[1]) * normal[1]
                       + (x0[2] - p0[2]) * normal[2];
    const T plane_dir = local_dot<T>(
        direction, std::span<const T>(normal.data(), 3));
    if (std::fabs(plane_seed) > tol * normal_norm
        || std::fabs(plane_dir) > tol * normal_norm * direction_norm)
    {
        return false;
    }

    std::array<T, 3> centroid = {T(0), T(0), T(0)};
    for (int i = 0; i < nverts; ++i)
    {
        const auto p = std::span<const T>(
            host.vertices.data() + static_cast<std::size_t>(3 * i), 3);
        for (int d = 0; d < 3; ++d)
            centroid[static_cast<std::size_t>(d)] += p[static_cast<std::size_t>(d)];
    }
    for (T& x : centroid)
        x /= static_cast<T>(nverts);

    for (int i = 0; i < nverts; ++i)
    {
        const auto a = std::span<const T>(
            host.vertices.data() + static_cast<std::size_t>(3 * i), 3);
        const auto b = std::span<const T>(
            host.vertices.data()
                + static_cast<std::size_t>(3 * ((i + 1) % nverts)), 3);
        std::array<T, 3> edge = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        auto inward = local_cross<T>(
            std::span<const T>(normal.data(), 3),
            std::span<const T>(edge.data(), 3));
        std::array<T, 3> ca = {
            centroid[0] - a[0], centroid[1] - a[1], centroid[2] - a[2]};
        if (local_dot<T>(
                std::span<const T>(inward.data(), 3),
                std::span<const T>(ca.data(), 3)) < T(0))
        {
            for (T& v : inward)
                v = -v;
        }
        const T offset = -local_dot<T>(
            std::span<const T>(inward.data(), 3), a);
        if (!clip_interval_by_halfspace<T>(
                x0,
                direction,
                std::span<const T>(inward.data(), 3),
                offset,
                tol,
                lo,
                hi))
        {
            return false;
        }
    }
    return lo <= hi;
}

template <std::floating_point T>
bool clip_line_interval_in_ordered_polygon_2d(const LocalHostDomain<T>& host,
                                              int tdim,
                                              std::span<const T> x0,
                                              std::span<const T> direction,
                                              T tol,
                                              T& lo,
                                              T& hi)
{
    if (tdim != 2 || !host.ordered_boundary)
        return false;
    const int nverts =
        static_cast<int>(host.vertices.size() / static_cast<std::size_t>(tdim));
    if (nverts < 3)
        return false;

    std::array<T, 2> centroid = {T(0), T(0)};
    for (int i = 0; i < nverts; ++i)
    {
        centroid[0] += host.vertices[static_cast<std::size_t>(2 * i)];
        centroid[1] += host.vertices[static_cast<std::size_t>(2 * i + 1)];
    }
    centroid[0] /= static_cast<T>(nverts);
    centroid[1] /= static_cast<T>(nverts);

    for (int i = 0; i < nverts; ++i)
    {
        const auto a = std::span<const T>(
            host.vertices.data() + static_cast<std::size_t>(2 * i), 2);
        const auto b = std::span<const T>(
            host.vertices.data()
                + static_cast<std::size_t>(2 * ((i + 1) % nverts)), 2);
        std::array<T, 2> edge = {b[0] - a[0], b[1] - a[1]};
        std::array<T, 2> inward = {-edge[1], edge[0]};
        std::array<T, 2> ca = {centroid[0] - a[0], centroid[1] - a[1]};
        if (local_dot<T>(
                std::span<const T>(inward.data(), 2),
                std::span<const T>(ca.data(), 2)) < T(0))
        {
            for (T& v : inward)
                v = -v;
        }
        const T offset = -local_dot<T>(
            std::span<const T>(inward.data(), 2), a);
        if (!clip_interval_by_halfspace<T>(
                x0,
                direction,
                std::span<const T>(inward.data(), 2),
                offset,
                tol,
                lo,
                hi))
        {
            return false;
        }
    }
    return lo <= hi;
}

template <std::floating_point T>
bool clip_line_interval_in_convex_hull(const LocalHostDomain<T>& host,
                                       int tdim,
                                       std::span<const T> x0,
                                       std::span<const T> direction,
                                       T tol,
                                       T& lo,
                                       T& hi)
{
    const int nverts =
        static_cast<int>(host.vertices.size() / static_cast<std::size_t>(tdim));
    if (host.dimension == 2)
    {
        if (tdim == 3)
            return clip_line_interval_in_ordered_face<T>(
                host, tdim, x0, direction, tol, lo, hi);
        if (tdim == 2)
            return clip_line_interval_in_ordered_polygon_2d<T>(
                host, tdim, x0, direction, tol, lo, hi);
        return false;
    }

    if (host.dimension != 3 || tdim != 3 || nverts < 4)
        return false;

    int accepted_planes = 0;
    for (int i = 0; i < nverts; ++i)
    {
        const auto pi = std::span<const T>(
            host.vertices.data() + static_cast<std::size_t>(3 * i), 3);
        for (int j = i + 1; j < nverts; ++j)
        {
            const auto pj = std::span<const T>(
                host.vertices.data() + static_cast<std::size_t>(3 * j), 3);
            for (int k = j + 1; k < nverts; ++k)
            {
                const auto pk = std::span<const T>(
                    host.vertices.data() + static_cast<std::size_t>(3 * k), 3);
                std::array<T, 3> e0 = {
                    pj[0] - pi[0], pj[1] - pi[1], pj[2] - pi[2]};
                std::array<T, 3> e1 = {
                    pk[0] - pi[0], pk[1] - pi[1], pk[2] - pi[2]};
                auto normal = local_cross<T>(
                    std::span<const T>(e0.data(), 3),
                    std::span<const T>(e1.data(), 3));
                const T nn = local_norm<T>(
                    std::span<const T>(normal.data(), 3));
                if (nn <= tol)
                    continue;

                bool has_pos = false;
                bool has_neg = false;
                for (int m = 0; m < nverts; ++m)
                {
                    const auto pm = std::span<const T>(
                        host.vertices.data() + static_cast<std::size_t>(3 * m), 3);
                    const T s = (pm[0] - pi[0]) * normal[0]
                              + (pm[1] - pi[1]) * normal[1]
                              + (pm[2] - pi[2]) * normal[2];
                    has_pos = has_pos || s > tol * nn;
                    has_neg = has_neg || s < -tol * nn;
                    if (has_pos && has_neg)
                        break;
                }
                if (has_pos && has_neg)
                    continue;
                if (!has_pos && !has_neg)
                    continue;
                if (has_neg)
                {
                    for (T& v : normal)
                        v = -v;
                }
                const T offset = -local_dot<T>(
                    std::span<const T>(normal.data(), 3), pi);
                if (!clip_interval_by_halfspace<T>(
                        x0,
                        direction,
                        std::span<const T>(normal.data(), 3),
                        offset,
                        tol,
                        lo,
                        hi))
                {
                    return false;
                }
                ++accepted_planes;
            }
        }
    }

    return accepted_planes > 0 && lo <= hi;
}

template <std::floating_point T>
std::vector<T> project_to_local_host_space(const AdaptCell<T>& ac,
                                           int local_zero_entity_id,
                                           std::span<const T> direction,
                                           T tol)
{
    std::vector<T> out(direction.begin(), direction.end());
    const auto host = local_host_domain_for_zero_entity<T>(
        ac, local_zero_entity_id);
    if (!host.valid(ac.tdim) || host.dimension == ac.tdim)
        return out;

    if (host.dimension == 2 && ac.tdim == 3)
    {
        const int nverts = static_cast<int>(
            host.vertices.size() / static_cast<std::size_t>(ac.tdim));
        if (nverts < 3)
            return out;
        const auto p0 = std::span<const T>(host.vertices.data(), 3);
        const auto p1 = std::span<const T>(host.vertices.data() + 3, 3);
        const auto p2 = std::span<const T>(host.vertices.data() + 6, 3);
        std::array<T, 3> e0 = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        std::array<T, 3> e1 = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
        auto normal = local_cross<T>(
            std::span<const T>(e0.data(), 3),
            std::span<const T>(e1.data(), 3));
        const T nn = local_dot<T>(
            std::span<const T>(normal.data(), 3),
            std::span<const T>(normal.data(), 3));
        if (nn <= tol * tol)
            return out;
        const T c = local_dot<T>(
            std::span<const T>(out.data(), out.size()),
            std::span<const T>(normal.data(), 3)) / nn;
        for (int d = 0; d < 3; ++d)
            out[static_cast<std::size_t>(d)] -= c * normal[static_cast<std::size_t>(d)];
    }

    return out;
}

template <std::floating_point T>
std::vector<T> local_host_face_normal_reference(const AdaptCell<T>& ac,
                                                int local_zero_entity_id,
                                                T tol)
{
    const auto host = local_host_domain_for_zero_entity<T>(
        ac, local_zero_entity_id);
    if (ac.tdim != 3 || !host.valid(ac.tdim) || host.dimension != 2)
        return {};

    const int nverts = static_cast<int>(
        host.vertices.size() / static_cast<std::size_t>(ac.tdim));
    if (nverts < 3)
        return {};

    const auto p0 = std::span<const T>(host.vertices.data(), 3);
    const auto p1 = std::span<const T>(host.vertices.data() + 3, 3);
    const auto p2 = std::span<const T>(host.vertices.data() + 6, 3);
    std::array<T, 3> e0 = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    std::array<T, 3> e1 = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
    auto normal = local_cross<T>(
        std::span<const T>(e0.data(), 3),
        std::span<const T>(e1.data(), 3));
    const T nn = local_norm<T>(std::span<const T>(normal.data(), 3));
    if (nn <= tol)
        return {};
    std::vector<T> out(3, T(0));
    for (int d = 0; d < 3; ++d)
        out[static_cast<std::size_t>(d)] = normal[static_cast<std::size_t>(d)] / nn;
    return out;
}

template <std::floating_point T>
struct BoundaryEdgeState
{
    std::array<int, 2> vertices = {-1, -1};
    const CurvedZeroEntityState<T>* state = nullptr;
    bool use_curved_state = true;
};

struct BoundaryEdgeRef
{
    int local_zero_entity_id = -1;
    int host_face_id = -1;
    bool use_curved_state = true;
};

template <std::floating_point T>
struct ProjectionStats
{
    int iterations = 0;
    CurvingStatus status = CurvingStatus::failed;
    CurvingFailureCode failure_code = CurvingFailureCode::projection_failed;
    T residual = std::numeric_limits<T>::infinity();
    std::uint32_t active_face_mask = 0;
    int closest_face_id = -1;
    int safe_subspace_dim = -1;
    CurvingProjectionMode projection_mode = CurvingProjectionMode::none;
    int retry_count = 0;
    std::vector<T> seed;
    std::vector<T> direction;
    T clip_lo = std::numeric_limits<T>::quiet_NaN();
    T clip_hi = std::numeric_limits<T>::quiet_NaN();
    T root_t = std::numeric_limits<T>::quiet_NaN();
};

template <std::floating_point T>
void append_node_stats(CurvedZeroEntityState<T>& state,
                       const ProjectionStats<T>& stats,
                       int tdim)
{
    state.node_iterations.push_back(static_cast<std::int32_t>(stats.iterations));
    state.node_status.push_back(static_cast<std::uint8_t>(stats.status));
    state.node_failure_code.push_back(static_cast<std::uint8_t>(stats.failure_code));
    state.node_residual.push_back(stats.residual);
    state.node_active_face_mask.push_back(stats.active_face_mask);
    state.node_closest_face_id.push_back(static_cast<std::int32_t>(stats.closest_face_id));
    state.node_safe_subspace_dim.push_back(static_cast<std::int32_t>(stats.safe_subspace_dim));
    state.node_projection_mode.push_back(static_cast<std::uint8_t>(stats.projection_mode));
    state.node_retry_count.push_back(static_cast<std::int32_t>(stats.retry_count));

    const T nan = std::numeric_limits<T>::quiet_NaN();
    auto append_vector = [&](std::vector<T>& dst, const std::vector<T>& src)
    {
        for (int d = 0; d < tdim; ++d)
        {
            dst.push_back(
                d < static_cast<int>(src.size())
                    ? src[static_cast<std::size_t>(d)]
                    : nan);
        }
    };
    append_vector(state.node_seed, stats.seed);
    append_vector(state.node_direction, stats.direction);
    state.node_clip_lo.push_back(stats.clip_lo);
    state.node_clip_hi.push_back(stats.clip_hi);
    state.node_root_t.push_back(stats.root_t);
}

template <std::floating_point T>
ProjectionStats<T> accepted_node_stats(CurvingFailureCode code, T residual = T(0))
{
    ProjectionStats<T> stats;
    stats.iterations = 0;
    stats.status = CurvingStatus::curved;
    stats.failure_code = code;
    stats.residual = residual;
    return stats;
}

inline bool same_unordered_vertices(std::span<const int> a, std::span<const int> b)
{
    if (a.size() != b.size())
        return false;
    for (const int av : a)
    {
        bool found = false;
        for (const int bv : b)
            found = found || (av == bv);
        if (!found)
            return false;
    }
    return true;
}

template <std::floating_point T>
std::vector<T> zero_entity_vertices(const AdaptCell<T>& ac, int local_zero_entity_id)
{
    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    const int zid = ac.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    std::vector<T> out;

    if (zdim == 0)
    {
        out.resize(static_cast<std::size_t>(ac.tdim));
        for (int d = 0; d < ac.tdim; ++d)
            out[static_cast<std::size_t>(d)] =
                ac.vertex_coords[static_cast<std::size_t>(zid * ac.tdim + d)];
        return out;
    }

    auto verts = ac.entity_to_vertex[zdim][static_cast<std::int32_t>(zid)];
    out.resize(static_cast<std::size_t>(verts.size() * ac.tdim));
    for (std::size_t i = 0; i < verts.size(); ++i)
    {
        const int v = static_cast<int>(verts[i]);
        for (int d = 0; d < ac.tdim; ++d)
            out[i * static_cast<std::size_t>(ac.tdim) + static_cast<std::size_t>(d)] =
                ac.vertex_coords[static_cast<std::size_t>(v * ac.tdim + d)];
    }
    return out;
}

template <std::floating_point T>
std::vector<int> zero_entity_vertex_ids(const AdaptCell<T>& ac, int local_zero_entity_id)
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

template <std::floating_point T>
T distance_between_vertices(const AdaptCell<T>& ac, int v0, int v1)
{
    T length2 = T(0);
    for (int d = 0; d < ac.tdim; ++d)
    {
        const T delta =
            ac.vertex_coords[static_cast<std::size_t>(v1 * ac.tdim + d)]
          - ac.vertex_coords[static_cast<std::size_t>(v0 * ac.tdim + d)];
        length2 += delta * delta;
    }
    return std::sqrt(length2);
}

template <std::floating_point T>
T zero_entity_edge_length(const AdaptCell<T>& ac, int local_zero_entity_id)
{
    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim != 1)
        return std::numeric_limits<T>::infinity();

    const int zid = ac.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    auto verts = ac.entity_to_vertex[1][static_cast<std::int32_t>(zid)];
    if (verts.size() != 2)
        return std::numeric_limits<T>::infinity();
    return distance_between_vertices<T>(
        ac, static_cast<int>(verts[0]), static_cast<int>(verts[1]));
}

template <std::floating_point T>
T zero_entity_face_max_edge_length(const AdaptCell<T>& ac, int local_zero_entity_id)
{
    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim != 2)
        return std::numeric_limits<T>::infinity();

    const int zid = ac.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    auto verts = ac.entity_to_vertex[2][static_cast<std::int32_t>(zid)];
    if (verts.size() < 2)
        return T(0);

    T max_length = T(0);
    for (std::size_t i = 0; i < verts.size(); ++i)
    {
        const int v0 = static_cast<int>(verts[i]);
        const int v1 = static_cast<int>(verts[(i + 1) % verts.size()]);
        max_length = std::max(
            max_length,
            distance_between_vertices<T>(
                ac,
                v0,
                v1));
    }
    return max_length;
}

template <std::floating_point T>
bool zero_entity_below_curving_threshold(const AdaptCell<T>& ac,
                                         int local_zero_entity_id,
                                         const CurvingOptions<T>& options)
{
    if (!(options.small_entity_tol > T(0)))
        return false;

    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim == 1)
        return zero_entity_edge_length<T>(ac, local_zero_entity_id)
            <= options.small_entity_tol;
    if (zdim == 2)
        return zero_entity_face_max_edge_length<T>(ac, local_zero_entity_id)
            <= options.small_entity_tol;
    return false;
}

template <std::floating_point T>
void append_straight_seed_node_stats(CurvedZeroEntityState<T>& state,
                                     const AdaptCell<T>& ac,
                                     std::span<const T> seeds)
{
    const int n_nodes = static_cast<int>(
        seeds.size() / static_cast<std::size_t>(std::max(ac.tdim, 1)));
    for (int i = 0; i < n_nodes; ++i)
        append_node_stats<T>(
            state,
            accepted_node_stats<T>(
                CurvingFailureCode::small_entity_kept_straight),
            ac.tdim);
}

template <std::floating_point T>
std::pair<T, T> legendre_value_and_derivative(int order, T x)
{
    if (order == 0)
        return {T(1), T(0)};

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
    const T derivative = T(order) * (pm2 - x * pm1) / denom;
    return {pm1, derivative};
}

template <std::floating_point T>
std::vector<T> gll_parameters(int order)
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
            const auto [p, dp] = legendre_value_and_derivative<T>(order, x);
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
std::vector<T> interpolation_parameters(int order, NodeFamily family)
{
    order = std::max(order, 1);
    std::vector<T> params(static_cast<std::size_t>(order + 1), T(0));

    if (family == NodeFamily::gll)
        return gll_parameters<T>(order);

    // The current reference Lagrange point generator is equispaced. Keep
    // lagrange distinct at the API/cache level so a later Basix-backed node
    // family can be added without changing cache semantics.
    for (int i = 0; i <= order; ++i)
        params[static_cast<std::size_t>(i)] = T(i) / T(order);
    return params;
}

template <std::floating_point T>
std::vector<std::array<T, 3>> triangle_interpolation_barycentric_nodes(
    int order,
    NodeFamily family);

template <std::floating_point T>
void append_edge_seed_nodes(const AdaptCell<T>& ac,
                            int local_zero_entity_id,
                            const CurvingOptions<T>& options,
                            std::vector<T>& seeds)
{
    const int zid = ac.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    auto verts = ac.entity_to_vertex[1][static_cast<std::int32_t>(zid)];
    if (verts.size() != 2)
        throw std::runtime_error("curving: zero edge does not have two vertices");

    const auto params = interpolation_parameters<T>(options.geometry_order, options.node_family);
    seeds.clear();
    seeds.reserve(params.size() * static_cast<std::size_t>(ac.tdim));
    for (const T s : params)
    {
        for (int d = 0; d < ac.tdim; ++d)
        {
            const T x0 = ac.vertex_coords[static_cast<std::size_t>(verts[0] * ac.tdim + d)];
            const T x1 = ac.vertex_coords[static_cast<std::size_t>(verts[1] * ac.tdim + d)];
            seeds.push_back((T(1) - s) * x0 + s * x1);
        }
    }
}

template <std::floating_point T>
void append_face_seed_nodes(const AdaptCell<T>& ac,
                            int local_zero_entity_id,
                            const CurvingOptions<T>& options,
                            std::vector<T>& seeds)
{
    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    const int zid = ac.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim != 2)
        throw std::runtime_error("curving: expected a zero face");

    auto verts = ac.entity_to_vertex[2][static_cast<std::int32_t>(zid)];
    const auto entity_type = ac.entity_types[2][static_cast<std::size_t>(zid)];
    const int order = std::max(options.geometry_order, 1);
    seeds.clear();

    auto append_combo = [&](std::span<const T> weights)
    {
        for (int d = 0; d < ac.tdim; ++d)
        {
            T x = T(0);
            for (std::size_t v = 0; v < weights.size(); ++v)
            {
                x += weights[v] * ac.vertex_coords[
                    static_cast<std::size_t>(verts[v] * ac.tdim + d)];
            }
            seeds.push_back(x);
        }
    };

    if (entity_type == cell::type::triangle)
    {
        const auto nodes =
            triangle_interpolation_barycentric_nodes<T>(order, options.node_family);
        for (const auto& w : nodes)
            append_combo(std::span<const T>(w.data(), 3));
        return;
    }

    if (entity_type == cell::type::quadrilateral)
    {
        const auto params = interpolation_parameters<T>(order, options.node_family);
        std::array<T, 4> w = {};
        for (const T v : params)
        {
            for (const T u : params)
            {
                w = {(T(1) - u) * (T(1) - v),
                     u * (T(1) - v),
                     (T(1) - u) * v,
                     u * v};
                append_combo(std::span<const T>(w.data(), 4));
            }
        }
        return;
    }

    throw std::runtime_error("curving: unsupported zero-face type");
}

template <std::floating_point T>
T lagrange_basis_1d(int i, std::span<const T> params, T x)
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
T warp_factor(int order, T r)
{
    if (order <= 1)
        return T(0);

    std::vector<T> equispaced(static_cast<std::size_t>(order + 1), T(0));
    std::vector<T> gll = gll_parameters<T>(order);
    for (int i = 0; i <= order; ++i)
    {
        equispaced[static_cast<std::size_t>(i)] = -T(1) + T(2 * i) / T(order);
        gll[static_cast<std::size_t>(i)] = T(2) * gll[static_cast<std::size_t>(i)] - T(1);
    }

    T warp = T(0);
    for (int i = 0; i <= order; ++i)
    {
        const T Li = lagrange_basis_1d<T>(
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
std::array<T, 3> equilateral_to_reference_barycentric(T x, T y)
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
std::vector<std::array<T, 3>> triangle_interpolation_barycentric_nodes(
    int order,
    NodeFamily family)
{
    order = std::max(order, 1);
    std::vector<std::array<T, 3>> nodes;
    nodes.reserve(static_cast<std::size_t>((order + 1) * (order + 2) / 2));

    if (family != NodeFamily::gll)
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
            const T warp1 = T(4) * L2 * L3 * warp_factor<T>(order, L3 - L2)
                          * (T(1) + (alpha * L1) * (alpha * L1));
            const T warp2 = T(4) * L1 * L3 * warp_factor<T>(order, L1 - L3)
                          * (T(1) + (alpha * L2) * (alpha * L2));
            const T warp3 = T(4) * L1 * L2 * warp_factor<T>(order, L2 - L1)
                          * (T(1) + (alpha * L3) * (alpha * L3));

            x += warp1 + cos120 * warp2 + cos240 * warp3;
            y += sin120 * warp2 + sin240 * warp3;
            nodes.push_back(equilateral_to_reference_barycentric<T>(x, y));
        }
    }
    return nodes;
}

template <std::floating_point T>
const BoundaryEdgeState<T>* find_boundary_edge_state(
    std::span<const BoundaryEdgeState<T>> boundary_edges,
    int v0,
    int v1)
{
    std::array<int, 2> query = {v0, v1};
    for (const auto& edge : boundary_edges)
    {
        if (same_unordered_vertices(
                std::span<const int>(query.data(), query.size()),
                std::span<const int>(edge.vertices.data(), edge.vertices.size())))
        {
            return &edge;
        }
    }
    return nullptr;
}

template <std::floating_point T>
std::vector<T> eval_edge_state_at(const BoundaryEdgeState<T>& edge,
                                  const CurvingOptions<T>& options,
                                  int from_vertex,
                                  int to_vertex,
                                  T t_from_to)
{
    if (edge.state == nullptr || edge.state->status != CurvingStatus::curved)
        throw std::runtime_error("curving: missing accepted curved boundary edge");

    const int order = std::max(options.geometry_order, 1);
    const int tdim = static_cast<int>(edge.state->ref_nodes.size()) / (order + 1);
    const bool same_orientation =
        edge.vertices[0] == from_vertex && edge.vertices[1] == to_vertex;
    const bool reverse_orientation =
        edge.vertices[0] == to_vertex && edge.vertices[1] == from_vertex;
    if (!same_orientation && !reverse_orientation)
        throw std::runtime_error("curving: boundary edge orientation mismatch");

    const T t = same_orientation ? t_from_to : T(1) - t_from_to;
    const auto params = interpolation_parameters<T>(order, options.node_family);
    std::vector<T> out(static_cast<std::size_t>(tdim), T(0));
    for (int i = 0; i <= order; ++i)
    {
        const T Li = lagrange_basis_1d<T>(
            i, std::span<const T>(params.data(), params.size()), t);
        for (int d = 0; d < tdim; ++d)
            out[static_cast<std::size_t>(d)] += Li * edge.state->ref_nodes[
                static_cast<std::size_t>(i * tdim + d)];
    }
    return out;
}

template <std::floating_point T>
void append_domain_inequality(std::vector<std::array<T, 4>>& ineq,
                              std::array<T, 3> a,
                              T b)
{
    ineq.push_back({a[0], a[1], a[2], b});
}

template <std::floating_point T>
std::vector<std::array<T, 4>> reference_domain_inequalities(cell::type cell_type)
{
    std::vector<std::array<T, 4>> ineq;
    switch (cell_type)
    {
    case cell::type::interval:
        append_domain_inequality<T>(ineq, {-T(1), T(0), T(0)}, T(0));
        append_domain_inequality<T>(ineq, { T(1), T(0), T(0)}, T(1));
        break;
    case cell::type::triangle:
        append_domain_inequality<T>(ineq, {-T(1), T(0), T(0)}, T(0));
        append_domain_inequality<T>(ineq, {T(0), -T(1), T(0)}, T(0));
        append_domain_inequality<T>(ineq, {T(1), T(1), T(0)}, T(1));
        break;
    case cell::type::tetrahedron:
        append_domain_inequality<T>(ineq, {-T(1), T(0), T(0)}, T(0));
        append_domain_inequality<T>(ineq, {T(0), -T(1), T(0)}, T(0));
        append_domain_inequality<T>(ineq, {T(0), T(0), -T(1)}, T(0));
        append_domain_inequality<T>(ineq, {T(1), T(1), T(1)}, T(1));
        break;
    case cell::type::quadrilateral:
        append_domain_inequality<T>(ineq, {-T(1), T(0), T(0)}, T(0));
        append_domain_inequality<T>(ineq, { T(1), T(0), T(0)}, T(1));
        append_domain_inequality<T>(ineq, {T(0), -T(1), T(0)}, T(0));
        append_domain_inequality<T>(ineq, {T(0),  T(1), T(0)}, T(1));
        break;
    case cell::type::hexahedron:
        for (int d = 0; d < 3; ++d)
        {
            std::array<T, 3> lo = {T(0), T(0), T(0)};
            std::array<T, 3> hi = {T(0), T(0), T(0)};
            lo[static_cast<std::size_t>(d)] = -T(1);
            hi[static_cast<std::size_t>(d)] = T(1);
            append_domain_inequality<T>(ineq, lo, T(0));
            append_domain_inequality<T>(ineq, hi, T(1));
        }
        break;
    case cell::type::prism:
        append_domain_inequality<T>(ineq, {-T(1), T(0), T(0)}, T(0));
        append_domain_inequality<T>(ineq, {T(0), -T(1), T(0)}, T(0));
        append_domain_inequality<T>(ineq, {T(1), T(1), T(0)}, T(1));
        append_domain_inequality<T>(ineq, {T(0), T(0), -T(1)}, T(0));
        append_domain_inequality<T>(ineq, {T(0), T(0), T(1)}, T(1));
        break;
    case cell::type::pyramid:
        append_domain_inequality<T>(ineq, {-T(1), T(0), T(0)}, T(0));
        append_domain_inequality<T>(ineq, {T(0), -T(1), T(0)}, T(0));
        append_domain_inequality<T>(ineq, {T(1), T(0), T(1)}, T(1));
        append_domain_inequality<T>(ineq, {T(0), T(1), T(1)}, T(1));
        append_domain_inequality<T>(ineq, {T(0), T(0), -T(1)}, T(0));
        append_domain_inequality<T>(ineq, {T(0), T(0), T(1)}, T(1));
        break;
    default:
        break;
    }
    return ineq;
}

template <std::floating_point T>
geom::ParentEntity zero_entity_parent_entity(const AdaptCell<T>& ac,
                                             int local_zero_entity_id)
{
    if (local_zero_entity_id < 0
        || local_zero_entity_id >= ac.n_zero_entities())
    {
        return {-1, -1};
    }

    const int parent_dim =
        ac.zero_entity_parent_dim[static_cast<std::size_t>(local_zero_entity_id)];
    int parent_id =
        ac.zero_entity_parent_id[static_cast<std::size_t>(local_zero_entity_id)];
    if (parent_dim == ac.tdim)
        parent_id = -1;
    return {parent_dim, parent_id};
}

template <std::floating_point T>
bool host_parameter_interval(const AdaptCell<T>& ac,
                             int local_zero_entity_id,
                             std::span<const T> x0,
                             std::span<const T> d,
                             T tol,
                             T& lo,
                             T& hi)
{
    const auto local_host =
        local_host_domain_for_zero_entity<T>(ac, local_zero_entity_id);
    const bool has_explicit_host =
        zero_entity_has_explicit_host_cell<T>(ac, local_zero_entity_id);
    if (local_host.valid(ac.tdim))
    {
        T local_lo = -std::numeric_limits<T>::infinity();
        T local_hi = std::numeric_limits<T>::infinity();
        if (clip_line_interval_in_convex_hull<T>(
                local_host, ac.tdim, x0, d, tol, local_lo, local_hi)
            && local_lo <= local_hi)
        {
            lo = local_lo;
            hi = local_hi;
            return true;
        }
        if (has_explicit_host)
            return false;
    }
    else if (has_explicit_host)
    {
        return false;
    }

    const auto host = zero_entity_parent_entity<T>(ac, local_zero_entity_id);
    const auto interval = geom::clip_line_interval_in_parent_entity<T>(
        ac.parent_cell_type,
        host,
        x0,
        d,
        -std::numeric_limits<T>::infinity(),
        std::numeric_limits<T>::infinity(),
        tol);
    lo = interval.t0;
    hi = interval.t1;
    return interval.valid && lo <= hi;
}

template <std::floating_point T>
T vec_dot(std::span<const T> a, std::span<const T> b)
{
    T out = T(0);
    for (std::size_t i = 0; i < a.size(); ++i)
        out += a[i] * b[i];
    return out;
}

template <std::floating_point T>
T vec_norm(std::span<const T> a)
{
    return std::sqrt(vec_dot<T>(a, a));
}

template <std::floating_point T>
std::vector<T> normalized_delta_direction(std::span<const T> from,
                                          std::span<const T> to,
                                          T tol)
{
    std::vector<T> direction(from.size(), T(0));
    for (std::size_t i = 0; i < from.size(); ++i)
        direction[i] = to[i] - from[i];

    const T norm = vec_norm<T>(
        std::span<const T>(direction.data(), direction.size()));
    if (norm <= tol)
        return {};

    for (T& value : direction)
        value /= norm;
    return direction;
}

template <std::floating_point T>
T face_normal_norm(const std::array<T, 4>& row, int tdim)
{
    T n2 = T(0);
    for (int d = 0; d < tdim; ++d)
        n2 += row[static_cast<std::size_t>(d)] * row[static_cast<std::size_t>(d)];
    return std::sqrt(n2);
}

template <std::floating_point T>
T face_value(const std::array<T, 4>& row, std::span<const T> x, int tdim)
{
    T ax = T(0);
    for (int d = 0; d < tdim; ++d)
        ax += row[static_cast<std::size_t>(d)] * x[static_cast<std::size_t>(d)];
    return ax;
}

template <std::floating_point T>
T normalized_face_slack(const std::array<T, 4>& row,
                        std::span<const T> x,
                        int tdim)
{
    const T n = face_normal_norm<T>(row, tdim);
    if (n <= T(0))
        return std::numeric_limits<T>::infinity();
    return (row[3] - face_value<T>(row, x, tdim)) / n;
}

template <std::floating_point T>
bool vertex_satisfies_face(const std::array<T, 4>& row,
                           std::span<const T> ref_vertices,
                           int vertex,
                           int tdim)
{
    std::array<T, 3> x = {T(0), T(0), T(0)};
    for (int d = 0; d < tdim; ++d)
        x[static_cast<std::size_t>(d)] =
            ref_vertices[static_cast<std::size_t>(vertex * tdim + d)];
    const T n = face_normal_norm<T>(row, tdim);
    const T tol = T(64) * std::numeric_limits<T>::epsilon() * std::max(T(1), n);
    return std::fabs(face_value<T>(row, std::span<const T>(x.data(), static_cast<std::size_t>(tdim)), tdim) - row[3]) <= tol;
}

template <std::floating_point T>
int face_inequality_index_from_vertices(cell::type cell_type,
                                        std::span<const int> vertices,
                                        int tdim)
{
    const auto ineq = reference_domain_inequalities<T>(cell_type);
    const auto ref_vertices = cell::reference_vertices<T>(cell_type);
    for (int i = 0; i < static_cast<int>(ineq.size()); ++i)
    {
        bool all_on_face = true;
        for (const int v : vertices)
        {
            all_on_face = all_on_face
                       && vertex_satisfies_face<T>(
                              ineq[static_cast<std::size_t>(i)],
                              std::span<const T>(ref_vertices.data(), ref_vertices.size()),
                              v,
                              tdim);
        }
        if (all_on_face)
            return i;
    }
    return -1;
}

template <std::floating_point T>
int parent_face_inequality_index(cell::type cell_type, int parent_face_id, int tdim)
{
    if (parent_face_id < 0 || parent_face_id >= cell::num_faces(cell_type))
        return -1;
    const auto fv = cell::face_vertices(cell_type, parent_face_id);
    return face_inequality_index_from_vertices<T>(cell_type, fv, tdim);
}

template <std::floating_point T>
std::uint32_t structural_active_face_mask(const AdaptCell<T>& ac,
                                          int local_zero_entity_id)
{
    if (ac.tdim != 3)
        return 0;

    const int parent_dim =
        ac.zero_entity_parent_dim[static_cast<std::size_t>(local_zero_entity_id)];
    const int parent_id =
        ac.zero_entity_parent_id[static_cast<std::size_t>(local_zero_entity_id)];

    std::uint32_t mask = 0;
    if (parent_dim == 2)
    {
        const int row = parent_face_inequality_index<T>(
            ac.parent_cell_type, parent_id, ac.tdim);
        if (row >= 0)
            mask |= std::uint32_t(1) << row;
        return mask;
    }

    const auto ineq = reference_domain_inequalities<T>(ac.parent_cell_type);
    const auto ref_vertices = cell::reference_vertices<T>(ac.parent_cell_type);
    if (parent_dim == 1)
    {
        const auto edges = cell::edges(ac.parent_cell_type);
        if (parent_id < 0 || parent_id >= static_cast<int>(edges.size()))
            return mask;
        const auto edge = edges[static_cast<std::size_t>(parent_id)];
        for (int i = 0; i < static_cast<int>(ineq.size()); ++i)
        {
            if (vertex_satisfies_face<T>(
                    ineq[static_cast<std::size_t>(i)],
                    std::span<const T>(ref_vertices.data(), ref_vertices.size()),
                    edge[0],
                    ac.tdim)
                && vertex_satisfies_face<T>(
                    ineq[static_cast<std::size_t>(i)],
                    std::span<const T>(ref_vertices.data(), ref_vertices.size()),
                    edge[1],
                    ac.tdim))
            {
                mask |= std::uint32_t(1) << i;
            }
        }
    }
    else if (parent_dim == 0)
    {
        for (int i = 0; i < static_cast<int>(ineq.size()); ++i)
        {
            if (vertex_satisfies_face<T>(
                    ineq[static_cast<std::size_t>(i)],
                    std::span<const T>(ref_vertices.data(), ref_vertices.size()),
                    parent_id,
                    ac.tdim))
            {
                mask |= std::uint32_t(1) << i;
            }
        }
    }
    return mask;
}

template <std::floating_point T>
std::uint32_t add_near_active_faces(const AdaptCell<T>& ac,
                                    std::span<const T> x,
                                    std::uint32_t mask,
                                    T tol)
{
    const auto ineq = reference_domain_inequalities<T>(ac.parent_cell_type);
    for (int i = 0; i < static_cast<int>(ineq.size()); ++i)
    {
        const T slack = normalized_face_slack<T>(
            ineq[static_cast<std::size_t>(i)], x, ac.tdim);
        if (slack <= tol)
            mask |= std::uint32_t(1) << i;
    }
    return mask;
}

template <std::floating_point T>
int closest_inactive_face(const AdaptCell<T>& ac,
                          std::span<const T> x,
                          std::uint32_t active_mask)
{
    const auto ineq = reference_domain_inequalities<T>(ac.parent_cell_type);
    int best = -1;
    T best_slack = std::numeric_limits<T>::infinity();
    for (int i = 0; i < static_cast<int>(ineq.size()); ++i)
    {
        if ((active_mask & (std::uint32_t(1) << i)) != 0)
            continue;
        const T slack = normalized_face_slack<T>(
            ineq[static_cast<std::size_t>(i)], x, ac.tdim);
        if (slack < best_slack)
        {
            best_slack = slack;
            best = i;
        }
    }
    return best;
}

template <std::floating_point T>
std::uint32_t add_closest_faces(const AdaptCell<T>& ac,
                                std::span<const T> x,
                                std::uint32_t active_mask,
                                int& closest_face_id)
{
    const auto ineq = reference_domain_inequalities<T>(ac.parent_cell_type);
    closest_face_id = closest_inactive_face<T>(ac, x, active_mask);
    if (closest_face_id < 0)
        return active_mask;

    const T best_slack = normalized_face_slack<T>(
        ineq[static_cast<std::size_t>(closest_face_id)], x, ac.tdim);
    const T tie_tol = std::max(T(128) * std::numeric_limits<T>::epsilon(),
                              T(10) * ac.tdim * std::numeric_limits<T>::epsilon());
    std::uint32_t out = active_mask;
    for (int i = 0; i < static_cast<int>(ineq.size()); ++i)
    {
        if ((active_mask & (std::uint32_t(1) << i)) != 0)
            continue;
        const T slack = normalized_face_slack<T>(
            ineq[static_cast<std::size_t>(i)], x, ac.tdim);
        if (std::fabs(slack - best_slack) <= tie_tol * std::max(T(1), std::fabs(best_slack)))
            out |= std::uint32_t(1) << i;
    }
    return out;
}

template <std::floating_point T>
std::vector<std::vector<T>> orthonormal_active_normals(const AdaptCell<T>& ac,
                                                       std::uint32_t active_mask)
{
    const auto ineq = reference_domain_inequalities<T>(ac.parent_cell_type);
    std::vector<std::vector<T>> normals;
    for (int i = 0; i < static_cast<int>(ineq.size()); ++i)
    {
        if ((active_mask & (std::uint32_t(1) << i)) == 0)
            continue;
        std::vector<T> q(static_cast<std::size_t>(ac.tdim), T(0));
        for (int d = 0; d < ac.tdim; ++d)
            q[static_cast<std::size_t>(d)] =
                ineq[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
        for (const auto& prev : normals)
        {
            const T c = vec_dot<T>(
                std::span<const T>(q.data(), q.size()),
                std::span<const T>(prev.data(), prev.size()));
            for (int d = 0; d < ac.tdim; ++d)
                q[static_cast<std::size_t>(d)] -= c * prev[static_cast<std::size_t>(d)];
        }
        const T n = vec_norm<T>(std::span<const T>(q.data(), q.size()));
        if (n <= T(64) * std::numeric_limits<T>::epsilon())
            continue;
        for (T& x : q)
            x /= n;
        normals.push_back(q);
    }
    return normals;
}

template <std::floating_point T>
std::vector<T> project_to_active_face_space(const AdaptCell<T>& ac,
                                            std::span<const T> direction,
                                            std::uint32_t active_mask)
{
    std::vector<T> out(direction.begin(), direction.end());
    const auto normals = orthonormal_active_normals<T>(ac, active_mask);
    for (const auto& n : normals)
    {
        const T c = vec_dot<T>(
            std::span<const T>(out.data(), out.size()),
            std::span<const T>(n.data(), n.size()));
        for (int d = 0; d < ac.tdim; ++d)
            out[static_cast<std::size_t>(d)] -= c * n[static_cast<std::size_t>(d)];
    }
    return out;
}

template <std::floating_point T>
std::vector<std::vector<T>> active_nullspace_basis(const AdaptCell<T>& ac,
                                                   std::uint32_t active_mask)
{
    std::vector<std::vector<T>> basis;
    const auto normals = orthonormal_active_normals<T>(ac, active_mask);
    for (int axis = 0; axis < ac.tdim; ++axis)
    {
        std::vector<T> q(static_cast<std::size_t>(ac.tdim), T(0));
        q[static_cast<std::size_t>(axis)] = T(1);
        for (const auto& n : normals)
        {
            const T c = vec_dot<T>(
                std::span<const T>(q.data(), q.size()),
                std::span<const T>(n.data(), n.size()));
            for (int d = 0; d < ac.tdim; ++d)
                q[static_cast<std::size_t>(d)] -= c * n[static_cast<std::size_t>(d)];
        }
        for (const auto& prev : basis)
        {
            const T c = vec_dot<T>(
                std::span<const T>(q.data(), q.size()),
                std::span<const T>(prev.data(), prev.size()));
            for (int d = 0; d < ac.tdim; ++d)
                q[static_cast<std::size_t>(d)] -= c * prev[static_cast<std::size_t>(d)];
        }
        const T n = vec_norm<T>(std::span<const T>(q.data(), q.size()));
        if (n <= T(64) * std::numeric_limits<T>::epsilon())
            continue;
        for (T& x : q)
            x /= n;
        basis.push_back(q);
    }
    return basis;
}

template <std::floating_point T>
int active_subspace_dim(const AdaptCell<T>& ac, std::uint32_t active_mask)
{
    return static_cast<int>(active_nullspace_basis<T>(ac, active_mask).size());
}

template <std::floating_point T>
void append_unique_direction(std::vector<std::vector<T>>& directions,
                             std::vector<T> dir,
                             T tol)
{
    const T n = vec_norm<T>(std::span<const T>(dir.data(), dir.size()));
    if (n <= tol)
        return;
    for (T& x : dir)
        x /= n;

    for (const auto& existing : directions)
    {
        const T c = std::fabs(vec_dot<T>(
            std::span<const T>(dir.data(), dir.size()),
            std::span<const T>(existing.data(), existing.size())));
        if (c >= T(1) - T(64) * std::numeric_limits<T>::epsilon())
            return;
    }
    directions.push_back(std::move(dir));
}

template <std::floating_point T>
void append_admissible_unique_direction(std::vector<std::vector<T>>& directions,
                                        const AdaptCell<T>& ac,
                                        int local_zero_entity_id,
                                        std::span<const T> raw_direction,
                                        T tol)
{
    const auto projected = geom::admissible_direction_in_parent_frame<T>(
        ac.parent_cell_type,
        zero_entity_parent_entity<T>(ac, local_zero_entity_id),
        raw_direction,
        tol);
    if (projected.degenerate())
        return;
    auto local_projected = project_to_local_host_space<T>(
        ac,
        local_zero_entity_id,
        std::span<const T>(projected.value.data(), projected.value.size()),
        tol);
    append_unique_direction<T>(directions, std::move(local_projected), tol);
}

template <std::floating_point T>
std::vector<T> zero_entity_edge_tangent(const AdaptCell<T>& ac,
                                        int local_zero_entity_id)
{
    std::vector<T> tangent(static_cast<std::size_t>(ac.tdim), T(0));
    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim != 1)
        return tangent;

    const int zid = ac.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    auto verts = ac.entity_to_vertex[1][static_cast<std::int32_t>(zid)];
    if (verts.size() != 2)
        return tangent;
    for (int d = 0; d < ac.tdim; ++d)
        tangent[static_cast<std::size_t>(d)] =
            ac.vertex_coords[static_cast<std::size_t>(verts[1] * ac.tdim + d)]
          - ac.vertex_coords[static_cast<std::size_t>(verts[0] * ac.tdim + d)];
    return tangent;
}

template <std::floating_point T>
std::vector<T> zero_entity_face_normal(const AdaptCell<T>& ac,
                                       int local_zero_entity_id)
{
    std::vector<T> normal(static_cast<std::size_t>(ac.tdim), T(0));
    if (ac.tdim != 3)
        return normal;

    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim != 2)
        return normal;

    const int zid = ac.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    auto verts = ac.entity_to_vertex[2][static_cast<std::int32_t>(zid)];
    if (verts.size() < 3)
        return normal;

    std::array<T, 3> e0 = {T(0), T(0), T(0)};
    std::array<T, 3> e1 = {T(0), T(0), T(0)};
    for (int d = 0; d < 3; ++d)
    {
        e0[static_cast<std::size_t>(d)] =
            ac.vertex_coords[static_cast<std::size_t>(verts[1] * ac.tdim + d)]
          - ac.vertex_coords[static_cast<std::size_t>(verts[0] * ac.tdim + d)];
        e1[static_cast<std::size_t>(d)] =
            ac.vertex_coords[static_cast<std::size_t>(verts[2] * ac.tdim + d)]
          - ac.vertex_coords[static_cast<std::size_t>(verts[0] * ac.tdim + d)];
    }
    normal[0] = e0[1] * e1[2] - e0[2] * e1[1];
    normal[1] = e0[2] * e1[0] - e0[0] * e1[2];
    normal[2] = e0[0] * e1[1] - e0[1] * e1[0];
    return normal;
}

template <std::floating_point T>
bool contains_vertex(std::span<const std::int32_t> vertices, std::int32_t query)
{
    for (const auto vertex : vertices)
        if (vertex == query)
            return true;
    return false;
}

template <std::floating_point T>
std::vector<T> adjacent_zero_face_normal_for_edge(const AdaptCell<T>& ac,
                                                  int local_zero_entity_id,
                                                  T tol)
{
    std::vector<T> best(static_cast<std::size_t>(ac.tdim), T(0));
    if (ac.tdim != 3)
        return best;

    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim != 1)
        return best;

    const int zid = ac.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    auto edge_vertices = ac.entity_to_vertex[1][static_cast<std::int32_t>(zid)];
    if (edge_vertices.size() != 2)
        return best;

    const std::uint64_t edge_mask =
        ac.zero_entity_zero_mask[static_cast<std::size_t>(local_zero_entity_id)];
    T best_norm = T(0);
    for (int z = 0; z < ac.n_zero_entities(); ++z)
    {
        if (ac.zero_entity_dim[static_cast<std::size_t>(z)] != 2)
            continue;

        const std::uint64_t face_mask =
            ac.zero_entity_zero_mask[static_cast<std::size_t>(z)];
        if ((face_mask & edge_mask) != face_mask)
            continue;

        const int face_id = ac.zero_entity_id[static_cast<std::size_t>(z)];
        auto face_vertices = ac.entity_to_vertex[2][static_cast<std::int32_t>(face_id)];
        if (!contains_vertex<T>(face_vertices, edge_vertices[0])
            || !contains_vertex<T>(face_vertices, edge_vertices[1]))
        {
            continue;
        }

        auto normal = zero_entity_face_normal<T>(ac, z);
        const T n = vec_norm<T>(std::span<const T>(normal.data(), normal.size()));
        if (n > best_norm)
        {
            best_norm = n;
            best = std::move(normal);
        }
    }

    if (best_norm <= tol)
        std::fill(best.begin(), best.end(), T(0));
    return best;
}

template <std::floating_point T>
geom::VectorQuantity<T> straight_zero_entity_normal_direction(
    const AdaptCell<T>& ac,
    int local_zero_entity_id,
    T tol)
{
    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    const int tdim = ac.tdim;
    const auto host = zero_entity_parent_entity<T>(ac, local_zero_entity_id);

    if (zdim == 1)
    {
        const int zid = ac.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
        auto verts = ac.entity_to_vertex[1][static_cast<std::int32_t>(zid)];
        if (verts.size() != 2)
            return {{}, T(0), geom::Degeneracy::invalid_parent_entity};

        std::vector<T> a(static_cast<std::size_t>(tdim), T(0));
        std::vector<T> b(static_cast<std::size_t>(tdim), T(0));
        for (int d = 0; d < tdim; ++d)
        {
            a[static_cast<std::size_t>(d)] =
                ac.vertex_coords[static_cast<std::size_t>(verts[0] * tdim + d)];
            b[static_cast<std::size_t>(d)] =
                ac.vertex_coords[static_cast<std::size_t>(verts[1] * tdim + d)];
        }

        geom::VectorQuantity<T> raw;
        if (tdim == 2)
        {
            raw = geom::segment_normal<T>(
                std::span<const T>(a.data(), a.size()),
                std::span<const T>(b.data(), b.size()),
                true,
                tol);
        }
        else if (tdim == 3)
        {
            auto explicit_host_normal =
                local_host_face_normal_reference<T>(
                    ac, local_zero_entity_id, tol);
            if (!explicit_host_normal.empty())
            {
                raw = geom::in_face_segment_normal<T>(
                    std::span<const T>(a.data(), a.size()),
                    std::span<const T>(b.data(), b.size()),
                    std::span<const T>(
                        explicit_host_normal.data(), explicit_host_normal.size()),
                    true,
                    tol);
            }
            else if (host.dim == 2)
            {
                const auto face_normal = geom::parent_face_normal<T>(
                    ac.parent_cell_type, host.id, true, tol);
                if (face_normal.degenerate())
                    return {{}, T(0), face_normal.degeneracy};
                raw = geom::in_face_segment_normal<T>(
                    std::span<const T>(a.data(), a.size()),
                    std::span<const T>(b.data(), b.size()),
                    std::span<const T>(face_normal.value.data(), face_normal.value.size()),
                    true,
                    tol);
            }
            else
            {
                auto face_normal = adjacent_zero_face_normal_for_edge<T>(
                    ac, local_zero_entity_id, tol);
                auto face_normal_quantity = geom::make_vector_quantity<T>(
                    std::move(face_normal),
                    tol,
                    geom::Degeneracy::zero_projection);
                if (!face_normal_quantity.degenerate())
                {
                    raw = geom::in_face_segment_normal<T>(
                        std::span<const T>(a.data(), a.size()),
                        std::span<const T>(b.data(), b.size()),
                        std::span<const T>(
                            face_normal_quantity.value.data(),
                            face_normal_quantity.value.size()),
                        true,
                        tol);
                }
                else
                {
                    raw = std::move(face_normal_quantity);
                }
            }
        }
        else
        {
            return {std::vector<T>(static_cast<std::size_t>(tdim), T(0)),
                    T(0),
                    geom::Degeneracy::zero_projection};
        }

        if (raw.degenerate())
            return raw;
        return geom::admissible_direction_in_parent_frame<T>(
            ac.parent_cell_type,
            host,
            std::span<const T>(raw.value.data(), raw.value.size()),
            tol);
    }

    if (zdim == 2 && tdim == 3)
    {
        auto normal = zero_entity_face_normal<T>(ac, local_zero_entity_id);
        return geom::admissible_direction_in_parent_frame<T>(
            ac.parent_cell_type,
            host,
            std::span<const T>(normal.data(), normal.size()),
            tol);
    }

    return {std::vector<T>(static_cast<std::size_t>(tdim), T(0)),
            T(0),
            geom::Degeneracy::zero_projection};
}

template <std::floating_point T>
std::vector<T> remove_component(std::vector<T> direction,
                                std::span<const T> tangent)
{
    const T tt = vec_dot<T>(tangent, tangent);
    if (tt <= T(0))
        return direction;
    const T c = vec_dot<T>(
        std::span<const T>(direction.data(), direction.size()), tangent) / tt;
    for (std::size_t i = 0; i < direction.size(); ++i)
        direction[i] -= c * tangent[i];
    return direction;
}

template <std::floating_point T>
void orient_with_level_set_gradient(std::vector<T>& direction,
                                    std::span<const T> grad)
{
    if (direction.size() != grad.size())
        return;
    if (vec_dot<T>(
            std::span<const T>(direction.data(), direction.size()),
            grad) < T(0))
    {
        for (T& value : direction)
            value = -value;
    }
}

template <std::floating_point T, std::integral I>
std::vector<T> physical_normal_reference_direction(
    const LevelSetCell<T, I>& ls_cell,
    std::span<const T> grad_ref)
{
    std::vector<T> out(grad_ref.begin(), grad_ref.end());
    const int tdim = ls_cell.tdim;
    const int gdim = ls_cell.gdim;
    if (tdim <= 0 || tdim > 3 || gdim <= 0
        || static_cast<int>(grad_ref.size()) != tdim
        || ls_cell.parent_vertex_coords.empty())
    {
        return out;
    }

    const auto cols = cell::jacobian_col_indices(ls_cell.cell_type);
    std::vector<T> gram(static_cast<std::size_t>(tdim * tdim), T(0));
    for (int a = 0; a < tdim; ++a)
    {
        const int va = cols[static_cast<std::size_t>(a)];
        if (va < 0)
            return out;
        for (int b = 0; b < tdim; ++b)
        {
            const int vb = cols[static_cast<std::size_t>(b)];
            if (vb < 0)
                return out;
            T value = T(0);
            for (int r = 0; r < gdim; ++r)
            {
                const T ja = ls_cell.parent_vertex_coords[
                    static_cast<std::size_t>(va * gdim + r)]
                           - ls_cell.parent_vertex_coords[
                    static_cast<std::size_t>(r)];
                const T jb = ls_cell.parent_vertex_coords[
                    static_cast<std::size_t>(vb * gdim + r)]
                           - ls_cell.parent_vertex_coords[
                    static_cast<std::size_t>(r)];
                value += ja * jb;
            }
            gram[static_cast<std::size_t>(a * tdim + b)] = value;
        }
    }

    std::vector<T> rhs(grad_ref.begin(), grad_ref.end());
    std::vector<T> metric_direction;
    if (solve_dense_small<T>(std::move(gram), std::move(rhs), tdim, metric_direction))
        return metric_direction;
    return out;
}

template <std::floating_point T, std::integral I>
bool affine_jacobian(const LevelSetCell<T, I>& ls_cell,
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

template <std::floating_point T>
std::vector<T> pull_back_physical_direction(std::span<const T> jacobian,
                                            int gdim,
                                            int tdim,
                                            std::span<const T> direction_phys,
                                            std::span<const T> fallback_ref)
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
    if (solve_dense_small<T>(std::move(gram), std::move(rhs), tdim, out))
        return out;
    return std::vector<T>(fallback_ref.begin(), fallback_ref.end());
}

template <std::floating_point T>
std::vector<T> solve_physical_gradient(std::span<const T> jacobian,
                                       int gdim,
                                       int tdim,
                                       std::span<const T> grad_ref,
                                       std::span<const T> fallback_ref)
{
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
        std::vector<T> grad_phys;
        if (solve_dense_small<T>(std::move(jt), std::move(rhs), tdim, grad_phys))
            return grad_phys;
    }

    std::vector<T> out(static_cast<std::size_t>(gdim), T(0));
    for (int r = 0; r < gdim; ++r)
    {
        for (int a = 0; a < tdim; ++a)
        {
            out[static_cast<std::size_t>(r)] +=
                jacobian[static_cast<std::size_t>(r * tdim + a)]
              * fallback_ref[static_cast<std::size_t>(a)];
        }
    }
    return out;
}

template <std::floating_point T, std::integral I>
std::vector<T> level_set_physical_point(const LevelSetCell<T, I>& ls_cell,
                                        std::span<const T> ref)
{
    return cell::push_forward_affine_map<T>(
        ls_cell.cell_type,
        ls_cell.parent_vertex_coords,
        ls_cell.gdim,
        ref);
}

template <std::floating_point T, std::integral I>
std::vector<T> local_host_face_normal_physical(
    const LevelSetCell<T, I>& ls_cell,
    const AdaptCell<T>& ac,
    int local_zero_entity_id,
    T tol)
{
    const auto host = local_host_domain_for_zero_entity<T>(
        ac, local_zero_entity_id);
    if (ac.tdim != 3 || ls_cell.gdim != 3
        || !host.valid(ac.tdim) || host.dimension != 2)
    {
        return {};
    }

    const int nverts = static_cast<int>(
        host.vertices.size() / static_cast<std::size_t>(ac.tdim));
    if (nverts < 3)
        return {};

    const auto p0_ref = std::span<const T>(host.vertices.data(), 3);
    const auto p1_ref = std::span<const T>(host.vertices.data() + 3, 3);
    const auto p2_ref = std::span<const T>(host.vertices.data() + 6, 3);
    const auto p0 = level_set_physical_point<T, I>(ls_cell, p0_ref);
    const auto p1 = level_set_physical_point<T, I>(ls_cell, p1_ref);
    const auto p2 = level_set_physical_point<T, I>(ls_cell, p2_ref);

    std::array<T, 3> e0 = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    std::array<T, 3> e1 = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
    auto normal = local_cross<T>(
        std::span<const T>(e0.data(), 3),
        std::span<const T>(e1.data(), 3));
    const T nn = local_norm<T>(std::span<const T>(normal.data(), 3));
    if (nn <= tol)
        return {};

    std::vector<T> out(3, T(0));
    for (int d = 0; d < 3; ++d)
        out[static_cast<std::size_t>(d)] = normal[static_cast<std::size_t>(d)] / nn;
    return out;
}

template <std::floating_point T, std::integral I>
T curving_level_set_value(const LevelSetCell<T, I>& ls_cell,
                          std::span<const T> ref)
{
    return ls_cell.value(ref);
}

template <std::floating_point T, std::integral I>
void curving_level_set_grad(const LevelSetCell<T, I>& ls_cell,
                            std::span<const T> ref,
                            std::span<T> grad_ref)
{
    ls_cell.grad(ref, grad_ref);
}

template <std::floating_point T, std::integral I>
std::vector<T> physical_host_gradient_reference_direction(
    const LevelSetCell<T, I>& ls_cell,
    const AdaptCell<T>& ac,
    int local_zero_entity_id,
    std::span<const T> grad_ref,
    std::span<const T> fallback_ref)
{
    std::vector<T> jacobian;
    if (!affine_jacobian<T, I>(ls_cell, jacobian))
        return std::vector<T>(fallback_ref.begin(), fallback_ref.end());

    const int tdim = ls_cell.tdim;
    const int gdim = ls_cell.gdim;
    auto direction_phys = solve_physical_gradient<T>(
        std::span<const T>(jacobian.data(), jacobian.size()),
        gdim,
        tdim,
        grad_ref,
        fallback_ref);

    bool restricted_to_explicit_host_face = false;
    const auto explicit_host_normal =
        local_host_face_normal_physical<T, I>(
            ls_cell,
            ac,
            local_zero_entity_id,
            std::numeric_limits<T>::epsilon() * T(128));
    if (!explicit_host_normal.empty()
        && direction_phys.size() == explicit_host_normal.size())
    {
        const T nn = vec_dot<T>(
            std::span<const T>(
                explicit_host_normal.data(), explicit_host_normal.size()),
            std::span<const T>(
                explicit_host_normal.data(), explicit_host_normal.size()));
        if (nn > T(0))
        {
            const T c = vec_dot<T>(
                std::span<const T>(direction_phys.data(), direction_phys.size()),
                std::span<const T>(
                    explicit_host_normal.data(), explicit_host_normal.size())) / nn;
            for (std::size_t d = 0; d < direction_phys.size(); ++d)
                direction_phys[d] -= c * explicit_host_normal[d];
            restricted_to_explicit_host_face = true;
        }
    }

    const int parent_dim =
        ac.zero_entity_parent_dim[static_cast<std::size_t>(local_zero_entity_id)];
    const int parent_id =
        ac.zero_entity_parent_id[static_cast<std::size_t>(local_zero_entity_id)];
    if (!restricted_to_explicit_host_face && parent_dim == 2 && gdim == 3)
    {
        auto face = cell::face_vertices(ac.parent_cell_type, parent_id);
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
            if (nn > T(0))
            {
                const T c = (direction_phys[0] * normal[0]
                           + direction_phys[1] * normal[1]
                           + direction_phys[2] * normal[2]) / nn;
                for (int d = 0; d < 3; ++d)
                    direction_phys[static_cast<std::size_t>(d)] -= c * normal[d];
            }
        }
    }
    else if (!restricted_to_explicit_host_face && parent_dim == 1)
    {
        auto edges = cell::edges(ac.parent_cell_type);
        if (parent_id >= 0 && parent_id < static_cast<int>(edges.size()))
        {
            const auto edge = edges[static_cast<std::size_t>(parent_id)];
            std::vector<T> tangent(static_cast<std::size_t>(gdim), T(0));
            for (int d = 0; d < gdim; ++d)
            {
                tangent[static_cast<std::size_t>(d)] =
                    ls_cell.parent_vertex_coords[
                        static_cast<std::size_t>(edge[1] * gdim + d)]
                  - ls_cell.parent_vertex_coords[
                        static_cast<std::size_t>(edge[0] * gdim + d)];
            }
            const T tt = vec_dot<T>(
                std::span<const T>(tangent.data(), tangent.size()),
                std::span<const T>(tangent.data(), tangent.size()));
            if (tt > T(0))
            {
                const T c = vec_dot<T>(
                    std::span<const T>(direction_phys.data(), direction_phys.size()),
                    std::span<const T>(tangent.data(), tangent.size())) / tt;
                for (int d = 0; d < gdim; ++d)
                    direction_phys[static_cast<std::size_t>(d)] =
                        c * tangent[static_cast<std::size_t>(d)];
            }
        }
    }
    else if (!restricted_to_explicit_host_face && parent_dim == 0)
    {
        return std::vector<T>(static_cast<std::size_t>(tdim), T(0));
    }

    auto pulled = pull_back_physical_direction<T>(
        std::span<const T>(jacobian.data(), jacobian.size()),
        gdim,
        tdim,
        std::span<const T>(direction_phys.data(), direction_phys.size()),
        fallback_ref);
    return project_to_local_host_space<T>(
        ac,
        local_zero_entity_id,
        std::span<const T>(pulled.data(), pulled.size()),
        std::numeric_limits<T>::epsilon() * T(128));
}

template <std::floating_point T, std::integral I>
bool accept_seed_if_level_sets_within_tolerance(
    const AdaptCell<T>& ac,
    int local_zero_entity_id,
    std::span<const LevelSetCell<T, I>* const> active_cells,
    std::span<const T> seed,
    const CurvingOptions<T>& options,
    std::vector<T>& out,
    ProjectionStats<T>& stats)
{
    if (active_cells.empty())
        return false;

    stats.seed.assign(seed.begin(), seed.end());
    T max_abs_value = T(0);
    for (const LevelSetCell<T, I>* ls_cell : active_cells)
    {
        if (ls_cell == nullptr)
            return false;
        const T value = curving_level_set_value<T, I>(*ls_cell, seed);
        max_abs_value = std::max(max_abs_value, std::fabs(value));
        if (max_abs_value > options.ftol)
            return false;
    }

    out.assign(seed.begin(), seed.end());
    stats = {};
    stats.status = CurvingStatus::curved;
    stats.failure_code = CurvingFailureCode::none;
    stats.residual = max_abs_value;
    stats.active_face_mask = add_near_active_faces<T>(
        ac,
        seed,
        structural_active_face_mask<T>(ac, local_zero_entity_id),
        options.active_face_tol);
    stats.safe_subspace_dim = active_subspace_dim<T>(ac, stats.active_face_mask);
    stats.projection_mode = CurvingProjectionMode::none;
    return true;
}

template <std::floating_point T>
std::vector<std::vector<T>> scalar_candidate_directions(const AdaptCell<T>& ac,
                                                        int local_zero_entity_id,
                                                        std::span<const T> seed,
                                                        std::span<const T> grad,
                                                        std::span<const T> metric_grad,
                                                        std::span<const T> host_metric_grad,
                                                        std::uint32_t active_mask,
                                                        const CurvingOptions<T>& options)
{
    std::vector<std::vector<T>> directions;
    const T tol = std::max(options.ftol, T(64) * std::numeric_limits<T>::epsilon());
    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    const std::uint32_t host_mask =
        structural_active_face_mask<T>(ac, local_zero_entity_id);

    auto append_host_gradient = [&]()
    {
        // The primary implicit-geometry direction is restricted only by the
        // structural parent host entity. Extra near-face constraints can
        // overrestrict seeds on warped cut quads and create large tangential
        // motion, so those are kept for fallback directions below.
        append_admissible_unique_direction<T>(
            directions, ac, local_zero_entity_id, host_metric_grad, tol);
    };

    auto append_straight_normal = [&]() -> bool
    {
        const auto normal = straight_zero_entity_normal_direction<T>(
            ac, local_zero_entity_id, tol);
        if (!normal.degenerate())
        {
            auto direction = normal.value;
            orient_with_level_set_gradient<T>(direction, grad);
            append_admissible_unique_direction<T>(
                directions,
                ac,
                local_zero_entity_id,
                std::span<const T>(direction.data(), direction.size()),
                tol);
            return true;
        }
        return false;
    };

    if (options.direction_mode == CurvingDirectionMode::straight_zero_entity_normal)
    {
        append_straight_normal();
        append_host_gradient();
    }
    else
    {
        append_host_gradient();
    }

    if (active_mask != host_mask)
    {
        auto d = project_to_active_face_space<T>(ac, metric_grad, active_mask);
        append_admissible_unique_direction<T>(
            directions,
            ac,
            local_zero_entity_id,
            std::span<const T>(d.data(), d.size()),
            tol);
    }

    if (zdim == 1)
    {
        const auto tangent = zero_entity_edge_tangent<T>(ac, local_zero_entity_id);
        auto d = project_to_active_face_space<T>(ac, metric_grad, active_mask);
        d = remove_component<T>(
            std::move(d), std::span<const T>(tangent.data(), tangent.size()));
        append_admissible_unique_direction<T>(
            directions,
            ac,
            local_zero_entity_id,
            std::span<const T>(d.data(), d.size()),
            tol);
    }
    else if (zdim == 2)
    {
        const auto normal = straight_zero_entity_normal_direction<T>(
            ac, local_zero_entity_id, tol);
        if (!normal.degenerate())
        {
            auto d = project_to_active_face_space<T>(
                ac,
                std::span<const T>(normal.value.data(), normal.value.size()),
                active_mask);
            orient_with_level_set_gradient<T>(d, grad);
            append_admissible_unique_direction<T>(
                directions,
                ac,
                local_zero_entity_id,
                std::span<const T>(d.data(), d.size()),
                tol);
        }
    }

    {
        auto d = project_to_active_face_space<T>(ac, grad, active_mask);
        append_admissible_unique_direction<T>(
            directions,
            ac,
            local_zero_entity_id,
            std::span<const T>(d.data(), d.size()),
            tol);
    }

    const auto parent_vertices = cell::reference_vertices<T>(ac.parent_cell_type);
    const int nv = cell::get_num_vertices(ac.parent_cell_type);
    for (int v = 0; v < nv; ++v)
    {
        std::vector<T> ray(static_cast<std::size_t>(ac.tdim), T(0));
        for (int d = 0; d < ac.tdim; ++d)
            ray[static_cast<std::size_t>(d)] =
                seed[static_cast<std::size_t>(d)]
              - parent_vertices[static_cast<std::size_t>(v * ac.tdim + d)];
        auto d = project_to_active_face_space<T>(
            ac, std::span<const T>(ray.data(), ray.size()), active_mask);
        append_admissible_unique_direction<T>(
            directions,
            ac,
            local_zero_entity_id,
            std::span<const T>(d.data(), d.size()),
            tol);
    }

    const auto basis = active_nullspace_basis<T>(ac, active_mask);
    for (const auto& b : basis)
        append_admissible_unique_direction<T>(
            directions,
            ac,
            local_zero_entity_id,
            std::span<const T>(b.data(), b.size()),
            tol);

    if (options.direction_mode == CurvingDirectionMode::level_set_gradient)
        append_straight_normal();

    return directions;
}

template <std::floating_point T, class Eval>
bool try_scalar_line_search(const AdaptCell<T>& ac,
                            int local_zero_entity_id,
                            std::span<const T> seed,
                            std::span<const std::vector<T>> directions,
                            std::span<const T> gradient_at_seed,
                            const CurvingOptions<T>& options,
                            Eval&& eval_on_point,
                            std::vector<T>& out,
                            ProjectionStats<T>& stats)
{
    const int tdim = ac.tdim;
    for (const auto& unit : directions)
    {
        std::vector<T> direction = unit;
        if (direction.size() != static_cast<std::size_t>(tdim))
        {
            stats.failure_code = CurvingFailureCode::singular_gradient_system;
            continue;
        }
        if (gradient_at_seed.size() == direction.size()
            && vec_dot<T>(gradient_at_seed, std::span<const T>(
                                         direction.data(), direction.size())) < T(0))
        {
            for (T& value : direction)
                value = -value;
        }
        stats.seed.assign(seed.begin(), seed.end());
        stats.direction = direction;
        stats.clip_lo = std::numeric_limits<T>::quiet_NaN();
        stats.clip_hi = std::numeric_limits<T>::quiet_NaN();
        stats.root_t = std::numeric_limits<T>::quiet_NaN();

        T lo = T(0), hi = T(0);
        if (!host_parameter_interval<T>(
                ac,
                local_zero_entity_id,
                seed,
                std::span<const T>(direction.data(), direction.size()),
                options.domain_tol, lo, hi))
        {
            stats.failure_code = CurvingFailureCode::no_host_interval;
            continue;
        }
        stats.clip_lo = lo;
        stats.clip_hi = hi;
        if (!(lo <= T(0) && T(0) <= hi))
        {
            stats.failure_code = CurvingFailureCode::no_host_interval;
            continue;
        }

        auto eval_line = [&](T t) -> T
        {
            std::array<T, 3> x = {T(0), T(0), T(0)};
            for (int d = 0; d < tdim; ++d)
                x[static_cast<std::size_t>(d)] =
                    seed[static_cast<std::size_t>(d)]
                  + t * direction[static_cast<std::size_t>(d)];
            return eval_on_point(std::span<const T>(x.data(), static_cast<std::size_t>(tdim)));
        };

        const T f_seed = eval_line(T(0));
        if (std::fabs(f_seed) <= options.ftol)
        {
            out.assign(seed.begin(), seed.end());
            stats.iterations = 0;
            stats.status = CurvingStatus::curved;
            stats.failure_code = CurvingFailureCode::none;
            stats.residual = std::fabs(f_seed);
            stats.root_t = T(0);
            return true;
        }

        T a = T(0);
        T b = T(0);
        if (f_seed > T(0))
        {
            a = lo;
            b = T(0);
        }
        else
        {
            a = T(0);
            b = hi;
        }
        if (!(a < b))
        {
            stats.failure_code = CurvingFailureCode::no_sign_changing_bracket;
            continue;
        }

        const T fa = eval_line(a);
        const T fb = (b == T(0)) ? f_seed : eval_line(b);
        if (!(std::isfinite(fa) && std::isfinite(fb)) || fa * fb > T(0))
        {
            stats.failure_code = CurvingFailureCode::no_sign_changing_bracket;
            continue;
        }

        int iterations = 0;
        bool converged = false;
        const T linear_guess = cell::edge_root::linear_root_parameter<T>(fa, fb);
        const T initial_guess = a + (b - a) * std::clamp(linear_guess, T(0), T(1));
        const T root_t = cell::edge_root::newton_parameter<T>(
            eval_line, a, b, fa, fb,
            initial_guess,
            options.max_iter, options.xtol, options.ftol,
            &iterations, &converged);
        T accepted_root_t = root_t;
        stats.iterations = iterations;
        T residual = std::fabs(eval_line(root_t));
        stats.residual = residual;
        if (!converged && residual > options.ftol)
        {
            int fallback_iterations = 0;
            bool fallback_converged = false;
            try
            {
                const T fallback_root_t = cell::edge_root::itp_parameter<T>(
                    eval_line,
                    a,
                    b,
                    fa,
                    fb,
                    initial_guess,
                    options.max_iter,
                    options.xtol,
                    options.ftol,
                    &fallback_iterations,
                    &fallback_converged);
                const T fallback_residual = std::fabs(eval_line(fallback_root_t));
                if (fallback_converged || fallback_residual <= options.ftol)
                {
                    accepted_root_t = fallback_root_t;
                    iterations += fallback_iterations;
                    residual = fallback_residual;
                    stats.iterations = iterations;
                    stats.residual = residual;
                    converged = true;
                }
            }
            catch (const std::exception&)
            {
            }
            if (!converged && residual > options.ftol)
            {
                stats.failure_code = CurvingFailureCode::brent_failed;
                stats.root_t = accepted_root_t;
                continue;
            }
        }

        out.assign(seed.begin(), seed.end());
        for (int d = 0; d < tdim; ++d)
            out[static_cast<std::size_t>(d)] += accepted_root_t * direction[static_cast<std::size_t>(d)];
        const bool inside = cell::edge_root::is_inside_reference_domain<T>(
            std::span<const T>(out.data(), out.size()),
            ac.parent_cell_type,
            options.domain_tol);
        if (!inside)
        {
            stats.status = CurvingStatus::failed;
            stats.failure_code = CurvingFailureCode::outside_host_domain;
            stats.root_t = accepted_root_t;
            return false;
        }
        stats.status = CurvingStatus::curved;
        stats.failure_code = CurvingFailureCode::none;
        stats.direction = direction;
        stats.root_t = accepted_root_t;
        stats.residual = residual;
        return true;
    }

    stats.status = CurvingStatus::failed;
    return false;
}

template <std::floating_point T, class Eval, class Grad>
bool constrained_scalar_newton(const AdaptCell<T>& ac,
                               int local_zero_entity_id,
                               std::span<const T> seed,
                               std::uint32_t active_mask,
                               const CurvingOptions<T>& options,
                               Eval&& eval_on_point,
                               Grad&& grad_on_point,
                               std::vector<T>& out,
                               ProjectionStats<T>& stats)
{
    const int tdim = ac.tdim;
    out.assign(seed.begin(), seed.end());
    stats.projection_mode = CurvingProjectionMode::constrained_newton;
    stats.active_face_mask = active_mask;
    stats.safe_subspace_dim = active_subspace_dim<T>(ac, active_mask);

    if (stats.safe_subspace_dim <= 0)
    {
        stats.failure_code = CurvingFailureCode::constrained_newton_failed;
        return false;
    }

    std::vector<T> grad(static_cast<std::size_t>(tdim), T(0));
    for (int iter = 0; iter < options.max_iter; ++iter)
    {
        const T f = eval_on_point(std::span<const T>(out.data(), out.size()));
        if (std::fabs(f) <= options.ftol)
        {
            stats.iterations = iter;
            stats.status = CurvingStatus::curved;
            stats.failure_code = CurvingFailureCode::none;
            stats.residual = std::fabs(f);
            stats.direction = normalized_delta_direction<T>(
                seed,
                std::span<const T>(out.data(), out.size()),
                options.xtol);
            return true;
        }

        grad_on_point(std::span<const T>(out.data(), out.size()),
                      std::span<T>(grad.data(), grad.size()));
        auto pg = project_to_active_face_space<T>(
            ac, std::span<const T>(grad.data(), grad.size()), active_mask);
        const T gg = vec_dot<T>(
            std::span<const T>(pg.data(), pg.size()),
            std::span<const T>(pg.data(), pg.size()));
        if (gg <= options.ftol * options.ftol)
        {
            stats.iterations = iter;
            stats.failure_code = CurvingFailureCode::singular_gradient_system;
            stats.residual = std::fabs(f);
            return false;
        }

        std::vector<T> delta(static_cast<std::size_t>(tdim), T(0));
        for (int d = 0; d < tdim; ++d)
            delta[static_cast<std::size_t>(d)] =
                -f * pg[static_cast<std::size_t>(d)] / gg;

        T lo = T(0), hi = T(0);
        if (!host_parameter_interval<T>(
                ac,
                local_zero_entity_id,
                std::span<const T>(out.data(), out.size()),
                std::span<const T>(delta.data(), delta.size()),
                options.domain_tol,
                lo,
                hi))
        {
            stats.failure_code = CurvingFailureCode::no_host_interval;
            stats.residual = std::fabs(f);
            return false;
        }

        T step = std::min(T(1), hi);
        if (!(step > T(0)))
        {
            stats.failure_code = CurvingFailureCode::line_search_failed;
            stats.residual = std::fabs(f);
            return false;
        }

        const T merit = T(0.5) * f * f;
        bool accepted = false;
        std::vector<T> candidate(static_cast<std::size_t>(tdim), T(0));
        for (int ls = 0; ls < 16; ++ls)
        {
            for (int d = 0; d < tdim; ++d)
                candidate[static_cast<std::size_t>(d)] =
                    out[static_cast<std::size_t>(d)]
                  + step * delta[static_cast<std::size_t>(d)];
            if (!cell::edge_root::is_inside_reference_domain<T>(
                    std::span<const T>(candidate.data(), candidate.size()),
                    ac.parent_cell_type,
                    options.domain_tol))
            {
                step *= T(0.5);
                continue;
            }
            const T cand_f = eval_on_point(
                std::span<const T>(candidate.data(), candidate.size()));
            const T cand_merit = T(0.5) * cand_f * cand_f;
            if (cand_merit < merit)
            {
                out = candidate;
                accepted = true;
                break;
            }
            step *= T(0.5);
        }
        if (!accepted)
        {
            stats.iterations = iter + 1;
            stats.failure_code = CurvingFailureCode::line_search_failed;
            stats.residual = std::fabs(f);
            return false;
        }
    }

    const T f = eval_on_point(std::span<const T>(out.data(), out.size()));
    stats.iterations = options.max_iter;
    stats.residual = std::fabs(f);
    if (stats.residual <= options.ftol)
    {
        stats.status = CurvingStatus::curved;
        stats.failure_code = CurvingFailureCode::none;
        stats.direction = normalized_delta_direction<T>(
            seed,
            std::span<const T>(out.data(), out.size()),
            options.xtol);
        return true;
    }
    stats.failure_code = CurvingFailureCode::constrained_newton_failed;
    return false;
}

template <std::floating_point T>
bool solve_dense_small(std::vector<T> A, std::vector<T> b, int n, std::vector<T>& x)
{
    x.assign(static_cast<std::size_t>(n), T(0));
    if (n <= 0 || n > 3)
        return false;

    for (int k = 0; k < n; ++k)
    {
        int pivot = k;
        T best = std::fabs(A[static_cast<std::size_t>(k * n + k)]);
        for (int r = k + 1; r < n; ++r)
        {
            const T val = std::fabs(A[static_cast<std::size_t>(r * n + k)]);
            if (val > best)
            {
                best = val;
                pivot = r;
            }
        }
        if (best <= T(64) * std::numeric_limits<T>::epsilon())
            return false;
        if (pivot != k)
        {
            for (int c = k; c < n; ++c)
                std::swap(A[static_cast<std::size_t>(k * n + c)],
                          A[static_cast<std::size_t>(pivot * n + c)]);
            std::swap(b[static_cast<std::size_t>(k)], b[static_cast<std::size_t>(pivot)]);
        }

        const T diag = A[static_cast<std::size_t>(k * n + k)];
        for (int c = k; c < n; ++c)
            A[static_cast<std::size_t>(k * n + c)] /= diag;
        b[static_cast<std::size_t>(k)] /= diag;

        for (int r = 0; r < n; ++r)
        {
            if (r == k)
                continue;
            const T factor = A[static_cast<std::size_t>(r * n + k)];
            if (factor == T(0))
                continue;
            for (int c = k; c < n; ++c)
                A[static_cast<std::size_t>(r * n + c)] -= factor * A[static_cast<std::size_t>(k * n + c)];
            b[static_cast<std::size_t>(r)] -= factor * b[static_cast<std::size_t>(k)];
        }
    }

    x = std::move(b);
    return true;
}

template <std::floating_point T>
std::vector<T> projected_direction_to_host(const AdaptCell<T>& ac,
                                           int local_zero_entity_id,
                                           std::span<const T> direction)
{
    const int tdim = ac.tdim;
    std::vector<T> out(direction.begin(), direction.end());
    const int parent_dim = ac.zero_entity_parent_dim[static_cast<std::size_t>(local_zero_entity_id)];
    const int parent_id = ac.zero_entity_parent_id[static_cast<std::size_t>(local_zero_entity_id)];
    if (parent_dim < 0 || parent_dim >= tdim)
        return out;

    std::vector<T> ref_vertices = cell::reference_vertices<T>(ac.parent_cell_type);
    std::vector<std::array<T, 3>> basis;

    auto add_basis_from_vertices = [&](int v0, int v1)
    {
        std::array<T, 3> b = {T(0), T(0), T(0)};
        for (int d = 0; d < tdim; ++d)
            b[static_cast<std::size_t>(d)] =
                ref_vertices[static_cast<std::size_t>(v1 * tdim + d)]
              - ref_vertices[static_cast<std::size_t>(v0 * tdim + d)];
        basis.push_back(b);
    };

    if (parent_dim == 1)
    {
        auto edges = cell::edges(ac.parent_cell_type);
        if (parent_id >= 0 && parent_id < static_cast<int>(edges.size()))
            add_basis_from_vertices(edges[static_cast<std::size_t>(parent_id)][0],
                                    edges[static_cast<std::size_t>(parent_id)][1]);
    }
    else if (parent_dim == 2 && ac.tdim == 3)
    {
        auto fv = cell::face_vertices(ac.parent_cell_type, parent_id);
        add_basis_from_vertices(fv[0], fv[1]);
        add_basis_from_vertices(fv[0], fv[2]);
    }

    if (basis.empty())
        return out;

    std::fill(out.begin(), out.end(), T(0));
    if (basis.size() == 1)
    {
        T bd = T(0);
        T bb = T(0);
        for (int d = 0; d < tdim; ++d)
        {
            bd += basis[0][static_cast<std::size_t>(d)] * direction[static_cast<std::size_t>(d)];
            bb += basis[0][static_cast<std::size_t>(d)] * basis[0][static_cast<std::size_t>(d)];
        }
        if (bb > T(0))
            for (int d = 0; d < tdim; ++d)
                out[static_cast<std::size_t>(d)] = (bd / bb) * basis[0][static_cast<std::size_t>(d)];
        return out;
    }

    T g00 = T(0), g01 = T(0), g11 = T(0), r0 = T(0), r1 = T(0);
    for (int d = 0; d < tdim; ++d)
    {
        g00 += basis[0][static_cast<std::size_t>(d)] * basis[0][static_cast<std::size_t>(d)];
        g01 += basis[0][static_cast<std::size_t>(d)] * basis[1][static_cast<std::size_t>(d)];
        g11 += basis[1][static_cast<std::size_t>(d)] * basis[1][static_cast<std::size_t>(d)];
        r0 += basis[0][static_cast<std::size_t>(d)] * direction[static_cast<std::size_t>(d)];
        r1 += basis[1][static_cast<std::size_t>(d)] * direction[static_cast<std::size_t>(d)];
    }
    const T det = g00 * g11 - g01 * g01;
    if (std::fabs(det) <= T(64) * std::numeric_limits<T>::epsilon())
        return out;
    const T a = ( r0 * g11 - r1 * g01) / det;
    const T b = (-r0 * g01 + r1 * g00) / det;
    for (int d = 0; d < tdim; ++d)
        out[static_cast<std::size_t>(d)] =
            a * basis[0][static_cast<std::size_t>(d)]
          + b * basis[1][static_cast<std::size_t>(d)];
    return out;
}

template <std::floating_point T>
std::vector<int> active_level_sets(std::uint64_t mask)
{
    std::vector<int> ids;
    for (int i = 0; i < 64; ++i)
        if ((mask & (std::uint64_t(1) << i)) != 0)
            ids.push_back(i);
    return ids;
}

template <std::floating_point T, std::integral I>
const LevelSetCell<T, I>& find_level_set_cell(std::span<const LevelSetCell<T, I>> cells,
                                              std::span<const int> offsets,
                                              int cut_cell_id,
                                              int level_set_id)
{
    for (int i = offsets[static_cast<std::size_t>(cut_cell_id)];
         i < offsets[static_cast<std::size_t>(cut_cell_id + 1)]; ++i)
    {
        if (cells[static_cast<std::size_t>(i)].level_set_id == level_set_id)
            return cells[static_cast<std::size_t>(i)];
    }
    throw std::runtime_error("curving: missing LevelSetCell for zero entity");
}

template <std::floating_point T, std::integral I>
bool scalar_project(const AdaptCell<T>& ac,
                    int local_zero_entity_id,
                    const LevelSetCell<T, I>& ls_cell,
                    std::span<const T> seed,
                    const CurvingOptions<T>& options,
                    std::vector<T>& out,
                    ProjectionStats<T>& stats)
{
    const int tdim = ac.tdim;
    stats = {};
    stats.failure_code = CurvingFailureCode::no_sign_changing_bracket;
    const T seed_value = curving_level_set_value<T, I>(ls_cell, seed);
    if (std::fabs(seed_value) <= options.ftol)
    {
        out.assign(seed.begin(), seed.end());
        stats.seed.assign(seed.begin(), seed.end());
        stats.status = CurvingStatus::curved;
        stats.failure_code = CurvingFailureCode::none;
        stats.residual = std::fabs(seed_value);
        stats.active_face_mask = add_near_active_faces<T>(
            ac,
            seed,
            structural_active_face_mask<T>(ac, local_zero_entity_id),
            options.active_face_tol);
        stats.safe_subspace_dim = active_subspace_dim<T>(ac, stats.active_face_mask);
        stats.projection_mode = CurvingProjectionMode::none;
        return true;
    }
    std::vector<T> grad(static_cast<std::size_t>(tdim), T(0));
    curving_level_set_grad<T, I>(
        ls_cell, seed, std::span<T>(grad.data(), grad.size()));
    const auto metric_grad = physical_normal_reference_direction<T, I>(
        ls_cell, std::span<const T>(grad.data(), grad.size()));
    const auto host_metric_grad = physical_host_gradient_reference_direction<T, I>(
        ls_cell,
        ac,
        local_zero_entity_id,
        std::span<const T>(grad.data(), grad.size()),
        std::span<const T>(metric_grad.data(), metric_grad.size()));

    auto eval_on_point = [&](std::span<const T> x) -> T
    {
        return curving_level_set_value<T, I>(ls_cell, x);
    };
    auto grad_on_point = [&](std::span<const T> x, std::span<T> out_grad)
    {
        curving_level_set_grad<T, I>(ls_cell, x, out_grad);
    };

    std::uint32_t active_mask = structural_active_face_mask<T>(ac, local_zero_entity_id);
    active_mask = add_near_active_faces<T>(ac, seed, active_mask, options.active_face_tol);

    stats.active_face_mask = active_mask;
    stats.safe_subspace_dim = active_subspace_dim<T>(ac, active_mask);
    stats.projection_mode = CurvingProjectionMode::safe_line;

    auto directions = scalar_candidate_directions<T>(
        ac,
        local_zero_entity_id,
        seed,
        std::span<const T>(grad.data(), grad.size()),
        std::span<const T>(metric_grad.data(), metric_grad.size()),
        std::span<const T>(host_metric_grad.data(), host_metric_grad.size()),
        active_mask,
        options);
    if (try_scalar_line_search<T>(
            ac,
            local_zero_entity_id,
            seed,
            std::span<const std::vector<T>>(directions.data(), directions.size()),
            std::span<const T>(grad.data(), grad.size()),
            options,
            eval_on_point,
            out,
            stats))
    {
        stats.active_face_mask = active_mask;
        stats.safe_subspace_dim = active_subspace_dim<T>(ac, active_mask);
        stats.projection_mode = CurvingProjectionMode::safe_line;
        return true;
    }

    int closest_face_id = -1;
    const std::uint32_t retry_mask =
        add_closest_faces<T>(ac, seed, active_mask, closest_face_id);
    if (retry_mask != active_mask)
    {
        ProjectionStats<T> retry_stats = stats;
        retry_stats.active_face_mask = retry_mask;
        retry_stats.closest_face_id = closest_face_id;
        retry_stats.safe_subspace_dim = active_subspace_dim<T>(ac, retry_mask);
        retry_stats.projection_mode = CurvingProjectionMode::closest_face_retry;
        retry_stats.retry_count = 1;

        directions = scalar_candidate_directions<T>(
            ac,
            local_zero_entity_id,
            seed,
            std::span<const T>(grad.data(), grad.size()),
            std::span<const T>(metric_grad.data(), metric_grad.size()),
            std::span<const T>(host_metric_grad.data(), host_metric_grad.size()),
            retry_mask,
            options);
        if (try_scalar_line_search<T>(
                ac,
                local_zero_entity_id,
                seed,
                std::span<const std::vector<T>>(directions.data(), directions.size()),
                std::span<const T>(grad.data(), grad.size()),
                options,
                eval_on_point,
                out,
                retry_stats))
        {
            retry_stats.active_face_mask = retry_mask;
            retry_stats.closest_face_id = closest_face_id;
            retry_stats.safe_subspace_dim = active_subspace_dim<T>(ac, retry_mask);
            retry_stats.projection_mode = CurvingProjectionMode::closest_face_retry;
            retry_stats.retry_count = 1;
            stats = retry_stats;
            return true;
        }
        stats = retry_stats;
        stats.failure_code = CurvingFailureCode::closest_face_retry_failed;
    }

    std::array<std::uint32_t, 2> newton_masks = {
        (retry_mask != active_mask) ? retry_mask : active_mask,
        active_mask
    };
    const int newton_attempts = (newton_masks[0] == newton_masks[1]) ? 1 : 2;
    ProjectionStats<T> best_newton_stats = stats;
    for (int attempt = 0; attempt < newton_attempts; ++attempt)
    {
        ProjectionStats<T> newton_stats = stats;
        const std::uint32_t newton_mask = newton_masks[static_cast<std::size_t>(attempt)];
        newton_stats.active_face_mask = newton_mask;
        newton_stats.closest_face_id = (newton_mask == retry_mask) ? closest_face_id : -1;
        newton_stats.retry_count = (retry_mask != active_mask) ? 2 + attempt : 1 + attempt;
        if (constrained_scalar_newton<T>(
                ac,
                local_zero_entity_id,
                seed,
                newton_mask,
                options,
                eval_on_point,
                grad_on_point,
                out,
                newton_stats))
        {
            stats = newton_stats;
            return true;
        }
        best_newton_stats = newton_stats;
    }

    stats = best_newton_stats;
    stats.status = CurvingStatus::failed;
    if (stats.failure_code == CurvingFailureCode::none)
        stats.failure_code = CurvingFailureCode::constrained_newton_failed;
    return false;
}

template <std::floating_point T, std::integral I>
bool fixed_ray_scalar_project(const AdaptCell<T>& ac,
                              int local_zero_entity_id,
                              const LevelSetCell<T, I>& ls_cell,
                              std::span<const T> seed,
                              const CurvingOptions<T>& options,
                              std::vector<T>& out,
                              ProjectionStats<T>& stats)
{
    const int tdim = ac.tdim;
    stats = {};
    const std::uint32_t host_mask =
        structural_active_face_mask<T>(ac, local_zero_entity_id);
    stats.active_face_mask = host_mask;
    stats.safe_subspace_dim = active_subspace_dim<T>(ac, host_mask);
    stats.projection_mode = CurvingProjectionMode::safe_line;
    stats.failure_code = CurvingFailureCode::no_sign_changing_bracket;

    const T seed_value = curving_level_set_value<T, I>(ls_cell, seed);
    if (std::fabs(seed_value) <= options.ftol)
    {
        out.assign(seed.begin(), seed.end());
        stats.seed.assign(seed.begin(), seed.end());
        stats.status = CurvingStatus::curved;
        stats.failure_code = CurvingFailureCode::none;
        stats.residual = std::fabs(seed_value);
        return true;
    }

    std::vector<T> grad(static_cast<std::size_t>(tdim), T(0));
    curving_level_set_grad<T, I>(
        ls_cell, seed, std::span<T>(grad.data(), grad.size()));
    const auto metric_grad = physical_normal_reference_direction<T, I>(
        ls_cell, std::span<const T>(grad.data(), grad.size()));
    auto direction = physical_host_gradient_reference_direction<T, I>(
        ls_cell,
        ac,
        local_zero_entity_id,
        std::span<const T>(grad.data(), grad.size()),
        std::span<const T>(metric_grad.data(), metric_grad.size()));
    direction = project_to_active_face_space<T>(
        ac, std::span<const T>(direction.data(), direction.size()), host_mask);

    std::vector<std::vector<T>> directions;
    const T direction_tol =
        std::max(options.ftol, T(64) * std::numeric_limits<T>::epsilon());

    auto append_gradient_direction = [&]()
    {
        append_admissible_unique_direction<T>(
            directions,
            ac,
            local_zero_entity_id,
            std::span<const T>(direction.data(), direction.size()),
            direction_tol);
    };
    auto append_normal_direction = [&]() -> bool
    {
        const auto normal = straight_zero_entity_normal_direction<T>(
            ac, local_zero_entity_id, direction_tol);
        if (!normal.degenerate())
        {
            auto oriented = normal.value;
            orient_with_level_set_gradient<T>(oriented, std::span<const T>(grad.data(), grad.size()));
            append_admissible_unique_direction<T>(
                directions,
                ac,
                local_zero_entity_id,
                std::span<const T>(oriented.data(), oriented.size()),
                direction_tol);
            return true;
        }
        return false;
    };
    auto append_gradient_fallback_directions = [&]()
    {
        {
            auto d = project_to_active_face_space<T>(
                ac, std::span<const T>(grad.data(), grad.size()), host_mask);
            append_admissible_unique_direction<T>(
                directions,
                ac,
                local_zero_entity_id,
                std::span<const T>(d.data(), d.size()),
                direction_tol);
        }
        {
            auto d = project_to_active_face_space<T>(
                ac, std::span<const T>(metric_grad.data(), metric_grad.size()), host_mask);
            append_admissible_unique_direction<T>(
                directions,
                ac,
                local_zero_entity_id,
                std::span<const T>(d.data(), d.size()),
                direction_tol);
        }

        const auto parent_vertices = cell::reference_vertices<T>(ac.parent_cell_type);
        const int nv = cell::get_num_vertices(ac.parent_cell_type);
        for (int v = 0; v < nv; ++v)
        {
            std::vector<T> ray(static_cast<std::size_t>(tdim), T(0));
            for (int d = 0; d < tdim; ++d)
                ray[static_cast<std::size_t>(d)] =
                    seed[static_cast<std::size_t>(d)]
                  - parent_vertices[static_cast<std::size_t>(v * tdim + d)];
            auto d = project_to_active_face_space<T>(
                ac, std::span<const T>(ray.data(), ray.size()), host_mask);
            append_admissible_unique_direction<T>(
                directions,
                ac,
                local_zero_entity_id,
                std::span<const T>(d.data(), d.size()),
                direction_tol);
        }

        const auto basis = active_nullspace_basis<T>(ac, host_mask);
        for (const auto& b : basis)
            append_admissible_unique_direction<T>(
                directions,
                ac,
                local_zero_entity_id,
                std::span<const T>(b.data(), b.size()),
                direction_tol);
    };

    if (options.direction_mode == CurvingDirectionMode::straight_zero_entity_normal)
    {
        append_normal_direction();
        append_gradient_direction();
        append_gradient_fallback_directions();
    }
    else
    {
        append_gradient_direction();
        append_gradient_fallback_directions();
        append_normal_direction();
    }

    if (directions.empty())
    {
        stats.status = CurvingStatus::failed;
        stats.failure_code = CurvingFailureCode::singular_gradient_system;
        return false;
    }

    auto eval_on_point = [&](std::span<const T> x) -> T
    {
        return curving_level_set_value<T, I>(ls_cell, x);
    };

    ProjectionStats<T> best_stats = stats;
    for (const auto& unit : directions)
    {
        std::vector<T> candidate;
        ProjectionStats<T> candidate_stats = stats;
        if (!try_scalar_line_search<T>(
                ac,
                local_zero_entity_id,
                seed,
                std::span<const std::vector<T>>(&unit, 1),
                std::span<const T>(grad.data(), grad.size()),
                options,
                eval_on_point,
                candidate,
                candidate_stats))
        {
            best_stats = candidate_stats;
            continue;
        }

        out = std::move(candidate);
        stats = candidate_stats;
        stats.active_face_mask = host_mask;
        stats.safe_subspace_dim = active_subspace_dim<T>(ac, host_mask);
        stats.projection_mode = CurvingProjectionMode::safe_line;
        return true;
    }

    stats = best_stats;
    stats.active_face_mask = host_mask;
    stats.safe_subspace_dim = active_subspace_dim<T>(ac, host_mask);
    stats.projection_mode = CurvingProjectionMode::safe_line;
    stats.status = CurvingStatus::failed;
    if (stats.failure_code == CurvingFailureCode::none)
        stats.failure_code = CurvingFailureCode::projection_failed;
    return false;
}

template <std::floating_point T, std::integral I>
bool vector_project(const AdaptCell<T>& ac,
                    int local_zero_entity_id,
                    std::span<const LevelSetCell<T, I>* const> ls_cells,
                    std::span<const T> seed,
                    const CurvingOptions<T>& options,
                    std::vector<T>& out,
                    ProjectionStats<T>& stats)
{
    const int tdim = ac.tdim;
    const int m = static_cast<int>(ls_cells.size());
    stats = {};
    stats.failure_code = CurvingFailureCode::projection_failed;
    out.assign(seed.begin(), seed.end());
    if (m <= 0 || m > tdim)
    {
        stats.failure_code = CurvingFailureCode::invalid_constraint_count;
        return false;
    }

    if (accept_seed_if_level_sets_within_tolerance<T, I>(
            ac, local_zero_entity_id, ls_cells, seed, options, out, stats))
        return true;

    std::vector<T> values(static_cast<std::size_t>(m), T(0));
    std::vector<T> grads(static_cast<std::size_t>(m * tdim), T(0));
    std::uint32_t active_mask = structural_active_face_mask<T>(ac, local_zero_entity_id);
    active_mask = add_near_active_faces<T>(ac, seed, active_mask, options.active_face_tol);
    int closest_face_id = -1;
    int retry_count = 0;

    for (int iter = 0; iter < options.max_iter; ++iter)
    {
        T norm_f = T(0);
        T max_f = T(0);
        for (int i = 0; i < m; ++i)
        {
            values[static_cast<std::size_t>(i)] =
                ls_cells[static_cast<std::size_t>(i)]->value(
                    std::span<const T>(out.data(), out.size()));
            norm_f += values[static_cast<std::size_t>(i)] * values[static_cast<std::size_t>(i)];
            max_f = std::max(max_f, std::fabs(values[static_cast<std::size_t>(i)]));
            ls_cells[static_cast<std::size_t>(i)]->grad(
                std::span<const T>(out.data(), out.size()),
                std::span<T>(grads.data() + static_cast<std::size_t>(i * tdim),
                             static_cast<std::size_t>(tdim)));
        }
        norm_f = std::sqrt(norm_f);
        if (max_f <= options.ftol)
        {
            stats.iterations = iter;
            stats.status = CurvingStatus::curved;
            stats.failure_code = CurvingFailureCode::none;
            stats.residual = max_f;
            stats.active_face_mask = active_mask;
            stats.closest_face_id = closest_face_id;
            stats.safe_subspace_dim = active_subspace_dim<T>(ac, active_mask);
            stats.projection_mode =
                (iter == 0) ? CurvingProjectionMode::none : CurvingProjectionMode::vector_newton;
            stats.retry_count = retry_count;
            stats.direction = normalized_delta_direction<T>(
                seed,
                std::span<const T>(out.data(), out.size()),
                options.xtol);
            return true;
        }

        const auto basis = active_nullspace_basis<T>(ac, active_mask);
        const int q = static_cast<int>(basis.size());
        if (q < m)
        {
            const std::uint32_t retry_mask =
                add_closest_faces<T>(
                    ac,
                    std::span<const T>(out.data(), out.size()),
                    active_mask,
                    closest_face_id);
            if (retry_mask != active_mask)
            {
                active_mask = retry_mask;
                ++retry_count;
                continue;
            }
            stats.iterations = iter;
            stats.failure_code = CurvingFailureCode::singular_gradient_system;
            stats.residual = norm_f;
            stats.active_face_mask = active_mask;
            stats.closest_face_id = closest_face_id;
            stats.safe_subspace_dim = q;
            stats.projection_mode = CurvingProjectionMode::vector_newton;
            stats.retry_count = retry_count;
            return false;
        }

        std::vector<T> jq(static_cast<std::size_t>(m * q), T(0));
        for (int i = 0; i < m; ++i)
        {
            std::span<const T> gi(
                grads.data() + static_cast<std::size_t>(i * tdim),
                static_cast<std::size_t>(tdim));
            for (int k = 0; k < q; ++k)
                jq[static_cast<std::size_t>(i * q + k)] =
                    vec_dot<T>(gi, std::span<const T>(basis[static_cast<std::size_t>(k)].data(),
                                                      basis[static_cast<std::size_t>(k)].size()));
        }

        std::vector<T> normal_matrix(static_cast<std::size_t>(m * m), T(0));
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                T a = T(0);
                for (int k = 0; k < q; ++k)
                    a += jq[static_cast<std::size_t>(i * q + k)]
                       * jq[static_cast<std::size_t>(j * q + k)];
                normal_matrix[static_cast<std::size_t>(i * m + j)] = a;
            }
        }
        std::vector<T> rhs(static_cast<std::size_t>(m), T(0));
        for (int i = 0; i < m; ++i)
            rhs[static_cast<std::size_t>(i)] = -values[static_cast<std::size_t>(i)];

        std::vector<T> lambda;
        if (!solve_dense_small<T>(normal_matrix, rhs, m, lambda))
        {
            const std::uint32_t retry_mask =
                add_closest_faces<T>(
                    ac,
                    std::span<const T>(out.data(), out.size()),
                    active_mask,
                    closest_face_id);
            if (retry_mask != active_mask)
            {
                active_mask = retry_mask;
                ++retry_count;
                continue;
            }
            stats.iterations = iter;
            stats.failure_code = CurvingFailureCode::singular_gradient_system;
            stats.residual = norm_f;
            stats.active_face_mask = active_mask;
            stats.closest_face_id = closest_face_id;
            stats.safe_subspace_dim = q;
            stats.projection_mode = CurvingProjectionMode::vector_newton;
            stats.retry_count = retry_count;
            return false;
        }

        std::vector<T> delta(static_cast<std::size_t>(tdim), T(0));
        for (int k = 0; k < q; ++k)
        {
            T dy = T(0);
            for (int i = 0; i < m; ++i)
                dy += jq[static_cast<std::size_t>(i * q + k)]
                    * lambda[static_cast<std::size_t>(i)];
            for (int d = 0; d < tdim; ++d)
                delta[static_cast<std::size_t>(d)] +=
                    dy * basis[static_cast<std::size_t>(k)][static_cast<std::size_t>(d)];
        }

        if (vec_norm<T>(std::span<const T>(delta.data(), delta.size())) <= options.xtol)
        {
            stats.iterations = iter;
            stats.failure_code = CurvingFailureCode::singular_gradient_system;
            stats.residual = norm_f;
            stats.active_face_mask = active_mask;
            stats.closest_face_id = closest_face_id;
            stats.safe_subspace_dim = q;
            stats.projection_mode = CurvingProjectionMode::vector_newton;
            stats.retry_count = retry_count;
            return false;
        }

        T lo = T(0), hi = T(0);
        if (!host_parameter_interval<T>(
                ac,
                local_zero_entity_id,
                std::span<const T>(out.data(), out.size()),
                std::span<const T>(delta.data(), delta.size()),
                options.domain_tol,
                lo,
                hi))
        {
            stats.iterations = iter;
            stats.failure_code = CurvingFailureCode::no_host_interval;
            stats.residual = norm_f;
            stats.active_face_mask = active_mask;
            stats.closest_face_id = closest_face_id;
            stats.safe_subspace_dim = q;
            stats.projection_mode = CurvingProjectionMode::vector_newton;
            stats.retry_count = retry_count;
            return false;
        }

        T step = T(1);
        if (hi > T(0) && hi < step)
            step = hi;
        bool accepted = false;
        std::vector<T> candidate(static_cast<std::size_t>(tdim), T(0));
        const T merit = T(0.5) * norm_f * norm_f;
        for (int ls = 0; ls < 12; ++ls)
        {
            for (int d = 0; d < tdim; ++d)
                candidate[static_cast<std::size_t>(d)] =
                    out[static_cast<std::size_t>(d)] + step * delta[static_cast<std::size_t>(d)];
            if (!cell::edge_root::is_inside_reference_domain<T>(
                    std::span<const T>(candidate.data(), candidate.size()),
                    ac.parent_cell_type,
                    options.domain_tol))
            {
                step *= T(0.5);
                continue;
            }
            T cand_norm = T(0);
            for (int i = 0; i < m; ++i)
            {
                const T f = ls_cells[static_cast<std::size_t>(i)]->value(
                    std::span<const T>(candidate.data(), candidate.size()));
                cand_norm += f * f;
            }
            const T cand_merit = T(0.5) * cand_norm;
            if (cand_merit < merit)
            {
                out = candidate;
                accepted = true;
                break;
            }
            step *= T(0.5);
        }
        if (!accepted)
        {
            const std::uint32_t retry_mask =
                add_closest_faces<T>(
                    ac,
                    std::span<const T>(out.data(), out.size()),
                    active_mask,
                    closest_face_id);
            if (retry_mask != active_mask)
            {
                active_mask = retry_mask;
                ++retry_count;
                continue;
            }
            stats.iterations = iter + 1;
            stats.failure_code = CurvingFailureCode::line_search_failed;
            stats.residual = norm_f;
            stats.active_face_mask = active_mask;
            stats.closest_face_id = closest_face_id;
            stats.safe_subspace_dim = q;
            stats.projection_mode = CurvingProjectionMode::vector_newton;
            stats.retry_count = retry_count;
            return false;
        }
    }

    T final_max = T(0);
    for (int i = 0; i < m; ++i)
    {
        const T f = ls_cells[static_cast<std::size_t>(i)]->value(
            std::span<const T>(out.data(), out.size()));
        final_max = std::max(final_max, std::fabs(f));
    }
    stats.iterations = options.max_iter;
    stats.residual = final_max;
    stats.active_face_mask = active_mask;
    stats.closest_face_id = closest_face_id;
    stats.safe_subspace_dim = active_subspace_dim<T>(ac, active_mask);
    stats.projection_mode = CurvingProjectionMode::vector_newton;
    stats.retry_count = retry_count;
    if (final_max <= options.ftol)
    {
        stats.status = CurvingStatus::curved;
        stats.failure_code = CurvingFailureCode::none;
        stats.direction = normalized_delta_direction<T>(
            seed,
            std::span<const T>(out.data(), out.size()),
            options.xtol);
        return true;
    }
    stats.status = CurvingStatus::failed;
    stats.failure_code = CurvingFailureCode::max_iterations;
    return false;
}

template <std::floating_point T, std::integral I>
bool project_seed_to_zero_entity(const AdaptCell<T>& ac,
                                 int local_zero_entity_id,
                                 std::span<const LevelSetCell<T, I>* const> active_cells,
                                 std::span<const T> seed,
                                 const CurvingOptions<T>& options,
                                 std::vector<T>& projected,
                                 ProjectionStats<T>& stats)
{
    if (active_cells.size() == 1)
    {
        return scalar_project<T, I>(
            ac, local_zero_entity_id, *active_cells[0], seed, options, projected, stats);
    }
    return vector_project<T, I>(
        ac, local_zero_entity_id, active_cells, seed, options, projected, stats);
}

template <std::floating_point T, std::integral I>
std::vector<T> preferred_scalar_projection_direction(
    const AdaptCell<T>& ac,
    int local_zero_entity_id,
    const LevelSetCell<T, I>& ls_cell,
    std::span<const T> seed,
    const CurvingOptions<T>& options)
{
    if (local_zero_entity_id < 0
        || local_zero_entity_id >= ac.n_zero_entities())
    {
        return {};
    }

    std::vector<T> grad(static_cast<std::size_t>(ac.tdim), T(0));
    curving_level_set_grad<T, I>(
        ls_cell, seed, std::span<T>(grad.data(), grad.size()));
    const auto metric_grad = physical_normal_reference_direction<T, I>(
        ls_cell, std::span<const T>(grad.data(), grad.size()));
    const auto host_metric_grad = physical_host_gradient_reference_direction<T, I>(
        ls_cell,
        ac,
        local_zero_entity_id,
        std::span<const T>(grad.data(), grad.size()),
        std::span<const T>(metric_grad.data(), metric_grad.size()));

    std::uint32_t active_mask = structural_active_face_mask<T>(
        ac, local_zero_entity_id);
    active_mask = add_near_active_faces<T>(
        ac, seed, active_mask, options.active_face_tol);

    auto directions = scalar_candidate_directions<T>(
        ac,
        local_zero_entity_id,
        seed,
        std::span<const T>(grad.data(), grad.size()),
        std::span<const T>(metric_grad.data(), metric_grad.size()),
        std::span<const T>(host_metric_grad.data(), host_metric_grad.size()),
        active_mask,
        options);
    if (directions.empty())
        return {};
    return directions.front();
}

template <std::floating_point T, std::integral I>
bool build_hierarchical_face_nodes(
    CurvedZeroEntityState<T>& state,
    const AdaptCell<T>& ac,
    int local_zero_entity_id,
    std::span<const LevelSetCell<T, I>* const> active_cells,
    std::span<const BoundaryEdgeState<T>> boundary_edges,
    const CurvingOptions<T>& options)
{
    const int zid = ac.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    const auto entity_type = ac.entity_types[2][static_cast<std::size_t>(zid)];
    auto verts_span = ac.entity_to_vertex[2][static_cast<std::int32_t>(zid)];
    std::vector<int> verts;
    verts.reserve(verts_span.size());
    for (const auto v : verts_span)
        verts.push_back(static_cast<int>(v));

    const int order = std::max(options.geometry_order, 1);
    const T eps = T(64) * std::numeric_limits<T>::epsilon();
    std::vector<T> projected;
    state.ref_nodes.clear();

    auto append_vertex = [&](int vertex_id)
    {
        for (int d = 0; d < ac.tdim; ++d)
            state.ref_nodes.push_back(
                ac.vertex_coords[static_cast<std::size_t>(vertex_id * ac.tdim + d)]);
        append_node_stats<T>(
            state, accepted_node_stats<T>(CurvingFailureCode::exact_vertex), ac.tdim);
    };

    auto straight_seed = [&](std::span<const T> weights)
    {
        std::array<T, 3> seed = {T(0), T(0), T(0)};
        for (std::size_t v = 0; v < weights.size(); ++v)
        {
            for (int d = 0; d < ac.tdim; ++d)
            {
                seed[static_cast<std::size_t>(d)] += weights[v] * ac.vertex_coords[
                    static_cast<std::size_t>(verts[v] * ac.tdim + d)];
            }
        }
        return seed;
    };

    auto project_face_interior_seed = [&](std::span<const T> seed,
                                          std::vector<T>& candidate,
                                          ProjectionStats<T>& stats)
    {
        candidate.clear();
        stats = {};
        if (active_cells.size() != 1)
        {
            stats.failure_code = CurvingFailureCode::invalid_constraint_count;
            stats.projection_mode = CurvingProjectionMode::safe_line;
            stats.active_face_mask =
                structural_active_face_mask<T>(ac, local_zero_entity_id);
            stats.safe_subspace_dim = active_subspace_dim<T>(
                ac, stats.active_face_mask);
            return false;
        }
        return fixed_ray_scalar_project<T, I>(
            ac,
            local_zero_entity_id,
            *active_cells[0],
            seed,
            options,
            candidate,
            stats);
    };

    auto append_straight_face_interior_seed = [&](std::span<const T> weights)
    {
        const auto seed = straight_seed(weights);
        ProjectionStats<T> stats;
        if (!project_face_interior_seed(
                std::span<const T>(seed.data(), static_cast<std::size_t>(ac.tdim)),
                projected,
                stats))
        {
            append_node_stats<T>(state, stats, ac.tdim);
            return false;
        }
        state.ref_nodes.insert(state.ref_nodes.end(), projected.begin(), projected.end());
        append_node_stats<T>(state, stats, ac.tdim);
        return true;
    };

    auto append_straight_edge_node = [&](int v0, int v1, T t)
    {
        for (int d = 0; d < ac.tdim; ++d)
        {
            const T x0 = ac.vertex_coords[static_cast<std::size_t>(v0 * ac.tdim + d)];
            const T x1 = ac.vertex_coords[static_cast<std::size_t>(v1 * ac.tdim + d)];
            state.ref_nodes.push_back((T(1) - t) * x0 + t * x1);
        }
        append_node_stats<T>(
            state, accepted_node_stats<T>(CurvingFailureCode::boundary_from_edge), ac.tdim);
    };

    auto eval_curved_edge = [&](int v0, int v1, T t, std::vector<T>& x) -> bool
    {
        const auto* edge = find_boundary_edge_state<T>(boundary_edges, v0, v1);
        if (edge == nullptr || !edge->use_curved_state)
            return false;
        x = eval_edge_state_at<T>(*edge, options, v0, v1, t);
        return true;
    };

    auto add_scaled = [&](std::array<T, 3>& out,
                          T scale,
                          std::span<const T> x)
    {
        for (int d = 0; d < ac.tdim; ++d)
            out[static_cast<std::size_t>(d)] += scale * x[static_cast<std::size_t>(d)];
    };

    auto sub_scaled_vertex = [&](std::array<T, 3>& out, T scale, int vertex_id)
    {
        for (int d = 0; d < ac.tdim; ++d)
            out[static_cast<std::size_t>(d)] -= scale * ac.vertex_coords[
                static_cast<std::size_t>(vertex_id * ac.tdim + d)];
    };

    auto append_transfinite_or_straight_seed =
        [&](std::array<T, 3> transfinite_seed,
            std::span<const T> straight_weights) -> bool
    {
        // Strict face-interior nodes are lifted by a scalar solve along one
        // frozen ray through the transfinite seed. This preserves the face
        // parameter location; the generic scalar projector can otherwise
        // select fallback directions that satisfy phi=0 but drift tangentially.
        ProjectionStats<T> transfinite_stats;
        std::vector<T> transfinite_projected;
        if (project_face_interior_seed(
                std::span<const T>(
                    transfinite_seed.data(), static_cast<std::size_t>(ac.tdim)),
                transfinite_projected,
                transfinite_stats))
        {
            state.ref_nodes.insert(
                state.ref_nodes.end(),
                transfinite_projected.begin(),
                transfinite_projected.end());
            append_node_stats<T>(state, transfinite_stats, ac.tdim);
            return true;
        }

        const auto fallback_seed = straight_seed(straight_weights);
        ProjectionStats<T> fallback_stats;
        std::vector<T> fallback_projected;
        if (!project_face_interior_seed(
                std::span<const T>(
                    fallback_seed.data(), static_cast<std::size_t>(ac.tdim)),
                fallback_projected,
                fallback_stats))
        {
            append_node_stats<T>(state, fallback_stats, ac.tdim);
            return false;
        }
        state.ref_nodes.insert(
            state.ref_nodes.end(),
            fallback_projected.begin(),
            fallback_projected.end());
        append_node_stats<T>(state, fallback_stats, ac.tdim);
        return true;
    };

    auto append_triangle_interior_node = [&](std::array<T, 3> w) -> bool
    {
        const int v0 = verts[0];
        const int v1 = verts[1];
        const int v2 = verts[2];
        std::vector<T> c01;
        std::vector<T> c12;
        std::vector<T> c02;

        const T s01 = w[0] + w[1];
        const T s12 = w[1] + w[2];
        const T s02 = w[0] + w[2];
        if (s01 <= eps || s12 <= eps || s02 <= eps
            || !eval_curved_edge(v0, v1, w[1] / s01, c01)
            || !eval_curved_edge(v1, v2, w[2] / s12, c12)
            || !eval_curved_edge(v0, v2, w[2] / s02, c02))
        {
            return append_straight_face_interior_seed(
                std::span<const T>(w.data(), w.size()));
        }

        std::array<T, 3> seed = {T(0), T(0), T(0)};
        add_scaled(seed, s01, std::span<const T>(c01.data(), c01.size()));
        add_scaled(seed, s12, std::span<const T>(c12.data(), c12.size()));
        add_scaled(seed, s02, std::span<const T>(c02.data(), c02.size()));
        sub_scaled_vertex(seed, w[0], v0);
        sub_scaled_vertex(seed, w[1], v1);
        sub_scaled_vertex(seed, w[2], v2);
        return append_transfinite_or_straight_seed(
            seed, std::span<const T>(w.data(), w.size()));
    };

    auto append_quad_interior_node = [&](T u, T v) -> bool
    {
        std::array<T, 4> w = {(T(1) - u) * (T(1) - v),
                              u * (T(1) - v),
                              (T(1) - u) * v,
                              u * v};
        std::vector<T> bottom;
        std::vector<T> top;
        std::vector<T> left;
        std::vector<T> right;
        if (!eval_curved_edge(verts[0], verts[1], u, bottom)
            || !eval_curved_edge(verts[2], verts[3], u, top)
            || !eval_curved_edge(verts[0], verts[2], v, left)
            || !eval_curved_edge(verts[1], verts[3], v, right))
        {
            return append_straight_face_interior_seed(
                std::span<const T>(w.data(), w.size()));
        }

        std::array<T, 3> seed = {T(0), T(0), T(0)};
        add_scaled(seed, T(1) - v, std::span<const T>(bottom.data(), bottom.size()));
        add_scaled(seed, v, std::span<const T>(top.data(), top.size()));
        add_scaled(seed, T(1) - u, std::span<const T>(left.data(), left.size()));
        add_scaled(seed, u, std::span<const T>(right.data(), right.size()));
        sub_scaled_vertex(seed, w[0], verts[0]);
        sub_scaled_vertex(seed, w[1], verts[1]);
        sub_scaled_vertex(seed, w[2], verts[2]);
        sub_scaled_vertex(seed, w[3], verts[3]);
        return append_transfinite_or_straight_seed(
            seed, std::span<const T>(w.data(), w.size()));
    };

    auto append_triangle_node = [&](std::array<T, 3> w) -> bool
    {
        int zero_count = 0;
        int one_index = -1;
        int edge_a = -1;
        int edge_b = -1;
        for (int k = 0; k < 3; ++k)
        {
            if (std::abs(w[static_cast<std::size_t>(k)]) <= eps)
                ++zero_count;
            if (std::abs(w[static_cast<std::size_t>(k)] - T(1)) <= eps)
                one_index = k;
            else if (w[static_cast<std::size_t>(k)] > eps)
            {
                if (edge_a < 0)
                    edge_a = k;
                else
                    edge_b = k;
            }
        }

        if (one_index >= 0)
        {
            append_vertex(verts[static_cast<std::size_t>(one_index)]);
            return true;
        }
        if (zero_count == 1 && edge_a >= 0 && edge_b >= 0)
        {
            const int v0 = verts[static_cast<std::size_t>(edge_a)];
            const int v1 = verts[static_cast<std::size_t>(edge_b)];
            const T denom = w[static_cast<std::size_t>(edge_a)]
                          + w[static_cast<std::size_t>(edge_b)];
            const T t = w[static_cast<std::size_t>(edge_b)] / denom;
            const auto* edge = find_boundary_edge_state<T>(boundary_edges, v0, v1);
            if (edge == nullptr)
            {
                ProjectionStats<T> stats;
                stats.failure_code = CurvingFailureCode::missing_boundary_edge;
                append_node_stats<T>(state, stats, ac.tdim);
                return false;
            }
            if (!edge->use_curved_state)
            {
                append_straight_edge_node(v0, v1, t);
                return true;
            }
            const auto x = eval_edge_state_at<T>(*edge, options, v0, v1, t);
            state.ref_nodes.insert(state.ref_nodes.end(), x.begin(), x.end());
            append_node_stats<T>(
                state, accepted_node_stats<T>(CurvingFailureCode::boundary_from_edge), ac.tdim);
            return true;
        }
        return append_triangle_interior_node(w);
    };

    if (entity_type == cell::type::triangle)
    {
        const auto nodes =
            triangle_interpolation_barycentric_nodes<T>(order, options.node_family);
        for (const auto& w : nodes)
        {
            if (!append_triangle_node(w))
            {
                state.failure_reason = "missing or failed hierarchical curved face boundary edge";
                return false;
            }
        }
        return true;
    }

    if (entity_type == cell::type::quadrilateral)
    {
        const auto params = interpolation_parameters<T>(order, options.node_family);
        for (const T v : params)
        {
            for (const T u : params)
            {
                const bool on_u0 = std::abs(u) <= eps;
                const bool on_u1 = std::abs(u - T(1)) <= eps;
                const bool on_v0 = std::abs(v) <= eps;
                const bool on_v1 = std::abs(v - T(1)) <= eps;

                if ((on_u0 || on_u1) && (on_v0 || on_v1))
                {
                    int vertex = -1;
                    if (on_u0 && on_v0)
                    {
                        vertex = verts[0];
                    }
                    else if (on_u1 && on_v0)
                    {
                        vertex = verts[1];
                    }
                    else if (on_u0 && on_v1)
                    {
                        vertex = verts[2];
                    }
                    else
                    {
                        vertex = verts[3];
                    }
                    append_vertex(vertex);
                    continue;
                }

                if (on_v0 || on_v1 || on_u0 || on_u1)
                {
                    int v0 = -1;
                    int v1 = -1;
                    T t = T(0);
                    if (on_v0)
                    {
                        v0 = verts[0]; v1 = verts[1]; t = u;
                    }
                    else if (on_v1)
                    {
                        v0 = verts[2]; v1 = verts[3]; t = u;
                    }
                    else if (on_u0)
                    {
                        v0 = verts[0]; v1 = verts[2]; t = v;
                    }
                    else
                    {
                        v0 = verts[1]; v1 = verts[3]; t = v;
                    }
                    const auto* edge = find_boundary_edge_state<T>(boundary_edges, v0, v1);
                    if (edge == nullptr)
                    {
                        ProjectionStats<T> stats;
                        stats.failure_code = CurvingFailureCode::missing_boundary_edge;
                        append_node_stats<T>(state, stats, ac.tdim);
                        state.failure_reason = "missing or failed hierarchical curved face boundary edge";
                        return false;
                    }
                    if (!edge->use_curved_state)
                    {
                        append_straight_edge_node(v0, v1, t);
                        continue;
                    }
                    const auto x = eval_edge_state_at<T>(*edge, options, v0, v1, t);
                    state.ref_nodes.insert(state.ref_nodes.end(), x.begin(), x.end());
                    append_node_stats<T>(
                        state, accepted_node_stats<T>(CurvingFailureCode::boundary_from_edge), ac.tdim);
                    continue;
                }

                if (!append_quad_interior_node(u, v))
                {
                    state.failure_reason = "projection failed inside zero-entity host domain";
                    return false;
                }
            }
        }
        return true;
    }

    state.failure_reason = "curving: unsupported zero-face type";
    return false;
}

template <std::floating_point T, std::integral I>
void build_curved_state(CurvedZeroEntityState<T>& state,
                        const AdaptCell<T>& ac,
                        std::span<const LevelSetCell<T, I>> level_set_cells,
                        std::span<const int> ls_offsets,
                        int cut_cell_id,
                        int local_zero_entity_id,
                        std::span<const BoundaryEdgeState<T>> boundary_edges,
                        const CurvingOptions<T>& options)
{
    state.status = CurvingStatus::in_progress;
    state.failure_reason.clear();
    state.geometry_order = options.geometry_order;
    state.node_family = options.node_family;
    state.direction_mode = options.direction_mode;
    state.small_entity_tol = options.small_entity_tol;
    state.zero_entity_version = ac.zero_entity_version;
    state.zero_mask = ac.zero_entity_zero_mask[static_cast<std::size_t>(local_zero_entity_id)];
    state.ref_nodes.clear();
    state.node_iterations.clear();
    state.node_status.clear();
    state.node_failure_code.clear();
    state.node_residual.clear();
    state.node_active_face_mask.clear();
    state.node_closest_face_id.clear();
    state.node_safe_subspace_dim.clear();
    state.node_projection_mode.clear();
    state.node_retry_count.clear();
    state.node_seed.clear();
    state.node_direction.clear();
    state.node_clip_lo.clear();
    state.node_clip_hi.clear();
    state.node_root_t.clear();

    std::vector<T> seeds;
    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim == 0)
    {
        state.ref_nodes = zero_entity_vertices<T>(ac, local_zero_entity_id);
        append_node_stats<T>(
            state, accepted_node_stats<T>(CurvingFailureCode::exact_vertex), ac.tdim);
        state.status = CurvingStatus::curved;
        return;
    }
    if (zdim == 1)
        append_edge_seed_nodes<T>(ac, local_zero_entity_id, options, seeds);
    else if (zdim == 2)
    {
        // Face states are built hierarchically below after active constraints
        // are available: vertices are kept exact, boundary nodes come from
        // already curved zero edges, and only strict face-interior nodes are
        // projected.
    }
    else
    {
        state.status = CurvingStatus::failed;
        state.failure_reason = "unsupported zero-entity dimension";
        ProjectionStats<T> stats;
        stats.failure_code = CurvingFailureCode::unsupported_entity;
        append_node_stats<T>(state, stats, ac.tdim);
        return;
    }

    if (zero_entity_below_curving_threshold<T>(ac, local_zero_entity_id, options))
    {
        if (zdim == 2)
            append_face_seed_nodes<T>(ac, local_zero_entity_id, options, seeds);
        state.ref_nodes = seeds;
        append_straight_seed_node_stats<T>(
            state,
            ac,
            std::span<const T>(state.ref_nodes.data(), state.ref_nodes.size()));
        state.status = CurvingStatus::curved;
        return;
    }

    const std::uint64_t effective_zero_mask = state.zero_mask & ac.active_level_set_mask;
    const auto ls_ids = active_level_sets<T>(effective_zero_mask);
    if (ls_ids.empty())
    {
        if (zdim == 2)
            append_face_seed_nodes<T>(ac, local_zero_entity_id, options, seeds);
        state.ref_nodes = seeds;
        const int n_nodes = static_cast<int>(
            seeds.size() / static_cast<std::size_t>(std::max(ac.tdim, 1)));
        for (int i = 0; i < n_nodes; ++i)
            append_node_stats<T>(
                state, accepted_node_stats<T>(CurvingFailureCode::none), ac.tdim);
        state.status = CurvingStatus::curved;
        return;
    }

    std::vector<const LevelSetCell<T, I>*> active_cells;
    active_cells.reserve(ls_ids.size());
    try
    {
        for (const int ls_id : ls_ids)
            active_cells.push_back(&find_level_set_cell<T, I>(
                level_set_cells, ls_offsets, cut_cell_id, ls_id));
    }
    catch (const std::exception& e)
    {
        state.status = CurvingStatus::failed;
        state.failure_reason = e.what();
        ProjectionStats<T> stats;
        stats.failure_code = CurvingFailureCode::missing_level_set_cell;
        append_node_stats<T>(state, stats, ac.tdim);
        return;
    }

    if (zdim == 2)
    {
        if (!build_hierarchical_face_nodes<T, I>(
                state,
                ac,
                local_zero_entity_id,
                std::span<const LevelSetCell<T, I>* const>(
                    active_cells.data(), active_cells.size()),
                boundary_edges,
                options))
        {
            state.status = CurvingStatus::failed;
            if (state.failure_reason.empty())
                state.failure_reason = "hierarchical zero-face curving failed";
            state.ref_nodes.clear();
            return;
        }
        state.status = CurvingStatus::curved;
        return;
    }

    const int nseeds = static_cast<int>(seeds.size() / static_cast<std::size_t>(ac.tdim));
    state.ref_nodes.reserve(seeds.size());
    std::vector<T> projected;
    for (int i = 0; i < nseeds; ++i)
    {
        std::span<const T> seed(
            seeds.data() + static_cast<std::size_t>(i * ac.tdim),
            static_cast<std::size_t>(ac.tdim));

        if (zdim == 1 && (i == 0 || i == nseeds - 1))
        {
            state.ref_nodes.insert(state.ref_nodes.end(), seed.begin(), seed.end());
            append_node_stats<T>(
                state,
                accepted_node_stats<T>(CurvingFailureCode::exact_vertex),
                ac.tdim);
            continue;
        }

        ProjectionStats<T> stats;
        const bool ok = project_seed_to_zero_entity<T, I>(
            ac,
            local_zero_entity_id,
            std::span<const LevelSetCell<T, I>* const>(
                active_cells.data(), active_cells.size()),
            seed,
            options,
            projected,
            stats);
        append_node_stats<T>(state, stats, ac.tdim);

        if (!ok)
        {
            state.status = CurvingStatus::failed;
            state.failure_reason = "projection failed inside zero-entity host domain";
            state.ref_nodes.clear();
            return;
        }
        state.ref_nodes.insert(state.ref_nodes.end(), projected.begin(), projected.end());
    }

    state.status = CurvingStatus::curved;
}

} // namespace

NodeFamily node_family_from_string(std::string_view name)
{
    if (name == "gll" || name == "lobatto"
        || name == "fekete" || name == "warp_blend" || name == "warp-blend")
        return NodeFamily::gll;
    if (name == "equispaced")
        return NodeFamily::equispaced;
    if (name == "lagrange")
        return NodeFamily::lagrange;
    throw std::invalid_argument("unknown curving node family");
}

std::string_view node_family_name(NodeFamily family)
{
    switch (family)
    {
    case NodeFamily::gll: return "gll";
    case NodeFamily::equispaced: return "equispaced";
    case NodeFamily::lagrange: return "lagrange";
    }
    return "unknown";
}

CurvingDirectionMode direction_mode_from_string(std::string_view name)
{
    if (name == "straight_zero_entity_normal"
        || name == "straight-zero-entity-normal"
        || name == "zero_entity_normal"
        || name == "zero-entity-normal"
        || name == "straight_normal"
        || name == "straight-normal"
        || name == "normal")
    {
        return CurvingDirectionMode::straight_zero_entity_normal;
    }
    if (name == "level_set_gradient"
        || name == "level-set-gradient"
        || name == "gradient")
    {
        return CurvingDirectionMode::level_set_gradient;
    }
    throw std::invalid_argument("unknown curving direction mode");
}

std::string_view direction_mode_name(CurvingDirectionMode mode)
{
    switch (mode)
    {
    case CurvingDirectionMode::straight_zero_entity_normal:
        return "straight_zero_entity_normal";
    case CurvingDirectionMode::level_set_gradient:
        return "level_set_gradient";
    }
    return "unknown";
}

template <std::floating_point T, std::integral I>
T reference_level_set_value(const LevelSetCell<T, I>& ls_cell,
                            std::span<const T> ref)
{
    return curving_level_set_value<T, I>(ls_cell, ref);
}

template <std::floating_point T, std::integral I>
void reference_level_set_gradient(const LevelSetCell<T, I>& ls_cell,
                                  std::span<const T> ref,
                                  std::span<T> grad_ref)
{
    curving_level_set_grad<T, I>(ls_cell, ref, grad_ref);
}

template <std::floating_point T>
ProjectionDiagnostic<T> make_projection_diagnostic(
    bool accepted,
    const std::vector<T>& projected,
    const ProjectionStats<T>& stats)
{
    ProjectionDiagnostic<T> diagnostic;
    diagnostic.accepted = accepted;
    diagnostic.status = stats.status;
    diagnostic.failure_code = stats.failure_code;
    diagnostic.projection_mode = stats.projection_mode;
    diagnostic.iterations = stats.iterations;
    diagnostic.residual = stats.residual;
    diagnostic.active_face_mask = stats.active_face_mask;
    diagnostic.closest_face_id = stats.closest_face_id;
    diagnostic.safe_subspace_dim = stats.safe_subspace_dim;
    diagnostic.retry_count = stats.retry_count;
    diagnostic.projected = projected;
    diagnostic.seed = stats.seed;
    diagnostic.direction = stats.direction;
    diagnostic.clip_lo = stats.clip_lo;
    diagnostic.clip_hi = stats.clip_hi;
    diagnostic.root_t = stats.root_t;
    return diagnostic;
}

template <std::floating_point T, std::integral I>
ProjectionDiagnostic<T> project_seed_to_zero_entity_diagnostic(
    const AdaptCell<T>& ac,
    int local_zero_entity_id,
    const LevelSetCell<T, I>& ls_cell,
    std::span<const T> seed,
    const CurvingOptions<T>& options)
{
    std::array<const LevelSetCell<T, I>*, 1> active_cells = {&ls_cell};
    std::vector<T> projected;
    ProjectionStats<T> stats;
    bool accepted = false;
    try
    {
        accepted = project_seed_to_zero_entity<T, I>(
            ac,
            local_zero_entity_id,
            std::span<const LevelSetCell<T, I>* const>(
                active_cells.data(), active_cells.size()),
            seed,
            options,
            projected,
            stats);
    }
    catch (const std::exception&)
    {
        accepted = false;
        if (stats.failure_code == CurvingFailureCode::none)
            stats.failure_code = CurvingFailureCode::projection_failed;
        stats.status = CurvingStatus::failed;
    }

    if (accepted && stats.direction.empty())
    {
        stats.direction = normalized_delta_direction<T>(
            seed,
            std::span<const T>(projected.data(), projected.size()),
            options.xtol);
    }
    if (stats.direction.empty())
    {
        stats.direction = preferred_scalar_projection_direction<T, I>(
            ac, local_zero_entity_id, ls_cell, seed, options);
    }

    return make_projection_diagnostic<T>(accepted, projected, stats);
}

template <std::floating_point T, std::integral I>
void rebuild_identity(CurvingData<T, I>& curving,
                      std::span<const I> parent_cell_ids,
                      std::span<const AdaptCell<T>> adapt_cells)
{
    curving.clear();
    curving.num_cut_cells = static_cast<int>(adapt_cells.size());
    curving.local_to_canonical_offsets.reserve(adapt_cells.size() + 1);
    curving.local_to_canonical_offsets.push_back(0);

    for (int c = 0; c < static_cast<int>(adapt_cells.size()); ++c)
    {
        const auto& ac = adapt_cells[static_cast<std::size_t>(c)];
        for (int z = 0; z < ac.n_zero_entities(); ++z)
        {
            CurvingIdentity id;
            id.cut_cell_id = c;
            id.local_zero_entity_id = z;
            id.dim = ac.zero_entity_dim[static_cast<std::size_t>(z)];
            id.parent_dim = ac.zero_entity_parent_dim[static_cast<std::size_t>(z)];
            id.parent_id = ac.zero_entity_parent_id[static_cast<std::size_t>(z)];
            if (id.parent_dim == ac.tdim && id.parent_id < 0)
                id.parent_id = static_cast<std::int32_t>(parent_cell_ids[static_cast<std::size_t>(c)]);
            id.zero_mask = ac.zero_entity_zero_mask[static_cast<std::size_t>(z)];
            curving.identities.push_back(id);
            curving.states.emplace_back();
            curving.local_to_canonical.push_back(
                static_cast<int>(curving.identities.size()) - 1);
        }
        curving.local_to_canonical_offsets.push_back(
            static_cast<int>(curving.local_to_canonical.size()));
    }

    curving.identity_valid = true;
}

template <std::floating_point T>
int host_face_for_zero_face_boundary_edge(const AdaptCell<T>& ac,
                                          int zero_face_id,
                                          int zero_edge_id,
                                          T tol);

template <std::floating_point T>
std::vector<BoundaryEdgeRef> face_boundary_zero_edges(const AdaptCell<T>& ac,
                                                      int local_zero_entity_id)
{
    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)];
    if (zdim != 2)
        return {};

    const std::uint64_t face_mask =
        ac.zero_entity_zero_mask[static_cast<std::size_t>(local_zero_entity_id)];
    const int zid = ac.zero_entity_id[static_cast<std::size_t>(local_zero_entity_id)];
    const auto face_type = ac.entity_types[2][static_cast<std::size_t>(zid)];
    auto face_verts_span = ac.entity_to_vertex[2][static_cast<std::int32_t>(zid)];
    std::vector<int> face_verts;
    face_verts.reserve(face_verts_span.size());
    for (const auto v : face_verts_span)
        face_verts.push_back(static_cast<int>(v));

    std::vector<BoundaryEdgeRef> edge_ids;
    for (const auto& edge : cell::edges(face_type))
    {
        std::array<int, 2> target = {
            face_verts[static_cast<std::size_t>(edge[0])],
            face_verts[static_cast<std::size_t>(edge[1])]
        };
        int match = -1;
        for (int z = 0; z < ac.n_zero_entities(); ++z)
        {
            if (ac.zero_entity_dim[static_cast<std::size_t>(z)] != 1)
                continue;
            const auto edge_mask = ac.zero_entity_zero_mask[static_cast<std::size_t>(z)];
            if ((edge_mask & face_mask) != face_mask)
                continue;
            const auto edge_verts = zero_entity_vertex_ids<T>(ac, z);
            if (same_unordered_vertices(
                    std::span<const int>(target.data(), target.size()),
                    std::span<const int>(edge_verts.data(), edge_verts.size())))
            {
                match = z;
                break;
            }
        }
        if (match < 0)
            return {};
        const int host_face_id =
            host_face_for_zero_face_boundary_edge<T>(
                ac,
                local_zero_entity_id,
                match,
                T(256) * std::numeric_limits<T>::epsilon());
        // Triangulated zero quadrilaterals introduce an artificial diagonal.
        // It is a zero subedge of the output triangle, but it is not on a
        // host-cell boundary face. Keep it as a curved zero edge; a negative
        // host_face_id makes its contextual host the uncut subcell volume.
        edge_ids.push_back({match, host_face_id, true});
    }
    return edge_ids;
}

template <std::floating_point T>
bool point_in_triangle_ref(const AdaptCell<T>& ac,
                           std::span<const int> tri_vertices,
                           std::span<const T> point,
                           T tol)
{
    if (ac.tdim != 3 || tri_vertices.size() != 3 || point.size() != 3)
        return false;
    const auto a = vertex_ref_point<T>(ac, tri_vertices[0]);
    const auto b = vertex_ref_point<T>(ac, tri_vertices[1]);
    const auto c = vertex_ref_point<T>(ac, tri_vertices[2]);
    std::array<T, 3> v0 = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    std::array<T, 3> v1 = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    std::array<T, 3> v2 = {point[0] - a[0], point[1] - a[1], point[2] - a[2]};
    const auto n = local_cross<T>(
        std::span<const T>(v0.data(), 3),
        std::span<const T>(v1.data(), 3));
    const T area2 = local_norm<T>(std::span<const T>(n.data(), 3));
    if (area2 <= tol)
        return false;
    const T plane =
        std::fabs(local_dot<T>(
            std::span<const T>(v2.data(), 3),
            std::span<const T>(n.data(), 3))) / area2;
    if (plane > tol)
        return false;

    const T d00 = local_dot<T>(
        std::span<const T>(v0.data(), 3),
        std::span<const T>(v0.data(), 3));
    const T d01 = local_dot<T>(
        std::span<const T>(v0.data(), 3),
        std::span<const T>(v1.data(), 3));
    const T d11 = local_dot<T>(
        std::span<const T>(v1.data(), 3),
        std::span<const T>(v1.data(), 3));
    const T d20 = local_dot<T>(
        std::span<const T>(v2.data(), 3),
        std::span<const T>(v0.data(), 3));
    const T d21 = local_dot<T>(
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
int host_face_for_zero_face_boundary_edge(const AdaptCell<T>& ac,
                                          int zero_face_id,
                                          int zero_edge_id,
                                          T tol)
{
    if (zero_face_id < 0 || zero_edge_id < 0
        || zero_face_id >= ac.n_zero_entities()
        || zero_edge_id >= ac.n_zero_entities())
    {
        return -1;
    }
    if (zero_face_id >= static_cast<int>(ac.zero_entity_host_cell_type.size()))
        return -1;
    const cell::type host_type =
        ac.zero_entity_host_cell_type[static_cast<std::size_t>(zero_face_id)];
    const auto host_cell_vertices =
        zero_entity_host_cell_vertex_ids<T>(ac, zero_face_id);
    const auto edge_vertices =
        local_zero_entity_vertex_ids<T>(ac, zero_edge_id);
    if (host_cell_vertices.empty() || edge_vertices.size() != 2
        || host_type == cell::type::point)
    {
        return -1;
    }

    const auto p0 = vertex_ref_point<T>(ac, edge_vertices[0]);
    const auto p1 = vertex_ref_point<T>(ac, edge_vertices[1]);
    for (int f = 0; f < cell::num_faces(host_type); ++f)
    {
        auto local_face = cell::face_vertices(host_type, f);
        if (local_face.size() != 3)
            continue;
        std::array<int, 3> face_vertices = {
            host_cell_vertices[static_cast<std::size_t>(local_face[0])],
            host_cell_vertices[static_cast<std::size_t>(local_face[1])],
            host_cell_vertices[static_cast<std::size_t>(local_face[2])]
        };
        if (point_in_triangle_ref<T>(
                ac,
                std::span<const int>(face_vertices.data(), face_vertices.size()),
                p0,
                tol)
            && point_in_triangle_ref<T>(
                ac,
                std::span<const int>(face_vertices.data(), face_vertices.size()),
                p1,
                tol))
        {
            return f;
        }
    }
    return -1;
}

template <std::floating_point T>
void override_zero_edge_host_from_face(AdaptCell<T>& ac,
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
    if (zero_face_id >= static_cast<int>(ac.zero_entity_host_cell_id.size())
        || zero_edge_id >= static_cast<int>(ac.zero_entity_host_cell_id.size()))
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
            const auto old = ac.zero_entity_host_cell_vertices[static_cast<std::int32_t>(z)];
            for (const auto v : old)
                updated.indices.push_back(v);
        }
        updated.offsets.push_back(
            static_cast<std::int32_t>(updated.indices.size()));
    }
    ac.zero_entity_host_cell_vertices = std::move(updated);
}

template <std::floating_point T, std::integral I>
const CurvedZeroEntityState<T>& ensure_curved(
    CurvingData<T, I>& curving,
    std::span<const I> parent_cell_ids,
    std::span<const AdaptCell<T>> adapt_cells,
    std::span<const LevelSetCell<T, I>> level_set_cells,
    std::span<const int> ls_offsets,
    int cut_cell_id,
    int local_zero_entity_id,
    const CurvingOptions<T>& options)
{
    if (!curving.identity_valid
        || curving.num_cut_cells != static_cast<int>(adapt_cells.size()))
    {
        rebuild_identity<T, I>(curving, parent_cell_ids, adapt_cells);
    }
    if (cut_cell_id < 0 || cut_cell_id >= static_cast<int>(adapt_cells.size()))
        throw std::out_of_range("curving: cut_cell_id out of range");

    const int begin = curving.local_to_canonical_offsets[static_cast<std::size_t>(cut_cell_id)];
    const int end = curving.local_to_canonical_offsets[static_cast<std::size_t>(cut_cell_id + 1)];
    if (local_zero_entity_id < 0 || begin + local_zero_entity_id >= end)
        throw std::out_of_range("curving: local_zero_entity_id out of range");

    const int canonical = curving.local_to_canonical[static_cast<std::size_t>(begin + local_zero_entity_id)];
    auto& state = curving.states[static_cast<std::size_t>(canonical)];
    const auto& ac = adapt_cells[static_cast<std::size_t>(cut_cell_id)];
    const bool valid =
        state.status == CurvingStatus::curved
        && state.geometry_order == options.geometry_order
        && state.node_family == options.node_family
        && state.direction_mode == options.direction_mode
        && state.small_entity_tol == options.small_entity_tol
        && state.zero_entity_version == ac.zero_entity_version;
    const bool failed =
        state.status == CurvingStatus::failed
        && state.geometry_order == options.geometry_order
        && state.node_family == options.node_family
        && state.direction_mode == options.direction_mode
        && state.small_entity_tol == options.small_entity_tol
        && state.zero_entity_version == ac.zero_entity_version;
    if (!valid && !failed)
    {
        std::vector<BoundaryEdgeState<T>> boundary_edges;
        std::vector<CurvedZeroEntityState<T>> contextual_boundary_edge_states;
        if (ac.zero_entity_dim[static_cast<std::size_t>(local_zero_entity_id)] == 2
            && !zero_entity_below_curving_threshold<T>(ac, local_zero_entity_id, options))
        {
            const auto boundary_edge_refs =
                face_boundary_zero_edges<T>(ac, local_zero_entity_id);
            if (boundary_edge_refs.empty())
            {
                state.status = CurvingStatus::failed;
                state.failure_reason = "hierarchical zero-face curving requires boundary zero edges";
                state.geometry_order = options.geometry_order;
                state.node_family = options.node_family;
                state.direction_mode = options.direction_mode;
                state.small_entity_tol = options.small_entity_tol;
                state.zero_entity_version = ac.zero_entity_version;
                state.zero_mask = ac.zero_entity_zero_mask[
                    static_cast<std::size_t>(local_zero_entity_id)];
                state.ref_nodes.clear();
                state.node_iterations.clear();
                state.node_status.clear();
                state.node_failure_code.clear();
                state.node_residual.clear();
                state.node_active_face_mask.clear();
                state.node_closest_face_id.clear();
                state.node_safe_subspace_dim.clear();
                state.node_projection_mode.clear();
                state.node_retry_count.clear();
                state.node_seed.clear();
                state.node_direction.clear();
                state.node_clip_lo.clear();
                state.node_clip_hi.clear();
                state.node_root_t.clear();
                ProjectionStats<T> stats;
                stats.failure_code = CurvingFailureCode::missing_boundary_edge;
                append_node_stats<T>(state, stats, ac.tdim);
                return state;
            }

            boundary_edges.reserve(boundary_edge_refs.size());
            contextual_boundary_edge_states.reserve(boundary_edge_refs.size());
            for (const auto& edge_ref : boundary_edge_refs)
            {
                const CurvedZeroEntityState<T>* edge_state = nullptr;
                if (edge_ref.use_curved_state)
                {
                    AdaptCell<T> edge_context_ac = ac;
                    override_zero_edge_host_from_face<T>(
                        edge_context_ac,
                        edge_ref.local_zero_entity_id,
                        local_zero_entity_id,
                        edge_ref.host_face_id);
                    contextual_boundary_edge_states.emplace_back();
                    auto& mutable_edge_state =
                        contextual_boundary_edge_states.back();
                    build_curved_state<T, I>(
                        mutable_edge_state,
                        edge_context_ac,
                        level_set_cells,
                        ls_offsets,
                        cut_cell_id,
                        edge_ref.local_zero_entity_id,
                        std::span<const BoundaryEdgeState<T>>(),
                        options);
                    edge_state = &mutable_edge_state;
                    if (edge_state->status != CurvingStatus::curved)
                    {
                        state.status = CurvingStatus::failed;
                        state.failure_reason = "hierarchical zero-face boundary edge curving failed";
                        state.geometry_order = options.geometry_order;
                        state.node_family = options.node_family;
                        state.direction_mode = options.direction_mode;
                        state.small_entity_tol = options.small_entity_tol;
                        state.zero_entity_version = ac.zero_entity_version;
                        state.zero_mask = ac.zero_entity_zero_mask[
                            static_cast<std::size_t>(local_zero_entity_id)];
                        state.ref_nodes.clear();
                        state.node_iterations.clear();
                        state.node_status.clear();
                        state.node_failure_code.clear();
                        state.node_residual.clear();
                        state.node_active_face_mask.clear();
                        state.node_closest_face_id.clear();
                        state.node_safe_subspace_dim.clear();
                        state.node_projection_mode.clear();
                        state.node_retry_count.clear();
                        state.node_seed.clear();
                        state.node_direction.clear();
                        state.node_clip_lo.clear();
                        state.node_clip_hi.clear();
                        state.node_root_t.clear();
                        ProjectionStats<T> stats;
                        stats.failure_code = CurvingFailureCode::boundary_edge_failed;
                        append_node_stats<T>(state, stats, ac.tdim);
                        return state;
                    }
                }
                const auto edge_verts =
                    zero_entity_vertex_ids<T>(ac, edge_ref.local_zero_entity_id);
                BoundaryEdgeState<T> boundary;
                boundary.vertices = {edge_verts[0], edge_verts[1]};
                boundary.state = edge_state;
                boundary.use_curved_state = edge_ref.use_curved_state;
                boundary_edges.push_back(boundary);
            }
        }

        build_curved_state<T, I>(
            state, ac, level_set_cells, ls_offsets,
            cut_cell_id, local_zero_entity_id,
            std::span<const BoundaryEdgeState<T>>(
                boundary_edges.data(), boundary_edges.size()),
            options);
    }
    return state;
}

template <std::floating_point T, std::integral I>
void ensure_all_curved(CurvingData<T, I>& curving,
                       std::span<const I> parent_cell_ids,
                       std::span<const AdaptCell<T>> adapt_cells,
                       std::span<const LevelSetCell<T, I>> level_set_cells,
                       std::span<const int> ls_offsets,
                       const CurvingOptions<T>& options)
{
    if (!curving.identity_valid
        || curving.num_cut_cells != static_cast<int>(adapt_cells.size()))
    {
        rebuild_identity<T, I>(curving, parent_cell_ids, adapt_cells);
    }
    for (int c = 0; c < static_cast<int>(adapt_cells.size()); ++c)
    {
        for (int z = 0; z < adapt_cells[static_cast<std::size_t>(c)].n_zero_entities(); ++z)
            (void)ensure_curved<T, I>(
                curving, parent_cell_ids, adapt_cells, level_set_cells, ls_offsets,
                c, z, options);
    }
}

template double reference_level_set_value(
    const LevelSetCell<double, int>&, std::span<const double>);
template float reference_level_set_value(
    const LevelSetCell<float, int>&, std::span<const float>);
template double reference_level_set_value(
    const LevelSetCell<double, long>&, std::span<const double>);
template float reference_level_set_value(
    const LevelSetCell<float, long>&, std::span<const float>);

template void reference_level_set_gradient(
    const LevelSetCell<double, int>&, std::span<const double>, std::span<double>);
template void reference_level_set_gradient(
    const LevelSetCell<float, int>&, std::span<const float>, std::span<float>);
template void reference_level_set_gradient(
    const LevelSetCell<double, long>&, std::span<const double>, std::span<double>);
template void reference_level_set_gradient(
    const LevelSetCell<float, long>&, std::span<const float>, std::span<float>);

template ProjectionDiagnostic<double> project_seed_to_zero_entity_diagnostic(
    const AdaptCell<double>&, int, const LevelSetCell<double, int>&,
    std::span<const double>, const CurvingOptions<double>&);
template ProjectionDiagnostic<float> project_seed_to_zero_entity_diagnostic(
    const AdaptCell<float>&, int, const LevelSetCell<float, int>&,
    std::span<const float>, const CurvingOptions<float>&);
template ProjectionDiagnostic<double> project_seed_to_zero_entity_diagnostic(
    const AdaptCell<double>&, int, const LevelSetCell<double, long>&,
    std::span<const double>, const CurvingOptions<double>&);
template ProjectionDiagnostic<float> project_seed_to_zero_entity_diagnostic(
    const AdaptCell<float>&, int, const LevelSetCell<float, long>&,
    std::span<const float>, const CurvingOptions<float>&);

template void rebuild_identity(CurvingData<double, int>&, std::span<const int>, std::span<const AdaptCell<double>>);
template void rebuild_identity(CurvingData<float, int>&, std::span<const int>, std::span<const AdaptCell<float>>);
template void rebuild_identity(CurvingData<double, long>&, std::span<const long>, std::span<const AdaptCell<double>>);
template void rebuild_identity(CurvingData<float, long>&, std::span<const long>, std::span<const AdaptCell<float>>);

template const CurvedZeroEntityState<double>& ensure_curved(
    CurvingData<double, int>&, std::span<const int>, std::span<const AdaptCell<double>>,
    std::span<const LevelSetCell<double, int>>, std::span<const int>, int, int,
    const CurvingOptions<double>&);
template const CurvedZeroEntityState<float>& ensure_curved(
    CurvingData<float, int>&, std::span<const int>, std::span<const AdaptCell<float>>,
    std::span<const LevelSetCell<float, int>>, std::span<const int>, int, int,
    const CurvingOptions<float>&);
template const CurvedZeroEntityState<double>& ensure_curved(
    CurvingData<double, long>&, std::span<const long>, std::span<const AdaptCell<double>>,
    std::span<const LevelSetCell<double, long>>, std::span<const int>, int, int,
    const CurvingOptions<double>&);
template const CurvedZeroEntityState<float>& ensure_curved(
    CurvingData<float, long>&, std::span<const long>, std::span<const AdaptCell<float>>,
    std::span<const LevelSetCell<float, long>>, std::span<const int>, int, int,
    const CurvingOptions<float>&);

template void ensure_all_curved(
    CurvingData<double, int>&, std::span<const int>, std::span<const AdaptCell<double>>,
    std::span<const LevelSetCell<double, int>>, std::span<const int>,
    const CurvingOptions<double>&);
template void ensure_all_curved(
    CurvingData<float, int>&, std::span<const int>, std::span<const AdaptCell<float>>,
    std::span<const LevelSetCell<float, int>>, std::span<const int>,
    const CurvingOptions<float>&);
template void ensure_all_curved(
    CurvingData<double, long>&, std::span<const long>, std::span<const AdaptCell<double>>,
    std::span<const LevelSetCell<double, long>>, std::span<const int>,
    const CurvingOptions<double>&);
template void ensure_all_curved(
    CurvingData<float, long>&, std::span<const long>, std::span<const AdaptCell<float>>,
    std::span<const LevelSetCell<float, long>>, std::span<const int>,
    const CurvingOptions<float>&);

} // namespace cutcells::curving
