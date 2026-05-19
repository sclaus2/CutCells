// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "implicit_quadrature.h"

#include "cell_certification.h"
#include "cell_flags.h"
#include "cell_topology.h"
#include "edge_certification.h"
#include "mapping.h"
#include "quadrature_tables.h"
#include "reference_cell.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cutcells::implicit_quadrature
{
namespace
{

enum class ChartKind : std::uint8_t
{
    triangle_apex = 0,
    tetra_apex_1_3 = 1,
    tetra_edge_edge_2_2 = 2,
    projected_height_direction = 3
};

enum DebugCandidateMask : std::int32_t
{
    debug_candidate_duffy = 1 << 0,
    debug_candidate_refined = 1 << 1,
    debug_candidate_height = 1 << 2,
    debug_candidate_linear_fallback = 1 << 3
};

enum DebugRejectionReason : std::int32_t
{
    debug_rejection_none = 0,
    debug_rejection_area_mismatch = 1 << 0,
    debug_rejection_q_error = 1 << 1,
    debug_rejection_hard_fallback = 1 << 2,
    debug_rejection_chart_line = 1 << 3,
    debug_rejection_height_failed = 1 << 4,
    debug_rejection_max_depth = 1 << 5,
    debug_rejection_topology = 1 << 6,
    debug_rejection_chart_endpoint = 1 << 7,
    debug_rejection_chart_root_structure = 1 << 8,
    debug_rejection_chart_root_solve = 1 << 9,
    debug_rejection_chart_transversality = 1 << 10
};

constexpr std::uint64_t fnv_offset_basis = 14695981039346656037ull;
constexpr std::uint64_t fnv_prime = 1099511628211ull;

inline void hash_u64(std::uint64_t& hash, std::uint64_t value)
{
    hash ^= value;
    hash *= fnv_prime;
}

inline void hash_i64(std::uint64_t& hash, std::int64_t value)
{
    hash_u64(hash, static_cast<std::uint64_t>(value));
}

template <std::floating_point T>
void hash_quantized(std::uint64_t& hash, T value)
{
    constexpr long double scale = 1000000000000.0L;
    const auto quantized =
        static_cast<std::int64_t>(std::llround(static_cast<long double>(value) * scale));
    hash_i64(hash, quantized);
}

template <std::floating_point T>
void hash_quantized_span(std::uint64_t& hash, std::span<const T> values)
{
    hash_i64(hash, static_cast<std::int64_t>(values.size()));
    for (const T value : values)
        hash_quantized<T>(hash, value);
}

inline std::int64_t signed_plan_hash(std::uint64_t hash)
{
    return static_cast<std::int64_t>(hash & 0x7fffffffffffffffULL);
}

inline std::uint64_t plan_hash_seed(int parent_cell_id,
                                    int local_cell_id,
                                    int chart_path,
                                    int refinement_depth)
{
    std::uint64_t hash = fnv_offset_basis;
    hash_i64(hash, parent_cell_id);
    hash_i64(hash, local_cell_id);
    hash_i64(hash, chart_path);
    hash_i64(hash, refinement_depth);
    return hash;
}

template <std::floating_point T>
struct LineClipResult
{
    bool valid = false;
    T r_lo = T(0);
    T r_hi = T(0);
    int lower_facet = -1;
    int upper_facet = -1;
};

template <std::floating_point T>
struct HeightDirectionFrame
{
    int tdim = 0;
    int base_dim = 0;
    std::array<T, 3> origin = {};
    std::array<T, 3> er = {};
    std::array<T, 6> base = {};
};

template <std::floating_point T>
struct BaseTriangle2
{
    std::array<T, 6> vertices = {};
    int depth = 0;
};

template <std::floating_point T>
struct ProjectedBaseSplitResult
{
    std::vector<BaseTriangle2<T>> cells;
    bool success = true;
};

template <std::floating_point T>
struct SimplexChart
{
    ChartKind kind = ChartKind::triangle_apex;
    int tdim = 0;
    int base_dim = 0;
    cell::type base_type = cell::type::interval;
    std::vector<T> r0_vertices;
    std::vector<T> r1_vertices;
    int r0_count = 0;
    int r1_count = 0;

    void edge_point(const std::vector<T>& verts,
                    int count,
                    T s,
                    std::span<T> out) const
    {
        if (count != 2)
            throw std::runtime_error("SimplexChart: expected edge entity");
        for (int d = 0; d < tdim; ++d)
        {
            out[static_cast<std::size_t>(d)] =
                (T(1) - s) * verts[static_cast<std::size_t>(d)]
              + s * verts[static_cast<std::size_t>(tdim + d)];
        }
    }

    void triangle_point(const std::vector<T>& verts,
                        int count,
                        std::span<const T> y,
                        std::span<T> out) const
    {
        if (count != 3)
            throw std::runtime_error("SimplexChart: expected triangle entity");
        const T s = y[0];
        const T t = y[1];
        for (int d = 0; d < tdim; ++d)
        {
            const T v0 = verts[static_cast<std::size_t>(d)];
            const T v1 = verts[static_cast<std::size_t>(tdim + d)];
            const T v2 = verts[static_cast<std::size_t>(2 * tdim + d)];
            out[static_cast<std::size_t>(d)] =
                v0 + s * (v1 - v0) + t * (v2 - v0);
        }
    }

    void entity_point(const std::vector<T>& verts,
                      int count,
                      std::span<const T> y,
                      std::span<T> out) const
    {
        if (count == 1)
        {
            std::copy_n(verts.data(), tdim, out.data());
        }
        else if (count == 2)
        {
            edge_point(verts, count, y[0], out);
        }
        else if (count == 3)
        {
            triangle_point(verts, count, y, out);
        }
        else
        {
            throw std::runtime_error("SimplexChart: unsupported entity size");
        }
    }

    void entity_derivative(const std::vector<T>& verts,
                           int count,
                           int axis,
                           std::span<T> out) const
    {
        std::fill(out.begin(), out.end(), T(0));
        if (count == 1)
            return;

        if (count == 2)
        {
            if (axis != 0)
                return;
            for (int d = 0; d < tdim; ++d)
            {
                out[static_cast<std::size_t>(d)] =
                    verts[static_cast<std::size_t>(tdim + d)]
                  - verts[static_cast<std::size_t>(d)];
            }
            return;
        }

        if (count == 3)
        {
            if (axis < 0 || axis > 1)
                return;
            const int vi = axis + 1;
            for (int d = 0; d < tdim; ++d)
            {
                out[static_cast<std::size_t>(d)] =
                    verts[static_cast<std::size_t>(vi * tdim + d)]
                  - verts[static_cast<std::size_t>(d)];
            }
            return;
        }

        throw std::runtime_error("SimplexChart: unsupported entity derivative");
    }

    void eval(T r, std::span<const T> y, std::span<T> out) const
    {
        std::array<T, 3> p0{};
        std::array<T, 3> p1{};
        if (kind == ChartKind::tetra_edge_edge_2_2)
        {
            edge_point(r0_vertices, r0_count, y[0],
                       std::span<T>(p0.data(), static_cast<std::size_t>(tdim)));
            edge_point(r1_vertices, r1_count, y[1],
                       std::span<T>(p1.data(), static_cast<std::size_t>(tdim)));
        }
        else
        {
            entity_point(r0_vertices, r0_count, y,
                         std::span<T>(p0.data(), static_cast<std::size_t>(tdim)));
            entity_point(r1_vertices, r1_count, y,
                         std::span<T>(p1.data(), static_cast<std::size_t>(tdim)));
        }

        for (int d = 0; d < tdim; ++d)
        {
            out[static_cast<std::size_t>(d)] =
                (T(1) - r) * p0[static_cast<std::size_t>(d)]
              + r * p1[static_cast<std::size_t>(d)];
        }
    }

    void derivatives(T r,
                     std::span<const T> y,
                     std::span<T> cr,
                     std::vector<T>& cy) const
    {
        std::array<T, 3> p0{};
        std::array<T, 3> p1{};
        cy.assign(static_cast<std::size_t>(base_dim * tdim), T(0));

        if (kind == ChartKind::tetra_edge_edge_2_2)
        {
            edge_point(r0_vertices, r0_count, y[0],
                       std::span<T>(p0.data(), static_cast<std::size_t>(tdim)));
            edge_point(r1_vertices, r1_count, y[1],
                       std::span<T>(p1.data(), static_cast<std::size_t>(tdim)));
            for (int d = 0; d < tdim; ++d)
            {
                cr[static_cast<std::size_t>(d)] =
                    p1[static_cast<std::size_t>(d)] - p0[static_cast<std::size_t>(d)];
                cy[static_cast<std::size_t>(d)] =
                    (T(1) - r) * (r0_vertices[static_cast<std::size_t>(tdim + d)]
                                  - r0_vertices[static_cast<std::size_t>(d)]);
                cy[static_cast<std::size_t>(tdim + d)] =
                    r * (r1_vertices[static_cast<std::size_t>(tdim + d)]
                         - r1_vertices[static_cast<std::size_t>(d)]);
            }
            return;
        }

        entity_point(r0_vertices, r0_count, y,
                     std::span<T>(p0.data(), static_cast<std::size_t>(tdim)));
        entity_point(r1_vertices, r1_count, y,
                     std::span<T>(p1.data(), static_cast<std::size_t>(tdim)));

        for (int d = 0; d < tdim; ++d)
        {
            cr[static_cast<std::size_t>(d)] =
                p1[static_cast<std::size_t>(d)] - p0[static_cast<std::size_t>(d)];
        }

        std::array<T, 3> deriv{};
        for (int a = 0; a < base_dim; ++a)
        {
            std::span<T> out(
                cy.data() + static_cast<std::size_t>(a * tdim),
                static_cast<std::size_t>(tdim));
            if (r0_count > 1)
            {
                entity_derivative(
                    r0_vertices, r0_count, a,
                    std::span<T>(deriv.data(), static_cast<std::size_t>(tdim)));
                for (int d = 0; d < tdim; ++d)
                    out[static_cast<std::size_t>(d)] =
                        (T(1) - r) * deriv[static_cast<std::size_t>(d)];
            }
            else if (r1_count > 1)
            {
                entity_derivative(
                    r1_vertices, r1_count, a,
                    std::span<T>(deriv.data(), static_cast<std::size_t>(tdim)));
                for (int d = 0; d < tdim; ++d)
                    out[static_cast<std::size_t>(d)] =
                        r * deriv[static_cast<std::size_t>(d)];
            }
        }
    }
};

template <std::floating_point T>
std::int64_t simplex_chart_plan_hash(const SimplexChart<T>& chart,
                                     int parent_cell_id,
                                     int local_cell_id,
                                     int refinement_depth)
{
    std::uint64_t hash =
        plan_hash_seed(parent_cell_id, local_cell_id, 1, refinement_depth);
    hash_i64(hash, static_cast<std::int64_t>(chart.kind));
    hash_i64(hash, static_cast<std::int64_t>(chart.tdim));
    hash_i64(hash, static_cast<std::int64_t>(chart.base_dim));
    hash_i64(hash, static_cast<std::int64_t>(chart.base_type));
    hash_i64(hash, static_cast<std::int64_t>(chart.r0_count));
    hash_i64(hash, static_cast<std::int64_t>(chart.r1_count));
    hash_quantized_span<T>(
        hash, std::span<const T>(chart.r0_vertices.data(), chart.r0_vertices.size()));
    hash_quantized_span<T>(
        hash, std::span<const T>(chart.r1_vertices.data(), chart.r1_vertices.size()));
    return signed_plan_hash(hash);
}

template <std::floating_point T>
std::int64_t straight_fallback_plan_hash(std::span<const T> roots_ref,
                                         int tdim,
                                         int parent_cell_id,
                                         int local_cell_id,
                                         int refinement_depth)
{
    std::uint64_t hash =
        plan_hash_seed(parent_cell_id, local_cell_id, 3, refinement_depth);
    hash_i64(hash, static_cast<std::int64_t>(tdim));
    hash_quantized_span<T>(hash, roots_ref);
    return signed_plan_hash(hash);
}

template <std::floating_point T>
std::int64_t height_frame_plan_component(const HeightDirectionFrame<T>& frame)
{
    std::uint64_t hash = fnv_offset_basis;
    hash_i64(hash, static_cast<std::int64_t>(ChartKind::projected_height_direction));
    hash_i64(hash, static_cast<std::int64_t>(frame.tdim));
    hash_i64(hash, static_cast<std::int64_t>(frame.base_dim));
    hash_quantized_span<T>(
        hash, std::span<const T>(frame.origin.data(), static_cast<std::size_t>(frame.tdim)));
    hash_quantized_span<T>(
        hash, std::span<const T>(frame.er.data(), static_cast<std::size_t>(frame.tdim)));
    hash_quantized_span<T>(
        hash, std::span<const T>(frame.base.data(),
                                 static_cast<std::size_t>(frame.base_dim * frame.tdim)));
    return signed_plan_hash(hash);
}

template <std::floating_point T>
std::int64_t triangle_height_plan_component(const HeightDirectionFrame<T>& frame,
                                           std::span<const T> split_points)
{
    std::uint64_t hash = static_cast<std::uint64_t>(
        height_frame_plan_component<T>(frame));
    hash_quantized_span<T>(hash, split_points);
    return signed_plan_hash(hash);
}

template <std::floating_point T>
std::int64_t tetra_height_plan_component(
    const HeightDirectionFrame<T>& frame,
    const std::vector<BaseTriangle2<T>>& base_cells)
{
    std::uint64_t hash = static_cast<std::uint64_t>(
        height_frame_plan_component<T>(frame));
    hash_i64(hash, static_cast<std::int64_t>(base_cells.size()));
    int max_depth = 0;
    for (const auto& cell : base_cells)
        max_depth = std::max(max_depth, cell.depth);
    hash_i64(hash, static_cast<std::int64_t>(max_depth));
    for (const auto& cell : base_cells)
    {
        hash_i64(hash, static_cast<std::int64_t>(cell.depth));
        hash_quantized_span<T>(
            hash, std::span<const T>(cell.vertices.data(), cell.vertices.size()));
    }
    return signed_plan_hash(hash);
}

inline std::int64_t accepted_height_plan_hash(std::int64_t component_hash,
                                              int parent_cell_id,
                                              int local_cell_id,
                                              int refinement_depth)
{
    std::uint64_t hash =
        plan_hash_seed(parent_cell_id, local_cell_id, 2, refinement_depth);
    hash_i64(hash, component_hash);
    return signed_plan_hash(hash);
}

template <std::floating_point T>
T dot(std::span<const T> a, std::span<const T> b)
{
    T value = T(0);
    for (std::size_t i = 0; i < a.size(); ++i)
        value += a[i] * b[i];
    return value;
}

template <std::floating_point T>
T norm(std::span<const T> a)
{
    return std::sqrt(dot<T>(a, a));
}

template <std::floating_point T>
std::array<T, 3> cross3(const std::array<T, 3>& a, const std::array<T, 3>& b)
{
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]};
}

template <std::floating_point T>
T dot3(const std::array<T, 3>& a, const std::array<T, 3>& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <std::floating_point T>
T norm3(const std::array<T, 3>& a)
{
    return std::sqrt(dot3<T>(a, a));
}

template <std::floating_point T>
bool normalize3(std::array<T, 3>& a)
{
    const T n = norm3<T>(a);
    if (n <= std::numeric_limits<T>::epsilon())
        return false;
    for (auto& value : a)
        value /= n;
    return true;
}

template <std::floating_point T>
T gram_measure(std::span<const T> tangents, int tangent_count, int gdim)
{
    if (tangent_count == 0)
        return T(1);

    if (tangent_count == 1)
        return norm<T>(tangents);

    T gram[9] = {};
    for (int i = 0; i < tangent_count; ++i)
    {
        for (int j = 0; j < tangent_count; ++j)
        {
            T value = T(0);
            for (int d = 0; d < gdim; ++d)
            {
                value += tangents[static_cast<std::size_t>(i * gdim + d)]
                       * tangents[static_cast<std::size_t>(j * gdim + d)];
            }
            gram[static_cast<std::size_t>(i * tangent_count + j)] = value;
        }
    }

    if (tangent_count == 2)
    {
        const T det = gram[0] * gram[3] - gram[1] * gram[2];
        return std::sqrt(std::max(det, T(0)));
    }

    const T det =
        gram[0] * (gram[4] * gram[8] - gram[5] * gram[7])
      - gram[1] * (gram[3] * gram[8] - gram[5] * gram[6])
      + gram[2] * (gram[3] * gram[7] - gram[4] * gram[6]);
    return std::sqrt(std::max(det, T(0)));
}

template <std::floating_point T, std::integral I>
std::vector<T> parent_cell_vertex_coords_vtk(const MeshView<T, I>& mesh,
                                             I cell_id)
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
void append_rule_points(quadrature::QuadratureRules<T>& rules,
                        std::span<const T> points,
                        std::span<const T> weights,
                        int tdim,
                        int parent_cell_id)
{
    if (rules._tdim == 0)
        rules._tdim = tdim;
    if (rules._offset.empty())
        rules._offset.push_back(0);
    if (rules._tdim != tdim)
        throw std::runtime_error("implicit quadrature: inconsistent point dimension");

    rules._points.insert(rules._points.end(), points.begin(), points.end());
    rules._weights.insert(rules._weights.end(), weights.begin(), weights.end());
    rules._parent_map.push_back(static_cast<std::int32_t>(parent_cell_id));
    rules._offset.push_back(static_cast<std::int32_t>(rules._weights.size()));
    rules._debug_local_cell_id.push_back(std::int32_t(-1));
    rules._debug_chart_path.push_back(std::int32_t(0));
    rules._debug_refinement_depth.push_back(std::int32_t(0));
    rules._debug_chart_plan_hash.push_back(std::int64_t(0));
    rules._debug_candidate_mask.push_back(std::int32_t(0));
    rules._debug_rejection_reason.push_back(std::int32_t(0));
    rules._debug_measure_probe.push_back(T(0));
    rules._debug_validation_weight_sum.push_back(T(0));
}

template <std::floating_point T>
void annotate_rule_debug(quadrature::QuadratureRules<T>& rules,
                         int local_cell_id,
                         int chart_path,
                         int refinement_depth,
                         std::int64_t chart_plan_hash,
                         std::int32_t candidate_mask,
                         std::int32_t rejection_reason,
                         T measure_probe,
                         T validation_weight_sum)
{
    const std::size_t n = rules._parent_map.size();
    rules._debug_local_cell_id.assign(n, static_cast<std::int32_t>(local_cell_id));
    rules._debug_chart_path.assign(n, static_cast<std::int32_t>(chart_path));
    rules._debug_refinement_depth.assign(
        n, static_cast<std::int32_t>(refinement_depth));
    rules._debug_chart_plan_hash.assign(n, chart_plan_hash);
    rules._debug_candidate_mask.assign(n, candidate_mask);
    rules._debug_rejection_reason.assign(n, rejection_reason);
    rules._debug_measure_probe.assign(n, measure_probe);
    rules._debug_validation_weight_sum.assign(n, validation_weight_sum);
}

template <std::floating_point T>
void map_canonical_simplex_to_vertices(const quadrature::ReferenceQuadratureRule<T>& ref_rule,
                                       std::span<const T> vertices,
                                       int vertex_dim,
                                       std::vector<T>& out_points)
{
    const int simplex_dim = ref_rule._tdim;
    const int nq = ref_rule._num_points;
    out_points.assign(static_cast<std::size_t>(nq * vertex_dim), T(0));
    const T* v0 = vertices.data();

    for (int q = 0; q < nq; ++q)
    {
        const T* y = ref_rule._points.data() + static_cast<std::size_t>(q * simplex_dim);
        T* x = out_points.data() + static_cast<std::size_t>(q * vertex_dim);

        for (int d = 0; d < vertex_dim; ++d)
            x[d] = v0[d];

        for (int a = 1; a <= simplex_dim; ++a)
        {
            const T* va = vertices.data() + static_cast<std::size_t>(a * vertex_dim);
            for (int d = 0; d < vertex_dim; ++d)
                x[d] += y[a - 1] * (va[d] - v0[d]);
        }
    }
}

template <std::floating_point T>
void append_simplex_entity_quadrature(quadrature::QuadratureRules<T>& rules,
                                      cell::type simplex_type,
                                      std::span<const T> ref_vertices,
                                      std::span<const T> physical_vertices,
                                      int parent_tdim,
                                      int gdim,
                                      int parent_cell_id,
                                      int order)
{
    const auto ref_rule = quadrature::get_reference_rule<T>(simplex_type, order);
    std::vector<T> points;
    map_canonical_simplex_to_vertices(ref_rule, ref_vertices, parent_tdim, points);

    std::vector<T> weights;
    weights.reserve(static_cast<std::size_t>(ref_rule._num_points));
    std::vector<T> physical_tangents(
        static_cast<std::size_t>(ref_rule._tdim * gdim), T(0));
    const T* v0 = physical_vertices.data();
    for (int a = 0; a < ref_rule._tdim; ++a)
    {
        const T* va = physical_vertices.data()
                    + static_cast<std::size_t>((a + 1) * gdim);
        for (int d = 0; d < gdim; ++d)
        {
            physical_tangents[static_cast<std::size_t>(a * gdim + d)] =
                va[d] - v0[d];
        }
    }
    const T measure = gram_measure<T>(
        std::span<const T>(physical_tangents.data(), physical_tangents.size()),
        ref_rule._tdim,
        gdim);

    for (int q = 0; q < ref_rule._num_points; ++q)
        weights.push_back(ref_rule._weights[static_cast<std::size_t>(q)] * measure);

    append_rule_points(
        rules,
        std::span<const T>(points.data(), points.size()),
        std::span<const T>(weights.data(), weights.size()),
        parent_tdim,
        parent_cell_id);
}

template <std::floating_point T>
void append_uncut_parent_cell(quadrature::QuadratureRules<T>& rules,
                              cell::type ctype,
                              std::span<const T> physical_vertices,
                              int parent_tdim,
                              int gdim,
                              int parent_cell_id,
                              int order)
{
    const auto ref_rule = quadrature::get_reference_rule<T>(ctype, order);
    std::vector<T> weights;
    weights.reserve(static_cast<std::size_t>(ref_rule._num_points));

    const T measure = cell::affine_volume_factor<T>(
        ctype, physical_vertices.data(), gdim);
    for (int q = 0; q < ref_rule._num_points; ++q)
        weights.push_back(ref_rule._weights[static_cast<std::size_t>(q)] * measure);

    append_rule_points(
        rules,
        std::span<const T>(ref_rule._points.data(), ref_rule._points.size()),
        std::span<const T>(weights.data(), weights.size()),
        parent_tdim,
        parent_cell_id);
}

template <std::floating_point T>
void gather_entity_vertices(const AdaptCell<T>& ac,
                            std::span<const std::int32_t> cell_vertices,
                            std::span<const int> local_ids,
                            std::vector<T>& out)
{
    out.assign(static_cast<std::size_t>(local_ids.size() * ac.tdim), T(0));
    for (std::size_t i = 0; i < local_ids.size(); ++i)
    {
        const int gv = static_cast<int>(
            cell_vertices[static_cast<std::size_t>(local_ids[i])]);
        for (int d = 0; d < ac.tdim; ++d)
        {
            out[static_cast<std::size_t>(i * ac.tdim + d)] =
                ac.vertex_coords[static_cast<std::size_t>(gv * ac.tdim + d)];
        }
    }
}

template <std::floating_point T>
void gather_leaf_ref_vertices(const AdaptCell<T>& ac,
                              int cell_id,
                              std::vector<T>& out)
{
    const auto verts = ac.entity_to_vertex[ac.tdim][static_cast<std::int32_t>(cell_id)];
    out.assign(static_cast<std::size_t>(verts.size() * ac.tdim), T(0));
    for (std::size_t i = 0; i < verts.size(); ++i)
    {
        const int gv = static_cast<int>(verts[i]);
        for (int d = 0; d < ac.tdim; ++d)
        {
            out[static_cast<std::size_t>(i * ac.tdim + d)] =
                ac.vertex_coords[static_cast<std::size_t>(gv * ac.tdim + d)];
        }
    }
}

template <std::floating_point T>
int vertex_sign(const AdaptCell<T>& ac, int vertex_id, int level_set_id)
{
    const std::uint64_t bit = std::uint64_t(1) << level_set_id;
    const auto zero = ac.zero_mask_per_vertex[static_cast<std::size_t>(vertex_id)] & bit;
    if (zero != 0)
        return 0;
    const auto neg = ac.negative_mask_per_vertex[static_cast<std::size_t>(vertex_id)] & bit;
    return neg != 0 ? -1 : 1;
}

template <std::floating_point T>
bool leaf_has_uniform_sign(const AdaptCell<T>& ac,
                           std::span<const std::int32_t> cell_vertices,
                           int level_set_id,
                           int wanted_sign)
{
    for (const auto gv : cell_vertices)
    {
        if (vertex_sign(ac, static_cast<int>(gv), level_set_id) != wanted_sign)
            return false;
    }
    return true;
}

template <std::floating_point T>
SimplexChart<T> make_chart_from_leaf(const AdaptCell<T>& ac,
                                     int cell_id,
                                     int level_set_id)
{
    const auto ctype = ac.entity_types[ac.tdim][static_cast<std::size_t>(cell_id)];
    if (ctype != cell::type::triangle && ctype != cell::type::tetrahedron)
        throw std::runtime_error(
            "implicit quadrature supports only triangle/tetrahedron simplex leaves in v1");

    const auto verts = ac.entity_to_vertex[ac.tdim][static_cast<std::int32_t>(cell_id)];
    std::vector<int> negative_local;
    std::vector<int> positive_local;
    for (int i = 0; i < static_cast<int>(verts.size()); ++i)
    {
        const int sign = vertex_sign(ac, static_cast<int>(verts[static_cast<std::size_t>(i)]),
                                     level_set_id);
        if (sign < 0)
            negative_local.push_back(i);
        else if (sign > 0)
            positive_local.push_back(i);
        else
            throw std::runtime_error(
                "implicit quadrature v1 does not accept zero vertices in chart leaves");
    }

    SimplexChart<T> chart;
    chart.tdim = ac.tdim;

    if (ctype == cell::type::triangle)
    {
        if (!((negative_local.size() == 1 && positive_local.size() == 2)
              || (negative_local.size() == 2 && positive_local.size() == 1)))
        {
            throw std::runtime_error("implicit quadrature: unsupported triangle sign split");
        }

        chart.kind = ChartKind::triangle_apex;
        chart.base_dim = 1;
        chart.base_type = cell::type::interval;
        chart.r0_count = static_cast<int>(negative_local.size());
        chart.r1_count = static_cast<int>(positive_local.size());
        gather_entity_vertices(
            ac, verts,
            std::span<const int>(negative_local.data(), negative_local.size()),
            chart.r0_vertices);
        gather_entity_vertices(
            ac, verts,
            std::span<const int>(positive_local.data(), positive_local.size()),
            chart.r1_vertices);
        return chart;
    }

    if (negative_local.size() == 2 && positive_local.size() == 2)
    {
        chart.kind = ChartKind::tetra_edge_edge_2_2;
        chart.base_dim = 2;
        chart.base_type = cell::type::quadrilateral;
    }
    else if ((negative_local.size() == 1 && positive_local.size() == 3)
             || (negative_local.size() == 3 && positive_local.size() == 1))
    {
        chart.kind = ChartKind::tetra_apex_1_3;
        chart.base_dim = 2;
        chart.base_type = cell::type::triangle;
    }
    else
    {
        throw std::runtime_error("implicit quadrature: unsupported tetrahedron sign split");
    }

    chart.r0_count = static_cast<int>(negative_local.size());
    chart.r1_count = static_cast<int>(positive_local.size());
    gather_entity_vertices(
        ac, verts,
        std::span<const int>(negative_local.data(), negative_local.size()),
        chart.r0_vertices);
    gather_entity_vertices(
        ac, verts,
        std::span<const int>(positive_local.data(), positive_local.size()),
        chart.r1_vertices);
    return chart;
}

template <std::floating_point T>
std::vector<T> push_forward_points(cell::type parent_type,
                                   std::span<const T> parent_vertices,
                                   int gdim,
                                   std::span<const T> ref_points)
{
    std::vector<T> parent_storage(parent_vertices.begin(), parent_vertices.end());
    return cell::push_forward_affine_map<T>(
        parent_type, parent_storage, gdim, ref_points);
}

template <std::floating_point T>
void push_forward_tangent(cell::type parent_type,
                          std::span<const T> parent_vertices,
                          int gdim,
                          std::span<const T> xi,
                          std::span<const T> tangent_ref,
                          std::span<T> tangent_phys)
{
    const int tdim = static_cast<int>(xi.size());
    std::vector<T> pts(static_cast<std::size_t>(2 * tdim), T(0));
    for (int d = 0; d < tdim; ++d)
    {
        pts[static_cast<std::size_t>(d)] = xi[static_cast<std::size_t>(d)];
        pts[static_cast<std::size_t>(tdim + d)] =
            xi[static_cast<std::size_t>(d)] + tangent_ref[static_cast<std::size_t>(d)];
    }

    const auto phys = push_forward_points<T>(
        parent_type, parent_vertices, gdim,
        std::span<const T>(pts.data(), pts.size()));
    for (int d = 0; d < gdim; ++d)
    {
        tangent_phys[static_cast<std::size_t>(d)] =
            phys[static_cast<std::size_t>(gdim + d)]
          - phys[static_cast<std::size_t>(d)];
    }
}

template <std::floating_point T>
T physical_measure_from_ref_tangents(cell::type parent_type,
                                     std::span<const T> parent_vertices,
                                     int gdim,
                                     std::span<const T> xi,
                                     std::span<const T> tangents_ref,
                                     int tangent_count)
{
    const int tdim = static_cast<int>(xi.size());
    std::vector<T> tangents_phys(
        static_cast<std::size_t>(tangent_count * gdim), T(0));
    for (int a = 0; a < tangent_count; ++a)
    {
        push_forward_tangent<T>(
            parent_type,
            parent_vertices,
            gdim,
            xi,
            std::span<const T>(
                tangents_ref.data() + static_cast<std::size_t>(a * tdim),
                static_cast<std::size_t>(tdim)),
            std::span<T>(
                tangents_phys.data() + static_cast<std::size_t>(a * gdim),
                static_cast<std::size_t>(gdim)));
    }
    return gram_measure<T>(
        std::span<const T>(tangents_phys.data(), tangents_phys.size()),
        tangent_count,
        gdim);
}

template <std::floating_point T>
T polygon_measure_from_physical_points(std::span<const T> points, int gdim)
{
    const int n = static_cast<int>(points.size()) / gdim;
    if (n < 2)
        return T(0);

    if (n == 2)
    {
        std::vector<T> tangent(static_cast<std::size_t>(gdim), T(0));
        for (int d = 0; d < gdim; ++d)
            tangent[static_cast<std::size_t>(d)] =
                points[static_cast<std::size_t>(gdim + d)] - points[static_cast<std::size_t>(d)];
        return norm<T>(std::span<const T>(tangent.data(), tangent.size()));
    }

    std::vector<int> order(static_cast<std::size_t>(n));
    std::iota(order.begin(), order.end(), 0);
    if (n > 3)
    {
        std::vector<T> centroid(static_cast<std::size_t>(gdim), T(0));
        for (int i = 0; i < n; ++i)
            for (int d = 0; d < gdim; ++d)
                centroid[static_cast<std::size_t>(d)]
                    += points[static_cast<std::size_t>(i * gdim + d)];
        for (auto& value : centroid)
            value /= static_cast<T>(n);

        if (gdim == 3)
        {
            std::array<T, 3> e0{};
            std::array<T, 3> e1{};
            for (int d = 0; d < 3; ++d)
            {
                e0[static_cast<std::size_t>(d)] =
                    points[static_cast<std::size_t>(gdim + d)] - points[static_cast<std::size_t>(d)];
                e1[static_cast<std::size_t>(d)] =
                    points[static_cast<std::size_t>(2 * gdim + d)] - points[static_cast<std::size_t>(d)];
            }
            const auto normal = cross3<T>(e0, e1);
            int drop = 0;
            if (std::abs(normal[1]) > std::abs(normal[drop]))
                drop = 1;
            if (std::abs(normal[2]) > std::abs(normal[drop]))
                drop = 2;
            const int a = (drop + 1) % 3;
            const int b = (drop + 2) % 3;
            std::sort(
                order.begin(),
                order.end(),
                [&](int lhs, int rhs)
                {
                    const T lx = points[static_cast<std::size_t>(lhs * gdim + a)]
                               - centroid[static_cast<std::size_t>(a)];
                    const T ly = points[static_cast<std::size_t>(lhs * gdim + b)]
                               - centroid[static_cast<std::size_t>(b)];
                    const T rx = points[static_cast<std::size_t>(rhs * gdim + a)]
                               - centroid[static_cast<std::size_t>(a)];
                    const T ry = points[static_cast<std::size_t>(rhs * gdim + b)]
                               - centroid[static_cast<std::size_t>(b)];
                    return std::atan2(ly, lx) < std::atan2(ry, rx);
                });
        }
    }

    T measure = T(0);
    for (int k = 1; k + 1 < n; ++k)
    {
        std::vector<T> tangents(static_cast<std::size_t>(2 * gdim), T(0));
        const int i0 = order[0];
        const int i1 = order[static_cast<std::size_t>(k)];
        const int i2 = order[static_cast<std::size_t>(k + 1)];
        for (int d = 0; d < gdim; ++d)
        {
            tangents[static_cast<std::size_t>(d)] =
                points[static_cast<std::size_t>(i1 * gdim + d)]
              - points[static_cast<std::size_t>(i0 * gdim + d)];
            tangents[static_cast<std::size_t>(gdim + d)] =
                points[static_cast<std::size_t>(i2 * gdim + d)]
              - points[static_cast<std::size_t>(i0 * gdim + d)];
        }
        measure += T(0.5) * gram_measure<T>(
            std::span<const T>(tangents.data(), tangents.size()), 2, gdim);
    }
    return measure;
}

template <std::floating_point T, std::integral I>
T linearized_interface_measure_for_leaf(const AdaptCell<T>& ac,
                                        int cell_id,
                                        const LevelSetCell<T, I>& ls_cell,
                                        int level_set_id,
                                        const ImplicitQuadratureOptions& opts)
{
    const auto verts = ac.entity_to_vertex[ac.tdim][static_cast<std::int32_t>(cell_id)];
    const auto ctype = ac.entity_types[ac.tdim][static_cast<std::size_t>(cell_id)];
    const auto edges = cell::edges(ctype);
    std::vector<T> roots_ref;
    roots_ref.reserve(static_cast<std::size_t>(4 * ac.tdim));

    for (const auto& edge : edges)
    {
        const int a = static_cast<int>(edge[0]);
        const int b = static_cast<int>(edge[1]);
        const int va = static_cast<int>(verts[static_cast<std::size_t>(a)]);
        const int vb = static_cast<int>(verts[static_cast<std::size_t>(b)]);
        const int sa = vertex_sign(ac, va, level_set_id);
        const int sb = vertex_sign(ac, vb, level_set_id);
        if (sa == 0 || sb == 0 || sa == sb)
            continue;

        std::array<T, 3> xa{};
        std::array<T, 3> xb{};
        for (int d = 0; d < ac.tdim; ++d)
        {
            xa[static_cast<std::size_t>(d)] =
                ac.vertex_coords[static_cast<std::size_t>(va * ac.tdim + d)];
            xb[static_cast<std::size_t>(d)] =
                ac.vertex_coords[static_cast<std::size_t>(vb * ac.tdim + d)];
        }
        const T fa = ls_cell.value(
            std::span<const T>(xa.data(), static_cast<std::size_t>(ac.tdim)));
        const T fb = ls_cell.value(
            std::span<const T>(xb.data(), static_cast<std::size_t>(ac.tdim)));
        if (std::abs(fb - fa) <= static_cast<T>(opts.root_tol))
            continue;
        const T t = std::clamp(-fa / (fb - fa), T(0), T(1));
        for (int d = 0; d < ac.tdim; ++d)
        {
            roots_ref.push_back(
                xa[static_cast<std::size_t>(d)]
              + t * (xb[static_cast<std::size_t>(d)] - xa[static_cast<std::size_t>(d)]));
        }
    }

    if (static_cast<int>(roots_ref.size()) / ac.tdim < ac.tdim)
        return T(0);
    std::vector<T> unique_roots;
    const T unique_tol = static_cast<T>(100) * static_cast<T>(opts.root_tol);
    const int n_roots = static_cast<int>(roots_ref.size()) / ac.tdim;
    for (int i = 0; i < n_roots; ++i)
    {
        bool duplicate = false;
        for (int j = 0; j < static_cast<int>(unique_roots.size()) / ac.tdim; ++j)
        {
            T distance2 = T(0);
            for (int d = 0; d < ac.tdim; ++d)
            {
                const T delta =
                    roots_ref[static_cast<std::size_t>(i * ac.tdim + d)]
                  - unique_roots[static_cast<std::size_t>(j * ac.tdim + d)];
                distance2 += delta * delta;
            }
            if (distance2 <= unique_tol * unique_tol)
            {
                duplicate = true;
                break;
            }
        }
        if (!duplicate)
        {
            for (int d = 0; d < ac.tdim; ++d)
                unique_roots.push_back(roots_ref[static_cast<std::size_t>(i * ac.tdim + d)]);
        }
    }
    roots_ref = std::move(unique_roots);
    const auto phys = push_forward_points<T>(
        ls_cell.cell_type,
        std::span<const T>(
            ls_cell.parent_vertex_coords.data(),
            ls_cell.parent_vertex_coords.size()),
        ls_cell.gdim,
        std::span<const T>(roots_ref.data(), roots_ref.size()));
    return polygon_measure_from_physical_points<T>(
        std::span<const T>(phys.data(), phys.size()), ls_cell.gdim);
}

template <std::floating_point T, std::integral I>
bool surface_rule_matches_linearized_leaf(
    const quadrature::QuadratureRules<T>& leaf_rules,
    const AdaptCell<T>& ac,
    int cell_id,
    const LevelSetCell<T, I>& ls_cell,
    int level_set_id,
    IntegrationDomain domain,
    const ImplicitQuadratureOptions& opts,
    T min_ratio = T(0.02),
    T max_ratio = T(50))
{
    if (domain != IntegrationDomain::interface_surface)
        return true;
    const T weight_sum = std::accumulate(
        leaf_rules._weights.begin(), leaf_rules._weights.end(), T(0));
    const T linear_measure = linearized_interface_measure_for_leaf<T, I>(
        ac, cell_id, ls_cell, level_set_id, opts);
    const T scale = std::max({T(1), std::abs(weight_sum), std::abs(linear_measure)});
    const T tiny = static_cast<T>(1000) * static_cast<T>(opts.root_tol) * scale;
    if (linear_measure <= tiny)
        return weight_sum <= tiny;
    const T ratio = weight_sum / linear_measure;
    return ratio >= min_ratio && ratio <= max_ratio;
}

template <std::floating_point T>
void append_unique_ref_point(std::vector<T>& points,
                             std::span<const T> candidate,
                             int tdim,
                             T tol)
{
    const int n = static_cast<int>(points.size()) / tdim;
    for (int i = 0; i < n; ++i)
    {
        T distance2 = T(0);
        for (int d = 0; d < tdim; ++d)
        {
            const T delta =
                points[static_cast<std::size_t>(i * tdim + d)]
              - candidate[static_cast<std::size_t>(d)];
            distance2 += delta * delta;
        }
        if (distance2 <= tol * tol)
            return;
    }
    for (int d = 0; d < tdim; ++d)
        points.push_back(candidate[static_cast<std::size_t>(d)]);
}

template <std::floating_point T>
void sort_coplanar_polygon_points(std::vector<T>& points, int tdim)
{
    const int n = static_cast<int>(points.size()) / tdim;
    if (tdim != 3 || n <= 3)
        return;

    std::array<T, 3> centroid{};
    for (int i = 0; i < n; ++i)
    {
        for (int d = 0; d < 3; ++d)
            centroid[static_cast<std::size_t>(d)] +=
                points[static_cast<std::size_t>(i * 3 + d)];
    }
    for (int d = 0; d < 3; ++d)
        centroid[static_cast<std::size_t>(d)] /= static_cast<T>(n);

    std::array<T, 3> normal{};
    bool have_normal = false;
    for (int i = 1; i + 1 < n && !have_normal; ++i)
    {
        std::array<T, 3> a{};
        std::array<T, 3> b{};
        for (int d = 0; d < 3; ++d)
        {
            a[static_cast<std::size_t>(d)] =
                points[static_cast<std::size_t>(i * 3 + d)] - points[static_cast<std::size_t>(d)];
            b[static_cast<std::size_t>(d)] =
                points[static_cast<std::size_t>((i + 1) * 3 + d)]
              - points[static_cast<std::size_t>(d)];
        }
        normal = cross3<T>(a, b);
        have_normal = normalize3<T>(normal);
    }
    if (!have_normal)
        return;

    int drop = 0;
    if (std::abs(normal[1]) > std::abs(normal[drop]))
        drop = 1;
    if (std::abs(normal[2]) > std::abs(normal[drop]))
        drop = 2;
    const int a = (drop + 1) % 3;
    const int b = (drop + 2) % 3;

    std::vector<int> order(static_cast<std::size_t>(n));
    std::iota(order.begin(), order.end(), 0);
    std::sort(
        order.begin(),
        order.end(),
        [&](int lhs, int rhs)
        {
            const T lx = points[static_cast<std::size_t>(lhs * 3 + a)] - centroid[a];
            const T ly = points[static_cast<std::size_t>(lhs * 3 + b)] - centroid[b];
            const T rx = points[static_cast<std::size_t>(rhs * 3 + a)] - centroid[a];
            const T ry = points[static_cast<std::size_t>(rhs * 3 + b)] - centroid[b];
            return std::atan2(ly, lx) < std::atan2(ry, rx);
        });

    std::vector<T> sorted;
    sorted.reserve(points.size());
    for (const int index : order)
    {
        for (int d = 0; d < 3; ++d)
            sorted.push_back(points[static_cast<std::size_t>(index * 3 + d)]);
    }
    points = std::move(sorted);
}

template <std::floating_point T, std::integral I>
bool append_straight_interface_fallback_quadrature(
    quadrature::QuadratureRules<T>& rules,
    const AdaptCell<T>& ac,
    int cell_id,
    const LevelSetCell<T, I>& ls_cell,
    const ImplicitQuadratureOptions& opts,
    int parent_cell_id,
    std::vector<T>& roots_ref)
{
    roots_ref.clear();
    const auto verts = ac.entity_to_vertex[ac.tdim][static_cast<std::int32_t>(cell_id)];
    const auto ctype = ac.entity_types[ac.tdim][static_cast<std::size_t>(cell_id)];
    const auto edges = cell::edges(ctype);
    const T point_tol = static_cast<T>(100) * static_cast<T>(opts.root_tol);

    for (const auto& edge : edges)
    {
        const int va = static_cast<int>(verts[static_cast<std::size_t>(edge[0])]);
        const int vb = static_cast<int>(verts[static_cast<std::size_t>(edge[1])]);
        std::array<T, 3> xa{};
        std::array<T, 3> xb{};
        for (int d = 0; d < ac.tdim; ++d)
        {
            xa[static_cast<std::size_t>(d)] =
                ac.vertex_coords[static_cast<std::size_t>(va * ac.tdim + d)];
            xb[static_cast<std::size_t>(d)] =
                ac.vertex_coords[static_cast<std::size_t>(vb * ac.tdim + d)];
        }
        const T fa = ls_cell.value(
            std::span<const T>(xa.data(), static_cast<std::size_t>(ac.tdim)));
        const T fb = ls_cell.value(
            std::span<const T>(xb.data(), static_cast<std::size_t>(ac.tdim)));

        std::array<T, 3> root{};
        bool have_root = false;
        if (std::abs(fa) <= static_cast<T>(opts.root_tol))
        {
            root = xa;
            have_root = true;
        }
        else if (std::abs(fb) <= static_cast<T>(opts.root_tol))
        {
            root = xb;
            have_root = true;
        }
        else if (fa * fb < T(0))
        {
            auto phi = [&](std::span<const T> xi) -> T
            {
                return ls_cell.value(xi);
            };
            const auto info = cell::edge_root::find_root_parameter_info<T>(
                std::span<const T>(xa.data(), static_cast<std::size_t>(ac.tdim)),
                std::span<const T>(xb.data(), static_cast<std::size_t>(ac.tdim)),
                phi,
                opts.root_method,
                T(0),
                /*max_iter=*/64,
                static_cast<T>(opts.root_tol),
                static_cast<T>(opts.root_tol));
            if (info.converged
                && info.residual <= static_cast<T>(T(10) * opts.root_tol))
            {
                const T t = std::clamp(info.t, T(0), T(1));
                for (int d = 0; d < ac.tdim; ++d)
                {
                    root[static_cast<std::size_t>(d)] =
                        xa[static_cast<std::size_t>(d)]
                      + t * (xb[static_cast<std::size_t>(d)]
                             - xa[static_cast<std::size_t>(d)]);
                }
                have_root = true;
            }
        }

        if (have_root)
        {
            append_unique_ref_point<T>(
                roots_ref,
                std::span<const T>(root.data(), static_cast<std::size_t>(ac.tdim)),
                ac.tdim,
                point_tol);
        }
    }

    const int n_roots = static_cast<int>(roots_ref.size()) / ac.tdim;
    if (n_roots < ac.tdim)
        return false;
    sort_coplanar_polygon_points<T>(roots_ref, ac.tdim);

    auto append_surface_simplex = [&](std::span<const T> simplex_ref,
                                      cell::type simplex_type)
    {
        const auto physical_vertices = push_forward_points<T>(
            ls_cell.cell_type,
            std::span<const T>(
                ls_cell.parent_vertex_coords.data(),
                ls_cell.parent_vertex_coords.size()),
            ls_cell.gdim,
            simplex_ref);
        append_simplex_entity_quadrature<T>(
            rules,
            simplex_type,
            simplex_ref,
            std::span<const T>(physical_vertices.data(), physical_vertices.size()),
            ac.tdim,
            ls_cell.gdim,
            parent_cell_id,
            opts.order);
    };

    if (ac.tdim == 2)
    {
        std::vector<T> segment_ref{
            roots_ref.begin(),
            roots_ref.begin() + static_cast<std::ptrdiff_t>(2 * ac.tdim)};
        append_surface_simplex(
            std::span<const T>(segment_ref.data(), segment_ref.size()),
            cell::type::interval);
        return true;
    }

    if (ac.tdim == 3)
    {
        bool appended = false;
        for (int k = 1; k + 1 < n_roots; ++k)
        {
            std::vector<T> tri_ref;
            tri_ref.reserve(9);
            for (const int index : {0, k, k + 1})
            {
                for (int d = 0; d < 3; ++d)
                {
                    tri_ref.push_back(
                        roots_ref[static_cast<std::size_t>(index * 3 + d)]);
                }
            }
            append_surface_simplex(
                std::span<const T>(tri_ref.data(), tri_ref.size()),
                cell::type::triangle);
            appended = true;
        }
        return appended;
    }
    return false;
}

template <std::floating_point T, std::integral I>
struct RootOnChartLine
{
    bool valid = false;
    int failure_code = 0;
    T r = T(0);
    T den = T(0);
    std::array<T, 3> xi = {};
    std::array<T, 3> grad = {};
};

inline const char* chart_root_failure_message(int code)
{
    switch (code)
    {
    case 1:
        return "endpoint_orientation";
    case 2:
        return "root_structure";
    case 3:
        return "root_solve";
    case 4:
        return "transversality";
    default:
        return "unknown";
    }
}

template <std::floating_point T, std::integral I>
RootOnChartLine<T, I> solve_chart_root(const SimplexChart<T>& chart,
                                       const LevelSetCell<T, I>& ls_cell,
                                       std::span<const T> y,
                                       const ImplicitQuadratureOptions& opts)
{
    RootOnChartLine<T, I> result;
    std::array<T, 3> x0{};
    std::array<T, 3> x1{};
    chart.eval(T(0), y, std::span<T>(x0.data(), static_cast<std::size_t>(chart.tdim)));
    chart.eval(T(1), y, std::span<T>(x1.data(), static_cast<std::size_t>(chart.tdim)));

    const T f0 = ls_cell.value(
        std::span<const T>(x0.data(), static_cast<std::size_t>(chart.tdim)));
    const T f1 = ls_cell.value(
        std::span<const T>(x1.data(), static_cast<std::size_t>(chart.tdim)));
    if (f0 > static_cast<T>(opts.sign_tol) || f1 < -static_cast<T>(opts.sign_tol))
    {
        result.failure_code = 1;
        return result;
    }

    if (!ls_cell.bernstein_coeffs.empty())
    {
        std::vector<T> edge_coeffs;
        restrict_edge_bernstein_exact<T>(
            ls_cell.cell_type,
            ls_cell.bernstein_order,
            std::span<const T>(ls_cell.bernstein_coeffs),
            std::span<const T>(x0.data(), static_cast<std::size_t>(chart.tdim)),
            std::span<const T>(x1.data(), static_cast<std::size_t>(chart.tdim)),
            edge_coeffs);

        T green_split_t = T(0);
        bool has_green = false;
        const auto tag = classify_edge_roots<T>(
            std::span<const T>(edge_coeffs.data(), edge_coeffs.size()),
            static_cast<T>(opts.zero_tol),
            static_cast<T>(opts.sign_tol),
            opts.edge_root_max_depth,
            green_split_t,
            has_green);
        (void)green_split_t;
        (void)has_green;
        if (tag != EdgeRootTag::one_root)
        {
            result.failure_code = 2;
            return result;
        }
    }
    else if (f0 * f1 > T(0))
    {
        result.failure_code = 2;
        return result;
    }

    auto phi = [&](std::span<const T> xi) -> T
    {
        return ls_cell.value(xi);
    };

    const auto info = cell::edge_root::find_root_parameter_info<T>(
        std::span<const T>(x0.data(), static_cast<std::size_t>(chart.tdim)),
        std::span<const T>(x1.data(), static_cast<std::size_t>(chart.tdim)),
        phi,
        opts.root_method,
        T(0),
        /*max_iter=*/64,
        static_cast<T>(opts.root_tol),
        static_cast<T>(opts.root_tol));
    if (!info.converged || info.residual > static_cast<T>(T(10) * opts.root_tol))
    {
        result.failure_code = 3;
        return result;
    }

    result.r = std::clamp(info.t, T(0), T(1));
    chart.eval(result.r, y,
               std::span<T>(result.xi.data(), static_cast<std::size_t>(chart.tdim)));
    ls_cell.grad(
        std::span<const T>(result.xi.data(), static_cast<std::size_t>(chart.tdim)),
        std::span<T>(result.grad.data(), static_cast<std::size_t>(chart.tdim)));

    std::array<T, 3> cr{};
    std::vector<T> cy;
    chart.derivatives(
        result.r,
        y,
        std::span<T>(cr.data(), static_cast<std::size_t>(chart.tdim)),
        cy);
    result.den = dot<T>(
        std::span<const T>(result.grad.data(), static_cast<std::size_t>(chart.tdim)),
        std::span<const T>(cr.data(), static_cast<std::size_t>(chart.tdim)));

    const T grad_norm = norm<T>(
        std::span<const T>(result.grad.data(), static_cast<std::size_t>(chart.tdim)));
    const T cr_norm = norm<T>(
        std::span<const T>(cr.data(), static_cast<std::size_t>(chart.tdim)));
    // This is a hard numerical validity check, not the soft chart-quality
    // transversality threshold. A low but nonzero denominator means the Duffy
    // chart may need recovery/splitting/refinement; it is not by itself an
    // invalid root branch.
    const T threshold = std::max(
        static_cast<T>(opts.root_tol),
        static_cast<T>(100) * std::numeric_limits<T>::epsilon() * grad_norm * cr_norm);
    if (!(result.den > threshold))
    {
        result.failure_code = 4;
        return result;
    }

    result.valid = true;
    result.failure_code = 0;
    return result;
}

template <std::floating_point T>
void sort_unique_with_tolerance(std::vector<T>& values, T tol)
{
    std::sort(values.begin(), values.end());
    auto it = std::unique(
        values.begin(),
        values.end(),
        [tol](T a, T b) { return std::abs(a - b) <= tol; });
    values.erase(it, values.end());
}

template <std::floating_point T, std::integral I>
bool choose_triangle_height_frame(std::span<const T> ref_vertices,
                                  const LevelSetCell<T, I>& ls_cell,
                                  HeightDirectionFrame<T>& frame)
{
    if (ref_vertices.size() != 6)
        return false;

    std::array<T, 3> centroid{};
    for (int v = 0; v < 3; ++v)
    {
        centroid[0] += ref_vertices[static_cast<std::size_t>(2 * v)];
        centroid[1] += ref_vertices[static_cast<std::size_t>(2 * v + 1)];
    }
    centroid[0] /= T(3);
    centroid[1] /= T(3);

    std::array<T, 3> grad{};
    ls_cell.grad(std::span<const T>(centroid.data(), 2),
                 std::span<T>(grad.data(), 2));
    const T grad_norm = norm<T>(std::span<const T>(grad.data(), 2));
    if (grad_norm <= std::numeric_limits<T>::epsilon())
        return false;

    frame.tdim = 2;
    frame.base_dim = 1;
    frame.origin = {T(0), T(0), T(0)};
    frame.er = {grad[0] / grad_norm, grad[1] / grad_norm, T(0)};
    frame.base = {-frame.er[1], frame.er[0], T(0), T(0), T(0), T(0)};
    return true;
}

template <std::floating_point T>
void triangle_barycentric_affine(std::span<const T> ref_vertices,
                                 std::span<const T> point,
                                 std::span<const T> direction,
                                 std::span<T> alpha,
                                 std::span<T> beta)
{
    const T x0 = ref_vertices[0];
    const T y0 = ref_vertices[1];
    const T ax = ref_vertices[2] - x0;
    const T ay = ref_vertices[3] - y0;
    const T bx = ref_vertices[4] - x0;
    const T by = ref_vertices[5] - y0;
    const T det = ax * by - ay * bx;
    if (std::abs(det) <= std::numeric_limits<T>::epsilon())
        throw std::runtime_error("implicit quadrature: degenerate triangle leaf");

    const auto eval = [&](T px, T py, std::span<T> out)
    {
        const T rx = px - x0;
        const T ry = py - y0;
        out[1] = (rx * by - ry * bx) / det;
        out[2] = (ax * ry - ay * rx) / det;
        out[0] = T(1) - out[1] - out[2];
    };

    eval(point[0], point[1], alpha);

    std::array<T, 3> lambda_at_direction_end{};
    eval(point[0] + direction[0], point[1] + direction[1],
         std::span<T>(lambda_at_direction_end.data(), 3));
    for (int i = 0; i < 3; ++i)
        beta[static_cast<std::size_t>(i)] =
            lambda_at_direction_end[static_cast<std::size_t>(i)] - alpha[static_cast<std::size_t>(i)];
}

template <std::floating_point T>
LineClipResult<T> clip_triangle_height_line(std::span<const T> ref_vertices,
                                            const HeightDirectionFrame<T>& frame,
                                            T s,
                                            T tol)
{
    LineClipResult<T> result;
    result.r_lo = -std::numeric_limits<T>::infinity();
    result.r_hi = std::numeric_limits<T>::infinity();

    std::array<T, 2> point{
        frame.origin[0] + frame.base[0] * s,
        frame.origin[1] + frame.base[1] * s};
    std::array<T, 2> er{frame.er[0], frame.er[1]};
    std::array<T, 3> alpha{};
    std::array<T, 3> beta{};
    triangle_barycentric_affine<T>(
        ref_vertices,
        std::span<const T>(point.data(), 2),
        std::span<const T>(er.data(), 2),
        std::span<T>(alpha.data(), 3),
        std::span<T>(beta.data(), 3));

    for (int i = 0; i < 3; ++i)
    {
        const T a = alpha[static_cast<std::size_t>(i)];
        const T b = beta[static_cast<std::size_t>(i)];
        if (b > tol)
        {
            const T r = -a / b;
            if (r > result.r_lo)
            {
                result.r_lo = r;
                result.lower_facet = i;
            }
        }
        else if (b < -tol)
        {
            const T r = -a / b;
            if (r < result.r_hi)
            {
                result.r_hi = r;
                result.upper_facet = i;
            }
        }
        else if (a < -tol)
        {
            return result;
        }
    }

    result.valid = std::isfinite(result.r_lo) && std::isfinite(result.r_hi)
                && result.r_hi >= result.r_lo - tol;
    return result;
}

template <std::floating_point T>
void eval_height_point(const HeightDirectionFrame<T>& frame, T s, T r, std::span<T> xi)
{
    xi[0] = frame.origin[0] + frame.base[0] * s + frame.er[0] * r;
    xi[1] = frame.origin[1] + frame.base[1] * s + frame.er[1] * r;
}

template <std::floating_point T, std::integral I>
T eval_height_phi(const HeightDirectionFrame<T>& frame,
                  const LevelSetCell<T, I>& ls_cell,
                  T s,
                  T r)
{
    std::array<T, 3> xi{};
    eval_height_point<T>(frame, s, r, std::span<T>(xi.data(), 2));
    return ls_cell.value(std::span<const T>(xi.data(), 2));
}

template <std::floating_point T, std::integral I>
void append_endpoint_roots_on_base_interval(std::span<const T> ref_vertices,
                                            const HeightDirectionFrame<T>& frame,
                                            const LevelSetCell<T, I>& ls_cell,
                                            T s0,
                                            T s1,
                                            bool lower_endpoint,
                                            const ImplicitQuadratureOptions& opts,
                                            std::vector<T>& split_points)
{
    const int sample_count = std::max(12, 2 * opts.inspection_order + 8);
    const T clip_tol = static_cast<T>(100) * static_cast<T>(opts.root_tol);
    auto endpoint_phi = [&](T s, bool& valid) -> T
    {
        const auto clip = clip_triangle_height_line<T>(
            ref_vertices, frame, s, clip_tol);
        valid = clip.valid;
        if (!clip.valid)
            return std::numeric_limits<T>::quiet_NaN();
        return eval_height_phi<T, I>(
            frame, ls_cell, s, lower_endpoint ? clip.r_lo : clip.r_hi);
    };

    T prev_s = s0;
    bool prev_valid = false;
    T prev_f = endpoint_phi(prev_s, prev_valid);
    if (prev_valid && std::abs(prev_f) <= static_cast<T>(opts.root_tol))
        split_points.push_back(prev_s);

    for (int i = 1; i <= sample_count; ++i)
    {
        const T s = s0 + (s1 - s0) * (static_cast<T>(i) / static_cast<T>(sample_count));
        bool valid = false;
        const T f = endpoint_phi(s, valid);
        if (valid && std::abs(f) <= static_cast<T>(opts.root_tol))
            split_points.push_back(s);

        if (prev_valid && valid && prev_f * f < T(0))
        {
            auto eval = [&](T x) -> T
            {
                bool ok = false;
                const T value = endpoint_phi(x, ok);
                return ok ? value : std::numeric_limits<T>::quiet_NaN();
            };
            int iterations = 0;
            bool converged = false;
            const T root = cell::edge_root::brent_solve<T>(
                eval,
                prev_s,
                s,
                prev_f,
                f,
                64,
                static_cast<T>(opts.root_tol),
                static_cast<T>(opts.root_tol),
                &iterations,
                &converged);
            (void)iterations;
            if (converged)
                split_points.push_back(std::clamp(root, s0, s1));
        }

        prev_s = s;
        prev_valid = valid;
        prev_f = f;
    }
}

template <std::floating_point T, std::integral I>
std::vector<T> height_line_roots(std::span<const T> ref_vertices,
                                 const HeightDirectionFrame<T>& frame,
                                 const LevelSetCell<T, I>& ls_cell,
                                 T s,
                                 const ImplicitQuadratureOptions& opts)
{
    std::vector<T> roots;
    const T clip_tol = static_cast<T>(100) * static_cast<T>(opts.root_tol);
    const auto clip = clip_triangle_height_line<T>(
        ref_vertices, frame, s, clip_tol);
    if (!clip.valid || clip.r_hi - clip.r_lo <= clip_tol)
        return roots;

    const int sample_count = std::max(16, 4 * opts.order + 8);
    auto eval = [&](T r) -> T
    {
        return eval_height_phi<T, I>(frame, ls_cell, s, r);
    };

    T prev_r = clip.r_lo;
    T prev_f = eval(prev_r);
    if (std::abs(prev_f) <= static_cast<T>(opts.root_tol))
        roots.push_back(prev_r);

    for (int i = 1; i <= sample_count; ++i)
    {
        const T r = clip.r_lo
                  + (clip.r_hi - clip.r_lo)
                        * (static_cast<T>(i) / static_cast<T>(sample_count));
        const T f = eval(r);
        if (std::abs(f) <= static_cast<T>(opts.root_tol))
            roots.push_back(r);

        if (prev_f * f < T(0))
        {
            int iterations = 0;
            bool converged = false;
            const T root = cell::edge_root::brent_solve<T>(
                eval,
                prev_r,
                r,
                prev_f,
                f,
                64,
                static_cast<T>(opts.root_tol),
                static_cast<T>(opts.root_tol),
                &iterations,
                &converged);
            (void)iterations;
            if (converged)
                roots.push_back(std::clamp(root, clip.r_lo, clip.r_hi));
        }

        prev_r = r;
        prev_f = f;
    }

    sort_unique_with_tolerance<T>(
        roots, static_cast<T>(100) * static_cast<T>(opts.root_tol));
    return roots;
}

template <std::floating_point T, std::integral I>
bool append_triangle_height_recovery_quadrature(
    quadrature::QuadratureRules<T>& rules,
    const AdaptCell<T>& ac,
    int cell_id,
    const LevelSetCell<T, I>& ls_cell,
    IntegrationDomain domain,
    const ImplicitQuadratureOptions& opts,
    int parent_cell_id,
    std::int64_t* plan_component_hash = nullptr)
{
    if (ac.tdim != 2)
        return false;
    const auto leaf_type = ac.entity_types[ac.tdim][static_cast<std::size_t>(cell_id)];
    if (leaf_type != cell::type::triangle)
        return false;
    if (ls_cell.parent_vertex_coords.empty())
        return false;

    std::vector<T> ref_vertices;
    gather_leaf_ref_vertices(ac, cell_id, ref_vertices);

    HeightDirectionFrame<T> frame;
    if (!choose_triangle_height_frame<T, I>(
            std::span<const T>(ref_vertices.data(), ref_vertices.size()),
            ls_cell,
            frame))
    {
        return false;
    }

    std::vector<T> base_points;
    base_points.reserve(6);
    for (int v = 0; v < 3; ++v)
    {
        const T x = ref_vertices[static_cast<std::size_t>(2 * v)];
        const T y = ref_vertices[static_cast<std::size_t>(2 * v + 1)];
        base_points.push_back(
            (x - frame.origin[0]) * frame.base[0]
          + (y - frame.origin[1]) * frame.base[1]);
    }
    sort_unique_with_tolerance<T>(
        base_points, static_cast<T>(100) * static_cast<T>(opts.root_tol));
    if (base_points.size() < 2)
        return false;

    std::vector<T> split_points = base_points;
    for (std::size_t i = 0; i + 1 < base_points.size(); ++i)
    {
        const T s0 = base_points[i];
        const T s1 = base_points[i + 1];
        if (s1 - s0 <= static_cast<T>(opts.root_tol))
            continue;
        append_endpoint_roots_on_base_interval<T, I>(
            std::span<const T>(ref_vertices.data(), ref_vertices.size()),
            frame,
            ls_cell,
            s0,
            s1,
            true,
            opts,
            split_points);
        append_endpoint_roots_on_base_interval<T, I>(
            std::span<const T>(ref_vertices.data(), ref_vertices.size()),
            frame,
            ls_cell,
            s0,
            s1,
            false,
            opts,
            split_points);
    }
    sort_unique_with_tolerance<T>(
        split_points, static_cast<T>(100) * static_cast<T>(opts.root_tol));

    const auto base_rule =
        quadrature::get_reference_rule<T>(cell::type::interval, opts.order);
    const auto radial_rule =
        quadrature::get_reference_rule<T>(cell::type::interval, opts.order);
    const int tdim = 2;
    const int gdim = ls_cell.gdim;
    const T clip_tol = static_cast<T>(100) * static_cast<T>(opts.root_tol);

    std::vector<T> points;
    std::vector<T> weights;
    std::array<T, 3> xi{};
    std::array<T, 3> grad{};
    std::array<T, 4> tangents_ref{};

    for (std::size_t interval_id = 0; interval_id + 1 < split_points.size(); ++interval_id)
    {
        const T s0 = split_points[interval_id];
        const T s1 = split_points[interval_id + 1];
        const T base_len = s1 - s0;
        if (base_len <= static_cast<T>(opts.root_tol))
            continue;

        if (domain == IntegrationDomain::interface_surface)
        {
            auto midpoint_roots = height_line_roots<T, I>(
                std::span<const T>(ref_vertices.data(), ref_vertices.size()),
                frame,
                ls_cell,
                T(0.5) * (s0 + s1),
                opts);
            if (midpoint_roots.empty())
                continue;
        }

        std::array<T, 2> base_tangent{
            frame.base[0] * base_len,
            frame.base[1] * base_len};

        for (int qb = 0; qb < base_rule._num_points; ++qb)
        {
            const T u = base_rule._points[static_cast<std::size_t>(qb)];
            const T s = s0 + base_len * u;
            const auto clip = clip_triangle_height_line<T>(
                std::span<const T>(ref_vertices.data(), ref_vertices.size()),
                frame,
                s,
                clip_tol);
            if (!clip.valid || clip.r_hi - clip.r_lo <= clip_tol)
                continue;

            auto roots = height_line_roots<T, I>(
                std::span<const T>(ref_vertices.data(), ref_vertices.size()),
                frame,
                ls_cell,
                s,
                opts);

            if (domain == IntegrationDomain::interface_surface)
            {
                for (const T root : roots)
                {
                    eval_height_point<T>(
                        frame, s, root, std::span<T>(xi.data(), 2));
                    ls_cell.grad(
                        std::span<const T>(xi.data(), 2),
                        std::span<T>(grad.data(), 2));
                    const T den = grad[0] * frame.er[0] + grad[1] * frame.er[1];
                    const T grad_norm =
                        norm<T>(std::span<const T>(grad.data(), 2));
                    const T threshold = std::max(
                        static_cast<T>(opts.root_tol),
                        static_cast<T>(opts.min_height_transversality) * grad_norm);
                    if (std::abs(den) <= threshold)
                        continue;

                    const T drdu =
                        -(grad[0] * base_tangent[0] + grad[1] * base_tangent[1]) / den;
                    tangents_ref[0] = base_tangent[0] + frame.er[0] * drdu;
                    tangents_ref[1] = base_tangent[1] + frame.er[1] * drdu;

                    const T measure = physical_measure_from_ref_tangents<T>(
                        ls_cell.cell_type,
                        std::span<const T>(
                            ls_cell.parent_vertex_coords.data(),
                            ls_cell.parent_vertex_coords.size()),
                        gdim,
                        std::span<const T>(xi.data(), 2),
                        std::span<const T>(tangents_ref.data(), 2),
                        1);
                    points.insert(points.end(), xi.begin(), xi.begin() + 2);
                    weights.push_back(
                        base_rule._weights[static_cast<std::size_t>(qb)] * measure);
                }
                continue;
            }

            std::vector<T> r_points;
            r_points.reserve(roots.size() + 2);
            r_points.push_back(clip.r_lo);
            r_points.insert(r_points.end(), roots.begin(), roots.end());
            r_points.push_back(clip.r_hi);
            sort_unique_with_tolerance<T>(
                r_points, static_cast<T>(100) * static_cast<T>(opts.root_tol));

            for (std::size_t j = 0; j + 1 < r_points.size(); ++j)
            {
                const T r0 = r_points[j];
                const T r1 = r_points[j + 1];
                const T radial_len = r1 - r0;
                if (radial_len <= static_cast<T>(opts.root_tol))
                    continue;

                const T r_mid = T(0.5) * (r0 + r1);
                const T f_mid = eval_height_phi<T, I>(frame, ls_cell, s, r_mid);
                const bool accept =
                    (domain == IntegrationDomain::negative_volume
                     && f_mid < static_cast<T>(opts.sign_tol))
                 || (domain == IntegrationDomain::positive_volume
                     && f_mid > -static_cast<T>(opts.sign_tol));
                if (!accept)
                    continue;

                for (int qr = 0; qr < radial_rule._num_points; ++qr)
                {
                    const T v = radial_rule._points[static_cast<std::size_t>(qr)];
                    const T r = r0 + radial_len * v;
                    eval_height_point<T>(
                        frame, s, r, std::span<T>(xi.data(), 2));

                    tangents_ref[0] = base_tangent[0];
                    tangents_ref[1] = base_tangent[1];
                    tangents_ref[2] = frame.er[0] * radial_len;
                    tangents_ref[3] = frame.er[1] * radial_len;

                    const T measure = physical_measure_from_ref_tangents<T>(
                        ls_cell.cell_type,
                        std::span<const T>(
                            ls_cell.parent_vertex_coords.data(),
                            ls_cell.parent_vertex_coords.size()),
                        gdim,
                        std::span<const T>(xi.data(), 2),
                        std::span<const T>(tangents_ref.data(), 4),
                        2);
                    points.insert(points.end(), xi.begin(), xi.begin() + 2);
                    weights.push_back(
                        base_rule._weights[static_cast<std::size_t>(qb)]
                      * radial_rule._weights[static_cast<std::size_t>(qr)]
                      * measure);
                }
            }
        }
    }

    if (weights.empty())
        return false;

    append_rule_points(
        rules,
        std::span<const T>(points.data(), points.size()),
        std::span<const T>(weights.data(), weights.size()),
        tdim,
        parent_cell_id);
    if (plan_component_hash != nullptr)
    {
        *plan_component_hash = triangle_height_plan_component<T>(
            frame, std::span<const T>(split_points.data(), split_points.size()));
    }
    return true;
}

template <std::floating_point T, std::integral I>
bool choose_tetra_height_frame(std::span<const T> ref_vertices,
                               const LevelSetCell<T, I>& ls_cell,
                               HeightDirectionFrame<T>& frame)
{
    if (ref_vertices.size() != 12)
        return false;

    std::array<T, 3> centroid{};
    for (int v = 0; v < 4; ++v)
    {
        for (int d = 0; d < 3; ++d)
            centroid[static_cast<std::size_t>(d)]
                += ref_vertices[static_cast<std::size_t>(3 * v + d)];
    }
    for (auto& value : centroid)
        value /= T(4);

    std::array<T, 3> grad{};
    ls_cell.grad(std::span<const T>(centroid.data(), 3),
                 std::span<T>(grad.data(), 3));
    if (!normalize3<T>(grad))
        return false;

    std::array<T, 3> helper =
        (std::abs(grad[2]) < T(0.85))
            ? std::array<T, 3>{T(0), T(0), T(1)}
            : std::array<T, 3>{T(0), T(1), T(0)};
    std::array<T, 3> es = cross3<T>(helper, grad);
    if (!normalize3<T>(es))
        return false;
    std::array<T, 3> et = cross3<T>(grad, es);
    if (!normalize3<T>(et))
        return false;

    frame.tdim = 3;
    frame.base_dim = 2;
    frame.origin = {T(0), T(0), T(0)};
    frame.er = grad;
    frame.base = {es[0], es[1], es[2], et[0], et[1], et[2]};
    return true;
}

template <std::floating_point T>
bool make_tetra_height_frame_from_direction(std::array<T, 3> er,
                                            HeightDirectionFrame<T>& frame)
{
    if (!normalize3<T>(er))
        return false;

    std::array<T, 3> helper{T(1), T(0), T(0)};
    if (std::abs(er[1]) < std::abs(er[0]) && std::abs(er[1]) <= std::abs(er[2]))
        helper = {T(0), T(1), T(0)};
    else if (std::abs(er[2]) < std::abs(er[0]) && std::abs(er[2]) < std::abs(er[1]))
        helper = {T(0), T(0), T(1)};

    std::array<T, 3> es = cross3<T>(helper, er);
    if (!normalize3<T>(es))
        return false;
    std::array<T, 3> et = cross3<T>(er, es);
    if (!normalize3<T>(et))
        return false;

    frame.tdim = 3;
    frame.base_dim = 2;
    frame.origin = {T(0), T(0), T(0)};
    frame.er = er;
    frame.base = {es[0], es[1], es[2], et[0], et[1], et[2]};
    return true;
}

template <std::floating_point T>
void append_unique_direction(std::vector<std::array<T, 3>>& directions,
                             std::array<T, 3> direction,
                             T cosine_tol = T(0.98))
{
    if (!normalize3<T>(direction))
        return;
    for (const auto& existing : directions)
    {
        if (std::abs(dot3<T>(existing, direction)) >= cosine_tol)
            return;
    }
    directions.push_back(direction);
}

template <std::floating_point T, std::integral I>
std::vector<HeightDirectionFrame<T>> enumerate_tetra_height_frames(
    std::span<const T> ref_vertices,
    const LevelSetCell<T, I>& ls_cell)
{
    std::vector<std::array<T, 3>> directions;
    std::array<T, 3> centroid{};
    for (int v = 0; v < 4; ++v)
    {
        for (int d = 0; d < 3; ++d)
            centroid[static_cast<std::size_t>(d)]
                += ref_vertices[static_cast<std::size_t>(3 * v + d)];
    }
    for (auto& value : centroid)
        value /= T(4);

    std::array<T, 3> grad{};
    ls_cell.grad(std::span<const T>(centroid.data(), 3),
                 std::span<T>(grad.data(), 3));
    append_unique_direction<T>(directions, grad);

    // Grazing cap cells often have one vertex inserted by multi-root edge
    // splitting. The gradient at that isolated vertex or at the nearby edge
    // roots is a better graph direction than the sign-oriented Duffy ray.
    std::array<std::array<T, 3>, 4> vertex_points{};
    std::array<T, 4> vertex_values{};
    for (int v = 0; v < 4; ++v)
    {
        for (int d = 0; d < 3; ++d)
        {
            vertex_points[static_cast<std::size_t>(v)][static_cast<std::size_t>(d)] =
                ref_vertices[static_cast<std::size_t>(3 * v + d)];
        }
        vertex_values[static_cast<std::size_t>(v)] = ls_cell.value(
            std::span<const T>(
                vertex_points[static_cast<std::size_t>(v)].data(), 3));
        std::array<T, 3> vertex_grad{};
        ls_cell.grad(
            std::span<const T>(
                vertex_points[static_cast<std::size_t>(v)].data(), 3),
            std::span<T>(vertex_grad.data(), 3));
        append_unique_direction<T>(directions, vertex_grad);
    }

    constexpr int edge_vertices[6][2] = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    for (const auto& edge : edge_vertices)
    {
        const int a = edge[0];
        const int b = edge[1];
        const T fa = vertex_values[static_cast<std::size_t>(a)];
        const T fb = vertex_values[static_cast<std::size_t>(b)];
        if (!(fa * fb < T(0)))
            continue;

        auto phi = [&](std::span<const T> xi) -> T
        {
            return ls_cell.value(xi);
        };
        const auto info = cell::edge_root::find_root_parameter_info<T>(
            std::span<const T>(
                vertex_points[static_cast<std::size_t>(a)].data(), 3),
            std::span<const T>(
                vertex_points[static_cast<std::size_t>(b)].data(), 3),
            phi,
            cell::edge_root::method::itp,
            T(0),
            64,
            T(1e-12),
            T(1e-12));
        if (!info.converged)
            continue;

        std::array<T, 3> root{};
        for (int d = 0; d < 3; ++d)
        {
            root[static_cast<std::size_t>(d)] =
                vertex_points[static_cast<std::size_t>(a)][static_cast<std::size_t>(d)]
              + info.t
                    * (vertex_points[static_cast<std::size_t>(b)][static_cast<std::size_t>(d)]
                       - vertex_points[static_cast<std::size_t>(a)][static_cast<std::size_t>(d)]);
        }
        std::array<T, 3> root_grad{};
        ls_cell.grad(
            std::span<const T>(root.data(), 3),
            std::span<T>(root_grad.data(), 3));
        append_unique_direction<T>(directions, root_grad);
    }

    std::array<T, 3> neg_centroid{};
    std::array<T, 3> pos_centroid{};
    int neg_count = 0;
    int pos_count = 0;
    for (int v = 0; v < 4; ++v)
    {
        std::array<T, 3> xi{
            ref_vertices[static_cast<std::size_t>(3 * v)],
            ref_vertices[static_cast<std::size_t>(3 * v + 1)],
            ref_vertices[static_cast<std::size_t>(3 * v + 2)]};
        const T value = ls_cell.value(std::span<const T>(xi.data(), 3));
        auto& target = value < T(0) ? neg_centroid : pos_centroid;
        int& count = value < T(0) ? neg_count : pos_count;
        for (int d = 0; d < 3; ++d)
            target[static_cast<std::size_t>(d)] += xi[static_cast<std::size_t>(d)];
        ++count;
    }
    if (neg_count > 0 && pos_count > 0)
    {
        for (int d = 0; d < 3; ++d)
        {
            neg_centroid[static_cast<std::size_t>(d)] /= static_cast<T>(neg_count);
            pos_centroid[static_cast<std::size_t>(d)] /= static_cast<T>(pos_count);
            pos_centroid[static_cast<std::size_t>(d)] -= neg_centroid[static_cast<std::size_t>(d)];
        }
        append_unique_direction<T>(directions, pos_centroid);
    }

    append_unique_direction<T>(directions, {T(1), T(0), T(0)}, T(0.999));
    append_unique_direction<T>(directions, {T(0), T(1), T(0)}, T(0.999));
    append_unique_direction<T>(directions, {T(0), T(0), T(1)}, T(0.999));
    append_unique_direction<T>(directions, {T(1), T(1), T(0)}, T(0.999));
    append_unique_direction<T>(directions, {T(1), T(0), T(1)}, T(0.999));
    append_unique_direction<T>(directions, {T(0), T(1), T(1)}, T(0.999));
    append_unique_direction<T>(directions, {T(1), T(1), T(1)}, T(0.999));

    constexpr int faces[4][3] = {
        {0, 1, 2},
        {0, 1, 3},
        {0, 2, 3},
        {1, 2, 3}};
    for (const auto& face : faces)
    {
        std::array<T, 3> p0{};
        std::array<T, 3> e0{};
        std::array<T, 3> e1{};
        for (int d = 0; d < 3; ++d)
        {
            p0[static_cast<std::size_t>(d)] =
                ref_vertices[static_cast<std::size_t>(3 * face[0] + d)];
            e0[static_cast<std::size_t>(d)] =
                ref_vertices[static_cast<std::size_t>(3 * face[1] + d)]
              - p0[static_cast<std::size_t>(d)];
            e1[static_cast<std::size_t>(d)] =
                ref_vertices[static_cast<std::size_t>(3 * face[2] + d)]
              - p0[static_cast<std::size_t>(d)];
        }
        append_unique_direction<T>(directions, cross3<T>(e0, e1), T(0.999));
    }

    std::vector<HeightDirectionFrame<T>> frames;
    frames.reserve(directions.size());
    for (const auto& direction : directions)
    {
        HeightDirectionFrame<T> frame;
        if (make_tetra_height_frame_from_direction<T>(direction, frame))
            frames.push_back(frame);
    }
    return frames;
}

template <std::floating_point T>
T det3_cols(const std::array<T, 3>& a,
            const std::array<T, 3>& b,
            const std::array<T, 3>& c)
{
    return dot3<T>(a, cross3<T>(b, c));
}

template <std::floating_point T>
std::array<T, 3> solve_3x3_cols(const std::array<T, 3>& a,
                                const std::array<T, 3>& b,
                                const std::array<T, 3>& c,
                                const std::array<T, 3>& rhs)
{
    const T det = det3_cols<T>(a, b, c);
    if (std::abs(det) <= std::numeric_limits<T>::epsilon())
        throw std::runtime_error("implicit quadrature: degenerate tetrahedron leaf");
    return {
        det3_cols<T>(rhs, b, c) / det,
        det3_cols<T>(a, rhs, c) / det,
        det3_cols<T>(a, b, rhs) / det};
}

template <std::floating_point T>
void tetra_barycentric_affine(std::span<const T> ref_vertices,
                              std::span<const T> point,
                              std::span<const T> direction,
                              std::span<T> alpha,
                              std::span<T> beta)
{
    std::array<T, 3> v0{
        ref_vertices[0], ref_vertices[1], ref_vertices[2]};
    std::array<T, 3> a{};
    std::array<T, 3> b{};
    std::array<T, 3> c{};
    std::array<T, 3> rhs{};
    std::array<T, 3> dir{};
    for (int d = 0; d < 3; ++d)
    {
        a[static_cast<std::size_t>(d)] =
            ref_vertices[static_cast<std::size_t>(3 + d)] - v0[static_cast<std::size_t>(d)];
        b[static_cast<std::size_t>(d)] =
            ref_vertices[static_cast<std::size_t>(6 + d)] - v0[static_cast<std::size_t>(d)];
        c[static_cast<std::size_t>(d)] =
            ref_vertices[static_cast<std::size_t>(9 + d)] - v0[static_cast<std::size_t>(d)];
        rhs[static_cast<std::size_t>(d)] = point[static_cast<std::size_t>(d)] - v0[static_cast<std::size_t>(d)];
        dir[static_cast<std::size_t>(d)] = direction[static_cast<std::size_t>(d)];
    }

    const auto lambda123 = solve_3x3_cols<T>(a, b, c, rhs);
    alpha[1] = lambda123[0];
    alpha[2] = lambda123[1];
    alpha[3] = lambda123[2];
    alpha[0] = T(1) - alpha[1] - alpha[2] - alpha[3];

    const auto dlambda123 = solve_3x3_cols<T>(a, b, c, dir);
    beta[1] = dlambda123[0];
    beta[2] = dlambda123[1];
    beta[3] = dlambda123[2];
    beta[0] = -beta[1] - beta[2] - beta[3];
}

template <std::floating_point T>
void eval_height_point(const HeightDirectionFrame<T>& frame,
                       std::span<const T> y,
                       T r,
                       std::span<T> xi)
{
    for (int d = 0; d < frame.tdim; ++d)
    {
        xi[static_cast<std::size_t>(d)] =
            frame.origin[static_cast<std::size_t>(d)]
          + frame.er[static_cast<std::size_t>(d)] * r;
    }
    for (int a = 0; a < frame.base_dim; ++a)
    {
        for (int d = 0; d < frame.tdim; ++d)
        {
            xi[static_cast<std::size_t>(d)] +=
                frame.base[static_cast<std::size_t>(a * frame.tdim + d)]
              * y[static_cast<std::size_t>(a)];
        }
    }
}

template <std::floating_point T>
LineClipResult<T> clip_tetra_height_line(std::span<const T> ref_vertices,
                                         const HeightDirectionFrame<T>& frame,
                                         std::span<const T> y,
                                         T tol)
{
    LineClipResult<T> result;
    result.r_lo = -std::numeric_limits<T>::infinity();
    result.r_hi = std::numeric_limits<T>::infinity();

    std::array<T, 3> point{};
    eval_height_point<T>(frame, y, T(0), std::span<T>(point.data(), 3));
    std::array<T, 4> alpha{};
    std::array<T, 4> beta{};
    tetra_barycentric_affine<T>(
        ref_vertices,
        std::span<const T>(point.data(), 3),
        std::span<const T>(frame.er.data(), 3),
        std::span<T>(alpha.data(), 4),
        std::span<T>(beta.data(), 4));

    for (int i = 0; i < 4; ++i)
    {
        const T a = alpha[static_cast<std::size_t>(i)];
        const T b = beta[static_cast<std::size_t>(i)];
        if (b > tol)
        {
            const T r = -a / b;
            if (r > result.r_lo)
            {
                result.r_lo = r;
                result.lower_facet = i;
            }
        }
        else if (b < -tol)
        {
            const T r = -a / b;
            if (r < result.r_hi)
            {
                result.r_hi = r;
                result.upper_facet = i;
            }
        }
        else if (a < -tol)
        {
            return result;
        }
    }

    result.valid = std::isfinite(result.r_lo) && std::isfinite(result.r_hi)
                && result.r_hi >= result.r_lo - tol;
    return result;
}

template <std::floating_point T, std::integral I>
T eval_height_phi(const HeightDirectionFrame<T>& frame,
                  const LevelSetCell<T, I>& ls_cell,
                  std::span<const T> y,
                  T r)
{
    std::array<T, 3> xi{};
    eval_height_point<T>(frame, y, r,
                         std::span<T>(xi.data(), static_cast<std::size_t>(frame.tdim)));
    return ls_cell.value(
        std::span<const T>(xi.data(), static_cast<std::size_t>(frame.tdim)));
}

template <std::floating_point T, std::integral I>
std::vector<T> height_line_roots_tetra(std::span<const T> ref_vertices,
                                       const HeightDirectionFrame<T>& frame,
                                       const LevelSetCell<T, I>& ls_cell,
                                       std::span<const T> y,
                                       const ImplicitQuadratureOptions& opts)
{
    std::vector<T> roots;
    const T clip_tol = static_cast<T>(100) * static_cast<T>(opts.root_tol);
    const auto clip = clip_tetra_height_line<T>(ref_vertices, frame, y, clip_tol);
    if (!clip.valid || clip.r_hi - clip.r_lo <= clip_tol)
        return roots;

    const int sample_count = std::max(20, 4 * opts.order + 12);
    auto eval = [&](T r) -> T
    {
        return eval_height_phi<T, I>(frame, ls_cell, y, r);
    };

    T prev_r = clip.r_lo;
    T prev_f = eval(prev_r);
    if (std::abs(prev_f) <= static_cast<T>(opts.root_tol))
        roots.push_back(prev_r);

    for (int i = 1; i <= sample_count; ++i)
    {
        const T r = clip.r_lo
                  + (clip.r_hi - clip.r_lo)
                        * (static_cast<T>(i) / static_cast<T>(sample_count));
        const T f = eval(r);
        if (std::abs(f) <= static_cast<T>(opts.root_tol))
            roots.push_back(r);

        if (prev_f * f < T(0))
        {
            int iterations = 0;
            bool converged = false;
            const T root = cell::edge_root::brent_solve<T>(
                eval,
                prev_r,
                r,
                prev_f,
                f,
                64,
                static_cast<T>(opts.root_tol),
                static_cast<T>(opts.root_tol),
                &iterations,
                &converged);
            (void)iterations;
            if (converged)
                roots.push_back(std::clamp(root, clip.r_lo, clip.r_hi));
        }

        prev_r = r;
        prev_f = f;
    }

    sort_unique_with_tolerance<T>(
        roots, static_cast<T>(100) * static_cast<T>(opts.root_tol));
    return roots;
}

template <std::floating_point T>
T orient2(const std::array<T, 2>& a,
          const std::array<T, 2>& b,
          const std::array<T, 2>& c)
{
    return (b[0] - a[0]) * (c[1] - a[1])
         - (b[1] - a[1]) * (c[0] - a[0]);
}

template <std::floating_point T>
std::vector<BaseTriangle2<T>> build_projected_tetra_base_cells(
    std::span<const T> ref_vertices,
    const HeightDirectionFrame<T>& frame,
    T tol)
{
    std::vector<std::array<T, 2>> projected;
    projected.reserve(4);
    for (int v = 0; v < 4; ++v)
    {
        std::array<T, 3> x{};
        for (int d = 0; d < 3; ++d)
        {
            x[static_cast<std::size_t>(d)] =
                ref_vertices[static_cast<std::size_t>(3 * v + d)]
              - frame.origin[static_cast<std::size_t>(d)];
        }
        projected.push_back({
            x[0] * frame.base[0] + x[1] * frame.base[1] + x[2] * frame.base[2],
            x[0] * frame.base[3] + x[1] * frame.base[4] + x[2] * frame.base[5]});
    }

    std::sort(projected.begin(), projected.end(),
              [](const auto& a, const auto& b)
              {
                  if (a[0] == b[0])
                      return a[1] < b[1];
                  return a[0] < b[0];
              });
    projected.erase(
        std::unique(projected.begin(), projected.end(),
                    [tol](const auto& a, const auto& b)
                    {
                        return std::abs(a[0] - b[0]) <= tol
                            && std::abs(a[1] - b[1]) <= tol;
                    }),
        projected.end());
    if (projected.size() < 3)
        return {};

    std::vector<std::array<T, 2>> lower;
    for (const auto& p : projected)
    {
        while (lower.size() >= 2
               && orient2<T>(lower[lower.size() - 2], lower.back(), p) <= tol)
            lower.pop_back();
        lower.push_back(p);
    }
    std::vector<std::array<T, 2>> upper;
    for (auto it = projected.rbegin(); it != projected.rend(); ++it)
    {
        while (upper.size() >= 2
               && orient2<T>(upper[upper.size() - 2], upper.back(), *it) <= tol)
            upper.pop_back();
        upper.push_back(*it);
    }
    lower.pop_back();
    upper.pop_back();
    lower.insert(lower.end(), upper.begin(), upper.end());
    if (lower.size() < 3)
        return {};

    std::vector<BaseTriangle2<T>> cells;
    cells.reserve(lower.size() - 2);
    const auto p0 = lower[0];
    for (std::size_t i = 1; i + 1 < lower.size(); ++i)
    {
        const auto p1 = lower[i];
        const auto p2 = lower[i + 1];
        if (std::abs(orient2<T>(p0, p1, p2)) <= tol)
            continue;
        BaseTriangle2<T> cell;
        cell.vertices = {p0[0], p0[1], p1[0], p1[1], p2[0], p2[1]};
        cell.depth = 0;
        cells.push_back(cell);
    }
    return cells;
}

template <std::floating_point T>
std::array<T, 2> base_triangle_point(const BaseTriangle2<T>& cell, T u, T v)
{
    const T* p0 = cell.vertices.data();
    const T* p1 = cell.vertices.data() + 2;
    const T* p2 = cell.vertices.data() + 4;
    return {
        p0[0] + u * (p1[0] - p0[0]) + v * (p2[0] - p0[0]),
        p0[1] + u * (p1[1] - p0[1]) + v * (p2[1] - p0[1])};
}

template <std::floating_point T>
std::array<T, 3> base_tangent_ref(const HeightDirectionFrame<T>& frame,
                                  const BaseTriangle2<T>& cell,
                                  int axis)
{
    const int p = axis == 0 ? 2 : 4;
    const T ds = cell.vertices[static_cast<std::size_t>(p)]
               - cell.vertices[0];
    const T dt = cell.vertices[static_cast<std::size_t>(p + 1)]
               - cell.vertices[1];
    return {
        frame.base[0] * ds + frame.base[3] * dt,
        frame.base[1] * ds + frame.base[4] * dt,
        frame.base[2] * ds + frame.base[5] * dt};
}

template <std::floating_point T>
std::array<BaseTriangle2<T>, 4> subdivide_base_triangle(const BaseTriangle2<T>& cell)
{
    std::array<T, 2> p0{cell.vertices[0], cell.vertices[1]};
    std::array<T, 2> p1{cell.vertices[2], cell.vertices[3]};
    std::array<T, 2> p2{cell.vertices[4], cell.vertices[5]};
    std::array<T, 2> m01{T(0.5) * (p0[0] + p1[0]), T(0.5) * (p0[1] + p1[1])};
    std::array<T, 2> m12{T(0.5) * (p1[0] + p2[0]), T(0.5) * (p1[1] + p2[1])};
    std::array<T, 2> m20{T(0.5) * (p2[0] + p0[0]), T(0.5) * (p2[1] + p0[1])};
    const int depth = cell.depth + 1;
    return {
        BaseTriangle2<T>{{p0[0], p0[1], m01[0], m01[1], m20[0], m20[1]}, depth},
        BaseTriangle2<T>{{m01[0], m01[1], p1[0], p1[1], m12[0], m12[1]}, depth},
        BaseTriangle2<T>{{m20[0], m20[1], m12[0], m12[1], p2[0], p2[1]}, depth},
        BaseTriangle2<T>{{m01[0], m01[1], m12[0], m12[1], m20[0], m20[1]}, depth}};
}

template <std::floating_point T>
void append_base_triangle_from_points(std::vector<BaseTriangle2<T>>& out,
                                      const std::vector<std::array<T, 2>>& poly,
                                      int depth,
                                      T tol)
{
    if (poly.size() < 3)
        return;
    const auto p0 = poly[0];
    for (std::size_t i = 1; i + 1 < poly.size(); ++i)
    {
        auto p1 = poly[i];
        auto p2 = poly[i + 1];
        if (std::abs(orient2<T>(p0, p1, p2)) <= tol)
            continue;
        if (orient2<T>(p0, p1, p2) < T(0))
            std::swap(p1, p2);
        out.push_back(BaseTriangle2<T>{
            {p0[0], p0[1], p1[0], p1[1], p2[0], p2[1]},
            depth});
    }
}

template <std::floating_point T>
void deduplicate_polygon(std::vector<std::array<T, 2>>& poly, T tol)
{
    std::vector<std::array<T, 2>> unique;
    unique.reserve(poly.size());
    for (const auto& p : poly)
    {
        if (unique.empty()
            || std::abs(unique.back()[0] - p[0]) > tol
            || std::abs(unique.back()[1] - p[1]) > tol)
        {
            unique.push_back(p);
        }
    }
    if (unique.size() > 1
        && std::abs(unique.front()[0] - unique.back()[0]) <= tol
        && std::abs(unique.front()[1] - unique.back()[1]) <= tol)
    {
        unique.pop_back();
    }
    poly = std::move(unique);
}

template <std::floating_point T>
std::vector<std::array<T, 2>> clip_polygon_by_endpoint_sign(
    const std::array<std::array<T, 2>, 3>& points,
    const std::array<T, 3>& values,
    bool keep_positive,
    T tol)
{
    std::vector<std::array<T, 2>> poly;
    poly.reserve(4);
    for (int i = 0; i < 3; ++i)
    {
        const int j = (i + 1) % 3;
        const T vi = values[static_cast<std::size_t>(i)];
        const T vj = values[static_cast<std::size_t>(j)];
        const bool inside_i = keep_positive ? vi >= -tol : vi <= tol;
        const bool inside_j = keep_positive ? vj >= -tol : vj <= tol;
        if (inside_i)
            poly.push_back(points[static_cast<std::size_t>(i)]);
        if (inside_i != inside_j)
        {
            const T denom = vi - vj;
            if (std::abs(denom) > tol)
            {
                const T a = std::clamp(vi / denom, T(0), T(1));
                const auto& pi = points[static_cast<std::size_t>(i)];
                const auto& pj = points[static_cast<std::size_t>(j)];
                poly.push_back({
                    pi[0] + a * (pj[0] - pi[0]),
                    pi[1] + a * (pj[1] - pi[1])});
            }
        }
    }
    deduplicate_polygon<T>(poly, tol);
    return poly;
}

template <std::floating_point T, std::integral I>
bool tetra_endpoint_phi_at_base_point(std::span<const T> ref_vertices,
                                      const HeightDirectionFrame<T>& frame,
                                      const LevelSetCell<T, I>& ls_cell,
                                      const std::array<T, 2>& y,
                                      bool lower_endpoint,
                                      const ImplicitQuadratureOptions& opts,
                                      T& value)
{
    const auto clip = clip_tetra_height_line<T>(
        ref_vertices,
        frame,
        std::span<const T>(y.data(), 2),
        static_cast<T>(100) * static_cast<T>(opts.root_tol));
    if (!clip.valid)
        return false;
    value = eval_height_phi<T, I>(
        frame,
        ls_cell,
        std::span<const T>(y.data(), 2),
        lower_endpoint ? clip.r_lo : clip.r_hi);
    return true;
}

template <std::floating_point T, std::integral I>
std::vector<BaseTriangle2<T>> split_base_triangle_by_endpoint_function(
    const BaseTriangle2<T>& cell,
    std::span<const T> ref_vertices,
    const HeightDirectionFrame<T>& frame,
    const LevelSetCell<T, I>& ls_cell,
    bool lower_endpoint,
    const ImplicitQuadratureOptions& opts)
{
    const T tol = static_cast<T>(100) * static_cast<T>(opts.root_tol);
    const std::array<std::array<T, 2>, 3> pts = {
        std::array<T, 2>{cell.vertices[0], cell.vertices[1]},
        std::array<T, 2>{cell.vertices[2], cell.vertices[3]},
        std::array<T, 2>{cell.vertices[4], cell.vertices[5]}};
    std::array<T, 3> values{};
    for (int i = 0; i < 3; ++i)
    {
        if (!tetra_endpoint_phi_at_base_point<T, I>(
                ref_vertices, frame, ls_cell, pts[static_cast<std::size_t>(i)],
                lower_endpoint, opts, values[static_cast<std::size_t>(i)]))
        {
            return {cell};
        }
    }

    bool has_negative = false;
    bool has_positive = false;
    for (const T value : values)
    {
        has_negative = has_negative || value < -static_cast<T>(opts.sign_tol);
        has_positive = has_positive || value > static_cast<T>(opts.sign_tol);
    }
    if (!(has_negative && has_positive))
        return {cell};

    std::vector<BaseTriangle2<T>> out;
    auto negative_poly = clip_polygon_by_endpoint_sign<T>(pts, values, false, tol);
    auto positive_poly = clip_polygon_by_endpoint_sign<T>(pts, values, true, tol);
    append_base_triangle_from_points<T>(out, negative_poly, cell.depth, tol);
    append_base_triangle_from_points<T>(out, positive_poly, cell.depth, tol);
    if (out.empty())
        return {cell};
    return out;
}

template <std::floating_point T, std::integral I>
std::vector<BaseTriangle2<T>> split_base_cells_by_endpoint_functions(
    std::vector<BaseTriangle2<T>> cells,
    std::span<const T> ref_vertices,
    const HeightDirectionFrame<T>& frame,
    const LevelSetCell<T, I>& ls_cell,
    const ImplicitQuadratureOptions& opts)
{
    for (const bool lower_endpoint : {true, false})
    {
        std::vector<BaseTriangle2<T>> next;
        for (const auto& cell : cells)
        {
            auto split = split_base_triangle_by_endpoint_function<T, I>(
                cell, ref_vertices, frame, ls_cell, lower_endpoint, opts);
            next.insert(next.end(), split.begin(), split.end());
        }
        cells = std::move(next);
    }
    return cells;
}

template <std::floating_point T, std::integral I>
bool tetra_base_cell_needs_split(const BaseTriangle2<T>& cell,
                                 std::span<const T> ref_vertices,
                                 const HeightDirectionFrame<T>& frame,
                                 const LevelSetCell<T, I>& ls_cell,
                                 const ImplicitQuadratureOptions& opts)
{
    const T clip_tol = static_cast<T>(100) * static_cast<T>(opts.root_tol);
    // Use strictly interior samples. Endpoint splitting deliberately puts root
    // branch entry/exit events on base-cell boundaries; sampling those
    // boundaries would report false root-count changes forever.
    const std::array<std::array<T, 2>, 4> samples = {
        base_triangle_point<T>(cell, T(1) / T(3), T(1) / T(3)),
        base_triangle_point<T>(cell, T(3) / T(5), T(1) / T(5)),
        base_triangle_point<T>(cell, T(1) / T(5), T(3) / T(5)),
        base_triangle_point<T>(cell, T(1) / T(5), T(1) / T(5))};

    bool have_reference = false;
    int lower = -1;
    int upper = -1;
    int root_count = -1;
    for (const auto& sample : samples)
    {
        const auto y = std::span<const T>(sample.data(), 2);
        const auto clip = clip_tetra_height_line<T>(
            ref_vertices, frame, y, clip_tol);
        if (!clip.valid || clip.r_hi - clip.r_lo <= clip_tol)
            continue;
        const int roots = static_cast<int>(
            height_line_roots_tetra<T, I>(ref_vertices, frame, ls_cell, y, opts).size());
        if (!have_reference)
        {
            have_reference = true;
            lower = clip.lower_facet;
            upper = clip.upper_facet;
            root_count = roots;
            continue;
        }
        if (clip.lower_facet != lower || clip.upper_facet != upper || roots != root_count)
            return true;
    }
    return false;
}

template <std::floating_point T, std::integral I>
ProjectedBaseSplitResult<T> split_projected_tetra_base_cells(
    std::vector<BaseTriangle2<T>> cells,
    std::span<const T> ref_vertices,
    const HeightDirectionFrame<T>& frame,
    const LevelSetCell<T, I>& ls_cell,
    const ImplicitQuadratureOptions& opts)
{
    ProjectedBaseSplitResult<T> result;
    std::vector<BaseTriangle2<T>> stack = std::move(cells);
    while (!stack.empty())
    {
        BaseTriangle2<T> cell = stack.back();
        stack.pop_back();
        const T area2 = orient2<T>(
            {cell.vertices[0], cell.vertices[1]},
            {cell.vertices[2], cell.vertices[3]},
            {cell.vertices[4], cell.vertices[5]});
        if (std::abs(area2) <= static_cast<T>(100) * static_cast<T>(opts.root_tol))
            continue;
        if (cell.depth < opts.max_base_split_depth
            && tetra_base_cell_needs_split<T, I>(
                cell, ref_vertices, frame, ls_cell, opts))
        {
            const auto children = subdivide_base_triangle<T>(cell);
            for (const auto& child : children)
                stack.push_back(child);
            continue;
        }
        if (cell.depth >= opts.max_base_split_depth
            && tetra_base_cell_needs_split<T, I>(
                cell, ref_vertices, frame, ls_cell, opts))
        {
            result.success = false;
            result.cells.clear();
            return result;
        }
        if (area2 < T(0))
            std::swap(cell.vertices[2], cell.vertices[4]),
            std::swap(cell.vertices[3], cell.vertices[5]);
        result.cells.push_back(cell);
    }
    return result;
}

template <std::floating_point T>
T tetra_leaf_reference_volume(std::span<const T> ref_vertices)
{
    std::array<T, 3> a{};
    std::array<T, 3> b{};
    std::array<T, 3> c{};
    for (int d = 0; d < 3; ++d)
    {
        a[static_cast<std::size_t>(d)] =
            ref_vertices[static_cast<std::size_t>(3 + d)] - ref_vertices[static_cast<std::size_t>(d)];
        b[static_cast<std::size_t>(d)] =
            ref_vertices[static_cast<std::size_t>(6 + d)] - ref_vertices[static_cast<std::size_t>(d)];
        c[static_cast<std::size_t>(d)] =
            ref_vertices[static_cast<std::size_t>(9 + d)] - ref_vertices[static_cast<std::size_t>(d)];
    }
    return std::abs(det3_cols<T>(a, b, c)) / T(6);
}

template <std::floating_point T>
bool projected_base_cells_cover_tetra_leaf(
    std::span<const T> ref_vertices,
    const HeightDirectionFrame<T>& frame,
    const std::vector<BaseTriangle2<T>>& base_cells,
    const ImplicitQuadratureOptions& opts)
{
    const T clip_tol = static_cast<T>(100) * static_cast<T>(opts.root_tol);
    T volume = T(0);
    for (const auto& cell : base_cells)
    {
        const T area2 = orient2<T>(
            {cell.vertices[0], cell.vertices[1]},
            {cell.vertices[2], cell.vertices[3]},
            {cell.vertices[4], cell.vertices[5]});
        if (std::abs(area2) <= clip_tol)
            continue;
        const auto centroid = base_triangle_point<T>(
            cell, T(1) / T(3), T(1) / T(3));
        const auto clip = clip_tetra_height_line<T>(
            ref_vertices,
            frame,
            std::span<const T>(centroid.data(), 2),
            clip_tol);
        if (!clip.valid || clip.r_hi < clip.r_lo)
            return false;
        volume += T(0.5) * std::abs(area2) * (clip.r_hi - clip.r_lo);
    }

    const T exact = tetra_leaf_reference_volume<T>(ref_vertices);
    const T tolerance = std::max(
        static_cast<T>(1e-8),
        static_cast<T>(opts.max_geometry_q_error))
        * std::max(T(1), exact);
    return std::abs(volume - exact) <= tolerance;
}

template <std::floating_point T>
bool tetra_height_reference_rules_available(int order)
{
    if (order <= 0)
        return false;
    try
    {
        (void)quadrature::get_reference_rule<T>(cell::type::triangle, order);
        (void)quadrature::get_reference_rule<T>(cell::type::interval, order);
    }
    catch (const std::exception&)
    {
        return false;
    }
    return true;
}

template <std::floating_point T>
bool reference_rule_available(cell::type cell_type, int order)
{
    if (order <= 0)
        return false;
    try
    {
        (void)quadrature::get_reference_rule<T>(cell_type, order);
    }
    catch (const std::exception&)
    {
        return false;
    }
    return true;
}

template <std::floating_point T>
int highest_available_common_order(std::span<const cell::type> cell_types,
                                   int requested_order)
{
    for (int order = requested_order; order >= 1; --order)
    {
        bool ok = true;
        for (const auto cell_type : cell_types)
            ok = ok && reference_rule_available<T>(cell_type, order);
        if (ok)
            return order;
    }
    return 0;
}

template <std::floating_point T>
int fixed_inspection_order(std::span<const cell::type> cell_types,
                           const ImplicitQuadratureOptions& opts)
{
    return highest_available_common_order<T>(
        cell_types, std::max(1, opts.inspection_order));
}

template <std::floating_point T>
int lower_fixed_inspection_order(std::span<const cell::type> cell_types,
                                 int high_order)
{
    if (high_order <= 1)
        return 0;
    return highest_available_common_order<T>(cell_types, high_order - 2);
}

template <std::floating_point T>
T sum_weights(std::span<const T> weights)
{
    return std::accumulate(weights.begin(), weights.end(), T(0));
}

template <std::floating_point T, std::integral I>
bool build_tetra_height_recovery_weights(
    std::span<const T> ref_vertices,
    const HeightDirectionFrame<T>& frame,
    const std::vector<BaseTriangle2<T>>& base_cells,
    const LevelSetCell<T, I>& ls_cell,
    IntegrationDomain domain,
    const ImplicitQuadratureOptions& opts,
    int quadrature_order,
    std::vector<T>& points,
    std::vector<T>& weights)
{
    points.clear();
    weights.clear();

    ImplicitQuadratureOptions local_opts = opts;
    local_opts.order = quadrature_order;

    const auto base_rule =
        quadrature::get_reference_rule<T>(cell::type::triangle, quadrature_order);
    const auto radial_rule =
        quadrature::get_reference_rule<T>(cell::type::interval, quadrature_order);
    const int gdim = ls_cell.gdim;
    const T clip_tol = static_cast<T>(100) * static_cast<T>(opts.root_tol);

    std::array<T, 3> xi{};
    std::array<T, 3> grad{};
    std::array<T, 9> tangents_ref{};

    for (const auto& base_cell : base_cells)
    {
        const auto base_t0 = base_tangent_ref<T>(frame, base_cell, 0);
        const auto base_t1 = base_tangent_ref<T>(frame, base_cell, 1);

        for (int qb = 0; qb < base_rule._num_points; ++qb)
        {
            const T u = base_rule._points[static_cast<std::size_t>(2 * qb)];
            const T v = base_rule._points[static_cast<std::size_t>(2 * qb + 1)];
            const auto base_y = base_triangle_point<T>(base_cell, u, v);
            const auto y = std::span<const T>(base_y.data(), 2);
            const auto clip = clip_tetra_height_line<T>(
                ref_vertices,
                frame,
                y,
                clip_tol);
            if (!clip.valid || clip.r_hi - clip.r_lo <= clip_tol)
                continue;

            auto roots = height_line_roots_tetra<T, I>(
                ref_vertices,
                frame,
                ls_cell,
                y,
                local_opts);

            if (domain == IntegrationDomain::interface_surface)
            {
                for (const T root : roots)
                {
                    eval_height_point<T>(frame, y, root, std::span<T>(xi.data(), 3));
                    ls_cell.grad(
                        std::span<const T>(xi.data(), 3),
                        std::span<T>(grad.data(), 3));
                    const T den = dot3<T>(grad, frame.er);
                    const T grad_norm = norm3<T>(grad);
                    const T threshold = std::max(
                        static_cast<T>(opts.root_tol),
                        static_cast<T>(opts.min_height_transversality) * grad_norm);
                    if (std::abs(den) <= threshold)
                        continue;

                    const T drdu = -dot3<T>(grad, base_t0) / den;
                    const T drdv = -dot3<T>(grad, base_t1) / den;
                    for (int d = 0; d < 3; ++d)
                    {
                        tangents_ref[static_cast<std::size_t>(d)] =
                            base_t0[static_cast<std::size_t>(d)]
                          + frame.er[static_cast<std::size_t>(d)] * drdu;
                        tangents_ref[static_cast<std::size_t>(3 + d)] =
                            base_t1[static_cast<std::size_t>(d)]
                          + frame.er[static_cast<std::size_t>(d)] * drdv;
                    }

                    const T measure = physical_measure_from_ref_tangents<T>(
                        ls_cell.cell_type,
                        std::span<const T>(
                            ls_cell.parent_vertex_coords.data(),
                            ls_cell.parent_vertex_coords.size()),
                        gdim,
                        std::span<const T>(xi.data(), 3),
                        std::span<const T>(tangents_ref.data(), 6),
                        2);
                    points.insert(points.end(), xi.begin(), xi.begin() + 3);
                    weights.push_back(
                        base_rule._weights[static_cast<std::size_t>(qb)] * measure);
                }
                continue;
            }

            std::vector<T> r_points;
            r_points.reserve(roots.size() + 2);
            r_points.push_back(clip.r_lo);
            r_points.insert(r_points.end(), roots.begin(), roots.end());
            r_points.push_back(clip.r_hi);
            sort_unique_with_tolerance<T>(
                r_points, static_cast<T>(100) * static_cast<T>(opts.root_tol));

            for (std::size_t j = 0; j + 1 < r_points.size(); ++j)
            {
                const T r0 = r_points[j];
                const T r1 = r_points[j + 1];
                const T radial_len = r1 - r0;
                if (radial_len <= static_cast<T>(opts.root_tol))
                    continue;

                const T r_mid = T(0.5) * (r0 + r1);
                const T f_mid = eval_height_phi<T, I>(frame, ls_cell, y, r_mid);
                const bool accept =
                    (domain == IntegrationDomain::negative_volume
                     && f_mid < static_cast<T>(opts.sign_tol))
                 || (domain == IntegrationDomain::positive_volume
                     && f_mid > -static_cast<T>(opts.sign_tol));
                if (!accept)
                    continue;

                for (int qr = 0; qr < radial_rule._num_points; ++qr)
                {
                    const T w = radial_rule._points[static_cast<std::size_t>(qr)];
                    const T r = r0 + radial_len * w;
                    eval_height_point<T>(frame, y, r, std::span<T>(xi.data(), 3));
                    for (int d = 0; d < 3; ++d)
                    {
                        tangents_ref[static_cast<std::size_t>(d)] =
                            base_t0[static_cast<std::size_t>(d)];
                        tangents_ref[static_cast<std::size_t>(3 + d)] =
                            base_t1[static_cast<std::size_t>(d)];
                        tangents_ref[static_cast<std::size_t>(6 + d)] =
                            frame.er[static_cast<std::size_t>(d)] * radial_len;
                    }

                    const T measure = physical_measure_from_ref_tangents<T>(
                        ls_cell.cell_type,
                        std::span<const T>(
                            ls_cell.parent_vertex_coords.data(),
                            ls_cell.parent_vertex_coords.size()),
                        gdim,
                        std::span<const T>(xi.data(), 3),
                        std::span<const T>(tangents_ref.data(), 9),
                        3);
                    points.insert(points.end(), xi.begin(), xi.begin() + 3);
                    weights.push_back(
                        base_rule._weights[static_cast<std::size_t>(qb)]
                      * radial_rule._weights[static_cast<std::size_t>(qr)]
                      * measure);
                }
            }
        }
    }

    return !weights.empty();
}

template <std::floating_point T, std::integral I>
bool append_tetra_height_recovery_quadrature(
    quadrature::QuadratureRules<T>& rules,
    const AdaptCell<T>& ac,
    int cell_id,
    const LevelSetCell<T, I>& ls_cell,
    IntegrationDomain domain,
    const ImplicitQuadratureOptions& opts,
    int parent_cell_id,
    std::int64_t* plan_component_hash = nullptr)
{
    if (ac.tdim != 3)
        return false;
    const auto leaf_type = ac.entity_types[ac.tdim][static_cast<std::size_t>(cell_id)];
    if (leaf_type != cell::type::tetrahedron)
        return false;
    if (ls_cell.parent_vertex_coords.empty())
        return false;

    std::vector<T> ref_vertices;
    gather_leaf_ref_vertices(ac, cell_id, ref_vertices);

    auto frames = enumerate_tetra_height_frames<T, I>(
            std::span<const T>(ref_vertices.data(), ref_vertices.size()),
            ls_cell);
    if (frames.empty())
        return false;

    const int tdim = 3;

    for (const auto& frame : frames)
    {
        auto base_cells = build_projected_tetra_base_cells<T>(
            std::span<const T>(ref_vertices.data(), ref_vertices.size()),
            frame,
            static_cast<T>(100) * static_cast<T>(opts.root_tol));
        if (base_cells.empty())
            continue;
        base_cells = split_base_cells_by_endpoint_functions<T, I>(
            std::move(base_cells),
            std::span<const T>(ref_vertices.data(), ref_vertices.size()),
            frame,
            ls_cell,
            opts);
        if (base_cells.empty())
            continue;
        auto split_result = split_projected_tetra_base_cells<T, I>(
            std::move(base_cells),
            std::span<const T>(ref_vertices.data(), ref_vertices.size()),
            frame,
            ls_cell,
            opts);
        if (!split_result.success || split_result.cells.empty())
            continue;

        if (!projected_base_cells_cover_tetra_leaf<T>(
                std::span<const T>(ref_vertices.data(), ref_vertices.size()),
                frame,
                split_result.cells,
                opts))
        {
            continue;
        }

        std::vector<T> points;
        std::vector<T> weights;
        if (!build_tetra_height_recovery_weights<T, I>(
                std::span<const T>(ref_vertices.data(), ref_vertices.size()),
                frame,
                split_result.cells,
                ls_cell,
                domain,
                opts,
                opts.order,
                points,
                weights))
        {
            continue;
        }
        if (weights.empty())
            continue;

        if (opts.enable_q_error_estimator)
        {
            const std::array<cell::type, 2> check_types{
                cell::type::triangle, cell::type::interval};
            const int high_order = fixed_inspection_order<T>(
                std::span<const cell::type>(check_types.data(), check_types.size()),
                opts);
            const int low_order = lower_fixed_inspection_order<T>(
                std::span<const cell::type>(check_types.data(), check_types.size()),
                high_order);

            if (high_order != 0 && low_order != 0)
            {
                std::vector<T> low_points;
                std::vector<T> low_weights;
                std::vector<T> high_points;
                std::vector<T> high_weights;
                if (!build_tetra_height_recovery_weights<T, I>(
                        std::span<const T>(ref_vertices.data(), ref_vertices.size()),
                        frame,
                        split_result.cells,
                        ls_cell,
                        domain,
                        opts,
                        low_order,
                        low_points,
                        low_weights))
                {
                    continue;
                }
                if (!build_tetra_height_recovery_weights<T, I>(
                        std::span<const T>(ref_vertices.data(), ref_vertices.size()),
                        frame,
                        split_result.cells,
                        ls_cell,
                        domain,
                        opts,
                        high_order,
                        high_points,
                        high_weights))
                {
                    continue;
                }
                const T low_sum = sum_weights<T>(
                    std::span<const T>(low_weights.data(), low_weights.size()));
                const T high_sum = sum_weights<T>(
                    std::span<const T>(high_weights.data(), high_weights.size()));
                const T scale = std::max({T(1), std::abs(low_sum), std::abs(high_sum)});
                const T tolerance =
                    std::max(static_cast<T>(opts.max_geometry_q_error),
                             static_cast<T>(1e-8))
                    * scale;
                if (std::abs(low_sum - high_sum) > tolerance)
                    continue;
            }
        }

        append_rule_points(
            rules,
            std::span<const T>(points.data(), points.size()),
            std::span<const T>(weights.data(), weights.size()),
            tdim,
            parent_cell_id);
        if (plan_component_hash != nullptr)
            *plan_component_hash =
                tetra_height_plan_component<T>(frame, split_result.cells);
        return true;
    }
    return false;
}

template <std::floating_point T, std::integral I>
void append_chart_quadrature(quadrature::QuadratureRules<T>& rules,
                             const SimplexChart<T>& chart,
                             const LevelSetCell<T, I>& ls_cell,
                             IntegrationDomain domain,
                             const ImplicitQuadratureOptions& opts,
                             int parent_cell_id)
{
    const int tdim = chart.tdim;
    const int gdim = ls_cell.gdim;
    if (ls_cell.parent_vertex_coords.empty())
        throw std::runtime_error(
            "implicit quadrature requires LevelSetCell parent physical coordinates");

    const std::array<cell::type, 1> inspection_types{chart.base_type};
    const int inspection_order = fixed_inspection_order<T>(
        std::span<const cell::type>(inspection_types.data(), inspection_types.size()),
        opts);
    if (inspection_order <= 0)
        throw std::runtime_error("implicit quadrature: no inspection rule available");
    const auto inspection_rule = quadrature::get_reference_rule<T>(
        chart.base_type, inspection_order);
    std::array<T, 3> ycheck{};
    for (int qb = 0; qb < inspection_rule._num_points; ++qb)
    {
        for (int a = 0; a < chart.base_dim; ++a)
        {
            ycheck[static_cast<std::size_t>(a)] =
                inspection_rule._points[static_cast<std::size_t>(
                    qb * chart.base_dim + a)];
        }
        const auto y = std::span<const T>(
            ycheck.data(), static_cast<std::size_t>(chart.base_dim));
        const auto root = solve_chart_root<T, I>(chart, ls_cell, y, opts);
        if (!root.valid)
        {
            throw std::runtime_error(
                std::string("implicit quadrature: fixed inspection rejected chart ray: ")
                + chart_root_failure_message(root.failure_code));
        }
    }

    const auto base_rule = quadrature::get_reference_rule<T>(
        chart.base_type, opts.order);
    const auto radial_rule = quadrature::get_reference_rule<T>(
        cell::type::interval, opts.order);

    std::vector<T> points;
    std::vector<T> weights;
    std::array<T, 3> ybuf{};
    std::array<T, 3> xi{};
    std::array<T, 3> cr{};
    std::vector<T> cy;
    std::vector<T> tangents_ref;

    for (int qb = 0; qb < base_rule._num_points; ++qb)
    {
        for (int a = 0; a < chart.base_dim; ++a)
        {
            ybuf[static_cast<std::size_t>(a)] =
                base_rule._points[static_cast<std::size_t>(
                    qb * chart.base_dim + a)];
        }
        const auto y = std::span<const T>(
            ybuf.data(), static_cast<std::size_t>(chart.base_dim));
        const auto root = solve_chart_root<T, I>(chart, ls_cell, y, opts);
        if (!root.valid)
        {
            throw std::runtime_error(
                std::string("implicit quadrature: failed to certify or solve a chart ray: ")
                + chart_root_failure_message(root.failure_code));
        }

        if (domain == IntegrationDomain::interface_surface)
        {
            chart.derivatives(
                root.r,
                y,
                std::span<T>(cr.data(), static_cast<std::size_t>(tdim)),
                cy);

            tangents_ref.assign(
                static_cast<std::size_t>(chart.base_dim * tdim), T(0));
            for (int a = 0; a < chart.base_dim; ++a)
            {
                const T drdy =
                    -dot<T>(
                        std::span<const T>(root.grad.data(), static_cast<std::size_t>(tdim)),
                        std::span<const T>(
                            cy.data() + static_cast<std::size_t>(a * tdim),
                            static_cast<std::size_t>(tdim)))
                    / root.den;
                for (int d = 0; d < tdim; ++d)
                {
                    tangents_ref[static_cast<std::size_t>(a * tdim + d)] =
                        cy[static_cast<std::size_t>(a * tdim + d)]
                      + cr[static_cast<std::size_t>(d)] * drdy;
                }
            }

            const T measure = physical_measure_from_ref_tangents<T>(
                ls_cell.cell_type,
                std::span<const T>(
                    ls_cell.parent_vertex_coords.data(),
                    ls_cell.parent_vertex_coords.size()),
                gdim,
                std::span<const T>(root.xi.data(), static_cast<std::size_t>(tdim)),
                std::span<const T>(tangents_ref.data(), tangents_ref.size()),
                chart.base_dim);
            points.insert(points.end(), root.xi.begin(), root.xi.begin() + tdim);
            weights.push_back(base_rule._weights[static_cast<std::size_t>(qb)] * measure);
            continue;
        }

        const bool negative = domain == IntegrationDomain::negative_volume;
        const T r0 = negative ? T(0) : root.r;
        const T r1 = negative ? root.r : T(1);
        const T length = r1 - r0;
        if (length <= static_cast<T>(opts.root_tol))
            continue;

        for (int qr = 0; qr < radial_rule._num_points; ++qr)
        {
            const T sr = radial_rule._points[static_cast<std::size_t>(qr)];
            const T r = r0 + length * sr;
            chart.eval(
                r,
                y,
                std::span<T>(xi.data(), static_cast<std::size_t>(tdim)));
            chart.derivatives(
                r,
                y,
                std::span<T>(cr.data(), static_cast<std::size_t>(tdim)),
                cy);

            tangents_ref.assign(static_cast<std::size_t>(tdim * tdim), T(0));
            for (int a = 0; a < chart.base_dim; ++a)
            {
                for (int d = 0; d < tdim; ++d)
                {
                    tangents_ref[static_cast<std::size_t>(a * tdim + d)] =
                        cy[static_cast<std::size_t>(a * tdim + d)];
                }
            }
            for (int d = 0; d < tdim; ++d)
            {
                tangents_ref[static_cast<std::size_t>(chart.base_dim * tdim + d)] =
                    cr[static_cast<std::size_t>(d)];
            }

            const T measure = physical_measure_from_ref_tangents<T>(
                ls_cell.cell_type,
                std::span<const T>(
                    ls_cell.parent_vertex_coords.data(),
                    ls_cell.parent_vertex_coords.size()),
                gdim,
                std::span<const T>(xi.data(), static_cast<std::size_t>(tdim)),
                std::span<const T>(tangents_ref.data(), tangents_ref.size()),
                tdim);
            points.insert(points.end(), xi.begin(), xi.begin() + tdim);
            weights.push_back(
                base_rule._weights[static_cast<std::size_t>(qb)]
              * radial_rule._weights[static_cast<std::size_t>(qr)]
              * length * measure);
        }
    }

    append_rule_points(
        rules,
        std::span<const T>(points.data(), points.size()),
        std::span<const T>(weights.data(), weights.size()),
        tdim,
        parent_cell_id);
}

template <std::floating_point T>
bool is_simplex(cell::type cell_type)
{
    return cell_type == cell::type::interval
        || cell_type == cell::type::triangle
        || cell_type == cell::type::tetrahedron;
}

template <std::floating_point T, std::integral I>
void append_full_leaf(quadrature::QuadratureRules<T>& rules,
                      const AdaptCell<T>& ac,
                      const LevelSetCell<T, I>& ls_cell,
                      int cell_id,
                      int order)
{
    const auto leaf_type = ac.entity_types[ac.tdim][static_cast<std::size_t>(cell_id)];
    if (!is_simplex<T>(leaf_type))
    {
        throw std::runtime_error(
            "implicit quadrature v1 supports full contributions only on simplex leaves");
    }

    std::vector<T> ref_vertices;
    gather_leaf_ref_vertices(ac, cell_id, ref_vertices);
    const auto physical_vertices = push_forward_points<T>(
        ls_cell.cell_type,
        std::span<const T>(
            ls_cell.parent_vertex_coords.data(),
            ls_cell.parent_vertex_coords.size()),
        ls_cell.gdim,
        std::span<const T>(ref_vertices.data(), ref_vertices.size()));

    append_simplex_entity_quadrature<T>(
        rules,
        leaf_type,
        std::span<const T>(ref_vertices.data(), ref_vertices.size()),
        std::span<const T>(physical_vertices.data(), physical_vertices.size()),
        ac.tdim,
        ls_cell.gdim,
        ac.parent_cell_id,
        order);
}

template <std::floating_point T>
quadrature::QuadratureRules<T> make_empty_rules(int tdim)
{
    quadrature::QuadratureRules<T> rules;
    rules._tdim = tdim;
    rules._offset.push_back(0);
    return rules;
}

template <std::floating_point T>
void merge_rules(quadrature::QuadratureRules<T>& out,
                 const quadrature::QuadratureRules<T>& in)
{
    if (in._parent_map.empty())
        return;
    if (out._tdim == 0)
        out._tdim = in._tdim;
    if (out._offset.empty())
        out._offset.push_back(0);
    if (out._tdim != in._tdim)
        throw std::runtime_error("implicit quadrature: cannot merge rules with different dimensions");

    for (std::size_t i = 0; i < in._parent_map.size(); ++i)
    {
        const int begin = in._offset[i];
        const int end = in._offset[i + 1];
        out._points.insert(
            out._points.end(),
            in._points.begin() + static_cast<std::ptrdiff_t>(begin * in._tdim),
            in._points.begin() + static_cast<std::ptrdiff_t>(end * in._tdim));
        out._weights.insert(
            out._weights.end(),
            in._weights.begin() + begin,
            in._weights.begin() + end);
        out._parent_map.push_back(in._parent_map[i]);
        out._debug_local_cell_id.push_back(
            i < in._debug_local_cell_id.size() ? in._debug_local_cell_id[i]
                                               : std::int32_t(-1));
        out._debug_chart_path.push_back(
            i < in._debug_chart_path.size() ? in._debug_chart_path[i]
                                            : std::int32_t(0));
        out._debug_refinement_depth.push_back(
            i < in._debug_refinement_depth.size()
                ? in._debug_refinement_depth[i]
                : std::int32_t(0));
        out._debug_chart_plan_hash.push_back(
            i < in._debug_chart_plan_hash.size()
                ? in._debug_chart_plan_hash[i]
                : std::int64_t(0));
        out._debug_candidate_mask.push_back(
            i < in._debug_candidate_mask.size()
                ? in._debug_candidate_mask[i]
                : std::int32_t(0));
        out._debug_rejection_reason.push_back(
            i < in._debug_rejection_reason.size()
                ? in._debug_rejection_reason[i]
                : std::int32_t(0));
        out._debug_measure_probe.push_back(
            i < in._debug_measure_probe.size() ? in._debug_measure_probe[i] : T(0));
        out._debug_validation_weight_sum.push_back(
            i < in._debug_validation_weight_sum.size()
                ? in._debug_validation_weight_sum[i]
                : T(0));
        out._offset.push_back(static_cast<std::int32_t>(out._weights.size()));
    }
}

template <std::floating_point T, std::integral I>
const LevelSetCell<T, I>& find_level_set_cell(const HOCutCells<T, I>& cut_cells,
                                              int cut_cell_id,
                                              int level_set_id)
{
    const int begin = cut_cells.ls_offsets[static_cast<std::size_t>(cut_cell_id)];
    const int end = cut_cells.ls_offsets[static_cast<std::size_t>(cut_cell_id + 1)];
    for (int i = begin; i < end; ++i)
    {
        const auto& ls_cell = cut_cells.level_set_cells[static_cast<std::size_t>(i)];
        if (ls_cell.level_set_id == level_set_id)
            return ls_cell;
    }
    throw std::runtime_error("implicit quadrature: selected level set is not active on cut cell");
}

template <std::floating_point T, std::integral I>
IntegrationDomain domain_from_relation(Relation relation)
{
    switch (relation)
    {
    case Relation::LessThan:
        return IntegrationDomain::negative_volume;
    case Relation::GreaterThan:
        return IntegrationDomain::positive_volume;
    case Relation::EqualTo:
        return IntegrationDomain::interface_surface;
    }
    throw std::runtime_error("implicit quadrature: invalid relation");
}

} // namespace

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> make_implicit_quadrature_impl(
    const AdaptCell<T>& base_adapt_cell,
    const LevelSetCell<T, I>& ls_cell,
    int level_set_id,
    IntegrationDomain domain,
    const ImplicitQuadratureOptions& opts,
    int chart_refinement_depth)
{
    auto rules = make_empty_rules<T>(base_adapt_cell.tdim);

    if (level_set_id < 0 || level_set_id >= 64)
        throw std::runtime_error("implicit quadrature: invalid level set id");
    if (base_adapt_cell.tdim != 2 && base_adapt_cell.tdim != 3)
        throw std::runtime_error("implicit quadrature v1 supports only triangles and tetrahedra");
    if (ls_cell.parent_vertex_coords.empty())
        throw std::runtime_error(
            "implicit quadrature requires polynomial LevelSetCell parent geometry");

    const int n_cells = base_adapt_cell.n_entities(base_adapt_cell.tdim);
    for (int c = 0; c < n_cells; ++c)
    {
        const auto verts = base_adapt_cell.entity_to_vertex[base_adapt_cell.tdim][
            static_cast<std::int32_t>(c)];

        const bool all_negative =
            leaf_has_uniform_sign(base_adapt_cell, verts, level_set_id, -1);
        const bool all_positive =
            leaf_has_uniform_sign(base_adapt_cell, verts, level_set_id, 1);

        if (domain == IntegrationDomain::negative_volume && all_negative)
        {
            append_full_leaf(rules, base_adapt_cell, ls_cell, c, opts.order);
            continue;
        }
        if (domain == IntegrationDomain::positive_volume && all_positive)
        {
            append_full_leaf(rules, base_adapt_cell, ls_cell, c, opts.order);
            continue;
        }
        if (all_negative || all_positive)
            continue;

        try
        {
            auto leaf_rules = make_empty_rules<T>(base_adapt_cell.tdim);
            std::int32_t sign_rejection_reason = debug_rejection_none;
            const auto chart = make_chart_from_leaf(base_adapt_cell, c, level_set_id);
            append_chart_quadrature(leaf_rules, chart, ls_cell, domain, opts,
                                    base_adapt_cell.parent_cell_id);
            auto sanity_rules = leaf_rules;
            const int diagnostic_refinement_limit = std::min(
                opts.max_refinement_iterations,
                std::max(0, opts.max_diagnostic_refinement_depth));
            const bool diagnostic_refinement_allowed =
                chart_refinement_depth < diagnostic_refinement_limit;
            if (domain == IntegrationDomain::interface_surface)
            {
                const std::array<cell::type, 1> check_types{chart.base_type};
                const int check_order = fixed_inspection_order<T>(
                    std::span<const cell::type>(check_types.data(), check_types.size()),
                    opts);
                if (check_order > 0 && check_order != opts.order)
                {
                    auto check_opts = opts;
                    check_opts.order = check_order;
                    sanity_rules = make_empty_rules<T>(base_adapt_cell.tdim);
                    append_chart_quadrature(
                        sanity_rules,
                        chart,
                        ls_cell,
                        domain,
                        check_opts,
                        base_adapt_cell.parent_cell_id);
                }
            }
            if (!surface_rule_matches_linearized_leaf<T, I>(
                    sanity_rules,
                    base_adapt_cell,
                    c,
                    ls_cell,
                    level_set_id,
                    domain,
                    opts,
                    static_cast<T>(opts.min_surface_area_ratio),
                    static_cast<T>(opts.max_surface_area_ratio)))
            {
                sign_rejection_reason |= debug_rejection_area_mismatch;
                if (diagnostic_refinement_allowed)
                {
                    throw std::runtime_error(
                        "implicit quadrature: surface rule failed linearized-measure sanity check");
                }
            }
            if (domain == IntegrationDomain::interface_surface
                && opts.enable_q_error_estimator)
            {
                const std::array<cell::type, 1> check_types{chart.base_type};
                const int high_order = fixed_inspection_order<T>(
                    std::span<const cell::type>(
                        check_types.data(), check_types.size()),
                    opts);
                const int low_order = lower_fixed_inspection_order<T>(
                    std::span<const cell::type>(
                        check_types.data(), check_types.size()),
                    high_order);
                if (high_order != 0 && low_order != 0)
                {
                    auto low_opts = opts;
                    low_opts.order = low_order;
                    auto low_rules = make_empty_rules<T>(base_adapt_cell.tdim);
                    append_chart_quadrature(
                        low_rules,
                        chart,
                        ls_cell,
                        domain,
                        low_opts,
                        base_adapt_cell.parent_cell_id);
                    const T low_sum = sum_weights<T>(
                        std::span<const T>(
                            low_rules._weights.data(), low_rules._weights.size()));
                    const T high_sum = sum_weights<T>(
                        std::span<const T>(
                            sanity_rules._weights.data(), sanity_rules._weights.size()));
                    const T scale = std::max(
                        static_cast<T>(opts.root_tol),
                        std::max(std::abs(low_sum), std::abs(high_sum)));
                    const T tolerance =
                        static_cast<T>(opts.max_surface_q_rel_error) * scale
                      + static_cast<T>(1000) * static_cast<T>(opts.root_tol);
                    if (std::abs(high_sum - low_sum) > tolerance)
                    {
                        sign_rejection_reason |= debug_rejection_q_error;
                        if (diagnostic_refinement_allowed)
                        {
                            throw std::runtime_error(
                                "implicit quadrature: surface chart failed fixed-order q agreement");
                        }
                    }
                }
            }
            const T measure_probe =
                domain == IntegrationDomain::interface_surface
                    ? linearized_interface_measure_for_leaf<T, I>(
                        base_adapt_cell, c, ls_cell, level_set_id, opts)
                    : T(0);
            const T validation_weight_sum =
                domain == IntegrationDomain::interface_surface
                    ? sum_weights<T>(
                        std::span<const T>(
                            sanity_rules._weights.data(), sanity_rules._weights.size()))
                    : sum_weights<T>(
                        std::span<const T>(
                            leaf_rules._weights.data(), leaf_rules._weights.size()));
            const auto plan_hash = simplex_chart_plan_hash<T>(
                chart,
                base_adapt_cell.parent_cell_id,
                c,
                chart_refinement_depth);
            const std::int32_t candidate_mask =
                debug_candidate_duffy
              | (chart_refinement_depth > 0 ? debug_candidate_refined : 0);
            annotate_rule_debug<T>(
                leaf_rules,
                c,
                1,
                chart_refinement_depth,
                plan_hash,
                candidate_mask,
                sign_rejection_reason,
                measure_probe,
                validation_weight_sum);
            merge_rules(rules, leaf_rules);
        }
        catch (const std::runtime_error& error)
        {
            std::int32_t failed_rejection_reason = debug_rejection_none;
            const std::string error_message = error.what();
            if (error_message.find("fixed inspection rejected chart ray")
                    != std::string::npos
                || error_message.find("failed to certify or solve a chart ray")
                    != std::string::npos)
            {
                failed_rejection_reason |= debug_rejection_chart_line;
            }
            if (error_message.find("endpoint_orientation") != std::string::npos)
                failed_rejection_reason |= debug_rejection_chart_endpoint;
            if (error_message.find("root_structure") != std::string::npos)
                failed_rejection_reason |= debug_rejection_chart_root_structure;
            if (error_message.find("root_solve") != std::string::npos)
                failed_rejection_reason |= debug_rejection_chart_root_solve;
            if (error_message.find("transversality") != std::string::npos)
                failed_rejection_reason |= debug_rejection_chart_transversality;
            if (error_message.find("unsupported triangle sign split")
                    != std::string::npos
                || error_message.find("unsupported tetrahedron sign split")
                    != std::string::npos
                || error_message.find("zero vertices in chart leaves")
                    != std::string::npos)
            {
                failed_rejection_reason |= debug_rejection_topology;
            }
            if (error_message.find("linearized-measure sanity check")
                    != std::string::npos)
            {
                failed_rejection_reason |= debug_rejection_area_mismatch;
            }
            if (error_message.find("fixed-order q agreement")
                    != std::string::npos)
            {
                failed_rejection_reason |= debug_rejection_q_error;
            }

            bool recovered = false;
            std::int64_t recovery_plan_component_hash = 0;
            std::int32_t recovery_rejection_reason = debug_rejection_none;
            T recovery_measure_probe = T(0);
            T recovery_validation_weight_sum = T(0);
            if (opts.enable_height_direction_recovery)
            {
                auto recovery_rules = make_empty_rules<T>(base_adapt_cell.tdim);
                recovered =
                    (base_adapt_cell.tdim == 2
                     && append_triangle_height_recovery_quadrature<T, I>(
                         recovery_rules,
                         base_adapt_cell,
                         c,
                         ls_cell,
                         domain,
                         opts,
                         base_adapt_cell.parent_cell_id,
                         &recovery_plan_component_hash))
                    || (base_adapt_cell.tdim == 3
                        && append_tetra_height_recovery_quadrature<T, I>(
                            recovery_rules,
                            base_adapt_cell,
                            c,
                            ls_cell,
                            domain,
                            opts,
                            base_adapt_cell.parent_cell_id,
                            &recovery_plan_component_hash));
                if (recovered
                    && domain == IntegrationDomain::interface_surface)
                {
                    auto recovery_sanity_rules = recovery_rules;
                    const std::array<cell::type, 2> check_types =
                        base_adapt_cell.tdim == 3
                            ? std::array<cell::type, 2>{
                                cell::type::triangle, cell::type::interval}
                            : std::array<cell::type, 2>{
                                cell::type::interval, cell::type::interval};
                    const int check_order = fixed_inspection_order<T>(
                        std::span<const cell::type>(check_types.data(), check_types.size()),
                        opts);
                    if (check_order > 0 && check_order != opts.order)
                    {
                        auto check_opts = opts;
                        check_opts.order = check_order;
                        recovery_sanity_rules = make_empty_rules<T>(base_adapt_cell.tdim);
                        recovered =
                            (base_adapt_cell.tdim == 2
                             && append_triangle_height_recovery_quadrature<T, I>(
                                     recovery_sanity_rules,
                                 base_adapt_cell,
                                 c,
                                 ls_cell,
                                 domain,
                                 check_opts,
                                 base_adapt_cell.parent_cell_id))
                             || (base_adapt_cell.tdim == 3
                                 && append_tetra_height_recovery_quadrature<T, I>(
                                    recovery_sanity_rules,
                                    base_adapt_cell,
                                    c,
                                    ls_cell,
                                    domain,
                                    check_opts,
                                    base_adapt_cell.parent_cell_id));
                    }
                    recovery_measure_probe = linearized_interface_measure_for_leaf<T, I>(
                        base_adapt_cell, c, ls_cell, level_set_id, opts);
                    recovery_validation_weight_sum = sum_weights<T>(
                        std::span<const T>(
                            recovery_sanity_rules._weights.data(),
                            recovery_sanity_rules._weights.size()));
                    if (recovered
                        && !surface_rule_matches_linearized_leaf<T, I>(
                            recovery_sanity_rules,
                            base_adapt_cell,
                            c,
                            ls_cell,
                            level_set_id,
                            domain,
                            opts,
                            static_cast<T>(opts.min_height_area_ratio),
                            static_cast<T>(opts.max_height_area_ratio)))
                    {
                        recovery_rejection_reason = debug_rejection_area_mismatch;
                        recovered = false;
                    }
                }
                else if (recovered
                         && !surface_rule_matches_linearized_leaf<T, I>(
                             recovery_rules,
                             base_adapt_cell,
                             c,
                             ls_cell,
                             level_set_id,
                             domain,
                             opts,
                             static_cast<T>(opts.min_height_area_ratio),
                             static_cast<T>(opts.max_height_area_ratio)))
                {
                    recovery_rejection_reason = debug_rejection_area_mismatch;
                    recovered = false;
                }

                if (recovered)
                {
                    if (domain != IntegrationDomain::interface_surface)
                    {
                        recovery_validation_weight_sum = sum_weights<T>(
                            std::span<const T>(
                                recovery_rules._weights.data(),
                                recovery_rules._weights.size()));
                    }
                    const auto plan_hash = accepted_height_plan_hash(
                        recovery_plan_component_hash,
                        base_adapt_cell.parent_cell_id,
                        c,
                        chart_refinement_depth);
                    const std::int32_t candidate_mask =
                        debug_candidate_height
                      | (chart_refinement_depth > 0 ? debug_candidate_refined : 0);
                    annotate_rule_debug<T>(
                        recovery_rules,
                        c,
                        2,
                        chart_refinement_depth,
                        plan_hash,
                        candidate_mask,
                        recovery_rejection_reason,
                        recovery_measure_probe,
                        recovery_validation_weight_sum);
                    merge_rules(rules, recovery_rules);
                }
                else
                    recovered = false;
            }
            if (!recovered && opts.enable_height_direction_recovery)
                failed_rejection_reason |= debug_rejection_height_failed
                                         | recovery_rejection_reason;

            if (!recovered)
            {
                if (opts.enable_chart_refinement
                    && chart_refinement_depth < opts.max_refinement_iterations)
                {
                    AdaptCell<T> refined = base_adapt_cell;
                    if (refine_ready_cell_on_largest_midpoint_value<T, I>(
                            refined, ls_cell, level_set_id, c))
                    {
                        fill_all_vertex_signs_from_level_set<T, I>(
                            refined,
                            ls_cell,
                            level_set_id,
                            static_cast<T>(opts.zero_tol));
                        certify_and_refine<T, I>(
                            refined,
                            ls_cell,
                            level_set_id,
                            std::max(1, opts.max_refinement_iterations - chart_refinement_depth - 1),
                            static_cast<T>(opts.zero_tol),
                            static_cast<T>(opts.sign_tol),
                            opts.edge_root_max_depth);
                        fill_all_vertex_signs_from_level_set<T, I>(
                            refined,
                            ls_cell,
                            level_set_id,
                            static_cast<T>(opts.zero_tol));
                        return make_implicit_quadrature_impl<T, I>(
                            refined,
                            ls_cell,
                            level_set_id,
                            domain,
                            opts,
                            chart_refinement_depth + 1);
                    }
                }
                if (chart_refinement_depth >= opts.max_refinement_iterations)
                    failed_rejection_reason |= debug_rejection_max_depth;
                if (domain == IntegrationDomain::interface_surface)
                {
                    auto fallback_rules = make_empty_rules<T>(base_adapt_cell.tdim);
                    std::vector<T> fallback_roots_ref;
                    if (append_straight_interface_fallback_quadrature<T, I>(
                            fallback_rules,
                            base_adapt_cell,
                            c,
                            ls_cell,
                            opts,
                            base_adapt_cell.parent_cell_id,
                            fallback_roots_ref))
                    {
                        const T measure_probe = linearized_interface_measure_for_leaf<T, I>(
                            base_adapt_cell, c, ls_cell, level_set_id, opts);
                        const T validation_weight_sum = sum_weights<T>(
                            std::span<const T>(
                                fallback_rules._weights.data(),
                                fallback_rules._weights.size()));
                        const auto plan_hash = straight_fallback_plan_hash<T>(
                            std::span<const T>(
                                fallback_roots_ref.data(), fallback_roots_ref.size()),
                            base_adapt_cell.tdim,
                            base_adapt_cell.parent_cell_id,
                            c,
                            chart_refinement_depth);
                        const std::int32_t candidate_mask =
                            debug_candidate_linear_fallback
                          | (chart_refinement_depth > 0 ? debug_candidate_refined : 0);
                        annotate_rule_debug<T>(
                            fallback_rules,
                            c,
                            3,
                            chart_refinement_depth,
                            plan_hash,
                            candidate_mask,
                            failed_rejection_reason | debug_rejection_hard_fallback,
                            measure_probe,
                            validation_weight_sum);
                        merge_rules(rules, fallback_rules);
                        continue;
                    }
                }
                throw std::runtime_error(
                    "implicit quadrature: failed on simplex leaf "
                    + std::to_string(c)
                    + " of parent cell "
                    + std::to_string(base_adapt_cell.parent_cell_id)
                    + "; sign chart error: "
                    + std::string(error.what()));
            }
        }
    }

    return rules;
}

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> make_implicit_quadrature(
    const AdaptCell<T>& base_adapt_cell,
    const LevelSetCell<T, I>& ls_cell,
    int level_set_id,
    IntegrationDomain domain,
    const ImplicitQuadratureOptions& opts)
{
    return make_implicit_quadrature_impl<T, I>(
        base_adapt_cell,
        ls_cell,
        level_set_id,
        domain,
        opts,
        0);
}

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> quadrature_rules_implicit(
    const HOMeshPart<T, I>& part,
    int order,
    bool include_uncut_cells,
    const ImplicitQuadratureOptions& input_opts)
{
    if (!part.cut_cells || !part.bg || !part.bg->mesh)
        throw std::runtime_error("HOMeshPart is not attached to cut-cell storage");
    if (part.expr.clauses.size() != 1)
    {
        throw std::runtime_error(
            "implicit_quadrature v1 supports exactly one level-set selection clause");
    }

    ImplicitQuadratureOptions opts = input_opts;
    opts.order = order;

    const auto& clause = part.expr.clauses.front();
    const int level_set_id = clause.level_set_index;
    const IntegrationDomain domain =
        domain_from_relation<T, I>(clause.relation);

    quadrature::QuadratureRules<T> rules;
    rules._tdim = part.bg->mesh->tdim;
    rules._offset.push_back(0);

    for (const std::int32_t cut_id32 : part.cut_cell_ids)
    {
        const int cut_id = static_cast<int>(cut_id32);
        const I parent_cell_id =
            part.cut_cells->parent_cell_ids[static_cast<std::size_t>(cut_id)];
        const auto& ls_cell =
            find_level_set_cell(*part.cut_cells, cut_id, level_set_id);

        AdaptCell<T> ac = make_adapt_cell(*part.bg->mesh, parent_cell_id);
        fill_all_vertex_signs_from_level_set(
            ac, ls_cell, level_set_id, static_cast<T>(opts.zero_tol));
        certify_and_refine(
            ac,
            ls_cell,
            level_set_id,
            opts.max_refinement_iterations,
            static_cast<T>(opts.zero_tol),
            static_cast<T>(opts.sign_tol),
            opts.edge_root_max_depth);
        fill_all_vertex_signs_from_level_set(
            ac, ls_cell, level_set_id, static_cast<T>(opts.zero_tol));
        build_edges(ac);
        if (ac.tdim == 3)
            build_faces(ac);
        recompute_active_level_set_masks(ac, level_set_id + 1);
        rebuild_zero_entity_inventory(ac);

        auto local = make_implicit_quadrature<T, I>(
            ac, ls_cell, level_set_id, domain, opts);
        merge_rules(rules, local);
    }

    if (include_uncut_cells && domain != IntegrationDomain::interface_surface)
    {
        const auto& mesh = *part.bg->mesh;
        for (const I cell_id : part.uncut_cell_ids)
        {
            const auto ctype = mesh.cell_type(cell_id);
            const auto physical_vertices = parent_cell_vertex_coords_vtk(mesh, cell_id);
            append_uncut_parent_cell<T>(
                rules,
                ctype,
                std::span<const T>(physical_vertices.data(), physical_vertices.size()),
                mesh.tdim,
                mesh.gdim,
                static_cast<int>(cell_id),
                opts.order);
        }
    }

    return rules;
}

template quadrature::QuadratureRules<double> make_implicit_quadrature(
    const AdaptCell<double>&,
    const LevelSetCell<double, int>&,
    int,
    IntegrationDomain,
    const ImplicitQuadratureOptions&);
template quadrature::QuadratureRules<float> make_implicit_quadrature(
    const AdaptCell<float>&,
    const LevelSetCell<float, int>&,
    int,
    IntegrationDomain,
    const ImplicitQuadratureOptions&);
template quadrature::QuadratureRules<double> make_implicit_quadrature(
    const AdaptCell<double>&,
    const LevelSetCell<double, long>&,
    int,
    IntegrationDomain,
    const ImplicitQuadratureOptions&);
template quadrature::QuadratureRules<float> make_implicit_quadrature(
    const AdaptCell<float>&,
    const LevelSetCell<float, long>&,
    int,
    IntegrationDomain,
    const ImplicitQuadratureOptions&);

template quadrature::QuadratureRules<double> quadrature_rules_implicit(
    const HOMeshPart<double, int>&,
    int,
    bool,
    const ImplicitQuadratureOptions&);
template quadrature::QuadratureRules<float> quadrature_rules_implicit(
    const HOMeshPart<float, int>&,
    int,
    bool,
    const ImplicitQuadratureOptions&);
template quadrature::QuadratureRules<double> quadrature_rules_implicit(
    const HOMeshPart<double, long>&,
    int,
    bool,
    const ImplicitQuadratureOptions&);
template quadrature::QuadratureRules<float> quadrature_rules_implicit(
    const HOMeshPart<float, long>&,
    int,
    bool,
    const ImplicitQuadratureOptions&);

} // namespace cutcells::implicit_quadrature
