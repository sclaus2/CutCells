// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include "geometric_quantity.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <span>
#include <string_view>
#include <vector>

namespace cutcells::graph_criteria
{

enum class HostDimension : std::uint8_t
{
    edge = 1,
    face = 2
};

enum class DirectionKind : std::uint8_t
{
    projected_straight_host_normal = 0,
    projected_level_set_gradient = 1
};

enum class FailureReason : std::uint8_t
{
    none = 0,
    invalid_input,
    invalid_parent_entity,
    host_point_outside_parent_entity,
    invalid_host_frame,
    weak_restricted_gradient,
    degenerate_direction,
    direction_not_admissible,
    tangent_to_zero_set,
    root_not_on_search_line,
    root_segment_leaves_parent_entity,
    excessive_correction_distance,
    excessive_tangential_shift,
    excessive_drift_amplification,
    direction_too_tangential_to_host,
    edge_ordering_fold,
    face_degenerate,
    surface_jacobian_not_positive,
    refinement_requested
};

inline std::string_view failure_reason_name(FailureReason reason)
{
    switch (reason)
    {
    case FailureReason::none:
        return "none";
    case FailureReason::invalid_input:
        return "invalid_input";
    case FailureReason::invalid_parent_entity:
        return "invalid_parent_entity";
    case FailureReason::host_point_outside_parent_entity:
        return "host_point_outside_parent_entity";
    case FailureReason::invalid_host_frame:
        return "invalid_host_frame";
    case FailureReason::weak_restricted_gradient:
        return "weak_restricted_gradient";
    case FailureReason::degenerate_direction:
        return "degenerate_direction";
    case FailureReason::direction_not_admissible:
        return "direction_not_admissible";
    case FailureReason::tangent_to_zero_set:
        return "tangent_to_zero_set";
    case FailureReason::root_not_on_search_line:
        return "root_not_on_search_line";
    case FailureReason::root_segment_leaves_parent_entity:
        return "root_segment_leaves_parent_entity";
    case FailureReason::excessive_correction_distance:
        return "excessive_correction_distance";
    case FailureReason::excessive_tangential_shift:
        return "excessive_tangential_shift";
    case FailureReason::excessive_drift_amplification:
        return "excessive_drift_amplification";
    case FailureReason::direction_too_tangential_to_host:
        return "direction_too_tangential_to_host";
    case FailureReason::edge_ordering_fold:
        return "edge_ordering_fold";
    case FailureReason::face_degenerate:
        return "face_degenerate";
    case FailureReason::surface_jacobian_not_positive:
        return "surface_jacobian_not_positive";
    case FailureReason::refinement_requested:
        return "refinement_requested";
    }
    return "unknown";
}

template <std::floating_point T>
struct Options
{
    T tolerance = geom::default_tolerance<T>();
    T min_restricted_gradient_strength = T(64) * std::numeric_limits<T>::epsilon();
    T min_transversality = T(1.0e-6);
    T min_host_normal_alignment = T(0.25);
    T max_drift_amplification = T(4);
    T max_relative_correction_distance = T(0.5);
    T max_relative_tangential_shift = T(0.25);
    T max_root_line_residual = std::sqrt(std::numeric_limits<T>::epsilon());
    T min_edge_ordering_fraction = std::sqrt(std::numeric_limits<T>::epsilon());
    T min_surface_jacobian_ratio = std::sqrt(std::numeric_limits<T>::epsilon());
};

template <std::floating_point T>
struct HostFrame
{
    HostDimension dimension = HostDimension::edge;
    std::vector<T> normal;
    std::vector<T> tangent;
    T h = T(1);
};

template <std::floating_point T>
struct Metrics
{
    T restricted_gradient_strength = T(0);
    T true_transversality = T(0);
    T host_normal_alignment = T(0);
    T drift_amplification = std::numeric_limits<T>::infinity();
    T root_alpha = T(0);
    T root_line_residual = T(0);
    bool root_segment_contained = false;
    T relative_correction_distance = T(0);
    T relative_tangential_shift = T(0);
};

template <std::floating_point T>
struct DirectionReport
{
    bool accepted = false;
    DirectionKind kind = DirectionKind::projected_straight_host_normal;
    FailureReason failure_reason = FailureReason::none;
    Metrics<T> metrics;
};

template <std::floating_point T>
struct SelectionReport
{
    bool accepted = false;
    bool request_refinement = false;
    DirectionKind selected_kind = DirectionKind::projected_straight_host_normal;
    FailureReason failure_reason = FailureReason::none;
    DirectionReport<T> straight_host_normal;
    DirectionReport<T> level_set_gradient;
    std::vector<T> selected_direction;
    std::vector<T> selected_root;
};

template <std::floating_point T>
struct EdgeOrderingReport
{
    bool accepted = false;
    FailureReason failure_reason = FailureReason::none;
    T minimum_gap_ratio = std::numeric_limits<T>::infinity();
};

template <std::floating_point T>
struct FaceQualityReport
{
    bool accepted = false;
    FailureReason failure_reason = FailureReason::none;
    T minimum_surface_jacobian_ratio = std::numeric_limits<T>::infinity();
    int failed_triangle_index = -1;
    T failed_surface_jacobian_ratio = std::numeric_limits<T>::quiet_NaN();
};

template <std::floating_point T>
inline std::span<const T> vector_span(const std::vector<T>& v)
{
    return std::span<const T>(v.data(), v.size());
}

template <std::floating_point T>
inline bool same_dimension(std::size_t n, std::span<const T> a)
{
    return a.size() == n;
}

template <std::floating_point T, typename... Spans>
inline bool same_dimension(std::size_t n, std::span<const T> a, Spans... rest)
{
    return a.size() == n && same_dimension<T>(n, rest...);
}

template <std::floating_point T>
inline T safe_scale(T a, T b, T c, T d)
{
    return std::max({T(1), std::fabs(a), std::fabs(b), std::fabs(c), std::fabs(d)});
}

template <std::floating_point T>
inline T absolute_cosine(std::span<const T> a,
                         std::span<const T> b,
                         T tol)
{
    const auto alignment = geom::alignment<T>(a, b, tol);
    if (alignment.degenerate)
        return T(0);
    return std::fabs(alignment.cosine);
}

template <std::floating_point T>
inline T drift_from_alignment(T abs_cosine, T tol)
{
    if (abs_cosine <= tol)
        return std::numeric_limits<T>::infinity();
    const T tangential2 = std::max(T(0), T(1) - abs_cosine * abs_cosine);
    return std::sqrt(tangential2) / abs_cosine;
}

template <std::floating_point T>
inline bool valid_host_frame(const HostFrame<T>& host,
                             std::size_t tdim,
                             T tol)
{
    if (!(host.h > tol))
        return false;
    if (host.normal.size() != tdim
        || geom::norm<T>(vector_span(host.normal)) <= tol)
    {
        return false;
    }
    if (host.dimension == HostDimension::edge)
    {
        if (host.tangent.size() != tdim
            || geom::norm<T>(vector_span(host.tangent)) <= tol)
        {
            return false;
        }
    }
    return host.dimension == HostDimension::edge
        || host.dimension == HostDimension::face;
}

template <std::floating_point T>
inline bool direction_admissible_in_parent_entity(
    cell::type parent_cell_type,
    geom::ParentEntity admissible_parent_entity,
    std::span<const T> x_h,
    std::span<const T> direction,
    T tol)
{
    if (!geom::point_in_parent_entity<T>(
            parent_cell_type, admissible_parent_entity, x_h, tol))
    {
        return false;
    }
    if (admissible_parent_entity.dim == 0
        && geom::norm<T>(direction) > tol)
    {
        return false;
    }

    const auto interval = geom::clip_line_interval_in_parent_entity<T>(
        parent_cell_type,
        admissible_parent_entity,
        x_h,
        direction,
        -std::numeric_limits<T>::infinity(),
        std::numeric_limits<T>::infinity(),
        tol);

    return interval.valid && interval.t0 <= tol && interval.t1 >= -tol;
}

template <std::floating_point T>
inline bool root_segment_contained_in_parent_entity(
    cell::type parent_cell_type,
    geom::ParentEntity admissible_parent_entity,
    std::span<const T> x_h,
    std::span<const T> direction,
    T alpha,
    T tol)
{
    const auto interval = geom::clip_line_interval_in_parent_entity<T>(
        parent_cell_type,
        admissible_parent_entity,
        x_h,
        direction,
        -std::numeric_limits<T>::infinity(),
        std::numeric_limits<T>::infinity(),
        tol);
    if (!interval.valid)
        return false;

    const T lo = std::min(T(0), alpha);
    const T hi = std::max(T(0), alpha);
    return lo >= interval.t0 - tol && hi <= interval.t1 + tol;
}

template <std::floating_point T>
inline T relative_tangential_shift(std::span<const T> delta,
                                   const HostFrame<T>& host,
                                   T tol)
{
    if (host.dimension == HostDimension::edge)
    {
        const T tangent_norm = geom::norm<T>(vector_span(host.tangent));
        if (tangent_norm <= tol)
            return std::numeric_limits<T>::infinity();
        return std::fabs(geom::dot<T>(delta, vector_span(host.tangent)))
             / (tangent_norm * host.h);
    }

    const auto tangential = geom::project_into_plane<T>(
        delta, vector_span(host.normal), tol);
    if (tangential.degenerate() && tangential.degeneracy == geom::Degeneracy::zero_frame)
        return std::numeric_limits<T>::infinity();
    return tangential.norm / host.h;
}

template <std::floating_point T>
inline DirectionReport<T> evaluate_direction(
    cell::type parent_cell_type,
    geom::ParentEntity admissible_parent_entity,
    const HostFrame<T>& host,
    std::span<const T> x_h,
    std::span<const T> raw_level_set_gradient,
    std::span<const T> candidate_direction,
    std::span<const T> x_c,
    DirectionKind kind,
    const Options<T>& options = {})
{
    DirectionReport<T> report;
    report.kind = kind;

    const int tdim = cell::get_tdim(parent_cell_type);
    if (tdim <= 0
        || !admissible_parent_entity.valid()
        || admissible_parent_entity.dim > tdim)
    {
        report.failure_reason = FailureReason::invalid_parent_entity;
        return report;
    }

    if (!same_dimension<T>(
            static_cast<std::size_t>(tdim),
            x_h,
            raw_level_set_gradient,
            candidate_direction,
            x_c))
    {
        report.failure_reason = FailureReason::invalid_input;
        return report;
    }

    const T tol = options.tolerance;
    if (!valid_host_frame<T>(host, static_cast<std::size_t>(tdim), tol))
    {
        report.failure_reason = FailureReason::invalid_host_frame;
        return report;
    }

    if (!geom::point_in_parent_entity<T>(
            parent_cell_type, admissible_parent_entity, x_h, tol))
    {
        report.failure_reason = FailureReason::host_point_outside_parent_entity;
        return report;
    }

    const auto restricted_gradient =
        geom::restricted_level_set_gradient_in_parent_frame<T>(
            parent_cell_type,
            admissible_parent_entity,
            raw_level_set_gradient,
            tol);
    report.metrics.restricted_gradient_strength = restricted_gradient.norm;
    if (restricted_gradient.degenerate()
        || restricted_gradient.norm < options.min_restricted_gradient_strength)
    {
        report.failure_reason = FailureReason::weak_restricted_gradient;
        return report;
    }

    const T direction_norm = geom::norm<T>(candidate_direction);
    if (direction_norm <= tol)
    {
        report.failure_reason = FailureReason::degenerate_direction;
        return report;
    }

    if (!direction_admissible_in_parent_entity<T>(
            parent_cell_type,
            admissible_parent_entity,
            x_h,
            candidate_direction,
            tol))
    {
        report.failure_reason = FailureReason::direction_not_admissible;
        return report;
    }

    report.metrics.true_transversality =
        std::fabs(geom::dot<T>(
            candidate_direction,
            vector_span(restricted_gradient.value)))
        / (direction_norm * restricted_gradient.norm);

    report.metrics.host_normal_alignment =
        absolute_cosine<T>(candidate_direction, vector_span(host.normal), tol);
    report.metrics.drift_amplification =
        drift_from_alignment<T>(report.metrics.host_normal_alignment, tol);

    const auto delta = geom::subtract<T>(x_c, x_h);
    const T delta_norm = geom::norm<T>(vector_span(delta));
    const T dd = direction_norm * direction_norm;
    report.metrics.root_alpha =
        geom::dot<T>(vector_span(delta), candidate_direction) / dd;

    std::vector<T> residual(delta.size(), T(0));
    for (std::size_t i = 0; i < delta.size(); ++i)
    {
        residual[i] = delta[i]
            - report.metrics.root_alpha * candidate_direction[i];
    }
    const T residual_norm = geom::norm<T>(vector_span(residual));
    const T residual_scale = safe_scale<T>(
        host.h,
        delta_norm,
        std::fabs(report.metrics.root_alpha) * direction_norm,
        direction_norm);
    report.metrics.root_line_residual = residual_norm / residual_scale;

    report.metrics.root_segment_contained =
        root_segment_contained_in_parent_entity<T>(
            parent_cell_type,
            admissible_parent_entity,
            x_h,
            candidate_direction,
            report.metrics.root_alpha,
            tol)
        && geom::point_in_parent_entity<T>(
            parent_cell_type,
            admissible_parent_entity,
            x_c,
            tol);

    report.metrics.relative_correction_distance = delta_norm / host.h;
    report.metrics.relative_tangential_shift =
        relative_tangential_shift<T>(vector_span(delta), host, tol);

    if (report.metrics.true_transversality < options.min_transversality)
    {
        report.failure_reason = FailureReason::tangent_to_zero_set;
        return report;
    }
    if (report.metrics.root_line_residual > options.max_root_line_residual)
    {
        report.failure_reason = FailureReason::root_not_on_search_line;
        return report;
    }
    if (!report.metrics.root_segment_contained)
    {
        report.failure_reason = FailureReason::root_segment_leaves_parent_entity;
        return report;
    }
    if (report.metrics.relative_correction_distance
        > options.max_relative_correction_distance)
    {
        report.failure_reason = FailureReason::excessive_correction_distance;
        return report;
    }
    if (report.metrics.relative_tangential_shift
        > options.max_relative_tangential_shift)
    {
        report.failure_reason = FailureReason::excessive_tangential_shift;
        return report;
    }
    if (report.metrics.drift_amplification > options.max_drift_amplification)
    {
        report.failure_reason = FailureReason::excessive_drift_amplification;
        return report;
    }
    if (report.metrics.host_normal_alignment < options.min_host_normal_alignment)
    {
        report.failure_reason = FailureReason::direction_too_tangential_to_host;
        return report;
    }

    report.accepted = true;
    report.failure_reason = FailureReason::none;
    return report;
}

template <std::floating_point T>
inline SelectionReport<T> select_preferred_direction(
    cell::type parent_cell_type,
    geom::ParentEntity admissible_parent_entity,
    const HostFrame<T>& host,
    std::span<const T> x_h,
    std::span<const T> raw_level_set_gradient,
    std::span<const T> projected_level_set_gradient_direction,
    std::span<const T> projected_level_set_gradient_root,
    std::span<const T> projected_straight_host_normal_direction,
    std::span<const T> projected_straight_host_normal_root,
    const Options<T>& options = {})
{
    SelectionReport<T> report;
    report.straight_host_normal = evaluate_direction<T>(
        parent_cell_type,
        admissible_parent_entity,
        host,
        x_h,
        raw_level_set_gradient,
        projected_straight_host_normal_direction,
        projected_straight_host_normal_root,
        DirectionKind::projected_straight_host_normal,
        options);
    report.level_set_gradient = evaluate_direction<T>(
        parent_cell_type,
        admissible_parent_entity,
        host,
        x_h,
        raw_level_set_gradient,
        projected_level_set_gradient_direction,
        projected_level_set_gradient_root,
        DirectionKind::projected_level_set_gradient,
        options);

    if (report.straight_host_normal.accepted)
    {
        report.accepted = true;
        report.selected_kind = DirectionKind::projected_straight_host_normal;
        report.failure_reason = FailureReason::none;
        report.selected_direction.assign(
            projected_straight_host_normal_direction.begin(),
            projected_straight_host_normal_direction.end());
        report.selected_root.assign(
            projected_straight_host_normal_root.begin(),
            projected_straight_host_normal_root.end());
        return report;
    }

    if (report.level_set_gradient.accepted
        && report.level_set_gradient.metrics.drift_amplification
            <= options.max_drift_amplification)
    {
        report.accepted = true;
        report.selected_kind = DirectionKind::projected_level_set_gradient;
        report.failure_reason = FailureReason::none;
        report.selected_direction.assign(
            projected_level_set_gradient_direction.begin(),
            projected_level_set_gradient_direction.end());
        report.selected_root.assign(
            projected_level_set_gradient_root.begin(),
            projected_level_set_gradient_root.end());
        return report;
    }

    report.accepted = false;
    report.request_refinement = true;
    report.failure_reason = FailureReason::refinement_requested;
    return report;
}

template <std::floating_point T>
inline EdgeOrderingReport<T> evaluate_projected_edge_ordering(
    std::span<const T> host_points,
    std::span<const T> corrected_points,
    int point_dim,
    std::span<const T> host_tangent,
    T h_edge,
    const Options<T>& options = {})
{
    EdgeOrderingReport<T> report;
    const T tol = options.tolerance;

    if (point_dim <= 0
        || host_points.size() != corrected_points.size()
        || host_points.size() % static_cast<std::size_t>(point_dim) != 0
        || host_tangent.size() != static_cast<std::size_t>(point_dim)
        || !(h_edge > tol))
    {
        report.failure_reason = FailureReason::invalid_input;
        return report;
    }

    const std::size_t n =
        host_points.size() / static_cast<std::size_t>(point_dim);
    if (n < 2 || geom::norm<T>(host_tangent) <= tol)
    {
        report.failure_reason = FailureReason::invalid_host_frame;
        return report;
    }

    std::vector<T> host_s(n, T(0));
    std::vector<T> corr_s(n, T(0));
    for (std::size_t i = 0; i < n; ++i)
    {
        const auto hp = host_points.subspan(
            i * static_cast<std::size_t>(point_dim),
            static_cast<std::size_t>(point_dim));
        const auto cp = corrected_points.subspan(
            i * static_cast<std::size_t>(point_dim),
            static_cast<std::size_t>(point_dim));
        host_s[i] = geom::dot<T>(hp, host_tangent);
        corr_s[i] = geom::dot<T>(cp, host_tangent);
    }

    const T orientation = (host_s.back() >= host_s.front()) ? T(1) : T(-1);
    report.minimum_gap_ratio = std::numeric_limits<T>::infinity();
    for (std::size_t i = 1; i < n; ++i)
    {
        const T host_gap = orientation * (host_s[i] - host_s[i - 1]);
        const T corr_gap = orientation * (corr_s[i] - corr_s[i - 1]);
        if (host_gap <= tol * h_edge)
        {
            report.failure_reason = FailureReason::invalid_host_frame;
            return report;
        }
        const T gap_ratio = corr_gap / host_gap;
        report.minimum_gap_ratio = std::min(report.minimum_gap_ratio, gap_ratio);
        if (gap_ratio <= options.min_edge_ordering_fraction)
        {
            report.failure_reason = FailureReason::edge_ordering_fold;
            return report;
        }
    }

    report.accepted = true;
    report.failure_reason = FailureReason::none;
    return report;
}

} // namespace cutcells::graph_criteria
