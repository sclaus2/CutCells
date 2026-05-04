// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include "geometric_quantity.h"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

namespace cutcells::tangent_shift
{

enum class FailureReason : std::uint8_t
{
    none = 0,
    invalid_input,
    invalid_parent_entity,
    invalid_curve,
    node_index_out_of_range,
    target_arclength_out_of_range,
    current_point_outside_parent_entity,
    correction_direction_not_admissible,
    correction_left_parent_entity,
    correction_failed,
    phi_residual_too_large,
    ordering_destroyed,
    spacing_too_small,
    arclength_error_too_large,
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
    case FailureReason::invalid_curve:
        return "invalid_curve";
    case FailureReason::node_index_out_of_range:
        return "node_index_out_of_range";
    case FailureReason::target_arclength_out_of_range:
        return "target_arclength_out_of_range";
    case FailureReason::current_point_outside_parent_entity:
        return "current_point_outside_parent_entity";
    case FailureReason::correction_direction_not_admissible:
        return "correction_direction_not_admissible";
    case FailureReason::correction_left_parent_entity:
        return "correction_left_parent_entity";
    case FailureReason::correction_failed:
        return "correction_failed";
    case FailureReason::phi_residual_too_large:
        return "phi_residual_too_large";
    case FailureReason::ordering_destroyed:
        return "ordering_destroyed";
    case FailureReason::spacing_too_small:
        return "spacing_too_small";
    case FailureReason::arclength_error_too_large:
        return "arclength_error_too_large";
    case FailureReason::refinement_requested:
        return "refinement_requested";
    }
    return "unknown";
}

template <std::floating_point T>
struct Options
{
    T tolerance = geom::default_tolerance<T>();
    T phi_tolerance = std::sqrt(std::numeric_limits<T>::epsilon());
    T min_direction_norm = T(64) * std::numeric_limits<T>::epsilon();

    T max_relative_arclength_drift = T(1.0e-3);
    T max_absolute_arclength_drift = std::sqrt(std::numeric_limits<T>::epsilon());
    T max_relative_final_arclength_error = T(1.0e-3);
    T max_absolute_final_arclength_error = std::sqrt(std::numeric_limits<T>::epsilon());

    T min_relative_spacing = std::sqrt(std::numeric_limits<T>::epsilon());
    T min_absolute_spacing = T(0);

    int max_correction_iterations = 16;
    int max_line_search_iterations = 16;
};

template <std::floating_point T>
struct Metrics
{
    T curve_length = T(0);
    T desired_arclength = T(0);
    T actual_arclength = T(0);
    T initial_arclength_error = T(0);
    T final_arclength = T(0);
    T final_arclength_error = T(0);
    T phi_residual = std::numeric_limits<T>::infinity();
    T minimum_spacing = std::numeric_limits<T>::infinity();
    T minimum_spacing_ratio = std::numeric_limits<T>::infinity();
};

template <std::floating_point T>
struct CorrectionReport
{
    bool accepted = false;
    bool left_parent_entity = false;
    FailureReason failure_reason = FailureReason::none;
    std::vector<T> point;
    T residual = std::numeric_limits<T>::infinity();
    int iterations = 0;
};

template <std::floating_point T>
struct ShiftReport
{
    bool accepted = false;
    bool corrected = false;
    bool request_refinement = false;
    FailureReason failure_reason = FailureReason::none;

    Metrics<T> metrics;
    int correction_iterations = 0;

    std::vector<T> shifted_seed;
    std::vector<T> corrected_node;
    std::vector<T> corrected_edge_points;
};

template <std::floating_point T>
struct ArclengthPoint
{
    bool valid = false;
    std::vector<T> point;
    std::size_t segment = 0;
    T segment_parameter = T(0);
};

template <std::floating_point T>
struct ClosestArclength
{
    bool valid = false;
    T arclength = T(0);
    T distance = std::numeric_limits<T>::infinity();
    std::size_t segment = 0;
    T segment_parameter = T(0);
};

template <std::floating_point T>
inline std::span<const T> vector_span(const std::vector<T>& v)
{
    return std::span<const T>(v.data(), v.size());
}

template <std::floating_point T>
inline std::span<const T> point_span(std::span<const T> points,
                                     int point_dim,
                                     std::size_t i)
{
    return points.subspan(
        i * static_cast<std::size_t>(point_dim),
        static_cast<std::size_t>(point_dim));
}

template <std::floating_point T>
inline std::span<T> mutable_point_span(std::span<T> points,
                                       int point_dim,
                                       std::size_t i)
{
    return points.subspan(
        i * static_cast<std::size_t>(point_dim),
        static_cast<std::size_t>(point_dim));
}

template <std::floating_point T>
inline bool flat_point_array_valid(std::span<const T> points, int point_dim)
{
    return point_dim > 0
        && !points.empty()
        && points.size() % static_cast<std::size_t>(point_dim) == 0;
}

template <std::floating_point T>
inline T point_distance(std::span<const T> a, std::span<const T> b)
{
    const auto delta = geom::subtract<T>(a, b);
    return geom::norm<T>(vector_span(delta));
}

template <std::floating_point T>
inline T arclength_tolerance(T length, T absolute_tolerance, T relative_tolerance)
{
    return absolute_tolerance + relative_tolerance * length;
}

template <std::floating_point T>
inline bool cumulative_arclength(std::span<const T> points,
                                 int point_dim,
                                 std::vector<T>& cumulative,
                                 T tol = geom::default_tolerance<T>())
{
    if (!flat_point_array_valid<T>(points, point_dim))
        return false;

    const std::size_t n = points.size() / static_cast<std::size_t>(point_dim);
    if (n < 2)
        return false;

    cumulative.assign(n, T(0));
    for (std::size_t i = 1; i < n; ++i)
    {
        const T ds = point_distance<T>(
            point_span<T>(points, point_dim, i - 1),
            point_span<T>(points, point_dim, i));
        if (!std::isfinite(ds))
            return false;
        cumulative[i] = cumulative[i - 1] + ds;
    }

    return cumulative.back() > tol;
}

template <std::floating_point T>
inline ArclengthPoint<T> point_at_arclength(std::span<const T> points,
                                            int point_dim,
                                            T target_arclength,
                                            T tol = geom::default_tolerance<T>())
{
    ArclengthPoint<T> out;
    std::vector<T> cumulative;
    if (!cumulative_arclength<T>(points, point_dim, cumulative, tol))
        return out;

    const T length = cumulative.back();
    if (target_arclength < -tol || target_arclength > length + tol)
        return out;

    const std::size_t n = cumulative.size();
    if (target_arclength <= tol)
    {
        out.valid = true;
        out.point.assign(points.begin(), points.begin() + point_dim);
        out.segment = 0;
        out.segment_parameter = T(0);
        return out;
    }
    if (length - target_arclength <= tol)
    {
        out.valid = true;
        out.point.assign(
            points.end() - static_cast<std::ptrdiff_t>(point_dim),
            points.end());
        out.segment = n - 2;
        out.segment_parameter = T(1);
        return out;
    }

    for (std::size_t i = 1; i < n; ++i)
    {
        if (target_arclength > cumulative[i] + tol)
            continue;

        const T segment_length = cumulative[i] - cumulative[i - 1];
        if (segment_length <= tol)
            continue;

        T u = (target_arclength - cumulative[i - 1]) / segment_length;
        if (u < -tol || u > T(1) + tol)
            return out;
        if (u < T(0))
            u = T(0);
        if (u > T(1))
            u = T(1);

        const auto a = point_span<T>(points, point_dim, i - 1);
        const auto b = point_span<T>(points, point_dim, i);
        out.point.assign(static_cast<std::size_t>(point_dim), T(0));
        for (int d = 0; d < point_dim; ++d)
            out.point[static_cast<std::size_t>(d)] =
                (T(1) - u) * a[static_cast<std::size_t>(d)]
              + u * b[static_cast<std::size_t>(d)];
        out.valid = true;
        out.segment = i - 1;
        out.segment_parameter = u;
        return out;
    }

    return out;
}

template <std::floating_point T>
inline ClosestArclength<T> closest_arclength_on_polyline(
    std::span<const T> points,
    int point_dim,
    std::span<const T> x,
    T tol = geom::default_tolerance<T>())
{
    ClosestArclength<T> out;
    if (!flat_point_array_valid<T>(points, point_dim)
        || x.size() != static_cast<std::size_t>(point_dim))
    {
        return out;
    }

    std::vector<T> cumulative;
    if (!cumulative_arclength<T>(points, point_dim, cumulative, tol))
        return out;

    const std::size_t n = cumulative.size();
    T best_distance2 = std::numeric_limits<T>::infinity();
    for (std::size_t i = 1; i < n; ++i)
    {
        const auto a = point_span<T>(points, point_dim, i - 1);
        const auto b = point_span<T>(points, point_dim, i);
        const auto ab = geom::subtract<T>(b, a);
        const T ab2 = geom::squared_norm<T>(vector_span(ab));
        if (ab2 <= tol * tol)
            continue;

        T u = T(0);
        for (int d = 0; d < point_dim; ++d)
            u += (x[static_cast<std::size_t>(d)] - a[static_cast<std::size_t>(d)])
               * ab[static_cast<std::size_t>(d)];
        u /= ab2;
        const T u_closest = std::clamp(u, T(0), T(1));

        T distance2 = T(0);
        for (int d = 0; d < point_dim; ++d)
        {
            const T closest =
                a[static_cast<std::size_t>(d)]
              + u_closest * ab[static_cast<std::size_t>(d)];
            const T r = x[static_cast<std::size_t>(d)] - closest;
            distance2 += r * r;
        }

        if (distance2 < best_distance2)
        {
            best_distance2 = distance2;
            out.valid = true;
            out.distance = std::sqrt(distance2);
            out.segment = i - 1;
            out.segment_parameter = u_closest;
            out.arclength =
                cumulative[i - 1]
              + u_closest * (cumulative[i] - cumulative[i - 1]);
        }
    }

    return out;
}

template <std::floating_point T>
struct LayoutReport
{
    bool accepted = false;
    FailureReason failure_reason = FailureReason::none;
    std::vector<T> arclengths;
    T minimum_spacing = std::numeric_limits<T>::infinity();
    T minimum_spacing_ratio = std::numeric_limits<T>::infinity();
};

template <std::floating_point T>
inline LayoutReport<T> evaluate_ordering_and_spacing(
    std::span<const T> provisional_curve_points,
    std::span<const T> edge_points,
    int point_dim,
    const Options<T>& options = {})
{
    LayoutReport<T> report;
    if (!flat_point_array_valid<T>(provisional_curve_points, point_dim)
        || !flat_point_array_valid<T>(edge_points, point_dim))
    {
        report.failure_reason = FailureReason::invalid_input;
        return report;
    }

    const std::size_t n =
        edge_points.size() / static_cast<std::size_t>(point_dim);
    if (n < 2)
    {
        report.failure_reason = FailureReason::invalid_input;
        return report;
    }

    std::vector<T> cumulative;
    if (!cumulative_arclength<T>(
            provisional_curve_points, point_dim, cumulative, options.tolerance))
    {
        report.failure_reason = FailureReason::invalid_curve;
        return report;
    }
    const T length = cumulative.back();
    const T min_spacing =
        options.min_absolute_spacing + options.min_relative_spacing * length;
    const T ordering_tol = options.tolerance * std::max(T(1), length);

    report.arclengths.assign(n, T(0));
    for (std::size_t i = 0; i < n; ++i)
    {
        const auto closest = closest_arclength_on_polyline<T>(
            provisional_curve_points,
            point_dim,
            point_span<T>(edge_points, point_dim, i),
            options.tolerance);
        if (!closest.valid)
        {
            report.failure_reason = FailureReason::invalid_curve;
            return report;
        }
        report.arclengths[i] = closest.arclength;
    }

    for (std::size_t i = 1; i < n; ++i)
    {
        const T gap = report.arclengths[i] - report.arclengths[i - 1];
        if (gap < -ordering_tol)
        {
            report.failure_reason = FailureReason::ordering_destroyed;
            return report;
        }

        report.minimum_spacing = std::min(report.minimum_spacing, gap);
        if (gap <= min_spacing)
        {
            report.minimum_spacing_ratio = length > options.tolerance
                ? gap / length
                : std::numeric_limits<T>::infinity();
            report.failure_reason = FailureReason::spacing_too_small;
            return report;
        }
    }

    report.accepted = true;
    report.failure_reason = FailureReason::none;
    report.minimum_spacing_ratio = length > options.tolerance
        ? report.minimum_spacing / length
        : std::numeric_limits<T>::infinity();
    return report;
}

template <std::floating_point T, class Phi, class Grad>
inline CorrectionReport<T> correct_seed_to_zero(
    cell::type parent_cell_type,
    geom::ParentEntity admissible_parent_entity,
    std::span<const T> seed,
    Phi&& phi,
    Grad&& grad,
    const Options<T>& options = {})
{
    CorrectionReport<T> report;
    const int tdim = cell::get_tdim(parent_cell_type);
    if (!admissible_parent_entity.valid()
        || admissible_parent_entity.dim > tdim
        || static_cast<int>(seed.size()) != tdim)
    {
        report.failure_reason = FailureReason::invalid_parent_entity;
        return report;
    }

    if (!geom::point_in_parent_entity<T>(
            parent_cell_type, admissible_parent_entity, seed, options.tolerance))
    {
        report.failure_reason = FailureReason::correction_left_parent_entity;
        report.left_parent_entity = true;
        return report;
    }

    report.point.assign(seed.begin(), seed.end());
    T f = static_cast<T>(std::invoke(
        phi,
        std::span<const T>(report.point.data(), report.point.size())));
    report.residual = std::fabs(f);
    if (!std::isfinite(report.residual))
    {
        report.failure_reason = FailureReason::correction_failed;
        return report;
    }
    if (report.residual <= options.phi_tolerance)
    {
        report.accepted = true;
        report.failure_reason = FailureReason::none;
        return report;
    }

    std::vector<T> g(static_cast<std::size_t>(tdim), T(0));
    bool full_step_left_parent = false;
    for (int iter = 0; iter < options.max_correction_iterations; ++iter)
    {
        std::fill(g.begin(), g.end(), T(0));
        std::invoke(
            grad,
            std::span<const T>(report.point.data(), report.point.size()),
            std::span<T>(g.data(), g.size()));

        const auto direction = geom::admissible_direction_in_parent_frame<T>(
            parent_cell_type,
            admissible_parent_entity,
            std::span<const T>(g.data(), g.size()),
            options.tolerance);
        if (direction.degenerate() || direction.norm <= options.min_direction_norm)
        {
            report.failure_reason = FailureReason::correction_direction_not_admissible;
            return report;
        }

        const T denom = geom::dot<T>(
            std::span<const T>(g.data(), g.size()),
            vector_span(direction.value));
        if (std::fabs(denom) <= options.min_direction_norm * direction.norm)
        {
            report.failure_reason = FailureReason::correction_direction_not_admissible;
            return report;
        }

        std::vector<T> step(static_cast<std::size_t>(tdim), T(0));
        const T scale = -f / denom;
        for (int d = 0; d < tdim; ++d)
            step[static_cast<std::size_t>(d)] =
                scale * direction.value[static_cast<std::size_t>(d)];
        if (geom::norm<T>(vector_span(step)) <= options.tolerance)
            break;

        std::vector<T> full_step_point(report.point);
        for (int d = 0; d < tdim; ++d)
            full_step_point[static_cast<std::size_t>(d)] +=
                step[static_cast<std::size_t>(d)];
        if (!geom::point_in_parent_entity<T>(
                parent_cell_type,
                admissible_parent_entity,
                vector_span(full_step_point),
                options.tolerance))
        {
            full_step_left_parent = true;
            report.left_parent_entity = true;
        }

        bool accepted_step = false;
        std::vector<T> candidate(static_cast<std::size_t>(tdim), T(0));
        T candidate_f = f;
        T alpha = T(1);
        for (int ls = 0; ls < options.max_line_search_iterations; ++ls)
        {
            for (int d = 0; d < tdim; ++d)
                candidate[static_cast<std::size_t>(d)] =
                    report.point[static_cast<std::size_t>(d)]
                  + alpha * step[static_cast<std::size_t>(d)];

            if (!geom::point_in_parent_entity<T>(
                    parent_cell_type,
                    admissible_parent_entity,
                    vector_span(candidate),
                    options.tolerance))
            {
                alpha *= T(0.5);
                continue;
            }

            candidate_f = static_cast<T>(std::invoke(
                phi,
                std::span<const T>(candidate.data(), candidate.size())));
            if (!std::isfinite(candidate_f))
            {
                alpha *= T(0.5);
                continue;
            }
            if (std::fabs(candidate_f) <= options.phi_tolerance
                || std::fabs(candidate_f) < report.residual)
            {
                accepted_step = true;
                break;
            }
            alpha *= T(0.5);
        }

        if (!accepted_step)
        {
            report.failure_reason = full_step_left_parent
                ? FailureReason::correction_left_parent_entity
                : FailureReason::correction_failed;
            return report;
        }

        report.point = std::move(candidate);
        f = candidate_f;
        report.residual = std::fabs(f);
        report.iterations = iter + 1;
        if (report.residual <= options.phi_tolerance)
        {
            report.accepted = true;
            report.failure_reason = FailureReason::none;
            return report;
        }
    }

    report.failure_reason = full_step_left_parent
        ? FailureReason::correction_left_parent_entity
        : FailureReason::phi_residual_too_large;
    return report;
}

template <std::floating_point T>
inline void set_failure(ShiftReport<T>& report,
                        FailureReason reason,
                        bool request_refinement = true)
{
    report.accepted = false;
    report.request_refinement = request_refinement;
    report.failure_reason = reason == FailureReason::none
        ? FailureReason::refinement_requested
        : reason;
}

template <std::floating_point T, class Phi>
inline bool validate_candidate(
    ShiftReport<T>& report,
    cell::type parent_cell_type,
    geom::ParentEntity admissible_parent_entity,
    std::span<const T> provisional_curve_points,
    int point_dim,
    int node_index,
    T desired_arclength,
    Phi&& phi,
    const Options<T>& options)
{
    if (!geom::point_in_parent_entity<T>(
            parent_cell_type,
            admissible_parent_entity,
            vector_span(report.corrected_node),
            options.tolerance))
    {
        set_failure<T>(report, FailureReason::correction_left_parent_entity);
        return false;
    }

    const T residual = std::fabs(static_cast<T>(std::invoke(
        phi,
        std::span<const T>(
            report.corrected_node.data(),
            report.corrected_node.size()))));
    report.metrics.phi_residual = residual;
    if (!std::isfinite(residual) || residual > options.phi_tolerance)
    {
        set_failure<T>(report, FailureReason::phi_residual_too_large);
        return false;
    }

    const auto layout = evaluate_ordering_and_spacing<T>(
        provisional_curve_points,
        std::span<const T>(
            report.corrected_edge_points.data(),
            report.corrected_edge_points.size()),
        point_dim,
        options);
    report.metrics.minimum_spacing = layout.minimum_spacing;
    report.metrics.minimum_spacing_ratio = layout.minimum_spacing_ratio;
    if (!layout.accepted)
    {
        set_failure<T>(report, layout.failure_reason);
        return false;
    }

    const auto final_closest = closest_arclength_on_polyline<T>(
        provisional_curve_points,
        point_dim,
        vector_span(report.corrected_node),
        options.tolerance);
    if (!final_closest.valid)
    {
        set_failure<T>(report, FailureReason::invalid_curve);
        return false;
    }

    report.metrics.final_arclength = final_closest.arclength;
    report.metrics.final_arclength_error =
        std::fabs(final_closest.arclength - desired_arclength);
    const T final_tol = arclength_tolerance<T>(
        report.metrics.curve_length,
        options.max_absolute_final_arclength_error,
        options.max_relative_final_arclength_error);
    if (report.metrics.final_arclength_error > final_tol)
    {
        set_failure<T>(report, FailureReason::arclength_error_too_large);
        return false;
    }

    if (node_index >= 0
        && static_cast<std::size_t>(node_index) < layout.arclengths.size())
    {
        report.metrics.final_arclength =
            layout.arclengths[static_cast<std::size_t>(node_index)];
    }

    report.accepted = true;
    report.request_refinement = false;
    report.failure_reason = FailureReason::none;
    return true;
}

template <std::floating_point T, class Phi, class Grad>
inline ShiftReport<T> correct_projected_edge_node(
    cell::type parent_cell_type,
    geom::ParentEntity admissible_parent_entity,
    std::span<const T> provisional_curve_points,
    std::span<const T> projected_edge_points,
    int point_dim,
    std::span<const T> desired_arclength_fractions,
    int node_index,
    Phi&& phi,
    Grad&& grad,
    const Options<T>& options = {})
{
    ShiftReport<T> report;
    if (!flat_point_array_valid<T>(provisional_curve_points, point_dim)
        || !flat_point_array_valid<T>(projected_edge_points, point_dim)
        || point_dim != cell::get_tdim(parent_cell_type))
    {
        set_failure<T>(report, FailureReason::invalid_input, false);
        return report;
    }

    const std::size_t n_nodes =
        projected_edge_points.size() / static_cast<std::size_t>(point_dim);
    if (desired_arclength_fractions.size() != n_nodes)
    {
        set_failure<T>(report, FailureReason::invalid_input, false);
        return report;
    }
    if (node_index < 0 || static_cast<std::size_t>(node_index) >= n_nodes)
    {
        set_failure<T>(report, FailureReason::node_index_out_of_range, false);
        return report;
    }

    const int tdim = cell::get_tdim(parent_cell_type);
    if (!admissible_parent_entity.valid()
        || admissible_parent_entity.dim > tdim)
    {
        set_failure<T>(report, FailureReason::invalid_parent_entity, false);
        return report;
    }

    std::vector<T> cumulative;
    if (!cumulative_arclength<T>(
            provisional_curve_points, point_dim, cumulative, options.tolerance))
    {
        set_failure<T>(report, FailureReason::invalid_curve);
        return report;
    }
    report.metrics.curve_length = cumulative.back();

    const auto current = point_span<T>(
        projected_edge_points,
        point_dim,
        static_cast<std::size_t>(node_index));
    if (!geom::point_in_parent_entity<T>(
            parent_cell_type, admissible_parent_entity, current, options.tolerance))
    {
        set_failure<T>(report, FailureReason::current_point_outside_parent_entity);
        return report;
    }

    const T fraction = desired_arclength_fractions[static_cast<std::size_t>(node_index)];
    if (fraction < -options.tolerance || fraction > T(1) + options.tolerance)
    {
        set_failure<T>(report, FailureReason::target_arclength_out_of_range);
        return report;
    }

    report.metrics.desired_arclength = fraction * report.metrics.curve_length;
    const auto actual = closest_arclength_on_polyline<T>(
        provisional_curve_points,
        point_dim,
        current,
        options.tolerance);
    if (!actual.valid)
    {
        set_failure<T>(report, FailureReason::invalid_curve);
        return report;
    }
    report.metrics.actual_arclength = actual.arclength;
    report.metrics.initial_arclength_error =
        std::fabs(actual.arclength - report.metrics.desired_arclength);

    report.corrected_edge_points.assign(
        projected_edge_points.begin(), projected_edge_points.end());
    report.corrected_node.assign(current.begin(), current.end());

    const T drift_tol = arclength_tolerance<T>(
        report.metrics.curve_length,
        options.max_absolute_arclength_drift,
        options.max_relative_arclength_drift);
    if (report.metrics.initial_arclength_error <= drift_tol)
    {
        validate_candidate<T>(
            report,
            parent_cell_type,
            admissible_parent_entity,
            provisional_curve_points,
            point_dim,
            node_index,
            report.metrics.desired_arclength,
            phi,
            options);
        return report;
    }

    const auto target = point_at_arclength<T>(
        provisional_curve_points,
        point_dim,
        report.metrics.desired_arclength,
        options.tolerance);
    if (!target.valid)
    {
        set_failure<T>(report, FailureReason::target_arclength_out_of_range);
        return report;
    }
    report.shifted_seed = target.point;

    auto correction = correct_seed_to_zero<T>(
        parent_cell_type,
        admissible_parent_entity,
        vector_span(report.shifted_seed),
        phi,
        grad,
        options);
    report.correction_iterations = correction.iterations;
    report.metrics.phi_residual = correction.residual;
    if (!correction.accepted)
    {
        set_failure<T>(report, correction.failure_reason);
        return report;
    }

    report.corrected = true;
    report.corrected_node = std::move(correction.point);
    auto corrected_span = std::span<T>(
        report.corrected_edge_points.data(),
        report.corrected_edge_points.size());
    const auto node = mutable_point_span<T>(
        corrected_span,
        point_dim,
        static_cast<std::size_t>(node_index));
    for (int d = 0; d < point_dim; ++d)
        node[static_cast<std::size_t>(d)] =
            report.corrected_node[static_cast<std::size_t>(d)];

    validate_candidate<T>(
        report,
        parent_cell_type,
        admissible_parent_entity,
        provisional_curve_points,
        point_dim,
        node_index,
        report.metrics.desired_arclength,
        phi,
        options);
    return report;
}

} // namespace cutcells::tangent_shift
