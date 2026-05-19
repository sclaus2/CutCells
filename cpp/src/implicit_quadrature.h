// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <cstdint>

#include "adapt_cell.h"
#include "edge_root.h"
#include "ho_cut_mesh.h"
#include "level_set_cell.h"
#include "quadrature.h"

namespace cutcells::implicit_quadrature
{

enum class IntegrationDomain : std::uint8_t
{
    negative_volume = 0,
    positive_volume = 1,
    interface_surface = 2
};

struct ImplicitQuadratureOptions
{
    int order = 6;
    int inspection_order = 10;
    int max_refinement_iterations = 12;
    int edge_root_max_depth = 24;

    double zero_tol = 1e-12;
    double sign_tol = 1e-12;
    double root_tol = 1e-13;

    double min_transversality = 0.05;
    double max_root_slope = 10.0;
    double max_geometry_q_error = 1e-10;
    double max_surface_q_rel_error = 0.25;
    double min_surface_area_ratio = 0.02;
    double max_surface_area_ratio = 50.0;
    int max_base_split_depth = 8;
    int max_diagnostic_refinement_depth = 4;
    int max_height_direction_candidates = 32;

    cell::edge_root::method root_method = cell::edge_root::method::itp;

    bool enable_chart_refinement = true;
    bool enable_q_error_estimator = true;
    bool enable_height_direction_recovery = true;

    double min_height_transversality = 0.05;
    double max_face_trace_alpha = 0.95;
    double min_height_area_ratio = 0.25;
    double max_height_area_ratio = 4.0;
};

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> make_implicit_quadrature(
    const AdaptCell<T>& base_adapt_cell,
    const LevelSetCell<T, I>& ls_cell,
    int level_set_id,
    IntegrationDomain domain,
    const ImplicitQuadratureOptions& opts);

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> quadrature_rules_implicit(
    const HOMeshPart<T, I>& part,
    int order,
    bool include_uncut_cells,
    const ImplicitQuadratureOptions& opts = {});

} // namespace cutcells::implicit_quadrature
