// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "rules.h"

#include <stdexcept>

namespace cutcells::quadrays
{

template <std::floating_point T>
void append_rules(const ClippedBox<T>& cell, const BoxBernstein<T>& phi, Part part, int q,
                  const Options& opt, std::int32_t parent_cell,
                  quadrature::QuadratureRules<T>& rules, Stats& stats)
{
    thread_local CellPoints<T> points;
    points.points.clear();
    points.weights.clear();
    points.tags.clear();
    integrate(cell, phi, part, q, opt, points, stats);
    if (points.n_points() == 0)
        return;

    if (rules._tdim == 0)
        rules._tdim = 3;
    if (rules._tdim != 3)
        throw std::invalid_argument("quadrays: rules of 3D cells need _tdim = 3");
    if (rules._offset.empty())
        rules._offset.push_back(0);
    rules._points.reserve(rules._points.size() + points.points.size());
    rules._weights.reserve(rules._weights.size() + points.weights.size());
    for (int i = 0; i < points.n_points(); ++i)
    {
        const Vec3<T> u = {points.points[3 * i], points.points[3 * i + 1], points.points[3 * i + 2]};
        const Vec3<T> xi = reference_point(cell, u);
        rules._points.insert(rules._points.end(), xi.begin(), xi.end());
        rules._weights.push_back(points.weights[i]);
    }
    rules._offset.push_back(static_cast<std::int32_t>(rules._weights.size()));
    rules._parent_map.push_back(parent_cell);
}

template <std::floating_point T>
void append_cell_rules(cell::type cell_type, std::span<const T> vertex_coords, int degree,
                       std::span<const T> coeffs, const SelectionTerm& term, int level_set, int q,
                       const Options& opt, std::int32_t parent_cell,
                       quadrature::QuadratureRules<T>& rules, Stats& stats)
{
    thread_local ClippedBox<T> cell;
    thread_local BoxBernstein<T> phi;
    make_clipped_box(cell_type, vertex_coords, 3, cell);
    cell_bernstein_on_box(cell_type, degree, coeffs, phi);
    append_rules(cell, phi, part_of(term, level_set), q, opt, parent_cell, rules, stats);
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template void append_rules<float>(const ClippedBox<float>&, const BoxBernstein<float>&, Part, int,
                                  const Options&, std::int32_t, quadrature::QuadratureRules<float>&,
                                  Stats&);
template void append_rules<double>(const ClippedBox<double>&, const BoxBernstein<double>&, Part, int,
                                   const Options&, std::int32_t, quadrature::QuadratureRules<double>&,
                                   Stats&);
template void append_cell_rules<float>(cell::type, std::span<const float>, int, std::span<const float>,
                                       const SelectionTerm&, int, int, const Options&, std::int32_t,
                                       quadrature::QuadratureRules<float>&, Stats&);
template void append_cell_rules<double>(cell::type, std::span<const double>, int,
                                        std::span<const double>, const SelectionTerm&, int, int,
                                        const Options&, std::int32_t,
                                        quadrature::QuadratureRules<double>&, Stats&);

} // namespace cutcells::quadrays
