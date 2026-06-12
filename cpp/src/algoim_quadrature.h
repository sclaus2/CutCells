// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <string>
#include <utility>
#include <vector>

#include "ho_cut_mesh.h"
#include "quadrature.h"

namespace cutcells::output
{

template <std::floating_point T, std::integral I = int>
quadrature::QuadratureRules<T> algoim_quadrature_rules(
    const HOMeshPart<T, I>& part,
    int order,
    bool include_uncut_cells);

template <std::floating_point T, std::integral I = int>
quadrature::QuadratureRules<T> algoim_general_quadrature_rules(
    const HOMeshPart<T, I>& part,
    int order,
    bool include_uncut_cells);

template <std::floating_point T, std::integral I = int>
std::vector<std::pair<std::string, quadrature::QuadratureRules<T>>>
algoim_paired_quadrature_rules(
    const std::vector<std::pair<std::string, HOMeshPart<T, I>>>& parts,
    int order,
    bool include_uncut_cells);

} // namespace cutcells::output
