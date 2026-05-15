// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <vector>

#include "cut_mesh.h"
#include "ho_cut_mesh.h"
#include "quadrature.h"

namespace cutcells::output
{

struct SelectedZeroEntityInfo
{
    std::int32_t cut_cell_id = -1;
    std::int32_t parent_cell_id = -1;
    std::int32_t local_zero_entity_id = -1;
    std::int32_t dimension = -1;
};

template <std::floating_point T, std::integral I = int>
std::vector<SelectedZeroEntityInfo> selected_zero_entity_infos(
    const HOMeshPart<T, I>& part);

template <std::floating_point T, std::integral I = int>
mesh::CutMesh<T> visualization_mesh(const HOMeshPart<T, I>& part,
                                    bool include_uncut_cells);

template <std::floating_point T, std::integral I = int>
quadrature::QuadratureRules<T> quadrature_rules(const HOMeshPart<T, I>& part,
                                                int order,
                                                bool include_uncut_cells);

} // namespace cutcells::output
