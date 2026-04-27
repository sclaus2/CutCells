// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <cstdint>
#include <vector>

#include "cut_mesh.h"
#include "curving.h"
#include "ho_cut_mesh.h"
#include "quadrature.h"

namespace cutcells::output
{

template <std::floating_point T>
struct CurvedVTUGrid
{
    int gdim = 0;
    int tdim = 0;
    int geometry_order = 1;

    std::vector<T> points;
    std::vector<int> connectivity;
    std::vector<int> offsets;
    std::vector<int> vtk_types;

    std::vector<std::int32_t> parent_map;
    std::vector<std::int32_t> curved_valid;
    std::vector<std::int32_t> subdivision_depth;
    std::vector<std::int32_t> curving_status;
};

template <std::floating_point T, std::integral I = int>
mesh::CutMesh<T> visualization_mesh(const HOMeshPart<T, I>& part,
                                    bool include_uncut_cells,
                                    int geometry_order = -1,
                                    curving::NodeFamily node_family = curving::NodeFamily::gll);

template <std::floating_point T, std::integral I = int>
quadrature::QuadratureRules<T> quadrature_rules(const HOMeshPart<T, I>& part,
                                                int order,
                                                bool include_uncut_cells,
                                                int geometry_order = -1,
                                                curving::NodeFamily node_family = curving::NodeFamily::gll);

template <std::floating_point T, std::integral I = int>
CurvedVTUGrid<T> curved_lagrange_grid(const HOMeshPart<T, I>& part,
                                      bool include_uncut_cells,
                                      int geometry_order,
                                      curving::NodeFamily node_family = curving::NodeFamily::lagrange);

} // namespace cutcells::output
