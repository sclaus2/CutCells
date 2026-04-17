// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>

#include "cut_mesh.h"
#include "ho_cut_mesh.h"
#include "quadrature.h"

namespace cutcells::output
{

template <std::floating_point T, std::integral I = int>
mesh::CutMesh<T> visualization_mesh(const HOMeshPart<T, I>& part,
                                    bool include_uncut_cells,
                                    bool triangulate);

template <std::floating_point T, std::integral I = int>
quadrature::QuadratureRules<T> quadrature_rules(const HOMeshPart<T, I>& part,
                                                int order,
                                                bool include_uncut_cells,
                                                bool triangulate);

} // namespace cutcells::output
