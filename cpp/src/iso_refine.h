// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <concepts>
#include <span>
#include <vector>

#include "adapt_cell.h"
#include "cell_types.h"

namespace cutcells
{

/// Topology-only description of one Pk-iso-P1 refinement template.
///
/// Reference coordinates live in the parent reference cell. Parent entity ids
/// are local to the parent cell for edges/faces and local corner ids for
/// vertices.
struct IsoRefineTemplate
{
    int n_vertices = 0;
    int n_cells = 0;
    int tdim = 0;
    int vertices_per_cell = 0;
    cell::type parent_cell_type = cell::type::point;
    cell::type child_cell_type = cell::type::point;
    std::vector<double> ref_vertex_coords;
    std::vector<int> vertex_parent_dim;
    std::vector<int> vertex_parent_id;
    std::vector<int> cell_connectivity;
};

/// Backward-compatible terminology from the older CutCells cutting path.
using RefinementTemplate = IsoRefineTemplate;

/// Trivial P1 template for a parent cell type.
const IsoRefineTemplate& p1_template(cell::type cell_type);

/// Select a Pk-iso-P1 template for k in {1, 2, 3, 4}.
///
/// V1 supports interval, triangle, tetrahedron, quadrilateral, and hexahedron.
/// Quadrilateral and hexahedron templates use classical tensor-product
/// subdivision into smaller quadrilaterals and hexahedra.
const IsoRefineTemplate& iso_p1_template(cell::type cell_type, int order);

/// Reference coordinates for the selected iso-P1 template.
std::span<const double> iso_p1_ref_coords(cell::type cell_type, int order);

/// Replace the initial leaf topology in an AdaptCell by an iso-P1 template.
///
/// This must be called before level-set sign stamping and certification.
template <std::floating_point T>
void apply_iso_refine(AdaptCell<T>& adapt_cell, const IsoRefineTemplate& tpl);

} // namespace cutcells
