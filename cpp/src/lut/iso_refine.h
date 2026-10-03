// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#pragma once

#include <span>
#include <vector>

#include "../cell_types.h"

namespace cutcells
{

/// Topology-only description of one Pk-iso-P1 refinement template.
///
/// Reference coordinates live in the parent reference cell. Parent entity ids
/// are local to the parent cell for edges/faces and local corner ids for
/// vertices. Children are in CSR layout; a template whose children differ in
/// type (the pyramid's: pyramids and tetrahedra) has vertices_per_cell 0 and
/// child_cell_type point.
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
    std::vector<int> cell_offsets;      ///< child c: cell_connectivity[cell_offsets[c]] to [cell_offsets[c + 1] - 1]
    std::vector<cell::type> cell_types; ///< per child
};

/// Backward-compatible name.
using RefinementTemplate = IsoRefineTemplate;

/// Trivial P1 template for a parent cell type.
const IsoRefineTemplate& p1_template(cell::type cell_type);

/// Select a Pk-iso-P1 template for k in {1, 2, 3, 4}.
///
/// Intervals, triangles, tetrahedra, quadrilaterals, hexahedra and prisms;
/// pyramids for k in {1, 2}. Quadrilateral and hexahedron templates use
/// classical tensor-product subdivision into smaller quadrilaterals and
/// hexahedra, prisms k^2 triangles in k layers; a pyramid of order 2 splits
/// into six pyramids and four tetrahedra on its 14 nodes.
const IsoRefineTemplate& iso_p1_template(cell::type cell_type, int order);

/// Reference coordinates for the selected iso-P1 template.
std::span<const double> iso_p1_ref_coords(cell::type cell_type, int order);

} // namespace cutcells
