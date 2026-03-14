// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include "cell_types.h"
#include <span>
#include <vector>

namespace cutcells
{

/// Topology-only description of one iso-refinement step.
///
/// Reference coordinates live in the parent reference cell (tdim dimensions).
/// All parent entity indices are local to the parent background cell.
struct RefinementTemplate
{
    int n_vertices;           ///< total vertices in the refined mesh
    int n_cells;              ///< total child cells in the refined mesh
    int tdim;                 ///< topological dimension of background / refined cells
    int vertices_per_cell;    ///< vertices per child cell (uniform cell type)

    cell::type bg_cell_type;    ///< background (parent) cell type
    cell::type child_cell_type; ///< child cell type

    /// Reference vertex coordinates, size = n_vertices * tdim.
    std::vector<double> ref_vertex_coords;

    /// Parent entity dimension per vertex (0 = vertex, 1 = edge, 2 = face/cell interior in 2D).
    std::vector<int> vertex_parent_dim;

    /// Local parent entity index per vertex.
    /// For dim=0: local corner id; for dim=1: local edge id; for dim=2: 0 (cell interior).
    std::vector<int> vertex_parent_id;

    /// Flattened child cell-vertex connectivity,
    /// size = n_cells * vertices_per_cell.
    std::vector<int> cell_connectivity;
};

/// Trivial (no-refinement) P1 template for a given cell type.
const RefinementTemplate& p1_template(cell::type ct);

/// Triangle Pk-iso-P1 templates (Basix point ordering, equispaced variant).
const RefinementTemplate& p2_iso_p1_triangle_template();
const RefinementTemplate& p3_iso_p1_triangle_template();
const RefinementTemplate& p4_iso_p1_triangle_template();

/// Select triangle Pk-iso-P1 template for k in {2, 3, 4}.
const RefinementTemplate& triangle_iso_p1_template(int order);

/// Tetrahedron P2-iso-P1 template (Basix point ordering).
const RefinementTemplate& p2_iso_p1_tetrahedron_template();
const RefinementTemplate& p3_iso_p1_tetrahedron_template();
const RefinementTemplate& p4_iso_p1_tetrahedron_template();
const RefinementTemplate& tetrahedron_iso_p1_template(int order);

/// Quadrilateral Pk-iso-P1 templates (Basix point ordering).
const RefinementTemplate& p2_iso_p1_quadrilateral_template();
const RefinementTemplate& p3_iso_p1_quadrilateral_template();
const RefinementTemplate& p4_iso_p1_quadrilateral_template();
const RefinementTemplate& quadrilateral_iso_p1_template(int order);

/// Hexahedron Pk-iso-P1 templates (Basix point ordering).
const RefinementTemplate& p2_iso_p1_hexahedron_template();
const RefinementTemplate& p3_iso_p1_hexahedron_template();
const RefinementTemplate& p4_iso_p1_hexahedron_template();
const RefinementTemplate& hexahedron_iso_p1_template(int order);

/// Interval Pk-iso-P1 templates.
const RefinementTemplate& p2_iso_p1_interval_template();
const RefinementTemplate& p3_iso_p1_interval_template();
const RefinementTemplate& p4_iso_p1_interval_template();
const RefinementTemplate& interval_iso_p1_template(int order);

/// Prism Pk-iso-P1 templates (Basix point ordering).
const RefinementTemplate& p2_iso_p1_prism_template();
const RefinementTemplate& p3_iso_p1_prism_template();
const RefinementTemplate& p4_iso_p1_prism_template();
const RefinementTemplate& prism_iso_p1_template(int order);

/// Pyramid Pk-iso-P1 templates (Basix point ordering).
const RefinementTemplate& p2_iso_p1_pyramid_template();
const RefinementTemplate& p3_iso_p1_pyramid_template();
const RefinementTemplate& p4_iso_p1_pyramid_template();
const RefinementTemplate& pyramid_iso_p1_template(int order);

/// Generic selector:
///   point      -> p1 only
///   interval, triangle, quadrilateral, tetrahedron, hexahedron, prism, pyramid
///               -> orders 2,3,4 (where implemented)
const RefinementTemplate& iso_p1_template(cell::type ct, int order);

/// Reference coordinates for iso-P1 interpolation points in Basix ordering.
/// Size is n_vertices * tdim for the corresponding template.
std::span<const double> iso_p1_ref_coords(cell::type ct, int order);

/// Reference coordinates for P1 interpolation points in Basix ordering.
std::span<const double> p1_ref_coords(cell::type ct);

} // namespace cutcells
