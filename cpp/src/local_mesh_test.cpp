// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "local_mesh.h"
#include "edge_classification.h"
#include "cell_flags.h"

#include <cmath>
#include <iostream>
#include <span>
#include <vector>

#include "iso_refine.h"

namespace
{

bool approx_zero(double x, double tol = 1e-12)
{
  return std::abs(x) < tol;
}

int run_local_mesh_triangle_test()
{
  using namespace cutcells;

  // One P2 background triangle in Basix ordering.
  const std::vector<double> coords = {
      0.0, 0.0, // 0
      1.0, 0.0, // 1
      0.0, 1.0, // 2
      0.5, 0.5, // 3 edge (1,2)
      0.0, 0.5, // 4 edge (0,2)
      0.5, 0.0  // 5 edge (0,1)
  };
  LocalMesh<double> mesh;
  init_local_mesh_from_template(
      mesh, p2_iso_p1_triangle_template(),
      std::span<const double>(coords),
      cell::type::triangle,
      0, 1);

  if (mesh.n_vertices() != 6 || mesh.n_cells() != 4 || mesh.n_edges() != 9)
  {
    std::cerr << "Unexpected local-mesh sizes: vertices=" << mesh.n_vertices()
              << " cells=" << mesh.n_cells() << " edges=" << mesh.n_edges()
              << "\n";
    return 1;
  }

  // phi(x,y) = x + y - 0.3 intersects the reference triangle without
  // passing through existing template vertices.
  LevelSetFunction<double> phi0;
  phi0.gdim = 2;
  phi0.value_fn = [](const double* x, int) { return x[0] + x[1] - 0.3; };
  std::vector<LevelSetFunction<double>> phi = {phi0};

  evaluate_levelsets_on_vertices(mesh, phi, 1, 1e-14);
  classify_local_edges(mesh, 0);

  int one_root_edges = 0;
  for (uint8_t st : mesh.edge_state)
  {
    if (st == static_cast<uint8_t>(EdgeState::one_root))
      ++one_root_edges;
  }

  if (one_root_edges == 0)
  {
    std::cerr << "Expected at least one edge with one_root classification\n";
    std::cerr << "vertex diagnostics (x, y, phi):\n";
    for (int i = 0; i < mesh.n_vertices(); ++i)
    {
      std::cerr << "  v" << i << "=(" << mesh.vertex_x[2 * i] << ", "
                << mesh.vertex_x[2 * i + 1] << ") phi=" << mesh.vertex_phi[i]
                << " zero_mask=" << mesh.vertex_zero_mask[i]
                << " inside_mask=" << mesh.vertex_inside_mask[i] << "\n";
    }
    std::cerr << "edge diagnostics (v0, v1, state):\n";
    for (int e = 0; e < mesh.n_edges(); ++e)
    {
      std::cerr << "  e" << e << "=(" << mesh.edge_vertices[2 * e] << ", "
                << mesh.edge_vertices[2 * e + 1]
                << ") state=" << static_cast<int>(mesh.edge_state[e]) << "\n";
    }
    return 1;
  }

  const int n_vertices_before_roots = mesh.n_vertices();
  compute_all_roots_linear(mesh, 0);

  int rooted_edges = 0;
  for (int root_vertex : mesh.edge_root_vertex)
  {
    if (root_vertex >= 0)
      ++rooted_edges;
  }

  if (rooted_edges != one_root_edges)
  {
    std::cerr << "Root count mismatch: rooted_edges=" << rooted_edges
              << " one_root_edges=" << one_root_edges << "\n";
    return 1;
  }

  const int n_new_vertices = mesh.n_vertices() - n_vertices_before_roots;
  if (n_new_vertices != one_root_edges)
  {
    std::cerr << "Unexpected number of new root vertices: " << n_new_vertices
              << ", expected " << one_root_edges << "\n";
    return 1;
  }

  // Root vertices should be zero in the level-set slot by construction.
  for (int e = 0; e < mesh.n_edges(); ++e)
  {
    const int rv = mesh.edge_root_vertex[e];
    if (rv < 0)
      continue;
    if (!approx_zero(mesh.vertex_phi[rv * mesh.n_level_sets]))
    {
      std::cerr << "Root vertex phi is not approximately zero on edge " << e
                << "\n";
      return 1;
    }
    if (mesh.vertex_root_edge_id[static_cast<std::size_t>(rv)] != e)
    {
      std::cerr << "Root vertex does not store its local edge id\n";
      return 1;
    }
  }

  // Decompose intersected cells into both inside and outside fragments using the
  // cached linear roots.
  decompose_local_mesh_linear(mesh, 0, true);

  if (mesh.n_cells() <= 0)
  {
    std::cerr << "Expected non-empty decomposed local mesh\n";
    return 1;
  }
  int n_inside = 0;
  int n_outside = 0;
  for (const auto dom_u8 : mesh.cell_domain)
  {
    const auto dom = static_cast<cell::domain>(dom_u8);
    if (dom == cell::domain::inside)
      ++n_inside;
    else if (dom == cell::domain::outside)
      ++n_outside;
  }
  if (n_inside == 0 || n_outside == 0)
  {
    std::cerr << "Decomposition should contain both inside and outside fragments\n";
    return 1;
  }

  return 0;
}

int run_edge_classification_refine_test()
{
  using namespace cutcells;

  const std::vector<double> coords = {
      0.0, 0.0, // 0
      1.0, 0.0, // 1
      0.0, 1.0, // 2
      0.5, 0.5, // 3 edge (1,2)
      0.0, 0.5, // 4 edge (0,2)
      0.5, 0.0  // 5 edge (0,1)
  };

  LocalMesh<double> mesh;
  init_local_mesh_from_cell(
      mesh,
      std::span<const double>(coords),
      cell::type::triangle,
      0,
      1);

  if (mesh.n_cells() != 1 || mesh.n_edges() != 3)
  {
    std::cerr << "Unexpected initial mesh sizes in edge classification test\n";
    return 1;
  }
  if (mesh.parent_vertex_to_local_vertex.size() != 3
      || mesh.parent_edge_to_local_edge.size() != 3)
  {
    std::cerr << "Unexpected parent-entity map sizes in local mesh\n";
    return 1;
  }
  if (mesh.parent_vertex_to_local_vertex[0] < 0
      || mesh.parent_vertex_to_local_vertex[1] < 0
      || mesh.parent_vertex_to_local_vertex[2] < 0)
  {
    std::cerr << "Parent vertex to local vertex map not populated\n";
    return 1;
  }
  if (mesh.parent_edge_to_local_edge[0] < 0
      || mesh.parent_edge_to_local_edge[1] < 0
      || mesh.parent_edge_to_local_edge[2] < 0)
  {
    std::cerr << "Parent edge to local edge map not populated\n";
    return 1;
  }

  LevelSetFunction<double> phi;
  phi.gdim = 2;
  // Endpoints on edge (1,2) are positive, edge-interior node is negative:
  // this must be classified as multiple_roots and trigger refinement.
  const std::vector<double> nodal_phi = {1.0, 1.0, 1.0, -1.0, 1.0, 1.0};
  phi.nodal_values = std::span<const double>(nodal_phi.data(), nodal_phi.size());

  const bool should_refine = classify_edges_and_mark_refine(mesh, phi, 0, 1e-14);
  if (!should_refine)
  {
    std::cerr << "Expected refinement trigger from multiple_roots edge classification\n";
    return 1;
  }
  if (mesh.edge_state.empty()
      || mesh.edge_state[0] != static_cast<uint8_t>(EdgeState::multiple_roots))
  {
    std::cerr << "Expected edge 0 to be classified as multiple_roots\n";
    return 1;
  }

  const RefinementTemplate& tpl = iso_p1_template(cell::type::triangle, 2);
  const bool did_refine = classify_edges_and_refine_cell(mesh, phi, tpl, 0, 1e-14);

  if (!did_refine)
  {
    std::cerr << "Expected refinement trigger in classify_edges_and_refine_cell\n";
    return 1;
  }
  if (mesh.n_cells() != tpl.n_cells)
  {
    std::cerr << "Refined mesh has unexpected number of cells: " << mesh.n_cells()
              << " expected " << tpl.n_cells << "\n";
    return 1;
  }

  // Iterative driver: should converge (no ambiguous edges) for smooth affine phi.
  LevelSetFunction<double> phi_affine;
  phi_affine.gdim = 2;
  phi_affine.value_fn = [](const double* x, int) { return x[0] + x[1] - 0.3; };
  LocalMesh<double> mesh_iter;
  init_local_mesh_from_cell(
      mesh_iter,
      std::span<const double>(coords),
      cell::type::triangle,
      0,
      1);
  const auto iter_res = refine_until_edges_single_root(
      mesh_iter, phi_affine, tpl, 0, 1e-14, 4, false);
  if (!iter_res.converged)
  {
    std::cerr << "Expected iterative edge refinement to converge\n";
    return 1;
  }

  return 0;
}

int run_iso_refine_scheme_test()
{
  using namespace cutcells;

  const RefinementTemplate& t2 = triangle_iso_p1_template(2);
  const RefinementTemplate& t3 = triangle_iso_p1_template(3);
  const RefinementTemplate& t4 = triangle_iso_p1_template(4);
  const auto t2_ref = iso_p1_ref_coords(cell::type::triangle, 2);
  const auto t3_ref = iso_p1_ref_coords(cell::type::triangle, 3);
  const auto t4_ref = iso_p1_ref_coords(cell::type::triangle, 4);

  if (t2.n_vertices != 6 || t2.n_cells != 4 || t2.vertices_per_cell != 3)
  {
    std::cerr << "Unexpected P2 iso template sizes\n";
    return 1;
  }
  if (t3.n_vertices != 10 || t3.n_cells != 9 || t3.vertices_per_cell != 3)
  {
    std::cerr << "Unexpected P3 iso template sizes\n";
    return 1;
  }
  if (t4.n_vertices != 15 || t4.n_cells != 16 || t4.vertices_per_cell != 3)
  {
    std::cerr << "Unexpected P4 iso template sizes\n";
    return 1;
  }

  // Basix ordering spot-checks (equispaced triangle):
  // p2: point 3 must be (0.5, 0.5), p3: point 9 must be (1/3, 1/3),
  // p4: points 12..14 must be (1/4,1/4), (1/2,1/4), (1/4,1/2).
  if (!approx_zero(t2_ref[2 * 3] - 0.5)
      || !approx_zero(t2_ref[2 * 3 + 1] - 0.5))
  {
    std::cerr << "P2 Basix point ordering mismatch at index 3\n";
    return 1;
  }

  if (!approx_zero(t3_ref[2 * 9] - (1.0 / 3.0))
      || !approx_zero(t3_ref[2 * 9 + 1] - (1.0 / 3.0)))
  {
    std::cerr << "P3 Basix point ordering mismatch at index 9\n";
    return 1;
  }

  if (!approx_zero(t4_ref[2 * 12] - 0.25)
      || !approx_zero(t4_ref[2 * 12 + 1] - 0.25)
      || !approx_zero(t4_ref[2 * 13] - 0.5)
      || !approx_zero(t4_ref[2 * 13 + 1] - 0.25)
      || !approx_zero(t4_ref[2 * 14] - 0.25)
      || !approx_zero(t4_ref[2 * 14 + 1] - 0.5))
  {
    std::cerr << "P4 Basix interior point ordering mismatch\n";
    return 1;
  }

  const RefinementTemplate& tet2 = iso_p1_template(cell::type::tetrahedron, 2);
  const RefinementTemplate& quad2 = iso_p1_template(cell::type::quadrilateral, 2);
  const RefinementTemplate& quad3 = iso_p1_template(cell::type::quadrilateral, 3);
  const RefinementTemplate& quad4 = iso_p1_template(cell::type::quadrilateral, 4);
  const RefinementTemplate& hex2 = iso_p1_template(cell::type::hexahedron, 2);
  const RefinementTemplate& hex3 = iso_p1_template(cell::type::hexahedron, 3);
  const RefinementTemplate& hex4 = iso_p1_template(cell::type::hexahedron, 4);
  const RefinementTemplate& prism2 = iso_p1_template(cell::type::prism, 2);
  const RefinementTemplate& prism3 = iso_p1_template(cell::type::prism, 3);
  const RefinementTemplate& prism4 = iso_p1_template(cell::type::prism, 4);
  const RefinementTemplate& pyr2 = iso_p1_template(cell::type::pyramid, 2);
  const RefinementTemplate& pyr3 = iso_p1_template(cell::type::pyramid, 3);
  const RefinementTemplate& pyr4 = iso_p1_template(cell::type::pyramid, 4);
  const RefinementTemplate& int2 = iso_p1_template(cell::type::interval, 2);
  const RefinementTemplate& int3 = iso_p1_template(cell::type::interval, 3);
  const RefinementTemplate& int4 = iso_p1_template(cell::type::interval, 4);

  if (tet2.n_vertices != 10 || tet2.n_cells != 8 || tet2.vertices_per_cell != 4)
  {
    std::cerr << "Unexpected tetrahedron P2 iso template sizes\n";
    return 1;
  }
  if (quad2.n_vertices != 9 || quad2.n_cells != 8 || quad2.vertices_per_cell != 3)
  {
    std::cerr << "Unexpected quadrilateral P2 iso template sizes\n";
    return 1;
  }
  if (quad3.n_vertices != 16 || quad3.n_cells != 18 || quad3.vertices_per_cell != 3)
  {
    std::cerr << "Unexpected quadrilateral P3 iso template sizes\n";
    return 1;
  }
  if (quad4.n_vertices != 25 || quad4.n_cells != 32 || quad4.vertices_per_cell != 3)
  {
    std::cerr << "Unexpected quadrilateral P4 iso template sizes\n";
    return 1;
  }
  if (hex2.n_vertices != 27 || hex2.n_cells != 40 || hex2.vertices_per_cell != 4)
  {
    std::cerr << "Unexpected hexahedron P2 iso template sizes\n";
    return 1;
  }
  if (hex3.n_vertices != 64 || hex3.n_cells != 135 || hex3.vertices_per_cell != 4)
  {
    std::cerr << "Unexpected hexahedron P3 iso template sizes\n";
    return 1;
  }
  if (hex4.n_vertices != 125 || hex4.n_cells != 320 || hex4.vertices_per_cell != 4)
  {
    std::cerr << "Unexpected hexahedron P4 iso template sizes\n";
    return 1;
  }
  if (prism2.n_vertices != 18 || prism2.n_cells != 24 || prism2.vertices_per_cell != 4)
  {
    std::cerr << "Unexpected prism P2 iso template sizes\n";
    return 1;
  }
  if (prism3.n_vertices != 40 || prism3.n_cells != 81 || prism3.vertices_per_cell != 4)
  {
    std::cerr << "Unexpected prism P3 iso template sizes\n";
    return 1;
  }
  if (prism4.n_vertices != 75 || prism4.n_cells != 192 || prism4.vertices_per_cell != 4)
  {
    std::cerr << "Unexpected prism P4 iso template sizes\n";
    return 1;
  }
  if (pyr2.n_vertices != 14 || pyr2.n_cells != 5 || pyr2.vertices_per_cell != 5)
  {
    std::cerr << "Unexpected pyramid P2 iso template sizes\n";
    return 1;
  }
  if (pyr3.n_vertices != 30 || pyr3.n_cells != 14 || pyr3.vertices_per_cell != 5)
  {
    std::cerr << "Unexpected pyramid P3 iso template sizes\n";
    return 1;
  }
  if (pyr4.n_vertices != 55 || pyr4.n_cells != 30 || pyr4.vertices_per_cell != 5)
  {
    std::cerr << "Unexpected pyramid P4 iso template sizes\n";
    return 1;
  }
  if (int2.n_vertices != 3 || int2.n_cells != 2 || int2.vertices_per_cell != 2)
  {
    std::cerr << "Unexpected interval P2 iso template sizes\n";
    return 1;
  }
  if (int3.n_vertices != 4 || int3.n_cells != 3 || int3.vertices_per_cell != 2)
  {
    std::cerr << "Unexpected interval P3 iso template sizes\n";
    return 1;
  }
  if (int4.n_vertices != 5 || int4.n_cells != 4 || int4.vertices_per_cell != 2)
  {
    std::cerr << "Unexpected interval P4 iso template sizes\n";
    return 1;
  }

  // Basix spot-checks on key entity points.
  const auto tet2_ref = iso_p1_ref_coords(cell::type::tetrahedron, 2);
  const auto quad2_ref = iso_p1_ref_coords(cell::type::quadrilateral, 2);
  const auto hex2_ref = iso_p1_ref_coords(cell::type::hexahedron, 2);
  const auto prism2_ref = iso_p1_ref_coords(cell::type::prism, 2);
  const auto pyr2_ref = iso_p1_ref_coords(cell::type::pyramid, 2);

  if (!approx_zero(tet2_ref[3 * 4 + 1] - 0.5)
      || !approx_zero(tet2_ref[3 * 4 + 2] - 0.5))
  {
    std::cerr << "Tetrahedron P2 edge midpoint ordering mismatch\n";
    return 1;
  }
  if (!approx_zero(quad2_ref[2 * 8] - 0.5)
      || !approx_zero(quad2_ref[2 * 8 + 1] - 0.5))
  {
    std::cerr << "Quadrilateral P2 center point ordering mismatch\n";
    return 1;
  }
  if (!approx_zero(hex2_ref[3 * 26] - 0.5)
      || !approx_zero(hex2_ref[3 * 26 + 1] - 0.5)
      || !approx_zero(hex2_ref[3 * 26 + 2] - 0.5))
  {
    std::cerr << "Hexahedron P2 interior point ordering mismatch\n";
    return 1;
  }
  if (!approx_zero(prism2_ref[3 * 17] - 0.5)
      || !approx_zero(prism2_ref[3 * 17 + 1] - 0.5)
      || !approx_zero(prism2_ref[3 * 17 + 2] - 0.5))
  {
    std::cerr << "Prism P2 interior point ordering mismatch\n";
    return 1;
  }
  if (!approx_zero(pyr2_ref[3 * 13] - 0.5)
      || !approx_zero(pyr2_ref[3 * 13 + 1] - 0.5)
      || !approx_zero(pyr2_ref[3 * 13 + 2]))
  {
    std::cerr << "Pyramid P2 base-face interior ordering mismatch\n";
    return 1;
  }

  return 0;
}

} // namespace

int main()
{
  if (run_iso_refine_scheme_test() != 0)
    return 1;
  if (run_local_mesh_triangle_test() != 0)
    return 1;
  return run_edge_classification_refine_test();
}
