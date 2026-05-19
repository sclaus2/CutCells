// Copyright (c) 2024 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include "quadrature.h"
#include "quadrature_tables.h"
#include "cut_cell.h"
#include "triangulation.h"
#include "cell_flags.h"
#include "mapping.h"

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace cutcells::quadrature
{

// =============================================================================
// Degree-1 shape functions and their Jacobians on canonical cut-subcell types.
//
// For each subcell type the canonical reference vertices are (Basix ordering):
//
//   interval:      v0=(0), v1=(1)
//   triangle:      v0=(0,0), v1=(1,0), v2=(0,1)
//   tetrahedron:   v0=(0,0,0), v1=(1,0,0), v2=(0,1,0), v3=(0,0,1)
//
// Degree-1 shape functions on the canonical cell:
//   interval:      N0=1-s, N1=s
//   triangle:      N0=1-s-t, N1=s, N2=t
//   tetrahedron:   N0=1-s-t-u, N1=s, N2=t, N3=u
//
// The mapped point at canonical coordinate (s,[t,[u]]) in an actual subcell
// with real vertices V0..Vn is:
//   x = sum_i  N_i(s,t,...) * Vi
//
// The Jacobian of that map is:
//   J = [V1-V0, V2-V0, ...]   (column-major, one column per coordinate)
//
// Note: quad, hex, prism, pyramid subcells do NOT appear in the CutCells
// output (the library always triangulates cut regions into simplices).
// We only need interval, triangle, tetrahedron.
// =============================================================================

namespace
{

// ---------------------------------------------------------------------------
// Map a batch of n canonical quadrature points to actual subcell coordinates.
//
// can_pts : n * tdim  canonical reference points (flat row-major)
// verts   : nv * gdim actual subcell vertex coordinates (flat row-major,
//                     positions in whichever space we are mapping in)
// cell_type : subcell type (determines basis)
// out_pts : n * gdim  output points (flat row-major)
//
// For interval subcells embedded in higher-dimensional space the "gdim" output
// dimension follows the vertex dimension; tdim == 1 but vertices live in gdim.
// ---------------------------------------------------------------------------
template <typename T>
void map_canonical_to_subcell(const T* can_pts,  int n,
                               int tdim,
                               const T* verts,   int gdim,
                               cell::type cell_type,
                               T* out_pts)
{
    // N_i are standard Lagrange-1 basis functions on the canonical simplex.
    // J = [v1-v0 | v2-v0 | ...] (each column = vertex difference)
    //
    // x = v0 + J * can_pt
    //
    // This is identical to the affine push-forward used in mapping.cpp but
    // here the "parent" is the canonical subcell itself.

    const T* v0 = verts;

    for (int q = 0; q < n; ++q)
    {
        const T* X = can_pts + q * tdim;
        T* x = out_pts + q * gdim;

        // Start at v0
        for (int d = 0; d < gdim; ++d)
            x[d] = v0[d];

        // Add contributions from each basis function N_i = X[i-1], i >= 1
        for (int i = 1; i <= tdim; ++i)
        {
            const T* vi = verts + i * gdim;
            for (int d = 0; d < gdim; ++d)
                x[d] += X[i - 1] * (vi[d] - v0[d]);
        }
    }
}

// ---------------------------------------------------------------------------
// Compute |det J| for the simplex spanned by subcell vertices.
//
// verts: nv * gdim  (row-major)
// cell_type: interval(1D→1D/2D/3D), triangle(2D→2D/3D), tetrahedron(3D→3D)
//
// For embedded simplices (tdim < gdim) we compute the Gramian determinant
// sqrt(det(J^T J)).
// ---------------------------------------------------------------------------
template <typename T>
T simplex_physical_volume_factor(const T* verts, int tdim, int gdim)
{
    // Build J: j[col*gdim + row] = v_{col+1} - v_0,  col in [0, tdim)
    T J[9] = {};
    const T* v0 = verts;
    for (int col = 0; col < tdim; ++col)
    {
        const T* vi = verts + (col + 1) * gdim;
        for (int row = 0; row < gdim; ++row)
            J[col * gdim + row] = vi[row] - v0[row];
    }

    if (tdim == gdim)
    {
        // Square Jacobian: use direct det
        if (tdim == 1)
        {
            return std::abs(J[0]);
        }
        else if (tdim == 2)
        {
            return std::abs(J[0] * J[3] - J[2] * J[1]);
        }
        else // tdim == 3
        {
            const T det =
                J[0] * (J[4] * J[8] - J[7] * J[5])
              - J[3] * (J[1] * J[8] - J[7] * J[2])
              + J[6] * (J[1] * J[5] - J[4] * J[2]);
            return std::abs(det);
        }
    }
    else
    {
        // Non-square: Gramian sqrt(det(J^T J))
        // G[i][j] = sum_k J[i*gdim+k] * J[j*gdim+k]
        T G[9] = {};
        for (int i = 0; i < tdim; ++i)
            for (int j = 0; j < tdim; ++j)
            {
                T s = 0;
                for (int k = 0; k < gdim; ++k)
                    s += J[i * gdim + k] * J[j * gdim + k];
                G[i * tdim + j] = s;
            }

        if (tdim == 1)
            return std::sqrt(G[0]);
        else if (tdim == 2)
            return std::sqrt(G[0] * G[3] - G[1] * G[2]);
        else
        {
            const T det =
                G[0] * (G[4] * G[8] - G[7] * G[5])
              - G[3] * (G[1] * G[8] - G[7] * G[2])
              + G[6] * (G[1] * G[5] - G[4] * G[2]);
            return std::sqrt(det);
        }
    }
}

// ---------------------------------------------------------------------------
// Gather vertex coordinates for subcell id from either the reference or
// physical vertex array of cut_cell.
// ---------------------------------------------------------------------------
template <typename T>
void gather_subcell_verts(const cutcells::cell::CutCell<T>& cut_cell,
                          int subcell_id,
                          const std::vector<T>& coord_array,
                          int gdim,
                          std::vector<T>& out)
{
    const auto verts = cutcells::cell::cell_vertices(cut_cell, subcell_id);
    const int nv = static_cast<int>(verts.size());
    out.resize(static_cast<std::size_t>(nv) * gdim);
    for (int j = 0; j < nv; ++j)
        for (int d = 0; d < gdim; ++d)
            out[static_cast<std::size_t>(j) * gdim + d] =
                coord_array[static_cast<std::size_t>(verts[j]) * gdim + d];
}

} // anonymous namespace

// =============================================================================
// append_quadrature
// =============================================================================
template <std::floating_point T>
void append_quadrature(const cutcells::cell::CutCell<T>& cut_cell,
                       int order,
                       QuadratureRules<T>& rules)
{
    const int gdim = cut_cell._gdim;
    const int tdim = cut_cell._tdim;      // topological dim of the parent cell
    const int nc   = cutcells::cell::num_cells(cut_cell);

    if (nc == 0)
        return;

    // Initialise _tdim on first call
    if (rules._tdim == 0)
        rules._tdim = tdim;

    // Ensure offset is initialised
    if (rules._offset.empty())
        rules._offset.push_back(0);

    const int32_t pts_before = static_cast<int32_t>(rules._weights.size());

    // Scratch buffers
    std::vector<T> ref_verts_sub;   // reference subcell vertices
    std::vector<T> phys_verts_sub;  // physical subcell vertices
    std::vector<T> mapped_ref_pts;  // canonical pts mapped to ref subcell

    // Helper: process one simplex sub-sub-cell (triangle or tetrahedron).
    // ref_verts_sub / phys_verts_sub must already contain the nv*gdim vertex
    // coordinates of the simplex. stype must be interval/triangle/tetrahedron.
    auto process_simplex = [&](cutcells::cell::type stype_)
    {
        const auto& ref_rule = get_reference_rule<T>(stype_, order);
        const int   nq       = ref_rule._num_points;
        const int   sdim     = ref_rule._tdim;   // canonical simplex tdim

        // ref_verts_sub contains vertices in the parent reference space, which
        // has dimension gdim (= tdim for non-embedded cells). We must output
        // gdim-dimensional parent reference coordinates, not sdim-dimensional
        // canonical coordinates. For interface subcells sdim < gdim (e.g.
        // interval sdim=1 from a triangle parent gdim=2), using sdim as the
        // output stride would only read the first component of each vertex
        // and produce wrong parent-reference points for physical_points.
        mapped_ref_pts.resize(static_cast<std::size_t>(nq) * gdim);
        map_canonical_to_subcell(ref_rule._points.data(), nq, sdim,
                                 ref_verts_sub.data(), gdim,
                                 stype_,
                                 mapped_ref_pts.data());

        rules._points.insert(rules._points.end(),
                             mapped_ref_pts.begin(), mapped_ref_pts.end());

        const T det_phys = simplex_physical_volume_factor(
            phys_verts_sub.data(), sdim, gdim);

        for (int q = 0; q < nq; ++q)
            rules._weights.push_back(ref_rule._weights[q] * det_phys);
    };

    // Helper: set ref_verts_sub and phys_verts_sub from a list of local vertex
    // indices into an already-filled parent buffer (ref_verts_parent /
    // phys_verts_parent), each of size nv_parent * gdim.
    auto set_tet_verts = [&](const std::vector<T>& ref_parent,
                              const std::vector<T>& phys_parent,
                              const std::vector<int>& tet_local_ids)
    {
        const int nv_tet = static_cast<int>(tet_local_ids.size());
        ref_verts_sub.resize(static_cast<std::size_t>(nv_tet) * gdim);
        phys_verts_sub.resize(static_cast<std::size_t>(nv_tet) * gdim);
        for (int j = 0; j < nv_tet; ++j)
        {
            const int vi = tet_local_ids[j];
            for (int d = 0; d < gdim; ++d)
            {
                ref_verts_sub[j * gdim + d]  = ref_parent[vi * gdim + d];
                phys_verts_sub[j * gdim + d] = phys_parent[vi * gdim + d];
            }
        }
    };

    for (int s = 0; s < nc; ++s)
    {
        const cutcells::cell::type stype = cut_cell._types[s];

        // For non-simplex types (prism, pyramid, quadrilateral) we
        // sub-triangulate into simplices to keep the integration exact.
        // The sub-triangulation uses local indices into the subcell vertex list.
        if (stype == cutcells::cell::type::prism   ||
            stype == cutcells::cell::type::pyramid  ||
            stype == cutcells::cell::type::quadrilateral)
        {
            // Gather all vertices of this subcell into temporary arrays.
            std::vector<T> ref_parent, phys_parent;
            gather_subcell_verts(cut_cell, s, cut_cell._vertex_coords,      gdim, ref_parent);
            gather_subcell_verts(cut_cell, s, cut_cell._vertex_coords_phys, gdim, phys_parent);

            const int nv_sub = cutcells::cell::num_cell_vertices(cut_cell, s);
            std::vector<int> local_ids(nv_sub);
            for (int j = 0; j < nv_sub; ++j) local_ids[j] = j;

            // Obtain tet/triangle decomposition using existing triangulation.h
            std::vector<std::vector<int>> sub_simplices;
            cutcells::cell::triangulation(stype, local_ids.data(), sub_simplices);

            // Determine simplex type based on topological dimension
            const cutcells::cell::type simplex_type =
                (tdim == 3) ? cutcells::cell::type::tetrahedron
                            : cutcells::cell::type::triangle;

            for (const auto& simplex : sub_simplices)
            {
                set_tet_verts(ref_parent, phys_parent, simplex);
                process_simplex(simplex_type);
            }
            continue;
        }

        // Simplex path (interval, triangle, tetrahedron):
        gather_subcell_verts(cut_cell, s, cut_cell._vertex_coords,      gdim, ref_verts_sub);
        gather_subcell_verts(cut_cell, s, cut_cell._vertex_coords_phys, gdim, phys_verts_sub);
        process_simplex(stype);
    }

    // One rule entry per CutCell: record how many new points were added
    rules._parent_map.push_back(cut_cell._parent_cell_index);
    rules._debug_local_cell_id.push_back(std::int32_t(-1));
    rules._debug_chart_path.push_back(std::int32_t(0));
    rules._debug_refinement_depth.push_back(std::int32_t(0));
    rules._debug_chart_plan_hash.push_back(std::int64_t(0));
    rules._debug_candidate_mask.push_back(std::int32_t(0));
    rules._debug_rejection_reason.push_back(std::int32_t(0));
    rules._debug_measure_probe.push_back(T(0));
    rules._debug_validation_weight_sum.push_back(T(0));
    rules._offset.push_back(static_cast<int32_t>(rules._weights.size()));
}

// =============================================================================
// make_quadrature (in-place)
// =============================================================================
template <std::floating_point T>
void make_quadrature(const std::vector<cutcells::cell::CutCell<T>>& cut_cells,
                     int order,
                     QuadratureRules<T>& rules)
{
    rules._tdim = 0;
    rules._points.clear();
    rules._weights.clear();
    rules._offset.clear();
    rules._parent_map.clear();
    rules._debug_local_cell_id.clear();
    rules._debug_chart_path.clear();
    rules._debug_refinement_depth.clear();
    rules._debug_chart_plan_hash.clear();
    rules._debug_candidate_mask.clear();
    rules._debug_rejection_reason.clear();
    rules._debug_measure_probe.clear();
    rules._debug_validation_weight_sum.clear();

    for (std::size_t i = 0; i < cut_cells.size(); ++i)
    {
        append_quadrature(cut_cells[i], order, rules);
        // Override the parent_map entry written by append_quadrature to store
        // the position in the input vector rather than _parent_cell_index
        // (which may not be set by the caller).
        if (!rules._parent_map.empty())
            rules._parent_map.back() = static_cast<int32_t>(i);
    }
}

// =============================================================================
// make_quadrature (returns new object)
// =============================================================================
template <std::floating_point T>
QuadratureRules<T> make_quadrature(
    const std::vector<cutcells::cell::CutCell<T>>& cut_cells, int order)
{
    QuadratureRules<T> rules;
    make_quadrature(cut_cells, order, rules);
    return rules;
}

// =============================================================================
// runtime_quadrature
// =============================================================================
template <std::floating_point T>
QuadratureRules<T> runtime_quadrature(
    std::span<const T>   ls_vals,
    std::span<const T>   points,
    std::span<const int> connectivity,
    std::span<const int> offset,
    std::span<const int> vtk_type,
    const std::string&   cut_type_str,
    bool triangulate,
    int  order)
{
    using namespace cutcells::cell;

    const int ncells     = static_cast<int>(vtk_type.size());
    const cut_type ctype_enum = string_to_cut_type(cut_type_str);

    // Thread-local scratch (reused across calls, avoids repeated allocations)
    thread_local std::vector<T>   tl_vtx_buf;
    thread_local std::vector<T>   tl_ls_buf;
    thread_local std::vector<T>   tl_vtx_buf_2d;   // 2D projection for embedded 2D cells
    thread_local CutCell<T>       tl_scratch;
    thread_local std::vector<T>   tl_can_verts;
    thread_local std::vector<T>   tl_ref_sub;
    thread_local std::vector<T>   tl_phys_sub;
    thread_local std::vector<T>   tl_mapped_ref;

    QuadratureRules<T> rules;
    rules._offset.push_back(0);

    for (int ci = 0; ci < ncells; ++ci)
    {
        const int coff = offset[ci];
        const type ctype = map_vtk_type_to_cell_type(
            static_cast<vtk_types>(vtk_type[ci]));
        const int nv   = get_num_vertices(ctype);
        const int tdim = get_tdim(ctype);

        // Gather level-set values at the cell vertices
        tl_ls_buf.resize(nv);
        for (int j = 0; j < nv; ++j)
            tl_ls_buf[j] = ls_vals[connectivity[coff + j]];

        const domain dom = classify_cell_domain<T>(
            std::span<const T>(tl_ls_buf.data(), nv));

        const bool is_full = (ctype_enum == cut_type::philt0 && dom == domain::inside)  ||
                             (ctype_enum == cut_type::phigt0 && dom == domain::outside);
        const bool is_cut  = (ctype_enum != cut_type::unset  && dom == domain::intersected);

        if (!is_full && !is_cut)
            continue;

        // Latch output dimension on first contributing cell
        if (rules._tdim == 0)
            rules._tdim = tdim;

        // Gather physical vertex coordinates (gdim = 3, VTK convention)
        tl_vtx_buf.resize(nv * 3);
        for (int j = 0; j < nv; ++j)
        {
            const int vid = connectivity[coff + j];
            tl_vtx_buf[j * 3 + 0] = points[vid * 3 + 0];
            tl_vtx_buf[j * 3 + 1] = points[vid * 3 + 1];
            tl_vtx_buf[j * 3 + 2] = points[vid * 3 + 2];
        }

        // ------------------------------------------------------------------
        // Full-cell path
        // ------------------------------------------------------------------
        if (is_full)
        {
            const bool is_simplex = (ctype == type::interval  ||
                                     ctype == type::triangle   ||
                                     ctype == type::tetrahedron);

            if (!triangulate || is_simplex)
            {
                // Direct rule: reference points are canonical quadrature points
                const auto& ref_rule = get_reference_rule<T>(ctype, order);
                const T det_J = affine_volume_factor<T>(ctype, tl_vtx_buf.data(), 3);
                const int nq  = ref_rule._num_points;

                rules._points.insert(rules._points.end(),
                    ref_rule._points.begin(), ref_rule._points.end());
                for (int q = 0; q < nq; ++q)
                    rules._weights.push_back(ref_rule._weights[q] * det_J);
            }
            else
            {
                // Triangulate into simplices; integrate over each sub-simplex
                tl_can_verts = canonical_vertices<T>(ctype);

                std::vector<int> local_ids(nv);
                for (int j = 0; j < nv; ++j) local_ids[j] = j;

                std::vector<std::vector<int>> sub_simplices;
                triangulation(ctype, local_ids.data(), sub_simplices);

                const type simplex_type =
                    (tdim == 3) ? type::tetrahedron : type::triangle;
                const int sv = (tdim == 3) ? 4 : 3;   // vertices per simplex

                const auto& ref_rule_s = get_reference_rule<T>(simplex_type, order);
                const int nq = ref_rule_s._num_points;
                tl_mapped_ref.resize(nq * tdim);

                for (const auto& simplex : sub_simplices)
                {
                    // Gather ref-space vertices of this sub-simplex (tdim-D coords)
                    tl_ref_sub.resize(sv * tdim);
                    for (int k = 0; k < sv; ++k)
                    {
                        const int vi = simplex[k];
                        for (int d = 0; d < tdim; ++d)
                            tl_ref_sub[k * tdim + d] = tl_can_verts[vi * tdim + d];
                    }

                    // Gather physical vertices of this sub-simplex (3-D coords)
                    tl_phys_sub.resize(sv * 3);
                    for (int k = 0; k < sv; ++k)
                    {
                        const int vi = simplex[k];
                        for (int d = 0; d < 3; ++d)
                            tl_phys_sub[k * 3 + d] = tl_vtx_buf[vi * 3 + d];
                    }

                    // Map canonical simplex quadrature points → parent ref space
                    map_canonical_to_subcell(
                        ref_rule_s._points.data(), nq,
                        tdim,
                        tl_ref_sub.data(), tdim,
                        simplex_type,
                        tl_mapped_ref.data());

                    // Physical volume factor for this sub-simplex
                    const T det = simplex_physical_volume_factor(
                        tl_phys_sub.data(), tdim, 3);

                    rules._points.insert(rules._points.end(),
                        tl_mapped_ref.begin(), tl_mapped_ref.end());
                    for (int q = 0; q < nq; ++q)
                        rules._weights.push_back(ref_rule_s._weights[q] * det);
                }
            }

            rules._parent_map.push_back(static_cast<int32_t>(ci));
            rules._debug_local_cell_id.push_back(std::int32_t(-1));
            rules._debug_chart_path.push_back(std::int32_t(0));
            rules._debug_refinement_depth.push_back(std::int32_t(0));
            rules._debug_chart_plan_hash.push_back(std::int64_t(0));
            rules._debug_candidate_mask.push_back(std::int32_t(0));
            rules._debug_rejection_reason.push_back(std::int32_t(0));
            rules._debug_measure_probe.push_back(T(0));
            rules._debug_validation_weight_sum.push_back(T(0));
            rules._offset.push_back(static_cast<int32_t>(rules._weights.size()));
        }
        // ------------------------------------------------------------------
        // Cut-cell path
        // ------------------------------------------------------------------
        else  // is_cut
        {
            if (tdim == 1)
                continue;  // embedded 1D edges cannot be cut

            if (tdim == 2)
            {
                // For 2D cells embedded in a 3D VTK mesh (z-plane), project
                // vertices to 2D (x,y) so the cutter and pull-back work in a
                // square (2×2) Jacobian system.
                tl_vtx_buf_2d.resize(nv * 2);
                for (int j = 0; j < nv; ++j)
                {
                    tl_vtx_buf_2d[j * 2 + 0] = tl_vtx_buf[j * 3 + 0];
                    tl_vtx_buf_2d[j * 2 + 1] = tl_vtx_buf[j * 3 + 1];
                }
                cell::cut<T>(ctype,
                             std::span<const T>(tl_vtx_buf_2d.data(), nv * 2),
                             2,
                             std::span<const T>(tl_ls_buf.data(), nv),
                             cut_type_str,
                             tl_scratch,
                             triangulate);
                if (cell::num_cells(tl_scratch) == 0)
                    continue;
                tl_scratch._parent_cell_type     = ctype;
                tl_scratch._parent_vertex_coords = tl_vtx_buf_2d;  // 2D parent
                cell::complete_from_physical(tl_scratch);
                append_quadrature(tl_scratch, order, rules);
            }
            else  // tdim == 3
            {
                cell::cut<T>(ctype,
                             std::span<const T>(tl_vtx_buf.data(), nv * 3),
                             3,
                             std::span<const T>(tl_ls_buf.data(), nv),
                             cut_type_str,
                             tl_scratch,
                             triangulate);
                if (cell::num_cells(tl_scratch) == 0)
                    continue;
                tl_scratch._parent_cell_type      = ctype;
                tl_scratch._parent_vertex_coords  = tl_vtx_buf;
                cell::complete_from_physical(tl_scratch);
                append_quadrature(tl_scratch, order, rules);
            }

            // append_quadrature records tl_scratch._parent_cell_index; override
            // with the global mesh index
            if (!rules._parent_map.empty())
                rules._parent_map.back() = static_cast<int32_t>(ci);
        }
    }

    return rules;
}

// =============================================================================
// physical_points
// =============================================================================
template <std::floating_point T>
std::vector<T> physical_points(
    const QuadratureRules<T>&  rules,
    std::span<const T>   points,
    std::span<const int> connectivity,
    std::span<const int> offset,
    std::span<const int> vtk_type)
{
    using namespace cutcells::cell;

    const int nrules    = static_cast<int>(rules._parent_map.size());
    const int total_pts = static_cast<int>(rules._weights.size());
    const int tdim      = rules._tdim;

    std::vector<T> out(static_cast<std::size_t>(total_pts) * 3, T(0));

    if (nrules == 0 || total_pts == 0 || tdim == 0)
        return out;

    thread_local std::vector<T> tl_phys_verts;

    for (int i = 0; i < nrules; ++i)
    {
        const int q_begin = rules._offset[i];
        const int q_end   = rules._offset[i + 1];
        const int nq      = q_end - q_begin;
        if (nq == 0)
            continue;

        const int ci   = rules._parent_map[i];
        const int coff = offset[ci];
        const type ctype = map_vtk_type_to_cell_type(
            static_cast<vtk_types>(vtk_type[ci]));
        if (get_tdim(ctype) != tdim)
            throw std::invalid_argument(
                "physical_points: quadrature point dimension does not match "
                "the supplied VTK parent cell type");
        const int nv = get_num_vertices(ctype);

        // Gather physical vertices for this cell
        tl_phys_verts.resize(nv * 3);
        for (int j = 0; j < nv; ++j)
        {
            const int vid = connectivity[coff + j];
            tl_phys_verts[j * 3 + 0] = points[vid * 3 + 0];
            tl_phys_verts[j * 3 + 1] = points[vid * 3 + 1];
            tl_phys_verts[j * 3 + 2] = points[vid * 3 + 2];
        }

        // Push parent reference points to physical coordinates with the same
        // affine column convention used by the cut-cell maps.
        {
            const auto cols  = jacobian_col_indices(ctype);
            const T*   v0    = tl_phys_verts.data();
            const T*   Xref  = rules._points.data() + q_begin * tdim;
            T*         xout  = out.data() + q_begin * 3;
            for (int q = 0; q < nq; ++q)
            {
                for (int d = 0; d < 3; ++d)
                    xout[q * 3 + d] = v0[d];
                for (int k = 0; k < tdim; ++k)
                {
                    const T*   vk   = tl_phys_verts.data() + cols[k] * 3;
                    const T    xi_k = Xref[q * tdim + k];
                    for (int d = 0; d < 3; ++d)
                        xout[q * 3 + d] += (vk[d] - v0[d]) * xi_k;
                }
            }
        }
    }

    return out;
}

// =============================================================================
// Explicit instantiations
// =============================================================================
template void append_quadrature(const cutcells::cell::CutCell<double>&, int, QuadratureRules<double>&);
template void append_quadrature(const cutcells::cell::CutCell<float>&,  int, QuadratureRules<float>&);

template void make_quadrature(const std::vector<cutcells::cell::CutCell<double>>&, int, QuadratureRules<double>&);
template void make_quadrature(const std::vector<cutcells::cell::CutCell<float>>&,  int, QuadratureRules<float>&);

template QuadratureRules<double> make_quadrature(const std::vector<cutcells::cell::CutCell<double>>&, int);
template QuadratureRules<float>  make_quadrature(const std::vector<cutcells::cell::CutCell<float>>&,  int);

template QuadratureRules<double> runtime_quadrature(
    std::span<const double>, std::span<const double>,
    std::span<const int>, std::span<const int>, std::span<const int>,
    const std::string&, bool, int);
template QuadratureRules<float>  runtime_quadrature(
    std::span<const float>,  std::span<const float>,
    std::span<const int>, std::span<const int>, std::span<const int>,
    const std::string&, bool, int);

template std::vector<double> physical_points(
    const QuadratureRules<double>&,
    std::span<const double>, std::span<const int>, std::span<const int>, std::span<const int>);
template std::vector<float>  physical_points(
    const QuadratureRules<float>&,
    std::span<const float>,  std::span<const int>, std::span<const int>, std::span<const int>);

} // namespace cutcells::quadrature
