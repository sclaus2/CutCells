// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "ho_cut_mesh.h"
#include "cell_certification.h"
#include "edge_certification.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cutcells
{
namespace
{

/// Gather vertex level-set values for the given cell.
template <std::floating_point T, std::integral I>
void gather_vertex_ls_values(const MeshView<T, I>& mesh,
                             const LevelSetFunction<T, I>& ls,
                             I cell_id, int nv,
                             std::vector<T>& out)
{
    out.resize(static_cast<std::size_t>(nv));

    if (ls.has_dof_values() && ls.has_mesh_data())
    {
        auto dofs = ls.mesh_data->cell_dofs_span(cell_id);
        for (int v = 0; v < nv; ++v)
            out[static_cast<std::size_t>(v)] =
                ls.dof_values[static_cast<std::size_t>(
                    dofs[static_cast<std::size_t>(v)])];
    }
    else if (ls.has_value())
    {
        for (int v = 0; v < nv; ++v)
        {
            const I node_id = mesh.cell_node(cell_id, static_cast<I>(v));
            out[static_cast<std::size_t>(v)] =
                ls.value(mesh.node(node_id), cell_id);
        }
    }
    else
    {
        throw std::runtime_error(
            "cut: LevelSetFunction has neither nodal_values, "
            "dof_values, nor an analytical value function");
    }
}

template <std::floating_point T>
T bernstein_cell_sign_tol(std::span<const T> coeffs)
{
    T max_abs = 0;
    for (const T c : coeffs)
        max_abs = std::max(max_abs, std::fabs(c));

    return std::max(T(64) * std::numeric_limits<T>::epsilon()
                        * std::max(T(1), max_abs),
                    T(1e-14));
}

template <std::floating_point T>
cell::domain classify_cell_domain_from_bernstein(std::span<const T> coeffs)
{
    if (coeffs.empty())
        return cell::domain::unset;

    const T sign_tol = bernstein_cell_sign_tol(coeffs);
    if (bernstein_all_positive(coeffs, sign_tol))
        return cell::domain::outside;
    if (bernstein_all_negative(coeffs, sign_tol))
        return cell::domain::inside;
    return cell::domain::intersected;
}

template <std::floating_point T, std::integral I>
cell::domain classify_cell_domain_fast(const MeshView<T, I>& mesh,
                                       const LevelSetFunction<T, I>& ls,
                                       I cell_id, int nv,
                                       bool use_bernstein_classification,
                                       std::vector<T>& vertex_ls_values,
                                       LevelSetCell<T, I>* intersected_ls_cell)
{
    if (use_bernstein_classification)
    {
        LevelSetCell<T, I> ls_cell = make_cell_level_set(ls, cell_id);
        const cell::domain dom = classify_cell_domain_from_bernstein<T>(
            std::span<const T>(ls_cell.bernstein_coeffs.data(),
                               ls_cell.bernstein_coeffs.size()));
        if (dom != cell::domain::unset)
        {
            if (dom == cell::domain::intersected
                && intersected_ls_cell != nullptr)
            {
                *intersected_ls_cell = std::move(ls_cell);
            }
            return dom;
        }
    }

    gather_vertex_ls_values(mesh, ls, cell_id, nv, vertex_ls_values);
    return cell::classify_cell_domain<T>(
        std::span<const T>(vertex_ls_values.data(),
                           static_cast<std::size_t>(nv)));
}

} // anonymous namespace

// =====================================================================
// cut() — single level set
// =====================================================================

template <std::floating_point T, std::integral I>
std::pair<HOCutCells<T, I>, BackgroundMeshData<T, I>>
cut(const MeshView<T, I>& mesh, const LevelSetFunction<T, I>& ls)
{
    if (!mesh.has_cell_types())
        throw std::runtime_error("cut: MeshView must have cell types");

    const I ncells = mesh.num_cells();

    // --- BackgroundMeshData ---
    BackgroundMeshData<T, I> bg;
    bg.mesh = &mesh;
    bg.level_set_names.push_back(ls.name);
    bg.num_cells = static_cast<int>(ncells);
    bg.num_level_sets = 1;
    bg.cell_domains.assign(static_cast<std::size_t>(ncells),
                            cell::domain::unset);
    bg.cell_to_cut_index.assign(static_cast<std::size_t>(ncells), -1);

    // --- HOCutCells ---
    HOCutCells<T, I> hc;
    hc.gdim = mesh.gdim;
    hc.tdim = cell::get_tdim(mesh.cell_type(I(0)));
    hc.ls_offsets.push_back(0);

    const bool use_bernstein_classification =
        ls.type == LevelSetType::Polynomial
        && ls.has_mesh_data()
        && ls.has_dof_values();
    std::vector<T> ls_vertex_vals;

    for (I ci = 0; ci < ncells; ++ci)
    {
        const cell::type ctype = mesh.cell_type(ci);
        const int nv = cell::get_num_vertices(ctype);
        LevelSetCell<T, I> ls_cell;
        const cell::domain dom = classify_cell_domain_fast(
            mesh, ls, ci, nv, use_bernstein_classification, ls_vertex_vals,
            &ls_cell);
        bg.cell_domains[static_cast<std::size_t>(ci)] = dom;

        if (dom != cell::domain::intersected)
            continue;

        const int cut_idx = hc.num_cut_cells();
        bg.cell_to_cut_index[static_cast<std::size_t>(ci)] = cut_idx;

        // LevelSetCell
        hc.level_set_cells.push_back(std::move(ls_cell));
        hc.ls_offsets.push_back(
            static_cast<int>(hc.level_set_cells.size()));

        // AdaptCell
        AdaptCell<T> ac = make_adapt_cell(mesh, ci);

        certify_refine_and_process_ready_cells(
            ac, hc.level_set_cells.back(), /*level_set_id=*/0,
            /*max_iterations=*/8, T(1e-12), T(1e-12), /*edge_max_depth=*/20);
        hc.adapt_cells.push_back(std::move(ac));

        // Single LS: bit 0 is always set.
        hc.active_level_set_mask.push_back(std::uint64_t(1));

        hc.parent_cell_ids.push_back(ci);
    }

    return {std::move(hc), std::move(bg)};
}

// =====================================================================
// cut() — multiple level sets
// =====================================================================

template <std::floating_point T, std::integral I>
std::pair<HOCutCells<T, I>, BackgroundMeshData<T, I>>
cut(const MeshView<T, I>& mesh,
    const std::vector<LevelSetFunction<T, I>>& level_sets)
{
    if (!mesh.has_cell_types())
        throw std::runtime_error("cut: MeshView must have cell types");
    if (level_sets.empty())
        throw std::runtime_error("cut: no level sets provided");
    if (level_sets.size() > 64)
        throw std::runtime_error("cut: more than 64 level sets not supported");

    const I ncells = mesh.num_cells();
    const int nls = static_cast<int>(level_sets.size());

    // --- BackgroundMeshData ---
    BackgroundMeshData<T, I> bg;
    bg.mesh = &mesh;
    bg.num_cells = static_cast<int>(ncells);
    bg.num_level_sets = nls;
    bg.cell_domains.assign(
        static_cast<std::size_t>(nls) * static_cast<std::size_t>(ncells),
        cell::domain::unset);
    bg.cell_to_cut_index.assign(static_cast<std::size_t>(ncells), -1);
    for (const auto& ls : level_sets)
        bg.level_set_names.push_back(ls.name);

    // --- HOCutCells ---
    HOCutCells<T, I> hc;
    hc.gdim = mesh.gdim;
    hc.tdim = cell::get_tdim(mesh.cell_type(I(0)));
    hc.ls_offsets.push_back(0);

    std::vector<bool> use_bernstein_classification(
        static_cast<std::size_t>(nls), false);
    for (int li = 0; li < nls; ++li)
    {
        const auto& ls = level_sets[static_cast<std::size_t>(li)];
        use_bernstein_classification[static_cast<std::size_t>(li)] =
            ls.type == LevelSetType::Polynomial
            && ls.has_mesh_data()
            && ls.has_dof_values();
    }
    std::vector<T> ls_vertex_vals;

    for (I ci = 0; ci < ncells; ++ci)
    {
        const cell::type ctype = mesh.cell_type(ci);
        const int nv = cell::get_num_vertices(ctype);

        // Classify each level set individually; track if any intersects.
        bool any_intersected = false;
        std::vector<int> intersected_ls_indices;
        std::vector<LevelSetCell<T, I>> intersected_ls_cells;
        intersected_ls_indices.reserve(static_cast<std::size_t>(nls));
        intersected_ls_cells.reserve(static_cast<std::size_t>(nls));
        for (int li = 0; li < nls; ++li)
        {
            LevelSetCell<T, I> ls_cell;
            const cell::domain dom = classify_cell_domain_fast(
                mesh, level_sets[static_cast<std::size_t>(li)], ci, nv,
                use_bernstein_classification[static_cast<std::size_t>(li)],
                ls_vertex_vals, &ls_cell);

            bg.cell_domains[static_cast<std::size_t>(
                li * static_cast<int>(ncells) + static_cast<int>(ci))] = dom;

            if (dom == cell::domain::intersected)
            {
                any_intersected = true;
                intersected_ls_indices.push_back(li);
                intersected_ls_cells.push_back(std::move(ls_cell));
            }
        }

        if (!any_intersected)
            continue;

        const int cut_idx = hc.num_cut_cells();
        bg.cell_to_cut_index[static_cast<std::size_t>(ci)] = cut_idx;

        // Build AdaptCell once per cell.
        AdaptCell<T> ac = make_adapt_cell(mesh, ci);

        // Build LevelSetCell and fill vertex signs for each intersecting LS.
        std::uint64_t cell_active_mask = 0;
        for (std::size_t k = 0; k < intersected_ls_indices.size(); ++k)
        {
            const int li = intersected_ls_indices[k];
            hc.level_set_cells.push_back(std::move(intersected_ls_cells[k]));
            certify_refine_and_process_ready_cells(
                ac, hc.level_set_cells.back(), li,
                /*max_iterations=*/8, T(1e-12), T(1e-12), /*edge_max_depth=*/20);

            cell_active_mask |= std::uint64_t(1) << li;
        }
        hc.ls_offsets.push_back(
            static_cast<int>(hc.level_set_cells.size()));
        hc.active_level_set_mask.push_back(cell_active_mask);

        hc.adapt_cells.push_back(std::move(ac));
        hc.parent_cell_ids.push_back(ci);
    }

    return {std::move(hc), std::move(bg)};
}

// =====================================================================
// select_part()
// =====================================================================

template <std::floating_point T, std::integral I>
HOMeshPart<T, I> select_part(const HOCutCells<T, I>& cut_cells,
                              const BackgroundMeshData<T, I>& bg,
                              std::string_view expr_str)
{
    HOMeshPart<T, I> part;
    part.cut_cells = &cut_cells;
    part.bg        = &bg;

    // Parse and compile the expression.
    part.expr = parse_selection_expr(expr_str);
    compile_selection_expr(part.expr, bg.level_set_names);
    part.dim = infer_selection_dim(part.expr, cut_cells.tdim);

    const bool is_volume = (part.dim == cut_cells.tdim);

    // --- Uncut cells: scan cell_domains ---
    const int ncells = bg.num_cells;
    for (I ci = 0; ci < static_cast<I>(ncells); ++ci)
    {
        if (bg.cell_to_cut_index[static_cast<std::size_t>(ci)] >= 0)
            continue; // cut cell, handled below

        // For an uncut cell, every vertex has the same sign for each LS.
        // Check if the cell domain is compatible with the expression.
        bool match = true;
        for (const auto& clause : part.expr.clauses)
        {
            const int li = clause.level_set_index;
            const cell::domain dom = bg.domain(li, ci);

            switch (clause.relation)
            {
                case Relation::LessThan:
                    if (dom != cell::domain::inside) match = false;
                    break;
                case Relation::GreaterThan:
                    if (dom != cell::domain::outside) match = false;
                    break;
                case Relation::EqualTo:
                    // An uncut (non-intersected) cell has no zero interface.
                    match = false;
                    break;
            }
            if (!match) break;
        }

        if (match)
            part.uncut_cell_ids.push_back(ci);
    }

    // --- Cut cells: check vertex bitmasks for volume selections,
    //     or entity inventory for interface selections ---
    const int ncut = cut_cells.num_cut_cells();
    for (int k = 0; k < ncut; ++k)
    {
        const auto& ac = cut_cells.adapt_cells[static_cast<std::size_t>(k)];
        const int nv = ac.n_vertices();

        if (is_volume)
        {
            // Volume selection: a cut cell contributes if at least one vertex
            // satisfies all sign constraints.
            // Actually, a cut cell always contributes to a volume part
            // (it has sub-entities on the requested side), so we include it.
            // More refined: check if the required bitmask pattern is present.
            bool has_matching_vertex = false;
            for (int v = 0; v < nv; ++v)
            {
                const auto zm = ac.zero_mask_per_vertex[static_cast<std::size_t>(v)];
                const auto nm = ac.negative_mask_per_vertex[static_cast<std::size_t>(v)];

                // Positive bit: not zero, not negative.
                const auto pm = ~zm & ~nm;

                bool vertex_ok = true;

                // Check if zero requirements are met (for EqualTo clauses).
                if ((zm & part.expr.zero_required) != part.expr.zero_required)
                    vertex_ok = false;

                // Check if negative requirements are met.
                if ((nm & part.expr.negative_required) != part.expr.negative_required)
                    vertex_ok = false;

                // Check if positive requirements are met.
                if ((pm & part.expr.positive_required) != part.expr.positive_required)
                    vertex_ok = false;

                if (vertex_ok)
                {
                    has_matching_vertex = true;
                    break;
                }
            }

            if (has_matching_vertex)
                part.cut_cell_ids.push_back(static_cast<std::int32_t>(k));
        }
        else
        {
            // Interface selection: for now include any cut cell that has
            // the required zero level set(s) intersecting it.
            // A cell is included if:
            //  - the zero_required LSs actually cut the cell
            //  - the sign constraints on the non-zero LSs are satisfiable
            //    (at least one vertex per non-zero constraint).
            bool has_zero_ls = true;
            for (int v = 0; v < nv && has_zero_ls; ++v)
            {
                // Check eventually; for now, if the LS is marked as
                // intersected in cell_domains, the zero surface exists.
            }

            // Check that the intersecting LS indices match zero_required.
            bool zero_ok = true;
            for (const auto& clause : part.expr.clauses)
            {
                if (clause.relation != Relation::EqualTo) continue;
                const int li = clause.level_set_index;
                const I bg_cell = cut_cells.parent_cell_ids[static_cast<std::size_t>(k)];
                if (bg.domain(li, bg_cell) != cell::domain::intersected)
                {
                    zero_ok = false;
                    break;
                }
            }

            if (zero_ok)
                part.cut_cell_ids.push_back(static_cast<std::int32_t>(k));
        }
    }

    part.cut_only = !is_volume;

    return part;
}

// ---------------------------------------------------------------------------
// Explicit template instantiations
// ---------------------------------------------------------------------------

// cut() single LS
template std::pair<HOCutCells<double, int>, BackgroundMeshData<double, int>>
cut(const MeshView<double, int>&, const LevelSetFunction<double, int>&);

template std::pair<HOCutCells<float, int>, BackgroundMeshData<float, int>>
cut(const MeshView<float, int>&, const LevelSetFunction<float, int>&);

template std::pair<HOCutCells<double, long>, BackgroundMeshData<double, long>>
cut(const MeshView<double, long>&, const LevelSetFunction<double, long>&);

template std::pair<HOCutCells<float, long>, BackgroundMeshData<float, long>>
cut(const MeshView<float, long>&, const LevelSetFunction<float, long>&);

// cut() multi LS
template std::pair<HOCutCells<double, int>, BackgroundMeshData<double, int>>
cut(const MeshView<double, int>&, const std::vector<LevelSetFunction<double, int>>&);

template std::pair<HOCutCells<float, int>, BackgroundMeshData<float, int>>
cut(const MeshView<float, int>&, const std::vector<LevelSetFunction<float, int>>&);

template std::pair<HOCutCells<double, long>, BackgroundMeshData<double, long>>
cut(const MeshView<double, long>&, const std::vector<LevelSetFunction<double, long>>&);

template std::pair<HOCutCells<float, long>, BackgroundMeshData<float, long>>
cut(const MeshView<float, long>&, const std::vector<LevelSetFunction<float, long>>&);

// select_part()
template HOMeshPart<double, int>
select_part(const HOCutCells<double, int>&, const BackgroundMeshData<double, int>&,
            std::string_view);

template HOMeshPart<float, int>
select_part(const HOCutCells<float, int>&, const BackgroundMeshData<float, int>&,
            std::string_view);

template HOMeshPart<double, long>
select_part(const HOCutCells<double, long>&, const BackgroundMeshData<double, long>&,
            std::string_view);

template HOMeshPart<float, long>
select_part(const HOCutCells<float, long>&, const BackgroundMeshData<float, long>&,
            std::string_view);

} // namespace cutcells
