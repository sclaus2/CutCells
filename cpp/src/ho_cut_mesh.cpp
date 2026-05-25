// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "ho_cut_mesh.h"
#include "cell_certification.h"
#include "edge_certification.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
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
                             std::vector<T>& out,
                             std::vector<I>& mesh_cell_scratch,
                             std::vector<I>& ls_dof_scratch)
{
    out.resize(static_cast<std::size_t>(nv));

    if (ls.has_dof_values() && ls.has_mesh_data())
    {
        auto dofs = ls.mesh_data.cell_dofs_span(cell_id, ls_dof_scratch);
        for (int v = 0; v < nv; ++v)
            out[static_cast<std::size_t>(v)] =
                ls.dof_values[static_cast<std::size_t>(
                    dofs[static_cast<std::size_t>(v)])];
    }
    else if (ls.has_value())
    {
        auto nodes = mesh.cell_nodes(cell_id, mesh_cell_scratch);
        for (int v = 0; v < nv; ++v)
        {
            const I node_id = nodes[static_cast<std::size_t>(v)];
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
                                       std::vector<I>& mesh_cell_scratch,
                                       std::vector<I>& ls_dof_scratch,
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

    gather_vertex_ls_values(mesh, ls, cell_id, nv, vertex_ls_values,
                            mesh_cell_scratch, ls_dof_scratch);
    return cell::classify_cell_domain<T>(
        std::span<const T>(vertex_ls_values.data(),
                           static_cast<std::size_t>(nv)));
}

template <std::floating_point T>
bool cell_contains_all_vertices(std::span<const std::int32_t> cell_verts,
                                std::span<const int> entity_verts)
{
    for (const int v : entity_verts)
    {
        bool found = false;
        for (const auto cv : cell_verts)
        {
            if (cv == v)
            {
                found = true;
                break;
            }
        }
        if (!found)
            return false;
    }
    return true;
}

template <std::floating_point T>
bool vertex_state_for_level_set(const AdaptCell<T>& ac,
                                int vertex_id,
                                int level_set_id,
                                bool& is_negative,
                                bool& is_positive,
                                bool& is_zero)
{
    if (level_set_id < 0 || level_set_id >= 64)
        return false;

    const std::uint64_t bit = std::uint64_t(1) << level_set_id;
    const auto zm = ac.zero_mask_per_vertex[static_cast<std::size_t>(vertex_id)];
    const auto nm = ac.negative_mask_per_vertex[static_cast<std::size_t>(vertex_id)];
    is_zero = (zm & bit) != 0;
    is_negative = !is_zero && ((nm & bit) != 0);
    is_positive = !is_zero && !is_negative;
    return true;
}

template <std::floating_point T, std::integral I>
bool leaf_cell_matches_sign_requirements(
    const AdaptCell<T>& ac,
    std::span<const std::int32_t> cell_verts,
    const SelectionTerm& term,
    std::uint64_t cut_cell_active_mask,
    const ParentCellClassification<T, I>& parent_cells,
    I parent_cell_id)
{
    const int nls = std::min(parent_cells.num_level_sets, 64);
    for (int li = 0; li < nls; ++li)
    {
        const std::uint64_t bit = std::uint64_t(1) << li;
        const bool require_neg = (term.negative_required & bit) != 0;
        const bool require_pos = (term.positive_required & bit) != 0;
        if (!require_neg && !require_pos)
            continue;

        const bool ls_is_active = (cut_cell_active_mask & bit) != 0;
        if (!ls_is_active)
        {
            const auto dom = parent_cells.domain(li, parent_cell_id);
            if (require_neg && dom != cell::domain::inside)
                return false;
            if (require_pos && dom != cell::domain::outside)
                return false;
            continue;
        }

        bool has_neg = false;
        bool has_pos = false;
        for (const auto cv : cell_verts)
        {
            bool is_neg = false;
            bool is_pos = false;
            bool is_zero = false;
            vertex_state_for_level_set(ac, static_cast<int>(cv), li, is_neg, is_pos, is_zero);
            has_neg = has_neg || is_neg;
            has_pos = has_pos || is_pos;
        }

        if (require_neg)
        {
            if (has_pos || !has_neg)
                return false;
        }
        if (require_pos)
        {
            if (has_neg || !has_pos)
                return false;
        }
    }

    return true;
}

template <std::floating_point T, std::integral I>
bool leaf_cell_matches_selection_expr(
    const AdaptCell<T>& ac,
    std::span<const std::int32_t> cell_verts,
    const SelectionExpr& expr,
    std::uint64_t cut_cell_active_mask,
    const ParentCellClassification<T, I>& parent_cells,
    I parent_cell_id)
{
    for (const auto& term : expr.terms)
    {
        if (leaf_cell_matches_sign_requirements(
                ac, cell_verts, term, cut_cell_active_mask, parent_cells, parent_cell_id))
        {
            return true;
        }
    }
    return false;
}

template <std::floating_point T, std::integral I>
bool zero_entity_matches(
    const AdaptCell<T>& ac,
    int zero_entity_index,
    int target_dim,
    const SelectionTerm& term,
    std::uint64_t cut_cell_active_mask,
    const ParentCellClassification<T, I>& parent_cells,
    I parent_cell_id)
{
    if (ac.zero_entity_dim[static_cast<std::size_t>(zero_entity_index)] != target_dim)
        return false;

    const auto zero_mask = ac.zero_entity_zero_mask[static_cast<std::size_t>(zero_entity_index)];
    if ((cut_cell_active_mask & term.zero_required) != term.zero_required)
        return false;

    if ((zero_mask & term.zero_required) != term.zero_required)
        return false;

    if (term.negative_required == 0 && term.positive_required == 0)
        return true;

    std::vector<int> zero_verts;
    const int zdim = ac.zero_entity_dim[static_cast<std::size_t>(zero_entity_index)];
    const int zid = ac.zero_entity_id[static_cast<std::size_t>(zero_entity_index)];
    if (zdim == 0)
    {
        zero_verts.push_back(zid);
    }
    else
    {
        auto verts = ac.entity_to_vertex[zdim][static_cast<std::int32_t>(zid)];
        zero_verts.reserve(verts.size());
        for (const auto v : verts)
            zero_verts.push_back(static_cast<int>(v));
    }

    const int tdim = ac.tdim;
    const int n_cells = ac.n_entities(tdim);
    for (int c = 0; c < n_cells; ++c)
    {
        auto cell_verts = ac.entity_to_vertex[tdim][static_cast<std::int32_t>(c)];
        if (!cell_contains_all_vertices<T>(cell_verts, std::span<const int>(zero_verts)))
            continue;

        if (leaf_cell_matches_sign_requirements(
                ac, cell_verts, term, cut_cell_active_mask, parent_cells, parent_cell_id))
        {
            return true;
        }
    }

    return false;
}

template <std::floating_point T, std::integral I>
bool zero_entity_matches_selection_expr(
    const AdaptCell<T>& ac,
    int zero_entity_index,
    int target_dim,
    const SelectionExpr& expr,
    std::uint64_t cut_cell_active_mask,
    const ParentCellClassification<T, I>& parent_cells,
    I parent_cell_id)
{
    for (const auto& term : expr.terms)
    {
        if (zero_entity_matches(
                ac, zero_entity_index, target_dim, term,
                cut_cell_active_mask, parent_cells, parent_cell_id))
        {
            return true;
        }
    }
    return false;
}

template <std::floating_point T, std::integral I>
void refresh_adapt_cell_semantics(
    AdaptCell<T>& ac,
    std::span<const int> processed_level_set_indices,
    std::span<const LevelSetCell<T, I>> processed_level_set_cells,
    int total_num_level_sets,
    T zero_tol)
{
    for (std::size_t i = 0; i < processed_level_set_indices.size(); ++i)
    {
        fill_all_vertex_signs_from_level_set(
            ac,
            processed_level_set_cells[i],
            processed_level_set_indices[i],
            zero_tol);
    }

    recompute_active_level_set_masks(ac, total_num_level_sets);
}

inline void validate_resolved_cut_options(const CutOptions& options)
{
    if (options.cut_approximation == "linear")
    {
        if (options.cut_approximation_order != 1)
        {
            throw std::invalid_argument(
                "cut: cut_approximation='linear' requires "
                "cut_approximation_order=1");
        }
        return;
    }

    if (options.cut_approximation != "iso_p1")
    {
        throw std::invalid_argument(
            "cut: cut_approximation must be 'auto', 'linear', or 'iso_p1'");
    }

    if (options.cut_approximation_order < 1
        || options.cut_approximation_order > 4)
    {
        throw std::invalid_argument(
            "cut: cut_approximation_order must be 1, 2, 3, or 4");
    }
}

inline void validate_cut_options(const CutOptions& options)
{
    if (options.cut_approximation == "auto")
        return;

    validate_resolved_cut_options(options);
}

template <std::floating_point T, std::integral I>
int level_set_space_order(const LevelSetFunction<T, I>& ls)
{
    return ls.has_mesh_data() ? ls.mesh_data.degree : 1;
}

template <std::floating_point T, std::integral I>
CutOptions resolve_cut_options(const CutOptions& options,
                               const LevelSetFunction<T, I>& ls)
{
    validate_cut_options(options);
    if (options.cut_approximation != "auto")
        return options;

    CutOptions resolved = options;
    const int order = level_set_space_order(ls);
    if (order > 1)
    {
        resolved.cut_approximation = "iso_p1";
        resolved.cut_approximation_order = order;
    }
    else
    {
        resolved.cut_approximation = "linear";
        resolved.cut_approximation_order = 1;
    }
    validate_resolved_cut_options(resolved);
    return resolved;
}

template <std::floating_point T, std::integral I>
CutOptions resolve_cut_options(
    const CutOptions& options,
    const std::vector<LevelSetFunction<T, I>>& level_sets)
{
    validate_cut_options(options);
    if (options.cut_approximation != "auto")
        return options;

    int order = 1;
    for (const auto& ls : level_sets)
        order = std::max(order, level_set_space_order(ls));

    CutOptions resolved = options;
    if (order > 1)
    {
        resolved.cut_approximation = "iso_p1";
        resolved.cut_approximation_order = order;
    }
    else
    {
        resolved.cut_approximation = "linear";
        resolved.cut_approximation_order = 1;
    }
    validate_resolved_cut_options(resolved);
    return resolved;
}

template <std::floating_point T>
void apply_cut_approximation(AdaptCell<T>& ac, const CutOptions& options)
{
    validate_resolved_cut_options(options);
    if (options.cut_approximation != "iso_p1"
        || options.cut_approximation_order == 1)
    {
        return;
    }

    apply_iso_refine(
        ac, iso_p1_template(ac.parent_cell_type, options.cut_approximation_order));
}

} // anonymous namespace

// =====================================================================
// cut() — single level set
// =====================================================================

template <std::floating_point T, std::integral I>
std::pair<HOCutCells<T, I>, ParentCellClassification<T, I>>
cut(const MeshView<T, I>& mesh,
    const LevelSetFunction<T, I>& ls,
    const CutOptions& options)
{
    const CutOptions resolved_options = resolve_cut_options(options, ls);
    if (!mesh.has_cell_types())
        throw std::runtime_error("cut: MeshView must have cell types");

    const I ncells = mesh.num_cells();

    // --- ParentCellClassification ---
    ParentCellClassification<T, I> parent_cells;
    parent_cells.level_set_names.push_back(ls.name);
    parent_cells.num_cells = static_cast<int>(ncells);
    parent_cells.num_level_sets = 1;
    parent_cells.cell_domains.assign(static_cast<std::size_t>(ncells),
                            cell::domain::unset);
    parent_cells.cell_to_cut_index.assign(static_cast<std::size_t>(ncells), -1);

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
    std::vector<I> mesh_cell_scratch;
    std::vector<I> ls_dof_scratch;

    for (I ci = 0; ci < ncells; ++ci)
    {
        const cell::type ctype = mesh.cell_type(ci);
        const int nv = cell::get_num_vertices(ctype);
        LevelSetCell<T, I> ls_cell;
        const cell::domain dom = classify_cell_domain_fast(
            mesh, ls, ci, nv, use_bernstein_classification, ls_vertex_vals,
            mesh_cell_scratch, ls_dof_scratch,
            &ls_cell);
        parent_cells.cell_domains[static_cast<std::size_t>(ci)] = dom;

        if (dom != cell::domain::intersected)
            continue;

        const int cut_idx = hc.num_cut_cells();
        parent_cells.cell_to_cut_index[static_cast<std::size_t>(ci)] = cut_idx;

        // LevelSetCell
        hc.level_set_cells.push_back(std::move(ls_cell));
        hc.ls_offsets.push_back(
            static_cast<int>(hc.level_set_cells.size()));

        // AdaptCell
        AdaptCell<T> ac = make_adapt_cell(mesh, ci);
        apply_cut_approximation(ac, resolved_options);

        certify_refine_and_process_ready_cells(
            ac, hc.level_set_cells.back(), /*level_set_id=*/0,
            /*max_iterations=*/8, T(1e-12), T(1e-12), /*edge_max_depth=*/20,
            resolved_options.triangulate_cut_parts
                ? resolved_options.triangulation_strategy
                : cell::TriangulationStrategy::none);
        {
            const std::array<int, 1> processed_ids = {0};
            const auto* processed_cell = &hc.level_set_cells.back();
            refresh_adapt_cell_semantics(
                ac,
                std::span<const int>(processed_ids.data(), processed_ids.size()),
                std::span<const LevelSetCell<T, I>>(processed_cell, std::size_t(1)),
                /*total_num_level_sets=*/1,
                T(1e-12));
        }

        // Finalize derived topology/semantic layers used by selection + output.
        build_edges(ac);
        if (ac.tdim == 3)
            build_faces(ac);
        rebuild_zero_entity_inventory(ac);
        hc.adapt_cells.push_back(std::move(ac));

        // Single LS: bit 0 is always set.
        hc.active_level_set_mask.push_back(std::uint64_t(1));

        hc.parent_cell_ids.push_back(ci);
    }

    return {std::move(hc), std::move(parent_cells)};
}

template <std::floating_point T, std::integral I>
std::pair<HOCutCells<T, I>, ParentCellClassification<T, I>>
cut(const MeshView<T, I>& mesh,
    const LevelSetFunction<T, I>& ls,
    bool triangulate_cut_parts)
{
    CutOptions options;
    options.triangulate_cut_parts = triangulate_cut_parts;
    options.triangulation_strategy = cell::TriangulationStrategy::classical;
    options.cut_approximation = "auto";
    options.cut_approximation_order = 1;
    return cut(mesh, ls, options);
}

// =====================================================================
// cut() — multiple level sets
// =====================================================================

template <std::floating_point T, std::integral I>
std::pair<HOCutCells<T, I>, ParentCellClassification<T, I>>
cut(const MeshView<T, I>& mesh,
    const std::vector<LevelSetFunction<T, I>>& level_sets,
    const CutOptions& options)
{
    const CutOptions resolved_options = resolve_cut_options(options, level_sets);
    if (!mesh.has_cell_types())
        throw std::runtime_error("cut: MeshView must have cell types");
    if (level_sets.empty())
        throw std::runtime_error("cut: no level sets provided");
    if (level_sets.size() > 64)
        throw std::runtime_error("cut: more than 64 level sets not supported");

    const I ncells = mesh.num_cells();
    const int nls = static_cast<int>(level_sets.size());

    // --- ParentCellClassification ---
    ParentCellClassification<T, I> parent_cells;
    parent_cells.num_cells = static_cast<int>(ncells);
    parent_cells.num_level_sets = nls;
    parent_cells.cell_domains.assign(
        static_cast<std::size_t>(nls) * static_cast<std::size_t>(ncells),
        cell::domain::unset);
    parent_cells.cell_to_cut_index.assign(static_cast<std::size_t>(ncells), -1);
    for (const auto& ls : level_sets)
        parent_cells.level_set_names.push_back(ls.name);

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
    std::vector<I> mesh_cell_scratch;
    std::vector<I> ls_dof_scratch;

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
                ls_vertex_vals, mesh_cell_scratch, ls_dof_scratch, &ls_cell);

            parent_cells.cell_domains[static_cast<std::size_t>(
                li * static_cast<int>(ncells) + static_cast<int>(ci))] = dom;

            if (dom == cell::domain::intersected)
            {
                any_intersected = true;
                ls_cell.level_set_id = li;
                intersected_ls_indices.push_back(li);
                intersected_ls_cells.push_back(std::move(ls_cell));
            }
        }

        if (!any_intersected)
            continue;

        std::vector<int> all_level_set_indices(static_cast<std::size_t>(nls));
        std::iota(all_level_set_indices.begin(), all_level_set_indices.end(), 0);
        std::vector<LevelSetCell<T, I>> all_level_set_cells;
        all_level_set_cells.reserve(static_cast<std::size_t>(nls));
        for (int li = 0; li < nls; ++li)
        {
            auto ls_cell = make_cell_level_set(level_sets[static_cast<std::size_t>(li)], ci);
            ls_cell.level_set_id = li;
            all_level_set_cells.push_back(std::move(ls_cell));
        }

        const int cut_idx = hc.num_cut_cells();
        parent_cells.cell_to_cut_index[static_cast<std::size_t>(ci)] = cut_idx;

        // Build AdaptCell once per cell.
        AdaptCell<T> ac = make_adapt_cell(mesh, ci);
        apply_cut_approximation(ac, resolved_options);

        // Process intersecting level sets recursively (input order).
        //
        std::uint64_t cell_active_mask = 0;
        for (std::size_t k = 0; k < intersected_ls_indices.size(); ++k)
        {
            const int li = intersected_ls_indices[k];
            certify_refine_and_process_ready_cells(
                    ac, intersected_ls_cells[k], li,
                    /*max_iterations=*/8, T(1e-12), T(1e-12),
                    /*edge_max_depth=*/20,
                    resolved_options.triangulate_cut_parts
                        ? resolved_options.triangulation_strategy
                        : cell::TriangulationStrategy::none);

            // New vertices created while processing level set li must be
            // reclassified for all already-processed level sets.
            refresh_adapt_cell_semantics(
                ac,
                std::span<const int>(intersected_ls_indices.data(), k + 1),
                std::span<const LevelSetCell<T, I>>(intersected_ls_cells.data(), k + 1),
                nls,
                T(1e-12));

            cell_active_mask |= std::uint64_t(1) << li;
        }

        // Finalize derived topology/semantic layers used by selection + output.
        build_edges(ac);
        if (ac.tdim == 3)
            build_faces(ac);
        recompute_active_level_set_masks(ac, nls);
        refresh_adapt_cell_semantics(
            ac,
            std::span<const int>(all_level_set_indices.data(), all_level_set_indices.size()),
            std::span<const LevelSetCell<T, I>>(all_level_set_cells.data(), all_level_set_cells.size()),
            nls,
            T(1e-12));
        rebuild_zero_entity_inventory(ac);

        // Persist only level sets actively changing sign in this parent cell.
        // Non-active level sets may still touch a vertex/face, but they do not
        // change sign in this cut-cell.
        for (auto& ls_cell : intersected_ls_cells)
            hc.level_set_cells.push_back(std::move(ls_cell));
        hc.ls_offsets.push_back(
            static_cast<int>(hc.level_set_cells.size()));
        hc.active_level_set_mask.push_back(cell_active_mask);

        hc.adapt_cells.push_back(std::move(ac));
        hc.parent_cell_ids.push_back(ci);
    }

    return {std::move(hc), std::move(parent_cells)};
}

template <std::floating_point T, std::integral I>
std::pair<HOCutCells<T, I>, ParentCellClassification<T, I>>
cut(const MeshView<T, I>& mesh,
    const std::vector<LevelSetFunction<T, I>>& level_sets,
    bool triangulate_cut_parts)
{
    CutOptions options;
    options.triangulate_cut_parts = triangulate_cut_parts;
    options.triangulation_strategy = cell::TriangulationStrategy::classical;
    options.cut_approximation = "auto";
    options.cut_approximation_order = 1;
    return cut(mesh, level_sets, options);
}

// =====================================================================
// select_part()
// =====================================================================

template <std::floating_point T, std::integral I>
HOMeshPart<T, I> select_part(const MeshView<T, I>& mesh,
                              const HOCutCells<T, I>& cut_cells,
                              const ParentCellClassification<T, I>& parent_cells,
                              std::string_view expr_str)
{
    HOMeshPart<T, I> part;
    part.mesh = &mesh;
    part.cut_cells = &cut_cells;
    part.parent_cells        = &parent_cells;

    // Parse and compile the expression.
    part.expr = parse_selection_expr(expr_str);
    compile_selection_expr(part.expr, parent_cells.level_set_names);
    part.dim = infer_selection_dim(part.expr, cut_cells.tdim);

    const bool is_volume = (part.dim == cut_cells.tdim);

    // --- Uncut cells: scan cell_domains ---
    const int ncells = parent_cells.num_cells;
    for (I ci = 0; ci < static_cast<I>(ncells); ++ci)
    {
        if (parent_cells.cell_to_cut_index[static_cast<std::size_t>(ci)] >= 0)
            continue; // cut cell, handled below

        // For an uncut cell, every vertex has the same sign for each LS.
        // Check if the cell domain is compatible with the expression.
        bool match = false;
        for (const auto& term : part.expr.terms)
        {
            bool term_match = true;
            for (const auto& clause : term.clauses)
            {
                const int li = clause.level_set_index;
                const cell::domain dom = parent_cells.domain(li, ci);

                switch (clause.relation)
                {
                    case Relation::LessThan:
                    case Relation::LessEqual:
                        if (dom != cell::domain::inside) term_match = false;
                        break;
                    case Relation::GreaterThan:
                    case Relation::GreaterEqual:
                        if (dom != cell::domain::outside) term_match = false;
                        break;
                    case Relation::EqualTo:
                        // An uncut (non-intersected) cell has no zero interface.
                        term_match = false;
                        break;
                }
                if (!term_match) break;
            }
            if (term_match)
            {
                match = true;
                break;
            }
        }

        if (match)
            part.uncut_cell_ids.push_back(ci);
    }

    // --- Cut cells ---
    const int ncut = cut_cells.num_cut_cells();
    for (int k = 0; k < ncut; ++k)
    {
        const auto& ac = cut_cells.adapt_cells[static_cast<std::size_t>(k)];
        const I bg_cell = cut_cells.parent_cell_ids[static_cast<std::size_t>(k)];
        const std::uint64_t cut_active_mask =
            cut_cells.active_level_set_mask[static_cast<std::size_t>(k)];

        if (is_volume)
        {
            bool has_matching_leaf = false;
            const int tdim = ac.tdim;
            const int n_leaf = ac.n_entities(tdim);
            for (int c = 0; c < n_leaf; ++c)
            {
                auto cell_verts = ac.entity_to_vertex[tdim][static_cast<std::int32_t>(c)];
                if (leaf_cell_matches_selection_expr(
                        ac, cell_verts, part.expr, cut_active_mask, parent_cells, bg_cell))
                {
                    has_matching_leaf = true;
                    break;
                }
            }

            if (has_matching_leaf)
                part.cut_cell_ids.push_back(static_cast<std::int32_t>(k));
        }
        else
        {
            bool has_matching_zero_entity = false;
            const int n_zero = ac.n_zero_entities();
            for (int z = 0; z < n_zero; ++z)
            {
                if (zero_entity_matches_selection_expr(
                        ac, z, part.dim, part.expr, cut_active_mask, parent_cells, bg_cell))
                {
                    has_matching_zero_entity = true;
                    break;
                }
            }

            if (has_matching_zero_entity)
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
template std::pair<HOCutCells<double, int>, ParentCellClassification<double, int>>
cut(const MeshView<double, int>&, const LevelSetFunction<double, int>&,
    const CutOptions&);

template std::pair<HOCutCells<float, int>, ParentCellClassification<float, int>>
cut(const MeshView<float, int>&, const LevelSetFunction<float, int>&,
    const CutOptions&);

template std::pair<HOCutCells<double, long>, ParentCellClassification<double, long>>
cut(const MeshView<double, long>&, const LevelSetFunction<double, long>&,
    const CutOptions&);

template std::pair<HOCutCells<float, long>, ParentCellClassification<float, long>>
cut(const MeshView<float, long>&, const LevelSetFunction<float, long>&,
    const CutOptions&);

template std::pair<HOCutCells<double, int>, ParentCellClassification<double, int>>
cut(const MeshView<double, int>&, const LevelSetFunction<double, int>&, bool);

template std::pair<HOCutCells<float, int>, ParentCellClassification<float, int>>
cut(const MeshView<float, int>&, const LevelSetFunction<float, int>&, bool);

template std::pair<HOCutCells<double, long>, ParentCellClassification<double, long>>
cut(const MeshView<double, long>&, const LevelSetFunction<double, long>&, bool);

template std::pair<HOCutCells<float, long>, ParentCellClassification<float, long>>
cut(const MeshView<float, long>&, const LevelSetFunction<float, long>&, bool);

// cut() multi LS
template std::pair<HOCutCells<double, int>, ParentCellClassification<double, int>>
cut(const MeshView<double, int>&, const std::vector<LevelSetFunction<double, int>>&,
    const CutOptions&);

template std::pair<HOCutCells<float, int>, ParentCellClassification<float, int>>
cut(const MeshView<float, int>&, const std::vector<LevelSetFunction<float, int>>&,
    const CutOptions&);

template std::pair<HOCutCells<double, long>, ParentCellClassification<double, long>>
cut(const MeshView<double, long>&, const std::vector<LevelSetFunction<double, long>>&,
    const CutOptions&);

template std::pair<HOCutCells<float, long>, ParentCellClassification<float, long>>
cut(const MeshView<float, long>&, const std::vector<LevelSetFunction<float, long>>&,
    const CutOptions&);

template std::pair<HOCutCells<double, int>, ParentCellClassification<double, int>>
cut(const MeshView<double, int>&, const std::vector<LevelSetFunction<double, int>>&,
    bool);

template std::pair<HOCutCells<float, int>, ParentCellClassification<float, int>>
cut(const MeshView<float, int>&, const std::vector<LevelSetFunction<float, int>>&,
    bool);

template std::pair<HOCutCells<double, long>, ParentCellClassification<double, long>>
cut(const MeshView<double, long>&, const std::vector<LevelSetFunction<double, long>>&,
    bool);

template std::pair<HOCutCells<float, long>, ParentCellClassification<float, long>>
cut(const MeshView<float, long>&, const std::vector<LevelSetFunction<float, long>>&,
    bool);

// select_part()
template HOMeshPart<double, int>
select_part(const MeshView<double, int>&, const HOCutCells<double, int>&,
            const ParentCellClassification<double, int>&, std::string_view);

template HOMeshPart<float, int>
select_part(const MeshView<float, int>&, const HOCutCells<float, int>&,
            const ParentCellClassification<float, int>&, std::string_view);

template HOMeshPart<double, long>
select_part(const MeshView<double, long>&, const HOCutCells<double, long>&,
            const ParentCellClassification<double, long>&, std::string_view);

template HOMeshPart<float, long>
select_part(const MeshView<float, long>&, const HOCutCells<float, long>&,
            const ParentCellClassification<float, long>&, std::string_view);

} // namespace cutcells
