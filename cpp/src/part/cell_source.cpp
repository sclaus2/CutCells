// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "cell_source.h"

#include <span>
#include <stdexcept>

namespace cutcells::part
{

namespace
{
/// The level set on the cell's box: the analytic level set, or the Bernstein
/// form of its cell polynomial, written to @p form.
template <std::floating_point T, std::integral I>
quadrays::Source<T> source_on_box(const LevelSetFunction<T, I>& ls, I cell_id, cell::type type,
                                  const quadrays::ClippedBox<T>& box, quadrays::BoxBernstein<T>& form,
                                  LevelSetCell<T, I>& ls_cell)
{
    if (ls.analytic)
        return quadrays::analytic_source(*ls.analytic, box);
    if (ls.type == LevelSetType::Polynomial && ls.has_mesh_data() && ls.has_dof_values())
    {
        make_cell_level_set(ls, cell_id, ls_cell);
        quadrays::cell_bernstein_on_box<T>(type, ls_cell.bernstein_order, std::span<const T>(ls_cell.bernstein_coeffs),
                                           box, form);
        return quadrays::bernstein_source(form);
    }
    throw std::invalid_argument("part: the level set '" + ls.name
                                + "' has neither dof values nor an analytic level set");
}
} // namespace

template <std::floating_point T, std::integral I>
bool cell_source(const MeshView<T, I>& mesh, const LevelSetFunction<T, I>& ls, I cell_id,
                 CellSource<T, I>& out, quadrays::BoxFrame frame)
{
    const cell::type type = mesh.cell_type(cell_id);
    if (!quadrays_takes(type, mesh.gdim, ls))
        return false;
    cell_vertex_coords_basix(mesh, cell_id, out.vertices, out.nodes);
    quadrays::make_clipped_box<T>(type, std::span<const T>(out.vertices), mesh.gdim, out.box, frame);
    out.source = source_on_box(ls, cell_id, type, out.box, out.form, out.ls_cell);
    return true;
}

template <std::floating_point T, std::integral I>
bool cell_sources(const MeshView<T, I>& mesh, std::span<const LevelSetFunction<T, I>* const> level_sets,
                  I cell_id, CellSources<T, I>& out)
{
    const cell::type type = mesh.cell_type(cell_id);
    for (const LevelSetFunction<T, I>* ls : level_sets)
        if (!quadrays_takes(type, mesh.gdim, *ls))
            return false;
    cell_vertex_coords_basix(mesh, cell_id, out.vertices, out.nodes);
    quadrays::make_clipped_box<T>(type, std::span<const T>(out.vertices), mesh.gdim, out.box);
    // the sources point into forms: size it first
    if (out.forms.size() < level_sets.size())
        out.forms.resize(level_sets.size());
    out.sources.clear();
    for (std::size_t l = 0; l < level_sets.size(); ++l)
        out.sources.push_back(source_on_box(*level_sets[l], cell_id, type, out.box, out.forms[l], out.ls_cell));
    return true;
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template bool cell_source<float, int>(const MeshView<float, int>&, const LevelSetFunction<float, int>&, int,
                                      CellSource<float, int>&, quadrays::BoxFrame);
template bool cell_source<double, int>(const MeshView<double, int>&, const LevelSetFunction<double, int>&, int,
                                       CellSource<double, int>&, quadrays::BoxFrame);
template bool cell_sources<float, int>(const MeshView<float, int>&,
                                       std::span<const LevelSetFunction<float, int>* const>, int,
                                       CellSources<float, int>&);
template bool cell_sources<double, int>(const MeshView<double, int>&,
                                        std::span<const LevelSetFunction<double, int>* const>, int,
                                        CellSources<double, int>&);

} // namespace cutcells::part
