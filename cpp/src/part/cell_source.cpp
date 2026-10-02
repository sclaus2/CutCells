// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "cell_source.h"

#include <span>
#include <stdexcept>

namespace cutcells::part
{

template <std::floating_point T, std::integral I>
bool cell_source(const MeshView<T, I>& mesh, const LevelSetFunction<T, I>& ls, I cell_id,
                 CellSource<T, I>& out)
{
    const cell::type type = mesh.cell_type(cell_id);
    if (mesh.gdim != 3 || (type != cell::type::tetrahedron && type != cell::type::hexahedron))
        return false;
    cell_vertex_coords_basix(mesh, cell_id, out.vertices, out.nodes);
    quadrays::make_clipped_box<T>(type, std::span<const T>(out.vertices), 3, out.box);
    if (ls.analytic)
    {
        out.source = quadrays::analytic_source(*ls.analytic, out.box);
        return true;
    }
    if (ls.type == LevelSetType::Polynomial && ls.has_mesh_data() && ls.has_dof_values())
    {
        make_cell_level_set(ls, cell_id, out.ls_cell);
        quadrays::cell_bernstein_on_box<T>(type, out.ls_cell.bernstein_order,
                                           std::span<const T>(out.ls_cell.bernstein_coeffs), out.form);
        out.source = quadrays::bernstein_source(out.form);
        return true;
    }
    throw std::invalid_argument("part: the level set '" + ls.name
                                + "' has neither dof values nor an analytic level set");
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template bool cell_source<float, int>(const MeshView<float, int>&, const LevelSetFunction<float, int>&, int,
                                      CellSource<float, int>&);
template bool cell_source<double, int>(const MeshView<double, int>&, const LevelSetFunction<double, int>&, int,
                                       CellSource<double, int>&);

} // namespace cutcells::part
