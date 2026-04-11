// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "level_set_cell.h"

#include <cassert>
#include <stdexcept>

#include "bernstein.h"
#include "mapping.h"

namespace cutcells
{
namespace
{

/// Determine the cell::type for a cell in the level-set mesh data.
/// Uses cell_types array if available, otherwise infers from tdim and
/// number of DOFs per cell (assuming equispaced Lagrange).
template <std::floating_point T, std::integral I>
cell::type infer_ls_cell_type(const LevelSetMeshData<T, I>& md, I cell_id)
{
    if (!md.cell_types.empty())
    {
        return md.cell_types[static_cast<std::size_t>(cell_id)];
    }

    // Infer from tdim and number of DOFs per cell
    const I ndofs = md.cell_num_dofs(cell_id);
    const int n = md.degree;

    switch (md.tdim)
    {
        case 1:
            return cell::type::interval;
        case 2:
        {
            // triangle: C(n+2,2) = (n+1)(n+2)/2
            // quad: (n+1)^2
            int tri_dofs = (n + 1) * (n + 2) / 2;
            if (ndofs == static_cast<I>(tri_dofs))
                return cell::type::triangle;
            else
                return cell::type::quadrilateral;
        }
        case 3:
        {
            // tet: C(n+3,3) = (n+1)(n+2)(n+3)/6
            // hex: (n+1)^3
            int tet_dofs = (n + 1) * (n + 2) * (n + 3) / 6;
            int hex_dofs = (n + 1) * (n + 1) * (n + 1);
            if (ndofs == static_cast<I>(tet_dofs))
                return cell::type::tetrahedron;
            else if (ndofs == static_cast<I>(hex_dofs))
                return cell::type::hexahedron;
            else
            {
                throw std::runtime_error(
                    "make_cell_level_set: cannot infer cell type from "
                    + std::to_string(ndofs) + " DOFs in 3D");
            }
        }
        default:
            throw std::runtime_error(
                "make_cell_level_set: unsupported tdim = "
                + std::to_string(md.tdim));
    }
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// make_cell_level_set
// ---------------------------------------------------------------------------

template <std::floating_point T, std::integral I>
LevelSetCell<T, I>
make_cell_level_set(const LevelSetFunction<T, I>& global_ls,
                    I cell_id)
{
    LevelSetCell<T, I> cell_ls;
    cell_ls.global_level_set = &global_ls;
    cell_ls.cell_id = cell_id;

    // --- Polynomial path: extract nodal values, convert to Bernstein ---
    if (global_ls.type == LevelSetType::Polynomial
        && global_ls.has_mesh_data()
        && global_ls.has_dof_values())
    {
        const auto& md = *global_ls.mesh_data;
        const int degree = md.degree;
        const int gdim = md.gdim;

        // Determine cell type
        cell::type ctype = infer_ls_cell_type(md, cell_id);
        cell_ls.cell_type = ctype;
        cell_ls.tdim = cell::get_tdim(ctype);

        // Get the cell's DOF indices
        auto cell_dofs = md.cell_dofs_span(cell_id);
        const int ndofs = static_cast<int>(cell_dofs.size());

        // Extract nodal values from the global DOF array
        cell_ls.nodal_values.resize(static_cast<std::size_t>(ndofs));
        for (int i = 0; i < ndofs; ++i)
            cell_ls.nodal_values[static_cast<std::size_t>(i)]
                = global_ls.dof_values[static_cast<std::size_t>(cell_dofs[static_cast<std::size_t>(i)])];
        cell_ls.nodal_order = degree;

        // Extract physical coordinates of all DOFs on this cell
        std::vector<T> dof_phys(static_cast<std::size_t>(ndofs * gdim));
        for (int i = 0; i < ndofs; ++i)
        {
            const T* x = md.dof_coordinate(cell_dofs[static_cast<std::size_t>(i)]);
            for (int d = 0; d < gdim; ++d)
                dof_phys[static_cast<std::size_t>(i * gdim + d)] = x[d];
        }

        // Extract physical vertex coordinates (first n_vertices DOFs in Basix ordering)
        const int nv = cell::get_num_vertices(ctype);
        std::vector<T> vertex_coords(static_cast<std::size_t>(nv * gdim));
        for (int v = 0; v < nv; ++v)
            for (int d = 0; d < gdim; ++d)
                vertex_coords[static_cast<std::size_t>(v * gdim + d)]
                    = dof_phys[static_cast<std::size_t>(v * gdim + d)];

        // Pull back all DOF coordinates from physical to reference space.
        // Note: pull_back_affine is supported for volume cells (gdim == tdim).
        assert(gdim == md.tdim && "pull_back_affine requires gdim == tdim");
        std::vector<T> dof_ref(static_cast<std::size_t>(ndofs * gdim));
        cell::pull_back_affine(ctype, vertex_coords, gdim,
                               std::span<const T>(dof_phys),
                               std::span<T>(dof_ref));

        // Convert Lagrange nodal values to Bernstein coefficients
        cell_ls.bernstein_order = degree;
        bernstein::lagrange_to_bernstein(
            ctype, degree,
            std::span<const T>(dof_ref),
            std::span<const T>(cell_ls.nodal_values),
            cell_ls.bernstein_coeffs);
    }

    return cell_ls;
}

// ---------------------------------------------------------------------------
// LevelSetCell::value
// ---------------------------------------------------------------------------

template <std::floating_point T, std::integral I>
T LevelSetCell<T, I>::value(std::span<const T> xi) const
{
    // Polynomial backend: evaluate the Bernstein expansion
    if (!bernstein_coeffs.empty())
    {
        return bernstein::evaluate(
            cell_type, bernstein_order,
            std::span<const T>(bernstein_coeffs), xi);
    }

    // Fallback: delegate to the global level set value function (in physical
    // coordinates). This path requires an external mapping from reference to
    // physical space, which is not yet wired up here.
    throw std::runtime_error(
        "LevelSetCell::value: no Bernstein coefficients available "
        "and analytical fallback not yet implemented");
}

// ---------------------------------------------------------------------------
// LevelSetCell::grad
// ---------------------------------------------------------------------------

template <std::floating_point T, std::integral I>
void LevelSetCell<T, I>::grad(std::span<const T> xi, std::span<T> g) const
{
    // Polynomial backend: gradient of the Bernstein expansion
    if (!bernstein_coeffs.empty())
    {
        bernstein::gradient(
            cell_type, bernstein_order,
            std::span<const T>(bernstein_coeffs), xi, g);
        return;
    }

    throw std::runtime_error(
        "LevelSetCell::grad: no Bernstein coefficients available "
        "and analytical fallback not yet implemented");
}

// ---------------------------------------------------------------------------
// Explicit template instantiations
// ---------------------------------------------------------------------------

template LevelSetCell<double, int>
make_cell_level_set(const LevelSetFunction<double, int>&, int);

template LevelSetCell<float, int>
make_cell_level_set(const LevelSetFunction<float, int>&, int);

template LevelSetCell<double, long>
make_cell_level_set(const LevelSetFunction<double, long>&, long);

template LevelSetCell<float, long>
make_cell_level_set(const LevelSetFunction<float, long>&, long);

template double LevelSetCell<double, int>::value(std::span<const double>) const;
template float  LevelSetCell<float,  int>::value(std::span<const float>)  const;
template double LevelSetCell<double, long>::value(std::span<const double>) const;
template float  LevelSetCell<float,  long>::value(std::span<const float>)  const;

template void LevelSetCell<double, int>::grad(std::span<const double>, std::span<double>) const;
template void LevelSetCell<float,  int>::grad(std::span<const float>,  std::span<float>)  const;
template void LevelSetCell<double, long>::grad(std::span<const double>, std::span<double>) const;
template void LevelSetCell<float,  long>::grad(std::span<const float>,  std::span<float>)  const;

} // namespace cutcells
