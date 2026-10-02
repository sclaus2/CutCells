// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "cut_result.h"

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <tuple>

#include "../cell_topology.h"
#include "../reference_cell.h"
#include "cell_source.h"
#include "classify.h"

namespace cutcells::part
{

namespace
{

/// Domain from the signs of a cell's Bernstein coefficients, for cells that
/// quadrays does not take.
template <std::floating_point T, std::integral I>
cell::domain coefficient_domain(const LevelSetFunction<T, I>& ls, I cell_id, LevelSetCell<T, I>& scratch)
{
    if (ls.analytic)
    {
        throw std::invalid_argument("part::cut: analytic level sets need tetrahedra or hexahedra in 3D "
                                    "(other cells come in phase 6)");
    }
    if (!(ls.type == LevelSetType::Polynomial && ls.has_mesh_data() && ls.has_dof_values()))
    {
        throw std::invalid_argument("part::cut: the level set '" + ls.name
                                    + "' has neither dof values nor an analytic level set");
    }
    make_cell_level_set(ls, cell_id, scratch);
    const int s = quadrays::coefficient_sign(std::span<const T>(scratch.bernstein_coeffs));
    return s > 0 ? cell::domain::outside : (s < 0 ? cell::domain::inside : cell::domain::intersected);
}

/// A face that lies in a zero set, before its owner is chosen.
template <std::integral I>
struct FaceEntry
{
    int level_set = 0;
    std::array<I, 4> key{}; ///< mesh nodes of the face, sorted, padded
    I cell = 0;
    int local = 0;
    int side = 0; ///< -1: phi < 0 next to the face in the cell, +1: > 0, 0: unknown
};

/// The faces of a cell on which level set l vanishes, appended to @p entries:
/// the level set is 0 up to the engine's tolerance at the face's vertices and,
/// by its bounds, on the whole face.
template <std::floating_point T, std::integral I>
void find_zero_faces(const MeshView<T, I>& mesh, CellSource<T, I>& cs, I cell_id, int l, cell::domain dom,
                     std::vector<FaceEntry<I>>& entries)
{
    const cell::type type = mesh.cell_type(cell_id);
    const T scale = quadrays::reference_magnitude(cs.source);
    if (!(scale > T(0)))
        return;
    const T tol = quadrays::scaled_tolerance<T>(1e-12) * scale;
    const int nv = cell::get_num_vertices(type);
    std::array<T, 8> values{};
    for (int v = 0; v < nv; ++v)
    {
        const quadrays::Vec3<T> u = box_vertex<T>(type, v);
        values[v] = quadrays::evaluate(cs.source, std::span<const T>(u));
    }
    thread_local quadrays::BoxBernstein<T> face_form;
    thread_local std::vector<T> work;
    const std::span<const I> nodes = mesh.cell_nodes(cell_id, cs.nodes);
    for (int f = 0; f < cell::num_faces(type); ++f)
    {
        const std::span<const int> fv = cell::face_vertices(type, f);
        bool zero = true;
        for (const int v : fv)
            zero &= std::abs(values[v]) <= tol;
        if (!zero)
            continue;
        // the face as the image of [0, 1]^2: u_a + s (u_b - u_a) + t (u_c - u_a)
        const quadrays::Vec3<T> ua = box_vertex<T>(type, fv[0]), ub = box_vertex<T>(type, fv[1]),
                                uc = box_vertex<T>(type, fv[2]);
        std::array<T, 6> matrix;
        for (int i = 0; i < 3; ++i)
        {
            matrix[2 * i] = ub[i] - ua[i];
            matrix[2 * i + 1] = uc[i] - ua[i];
        }
        if (cs.source.bernstein != nullptr)
        {
            quadrays::restrict_affine(*cs.source.bernstein, std::span<const T>(ua), std::span<const T>(matrix), 2,
                                      face_form, work);
            zero = quadrays::max_abs(std::span<const T>(face_form.coeffs)) <= tol;
        }
        else
        {
            quadrays::AffineBounds<T> b;
            zero = quadrays::affine_bounds(cs.source, std::span<const T>(ua), std::span<const T>(matrix), 2, b)
                   && b.magnitude <= tol;
        }
        if (!zero)
            continue;

        FaceEntry<I> e;
        e.level_set = l;
        e.cell = cell_id;
        e.local = f;
        if (dom == cell::domain::inside)
            e.side = -1;
        else if (dom == cell::domain::outside)
            e.side = 1;
        else
        {
            // a cut cell: the derivative into the cell at the face's centroid
            quadrays::Vec3<T> centroid = {0, 0, 0}, inner = {0, 0, 0};
            for (const int v : fv)
            {
                const quadrays::Vec3<T> u = box_vertex<T>(type, v);
                for (int i = 0; i < 3; ++i)
                    centroid[i] += u[i] / T(fv.size());
            }
            for (int v = 0; v < nv; ++v)
            {
                const quadrays::Vec3<T> u = box_vertex<T>(type, v);
                for (int i = 0; i < 3; ++i)
                    inner[i] += u[i] / T(nv);
            }
            quadrays::Vec3<T> g;
            quadrays::gradient(cs.source, std::span<const T>(centroid), std::span<T>(g));
            T d = T(0);
            for (int i = 0; i < 3; ++i)
                d += g[i] * (inner[i] - centroid[i]);
            e.side = d < -tol ? -1 : (d > tol ? 1 : 0);
        }
        e.key.fill(std::numeric_limits<I>::max());
        for (std::size_t k = 0; k < fv.size(); ++k)
        {
            const int local_vertex = mesh.vtk_vertex_order ? cell::basix_to_vtk_vertex(type, fv[k]) : fv[k];
            e.key[k] = nodes[static_cast<std::size_t>(local_vertex)];
        }
        std::sort(e.key.begin(), e.key.end());
        entries.push_back(e);
    }
}

} // namespace

template <std::floating_point T, std::integral I>
CutResult<T, I> cut(const MeshView<T, I>& mesh, std::span<const LevelSetFunction<T, I>> level_sets,
                    const ClassifyOptions& options)
{
    if (!mesh.has_cell_types())
        throw std::invalid_argument("part::cut: the mesh needs cell types");
    if (level_sets.empty() || level_sets.size() > 64)
        throw std::invalid_argument("part::cut: give 1 to 64 level sets");

    CutResult<T, I> r;
    r.mesh = &mesh;
    r.num_cells = static_cast<int>(mesh.num_cells());
    for (const LevelSetFunction<T, I>& ls : level_sets)
    {
        r.level_sets.push_back(&ls);
        r.level_set_names.push_back(ls.name);
    }
    const int nls = r.n_level_sets();
    r.domains.assign(static_cast<std::size_t>(nls) * static_cast<std::size_t>(r.num_cells), cell::domain::unset);

    std::vector<FaceEntry<I>> entries;
    CellSource<T, I> cs;
    LevelSetCell<T, I> scratch;
    for (I c = 0; c < static_cast<I>(r.num_cells); ++c)
    {
        bool any_cut = false;
        for (int l = 0; l < nls; ++l)
        {
            const LevelSetFunction<T, I>& ls = level_sets[static_cast<std::size_t>(l)];
            cell::domain dom;
            if (cell_source(mesh, ls, c, cs))
            {
                const int s = cell_sign(cs, mesh.cell_type(c), options.max_depth);
                dom = s > 0 ? cell::domain::outside : (s < 0 ? cell::domain::inside : cell::domain::intersected);
                find_zero_faces(mesh, cs, c, l, dom, entries);
            }
            else
                dom = coefficient_domain(ls, c, scratch);
            r.domains[static_cast<std::size_t>(l) * static_cast<std::size_t>(r.num_cells)
                      + static_cast<std::size_t>(c)]
                = dom;
            any_cut |= dom == cell::domain::intersected;
        }
        if (any_cut)
            r.cut_cells.push_back(c);
    }

    // one owner per face: the cell on its negative side, else the lower index
    std::sort(entries.begin(), entries.end(),
              [](const FaceEntry<I>& a, const FaceEntry<I>& b)
              { return std::tie(a.level_set, a.key, a.cell) < std::tie(b.level_set, b.key, b.cell); });
    for (std::size_t i = 0; i < entries.size();)
    {
        std::size_t j = i;
        while (j < entries.size() && entries[j].level_set == entries[i].level_set && entries[j].key == entries[i].key)
            ++j;
        std::size_t owner = i;
        int negative = 0;
        for (std::size_t k = i; k < j; ++k)
            if (entries[k].side < 0)
            {
                ++negative;
                owner = k;
            }
        if (negative != 1)
            owner = i;
        r.zero_face_level_sets.push_back(entries[owner].level_set);
        r.zero_face_cells.push_back(entries[owner].cell);
        r.zero_face_local.push_back(static_cast<std::int8_t>(entries[owner].local));
        i = j;
    }
    return r;
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template CutResult<float, int> cut<float, int>(const MeshView<float, int>&,
                                               std::span<const LevelSetFunction<float, int>>,
                                               const ClassifyOptions&);
template CutResult<double, int> cut<double, int>(const MeshView<double, int>&,
                                                 std::span<const LevelSetFunction<double, int>>,
                                                 const ClassifyOptions&);

} // namespace cutcells::part
