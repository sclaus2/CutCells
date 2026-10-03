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

#include "../bernstein.h"
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
        throw std::invalid_argument("part::cut: analytic level sets need triangles or quadrilaterals in 2D, or "
                                    "tetrahedra, hexahedra, prisms or pyramids in 3D");
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

/// The facets of a cell (faces in 3D, edges in 2D) on which level set l
/// vanishes, appended to @p entries: the level set is 0 up to the engine's
/// tolerance at the facet's vertices and, by its bounds, on the whole facet.
template <std::floating_point T, std::integral I>
void find_zero_faces(const MeshView<T, I>& mesh, CellSource<T, I>& cs, I cell_id, int l, cell::domain dom,
                     std::vector<FaceEntry<I>>& entries)
{
    const cell::type type = mesh.cell_type(cell_id);
    const int tdim = cell::get_tdim(type);
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
    const int m = tdim - 1; // dimension of the facets
    for (int f = 0; f < num_facets(type); ++f)
    {
        const std::span<const int> fv = facet_vertices(type, f);
        bool zero = true;
        for (const int v : fv)
            zero &= std::abs(values[v]) <= tol;
        if (!zero)
            continue;
        // the facet as the image of [0, 1]^m: u_a + s (u_b - u_a) (+ t (u_c - u_a))
        const quadrays::Vec3<T> ua = box_vertex<T>(type, fv[0]), ub = box_vertex<T>(type, fv[1]),
                                uc = box_vertex<T>(type, fv[m > 1 ? 2 : 1]);
        std::array<T, 6> matrix{};
        for (int i = 0; i < 3; ++i)
        {
            matrix[i * m] = ub[i] - ua[i];
            if (m > 1)
                matrix[i * m + 1] = uc[i] - ua[i];
        }
        if (cs.source.bernstein != nullptr)
        {
            const std::size_t dim = static_cast<std::size_t>(cs.source.bernstein->dim);
            quadrays::restrict_affine(*cs.source.bernstein, std::span<const T>(ua.data(), dim),
                                      std::span<const T>(matrix.data(), dim * static_cast<std::size_t>(m)), m,
                                      face_form, work);
            zero = quadrays::max_abs(std::span<const T>(face_form.coeffs)) <= tol;
        }
        else
        {
            quadrays::AffineBounds<T> b;
            zero = quadrays::affine_bounds(cs.source, std::span<const T>(ua), std::span<const T>(matrix.data(), 3 * m),
                                           m, b)
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
            // a cut cell: the derivative into the cell at the facet's centroid
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
            quadrays::Vec3<T> g = {0, 0, 0};
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

/// The facets of a cell classified by its Bernstein coefficients (cells
/// quadrays does not take) on which level set l vanishes, appended to @p entries:
/// the level set is 0 up to the tolerance on the facet's lattice of its
/// degree, which fixes the polynomial on the facet.
template <std::floating_point T, std::integral I>
void find_zero_facets(const MeshView<T, I>& mesh, const LevelSetCell<T, I>& ls_cell, I cell_id, int l,
                      cell::domain dom, std::vector<I>& node_scratch, std::vector<FaceEntry<I>>& entries)
{
    const cell::type type = mesh.cell_type(cell_id);
    const int tdim = cell::get_tdim(type);
    const std::span<const T> c(ls_cell.bernstein_coeffs);
    T scale = T(0);
    for (const T x : c)
        scale = std::max(scale, std::abs(x));
    if (tdim < 2 || !(scale > T(0)))
        return;
    const T tol = quadrays::scaled_tolerance<T>(1e-12) * scale;
    const int degree = ls_cell.bernstein_order, p = std::max(degree, 1);
    const std::vector<T> ref = cell::reference_vertices<T>(type);
    auto value = [&](const T* x) { return bernstein::evaluate<T>(type, degree, c, std::span<const T>(x, tdim)); };
    const std::span<const I> nodes = mesh.cell_nodes(cell_id, node_scratch);
    std::array<T, 3> x{}, centroid{}, inner{}, g{};
    for (int f = 0; f < num_facets(type); ++f)
    {
        const std::span<const int> fv = facet_vertices(type, f);
        const cell::type ft = facet_type(type, f);
        // the lattice a + i/p (b - a) + j/p (c - a)
        const T* a = ref.data() + fv[0] * tdim;
        const T* b = ref.data() + fv[1] * tdim;
        const T* e = fv.size() > 2 ? ref.data() + fv[2] * tdim : a;
        const int nj = ft == cell::type::interval ? 0 : p;
        bool zero = true;
        for (int j = 0; zero && j <= nj; ++j)
            for (int i = 0; zero && i <= p; ++i)
            {
                if (ft == cell::type::triangle && i + j > p)
                    continue;
                for (int d = 0; d < tdim; ++d)
                    x[d] = a[d] + T(i) / T(p) * (b[d] - a[d]) + T(j) / T(p) * (e[d] - a[d]);
                zero = std::abs(value(x.data())) <= tol;
            }
        if (!zero)
            continue;

        FaceEntry<I> entry;
        entry.level_set = l;
        entry.cell = cell_id;
        entry.local = f;
        if (dom == cell::domain::inside)
            entry.side = -1;
        else if (dom == cell::domain::outside)
            entry.side = 1;
        else
        {
            // a cut cell: the derivative into the cell at the facet's centroid
            centroid.fill(T(0));
            inner.fill(T(0));
            const int nv = cell::get_num_vertices(type);
            for (const int v : fv)
                for (int d = 0; d < tdim; ++d)
                    centroid[d] += ref[static_cast<std::size_t>(v * tdim + d)] / T(fv.size());
            for (int v = 0; v < nv; ++v)
                for (int d = 0; d < tdim; ++d)
                    inner[d] += ref[static_cast<std::size_t>(v * tdim + d)] / T(nv);
            bernstein::gradient<T>(type, degree, c, std::span<const T>(centroid.data(), tdim),
                                   std::span<T>(g.data(), tdim));
            T dd = T(0);
            for (int d = 0; d < tdim; ++d)
                dd += g[d] * (inner[d] - centroid[d]);
            entry.side = dd < -tol ? -1 : (dd > tol ? 1 : 0);
        }
        entry.key.fill(std::numeric_limits<I>::max());
        for (std::size_t k = 0; k < fv.size(); ++k)
        {
            const int local_vertex = mesh.vtk_vertex_order ? cell::basix_to_vtk_vertex(type, fv[k]) : fv[k];
            entry.key[k] = nodes[static_cast<std::size_t>(local_vertex)];
        }
        std::sort(entry.key.begin(), entry.key.end());
        entries.push_back(entry);
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
    std::vector<I> nodes;
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
            {
                dom = coefficient_domain(ls, c, scratch);
                find_zero_facets(mesh, scratch, c, l, dom, nodes, entries);
            }
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
