// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "cell_pieces.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include "cut_cell.h"
#include "iso_refine.h"
#include "triangulation.h"
#include "piece_rules.h"

namespace cutcells::lut
{

namespace
{

bool is_simplex(cell::type type)
{
    return type == cell::type::interval || type == cell::type::triangle || type == cell::type::tetrahedron;
}

/// A template sub-cell and its interpolant: P1 on simplices, multilinear on
/// the axis-aligned quadrilaterals and hexahedra of tensor templates.
template <std::floating_point T>
struct SubCell
{
    cell::type type = cell::type::point;
    int tdim = 0;
    std::array<T, 3> origin{};
    std::array<T, 9> inverse{}; ///< simplices: inverse of the edge matrix, row-major
    std::array<T, 3> scale{};   ///< boxes: 1 / edge lengths
};

template <std::floating_point T>
void make_sub_cell(cell::type type, int tdim, const T* x, SubCell<T>& s)
{
    s.type = type;
    s.tdim = tdim;
    for (int d = 0; d < tdim; ++d)
        s.origin[d] = x[d];
    if (!is_simplex(type))
    {
        // vertex 0 and the last vertex are opposite corners of the box
        const int last = cell::get_num_vertices(type) - 1;
        for (int d = 0; d < tdim; ++d)
            s.scale[d] = T(1) / (x[last * tdim + d] - x[d]);
        return;
    }
    // edge matrix e[d][j] = x_{j+1}[d] - x_0[d] and its inverse
    std::array<T, 9> e{};
    for (int d = 0; d < tdim; ++d)
        for (int j = 0; j < tdim; ++j)
            e[d * tdim + j] = x[(j + 1) * tdim + d] - x[d];
    if (tdim == 1)
        s.inverse[0] = T(1) / e[0];
    else if (tdim == 2)
    {
        const T det = e[0] * e[3] - e[1] * e[2];
        s.inverse = {e[3] / det, -e[1] / det, -e[2] / det, e[0] / det};
    }
    else
    {
        const T det = e[0] * (e[4] * e[8] - e[5] * e[7]) - e[1] * (e[3] * e[8] - e[5] * e[6])
                      + e[2] * (e[3] * e[7] - e[4] * e[6]);
        s.inverse = {(e[4] * e[8] - e[5] * e[7]) / det, (e[2] * e[7] - e[1] * e[8]) / det,
                     (e[1] * e[5] - e[2] * e[4]) / det, (e[5] * e[6] - e[3] * e[8]) / det,
                     (e[0] * e[8] - e[2] * e[6]) / det, (e[2] * e[3] - e[0] * e[5]) / det,
                     (e[3] * e[7] - e[4] * e[6]) / det, (e[1] * e[6] - e[0] * e[7]) / det,
                     (e[0] * e[4] - e[1] * e[3]) / det};
    }
}

/// The sub-cell's interpolant of @p values (at its vertices) at point @p p.
template <std::floating_point T>
T interpolate(const SubCell<T>& s, const T* values, const T* p)
{
    const int tdim = s.tdim;
    if (is_simplex(s.type))
    {
        T sum = 0, result = 0;
        for (int j = 0; j < tdim; ++j)
        {
            T lambda = 0;
            for (int d = 0; d < tdim; ++d)
                lambda += s.inverse[j * tdim + d] * (p[d] - s.origin[d]);
            sum += lambda;
            result += values[j + 1] * lambda;
        }
        return result + values[0] * (T(1) - sum);
    }
    std::array<T, 3> xi{};
    for (int d = 0; d < tdim; ++d)
        xi[d] = (p[d] - s.origin[d]) * s.scale[d];
    T result = 0;
    for (int v = 0; v < (1 << tdim); ++v)
    {
        T w = values[v];
        for (int d = 0; d < tdim; ++d)
            w *= ((v >> d) & 1) ? xi[d] : T(1) - xi[d];
        result += w;
    }
    return result;
}

/// -1 if all values are negative, +1 if none is (0 counts as positive), 0 if both occur.
template <std::floating_point T>
int sign_class(std::span<const T> v)
{
    bool negative = false, positive = false;
    for (const T x : v)
    {
        negative |= x < T(0);
        positive |= x >= T(0);
    }
    return negative && positive ? 0 : (negative ? -1 : 1);
}

/// One straight cell of the arrangement, its sides and the zero sets it lies in.
template <std::floating_point T>
struct Leaf
{
    cell::type type = cell::type::point;
    std::vector<T> x; ///< reference coordinates, tdim per vertex
    std::uint64_t negative = 0, positive = 0, zero = 0;
};

/// The cells of a part that the lookup tables cut from a cell of type
/// @p cut_type, as leaves in Basix vertex order. The tables of quadrilaterals,
/// hexahedra, prisms and pyramids take and give Basix order; those of
/// triangles and tetrahedra give their quadrilaterals in cyclic (VTK) order.
template <std::floating_point T>
void append_part(const cell::CutCell<T>& part, cell::type cut_type, int tdim, std::uint64_t negative,
                 std::uint64_t positive, std::uint64_t zero, std::vector<Leaf<T>>& out)
{
    const bool cyclic_quadrilaterals = cut_type == cell::type::triangle || cut_type == cell::type::tetrahedron;
    for (int q = 0; q < cell::num_cells(part); ++q)
    {
        Leaf<T> leaf;
        leaf.type = part._types[static_cast<std::size_t>(q)];
        const std::span<const int> vertices = cell::cell_vertices(part, q);
        const bool swap = cyclic_quadrilaterals && leaf.type == cell::type::quadrilateral;
        for (std::size_t j = 0; j < vertices.size(); ++j)
        {
            const int v = vertices[swap && j >= 2 ? 5 - j : j];
            leaf.x.insert(leaf.x.end(), part._vertex_coords.begin() + v * tdim,
                          part._vertex_coords.begin() + (v + 1) * tdim);
        }
        leaf.negative = negative;
        leaf.positive = positive;
        leaf.zero = zero;
        out.push_back(std::move(leaf));
    }
}

/// Values the lookup tables read as we classify them: values within
/// 64 eps max|v| of 0 become that bound, so they count as positive here and in
/// both of the tables' case masks (which treat values within 2 eps |v_0| as 0),
/// and intersections next to them land on the vertex.
template <std::floating_point T>
void snap_values(std::vector<T>& v)
{
    T scale = T(0);
    for (const T x : v)
        scale = std::max(scale, std::abs(x));
    const T zero = T(64) * std::numeric_limits<T>::epsilon() * scale;
    for (T& x : v)
        if (std::abs(x) <= zero)
            x = zero;
}

/// Split a leaf by level set i (its interpolant on the sub-cell), both sides
/// kept and tagged; with @p add_zero also its zero set in the leaf. The tables
/// cut parallelograms and parallelepipeds exactly; other quadrilaterals and
/// hexahedra, which earlier cuts leave, are cut as simplices.
template <std::floating_point T>
void cut_leaf(Leaf<T>& leaf, const SubCell<T>& s, const T* sub_values, int i, bool add_zero,
              cell::TriangulationStrategy triangulation, std::vector<Leaf<T>>& out)
{
    const int dim = s.tdim;
    const std::uint64_t bit = std::uint64_t(1) << i;
    const int nv = cell::get_num_vertices(leaf.type);
    std::vector<T> v(static_cast<std::size_t>(nv));
    for (int j = 0; j < nv; ++j)
        v[static_cast<std::size_t>(j)] = interpolate(s, sub_values, leaf.x.data() + j * dim);
    snap_values(v);
    const int c = sign_class(std::span<const T>(v));
    if (c != 0)
    {
        (c < 0 ? leaf.negative : leaf.positive) |= bit;
        out.push_back(std::move(leaf));
        return;
    }
    if ((leaf.type == cell::type::quadrilateral || leaf.type == cell::type::hexahedron)
        && !is_parallelotope(std::span<const T>(leaf.x), cell::get_tdim(leaf.type), dim))
    {
        int ids[8] = {0, 1, 2, 3, 4, 5, 6, 7};
        std::vector<std::vector<int>> simplices;
        cell::triangulation(leaf.type, ids, simplices);
        const cell::type simplex
            = leaf.type == cell::type::hexahedron ? cell::type::tetrahedron : cell::type::triangle;
        for (const std::vector<int>& t : simplices)
        {
            Leaf<T> piece{simplex, {}, leaf.negative, leaf.positive, leaf.zero};
            for (const int j : t)
                piece.x.insert(piece.x.end(), leaf.x.begin() + j * dim, leaf.x.begin() + (j + 1) * dim);
            cut_leaf(piece, s, sub_values, i, add_zero, triangulation, out);
        }
        return;
    }
    for (const bool below : {true, false})
    {
        cell::CutCell<T> part;
        cell::cut<T>(leaf.type, std::span<const T>(leaf.x), dim, std::span<const T>(v), below ? "phi<0" : "phi>0",
                     part, triangulation);
        append_part(part, leaf.type, dim, leaf.negative | (below ? bit : 0), leaf.positive | (below ? 0 : bit),
                    leaf.zero, out);
    }
    if (add_zero)
    {
        cell::CutCell<T> part;
        cell::cut<T>(leaf.type, std::span<const T>(leaf.x), dim, std::span<const T>(v), "phi=0", part, triangulation);
        append_part(part, leaf.type, dim, leaf.negative, leaf.positive, leaf.zero | bit, out);
    }
}

/// Split every leaf by level set i; the leaves in one zero set of
/// @p curve_sets also give their part in i's zero set if i is in it too.
template <std::floating_point T>
void cut_leaves(std::vector<Leaf<T>>& leaves, const SubCell<T>& s, const T* sub_values, int i,
                std::uint64_t curve_sets, cell::TriangulationStrategy triangulation)
{
    std::vector<Leaf<T>> out;
    for (Leaf<T>& leaf : leaves)
    {
        const bool add_zero = ((curve_sets >> i) & 1) && std::popcount(leaf.zero) == 1
                              && i > std::countr_zero(leaf.zero);
        cut_leaf(leaf, s, sub_values, i, add_zero, triangulation, out);
    }
    leaves = std::move(out);
}

template <std::floating_point T>
void append_leaves(const std::vector<Leaf<T>>& leaves, Pieces<T>& out)
{
    for (const Leaf<T>& leaf : leaves)
    {
        out.vertices.insert(out.vertices.end(), leaf.x.begin(), leaf.x.end());
        out.offsets.push_back(out.offsets.back() + cell::get_num_vertices(leaf.type));
        out.types.push_back(leaf.type);
        out.negative.push_back(leaf.negative);
        out.positive.push_back(leaf.positive);
        out.zero.push_back(leaf.zero);
    }
}

/// Whether values at the vertices of a hexahedron (or quadrilateral, Basix
/// order) are an affine function's, up to rounding: their multilinear
/// interpolant has no bilinear or trilinear part.
template <std::floating_point T>
bool affine_values(const T* v, int tdim)
{
    T scale = T(0), gap = T(0);
    for (int b = 0; b < (1 << tdim); ++b)
        scale = std::max(scale, std::abs(v[b]));
    for (int b = 1; b < (1 << tdim); ++b)
    {
        T p = v[0];
        for (int d = 0; d < tdim; ++d)
            if ((b >> d) & 1)
                p += v[1 << d] - v[0];
        gap = std::max(gap, std::abs(v[b] - p));
    }
    return gap <= T(1024) * std::numeric_limits<T>::epsilon() * scale;
}

/// The pieces of one sub-cell: its volume cut by one level set after the
/// other, and the zero sets asked for, cut by the other level sets (with
/// @p curves also where two of them vanish).
/// @param sv  the level sets at the sub-cell's vertices, level set by level set
template <std::floating_point T>
void cut_sub_cell(cell::type type, int tdim, const std::vector<T>& x, const std::vector<T>& sv, int n_level_sets,
                  std::uint64_t zero_sets, bool curves, cell::TriangulationStrategy triangulation, Pieces<T>& out)
{
    const int nv = cell::get_num_vertices(type);
    const std::uint64_t curve_sets = curves ? zero_sets : 0;
    SubCell<T> s;
    make_sub_cell(type, tdim, x.data(), s);
    std::vector<Leaf<T>> leaves(1, Leaf<T>{type, x, 0, 0, 0});
    for (int i = 0; i < n_level_sets; ++i)
        cut_leaves(leaves, s, sv.data() + i * nv, i, 0, triangulation);
    append_leaves(leaves, out);

    for (int i = 0; i < n_level_sets; ++i)
    {
        if (!((zero_sets >> i) & 1))
            continue;
        std::vector<T> v(sv.begin() + i * nv, sv.begin() + (i + 1) * nv);
        snap_values(v);
        if (sign_class(std::span<const T>(v)) != 0)
            continue;
        cell::CutCell<T> part;
        cell::cut<T>(type, std::span<const T>(x), tdim, std::span<const T>(v), "phi=0", part, triangulation);
        leaves.clear();
        append_part(part, type, tdim, 0, 0, std::uint64_t(1) << i, leaves);
        for (int j = 0; j < n_level_sets; ++j)
            if (j != i)
                cut_leaves(leaves, s, sv.data() + j * nv, j, curve_sets, triangulation);
        append_leaves(leaves, out);
    }
}

} // namespace

std::span<const double> template_vertices(cell::type cell_type, int template_order)
{
    return std::span<const double>(iso_p1_template(cell_type, template_order).ref_vertex_coords);
}

template <std::floating_point T>
void cut_cell(cell::type cell_type, int template_order, std::span<const T> values, int n_level_sets,
              std::uint64_t zero_sets, bool curves, cell::TriangulationStrategy triangulation, Pieces<T>& out)
{
    if (n_level_sets < 1 || n_level_sets > 64)
        throw std::invalid_argument("lut: give 1 to 64 level sets");
    const IsoRefineTemplate& tpl = iso_p1_template(cell_type, template_order);
    const int tdim = tpl.tdim, nvt = tpl.n_vertices, vpc = tpl.vertices_per_cell;
    if (static_cast<int>(values.size()) != n_level_sets * nvt)
        throw std::invalid_argument("lut: the level sets need one value per template vertex");

    out = Pieces<T>{};
    out.tdim = tdim;
    const cell::type type = tpl.child_cell_type;
    std::vector<T> x(static_cast<std::size_t>(vpc * tdim)), sv(static_cast<std::size_t>(n_level_sets * vpc));
    std::vector<T> xs, svs;
    for (int c = 0; c < tpl.n_cells; ++c)
    {
        const int* ids = tpl.cell_connectivity.data() + c * vpc;
        for (int j = 0; j < vpc; ++j)
        {
            for (int d = 0; d < tdim; ++d)
                x[static_cast<std::size_t>(j * tdim + d)]
                    = static_cast<T>(tpl.ref_vertex_coords[static_cast<std::size_t>(ids[j] * tdim + d)]);
            for (int i = 0; i < n_level_sets; ++i)
                sv[static_cast<std::size_t>(i * vpc + j)] = values[static_cast<std::size_t>(i * nvt + ids[j])];
        }
        bool affine = true;
        if (type == cell::type::hexahedron)
            for (int i = 0; affine && i < n_level_sets; ++i)
                affine = affine_values(sv.data() + i * vpc, tdim);
        if (affine)
        {
            cut_sub_cell(type, tdim, x, sv, n_level_sets, zero_sets, curves, triangulation, out);
            continue;
        }
        // multilinear values: the hexahedron's tables for both sides do not fit
        // together there (up to 0.8% of a cell), P1 on its Kuhn tetrahedra does
        int local[8] = {0, 1, 2, 3, 4, 5, 6, 7};
        std::vector<std::vector<int>> simplices;
        cell::triangulation(type, local, simplices);
        const cell::type simplex = cell::type::tetrahedron;
        for (const std::vector<int>& t : simplices)
        {
            xs.clear();
            svs.clear();
            for (const int j : t)
                xs.insert(xs.end(), x.begin() + j * tdim, x.begin() + (j + 1) * tdim);
            for (int i = 0; i < n_level_sets; ++i)
                for (const int j : t)
                    svs.push_back(sv[static_cast<std::size_t>(i * vpc + j)]);
            cut_sub_cell(simplex, tdim, xs, svs, n_level_sets, zero_sets, curves, triangulation, out);
        }
    }
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template void cut_cell<float>(cell::type, int, std::span<const float>, int, std::uint64_t, bool,
                              cell::TriangulationStrategy, Pieces<float>&);
template void cut_cell<double>(cell::type, int, std::span<const double>, int, std::uint64_t, bool,
                               cell::TriangulationStrategy, Pieces<double>&);

} // namespace cutcells::lut
