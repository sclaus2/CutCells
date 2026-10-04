// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "leaves.h"

#include <algorithm>
#include <array>
#include <bit>
#include <map>
#include <stdexcept>

#include "../reference_cell.h"
#include "../write_vtk.h"

namespace cutcells::quadrays
{

namespace
{
/// Position of node i, 0 <= i <= p, of a degree-p Lagrange curve in VTK's node
/// order: the ends, then the interior.
int vtk_lagrange_curve_index(int i, int p)
{
    return i == 0 ? 0 : (i == p ? 1 : i + 1);
}
} // namespace

template <std::floating_point T>
void append_leaves(const ClippedBox<T>& cell, std::span<const Source<T>> phis, std::span<const SelectionTerm> terms,
                   int degree, const Options& opt, std::int32_t parent_cell, LeafMesh<T>& mesh, Stats& stats)
{
    thread_local CellPoints<T> nodes;
    nodes.points.clear();
    nodes.weights.clear();
    nodes.tags.clear();
    leaf_nodes(cell, phis, terms, degree, opt, nodes, stats);

    std::uint64_t zero = 0, signed_sets = 0;
    for (const SelectionTerm& t : terms)
    {
        zero |= t.zero_required;
        signed_sets |= t.negative_required | t.positive_required;
    }
    const bool surface = zero != 0;
    const int surface_ls = surface ? std::countr_zero(zero) : -1;
    const int partner_ls = std::popcount(zero) == 2 ? 63 - std::countl_zero(zero) : -1; // a curve's second
    const int dim = cell.tdim - std::popcount(zero);                                  // of the leaves
    const int p1 = degree + 1;
    int n_nodes = 1;
    for (int d = 0; d < dim; ++d)
        n_nodes *= p1;

    struct Leaf
    {
        std::vector<Vec3<T>> u;
        std::vector<char> set;
    };
    // a leaf is one segment per level of one certified box per level
    std::map<std::array<int, 6>, Leaf> leaves;
    for (int i = 0; i < nodes.n_points(); ++i)
    {
        const NodeTag& tag = nodes.tags[i];
        const std::array<int, 6> key
            = {tag.box[0], tag.box[1], tag.box[2], tag.segment[0], tag.segment[1], tag.segment[2]};
        Leaf& leaf = leaves[key];
        if (leaf.u.empty())
        {
            leaf.u.resize(n_nodes);
            leaf.set.assign(n_nodes, 0);
        }
        int index = 0;
        for (int d = 0, stride = 1; d < dim; ++d, stride *= p1)
            index += stride * tag.node[d];
        if (index < 0 || index >= n_nodes || leaf.set[index])
            continue;
        leaf.u[index] = {nodes.points[3 * i], nodes.points[3 * i + 1], nodes.points[3 * i + 2]};
        leaf.set[index] = 1;
    }

    const Mat3<T> inv = inverse_jacobian(cell);
    std::vector<std::int32_t> conn(n_nodes);
    std::vector<Vec3<T>> x(n_nodes);
    std::vector<int> sign(phis.size(), 0);
    auto sub = [](const Vec3<T>& a, const Vec3<T>& b) { return Vec3<T>{a[0] - b[0], a[1] - b[1], a[2] - b[2]}; };
    auto cross = [](const Vec3<T>& a, const Vec3<T>& b)
    { return Vec3<T>{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; };
    // physical gradient of a level set: J^-T grad_u phi
    auto physical_gradient = [&](int ls, const Vec3<T>& u)
    {
        Vec3<T> gu = {0, 0, 0}, gx = {0, 0, 0};
        gradient(phis[static_cast<std::size_t>(ls)], std::span<const T>(u), std::span<T>(gu));
        for (int d = 0; d < 3; ++d)
            for (int e = 0; e < 3; ++e)
                gx[d] += inv[e][d] * gu[e];
        return gx;
    };
    for (const auto& [key, leaf] : leaves)
    {
        // A leaf lies on one side of each level set that splits it; nodes on
        // a zero set are ~0. The terms decide by the sign of the sums, over
        // the nodes a leaf has.
        if (signed_sets != 0)
        {
            for (std::uint64_t bits = signed_sets; bits != 0; bits &= bits - 1)
            {
                const int l = std::countr_zero(bits);
                T sum = T(0);
                for (int index = 0; index < n_nodes; ++index)
                    if (leaf.set[index])
                        sum += evaluate(phis[static_cast<std::size_t>(l)], std::span<const T>(leaf.u[index]));
                sign[static_cast<std::size_t>(l)] = sum < T(0) ? -1 : (sum > T(0) ? 1 : 0);
            }
            bool selected = false;
            for (const SelectionTerm& t : terms)
            {
                bool ok = true;
                for (std::uint64_t bits = t.negative_required | t.positive_required; bits != 0 && ok; bits &= bits - 1)
                {
                    const int l = std::countr_zero(bits);
                    ok = sign[static_cast<std::size_t>(l)] == ((t.negative_required >> l) & 1 ? -1 : 1);
                }
                selected |= ok;
            }
            if (!selected)
                continue;
        }
        if (std::count(leaf.set.begin(), leaf.set.end(), char(1)) != n_nodes)
        {
            ++stats.incomplete_leaves; // node counts differed across the leaf
            continue;
        }

        const std::int32_t first = mesh.n_points();
        for (int index = 0; index < n_nodes; ++index)
        {
            x[index] = physical_point(cell, leaf.u[index]);
            mesh.points.insert(mesh.points.end(), x[index].begin(), x[index].end());
        }
        // Orientation: volume leaves with a positive Jacobian, interface leaves
        // with their normal along the gradient. Summing over all corners keeps the
        // sign reliable for collapsed leaves.
        auto node = [&](int i, int j, int k) -> const Vec3<T>& { return x[i + p1 * (j + p1 * k)]; };
        T orientation = T(0);
        if (dim == 3)
        {
            for (int c = 0; c < 8; ++c)
            {
                const int i0 = (c & 1) * degree, j0 = ((c >> 1) & 1) * degree, k0 = ((c >> 2) & 1) * degree;
                const T s = (i0 ? T(-1) : T(1)) * (j0 ? T(-1) : T(1)) * (k0 ? T(-1) : T(1));
                const Vec3<T>& o = node(i0, j0, k0);
                const Vec3<T> n = cross(sub(node(degree - i0, j0, k0), o), sub(node(i0, degree - j0, k0), o));
                const Vec3<T> e = sub(node(i0, j0, degree - k0), o);
                orientation += s * (n[0] * e[0] + n[1] * e[1] + n[2] * e[2]);
            }
        }
        else if (dim == 2)
        {
            Vec3<T> normal = {0, 0, 0};
            for (int c = 0; c < 4; ++c)
            {
                const int i0 = (c & 1) * degree, j0 = ((c >> 1) & 1) * degree;
                const T s = (i0 ? T(-1) : T(1)) * (j0 ? T(-1) : T(1));
                const Vec3<T> n = cross(sub(node(degree - i0, j0, 0), node(i0, j0, 0)),
                                        sub(node(i0, degree - j0, 0), node(i0, j0, 0)));
                for (int d = 0; d < 3; ++d)
                    normal[d] += s * n[d];
            }
            if (surface)
            {
                const Vec3<T> gx = physical_gradient(surface_ls, leaf.u[n_nodes / 2]);
                for (int d = 0; d < 3; ++d)
                    orientation += normal[d] * gx[d];
            }
            else
                orientation = normal[2]; // a 2D cell's leaf, in the plane z = 0
        }
        else if (dim == 1)
        {
            const Vec3<T> t = sub(x[degree], x[0]);
            const Vec3<T> gx = physical_gradient(surface_ls, leaf.u[n_nodes / 2]);
            if (partner_ls >= 0)
            {
                // where two level sets vanish: along grad a x grad b
                const Vec3<T> n = cross(gx, physical_gradient(partner_ls, leaf.u[n_nodes / 2]));
                orientation = t[0] * n[0] + t[1] * n[1] + t[2] * n[2];
            }
            else
            {
                // a curve in the plane: its normal (t_y, -t_x) along the gradient
                orientation = t[1] * gx[0] - t[0] * gx[1];
            }
        }
        const bool flip = orientation < T(0);
        if (dim == 0)
            conn[0] = first; // a point where two curves cross
        else if (dim == 1)
        {
            for (int i = 0; i < p1; ++i)
                conn[vtk_lagrange_curve_index(flip ? degree - i : i, degree)] = first + i;
        }
        else if (dim == 2)
        {
            for (int j = 0; j < p1; ++j)
                for (int i = 0; i < p1; ++i)
                    conn[io::vtk_lagrange_quadrilateral_index(flip ? degree - i : i, j, degree)]
                        = first + i + p1 * j;
        }
        else
        {
            for (int k = 0; k < p1; ++k)
                for (int j = 0; j < p1; ++j)
                    for (int i = 0; i < p1; ++i)
                        conn[io::vtk_lagrange_hexahedron_index(flip ? degree - i : i, j, k, degree)]
                            = first + i + p1 * (j + p1 * k);
        }
        mesh.connectivity.insert(mesh.connectivity.end(), conn.begin(), conn.end());
        mesh.offsets.push_back(static_cast<std::int32_t>(mesh.connectivity.size()));
        mesh.vtk_types.push_back(dim == 0   ? vtk_vertex
                                 : dim == 1 ? vtk_lagrange_curve
                                 : dim == 2 ? vtk_lagrange_quadrilateral
                                            : vtk_lagrange_hexahedron);
        mesh.parent.push_back(parent_cell);
        mesh.degree.push_back(dim == 0 ? 0 : degree);
    }
}

template <std::floating_point T>
void append_leaves(const ClippedBox<T>& cell, const Source<T>& phi, Part part, int degree, const Options& opt,
                   std::int32_t parent_cell, LeafMesh<T>& mesh, Stats& stats)
{
    const SelectionTerm term = term_of(part);
    append_leaves(cell, std::span<const Source<T>>(&phi, 1), std::span<const SelectionTerm>(&term, 1), degree, opt,
                  parent_cell, mesh, stats);
}

template <std::floating_point T>
void append_leaves(const ClippedBox<T>& cell, const BoxBernstein<T>& phi, Part part, int degree,
                   const Options& opt, std::int32_t parent_cell, LeafMesh<T>& mesh, Stats& stats)
{
    append_leaves(cell, bernstein_source(phi), part, degree, opt, parent_cell, mesh, stats);
}

template <std::floating_point T>
void append_linear_cell(cell::type cell_type, std::span<const T> vertex_coords, int gdim,
                        std::int32_t parent_cell, LeafMesh<T>& mesh)
{
    std::uint8_t vtk_type = 0;
    switch (cell_type)
    {
    case cell::type::interval:
        vtk_type = vtk_line;
        break;
    case cell::type::triangle:
        vtk_type = vtk_triangle;
        break;
    case cell::type::quadrilateral:
        vtk_type = vtk_quad;
        break;
    case cell::type::tetrahedron:
        vtk_type = vtk_tetra;
        break;
    case cell::type::hexahedron:
        vtk_type = vtk_hexahedron;
        break;
    case cell::type::prism:
        vtk_type = vtk_wedge;
        break;
    case cell::type::pyramid:
        vtk_type = vtk_pyramid;
        break;
    default:
        throw std::invalid_argument("quadrays: unsupported cell type " + cell::cell_type_to_str(cell_type));
    }
    if (gdim < 1 || gdim > 3)
        throw std::invalid_argument("quadrays: points have 1 to 3 coordinates");
    const std::span<const int> perm = cell::basix_to_vtk_vertex_permutation(cell_type);
    const std::int32_t first = mesh.n_points();
    for (std::size_t j = 0; j < perm.size(); ++j)
    {
        const T* v = vertex_coords.data() + gdim * perm[j];
        for (int d = 0; d < 3; ++d)
            mesh.points.push_back(d < gdim ? v[d] : T(0));
        mesh.connectivity.push_back(first + static_cast<std::int32_t>(j));
    }
    if (cell_type == cell::type::tetrahedron)
    {
        // VTK tetrahedra are positively oriented: the sign of the reference map
        ClippedBox<T> box;
        make_clipped_box(cell_type, vertex_coords, gdim, box, BoxFrame::reference);
        if (jacobian_determinant(box) < T(0))
            std::swap(mesh.connectivity[first + 1], mesh.connectivity[first + 2]);
    }
    mesh.offsets.push_back(static_cast<std::int32_t>(mesh.connectivity.size()));
    mesh.vtk_types.push_back(vtk_type);
    mesh.parent.push_back(parent_cell);
    mesh.degree.push_back(1);
}

template <std::floating_point T>
void write_leaves(const std::string& filename, const LeafMesh<T>& mesh)
{
    const std::vector<int> connectivity(mesh.connectivity.begin(), mesh.connectivity.end());
    const std::vector<int> offsets(mesh.offsets.begin(), mesh.offsets.end());
    const std::vector<int> types(mesh.vtk_types.begin(), mesh.vtk_types.end());
    io::write_lagrange_vtk<T>(filename, std::span<const T>(mesh.points), connectivity, offsets, types, 3,
                              std::span<const std::int32_t>(mesh.parent), {},
                              std::span<const std::int32_t>(mesh.degree));
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template void append_leaves<float>(const ClippedBox<float>&, std::span<const Source<float>>,
                                   std::span<const SelectionTerm>, int, const Options&, std::int32_t,
                                   LeafMesh<float>&, Stats&);
template void append_leaves<double>(const ClippedBox<double>&, std::span<const Source<double>>,
                                    std::span<const SelectionTerm>, int, const Options&, std::int32_t,
                                    LeafMesh<double>&, Stats&);
template void append_leaves<float>(const ClippedBox<float>&, const Source<float>&, Part, int, const Options&,
                                   std::int32_t, LeafMesh<float>&, Stats&);
template void append_leaves<double>(const ClippedBox<double>&, const Source<double>&, Part, int, const Options&,
                                    std::int32_t, LeafMesh<double>&, Stats&);
template void append_leaves<float>(const ClippedBox<float>&, const BoxBernstein<float>&, Part, int,
                                   const Options&, std::int32_t, LeafMesh<float>&, Stats&);
template void append_leaves<double>(const ClippedBox<double>&, const BoxBernstein<double>&, Part, int,
                                    const Options&, std::int32_t, LeafMesh<double>&, Stats&);
template void append_linear_cell<float>(cell::type, std::span<const float>, int, std::int32_t, LeafMesh<float>&);
template void append_linear_cell<double>(cell::type, std::span<const double>, int, std::int32_t,
                                         LeafMesh<double>&);
template void write_leaves<float>(const std::string&, const LeafMesh<float>&);
template void write_leaves<double>(const std::string&, const LeafMesh<double>&);

} // namespace cutcells::quadrays
