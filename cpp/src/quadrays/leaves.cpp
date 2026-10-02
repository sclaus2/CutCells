// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "leaves.h"

#include <algorithm>
#include <array>
#include <map>
#include <stdexcept>

#include "../reference_cell.h"
#include "../write_vtk.h"

namespace cutcells::quadrays
{

template <std::floating_point T>
void append_leaves(const ClippedBox<T>& cell, const Source<T>& phi, Part part, int degree, const Options& opt,
                   std::int32_t parent_cell, LeafMesh<T>& mesh, Stats& stats)
{
    thread_local CellPoints<T> nodes;
    nodes.points.clear();
    nodes.weights.clear();
    nodes.tags.clear();
    leaf_nodes(cell, phi, part, degree, opt, nodes, stats);

    const bool surface = part == Part::interface;
    const int p1 = degree + 1;
    const int n_nodes = surface ? p1 * p1 : p1 * p1 * p1;

    struct Leaf
    {
        std::vector<Vec3<T>> u;
        std::vector<char> set;
        T phi_sum = 0;
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
        const int index = surface ? tag.node[0] + p1 * tag.node[1]
                                  : tag.node[0] + p1 * (tag.node[1] + p1 * tag.node[2]);
        if (index < 0 || index >= n_nodes || leaf.set[index])
            continue;
        const Vec3<T> u = {nodes.points[3 * i], nodes.points[3 * i + 1], nodes.points[3 * i + 2]};
        leaf.u[index] = u;
        leaf.set[index] = 1;
        leaf.phi_sum += evaluate(phi, std::span<const T>(u));
    }

    const Mat3<T> inv = surface ? inverse_jacobian(cell) : Mat3<T>{};
    std::vector<std::int32_t> conn(n_nodes);
    std::vector<Vec3<T>> x(n_nodes);
    for (const auto& [key, leaf] : leaves)
    {
        if (std::count(leaf.set.begin(), leaf.set.end(), char(1)) != n_nodes)
        {
            ++stats.incomplete_leaves; // node counts differed across the leaf
            continue;
        }
        // a volume leaf lies on one side of the level set; nodes on the interface are ~0
        if (!surface && part == Part::negative && !(leaf.phi_sum < T(0)))
            continue;
        if (!surface && part == Part::positive && !(leaf.phi_sum > T(0)))
            continue;

        const std::int32_t first = mesh.n_points();
        for (int index = 0; index < n_nodes; ++index)
        {
            x[index] = physical_point(cell, leaf.u[index]);
            mesh.points.insert(mesh.points.end(), x[index].begin(), x[index].end());
        }
        // Orientation: hexahedra with a positive Jacobian, interface quadrilaterals
        // with their normal along grad phi. Summing over all corners keeps the sign
        // reliable for collapsed leaves.
        auto node = [&](int i, int j, int k) -> const Vec3<T>&
        { return x[surface ? i + p1 * j : i + p1 * (j + p1 * k)]; };
        auto sub = [](const Vec3<T>& a, const Vec3<T>& b) { return Vec3<T>{a[0] - b[0], a[1] - b[1], a[2] - b[2]}; };
        auto cross = [](const Vec3<T>& a, const Vec3<T>& b)
        { return Vec3<T>{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; };
        T orientation = T(0);
        if (surface)
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
            // physical gradient of phi at the leaf centre: J^-T grad_u phi
            Vec3<T> gu;
            gradient(phi, std::span<const T>(leaf.u[n_nodes / 2]), std::span<T>(gu));
            for (int d = 0; d < 3; ++d)
            {
                T gx = T(0);
                for (int e = 0; e < 3; ++e)
                    gx += inv[e][d] * gu[e];
                orientation += normal[d] * gx;
            }
        }
        else
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
        const bool flip = orientation < T(0);
        if (surface)
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
        mesh.vtk_types.push_back(surface ? vtk_lagrange_quadrilateral : vtk_lagrange_hexahedron);
        mesh.parent.push_back(parent_cell);
        mesh.degree.push_back(degree);
    }
}

template <std::floating_point T>
void append_leaves(const ClippedBox<T>& cell, const BoxBernstein<T>& phi, Part part, int degree,
                   const Options& opt, std::int32_t parent_cell, LeafMesh<T>& mesh, Stats& stats)
{
    append_leaves(cell, bernstein_source(phi), part, degree, opt, parent_cell, mesh, stats);
}

template <std::floating_point T>
void append_linear_cell(cell::type cell_type, std::span<const T> vertex_coords,
                        std::int32_t parent_cell, LeafMesh<T>& mesh)
{
    if (cell_type != cell::type::tetrahedron && cell_type != cell::type::hexahedron)
    {
        throw std::invalid_argument("quadrays: unsupported cell type "
                                    + cell::cell_type_to_str(cell_type));
    }
    const std::span<const int> perm = cell::basix_to_vtk_vertex_permutation(cell_type);
    const std::int32_t first = mesh.n_points();
    for (std::size_t j = 0; j < perm.size(); ++j)
    {
        const T* v = vertex_coords.data() + 3 * perm[j];
        mesh.points.insert(mesh.points.end(), v, v + 3);
        mesh.connectivity.push_back(first + static_cast<std::int32_t>(j));
    }
    if (cell_type == cell::type::tetrahedron)
    {
        // VTK tetrahedra are positively oriented
        ClippedBox<T> box;
        make_clipped_box(cell_type, vertex_coords, 3, box);
        if (jacobian_determinant(box) < T(0))
            std::swap(mesh.connectivity[first + 1], mesh.connectivity[first + 2]);
    }
    mesh.offsets.push_back(static_cast<std::int32_t>(mesh.connectivity.size()));
    mesh.vtk_types.push_back(cell_type == cell::type::tetrahedron ? vtk_tetra : vtk_hexahedron);
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

template void append_leaves<float>(const ClippedBox<float>&, const Source<float>&, Part, int, const Options&,
                                   std::int32_t, LeafMesh<float>&, Stats&);
template void append_leaves<double>(const ClippedBox<double>&, const Source<double>&, Part, int, const Options&,
                                    std::int32_t, LeafMesh<double>&, Stats&);
template void append_leaves<float>(const ClippedBox<float>&, const BoxBernstein<float>&, Part, int,
                                   const Options&, std::int32_t, LeafMesh<float>&, Stats&);
template void append_leaves<double>(const ClippedBox<double>&, const BoxBernstein<double>&, Part, int,
                                    const Options&, std::int32_t, LeafMesh<double>&, Stats&);
template void append_linear_cell<float>(cell::type, std::span<const float>, std::int32_t, LeafMesh<float>&);
template void append_linear_cell<double>(cell::type, std::span<const double>, std::int32_t,
                                         LeafMesh<double>&);
template void write_leaves<float>(const std::string&, const LeafMesh<float>&);
template void write_leaves<double>(const std::string&, const LeafMesh<double>&);

} // namespace cutcells::quadrays
