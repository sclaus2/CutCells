// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "output.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "../cell_topology.h"
#include "../quadrature_tables.h"
#include "cell_source.h"

namespace cutcells::part
{

namespace
{

void check_backend(const std::string& backend)
{
    if (backend == "quadrays")
        return;
    if (backend == "lut")
        throw std::invalid_argument("part: the backend 'lut' comes in phase 4");
    throw std::invalid_argument("part: unknown backend '" + backend + "'; expected 'quadrays'");
}

/// The level set bounding a part's piece of a cut cell, and the parts of it
/// that the expression's terms select there.
struct Pieces
{
    int level_set = -1;
    std::vector<quadrays::Part> parts;
};

template <std::floating_point T, std::integral I>
Pieces pieces_in_cell(const MeshPart<T, I>& part, I cell_id)
{
    const CutResult<T, I>& r = *part.result;
    const std::uint64_t cut = cut_mask(r, cell_id);
    Pieces p;
    for (const SelectionTerm& term : part.expr.terms)
    {
        if (term_on_cell(term, r, cell_id) != TermCell::piece)
            continue;
        const std::uint64_t bounding
            = (term.negative_required | term.positive_required | term.zero_required) & cut;
        if (std::popcount(bounding) != 1)
        {
            throw std::runtime_error("part: quadrays takes one level set per cell, and cell "
                                     + std::to_string(cell_id) + " is cut by "
                                     + std::to_string(std::popcount(bounding))
                                     + " level sets of a term (several level sets per cell come in phase 6)");
        }
        const int l = std::countr_zero(bounding);
        if (p.level_set >= 0 && p.level_set != l)
        {
            throw std::runtime_error("part: terms bounded by different level sets select pieces of cell "
                                     + std::to_string(cell_id)
                                     + " (several level sets per cell come in phase 6)");
        }
        p.level_set = l;
        const std::uint64_t bit = std::uint64_t(1) << l;
        const quadrays::Part q = (term.zero_required & bit)       ? quadrays::Part::interface
                                 : (term.negative_required & bit) ? quadrays::Part::negative
                                                                  : quadrays::Part::positive;
        if (std::find(p.parts.begin(), p.parts.end(), q) == p.parts.end())
            p.parts.push_back(q);
    }
    return p;
}

/// The reference rule of a whole cell, in box coordinates with physical weights.
template <std::floating_point T>
void append_cell_points(const quadrays::ClippedBox<T>& box, cell::type type, int degree,
                        quadrays::CellPoints<T>& out)
{
    const auto rule = quadrature::get_reference_rule<T>(type, degree);
    const T detj = std::abs(quadrays::jacobian_determinant(box));
    out.points.insert(out.points.end(), rule._points.begin(), rule._points.end());
    for (const T w : rule._weights)
        out.weights.push_back(w * detj);
}

/// The reference rule of face @p f of a cell, in box coordinates with
/// physical surface weights.
template <std::floating_point T>
void append_face_points(const quadrays::ClippedBox<T>& box, cell::type type, int f, int degree,
                        quadrays::CellPoints<T>& out)
{
    const std::span<const int> fv = cell::face_vertices(type, f);
    const quadrays::Vec3<T> ua = box_vertex<T>(type, fv[0]), ub = box_vertex<T>(type, fv[1]),
                            uc = box_vertex<T>(type, fv[2]);
    // physical edges of the face and its area factor
    std::array<T, 3> e1{}, e2{};
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
        {
            e1[i] += box.jacobian[i][k] * (ub[k] - ua[k]);
            e2[i] += box.jacobian[i][k] * (uc[k] - ua[k]);
        }
    const T area = std::sqrt((e1[1] * e2[2] - e1[2] * e2[1]) * (e1[1] * e2[2] - e1[2] * e2[1])
                             + (e1[2] * e2[0] - e1[0] * e2[2]) * (e1[2] * e2[0] - e1[0] * e2[2])
                             + (e1[0] * e2[1] - e1[1] * e2[0]) * (e1[0] * e2[1] - e1[1] * e2[0]));
    const auto rule = quadrature::get_reference_rule<T>(cell::face_type(type, f), degree);
    for (int q = 0; q < rule._num_points; ++q)
    {
        const T s = rule._points[2 * q], t = rule._points[2 * q + 1];
        for (int i = 0; i < 3; ++i)
            out.points.push_back(ua[i] + s * (ub[i] - ua[i]) + t * (uc[i] - ua[i]));
        out.weights.push_back(rule._weights[q] * area);
    }
}

/// What a part has in one cell.
template <std::integral I>
struct Entry
{
    I cell;
    int kind;  ///< 0: a piece, 1: a zero face, 2: the whole cell
    int index; ///< the zero face
};

template <std::floating_point T, std::integral I>
std::vector<Entry<I>> part_entries(const MeshPart<T, I>& part, bool include_uncut_cells)
{
    const CutResult<T, I>& r = *part.result;
    std::vector<Entry<I>> entries;
    for (const I c : part.cut_cells)
        entries.push_back({c, 0, -1});
    for (const int z : part.zero_faces)
        entries.push_back({r.zero_face_cells[static_cast<std::size_t>(z)], 1, z});
    if (include_uncut_cells)
        for (const I c : part.uncut_cells)
            entries.push_back({c, 2, -1});
    std::stable_sort(entries.begin(), entries.end(),
                     [](const Entry<I>& a, const Entry<I>& b) { return a.cell < b.cell; });
    return entries;
}

} // namespace

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> quadrature_rules(const MeshPart<T, I>& part, int order, bool include_uncut_cells,
                                                const std::string& backend, const quadrays::Options& options)
{
    check_backend(backend);
    if (order < 1)
        throw std::invalid_argument("part: the quadrature order must be at least 1");
    const CutResult<T, I>& r = *part.result;
    const MeshView<T, I>& mesh = *r.mesh;
    const int degree = std::min(2 * order - 1, 10);

    quadrature::QuadratureRules<T> rules;
    rules._tdim = r.num_cells > 0 ? cell::get_tdim(mesh.cell_type(I(0))) : 3;
    rules._offset.push_back(0);
    const std::vector<Entry<I>> entries = part_entries(part, include_uncut_cells);
    CellSource<T, I> cs;
    quadrays::ClippedBox<T> box;
    quadrays::CellPoints<T> points;
    quadrays::Stats stats;
    for (std::size_t i = 0; i < entries.size();)
    {
        const I c = entries[i].cell;
        const cell::type type = mesh.cell_type(c);
        if (mesh.gdim != 3 || (type != cell::type::tetrahedron && type != cell::type::hexahedron))
            throw std::invalid_argument("part: quadrays takes tetrahedra and hexahedra in 3D");
        cell_vertex_coords_basix(mesh, c, cs.vertices, cs.nodes);
        quadrays::make_clipped_box<T>(type, std::span<const T>(cs.vertices), 3, box);
        points.points.clear();
        points.weights.clear();
        for (; i < entries.size() && entries[i].cell == c; ++i)
        {
            if (entries[i].kind == 0)
            {
                const Pieces p = pieces_in_cell(part, c);
                cell_source(mesh, *r.level_sets[static_cast<std::size_t>(p.level_set)], c, cs);
                for (const quadrays::Part q : p.parts)
                    quadrays::integrate(cs.box, cs.source, q, order, options, points, stats);
            }
            else if (entries[i].kind == 1)
            {
                const int f = r.zero_face_local[static_cast<std::size_t>(entries[i].index)];
                append_face_points(box, type, f, degree, points);
            }
            else
                append_cell_points(box, type, degree, points);
        }
        if (points.n_points() == 0)
            continue;
        for (int k = 0; k < points.n_points(); ++k)
        {
            const quadrays::Vec3<T> u = {points.points[3 * k], points.points[3 * k + 1], points.points[3 * k + 2]};
            const quadrays::Vec3<T> xi = quadrays::reference_point(box, u);
            rules._points.insert(rules._points.end(), xi.begin(), xi.end());
            rules._weights.push_back(points.weights[k]);
        }
        rules._offset.push_back(static_cast<std::int32_t>(rules._weights.size()));
        rules._parent_map.push_back(static_cast<std::int32_t>(c));
    }
    return rules;
}

template <std::floating_point T, std::integral I>
quadrays::LeafMesh<T> visualization_mesh(const MeshPart<T, I>& part, int degree, bool include_uncut_cells,
                                         const std::string& backend, const quadrays::Options& options)
{
    check_backend(backend);
    const CutResult<T, I>& r = *part.result;
    const MeshView<T, I>& mesh = *r.mesh;
    quadrays::LeafMesh<T> leaves;
    CellSource<T, I> cs;
    quadrays::Stats stats;
    std::vector<T> face;
    for (const Entry<I>& e : part_entries(part, include_uncut_cells))
    {
        const cell::type type = mesh.cell_type(e.cell);
        if (e.kind == 0)
        {
            const Pieces p = pieces_in_cell(part, e.cell);
            if (!cell_source(mesh, *r.level_sets[static_cast<std::size_t>(p.level_set)], e.cell, cs))
                throw std::invalid_argument("part: quadrays takes tetrahedra and hexahedra in 3D");
            for (const quadrays::Part q : p.parts)
                quadrays::append_leaves(cs.box, cs.source, q, degree, options, static_cast<std::int32_t>(e.cell),
                                        leaves, stats);
            continue;
        }
        cell_vertex_coords_basix(mesh, e.cell, cs.vertices, cs.nodes);
        if (e.kind == 1)
        {
            const int f = r.zero_face_local[static_cast<std::size_t>(e.index)];
            const std::span<const int> fv = cell::face_vertices(type, f);
            face.clear();
            for (const int v : fv)
                face.insert(face.end(), cs.vertices.begin() + 3 * v, cs.vertices.begin() + 3 * v + 3);
            quadrays::append_linear_cell<T>(cell::face_type(type, f), std::span<const T>(face),
                                            static_cast<std::int32_t>(e.cell), leaves);
        }
        else
            quadrays::append_linear_cell<T>(type, std::span<const T>(cs.vertices), static_cast<std::int32_t>(e.cell),
                                            leaves);
    }
    return leaves;
}

template <std::floating_point T, std::integral I>
void write_vtu(const std::string& filename, const MeshPart<T, I>& part, int degree, bool include_uncut_cells,
               const std::string& backend, const quadrays::Options& options)
{
    quadrays::write_leaves(filename, visualization_mesh(part, degree, include_uncut_cells, backend, options));
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template quadrature::QuadratureRules<float> quadrature_rules<float, int>(const MeshPart<float, int>&, int, bool,
                                                                         const std::string&,
                                                                         const quadrays::Options&);
template quadrature::QuadratureRules<double> quadrature_rules<double, int>(const MeshPart<double, int>&, int, bool,
                                                                           const std::string&,
                                                                           const quadrays::Options&);
template quadrays::LeafMesh<float> visualization_mesh<float, int>(const MeshPart<float, int>&, int, bool,
                                                                  const std::string&, const quadrays::Options&);
template quadrays::LeafMesh<double> visualization_mesh<double, int>(const MeshPart<double, int>&, int, bool,
                                                                    const std::string&, const quadrays::Options&);
template void write_vtu<float, int>(const std::string&, const MeshPart<float, int>&, int, bool, const std::string&,
                                    const quadrays::Options&);
template void write_vtu<double, int>(const std::string&, const MeshPart<double, int>&, int, bool,
                                     const std::string&, const quadrays::Options&);

} // namespace cutcells::part
