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

#include "../bernstein.h"
#include "../cell_topology.h"
#include "../lut/piece_rules.h"
#include "../quadrature_tables.h"
#include "../quadrays/analytic.h"
#include "../reference_cell.h"
#include "../write_vtk.h"
#include "cell_source.h"

namespace cutcells::part
{

namespace
{

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

// ============================================================================
// The lookup-table backend
// ============================================================================

template <std::floating_point T, std::integral I>
void check_lut(const MeshPart<T, I>& part, const lut::Options& options)
{
    if (options.template_order < 0 || options.template_order > 4)
        throw std::invalid_argument("part: the template order of the lookup tables goes from 1 to 4 "
                                    "(0: from the level sets)");
    const CutResult<T, I>& r = *part.result;
    if (r.num_cells > 0 && part.dim < cell::get_tdim(r.mesh->cell_type(I(0))) - 1)
        throw std::invalid_argument("part: the lookup tables give volumes and interfaces, not the sets "
                                    "where two level sets vanish");
}

/// The map of a cell from its vertices in Basix order.
template <std::floating_point T, std::integral I>
void cell_map(const MeshView<T, I>& mesh, I cell_id, lut::CellMap<T>& map, std::vector<I>& scratch)
{
    map.type = mesh.cell_type(cell_id);
    if (map.type != cell::type::triangle && map.type != cell::type::quadrilateral
        && map.type != cell::type::tetrahedron && map.type != cell::type::hexahedron)
    {
        throw std::invalid_argument("part: the lookup tables take triangles, quadrilaterals, tetrahedra and "
                                    "hexahedra");
    }
    map.gdim = mesh.gdim;
    cell_vertex_coords_basix(mesh, cell_id, map.vertices, scratch);
}

/// A cut cell of a part as the lookup tables see it.
template <std::floating_point T>
struct LutCell
{
    std::vector<int> level_sets; ///< the level sets that cut the cell and that the expression names
    std::vector<T> xi, x;        ///< the template's vertices, reference and physical
    std::vector<T> values;       ///< the level sets there
    lut::Pieces<T> pieces;
    std::vector<char> selected;  ///< per piece: in the part
};

template <std::floating_point T>
std::span<const T> piece_vertices(const lut::Pieces<T>& pieces, int p)
{
    const std::size_t tdim = static_cast<std::size_t>(pieces.tdim);
    return std::span<const T>(pieces.vertices)
        .subspan(static_cast<std::size_t>(pieces.offsets[p]) * tdim,
                 static_cast<std::size_t>(pieces.offsets[p + 1] - pieces.offsets[p]) * tdim);
}

/// Whether a term holds on piece @p p: by the piece's sides for the level sets
/// cutting the cell, by the cell's domains for the others.
template <std::floating_point T, std::integral I>
bool term_on_piece(const SelectionTerm& term, const CutResult<T, I>& r, I cell_id, const LutCell<T>& lc, int p)
{
    const int zero = lc.pieces.zero[static_cast<std::size_t>(p)];
    // volume pieces answer the terms without zero clause, a zero piece those with its own
    if ((term.zero_required == 0) != (zero < 0))
        return false;
    const std::uint64_t all = term.negative_required | term.positive_required | term.zero_required;
    for (int l = 0; l < r.n_level_sets(); ++l)
    {
        const std::uint64_t bit = std::uint64_t(1) << l;
        if (!(all & bit))
            continue;
        if ((term.negative_required & bit) && (term.positive_required & bit))
            return false;
        const auto it = std::find(lc.level_sets.begin(), lc.level_sets.end(), l);
        if (it == lc.level_sets.end())
        {
            const cell::domain side = (term.negative_required & bit) ? cell::domain::inside : cell::domain::outside;
            if ((term.zero_required & bit) || r.domain(l, cell_id) != side)
                return false;
            continue;
        }
        const int j = static_cast<int>(it - lc.level_sets.begin());
        const std::uint64_t local = std::uint64_t(1) << j;
        if (term.zero_required & bit)
        {
            if (zero != j)
                return false;
        }
        else if (term.negative_required & bit)
        {
            if (!(lc.pieces.negative[static_cast<std::size_t>(p)] & local))
                return false;
        }
        else if (!(lc.pieces.positive[static_cast<std::size_t>(p)] & local))
            return false;
    }
    return true;
}

/// The pieces of a cut cell: the template of the options' order (by default
/// the highest degree of the level sets, 2 for analytic ones), the values of
/// the level sets there (exact for analytic level sets), the lookup tables,
/// and which pieces the part selects.
template <std::floating_point T, std::integral I>
void lut_cell(const MeshPart<T, I>& part, I cell_id, const lut::CellMap<T>& map, const lut::Options& options,
              LevelSetCell<T, I>& scratch, LutCell<T>& lc)
{
    const CutResult<T, I>& r = *part.result;
    std::uint64_t named = 0, zero_named = 0;
    for (const SelectionTerm& term : part.expr.terms)
    {
        named |= term.negative_required | term.positive_required | term.zero_required;
        zero_named |= term.zero_required;
    }
    lc.level_sets.clear();
    std::uint64_t zero_sets = 0;
    int degree = 1;
    bool analytic = false;
    for (int l = 0; l < r.n_level_sets(); ++l)
    {
        if (!((named >> l) & 1) || r.domain(l, cell_id) != cell::domain::intersected)
            continue;
        if ((zero_named >> l) & 1)
            zero_sets |= std::uint64_t(1) << lc.level_sets.size();
        lc.level_sets.push_back(l);
        const LevelSetFunction<T, I>& ls = *r.level_sets[static_cast<std::size_t>(l)];
        analytic |= static_cast<bool>(ls.analytic);
        degree = std::max(degree, ls.analytic ? 2 : ls.mesh_data.degree);
    }
    const int k = options.template_order > 0 ? options.template_order : std::min(degree, 4);
    const int tdim = map.tdim();
    const std::span<const double> tv = lut::template_vertices(map.type, k);
    lc.xi.assign(tv.begin(), tv.end());
    const std::size_t nv = lc.xi.size() / static_cast<std::size_t>(tdim);
    if (analytic)
        lut::push_forward(map, std::span<const T>(lc.xi), lc.x);

    lc.values.clear();
    for (const int l : lc.level_sets)
    {
        const LevelSetFunction<T, I>& ls = *r.level_sets[static_cast<std::size_t>(l)];
        if (ls.analytic)
        {
            const quadrays::AnalyticLevelSet& phi = *ls.analytic;
            for (std::size_t v = 0; v < nv; ++v)
            {
                double x[3] = {0, 0, 0};
                for (int d = 0; d < map.gdim; ++d)
                    x[d] = static_cast<double>(lc.x[v * static_cast<std::size_t>(map.gdim) + d]);
                lc.values.push_back(static_cast<T>(phi.value(x, phi.context)));
            }
            continue;
        }
        make_cell_level_set(ls, cell_id, scratch);
        for (std::size_t v = 0; v < nv; ++v)
        {
            lc.values.push_back(bernstein::evaluate<T>(map.type, scratch.bernstein_order,
                                                       std::span<const T>(scratch.bernstein_coeffs),
                                                       std::span<const T>(lc.xi).subspan(v * tdim, tdim)));
        }
    }
    lc.selected.clear();
    if (lc.level_sets.empty())
    {
        lc.pieces = lut::Pieces<T>{};
        return;
    }
    lut::cut_cell<T>(map.type, k, std::span<const T>(lc.values), static_cast<int>(lc.level_sets.size()),
                     part.dim < tdim ? zero_sets : 0, options.triangulate, lc.pieces);
    lc.selected.assign(static_cast<std::size_t>(lc.pieces.n_pieces()), 0);
    for (int p = 0; p < lc.pieces.n_pieces(); ++p)
        for (const SelectionTerm& term : part.expr.terms)
            if (term_on_piece(term, r, cell_id, lc, p))
            {
                lc.selected[static_cast<std::size_t>(p)] = 1;
                break;
            }
}

/// Face @p f of a cell in the cell's reference coordinates.
template <std::floating_point T>
void reference_face(cell::type type, int f, std::vector<T>& face)
{
    const std::vector<T> ref = cell::reference_vertices<T>(type);
    const int tdim = cell::get_tdim(type);
    face.clear();
    for (const int v : cell::face_vertices(type, f))
        face.insert(face.end(), ref.begin() + v * tdim, ref.begin() + (v + 1) * tdim);
}

template <std::floating_point T>
void append_cell(mesh::CutMesh<T>& out, cell::type type, std::span<const T> x, std::int32_t parent)
{
    const int nv = cell::get_num_vertices(type);
    for (int v = 0; v < nv; ++v)
        out._connectivity.push_back(out._num_vertices + v);
    out._vertex_coords.insert(out._vertex_coords.end(), x.begin(), x.begin() + nv * out._gdim);
    out._num_vertices += nv;
    out._offset.push_back(static_cast<int>(out._connectivity.size()));
    out._types.push_back(type);
    out._parent_map.push_back(parent);
    out._num_cells += 1;
}

} // namespace

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> quadrature_rules(const MeshPart<T, I>& part, int order, bool include_uncut_cells,
                                                const quadrays::Options& options)
{
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
                                         const quadrays::Options& options)
{
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
               const quadrays::Options& options)
{
    quadrays::write_leaves(filename, visualization_mesh(part, degree, include_uncut_cells, options));
}

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> quadrature_rules(const MeshPart<T, I>& part, int order, bool include_uncut_cells,
                                                const lut::Options& options)
{
    check_lut(part, options);
    if (order < 1)
        throw std::invalid_argument("part: the quadrature order must be at least 1");
    const CutResult<T, I>& r = *part.result;
    const MeshView<T, I>& mesh = *r.mesh;
    const int degree = std::min(2 * order - 1, 10);

    quadrature::QuadratureRules<T> rules;
    rules._tdim = r.num_cells > 0 ? cell::get_tdim(mesh.cell_type(I(0))) : 3;
    rules._offset.push_back(0);
    const std::vector<Entry<I>> entries = part_entries(part, include_uncut_cells);
    lut::CellMap<T> map;
    std::vector<I> nodes;
    LutCell<T> lc;
    LevelSetCell<T, I> scratch;
    std::vector<T> face;
    for (std::size_t i = 0; i < entries.size();)
    {
        const I c = entries[i].cell;
        cell_map(mesh, c, map, nodes);
        const std::size_t first = rules._weights.size();
        for (; i < entries.size() && entries[i].cell == c; ++i)
        {
            if (entries[i].kind == 0)
            {
                lut_cell(part, c, map, options, scratch, lc);
                for (int p = 0; p < lc.pieces.n_pieces(); ++p)
                    if (lc.selected[static_cast<std::size_t>(p)])
                        lut::append_piece_rule(map, lc.pieces.types[static_cast<std::size_t>(p)],
                                               piece_vertices(lc.pieces, p), degree, rules._points, rules._weights);
            }
            else if (entries[i].kind == 1)
            {
                const int f = r.zero_face_local[static_cast<std::size_t>(entries[i].index)];
                reference_face(map.type, f, face);
                lut::append_piece_rule(map, cell::face_type(map.type, f), std::span<const T>(face), degree,
                                       rules._points, rules._weights);
            }
            else
            {
                const std::vector<T> ref = cell::reference_vertices<T>(map.type);
                lut::append_piece_rule(map, map.type, std::span<const T>(ref), degree, rules._points, rules._weights);
            }
        }
        if (rules._weights.size() == first)
            continue;
        rules._offset.push_back(static_cast<std::int32_t>(rules._weights.size()));
        rules._parent_map.push_back(static_cast<std::int32_t>(c));
    }
    return rules;
}

template <std::floating_point T, std::integral I>
mesh::CutMesh<T> visualization_mesh(const MeshPart<T, I>& part, bool include_uncut_cells,
                                    const lut::Options& options)
{
    check_lut(part, options);
    const CutResult<T, I>& r = *part.result;
    const MeshView<T, I>& mesh = *r.mesh;
    mesh::CutMesh<T> out;
    out._gdim = mesh.gdim;
    out._tdim = part.dim;
    out._offset.push_back(0);
    lut::CellMap<T> map;
    std::vector<I> nodes;
    LutCell<T> lc;
    LevelSetCell<T, I> scratch;
    std::vector<T> x;
    for (const Entry<I>& e : part_entries(part, include_uncut_cells))
    {
        cell_map(mesh, e.cell, map, nodes);
        const std::int32_t parent = static_cast<std::int32_t>(e.cell);
        if (e.kind == 0)
        {
            lut_cell(part, e.cell, map, options, scratch, lc);
            for (int p = 0; p < lc.pieces.n_pieces(); ++p)
            {
                if (!lc.selected[static_cast<std::size_t>(p)])
                    continue;
                lut::push_forward(map, piece_vertices(lc.pieces, p), x);
                append_cell(out, lc.pieces.types[static_cast<std::size_t>(p)], std::span<const T>(x), parent);
            }
        }
        else if (e.kind == 1)
        {
            const int f = r.zero_face_local[static_cast<std::size_t>(e.index)];
            x.clear();
            for (const int v : cell::face_vertices(map.type, f))
                x.insert(x.end(), map.vertices.begin() + v * map.gdim, map.vertices.begin() + (v + 1) * map.gdim);
            append_cell(out, cell::face_type(map.type, f), std::span<const T>(x), parent);
        }
        else
            append_cell(out, map.type, std::span<const T>(map.vertices), parent);
    }
    return out;
}

template <std::floating_point T, std::integral I>
void write_vtu(const std::string& filename, const MeshPart<T, I>& part, bool include_uncut_cells,
               const lut::Options& options)
{
    io::write_vtk(filename, visualization_mesh(part, include_uncut_cells, options));
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template quadrature::QuadratureRules<float> quadrature_rules<float, int>(const MeshPart<float, int>&, int, bool,
                                                                         const quadrays::Options&);
template quadrature::QuadratureRules<double> quadrature_rules<double, int>(const MeshPart<double, int>&, int, bool,
                                                                           const quadrays::Options&);
template quadrature::QuadratureRules<float> quadrature_rules<float, int>(const MeshPart<float, int>&, int, bool,
                                                                         const lut::Options&);
template quadrature::QuadratureRules<double> quadrature_rules<double, int>(const MeshPart<double, int>&, int, bool,
                                                                           const lut::Options&);
template quadrays::LeafMesh<float> visualization_mesh<float, int>(const MeshPart<float, int>&, int, bool,
                                                                  const quadrays::Options&);
template quadrays::LeafMesh<double> visualization_mesh<double, int>(const MeshPart<double, int>&, int, bool,
                                                                    const quadrays::Options&);
template mesh::CutMesh<float> visualization_mesh<float, int>(const MeshPart<float, int>&, bool, const lut::Options&);
template mesh::CutMesh<double> visualization_mesh<double, int>(const MeshPart<double, int>&, bool,
                                                               const lut::Options&);
template void write_vtu<float, int>(const std::string&, const MeshPart<float, int>&, int, bool,
                                    const quadrays::Options&);
template void write_vtu<double, int>(const std::string&, const MeshPart<double, int>&, int, bool,
                                     const quadrays::Options&);
template void write_vtu<float, int>(const std::string&, const MeshPart<float, int>&, bool, const lut::Options&);
template void write_vtu<double, int>(const std::string&, const MeshPart<double, int>&, bool, const lut::Options&);

} // namespace cutcells::part
