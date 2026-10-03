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
#include "../quadrays/rules.h"
#include "cell_source.h"

namespace cutcells::part
{

namespace
{

/// What a part asks of one cut cell, for quadrays: the level sets that cut
/// the cell and that some term holding on a piece of it names, and those
/// terms with bits over them (local indices), grouped by the zero set they
/// ask for: volume terms first (zero 0), then one group per level set.
/// Conditions on level sets that do not cut the cell hold on all of it.
struct CellTerms
{
    std::vector<int> level_sets;
    std::vector<SelectionTerm> terms;
    std::vector<int> offsets = {0}; ///< groups of terms with the same zero mask
};

template <std::floating_point T, std::integral I>
void cell_terms(const MeshPart<T, I>& part, I cell_id, CellTerms& out)
{
    const CutResult<T, I>& r = *part.result;
    const std::uint64_t cut = cut_mask(r, cell_id);
    out.level_sets.clear();
    out.terms.clear();
    out.offsets.assign(1, 0);
    std::uint64_t named = 0;
    for (const SelectionTerm& term : part.expr.terms)
        if (term_on_cell(term, r, cell_id) == TermCell::piece)
            named |= (term.negative_required | term.positive_required | term.zero_required) & cut;
    std::array<int, 64> local;
    local.fill(-1);
    for (std::uint64_t bits = named; bits != 0; bits &= bits - 1)
    {
        local[static_cast<std::size_t>(std::countr_zero(bits))] = static_cast<int>(out.level_sets.size());
        out.level_sets.push_back(std::countr_zero(bits));
    }
    auto to_local = [&](std::uint64_t mask)
    {
        std::uint64_t m = 0;
        for (std::uint64_t bits = mask & cut; bits != 0; bits &= bits - 1)
            m |= std::uint64_t(1) << local[static_cast<std::size_t>(std::countr_zero(bits))];
        return m;
    };
    for (const SelectionTerm& term : part.expr.terms)
    {
        if (term_on_cell(term, r, cell_id) != TermCell::piece)
            continue;
        SelectionTerm t;
        t.negative_required = to_local(term.negative_required);
        t.positive_required = to_local(term.positive_required);
        t.zero_required = to_local(term.zero_required);
        out.terms.push_back(t);
    }
    std::stable_sort(out.terms.begin(), out.terms.end(), [](const SelectionTerm& a, const SelectionTerm& b)
                     { return a.zero_required < b.zero_required; });
    for (std::size_t i = 1; i <= out.terms.size(); ++i)
        if (i == out.terms.size() || out.terms[i].zero_required != out.terms[i - 1].zero_required)
            out.offsets.push_back(static_cast<int>(i));
}

/// quadrays integrates volumes and interfaces.
template <std::floating_point T, std::integral I>
void check_quadrays(const MeshPart<T, I>& part)
{
    const CutResult<T, I>& r = *part.result;
    if (r.num_cells > 0 && part.dim < cell::get_tdim(r.mesh->cell_type(I(0))) - 1)
        throw std::invalid_argument("part: quadrays integrates volumes and interfaces; the curves where two level "
                                    "sets vanish come from the lookup tables (backend 'lut')");
}

/// The cell's level sets the terms name, as quadrays reads them; throws if
/// quadrays does not take the cell.
template <std::floating_point T, std::integral I>
void sources_of(const MeshPart<T, I>& part, I cell_id, const CellTerms& ct, CellSources<T, I>& cs)
{
    const CutResult<T, I>& r = *part.result;
    thread_local std::vector<const LevelSetFunction<T, I>*> level_sets;
    level_sets.clear();
    for (const int l : ct.level_sets)
        level_sets.push_back(r.level_sets[static_cast<std::size_t>(l)]);
    if (!cell_sources(*r.mesh, std::span<const LevelSetFunction<T, I>* const>(level_sets), cell_id, cs))
    {
        throw std::invalid_argument("part: quadrays takes triangles and quadrilaterals in 2D, tetrahedra, "
                                    "hexahedra, prisms and pyramids in 3D (Pk level sets not on prisms and "
                                    "pyramids)");
    }
}

/// The reference rule of a whole cell, in box coordinates with physical weights.
template <std::floating_point T>
void append_cell_points(const quadrays::ClippedBox<T>& box, cell::type type, int degree,
                        quadrays::CellPoints<T>& out)
{
    const auto rule = quadrature::get_reference_rule<T>(type, degree);
    const T detj = std::abs(quadrays::jacobian_determinant(box));
    const int tdim = rule._tdim;
    for (int q = 0; q < rule._num_points; ++q)
    {
        for (int i = 0; i < 3; ++i)
            out.points.push_back(i < tdim ? rule._points[static_cast<std::size_t>(q * tdim + i)] : T(0));
        out.weights.push_back(rule._weights[static_cast<std::size_t>(q)] * detj);
    }
}

/// The reference rule of facet @p f of a cell (a face in 3D, an edge in 2D),
/// in box coordinates with physical weights.
template <std::floating_point T>
void append_face_points(const quadrays::ClippedBox<T>& box, cell::type type, int f, int degree,
                        quadrays::CellPoints<T>& out)
{
    const std::span<const int> fv = facet_vertices(type, f);
    const bool edge = cell::get_tdim(type) == 2;
    const quadrays::Vec3<T> ua = box_vertex<T>(type, fv[0]), ub = box_vertex<T>(type, fv[1]),
                            uc = box_vertex<T>(type, fv[edge ? 1 : 2]);
    // physical edges of the facet and its measure factor
    std::array<T, 3> e1{}, e2{};
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
        {
            e1[i] += box.jacobian[i][k] * (ub[k] - ua[k]);
            e2[i] += box.jacobian[i][k] * (uc[k] - ua[k]);
        }
    const T measure
        = edge ? std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2])
               : std::sqrt((e1[1] * e2[2] - e1[2] * e2[1]) * (e1[1] * e2[2] - e1[2] * e2[1])
                           + (e1[2] * e2[0] - e1[0] * e2[2]) * (e1[2] * e2[0] - e1[0] * e2[2])
                           + (e1[0] * e2[1] - e1[1] * e2[0]) * (e1[0] * e2[1] - e1[1] * e2[0]));
    const auto rule = quadrature::get_reference_rule<T>(facet_type(type, f), degree);
    for (int q = 0; q < rule._num_points; ++q)
    {
        const T s = rule._points[static_cast<std::size_t>(rule._tdim * q)];
        const T t = edge ? T(0) : rule._points[static_cast<std::size_t>(rule._tdim * q + 1)];
        for (int i = 0; i < 3; ++i)
            out.points.push_back(ua[i] + s * (ub[i] - ua[i]) + t * (uc[i] - ua[i]));
        out.weights.push_back(rule._weights[static_cast<std::size_t>(q)] * measure);
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
    if (r.num_cells > 0 && part.dim < cell::get_tdim(r.mesh->cell_type(I(0))) - 2)
        throw std::invalid_argument("part: the lookup tables give volumes, interfaces and the curves where two "
                                    "level sets vanish, not where three do");
}

/// The map of a cell from its vertices in Basix order.
template <std::floating_point T, std::integral I>
void cell_map(const MeshView<T, I>& mesh, I cell_id, lut::CellMap<T>& map, std::vector<I>& scratch)
{
    map.type = mesh.cell_type(cell_id);
    if (map.type != cell::type::interval && map.type != cell::type::triangle
        && map.type != cell::type::quadrilateral && map.type != cell::type::tetrahedron
        && map.type != cell::type::hexahedron)
    {
        throw std::invalid_argument("part: the lookup tables take intervals, triangles, quadrilaterals, "
                                    "tetrahedra and hexahedra");
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

/// Whether a term holds on piece @p p: the piece lies in exactly the zero sets
/// the term names, and on its sides for the other level sets that cut the
/// cell; the cell's domains decide for the level sets that do not.
template <std::floating_point T, std::integral I>
bool term_on_piece(const SelectionTerm& term, const CutResult<T, I>& r, I cell_id, const LutCell<T>& lc, int p)
{
    const std::uint64_t all = term.negative_required | term.positive_required | term.zero_required;
    std::uint64_t zero = 0;
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
        const std::uint64_t local = std::uint64_t(1) << (it - lc.level_sets.begin());
        if (term.zero_required & bit)
            zero |= local;
        else if (term.negative_required & bit)
        {
            if (!(lc.pieces.negative[static_cast<std::size_t>(p)] & local))
                return false;
        }
        else if (!(lc.pieces.positive[static_cast<std::size_t>(p)] & local))
            return false;
    }
    return lc.pieces.zero[static_cast<std::size_t>(p)] == zero;
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
                     part.dim < tdim ? zero_sets : 0, part.dim == tdim - 2,
                     options.triangulate ? options.triangulation : cell::TriangulationStrategy::none, lc.pieces);
    lc.selected.assign(static_cast<std::size_t>(lc.pieces.n_pieces()), 0);
    for (int p = 0; p < lc.pieces.n_pieces(); ++p)
        for (const SelectionTerm& term : part.expr.terms)
            if (term_on_piece(term, r, cell_id, lc, p))
            {
                lc.selected[static_cast<std::size_t>(p)] = 1;
                break;
            }
}

/// Facet @p f of a cell (a face in 3D, an edge in 2D) in the cell's reference
/// coordinates.
template <std::floating_point T>
void reference_facet(cell::type type, int f, std::vector<T>& facet)
{
    const std::vector<T> ref = cell::reference_vertices<T>(type);
    const int tdim = cell::get_tdim(type);
    facet.clear();
    for (const int v : facet_vertices(type, f))
        facet.insert(facet.end(), ref.begin() + v * tdim, ref.begin() + (v + 1) * tdim);
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
    check_quadrays(part);
    const CutResult<T, I>& r = *part.result;
    const MeshView<T, I>& mesh = *r.mesh;
    const int degree = std::min(2 * order - 1, 10);

    quadrature::QuadratureRules<T> rules;
    rules._tdim = r.num_cells > 0 ? cell::get_tdim(mesh.cell_type(I(0))) : 3;
    rules._offset.push_back(0);
    const std::vector<Entry<I>> entries = part_entries(part, include_uncut_cells);
    CellSources<T, I> cs;
    CellTerms ct;
    quadrays::ClippedBox<T> box;
    quadrays::CellPoints<T> points;
    quadrays::Stats stats;
    for (std::size_t i = 0; i < entries.size();)
    {
        const I c = entries[i].cell;
        const cell::type type = mesh.cell_type(c);
        const int tdim = cell::get_tdim(type);
        if (!quadrays::supported_cell(type) || mesh.gdim != tdim)
        {
            throw std::invalid_argument("part: quadrays takes triangles and quadrilaterals in 2D, tetrahedra, "
                                        "hexahedra, prisms and pyramids in 3D");
        }
        cell_vertex_coords_basix(mesh, c, cs.vertices, cs.nodes);
        quadrays::make_clipped_box<T>(type, std::span<const T>(cs.vertices), mesh.gdim, box);
        points.points.clear();
        points.weights.clear();
        for (; i < entries.size() && entries[i].cell == c; ++i)
        {
            if (entries[i].kind == 0)
            {
                cell_terms(part, c, ct);
                if (ct.terms.empty())
                    continue;
                sources_of(part, c, ct, cs);
                const std::span<const quadrays::Source<T>> sources(cs.sources);
                for (std::size_t g = 0; g + 1 < ct.offsets.size(); ++g)
                    quadrays::integrate(cs.box, sources,
                                        std::span<const SelectionTerm>(ct.terms).subspan(
                                            static_cast<std::size_t>(ct.offsets[g]),
                                            static_cast<std::size_t>(ct.offsets[g + 1] - ct.offsets[g])),
                                        order, options, points, stats);
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
            rules._points.insert(rules._points.end(), xi.begin(), xi.begin() + tdim);
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
    check_quadrays(part);
    const CutResult<T, I>& r = *part.result;
    const MeshView<T, I>& mesh = *r.mesh;
    quadrays::LeafMesh<T> leaves;
    CellSources<T, I> cs;
    CellTerms ct;
    quadrays::Stats stats;
    std::vector<T> face;
    for (const Entry<I>& e : part_entries(part, include_uncut_cells))
    {
        const cell::type type = mesh.cell_type(e.cell);
        if (e.kind == 0)
        {
            cell_terms(part, e.cell, ct);
            if (ct.terms.empty())
                continue;
            sources_of(part, e.cell, ct, cs);
            for (std::size_t g = 0; g + 1 < ct.offsets.size(); ++g)
                quadrays::append_leaves(cs.box, std::span<const quadrays::Source<T>>(cs.sources),
                                        std::span<const SelectionTerm>(ct.terms).subspan(
                                            static_cast<std::size_t>(ct.offsets[g]),
                                            static_cast<std::size_t>(ct.offsets[g + 1] - ct.offsets[g])),
                                        degree, options, static_cast<std::int32_t>(e.cell), leaves, stats);
            continue;
        }
        cell_vertex_coords_basix(mesh, e.cell, cs.vertices, cs.nodes);
        if (e.kind == 1)
        {
            const int f = r.zero_face_local[static_cast<std::size_t>(e.index)];
            face.clear();
            for (const int v : facet_vertices(type, f))
                face.insert(face.end(), cs.vertices.begin() + mesh.gdim * v,
                            cs.vertices.begin() + mesh.gdim * (v + 1));
            quadrays::append_linear_cell<T>(facet_type(type, f), std::span<const T>(face), mesh.gdim,
                                            static_cast<std::int32_t>(e.cell), leaves);
        }
        else
            quadrays::append_linear_cell<T>(type, std::span<const T>(cs.vertices), mesh.gdim,
                                            static_cast<std::int32_t>(e.cell), leaves);
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
                reference_facet(map.type, f, face);
                lut::append_piece_rule(map, facet_type(map.type, f), std::span<const T>(face), degree, rules._points,
                                       rules._weights);
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
            for (const int v : facet_vertices(map.type, f))
                x.insert(x.end(), map.vertices.begin() + v * map.gdim, map.vertices.begin() + (v + 1) * map.gdim);
            append_cell(out, facet_type(map.type, f), std::span<const T>(x), parent);
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
