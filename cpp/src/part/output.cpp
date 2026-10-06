// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "output.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>
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
/// ask for and the signs of its normals: volume terms first (zero 0), then
/// one group per level set and sign, so that each run of the engine lies on
/// one zero set with one normal. Conditions on level sets that do not cut
/// the cell hold on all of it.
struct CellTerms
{
    std::vector<int> level_sets;
    std::vector<SelectionTerm> terms;
    std::vector<int> offsets = {0}; ///< groups of terms with the same zero and normal masks
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
        t.normal_flipped = to_local(term.normal_flipped);
        out.terms.push_back(t);
    }
    auto group = [](const SelectionTerm& t) { return std::make_pair(t.zero_required, t.normal_flipped); };
    std::stable_sort(out.terms.begin(), out.terms.end(),
                     [&group](const SelectionTerm& a, const SelectionTerm& b) { return group(a) < group(b); });
    for (std::size_t i = 1; i <= out.terms.size(); ++i)
        if (i == out.terms.size() || group(out.terms[i]) != group(out.terms[i - 1]))
            out.offsets.push_back(static_cast<int>(i));
}

/// quadrays integrates volumes, interfaces and the curves where two level sets
/// vanish (points on 2D cells).
template <std::floating_point T, std::integral I>
void check_quadrays(const MeshPart<T, I>& part)
{
    const CutResult<T, I>& r = *part.result;
    if (r.num_cells > 0 && part.dim < cell::get_tdim(r.mesh->cell_type(I(0))) - 2)
        throw std::invalid_argument("part: quadrays integrates volumes, interfaces and the curves where two level "
                                    "sets vanish, not the points where three do");
}

/// The cells quadrays takes, for errors.
constexpr const char* quadrays_cells
    = "part: quadrays takes intervals, triangles, quadrilaterals, tetrahedra, hexahedra, prisms and pyramids, in "
      "a space of their dimension or a higher one up to 3D (facets, surface meshes)";

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
        throw std::invalid_argument(quadrays_cells);
}

/// The reference rule of a cell type exact for @p degree: the lowest at least
/// as high with positive weights, made once per thread.
template <std::floating_point T>
const quadrature::ReferenceQuadratureRule<T>& reference_rule(cell::type type, int degree)
{
    thread_local std::array<std::array<quadrature::ReferenceQuadratureRule<T>, 11>, 8> cache;
    quadrature::ReferenceQuadratureRule<T>& rule = cache[static_cast<int>(type)][degree];
    if (rule._weights.empty())
        rule = quadrature::get_reference_rule<T>(type, quadrature::positive_rule_order(type, degree));
    return rule;
}

/// The reference rule of a whole cell, in box coordinates with physical weights.
template <std::floating_point T>
void append_cell_points(const quadrays::ClippedBox<T>& box, cell::type type, int degree,
                        quadrays::CellPoints<T>& out)
{
    const auto& rule = reference_rule<T>(type, degree);
    const T detj = std::abs(quadrays::reference_jacobian_determinant(box));
    const int tdim = rule._tdim;
    for (int q = 0; q < rule._num_points; ++q)
    {
        quadrays::Vec3<T> xi = {0, 0, 0};
        for (int i = 0; i < tdim; ++i)
            xi[i] = rule._points[static_cast<std::size_t>(q * tdim + i)];
        const quadrays::Vec3<T> u = quadrays::box_point(box, xi);
        out.points.insert(out.points.end(), u.begin(), u.end());
        out.weights.push_back(rule._weights[static_cast<std::size_t>(q)] * detj);
    }
}

/// What a part has in one cell: a piece the backend integrates, or the whole cell.
template <std::integral I>
struct Entry
{
    I cell;
    bool whole;
};

template <std::floating_point T, std::integral I>
std::vector<Entry<I>> part_entries(const MeshPart<T, I>& part, bool include_uncut_cells)
{
    std::vector<Entry<I>> entries;
    for (const I c : part.cut_cells)
        entries.push_back({c, false});
    if (include_uncut_cells)
        for (const I c : part.uncut_cells)
            entries.push_back({c, true});
    std::stable_sort(entries.begin(), entries.end(),
                     [](const Entry<I>& a, const Entry<I>& b) { return a.cell < b.cell; });
    return entries;
}

// ============================================================================
// Where points and pieces lie in their cell
// ============================================================================

/// Points within this distance of a facet, edge or vertex of their cell, in
/// reference coordinates, lie on it: twice the Newton distance within which
/// quadrays counts a root at an end of a line (engine.cpp, point_tol).
template <std::floating_point T>
constexpr T boundary_tol = T(2) * std::max(static_cast<T>(1e-10), T(100) * std::numeric_limits<T>::epsilon());

/// The facets of a reference cell as planes n . x = c with unit n (n has tdim
/// entries), made once per thread.
template <std::floating_point T>
const std::vector<std::array<T, 4>>& facet_planes(cell::type type)
{
    thread_local std::array<std::vector<std::array<T, 4>>, 8> cache;
    std::vector<std::array<T, 4>>& planes = cache[static_cast<int>(type)];
    if (!planes.empty())
        return planes;
    const int tdim = cell::get_tdim(type);
    const std::vector<T> ref = cell::reference_vertices<T>(type);
    for (int f = 0; f < cell::num_sub_entities(type, tdim - 1); ++f)
    {
        const std::span<const int> fv = cell::sub_entity_vertices(type, tdim - 1, f);
        const T* a = ref.data() + fv[0] * tdim;
        std::array<T, 4> p{};
        if (tdim == 1)
            p[0] = T(1);
        else if (tdim == 2)
        {
            const T* b = ref.data() + fv[1] * tdim;
            p[0] = -(b[1] - a[1]);
            p[1] = b[0] - a[0];
        }
        else
        {
            const T* b = ref.data() + fv[1] * tdim;
            const T* c = ref.data() + fv[2] * tdim;
            const T u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
            p[0] = u[1] * v[2] - u[2] * v[1];
            p[1] = u[2] * v[0] - u[0] * v[2];
            p[2] = u[0] * v[1] - u[1] * v[0];
        }
        T norm = T(0);
        for (int d = 0; d < tdim; ++d)
            norm += p[d] * p[d];
        norm = std::sqrt(norm);
        for (int d = 0; d < tdim; ++d)
        {
            p[d] /= norm;
            p[3] += p[d] * a[d];
        }
        planes.push_back(p);
    }
    return planes;
}

/// The facets of a cell (bits) whose planes hold the reference point @p xi.
template <std::floating_point T>
unsigned facets_holding(cell::type type, const T* xi, T tol)
{
    const int tdim = cell::get_tdim(type);
    unsigned bits = 0;
    const std::vector<std::array<T, 4>>& planes = facet_planes<T>(type);
    for (std::size_t f = 0; f < planes.size(); ++f)
    {
        T s = -planes[f][3];
        for (int d = 0; d < tdim; ++d)
            s += planes[f][d] * xi[d];
        if (std::abs(s) <= tol)
            bits |= 1u << f;
    }
    return bits;
}

/// The sub-entity of a cell that the facets @p bits (non-empty) share, as
/// (dimension, local number), or (-1, -1) if they share no vertex.
inline std::pair<int, int> shared_entity(cell::type type, unsigned bits)
{
    const int tdim = cell::get_tdim(type);
    unsigned vertices = (1u << cell::get_num_vertices(type)) - 1;
    for (unsigned b = bits; b != 0; b &= b - 1)
    {
        unsigned fv = 0;
        for (const int v : cell::sub_entity_vertices(type, tdim - 1, std::countr_zero(b)))
            fv |= 1u << v;
        vertices &= fv;
    }
    for (int dim = tdim - 1; dim >= 0; --dim)
        for (int i = 0; i < cell::num_sub_entities(type, dim); ++i)
        {
            unsigned ev = 0;
            for (const int v : cell::sub_entity_vertices(type, dim, i))
                ev |= 1u << v;
            if (ev == vertices)
                return {dim, i};
        }
    return {-1, -1};
}

/// The key of sub-entity @p i of dimension @p dim of a cell with vertex nodes @p nodes.
template <std::integral I>
std::array<I, 4> sub_entity_key(cell::type type, int dim, int i, std::span<const I> nodes)
{
    std::array<I, 4> local{};
    const std::span<const int> sv = cell::sub_entity_vertices(type, dim, i);
    for (std::size_t k = 0; k < sv.size(); ++k)
        local[k] = nodes[static_cast<std::size_t>(sv[k])];
    return entity_key(std::span<const I>(local.data(), sv.size()));
}

// ============================================================================
// Rules on entities
// ============================================================================

/// The entity rules of a part, by dimension, with the whole entities of its
/// own dimension.
template <std::floating_point T, std::integral I>
std::vector<EntityRules<T>> whole_entities(const MeshPart<T, I>& part, int order)
{
    const int tdim = part.result->tdim();
    std::vector<EntityRules<T>> rules(static_cast<std::size_t>(std::max(tdim, 0)));
    for (int d = 0; d < tdim; ++d)
    {
        rules[static_cast<std::size_t>(d)].tdim = tdim;
        rules[static_cast<std::size_t>(d)].entity_dim = d;
    }
    if (part.dim >= 0 && part.dim < tdim)
        rules[static_cast<std::size_t>(part.dim)] = whole_entity_rules(part, order);
    return rules;
}

/// A part of dimension 0 on 2D cells, the points where two curves cross: the
/// points of the cell rules @p rules within boundary_tol of an edge or a
/// vertex of their cell move to rules on that entity, one point per position
/// (weight 1, all sharing cells); the vertices the part holds whole keep their
/// own rule. The entities on which no level set vanishes, which the result
/// does not list, get their sharing cells here. Rebuilds @p entities[0] and [1].
template <std::floating_point T, std::integral I>
void points_to_entities(const MeshPart<T, I>& part, quadrature::QuadratureRules<T>& rules,
                        std::vector<EntityRules<T>>& entities)
{
    const CutResult<T, I>& r = *part.result;
    const MeshView<T, I>& mesh = *r.mesh;
    const T tol = boundary_tol<T>;
    struct Found
    {
        int dim;
        std::array<I, 4> key;
        T t; ///< along an edge's chart
        auto order() const { return std::tie(dim, key, t); }
    };
    std::vector<Found> found;
    quadrature::QuadratureRules<T> kept;
    kept._tdim = rules._tdim;
    kept._offset.push_back(0);
    std::vector<I> nodes, scratch;
    for (std::size_t k = 0; k + 1 < rules._offset.size(); ++k)
    {
        const I c = static_cast<I>(rules._parent_map[k]);
        const cell::type type = mesh.cell_type(c);
        const std::vector<T> ref = cell::reference_vertices<T>(type);
        cell_vertex_nodes(mesh, c, nodes, scratch);
        for (std::int32_t p = rules._offset[k]; p < rules._offset[k + 1]; ++p)
        {
            const T* xi = rules._points.data() + 2 * p;
            const unsigned bits = facets_holding<T>(type, xi, tol);
            const auto [dim, local] = bits != 0 ? shared_entity(type, bits) : std::pair<int, int>{-1, -1};
            if (dim < 0)
            {
                kept._points.insert(kept._points.end(), xi, xi + 2);
                kept._weights.push_back(rules._weights[static_cast<std::size_t>(p)]);
                continue;
            }
            const std::array<I, 4> key = sub_entity_key(type, dim, local, std::span<const I>(nodes));
            T t = T(0);
            if (dim == 1)
            {
                // the parameter along the edge from its first vertex in the
                // cell, then in its chart, which starts at its lower node
                const std::span<const int> sv = cell::sub_entity_vertices(type, 1, local);
                const T* a = ref.data() + 2 * sv[0];
                const T* b = ref.data() + 2 * sv[1];
                const T s = ((xi[0] - a[0]) * (b[0] - a[0]) + (xi[1] - a[1]) * (b[1] - a[1]))
                            / ((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]));
                t = std::clamp(nodes[static_cast<std::size_t>(sv[0])] == key[0] ? s : T(1) - s, T(0), T(1));
            }
            found.push_back({dim, key, t});
        }
        if (static_cast<std::size_t>(kept._offset.back()) == kept._weights.size())
            continue;
        kept._offset.push_back(static_cast<std::int32_t>(kept._weights.size()));
        kept._parent_map.push_back(rules._parent_map[k]);
    }
    rules = std::move(kept);
    std::sort(found.begin(), found.end(), [](const Found& a, const Found& b) { return a.order() < b.order(); });

    // the entities the points lie on: the result's, or tables of those it does not list
    std::array<EntityTable<I>, 2> others;
    {
        std::array<std::vector<std::array<I, 4>>, 2> keys;
        std::vector<I> star_nodes;
        for (const Found& f : found)
            if (find_entity(r.entities[static_cast<std::size_t>(f.dim)], f.key) < 0)
            {
                keys[static_cast<std::size_t>(f.dim)].push_back(f.key);
                for (int k = 0; k <= f.dim; ++k)
                    star_nodes.push_back(f.key[static_cast<std::size_t>(k)]);
            }
        std::sort(star_nodes.begin(), star_nodes.end());
        star_nodes.erase(std::unique(star_nodes.begin(), star_nodes.end()), star_nodes.end());
        NodeStars<I> stars;
        if (!star_nodes.empty())
            node_stars(mesh, star_nodes, std::span<const I>(star_nodes), stars);
        for (int dim = 0; dim < 2; ++dim)
        {
            std::vector<std::array<I, 4>>& k = keys[static_cast<std::size_t>(dim)];
            std::sort(k.begin(), k.end());
            k.erase(std::unique(k.begin(), k.end()), k.end());
            others[static_cast<std::size_t>(dim)] = make_entity_table(mesh, dim, std::move(k), stars);
        }
    }
    auto entity = [&](int dim, const std::array<I, 4>& key) -> std::pair<const EntityTable<I>*, std::int32_t>
    {
        const std::int32_t e = find_entity(r.entities[static_cast<std::size_t>(dim)], key);
        if (e >= 0)
            return {&r.entities[static_cast<std::size_t>(dim)], e};
        return {&others[static_cast<std::size_t>(dim)], find_entity(others[static_cast<std::size_t>(dim)], key)};
    };
    std::uint64_t zero = 0;
    for (const SelectionTerm& term : part.expr.terms)
        zero |= term.zero_required;

    // vertices: those the part holds whole, and those cells found a point at
    std::vector<std::pair<std::array<I, 4>, std::uint64_t>> vertices;
    for (std::size_t k = 0; k < part.entities[0].size(); ++k)
        vertices.emplace_back(r.entities[0].keys[static_cast<std::size_t>(part.entities[0][k])], part.entity_zero[0][k]);
    for (const Found& f : found)
        if (f.dim == 0)
            vertices.emplace_back(f.key, zero);
    std::stable_sort(vertices.begin(), vertices.end(),
                     [](const auto& a, const auto& b) { return a.first < b.first; });
    EntityRules<T>& points = entities[0];
    points = EntityRules<T>{};
    points.tdim = 2;
    points.entity_dim = 0;
    for (std::size_t k = 0; k < vertices.size(); ++k)
    {
        if (k > 0 && vertices[k].first == vertices[k - 1].first)
            continue;
        const auto [table, e] = entity(0, vertices[k].first);
        points.entity_types.push_back(cell::type::point);
        points.entity_nodes.push_back(static_cast<std::int32_t>(vertices[k].first[0]));
        points.node_offsets.push_back(static_cast<std::int32_t>(points.entity_nodes.size()));
        points.zero.push_back(vertices[k].second);
        points.weights.push_back(T(1));
        points.offsets.push_back(static_cast<std::int32_t>(points.weights.size()));
        append_sharing_cells(mesh, *table, e, points);
    }

    // edges: the points cells found on them, one per position
    EntityRules<T>& edges = entities[1];
    edges = EntityRules<T>{};
    edges.tdim = 2;
    edges.entity_dim = 1;
    for (std::size_t i = 0; i < found.size();)
    {
        if (found[i].dim != 1)
        {
            ++i;
            continue;
        }
        const std::array<I, 4> key = found[i].key;
        std::size_t j = i;
        T last = -std::numeric_limits<T>::infinity();
        for (; j < found.size() && found[j].dim == 1 && found[j].key == key; ++j)
        {
            if (found[j].t - last <= tol)
                continue;
            last = found[j].t;
            edges.points.push_back(last);
            edges.weights.push_back(T(1));
        }
        const auto [table, e] = entity(1, key);
        edges.entity_types.push_back(cell::type::interval);
        for (const I n : table->entity_nodes(e))
            edges.entity_nodes.push_back(static_cast<std::int32_t>(n));
        edges.node_offsets.push_back(static_cast<std::int32_t>(edges.entity_nodes.size()));
        edges.zero.push_back(zero);
        edges.offsets.push_back(static_cast<std::int32_t>(edges.weights.size()));
        append_sharing_cells(mesh, *table, e, edges);
        i = j;
    }
}

/// Cell rules with the entity rules taken in as @p policy says: one rule per
/// cell, cells ascending; each cell's own points first, then the entities'.
template <std::floating_point T>
quadrature::QuadratureRules<T> fold(const quadrature::QuadratureRules<T>& cells,
                                    const std::vector<EntityRules<T>>& entities, EntityPolicy policy)
{
    struct Contribution
    {
        std::int32_t cell;
        int dim;
        int rule;
        std::int32_t slot;
        T factor;
    };
    std::vector<Contribution> contributions;
    std::vector<std::int32_t> slots;
    for (std::size_t d = 0; d < entities.size(); ++d)
        for (int r = 0; r < entities[d].n_rules(); ++r)
        {
            T factor;
            policy_cells(entities[d], r, policy, slots, factor);
            for (const std::int32_t s : slots)
                contributions.push_back(
                    {entities[d].cells[static_cast<std::size_t>(s)], static_cast<int>(d), r, s, factor});
        }
    std::stable_sort(contributions.begin(), contributions.end(),
                     [](const Contribution& a, const Contribution& b) { return a.cell < b.cell; });
    quadrature::QuadratureRules<T> out;
    out._tdim = cells._tdim;
    out._offset.push_back(0);
    const std::size_t tdim = static_cast<std::size_t>(cells._tdim);
    const std::size_t n = cells._parent_map.size();
    std::size_t i = 0, j = 0;
    while (i < n || j < contributions.size())
    {
        std::int32_t c = std::numeric_limits<std::int32_t>::max();
        if (i < n)
            c = cells._parent_map[i];
        if (j < contributions.size())
            c = std::min(c, contributions[j].cell);
        if (i < n && cells._parent_map[i] == c)
        {
            out._points.insert(out._points.end(),
                               cells._points.begin() + static_cast<std::ptrdiff_t>(cells._offset[i] * tdim),
                               cells._points.begin() + static_cast<std::ptrdiff_t>(cells._offset[i + 1] * tdim));
            out._weights.insert(out._weights.end(), cells._weights.begin() + cells._offset[i],
                                cells._weights.begin() + cells._offset[i + 1]);
            ++i;
        }
        for (; j < contributions.size() && contributions[j].cell == c; ++j)
        {
            const Contribution& k = contributions[j];
            append_cell_contribution(entities[static_cast<std::size_t>(k.dim)], k.rule, k.slot, k.factor, out._points,
                                     out._weights);
        }
        if (static_cast<std::size_t>(out._offset.back()) == out._weights.size())
            continue;
        out._offset.push_back(static_cast<std::int32_t>(out._weights.size()));
        out._parent_map.push_back(c);
    }
    return out;
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
        && map.type != cell::type::hexahedron && map.type != cell::type::prism && map.type != cell::type::pyramid)
    {
        throw std::invalid_argument("part: the lookup tables take intervals, triangles, quadrilaterals, "
                                    "tetrahedra, hexahedra, prisms and pyramids");
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
/// the highest degree of the level sets, 2 for analytic ones, at most 4, on
/// pyramids 2), the values of
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
    // pyramids have templates of orders 1 and 2
    const int k = options.template_order > 0 ? options.template_order
                                             : std::min(degree, map.type == cell::type::pyramid ? 2 : 4);
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

/// Whether piece @p p of a cut cell is an entity's: it lies in a facet or an
/// edge of the cell (or at a vertex) on which the zero sets it lies in all
/// vanish. The entity's rule then has it, or it is no surface (the zero set
/// does not change sign there). Cells keep to their interior.
template <std::floating_point T, std::integral I>
bool on_zero_entity(const CutResult<T, I>& r, cell::type type, std::span<const I> nodes, const LutCell<T>& lc, int p)
{
    std::uint64_t zero = 0;
    for (std::uint64_t bits = lc.pieces.zero[static_cast<std::size_t>(p)]; bits != 0; bits &= bits - 1)
        zero |= std::uint64_t(1) << lc.level_sets[static_cast<std::size_t>(std::countr_zero(bits))];
    if (zero == 0)
        return false;
    const std::span<const T> x = piece_vertices(lc.pieces, p);
    const std::size_t tdim = static_cast<std::size_t>(cell::get_tdim(type));
    unsigned bits = ~0u;
    for (std::size_t v = 0; v < x.size() / tdim && bits != 0; ++v)
        bits &= facets_holding<T>(type, x.data() + v * tdim, boundary_tol<T>);
    if (bits == 0)
        return false;
    const auto [dim, local] = shared_entity(type, bits);
    if (dim < 0)
        return false;
    const EntityTable<I>& table = r.entities[static_cast<std::size_t>(dim)];
    const std::int32_t e = find_entity(table, sub_entity_key(type, dim, local, nodes));
    return e >= 0 && (zero & ~table.vanishing[static_cast<std::size_t>(e)]) == 0;
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

/// The polynomial degree a quadrature order integrates exactly on flat pieces,
/// whole cells and entities: the order itself, within the reference tables.
int exact_degree(int order)
{
    if (order < 1 || order > 10)
        throw std::invalid_argument("part: the quadrature order, the polynomial degree integrated exactly, goes "
                                    "from 1 to 10");
    return order;
}

/// quadrays' Gauss-Legendre points per segment for a quadrature order. On a
/// flat piece the m nested levels of a part of dimension m have linear bounds,
/// each of which raises the degree of the integrand by one: degree order needs
/// 2 q - 1 >= order + m - 1.
int points_per_segment(int order, int dim, const quadrays::Options& options)
{
    if (options.points_per_segment > 0)
        return options.points_per_segment;
    return std::max(1, (order + std::max(dim, 0) + 1) / 2);
}

// ============================================================================
// Cell rules and entity rules of a part
// ============================================================================

/// The rules of a part before the cells take in the entities' (quadrays): the
/// cells' own, inside them, and the entities' by dimension.
template <std::floating_point T, std::integral I>
void part_rules(const MeshPart<T, I>& part, int order, bool include_uncut_cells, const quadrays::Options& options,
                quadrature::QuadratureRules<T>& rules, std::vector<EntityRules<T>>& entities)
{
    const int degree = exact_degree(order);
    const int q = points_per_segment(order, part.dim, options);
    check_quadrays(part);
    const CutResult<T, I>& r = *part.result;
    const MeshView<T, I>& mesh = *r.mesh;

    rules = quadrature::QuadratureRules<T>{};
    rules._tdim = r.num_cells > 0 ? r.tdim() : 3;
    rules._offset.push_back(0);
    entities = whole_entities(part, degree);
    // a curve in 3D lying in a face or an edge of a cell is that entity's
    const bool curves = part.dim == 1 && rules._tdim == 3;
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
        if (!quadrays::supported_cell(type) || mesh.gdim < tdim || mesh.gdim > 3)
            throw std::invalid_argument(quadrays_cells);
        cell_vertex_coords_basix(mesh, c, cs.vertices, cs.nodes);
        // a cut cell's points come in the frame of its sources; whole cells
        // take the reference box, which costs nothing to build
        bool cut = false;
        for (std::size_t j = i; j < entries.size() && entries[j].cell == c; ++j)
            cut = cut || !entries[j].whole;
        quadrays::make_clipped_box<T>(type, std::span<const T>(cs.vertices), mesh.gdim, box,
                                      cut ? quadrays::BoxFrame::orthogonal : quadrays::BoxFrame::reference);
        points.points.clear();
        points.weights.clear();
        for (; i < entries.size() && entries[i].cell == c; ++i)
        {
            if (!entries[i].whole)
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
                                        q, options, points, stats);
            }
            else
                append_cell_points(box, type, degree, points);
        }
        const std::size_t first = rules._weights.size();
        for (int k = 0; k < points.n_points(); ++k)
        {
            const quadrays::Vec3<T> u = {points.points[3 * k], points.points[3 * k + 1], points.points[3 * k + 2]};
            const quadrays::Vec3<T> xi = quadrays::reference_point(box, u);
            if (curves && facets_holding<T>(type, xi.data(), boundary_tol<T>) != 0)
                continue;
            rules._points.insert(rules._points.end(), xi.begin(), xi.begin() + tdim);
            rules._weights.push_back(points.weights[k]);
        }
        if (rules._weights.size() == first)
            continue;
        rules._offset.push_back(static_cast<std::int32_t>(rules._weights.size()));
        rules._parent_map.push_back(static_cast<std::int32_t>(c));
    }
    if (part.dim == 0 && rules._tdim == 2)
        points_to_entities(part, rules, entities);
}

/// The rules of a part before the cells take in the entities' (lookup tables).
template <std::floating_point T, std::integral I>
void part_rules(const MeshPart<T, I>& part, int order, bool include_uncut_cells, const lut::Options& options,
                quadrature::QuadratureRules<T>& rules, std::vector<EntityRules<T>>& entities)
{
    check_lut(part, options);
    const int degree = exact_degree(order);
    const CutResult<T, I>& r = *part.result;
    const MeshView<T, I>& mesh = *r.mesh;

    rules = quadrature::QuadratureRules<T>{};
    rules._tdim = r.num_cells > 0 ? r.tdim() : 3;
    rules._offset.push_back(0);
    entities = whole_entities(part, degree);
    const std::vector<Entry<I>> entries = part_entries(part, include_uncut_cells);
    lut::CellMap<T> map;
    std::vector<I> scratch, nodes;
    LutCell<T> lc;
    LevelSetCell<T, I> ls_scratch;
    for (std::size_t i = 0; i < entries.size();)
    {
        const I c = entries[i].cell;
        cell_map(mesh, c, map, scratch);
        cell_vertex_nodes(mesh, c, nodes, scratch);
        const std::size_t first = rules._weights.size();
        for (; i < entries.size() && entries[i].cell == c; ++i)
        {
            if (!entries[i].whole)
            {
                lut_cell(part, c, map, options, ls_scratch, lc);
                for (int p = 0; p < lc.pieces.n_pieces(); ++p)
                {
                    if (!lc.selected[static_cast<std::size_t>(p)]
                        || (part.dim > 0 && on_zero_entity(r, map.type, std::span<const I>(nodes), lc, p)))
                        continue;
                    lut::append_piece_rule(map, lc.pieces.types[static_cast<std::size_t>(p)],
                                           piece_vertices(lc.pieces, p), degree, rules._points, rules._weights);
                }
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
    if (part.dim == 0 && rules._tdim == 2)
        points_to_entities(part, rules, entities);
}

/// The points of a part of dimension 0 on 2D cells, each once, for
/// visualisation: inside cells with their cell as parent, on entities with
/// their lowest sharing cell.
template <std::floating_point T, std::integral I, typename Emit>
void emit_points(const MeshPart<T, I>& part, const quadrature::QuadratureRules<T>& rules,
                 const std::vector<EntityRules<T>>& entities, Emit&& emit)
{
    const MeshView<T, I>& mesh = *part.result->mesh;
    lut::CellMap<T> map;
    std::vector<I> scratch;
    std::vector<T> x;
    for (std::size_t k = 0; k + 1 < rules._offset.size(); ++k)
    {
        const I c = static_cast<I>(rules._parent_map[k]);
        cell_map(mesh, c, map, scratch);
        lut::push_forward(map,
                          std::span<const T>(rules._points)
                              .subspan(static_cast<std::size_t>(rules._offset[k]) * 2,
                                       static_cast<std::size_t>(rules._offset[k + 1] - rules._offset[k]) * 2),
                          x);
        for (std::size_t p = 0; p < x.size() / static_cast<std::size_t>(mesh.gdim); ++p)
            emit(std::span<const T>(x).subspan(p * static_cast<std::size_t>(mesh.gdim),
                                              static_cast<std::size_t>(mesh.gdim)),
                 static_cast<std::int32_t>(c));
    }
    for (const EntityRules<T>& rules_d : entities)
        for (int r = 0; r < rules_d.n_rules(); ++r)
        {
            const std::int32_t s = rules_d.cell_offsets[r];
            const I c = static_cast<I>(rules_d.cells[static_cast<std::size_t>(s)]);
            cell_map(mesh, c, map, scratch);
            lut::push_forward(map,
                              std::span<const T>(rules_d.cell_points)
                                  .subspan(static_cast<std::size_t>(rules_d.cell_point_offsets[s]) * 2,
                                           static_cast<std::size_t>(rules_d.cell_point_offsets[s + 1]
                                                                    - rules_d.cell_point_offsets[s])
                                               * 2),
                              x);
            for (std::size_t p = 0; p < x.size() / static_cast<std::size_t>(mesh.gdim); ++p)
                emit(std::span<const T>(x).subspan(p * static_cast<std::size_t>(mesh.gdim),
                                                  static_cast<std::size_t>(mesh.gdim)),
                     static_cast<std::int32_t>(c));
        }
}

/// The entities a part holds whole, once each, for visualisation: their type,
/// physical vertices and lowest sharing cell.
template <std::floating_point T, std::integral I, typename Emit>
void emit_entities(const MeshPart<T, I>& part, Emit&& emit)
{
    const CutResult<T, I>& r = *part.result;
    const MeshView<T, I>& mesh = *r.mesh;
    if (part.dim < 0 || part.dim >= r.tdim())
        return;
    const EntityTable<I>& table = r.entities[static_cast<std::size_t>(part.dim)];
    std::vector<T> x;
    for (const std::int32_t e : part.entities[static_cast<std::size_t>(part.dim)])
    {
        x.clear();
        for (const I n : table.entity_nodes(e))
            x.insert(x.end(), mesh.node(n), mesh.node(n) + mesh.gdim);
        emit(table.types[static_cast<std::size_t>(e)], std::span<const T>(x),
             static_cast<std::int32_t>(table.entity_cells(e)[0]));
    }
}

} // namespace

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> quadrature_rules(const MeshPart<T, I>& part, int order, bool include_uncut_cells,
                                                const quadrays::Options& options, EntityPolicy entities)
{
    quadrature::QuadratureRules<T> cells;
    std::vector<EntityRules<T>> on_entities;
    part_rules(part, order, include_uncut_cells, options, cells, on_entities);
    return fold(cells, on_entities, entities);
}

template <std::floating_point T, std::integral I>
quadrature::QuadratureRules<T> quadrature_rules(const MeshPart<T, I>& part, int order, bool include_uncut_cells,
                                                const lut::Options& options, EntityPolicy entities)
{
    quadrature::QuadratureRules<T> cells;
    std::vector<EntityRules<T>> on_entities;
    part_rules(part, order, include_uncut_cells, options, cells, on_entities);
    return fold(cells, on_entities, entities);
}

template <std::floating_point T, std::integral I>
std::vector<EntityRules<T>> entity_rules(const MeshPart<T, I>& part, int order, const quadrays::Options& options)
{
    // only points on 2D cells need the backend to find them
    if (!(part.dim == 0 && part.result->tdim() == 2))
    {
        check_quadrays(part);
        return whole_entities(part, exact_degree(order));
    }
    quadrature::QuadratureRules<T> cells;
    std::vector<EntityRules<T>> on_entities;
    part_rules(part, order, false, options, cells, on_entities);
    return on_entities;
}

template <std::floating_point T, std::integral I>
std::vector<EntityRules<T>> entity_rules(const MeshPart<T, I>& part, int order, const lut::Options& options)
{
    if (!(part.dim == 0 && part.result->tdim() == 2))
    {
        check_lut(part, options);
        return whole_entities(part, exact_degree(order));
    }
    quadrature::QuadratureRules<T> cells;
    std::vector<EntityRules<T>> on_entities;
    part_rules(part, order, false, options, cells, on_entities);
    return on_entities;
}

template <std::floating_point T, std::integral I>
quadrays::LeafMesh<T> visualization_mesh(const MeshPart<T, I>& part, int degree, bool include_uncut_cells,
                                         const quadrays::Options& options)
{
    check_quadrays(part);
    const CutResult<T, I>& r = *part.result;
    const MeshView<T, I>& mesh = *r.mesh;
    quadrays::LeafMesh<T> leaves;
    if (part.dim == 0 && r.tdim() == 2)
    {
        // points, each once
        quadrature::QuadratureRules<T> cells;
        std::vector<EntityRules<T>> entities;
        part_rules(part, 1, false, options, cells, entities);
        emit_points(part, cells, entities, [&](std::span<const T> x, std::int32_t parent)
                    { quadrays::append_linear_cell<T>(cell::type::point, x, mesh.gdim, parent, leaves); });
        return leaves;
    }
    CellSources<T, I> cs;
    CellTerms ct;
    quadrays::Stats stats;
    for (const Entry<I>& e : part_entries(part, include_uncut_cells))
    {
        const cell::type type = mesh.cell_type(e.cell);
        if (!e.whole)
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
        quadrays::append_linear_cell<T>(type, std::span<const T>(cs.vertices), mesh.gdim,
                                        static_cast<std::int32_t>(e.cell), leaves);
    }
    emit_entities(part, [&](cell::type type, std::span<const T> x, std::int32_t parent)
                  { quadrays::append_linear_cell<T>(type, x, mesh.gdim, parent, leaves); });
    return leaves;
}

template <std::floating_point T, std::integral I>
void write_vtu(const std::string& filename, const MeshPart<T, I>& part, int degree, bool include_uncut_cells,
               const quadrays::Options& options)
{
    quadrays::write_leaves(filename, visualization_mesh(part, degree, include_uncut_cells, options));
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
    if (part.dim == 0 && r.tdim() == 2)
    {
        // points, each once
        quadrature::QuadratureRules<T> cells;
        std::vector<EntityRules<T>> entities;
        part_rules(part, 1, false, options, cells, entities);
        emit_points(part, cells, entities, [&](std::span<const T> x, std::int32_t parent)
                    { append_cell(out, cell::type::point, x, parent); });
        return out;
    }
    lut::CellMap<T> map;
    std::vector<I> scratch, nodes;
    LutCell<T> lc;
    LevelSetCell<T, I> ls_scratch;
    std::vector<T> x;
    for (const Entry<I>& e : part_entries(part, include_uncut_cells))
    {
        cell_map(mesh, e.cell, map, scratch);
        const std::int32_t parent = static_cast<std::int32_t>(e.cell);
        if (!e.whole)
        {
            cell_vertex_nodes(mesh, e.cell, nodes, scratch);
            lut_cell(part, e.cell, map, options, ls_scratch, lc);
            for (int p = 0; p < lc.pieces.n_pieces(); ++p)
            {
                if (!lc.selected[static_cast<std::size_t>(p)]
                    || (part.dim > 0 && on_zero_entity(r, map.type, std::span<const I>(nodes), lc, p)))
                    continue;
                lut::push_forward(map, piece_vertices(lc.pieces, p), x);
                append_cell(out, lc.pieces.types[static_cast<std::size_t>(p)], std::span<const T>(x), parent);
            }
        }
        else
            append_cell(out, map.type, std::span<const T>(map.vertices), parent);
    }
    emit_entities(part, [&](cell::type type, std::span<const T> v, std::int32_t parent)
                  { append_cell(out, type, v, parent); });
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
                                                                         const quadrays::Options&, EntityPolicy);
template quadrature::QuadratureRules<double> quadrature_rules<double, int>(const MeshPart<double, int>&, int, bool,
                                                                           const quadrays::Options&, EntityPolicy);
template quadrature::QuadratureRules<float> quadrature_rules<float, int>(const MeshPart<float, int>&, int, bool,
                                                                         const lut::Options&, EntityPolicy);
template quadrature::QuadratureRules<double> quadrature_rules<double, int>(const MeshPart<double, int>&, int, bool,
                                                                           const lut::Options&, EntityPolicy);
template std::vector<EntityRules<float>> entity_rules<float, int>(const MeshPart<float, int>&, int,
                                                                  const quadrays::Options&);
template std::vector<EntityRules<double>> entity_rules<double, int>(const MeshPart<double, int>&, int,
                                                                    const quadrays::Options&);
template std::vector<EntityRules<float>> entity_rules<float, int>(const MeshPart<float, int>&, int,
                                                                  const lut::Options&);
template std::vector<EntityRules<double>> entity_rules<double, int>(const MeshPart<double, int>&, int,
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
