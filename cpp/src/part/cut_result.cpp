// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "cut_result.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>

#include "../bernstein.h"
#include "../cell_topology.h"
#include "../reference_cell.h"
#include "../selection_expr.h"
#include "../tolerances.h"
#include "cell_source.h"
#include "classify.h"

namespace cutcells::part
{

namespace
{

// ============================================================================
// One cell and one level set
// ============================================================================

/// One cell and one level set as cut reads them: quadrays' box of the cell in
/// its reference frame and the level set on it (the cells quadrays takes), or
/// the Bernstein coefficients of the cell's polynomial on its reference cell
/// (other cells, Pk level sets).
template <std::floating_point T, std::integral I>
struct CellReading
{
    cell::type type = cell::type::point;
    bool box = false;
    CellSource<T, I> cs;
    LevelSetCell<T, I> lc;
};

/// Reads @p ls on cell @p c.
template <std::floating_point T, std::integral I>
void read(const MeshView<T, I>& mesh, const LevelSetFunction<T, I>& ls, I c, CellReading<T, I>& rd)
{
    rd.type = mesh.cell_type(c);
    // the reference frame: simplices are classified on their own Bernstein
    // coefficients, and it is the cheapest box to make
    rd.box = cell_source(mesh, ls, c, rd.cs, quadrays::BoxFrame::reference);
    if (rd.box)
        return;
    // cells quadrays does not take: the Bernstein coefficients of their polynomial
    if (ls.analytic)
    {
        throw std::invalid_argument("part::cut: analytic level sets need intervals, triangles, quadrilaterals, "
                                    "tetrahedra, hexahedra, prisms or pyramids, in a space of their dimension or "
                                    "a higher one up to 3D");
    }
    if (!(ls.type == LevelSetType::Polynomial && ls.has_mesh_data() && ls.has_dof_values()))
    {
        throw std::invalid_argument("part::cut: the level set '" + ls.name
                                    + "' has neither dof values nor an analytic level set");
    }
    make_cell_level_set(ls, c, rd.lc);
}

/// The domain of a cell read by read: values within @p touch of 0 count as
/// touching zero, beyond the bounds' own rounding allowance.
template <std::floating_point T, std::integral I>
cell::domain classify(const CellReading<T, I>& rd, int max_depth, T touch)
{
    const int s = rd.box ? cell_sign(rd.cs, rd.type, max_depth, touch)
                         : quadrays::coefficient_sign(std::span<const T>(rd.lc.bernstein_coeffs), touch);
    return s > 0 ? cell::domain::outside : (s < 0 ? cell::domain::inside : cell::domain::intersected);
}

/// The size of the level set on the cell, which its tolerances are relative
/// to: its largest Bernstein coefficient (on quadrays' box), or for an
/// analytic level set its largest value at the corners of the box and at the
/// cell's centroid (a cell may have all its corners on the zero set).
template <std::floating_point T, std::integral I>
T magnitude(const CellReading<T, I>& rd)
{
    if (!rd.box)
        return quadrays::max_abs(std::span<const T>(rd.lc.bernstein_coeffs));
    T m = quadrays::reference_magnitude(rd.cs.source);
    if (rd.cs.source.is_analytic())
    {
        const int nv = cell::get_num_vertices(rd.type);
        quadrays::Vec3<T> centroid = {0, 0, 0};
        for (int v = 0; v < nv; ++v)
        {
            const quadrays::Vec3<T> u = box_vertex<T>(rd.cs.box, rd.type, v);
            for (int k = 0; k < 3; ++k)
                centroid[k] += u[k] / T(nv);
        }
        m = std::max(m, std::abs(quadrays::evaluate(rd.cs.source, std::span<const T>(centroid))));
    }
    return m;
}

/// The level set at vertex @p v (Basix numbering) of the cell: for a Pk level
/// set its nodal value there (the vertex dofs come first; quadrays reads a
/// pyramid's as (1 - z)^n phi, 0 at the apex).
template <std::floating_point T, std::integral I>
T vertex_value(const CellReading<T, I>& rd, int v)
{
    if (rd.box && rd.cs.source.is_analytic())
    {
        const quadrays::Vec3<T> u = box_vertex<T>(rd.cs.box, rd.type, v);
        return quadrays::evaluate(rd.cs.source, std::span<const T>(u));
    }
    return (rd.box ? rd.cs.ls_cell : rd.lc).nodal_values[static_cast<std::size_t>(v)];
}

/// A bound of |phi| over the triangle with physical vertices @p xa, @p xb,
/// @p xc for an analytic level set: Taylor models over the parallelogram of
/// its edges from xa, whose linear part takes its extremes over the triangle at
/// its vertices, or intervals over the triangle's bounding box; infinite if no
/// bound holds.
template <std::floating_point T>
T triangle_bound(const quadrays::AnalyticLevelSet& phi, const quadrays::Vec3<T>& xa, const quadrays::Vec3<T>& xb,
                 const quadrays::Vec3<T>& xc)
{
    if (phi.taylor_bounds != nullptr)
    {
        double centre[3], axes[6], models[3 * 4];
        for (int i = 0; i < 3; ++i)
        {
            axes[2 * i] = 0.5 * (static_cast<double>(xb[i]) - xa[i]);
            axes[2 * i + 1] = 0.5 * (static_cast<double>(xc[i]) - xa[i]);
            centre[i] = xa[i] + axes[2 * i] + axes[2 * i + 1];
        }
        if (quadrays::parallelepiped_bounds(phi, centre, axes, 2, models) == 0)
            return std::numeric_limits<T>::infinity();
        // the triangle's vertices at t = (-1, -1), (1, -1), (-1, 1)
        const double alpha = models[0], b0 = models[1], b1 = models[2], eps = models[3];
        const double linear = std::max({std::abs(alpha - b0 - b1), std::abs(alpha + b0 - b1),
                                        std::abs(alpha - b0 + b1)});
        return static_cast<T>(linear + eps);
    }
    double lo[3], hi[3], b[8];
    for (int i = 0; i < 3; ++i)
    {
        lo[i] = std::min({xa[i], xb[i], xc[i]});
        hi[i] = std::max({xa[i], xb[i], xc[i]});
    }
    if (phi.box_bounds(lo, hi, b, phi.context) == 0)
        return std::numeric_limits<T>::infinity();
    return static_cast<T>(std::max(std::abs(b[0]), std::abs(b[1])));
}

/// A bound of |phi| over sub-entity @p i of dimension @p dim of the cell
/// (Basix numbering), from the cell's reading: the largest coefficient of the
/// restriction of its Bernstein form (exact), Taylor bounds over the entity
/// itself for an analytic level set, the value at a vertex; or, for cells
/// quadrays does not take, the largest value on the entity's lattice of the
/// level set's degree.
template <std::floating_point T, std::integral I>
T entity_bound(const CellReading<T, I>& rd, int dim, int i)
{
    const std::span<const int> sv = cell::sub_entity_vertices(rd.type, dim, i);
    if (dim == 0)
        return std::abs(vertex_value(rd, sv[0]));
    if (rd.box)
    {
        const quadrays::ClippedBox<T>& box = rd.cs.box;
        const quadrays::Vec3<T> ua = box_vertex<T>(box, rd.type, sv[0]), ub = box_vertex<T>(box, rd.type, sv[1]),
                                uc = box_vertex<T>(box, rd.type, sv[dim > 1 ? 2 : 1]);
        // the entity as the image of [0, 1]^dim: ua + s (ub - ua) (+ t (uc - ua))
        std::array<T, 6> matrix{};
        for (int k = 0; k < 3; ++k)
        {
            matrix[k * dim] = ub[k] - ua[k];
            if (dim > 1)
                matrix[k * dim + 1] = uc[k] - ua[k];
        }
        if (rd.cs.source.bernstein != nullptr)
        {
            thread_local quadrays::BoxBernstein<T> restricted;
            thread_local std::vector<T> work;
            const std::size_t d = static_cast<std::size_t>(rd.cs.source.bernstein->dim);
            quadrays::restrict_affine(*rd.cs.source.bernstein, std::span<const T>(ua.data(), d),
                                      std::span<const T>(matrix.data(), d * static_cast<std::size_t>(dim)), dim,
                                      restricted, work);
            return quadrays::max_abs(std::span<const T>(restricted.coeffs));
        }
        if (sv.size() == 3)
            return triangle_bound<T>(*rd.cs.source.analytic, quadrays::physical_point(box, ua),
                                     quadrays::physical_point(box, ub), quadrays::physical_point(box, uc));
        quadrays::AffineBounds<T> b;
        if (!quadrays::affine_bounds(rd.cs.source, std::span<const T>(ua),
                                     std::span<const T>(matrix.data(), static_cast<std::size_t>(3 * dim)), dim, b))
            return std::numeric_limits<T>::infinity();
        return b.magnitude;
    }
    // the lattice a + i/p (b - a) + j/p (c - a), which fixes the polynomial on the entity
    const LevelSetCell<T, I>& lc = rd.lc;
    const int tdim = cell::get_tdim(rd.type), p = std::max(lc.bernstein_order, 1);
    const std::vector<T> ref = cell::reference_vertices<T>(rd.type);
    const T* a = ref.data() + sv[0] * tdim;
    const T* b = ref.data() + sv[1] * tdim;
    const T* e = sv.size() > 2 ? ref.data() + sv[2] * tdim : a;
    const int nj = dim > 1 ? p : 0;
    std::array<T, 3> x{};
    T bound = T(0);
    for (int j = 0; j <= nj; ++j)
        for (int k = 0; k <= p; ++k)
        {
            if (sv.size() == 3 && k + j > p)
                continue;
            for (int d = 0; d < tdim; ++d)
                x[d] = a[d] + T(k) / T(p) * (b[d] - a[d]) + T(j) / T(p) * (e[d] - a[d]);
            bound = std::max(bound, std::abs(bernstein::evaluate<T>(rd.type, lc.bernstein_order,
                                                                    std::span<const T>(lc.bernstein_coeffs),
                                                                    std::span<const T>(x.data(), tdim))));
        }
    return bound;
}

/// The cell's side of the level set next to sub-entity @p i of dimension
/// @p dim: the sign of the level set a thousandth of the way from the
/// entity's centroid towards the cell's, -1, +1, or 0 within @p tol (scaled
/// to that step). Where the level set is smooth this is the sign of its
/// derivative into the cell; it also holds at a kink on the entity, such as
/// that of min(box1, box2) on the face two boxes share.
template <std::floating_point T, std::integral I>
int entity_side(const CellReading<T, I>& rd, int dim, int i, T tol)
{
    constexpr T step = T(1e-3);
    const std::span<const int> sv = cell::sub_entity_vertices(rd.type, dim, i);
    const int nv = cell::get_num_vertices(rd.type), tdim = cell::get_tdim(rd.type);
    quadrays::Vec3<T> centroid = {0, 0, 0}, inner = {0, 0, 0};
    auto point = [&](int v)
    { return rd.box ? box_vertex<T>(rd.cs.box, rd.type, v) : reference_vertex<T>(rd.type, v); };
    for (const int v : sv)
    {
        const quadrays::Vec3<T> u = point(v);
        for (int k = 0; k < 3; ++k)
            centroid[k] += u[k] / T(sv.size());
    }
    for (int v = 0; v < nv; ++v)
    {
        const quadrays::Vec3<T> u = point(v);
        for (int k = 0; k < 3; ++k)
            inner[k] += u[k] / T(nv);
    }
    quadrays::Vec3<T> x;
    for (int k = 0; k < 3; ++k)
        x[k] = centroid[k] + step * (inner[k] - centroid[k]);
    const T value = rd.box ? quadrays::evaluate(rd.cs.source, std::span<const T>(x))
                           : bernstein::evaluate<T>(rd.type, rd.lc.bernstein_order,
                                                    std::span<const T>(rd.lc.bernstein_coeffs),
                                                    std::span<const T>(x.data(), static_cast<std::size_t>(tdim)));
    return value < -step * tol ? -1 : (value > step * tol ? 1 : 0);
}

// ============================================================================
// Entities where a level set may vanish
// ============================================================================

/// A vertex at which a level set is about 0, as a cell sharing it reads it.
template <std::floating_point T, std::integral I>
struct VertexValue
{
    I node = 0;
    int level_set = 0;
    T value = 0; ///< |phi| there

    auto order() const { return std::tie(node, level_set); }
};

/// What the entities need of cell @p c and level set @p l, if the cell lies at
/// the level set's zero set (the level set cuts it, or is within twice
/// touch_tolerance of 0 at a vertex, relative to its magnitude on the cell):
/// the level set's magnitude on the cell, kept in @p magnitude_at_zero if it
/// is the largest so far, and its values at the cell's vertices.
template <std::floating_point T, std::integral I>
void record_vertices(const MeshView<T, I>& mesh, I c, const CellReading<T, I>& rd, int l, cell::domain dom,
                     T& magnitude_at_zero, std::vector<VertexValue<T, I>>& out)
{
    const T m = magnitude(rd);
    if (!(m > T(0)))
        return;
    const int nv = cell::get_num_vertices(rd.type);
    std::array<T, 8> values;
    bool near = dom == cell::domain::intersected;
    for (int v = 0; v < nv; ++v)
    {
        values[static_cast<std::size_t>(v)] = std::abs(vertex_value(rd, v));
        near |= values[static_cast<std::size_t>(v)] <= T(2) * touch_tolerance<T> * m;
    }
    if (!near)
        return;
    magnitude_at_zero = std::max(magnitude_at_zero, m);
    thread_local std::vector<I> nodes, scratch;
    cell_vertex_nodes(mesh, c, nodes, scratch);
    for (int v = 0; v < nv; ++v)
        out.push_back({nodes[static_cast<std::size_t>(v)], l, values[static_cast<std::size_t>(v)]});
}

/// Adds a level set to the result under @p name: itself, or, if it is made of
/// level sets, its own in its place, named "name.part", depth first, and its
/// inside. @p names gets every name in order, a level set made of level sets
/// before its own, and @p owned, for each of these, the range of @p names its
/// own take.
template <std::floating_point T, std::integral I>
void add_level_set(const LevelSetFunction<T, I>& ls, const std::string& name, CutResult<T, I>& r,
                   std::vector<std::string>& names, std::vector<std::array<std::size_t, 2>>& owned)
{
    check_level_sets(ls);
    names.push_back(name);
    if (!ls.has_level_sets())
    {
        r.level_sets.push_back(&ls);
        r.level_set_names.push_back(name);
        return;
    }
    if (!is_selection_name(ls.name))
        throw std::invalid_argument("part::cut: a level set made of level sets needs a name that selections can "
                                    "use, a word without spaces, parentheses, <, > or =; got '" + ls.name + "'");
    r.insides.emplace(name, ls.inside);
    const std::size_t first = names.size();
    for (const LevelSetFunction<T, I>& part : ls.level_sets)
        add_level_set(part, name + "." + part.name, r, names, owned);
    owned.push_back({first, names.size()});
}

/// The cells that some level set cuts, ascending.
template <std::floating_point T, std::integral I>
std::vector<I> cut_cells_of(const CutResult<T, I>& r)
{
    std::vector<I> cells;
    for (I c = 0; c < static_cast<I>(r.num_cells); ++c)
        if (cut_mask(r, c) != 0)
            cells.push_back(c);
    return cells;
}

/// The facets, edges and vertices on which a level set vanishes (CutResult::entities):
/// which level sets vanish on them, the cells beside them classified again
/// with the same threshold, the cells' sides, and the zero faces.
///
/// One threshold per level set, for a whole zero set at once: touch_tolerance
/// times the level set's largest magnitude over the cells at its zero set
/// (@p magnitude_at_zero). A level set vanishes on an entity if the bound of
/// its restriction there is within it, so it vanishes on the edges and
/// vertices of a facet it vanishes on, and a cut cell beside such an entity
/// is classified again with values within it counted as touching zero.
template <std::floating_point T, std::integral I>
void find_entities(CutResult<T, I>& r, std::vector<VertexValue<T, I>>& values,
                   const std::vector<T>& magnitude_at_zero, const ClassifyOptions& options)
{
    const MeshView<T, I>& mesh = *r.mesh;
    const int tdim = r.tdim(), nls = r.n_level_sets();
    if (tdim < 1)
        return;
    auto threshold = [&](int l) { return touch_tolerance<T> * magnitude_at_zero[static_cast<std::size_t>(l)]; };
    struct Vanishing
    {
        int dim;
        std::array<I, 4> key;
        int level_set;
    };
    std::vector<Vanishing> vanishing;

    // vertices: a level set vanishes at one if its value there is within the threshold
    std::sort(values.begin(), values.end(),
              [](const VertexValue<T, I>& a, const VertexValue<T, I>& b) { return a.order() < b.order(); });
    std::vector<std::pair<I, int>> zero_vertices; // (node, level set)
    for (std::size_t i = 0; i < values.size();)
    {
        std::size_t j = i;
        T value = T(0);
        for (; j < values.size() && values[j].order() == values[i].order(); ++j)
            value = std::max(value, values[j].value);
        const I n = values[i].node;
        const int l = values[i].level_set;
        if (value <= threshold(l))
        {
            zero_vertices.emplace_back(n, l);
            vanishing.push_back({0, entity_key(std::span<const I>(&n, 1)), l});
        }
        i = j;
    }

    std::sort(zero_vertices.begin(), zero_vertices.end(),
              [](const auto& a, const auto& b) { return std::tie(a.second, a.first) < std::tie(b.second, b.first); });
    // the cells around those vertices
    std::vector<I> star_nodes;
    for (const auto& [n, l] : zero_vertices)
        star_nodes.push_back(n);
    std::sort(star_nodes.begin(), star_nodes.end());
    star_nodes.erase(std::unique(star_nodes.begin(), star_nodes.end()), star_nodes.end());
    NodeStars<I> stars;
    node_stars(mesh, star_nodes, std::span<const I>(star_nodes), stars);

    // edges and facets whose vertices are all such vertices: decided by the
    // restriction from their lowest sharing cell
    struct Candidate
    {
        int dim;
        std::array<I, 4> key;
        int level_set;
        I cell;
        int local;
        auto order() const { return std::tie(dim, key, level_set, cell); }
    };
    std::vector<Candidate> candidates;
    std::array<I, 4> local;
    std::vector<char> vanishes(static_cast<std::size_t>(mesh.num_nodes()), 0);
    std::vector<I> around;
    for (std::size_t i = 0; i < zero_vertices.size();)
    {
        // the vertices of level set l, and the cells around them, once each
        const int l = zero_vertices[i].second;
        std::size_t j = i;
        around.clear();
        for (; j < zero_vertices.size() && zero_vertices[j].second == l; ++j)
        {
            vanishes[static_cast<std::size_t>(zero_vertices[j].first)] = 1;
            const std::span<const I> cells = stars.cells_around(zero_vertices[j].first);
            around.insert(around.end(), cells.begin(), cells.end());
        }
        std::sort(around.begin(), around.end());
        around.erase(std::unique(around.begin(), around.end()), around.end());
        for (const I c : around)
        {
            const cell::type type = mesh.cell_type(c);
            const std::span<const I> nodes = stars.vertex_nodes(c);
            unsigned mask = 0;
            for (std::size_t v = 0; v < nodes.size(); ++v)
                if (vanishes[static_cast<std::size_t>(nodes[v])])
                    mask |= 1u << v;
            for (int dim = 1; dim < tdim; ++dim)
                for (int e = 0; e < cell::num_sub_entities(type, dim); ++e)
                {
                    const std::span<const int> sv = cell::sub_entity_vertices(type, dim, e);
                    bool all = true;
                    for (std::size_t v = 0; v < sv.size(); ++v)
                    {
                        local[v] = nodes[static_cast<std::size_t>(sv[v])];
                        all &= ((mask >> sv[v]) & 1u) != 0;
                    }
                    if (all)
                        candidates.push_back({dim, entity_key(std::span<const I>(local.data(), sv.size())), l, c, e});
                }
        }
        for (std::size_t k = i; k < j; ++k)
            vanishes[static_cast<std::size_t>(zero_vertices[k].first)] = 0;
        i = j;
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& a, const Candidate& b) { return a.order() < b.order(); });
    // each entity and level set once, with its lowest sharing cell, read once per cell
    std::vector<Candidate> lowest;
    for (std::size_t i = 0; i < candidates.size(); ++i)
        if (i == 0 || candidates[i].dim != candidates[i - 1].dim || candidates[i].key != candidates[i - 1].key
            || candidates[i].level_set != candidates[i - 1].level_set)
            lowest.push_back(candidates[i]);
    std::sort(lowest.begin(), lowest.end(), [](const Candidate& a, const Candidate& b)
              { return std::tie(a.cell, a.level_set) < std::tie(b.cell, b.level_set); });
    CellReading<T, I> rd;
    for (std::size_t i = 0; i < lowest.size(); ++i)
    {
        const Candidate& e = lowest[i];
        if (i == 0 || e.cell != lowest[i - 1].cell || e.level_set != lowest[i - 1].level_set)
            read(mesh, *r.level_sets[static_cast<std::size_t>(e.level_set)], e.cell, rd);
        if (entity_bound(rd, e.dim, e.local) <= threshold(e.level_set))
            vanishing.push_back({e.dim, e.key, e.level_set});
    }

    // the entities: those where a level set vanishes
    std::array<std::vector<std::array<I, 4>>, 3> keys;
    for (const Vanishing& v : vanishing)
        keys[static_cast<std::size_t>(v.dim)].push_back(v.key);
    for (int dim = 0; dim < tdim; ++dim)
    {
        std::vector<std::array<I, 4>>& k = keys[static_cast<std::size_t>(dim)];
        std::sort(k.begin(), k.end());
        k.erase(std::unique(k.begin(), k.end()), k.end());
        r.entities[static_cast<std::size_t>(dim)] = make_entity_table(mesh, dim, std::move(k), stars);
    }

    // one threshold: the cut cells beside an entity a level set vanishes on are
    // classified again, with values within the threshold counted as touching
    std::vector<std::pair<I, int>> touching; // (cell, level set)
    for (const Vanishing& v : vanishing)
    {
        EntityTable<I>& table = r.entities[static_cast<std::size_t>(v.dim)];
        const std::int32_t e = find_entity(table, v.key);
        table.vanishing[static_cast<std::size_t>(e)] |= std::uint64_t(1) << v.level_set;
        for (const I c : table.entity_cells(e))
            if (r.domain(v.level_set, c) == cell::domain::intersected)
                touching.emplace_back(c, v.level_set);
    }
    std::sort(touching.begin(), touching.end());
    touching.erase(std::unique(touching.begin(), touching.end()), touching.end());
    for (const auto& [c, l] : touching)
    {
        read(mesh, *r.level_sets[static_cast<std::size_t>(l)], c, rd);
        r.domains[static_cast<std::size_t>(l) * static_cast<std::size_t>(r.num_cells) + static_cast<std::size_t>(c)]
            = classify(rd, options.max_depth, threshold(l));
    }

    // the sides of the vanishing level sets, and those that count: where they
    // change sign across the entity, or on the mesh boundary
    for (int dim = 0; dim < tdim; ++dim)
    {
        EntityTable<I>& table = r.entities[static_cast<std::size_t>(dim)];
        for (int e = 0; e < table.n_entities(); ++e)
        {
            const std::uint64_t vanish = table.vanishing[static_cast<std::size_t>(e)];
            if (vanish == 0)
                continue;
            std::uint64_t negative = 0, positive = 0;
            for (std::int32_t s = table.cell_offsets[e]; s < table.cell_offsets[e + 1]; ++s)
            {
                const I c = table.cells[static_cast<std::size_t>(s)];
                for (std::uint64_t bits = vanish; bits != 0; bits &= bits - 1)
                {
                    const int l = std::countr_zero(bits);
                    const cell::domain dom = r.domain(l, c);
                    int side = dom == cell::domain::inside ? -1 : (dom == cell::domain::outside ? 1 : 0);
                    if (dom == cell::domain::intersected)
                    {
                        // a cut cell: the level set just inside it, next to the entity
                        read(mesh, *r.level_sets[static_cast<std::size_t>(l)], c, rd);
                        side = entity_side(rd, dim, table.local_index[static_cast<std::size_t>(s)],
                                           touch_tolerance<T> * magnitude(rd));
                    }
                    if (side < 0)
                        table.negative[static_cast<std::size_t>(s)] |= std::uint64_t(1) << l;
                    else if (side > 0)
                        table.positive[static_cast<std::size_t>(s)] |= std::uint64_t(1) << l;
                }
                negative |= table.negative[static_cast<std::size_t>(s)];
                positive |= table.positive[static_cast<std::size_t>(s)];
            }
            table.zero[static_cast<std::size_t>(e)]
                = table.boundary[static_cast<std::size_t>(e)] ? vanish : vanish & negative & positive;
        }
    }

    // the zero faces, each facet once per level set that counts there: owned by
    // the cell on its negative side, else by the lowest cell
    const EntityTable<I>& facets = r.entities[static_cast<std::size_t>(tdim - 1)];
    for (int l = 0; l < nls; ++l)
    {
        const std::uint64_t bit = std::uint64_t(1) << l;
        for (int e = 0; e < facets.n_entities(); ++e)
        {
            if (!(facets.zero[static_cast<std::size_t>(e)] & bit))
                continue;
            std::int32_t owner = facets.cell_offsets[e];
            int negative = 0;
            for (std::int32_t s = facets.cell_offsets[e]; s < facets.cell_offsets[e + 1]; ++s)
                if (facets.negative[static_cast<std::size_t>(s)] & bit)
                {
                    ++negative;
                    owner = s;
                }
            if (negative != 1)
                owner = facets.cell_offsets[e];
            r.zero_face_level_sets.push_back(l);
            r.zero_face_cells.push_back(facets.cells[static_cast<std::size_t>(owner)]);
            r.zero_face_local.push_back(facets.local_index[static_cast<std::size_t>(owner)]);
        }
    }
}

} // namespace

template <std::floating_point T, std::integral I>
CutResult<T, I> cut(const MeshView<T, I>& mesh, std::span<const LevelSetFunction<T, I>> level_sets,
                    const ClassifyOptions& options)
{
    if (!mesh.has_cell_types())
        throw std::invalid_argument("part::cut: the mesh needs cell types");
    if (level_sets.empty())
        throw std::invalid_argument("part::cut: give 1 to 64 level sets");

    CutResult<T, I> r;
    r.mesh = &mesh;
    r.num_cells = static_cast<int>(mesh.num_cells());
    std::vector<std::string> names;
    std::vector<std::array<std::size_t, 2>> owned;
    for (const LevelSetFunction<T, I>& ls : level_sets)
        add_level_set(ls, ls.name, r, names, owned);
    if (r.level_sets.size() > 64)
        throw std::invalid_argument("part::cut: give 1 to 64 level sets, counting those that level sets made of "
                                    "level sets are made of ("
                                    + std::to_string(r.level_sets.size()) + ")");
    std::vector<std::string> sorted = names;
    std::sort(sorted.begin(), sorted.end());
    if (const auto twice = std::adjacent_find(sorted.begin(), sorted.end()); twice != sorted.end())
        throw std::invalid_argument("part::cut: two level sets are named '" + *twice + "'");
    // a level set made of level sets owns the names that begin with its own
    for (const auto& [first, end] : owned)
    {
        const std::string prefix = names[first - 1] + ".";
        for (std::size_t k = 0; k < names.size(); ++k)
            if ((k < first || k >= end) && names[k].starts_with(prefix))
                throw std::invalid_argument("part::cut: the level set named '" + names[k] + "' would read as one of '"
                                            + names[first - 1] + "''s own; rename it");
    }
    const int nls = r.n_level_sets();
    r.domains.assign(static_cast<std::size_t>(nls) * static_cast<std::size_t>(r.num_cells), cell::domain::unset);

    std::vector<VertexValue<T, I>> values;
    std::vector<T> magnitude_at_zero(static_cast<std::size_t>(nls), T(0));
    CellReading<T, I> rd;
    for (I c = 0; c < static_cast<I>(r.num_cells); ++c)
        for (int l = 0; l < nls; ++l)
        {
            read(mesh, *r.level_sets[static_cast<std::size_t>(l)], c, rd);
            const cell::domain dom = classify(rd, options.max_depth, T(0));
            r.domains[static_cast<std::size_t>(l) * static_cast<std::size_t>(r.num_cells)
                      + static_cast<std::size_t>(c)]
                = dom;
            record_vertices(mesh, c, rd, l, dom, magnitude_at_zero[static_cast<std::size_t>(l)], values);
        }
    find_entities(r, values, magnitude_at_zero, options);
    r.cut_cells = cut_cells_of(r);
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
