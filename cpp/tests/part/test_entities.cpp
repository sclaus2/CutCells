// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Pieces lying on facets, edges and vertices of the mesh (part/entities.h,
// part/entity_rules.h), on [-1, 1]^3 and [-1, 1]^2:
//  - the plane x = 0.25 on grid faces of hexahedra, tetrahedra, prisms and
//    pyramids: each facet once, with both cells, their sides, and local
//    numbers that name the facet's nodes in Basix's numbering; its area 4 once
//    in both backends with every option, none with exclude;
//  - planes 1e-16 to 1e-6 off grid faces, Pk and analytic: one threshold, so
//    a facet is in the zero set or a cell beside it is cut, never both, and
//    the area is 4 throughout;
//  - min(box1, box2) of two boxes side by side, an analytic level set with
//    Taylor models: no surface on the shared face, the union's 6;
//  - a P2 level set vanishing on grid faces and cutting the cells beside
//    them: the lookup tables count those facets once;
//  - x = 0.1 and y = 0.1 along mesh edges of Kuhn tetrahedra: length 2 once;
//  - two lines crossing at a grid vertex and on a grid edge (quadrilaterals,
//    triangles): one point of weight 1 with all the cells sharing it.
// Exits non-zero on failure.

#include <cutcells/cell_topology.h>
#include <cutcells/level_set.h>
#include <cutcells/part/cut_result.h>
#include <cutcells/part/entity_rules.h>
#include <cutcells/part/mesh_part.h>
#include <cutcells/part/output.h>
#include <cutcells/quadrays/analytic.h>
#include <cutcells/reference_cell.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "support/box_mesh.h"

using namespace cutcells;
using namespace cutcells::part;
using namespace cutcells::part::support;

namespace
{

int failures = 0;

void check(bool ok, const std::string& what)
{
    if (!ok)
    {
        std::printf("FAILED: %s\n", what.c_str());
        ++failures;
    }
}

bool close(double a, double b, double rel)
{
    return std::abs(a - b) <= rel * std::max(std::abs(b), 1.0);
}

/// The plane x_axis = offset.
struct AxisPlane
{
    int axis = 0;
    double offset = 0;

    template <typename V>
    V operator()(const std::array<V, 3>& x) const
    {
        return x[axis] - offset;
    }
};

/// min(box1, box2) of the boxes [-0.5, 0] x [-0.5, 0.5]^2 and [0, 0.5] x [-0.5, 0.5]^2,
/// each the max of its six half-spaces.
struct BoxUnion
{
    template <typename V>
    V operator()(const std::array<V, 3>& x) const
    {
        using std::max;
        using std::min;
        const double lo[2][3] = {{-0.5, -0.5, -0.5}, {0.0, -0.5, -0.5}};
        const double hi[2][3] = {{0.0, 0.5, 0.5}, {0.5, 0.5, 0.5}};
        V box[2];
        for (int b = 0; b < 2; ++b)
        {
            V m = lo[b][0] - x[0];
            for (int i = 0; i < 3; ++i)
            {
                if (i > 0)
                    m = max(m, lo[b][i] - x[i]);
                m = max(m, x[i] - hi[b][i]);
            }
            box[b] = m;
        }
        return min(box[0], box[1]);
    }
};

using Point = std::array<double, 3>;

LevelSetFunction<double, int> pk_level_set(const MeshView<double, int>& mesh, int degree,
                                           const std::function<double(const Point&)>& f, const std::string& name)
{
    LevelSetMeshData<double, int> data = create_level_set_mesh_data<double, int>(mesh, degree);
    std::vector<double> values(static_cast<std::size_t>(data.num_dofs()));
    for (int d = 0; d < data.num_dofs(); ++d)
    {
        const double* x = data.dof_coordinate(d);
        values[static_cast<std::size_t>(d)] = f({x[0], x[1], mesh.gdim > 2 ? x[2] : 0.0});
    }
    return create_level_set_function<double, int>(std::move(data), std::span<const double>(values), name);
}

/// @p functor must outlive the level set: the analytic level set points to it.
template <typename F>
LevelSetFunction<double, int> analytic_ls(const F& functor, const std::string& name, int gdim = 3)
{
    auto phi = std::make_shared<const quadrays::AnalyticLevelSet>(quadrays::analytic_level_set(functor));
    return create_level_set_function<double, int>(phi, gdim, name);
}

double total(const quadrature::QuadratureRules<double>& rules)
{
    double s = 0;
    for (const double w : rules._weights)
        s += w;
    return s;
}

double total(const EntityRules<double>& rules)
{
    double s = 0;
    for (const double w : rules.weights)
        s += w;
    return s;
}

/// A mesh with its storage; build in place, do not copy.
struct Mesh
{
    std::vector<double> coordinates;
    std::vector<int> connectivity, offsets = {0};
    std::vector<cell::type> types;
    MeshView<double, int> view;
};

/// [-1, 1]^2 with n^2 quadrilaterals or 2 n^2 triangles (split along x = y),
/// Basix vertex order.
void make_square_mesh(const std::string& kind, int n, Mesh& mesh)
{
    const double h = 2.0 / n;
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i)
            mesh.coordinates.insert(mesh.coordinates.end(), {-1 + h * i, -1 + h * j});
    auto node = [n](int i, int j) { return i + (n + 1) * j; };
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
        {
            const std::vector<std::vector<int>> cells
                = kind == "quad"
                      ? std::vector<std::vector<int>>{{node(i, j), node(i + 1, j), node(i, j + 1), node(i + 1, j + 1)}}
                      : std::vector<std::vector<int>>{{node(i, j), node(i + 1, j), node(i + 1, j + 1)},
                                                      {node(i, j), node(i + 1, j + 1), node(i, j + 1)}};
            for (const std::vector<int>& c : cells)
            {
                mesh.connectivity.insert(mesh.connectivity.end(), c.begin(), c.end());
                mesh.offsets.push_back(static_cast<int>(mesh.connectivity.size()));
                mesh.types.push_back(kind == "quad" ? cell::type::quadrilateral : cell::type::triangle);
            }
        }
    mesh.view.gdim = 2;
    mesh.view.tdim = 2;
    mesh.view.coordinates = mesh.coordinates;
    mesh.view.connectivity = mesh.connectivity;
    mesh.view.offsets = mesh.offsets;
    mesh.view.cell_types = mesh.types;
}

/// Kuhn tetrahedra of the tensor grid g x g x gz, Basix vertex order.
void make_tensor_tet_mesh(const std::vector<double>& g, const std::vector<double>& gz, Mesh& mesh)
{
    const int n = static_cast<int>(g.size());
    for (const double z : gz)
        for (const double y : g)
            for (const double x : g)
                mesh.coordinates.insert(mesh.coordinates.end(), {x, y, z});
    auto node = [n](const std::array<int, 3>& v) { return v[0] + n * (v[1] + n * v[2]); };
    for (int k = 0; k + 1 < static_cast<int>(gz.size()); ++k)
        for (int j = 0; j + 1 < n; ++j)
            for (int i = 0; i + 1 < n; ++i)
            {
                std::array<int, 3> perm = {0, 1, 2};
                do
                {
                    std::array<int, 3> v = {i, j, k};
                    mesh.connectivity.push_back(node(v));
                    for (const int axis : perm)
                    {
                        ++v[static_cast<std::size_t>(axis)];
                        mesh.connectivity.push_back(node(v));
                    }
                    mesh.offsets.push_back(static_cast<int>(mesh.connectivity.size()));
                    mesh.types.push_back(cell::type::tetrahedron);
                } while (std::next_permutation(perm.begin(), perm.end()));
            }
    mesh.view.gdim = 3;
    mesh.view.tdim = 3;
    mesh.view.coordinates = mesh.coordinates;
    mesh.view.connectivity = mesh.connectivity;
    mesh.view.offsets = mesh.offsets;
    mesh.view.cell_types = mesh.types;
}

const EntityPolicy policies[] = {EntityPolicy::split, EntityPolicy::negative, EntityPolicy::positive,
                                 EntityPolicy::owner};

/// The rules of @p part with every option and both backends, against @p exact.
bool totals_with_every_option(const MeshPart<double, int>& part, double exact, double rel, bool full = false)
{
    bool ok = true;
    for (const EntityPolicy p : policies)
    {
        ok &= close(total(quadrature_rules(part, 3, full, quadrays::Options{}, p)), exact, rel);
        ok &= close(total(quadrature_rules(part, 3, full, lut::Options{}, p)), exact, rel);
    }
    return ok;
}

/// The plane x = 0.25 on grid faces.
void test_plane_on_grid_faces()
{
    for (const char* kind : {"hex", "tet", "prism", "pyramid"})
    {
        BoxMesh mesh;
        make_box_mesh(kind, 8, {0, 0, 0}, mesh);
        const AxisPlane plane = {0, 0.25};
        for (int source = 0; source < 2; ++source)
        {
            const std::vector<LevelSetFunction<double, int>> ls
                = {source == 0 ? pk_level_set(mesh.view, 1, [](const Point& x) { return x[0] - 0.25; }, "phi")
                               : analytic_ls(plane, "phi")};
            const CutResult<double, int> r = cut<double, int>(mesh.view, ls);
            const EntityTable<int>& facets = r.entities[2];
            int zero = 0;
            bool sides = true, local = true;
            std::vector<int> nodes, scratch;
            for (int e = 0; e < facets.n_entities(); ++e)
            {
                const std::span<const int> en = facets.entity_nodes(e);
                // every entity's local numbers name its nodes
                for (std::int32_t s = facets.cell_offsets[e]; s < facets.cell_offsets[e + 1]; ++s)
                {
                    const int c = facets.cells[static_cast<std::size_t>(s)];
                    cell_vertex_nodes(mesh.view, c, nodes, scratch);
                    std::vector<int> held;
                    for (const int v : cell::sub_entity_vertices(mesh.view.cell_type(c), 2,
                                                                  facets.local_index[static_cast<std::size_t>(s)]))
                        held.push_back(nodes[static_cast<std::size_t>(v)]);
                    local &= std::is_permutation(held.begin(), held.end(), en.begin(), en.end());
                }
                if (!(facets.zero[static_cast<std::size_t>(e)] & 1))
                    continue;
                ++zero;
                sides &= facets.cell_offsets[e + 1] - facets.cell_offsets[e] == 2;
                for (std::int32_t s = facets.cell_offsets[e]; s < facets.cell_offsets[e + 1]; ++s)
                {
                    const int c = facets.cells[static_cast<std::size_t>(s)];
                    double cx = 0;
                    cell_vertex_nodes(mesh.view, c, nodes, scratch);
                    for (const int n : nodes)
                        cx += mesh.view.node(n)[0] / static_cast<double>(nodes.size());
                    const bool below = cx < 0.25;
                    sides &= facets.negative[static_cast<std::size_t>(s)] == (below ? 1u : 0u)
                             && facets.positive[static_cast<std::size_t>(s)] == (below ? 0u : 1u);
                }
            }
            const int expected = std::string(kind) == "tet" ? 128 : 64;
            const MeshPart<double, int> part = select(r, "phi = 0");
            const std::vector<EntityRules<double>> on = entity_rules(part, 3, quadrays::Options{});
            const double excluded = total(quadrature_rules(part, 3, false, quadrays::Options{}, EntityPolicy::exclude));
            const double below = total(quadrature_rules(select(r, "phi < 0"), 3, true));
            const bool totals = totals_with_every_option(part, 4.0, 1e-13) && close(total(on[2]), 4.0, 1e-13)
                                && on[2].n_rules() == expected && std::abs(excluded) < 1e-13
                                && close(below, 5.0, 1e-12);
            if (!totals)
                std::printf("  totals: every option %d, entity rules %.17g (%d), excluded %.3g, below %.17g\n",
                            totals_with_every_option(part, 4.0, 1e-13), total(on[2]), on[2].n_rules(), excluded,
                            below);
            std::printf("plane in faces %s %s: %d zero facets (%d expected), sides %d, local numbers %d, totals %d\n",
                        kind, source == 0 ? "P1" : "analytic", zero, expected, sides, local, totals);
            check(zero == expected && sides && local && totals && r.n_zero_faces() == expected,
                  std::string("plane on grid faces, ") + kind);
        }
    }
}

/// Planes off grid faces by 1e-16 to 1e-6: one threshold, the area 4.
void test_planes_off_grid_faces()
{
    for (const char* kind : {"hex", "tet"})
    {
        BoxMesh mesh;
        make_box_mesh(kind, 8, {0, 0, 0}, mesh);
        for (int source = 0; source < 2; ++source)
        {
            double worst = 0;
            bool one_threshold = true;
            for (const double sign : {1.0, -1.0})
                for (const double d : {1e-16, 1e-15, 3e-15, 1e-14, 2e-14, 1e-13, 1e-12, 1e-9, 1e-6})
                {
                    const double offset = 0.25 + sign * d;
                    const AxisPlane plane = {0, offset};
                    const std::vector<LevelSetFunction<double, int>> ls
                        = {source == 0
                               ? pk_level_set(mesh.view, 1, [offset](const Point& x) { return x[0] - offset; }, "phi")
                               : analytic_ls(plane, "phi")};
                    const CutResult<double, int> r = cut<double, int>(mesh.view, ls);
                    one_threshold &= r.cut_cells.empty() || r.n_zero_faces() == 0;
                    const MeshPart<double, int> part = select(r, "phi = 0");
                    for (const EntityPolicy p : policies)
                    {
                        worst = std::max(worst, std::abs(total(quadrature_rules(part, 2, false, quadrays::Options{}, p)) - 4));
                        worst = std::max(worst, std::abs(total(quadrature_rules(part, 2, false, lut::Options{}, p)) - 4));
                    }
                }
            std::printf("planes off faces %s %s: largest error of the area %.3g, one threshold %d\n", kind,
                        source == 0 ? "P1" : "analytic", worst, one_threshold);
            check(worst < 1e-12 && one_threshold, std::string("planes off grid faces, ") + kind);
        }
    }
}

/// min(box1, box2) as one level set: it vanishes on the shared face without
/// changing sign there, so that face is no surface. This checks the sign-change
/// rule; a shape with kinks belongs decomposed into its smooth level sets (a
/// level set made of level sets), whose edges and corners come from
/// intersecting them (test_level_sets_made_of). Tetrahedra take intervals over
/// boxes here: Taylor models of min and max over a tetrahedron's
/// parallelepiped cross the kink on the shared face and leave the cells
/// around it cut, where quadrays does not integrate a kink.
void test_two_boxes()
{
    for (const char* kind : {"hex", "tet"})
    {
        BoxMesh mesh;
        make_box_mesh(kind, 4, {0, 0, 0}, mesh);
        const BoxUnion boxes;
        auto phi = std::make_shared<quadrays::AnalyticLevelSet>(quadrays::analytic_level_set(boxes));
        if (std::string(kind) == "tet")
        {
            phi->taylor_bounds = nullptr;
            phi->hessian_bounds = nullptr;
        }
        const std::vector<LevelSetFunction<double, int>> ls
            = {create_level_set_function<double, int>(std::shared_ptr<const quadrays::AnalyticLevelSet>(phi), 3,
                                                      "phi")};
        const CutResult<double, int> r = cut<double, int>(mesh.view, ls);
        const EntityTable<int>& facets = r.entities[2];
        int shared = 0, counted = 0;
        for (int e = 0; e < facets.n_entities(); ++e)
        {
            bool on_plane = true;
            for (const int n : facets.entity_nodes(e))
                on_plane &= std::abs(mesh.view.node(n)[0]) < 1e-14;
            if (!on_plane)
                continue;
            ++shared;
            counted += facets.zero[static_cast<std::size_t>(e)] != 0;
        }
        const MeshPart<double, int> part = select(r, "phi = 0");
        const double inside = total(quadrature_rules(select(r, "phi < 0"), 3, true));
        const bool totals = totals_with_every_option(part, 6.0, 1e-13) && close(inside, 1.0, 1e-12);
        if (!totals)
            for (const EntityPolicy p : policies)
                std::printf("  %s: quadrays %.17g, lookup tables %.17g, volume %.17g\n", entity_policy_to_str(p).c_str(),
                            total(quadrature_rules(part, 3, false, quadrays::Options{}, p)),
                            total(quadrature_rules(part, 3, false, lut::Options{}, p)), inside);
        std::printf("two boxes %s: %zu cut cells, %d facets on the shared face (%d counted), totals %d\n", kind,
                    r.cut_cells.size(), shared, counted, totals);
        check(r.cut_cells.empty() && shared == (std::string(kind) == "tet" ? 8 : 4) && counted == 0 && totals,
              std::string("two boxes as one level set, ") + kind);
    }
}

/// (x - 0.25)(x - 0.45): its facets at x = 0.25 once in the lookup tables.
void test_lut_facets_of_cut_cells()
{
    for (const char* kind : {"hex", "tet"})
    {
        BoxMesh mesh;
        make_box_mesh(kind, 8, {0, 0, 0}, mesh);
        const std::vector<LevelSetFunction<double, int>> ls
            = {pk_level_set(mesh.view, 2, [](const Point& x) { return (x[0] - 0.25) * (x[0] - 0.45); }, "phi")};
        const CutResult<double, int> r = cut<double, int>(mesh.view, ls);
        const MeshPart<double, int> part = select(r, "phi = 0");
        const double lut = total(quadrature_rules(part, 2, false, lut::Options{}));
        const double qr = total(quadrature_rules(part, 2, false, quadrays::Options{}));
        std::printf("facets of cut cells %s: lookup tables %.15f, quadrays %.15f (8)\n", kind, lut, qr);
        check(close(lut, 8.0, 1e-12) && close(qr, 8.0, 1e-12) && totals_with_every_option(part, 8.0, 1e-12),
              std::string("lookup tables count facets of cut cells once, ") + kind);
    }
}

/// x = 0.1 and y = 0.1 along mesh edges.
void test_line_along_edges()
{
    Mesh mesh;
    make_tensor_tet_mesh({-1.0, -0.5, 0.0, 0.1, 0.5, 1.0}, {-1.0, -0.5, 0.0, 0.5, 1.0}, mesh);
    const AxisPlane plane = {1, 0.1};
    const std::vector<LevelSetFunction<double, int>> ls
        = {pk_level_set(mesh.view, 1, [](const Point& x) { return x[0] - 0.1; }, "a"), analytic_ls(plane, "b")};
    const CutResult<double, int> r = cut<double, int>(mesh.view, ls);
    const MeshPart<double, int> part = select(r, "a = 0 and b = 0");
    const std::vector<EntityRules<double>> on = entity_rules(part, 2, quadrays::Options{});
    bool cells = true;
    std::vector<int> nodes, scratch;
    for (int k = 0; k < on[1].n_rules(); ++k)
    {
        // every tetrahedron holding both nodes shares the edge
        const int a = on[1].entity_nodes[static_cast<std::size_t>(on[1].node_offsets[k])],
                  b = on[1].entity_nodes[static_cast<std::size_t>(on[1].node_offsets[k]) + 1];
        int holding = 0;
        for (int c = 0; c < static_cast<int>(mesh.types.size()); ++c)
        {
            cell_vertex_nodes(mesh.view, c, nodes, scratch);
            holding += std::count(nodes.begin(), nodes.end(), a) + std::count(nodes.begin(), nodes.end(), b) == 2;
        }
        cells &= holding == on[1].cell_offsets[k + 1] - on[1].cell_offsets[k];
    }
    const bool totals = totals_with_every_option(part, 2.0, 1e-14);
    std::printf("line along edges: %d edge rules, length %.15f, sharing cells %d, totals %d\n", on[1].n_rules(),
                total(on[1]), cells, totals);
    check(on[1].n_rules() == 4 && close(total(on[1]), 2.0, 1e-14) && cells && totals && r.cut_cells.empty(),
          "x = 0.1 and y = 0.1 along mesh edges");
}

/// Two lines crossing at a grid vertex and on a grid edge.
void test_points_on_entities()
{
    for (const char* kind : {"quad", "tri"})
        for (const bool vertex : {true, false})
        {
            Mesh mesh;
            make_square_mesh(kind, 8, mesh);
            const Point p = {0.25, vertex ? 0.5 : 0.4, 0};
            const std::vector<LevelSetFunction<double, int>> ls
                = {pk_level_set(mesh.view, 1, [p](const Point& x) { return 0.37 * (x[0] - p[0]) - (x[1] - p[1]); },
                                "a"),
                   pk_level_set(mesh.view, 1, [p](const Point& x) { return -3.45 * (x[0] - p[0]) - (x[1] - p[1]); },
                                "b")};
            const CutResult<double, int> r = cut<double, int>(mesh.view, ls);
            const MeshPart<double, int> part = select(r, "a = 0 and b = 0");
            bool ok = true;
            for (const bool lut_backend : {false, true})
            {
                const std::vector<EntityRules<double>> on = lut_backend ? entity_rules(part, 2, lut::Options{})
                                                                        : entity_rules(part, 2, quadrays::Options{});
                const EntityRules<double>& rules = on[vertex ? 0 : 1];
                const int n = vertex ? (std::string(kind) == "quad" ? 4 : 6) : 2;
                ok &= rules.n_rules() == 1 && on[vertex ? 1 : 0].n_rules() == 0 && rules.weights.size() == 1
                      && rules.weights[0] == 1.0 && static_cast<int>(rules.cells.size()) == n;
            }
            ok &= totals_with_every_option(part, 1.0, 1e-14);
            std::printf("crossing %s %s: %d\n", vertex ? "at a vertex" : "on an edge", kind, ok);
            check(ok, std::string("points on entities, ") + kind);
        }
}

} // namespace

int main()
{
    test_plane_on_grid_faces();
    test_planes_off_grid_faces();
    test_two_boxes();
    test_lut_facets_of_cut_cells();
    test_line_along_edges();
    test_points_on_entities();
    return failures == 0 ? 0 : 1;
}
