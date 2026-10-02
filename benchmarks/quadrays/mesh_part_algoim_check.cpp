// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Checks of algoim_rules (mesh_part_algoim.h) on single cells, for both of
// algoim's engines: linear cuts of a quadrilateral and a hexahedron against
// exact measures, the root of an interval embedded in 2D, and the interface
// of a quadrilateral embedded in 3D, whose weights see the induced metric.
// Exits non-zero on failure; ctest runs it with CUTCELLS_WITH_ALGOIM.

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include <cutcells/level_set.h>
#include <cutcells/mesh_view.h>
#include <cutcells/part/cut_result.h>
#include <cutcells/part/mesh_part.h>

#include "mesh_part_algoim.h"

using namespace cutcells;

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

/// One cell and the storage its view points to: build in place, do not copy.
struct OneCell
{
    std::vector<double> coordinates;
    std::vector<int> connectivity, offsets;
    std::vector<cell::type> types;
    MeshView<double, int> view;
};

void make_cell(cell::type type, int gdim, std::vector<double> coordinates, OneCell& m)
{
    m.coordinates = std::move(coordinates);
    const int nv = cell::get_num_vertices(type);
    m.connectivity.resize(static_cast<std::size_t>(nv));
    for (int v = 0; v < nv; ++v)
        m.connectivity[static_cast<std::size_t>(v)] = v;
    m.offsets = {0, nv};
    m.types = {type};
    m.view.gdim = gdim;
    m.view.tdim = cell::get_tdim(type);
    m.view.coordinates = m.coordinates;
    m.view.connectivity = m.connectivity;
    m.view.offsets = m.offsets;
    m.view.cell_types = m.types;
}

/// The P1 level set of @p f on @p mesh.
LevelSetFunction<double, int> p1(const MeshView<double, int>& mesh, const std::function<double(const double*)>& f)
{
    LevelSetMeshData<double, int> data = create_level_set_mesh_data<double, int>(mesh, 1);
    std::vector<double> values(static_cast<std::size_t>(data.num_dofs()));
    for (int d = 0; d < data.num_dofs(); ++d)
        values[static_cast<std::size_t>(d)] = f(data.dof_coordinate(d));
    return create_level_set_function<double, int>(std::move(data), std::span<const double>(values), "phi");
}

double total(const quadrature::QuadratureRules<double>& rules)
{
    double s = 0;
    for (const double w : rules._weights)
        s += w;
    return s;
}

void check_cut(const OneCell& m, const std::function<double(const double*)>& f, double volume, double interface,
               const std::string& what)
{
    const std::vector<LevelSetFunction<double, int>> ls = {p1(m.view, f)};
    const part::CutResult<double, int> r = part::cut<double, int>(m.view, ls);
    for (const bool general : {false, true})
    {
        const double v = total(benchmarks::algoim_rules(part::select(r, "phi < 0"), 4, false, general));
        const double a = total(benchmarks::algoim_rules(part::select(r, "phi = 0"), 4, false, general));
        std::printf("%s, %s engine: volume %.1e, interface %.1e\n", what.c_str(), general ? "general" : "Bernstein",
                    v - volume, a - interface);
        check(std::abs(v - volume) < 1e-12 && std::abs(a - interface) < 1e-12, what);
    }
}

} // namespace

int main()
{
    OneCell quad, hex, interval, embedded;
    make_cell(cell::type::quadrilateral, 2, {0, 0, 1, 0, 0, 1, 1, 1}, quad);
    check_cut(quad, [](const double* x) { return x[0] + x[1] - 0.8; }, 0.5 * 0.8 * 0.8, std::sqrt(2.0) * 0.8,
              "quadrilateral");
    make_cell(cell::type::hexahedron, 3, {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1}, hex);
    check_cut(hex, [](const double* x) { return x[0] + x[1] + x[2] - 0.8; }, 0.8 * 0.8 * 0.8 / 6.0,
              std::sqrt(3.0) * 0.8 * 0.8 / 2.0, "hexahedron");

    // the root of an interval embedded in 2D, at 0.3
    make_cell(cell::type::interval, 2, {0, 0, 1, 0}, interval);
    {
        const std::vector<LevelSetFunction<double, int>> ls = {p1(interval.view, [](const double* x) { return x[0] - 0.3; })};
        const part::CutResult<double, int> r = part::cut<double, int>(interval.view, ls);
        for (const bool general : {false, true})
        {
            const auto rules = benchmarks::algoim_rules(part::select(r, "phi = 0"), 4, false, general);
            check(rules._weights.size() == 1 && std::abs(rules._weights[0] - 1) < 1e-14
                      && std::abs(rules._points[0] - 0.3) < 1e-12,
                  "the root of an interval");
        }
    }

    // a 2 x 3 quadrilateral in the plane y = 0 of 3D, cut at x = 0.8: the
    // interface is 3 long and lies at the reference coordinate 0.4
    make_cell(cell::type::quadrilateral, 3, {0, 0, 0, 2, 0, 0, 0, 0, 3, 2, 0, 3}, embedded);
    {
        const std::vector<LevelSetFunction<double, int>> ls = {p1(embedded.view, [](const double* x) { return x[0] - 0.8; })};
        const part::CutResult<double, int> r = part::cut<double, int>(embedded.view, ls);
        for (const bool general : {false, true})
        {
            const auto rules = benchmarks::algoim_rules(part::select(r, "phi = 0"), 4, false, general);
            bool on_line = !rules._weights.empty();
            for (std::size_t q = 0; q < rules._weights.size(); ++q)
                on_line &= std::abs(rules._points[2 * q] - 0.4) < 1e-12;
            std::printf("embedded quadrilateral, %s engine: length %.1e\n", general ? "general" : "Bernstein",
                        total(rules) - 3);
            check(std::abs(total(rules) - 3) < 1e-12 && on_line, "the interface of an embedded quadrilateral");
        }
    }
    return failures == 0 ? 0 : 1;
}
