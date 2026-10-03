// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Sweeping circle and sphere: a circle (sphere) of radius r whose centre
// sweeps a cell, c = c0 + s h (1, 0.618..., 0.414...) for s = k / steps, so
// that the cuts take many positions relative to the grid, from the centre on
// a grid vertex (s = 0) on. The level set is the P2 interpolant of
// |x - c|^2 - r^2, exact on every cell, so that errors are the quadrature's;
// with --degree 1 its P1 interpolant (multilinear on quadrilaterals and
// hexahedra). The parts of the front end go to quadrays (part::quadrature_rules),
// to the lookup tables (lut: the Pk-iso-P1 template of the level set's degree,
// lut:k: of order k) and, with CUTCELLS_WITH_ALGOIM on quadrilaterals and
// hexahedra, to algoim's multi-polynomial engine on the same Bernstein forms
// (algoim_rules, its AutoMixed strategy).
//
// Per step and part, against exact per-cell values (exact_reference.h): the
// error of the total (cut cells plus whole cells inside), and the per-cell L1
// over the cut cells, relative to the measure of the disk (ball) or circle
// (sphere); the worst cell error relative to h^d (volume) or h^(d-1)
// (interface); robustness counters: negative or non-finite weights, points
// outside their cell or outside the part beyond 1e-8 h, generators that threw;
// points and microseconds per cut cell. With --degree 1 the reference is the
// part of the P1 level set, by quadrays with --ref-q points per segment (the
// lookup tables are exact on triangles and tetrahedra, which checks it), and
// the measures of the disk (ball) and circle (sphere) are those of that part.
// The P1 interpolant of |x - c|^2 - r^2 is affine on the axis-aligned
// quadrilaterals, hexahedra and prisms (it has no cross terms); --tilt a adds
// a sum_{i<j} (x_i - c_i)(x_j - c_j), an ellipse (ellipsoid) with tilted axes
// for |a| < 1, whose interpolant cuts those cells along curves. With --compress p, every rule is
// compressed onto P_p (triangles, tetrahedra) or Q_p (the other cells), and
// the points and microseconds per cut cell of the compressed rules are given.
// Run with OMP_NUM_THREADS=1: compression uses OpenMP.
//
// Usage:
//   quadrays_sweep [--mesh quad|tri|hex|tet|prism|prism-diag|pyramid] [--n 16] [--q 3,5] [--steps 64]
//                  [--radius r] [--centre x,y,z] [--gen quadrays,lut,lut:k,algoim]
//                  [--degree 1|2] [--ref-q 16] [--tilt a] [--compress p] [--csv file]

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <cutcells/compression/compress.h>
#include <cutcells/level_set.h>
#include <cutcells/level_set_cell.h>
#include <cutcells/lut/cell_pieces.h>
#include <cutcells/mesh_view.h>
#include <cutcells/part/cut_result.h>
#include <cutcells/part/mesh_part.h>
#include <cutcells/part/output.h>

#include "exact_reference.h"
#include "test_mesh.h"

#ifdef CUTCELLS_WITH_ALGOIM
#include "mesh_part_algoim.h"
#endif

using namespace cutcells;
using namespace cutcells::quadrays::support;

namespace
{
struct Config
{
    std::string mesh = "quad";
    int n = 16;
    std::vector<int> qs = {3, 5};
    int steps = 64;
    double radius = 0.7;
    V3 centre = {0, 0, 0};
    std::vector<std::string> generators = {"quadrays", "algoim"};
    int degree = 2;
    int ref_q = 16;
    double tilt = 0;
    int compress = -1;
    std::string csv;
};

std::vector<std::string> split_list(const std::string& s)
{
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ','))
        if (!item.empty())
            out.push_back(item);
    return out;
}

Config parse_args(int argc, char** argv)
{
    Config cfg;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        auto next = [&]() -> std::string
        {
            if (i + 1 >= argc)
                throw std::runtime_error("missing value after " + a);
            return argv[++i];
        };
        if (a == "--mesh")
            cfg.mesh = next();
        else if (a == "--n")
            cfg.n = std::stoi(next());
        else if (a == "--q")
        {
            cfg.qs.clear();
            for (const std::string& s : split_list(next()))
                cfg.qs.push_back(std::stoi(s));
        }
        else if (a == "--steps")
            cfg.steps = std::stoi(next());
        else if (a == "--radius")
            cfg.radius = std::stod(next());
        else if (a == "--centre")
        {
            const std::vector<std::string> c = split_list(next());
            for (std::size_t d = 0; d < c.size() && d < 3; ++d)
                cfg.centre[d] = std::stod(c[d]);
        }
        else if (a == "--gen")
            cfg.generators = split_list(next());
        else if (a == "--degree")
            cfg.degree = std::stoi(next());
        else if (a == "--ref-q")
            cfg.ref_q = std::stoi(next());
        else if (a == "--tilt")
            cfg.tilt = std::stod(next());
        else if (a == "--compress")
            cfg.compress = std::stoi(next());
        else if (a == "--csv")
            cfg.csv = next();
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    if (cfg.mesh != "quad" && cfg.mesh != "tri" && cfg.mesh != "hex" && cfg.mesh != "tet" && cfg.mesh != "prism"
        && cfg.mesh != "prism-diag" && cfg.mesh != "pyramid")
        throw std::runtime_error("--mesh takes quad, tri, hex, tet, prism, prism-diag or pyramid");
    if (cfg.degree != 1 && cfg.degree != 2)
        throw std::runtime_error("--degree takes 1 or 2");
    if (cfg.tilt != 0 && cfg.degree != 1)
        throw std::runtime_error("--tilt needs --degree 1: the exact references are those of the circle (sphere)");
    return cfg;
}

/// The test cells of [-1, 1]^d in grid order, as a mesh view over copies of
/// their vertices (Basix order).
struct SweepMesh
{
    std::string kind;
    int n = 0;
    int tdim = 2;
    double h = 0;
    int per_cube = 1;        ///< cells per grid square or cube
    double cell_measure = 0; ///< area or volume of every cell
    std::vector<double> coordinates;
    std::vector<int> connectivity, offsets = {0};
    std::vector<cell::type> types;
    MeshView<double, int> view;
};

std::vector<TestCell2D> square_cells(const SweepMesh& m, int square)
{
    const int i0 = square / m.n, i1 = square % m.n;
    return grid_cells_2d(m.kind, {-1 + m.h * i0, -1 + m.h * i1, 0}, m.h);
}

std::vector<TestCell> cube_cells(const SweepMesh& m, int cube, const V3& origin)
{
    const int i0 = cube / (m.n * m.n), i1 = (cube / m.n) % m.n, i2 = cube % m.n;
    return grid_cells(m.kind, {-1 + m.h * i0, -1 + m.h * i1, -1 + m.h * i2}, m.h, origin);
}

void build_mesh(const Config& cfg, SweepMesh& m)
{
    m.kind = cfg.mesh;
    m.n = cfg.n;
    m.tdim = cfg.mesh == "quad" || cfg.mesh == "tri" ? 2 : 3;
    m.h = 2.0 / cfg.n;
    m.per_cube = cfg.mesh == "tri" || cfg.mesh == "prism" || cfg.mesh == "prism-diag"
                     ? 2
                     : (cfg.mesh == "tet" || cfg.mesh == "pyramid" ? 6 : 1);
    m.cell_measure = std::pow(m.h, m.tdim) / m.per_cube;
    const int cubes = m.tdim == 2 ? m.n * m.n : m.n * m.n * m.n;
    for (int k = 0; k < cubes; ++k)
    {
        auto add = [&](cell::type type, const std::vector<double>& vertices)
        {
            const int first = static_cast<int>(m.coordinates.size()) / m.tdim;
            const int nv = static_cast<int>(vertices.size()) / m.tdim;
            m.coordinates.insert(m.coordinates.end(), vertices.begin(), vertices.end());
            for (int v = 0; v < nv; ++v)
                m.connectivity.push_back(first + v);
            m.offsets.push_back(static_cast<int>(m.connectivity.size()));
            m.types.push_back(type);
        };
        if (m.tdim == 2)
            for (const TestCell2D& c : square_cells(m, k))
                add(c.type, c.vertices);
        else
            for (const TestCell& c : cube_cells(m, k, {0, 0, 0}))
                add(c.type, c.vertices);
    }
    m.view.gdim = m.tdim;
    m.view.tdim = m.tdim;
    m.view.coordinates = m.coordinates;
    m.view.connectivity = m.connectivity;
    m.view.offsets = m.offsets;
    m.view.cell_types = m.types;
}

/// One cut cell: its exact shares of the disk (ball) and the circle (sphere),
/// its type, and the physical points of its reference coordinates.
struct CellInfo
{
    cell::type type = cell::type::point;
    double exact_volume = 0, exact_surface = 0;
    TestCell2D c2;
    TestCell c3;

    V3 physical_point(std::span<const double> xi) const
    {
        if (c2.loop.empty())
            return physical(c3, V3{xi[0], xi[1], xi[2]});
        return physical(c2, std::array<double, 2>{xi[0], xi[1]});
    }
};

CellInfo cell_info(const SweepMesh& m, int cell, const V3& c, double r)
{
    CellInfo info;
    const int cube = cell / m.per_cube, sub = cell % m.per_cube;
    if (m.tdim == 2)
    {
        info.c2 = square_cells(m, cube)[static_cast<std::size_t>(sub)];
        info.type = info.c2.type;
        info.exact_surface = exact::circle_polygon_length(info.c2.loop, c, r);
        info.exact_volume = exact::disk_polygon_area(info.c2.loop, c, r);
    }
    else
    {
        info.c3 = cube_cells(m, cube, c)[static_cast<std::size_t>(sub)]; // faces relative to c
        info.type = info.c3.type;
        info.exact_surface = exact::sphere_area(info.c3.faces, r);
        info.exact_volume = exact::ball_volume(info.c3.faces, r, info.exact_surface);
    }
    return info;
}

/// The cut cells of one centre: their CellInfo, and where each cell is among
/// them (-1: not cut); the measures of the disk (ball) and circle (sphere).
struct CutCells
{
    std::vector<CellInfo> infos;
    std::vector<int> where;
    double whole_volume = 0, whole_surface = 0;
};

quadrature::QuadratureRules<double> rules_of(const std::string& gen, const part::MeshPart<double, int>& part, int q)
{
    if (gen == "quadrays")
        return part::quadrature_rules(part, q, false, quadrays::Options{});
    if (gen == "lut" || gen.starts_with("lut:"))
    {
        lut::Options options;
        if (gen != "lut")
            options.template_order = std::stoi(gen.substr(4));
        return part::quadrature_rules(part, q, false, options);
    }
#ifdef CUTCELLS_WITH_ALGOIM
    if (gen == "algoim")
        return benchmarks::algoim_rules(part, q, false, false);
#endif
    throw std::runtime_error("generator " + gen + " is not available (algoim needs CUTCELLS_WITH_ALGOIM)");
}

/// With a P1 level set: the cut cells' measures of the parts, and of the disk
/// (ball) and circle (sphere), from quadrays with @p q points per segment.
void p1_reference(const SweepMesh& m, const part::CutResult<double, int>& result, int q, CutCells& cut)
{
    for (const bool surface : {false, true})
    {
        const part::MeshPart<double, int> part = part::select(result, surface ? "phi = 0" : "phi < 0");
        const quadrature::QuadratureRules<double> rules = part::quadrature_rules(part, q, false, quadrays::Options{});
        for (CellInfo& info : cut.infos)
            (surface ? info.exact_surface : info.exact_volume) = 0;
        double total = 0;
        for (std::size_t k = 0; k + 1 < rules._offset.size(); ++k)
        {
            CellInfo& info = cut.infos[static_cast<std::size_t>(cut.where[static_cast<std::size_t>(rules._parent_map[k])])];
            for (int p = rules._offset[k]; p < rules._offset[k + 1]; ++p)
            {
                (surface ? info.exact_surface : info.exact_volume) += rules._weights[static_cast<std::size_t>(p)];
                total += rules._weights[static_cast<std::size_t>(p)];
            }
        }
        if (!surface)
            for (int cell = 0; cell < result.num_cells; ++cell)
                if (result.domains[static_cast<std::size_t>(cell)] == cell::domain::inside)
                    total += m.cell_measure;
        (surface ? cut.whole_surface : cut.whole_volume) = total;
    }
}

/// Counters of one step, and their sums and maxima over the sweep.
struct Metrics
{
    long steps = 0, cut = 0, points = 0;
    long negative = 0, nonfinite = 0, outside_cell = 0, outside_part = 0, failed = 0;
    double total_error = 0; ///< max over steps, relative
    double l1 = 0;          ///< max over steps, relative
    double worst = 0;       ///< max over steps and cells, relative to h^d or h^(d-1)
    double seconds = 0;
    long compressed_points = 0;
    double compress_seconds = 0, residual = 0;
};

void merge(Metrics& sum, const Metrics& s)
{
    sum.steps += s.steps;
    sum.cut += s.cut;
    sum.points += s.points;
    sum.negative += s.negative;
    sum.nonfinite += s.nonfinite;
    sum.outside_cell += s.outside_cell;
    sum.outside_part += s.outside_part;
    sum.failed += s.failed;
    sum.total_error = std::max(sum.total_error, s.total_error);
    sum.l1 = std::max(sum.l1, s.l1);
    sum.worst = std::max(sum.worst, s.worst);
    sum.seconds += s.seconds;
    sum.compressed_points += s.compressed_points;
    sum.compress_seconds += s.compress_seconds;
    sum.residual = std::max(sum.residual, s.residual);
}

/// One generator on one part at one centre.
Metrics run_step(const Config& cfg, const SweepMesh& m, const part::CutResult<double, int>& result,
                 const CutCells& cut, const std::string& part_text, const std::string& gen, int q, const V3& c)
{
    const double r = cfg.radius;
    const bool surface = part_text.find('=') != std::string::npos;
    const int d = m.tdim;
    const double whole = surface ? cut.whole_surface : cut.whole_volume;
    const double unit = std::pow(m.h, surface ? d - 1 : d);
    Metrics s;
    s.steps = 1;
    s.cut = static_cast<long>(result.cut_cells.size());

    const part::MeshPart<double, int> part = part::select(result, part_text);
    quadrature::QuadratureRules<double> rules;
    const auto t0 = std::chrono::steady_clock::now();
    try
    {
        rules = rules_of(gen, part, q);
    }
    catch (const std::exception& e)
    {
        s.failed = 1;
        std::fprintf(stderr, "%s q = %d %s at c = (%g, %g, %g): %s\n", gen.c_str(), q, part_text.c_str(), c[0],
                     c[1], c[2], e.what());
        return s;
    }
    s.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    s.points = static_cast<long>(rules._weights.size());

    // per cut cell: the generated measure against the exact one
    std::vector<double> generated(static_cast<std::size_t>(m.view.num_cells()), 0.0);
    const int tdim = rules._tdim > 0 ? rules._tdim : d;
    LevelSetCell<double, int> ls_cell;
    for (std::size_t k = 0; k + 1 < rules._offset.size(); ++k)
    {
        const int cell = rules._parent_map[k];
        const int at = cut.where[static_cast<std::size_t>(cell)];
        if (at < 0)
            throw std::runtime_error("a rule on a cell that is not cut");
        const CellInfo& info = cut.infos[static_cast<std::size_t>(at)];
        if (cfg.degree == 1)
            make_cell_level_set(*result.level_sets[0], cell, ls_cell);
        for (int p = rules._offset[k]; p < rules._offset[k + 1]; ++p)
        {
            const double w = rules._weights[static_cast<std::size_t>(p)];
            const std::span<const double> xi(rules._points.data() + static_cast<std::size_t>(tdim) * p,
                                              static_cast<std::size_t>(tdim));
            bool finite = std::isfinite(w);
            for (const double x : xi)
                finite &= std::isfinite(x);
            if (!finite)
            {
                ++s.nonfinite;
                continue;
            }
            generated[static_cast<std::size_t>(cell)] += w;
            s.negative += w < 0;
            s.outside_cell += !in_reference_cell(info.type, xi, 1e-12);
            const V3 x = info.physical_point(xi);
            double dist2 = 0;
            for (int j = 0; j < d; ++j)
                dist2 += (x[j] - c[j]) * (x[j] - c[j]);
            // the P1 level set: its gradient is about that of |x - c|^2 - r^2 on the circle (sphere)
            const double phi = cfg.degree == 1 ? ls_cell.value(xi) : dist2 - r * r;
            const double reach = 1e-8 * m.h * 2 * (cfg.degree == 1 ? r : std::sqrt(dist2));
            s.outside_part += surface ? std::abs(phi) > reach : phi > reach;
        }
    }
    double total = 0, cut_error = 0;
    for (std::size_t i = 0; i < result.cut_cells.size(); ++i)
    {
        const int cell = result.cut_cells[i];
        const CellInfo& info = cut.infos[i];
        const double e = generated[static_cast<std::size_t>(cell)] - (surface ? info.exact_surface : info.exact_volume);
        total += generated[static_cast<std::size_t>(cell)];
        cut_error += std::abs(e);
        s.worst = std::max(s.worst, std::abs(e) / unit);
    }
    // whole cells inside the disk (ball)
    if (!surface)
        for (int cell = 0; cell < result.num_cells; ++cell)
            if (result.domains[static_cast<std::size_t>(cell)] == cell::domain::inside)
                total += m.cell_measure;
    s.total_error = std::abs(total - whole) / whole;
    s.l1 = cut_error / whole;

    if (cfg.compress >= 0)
    {
        // Q_p holds a prism's P_p x P_p; pyramids take Q_p too
        const compression::MomentSpace space = m.kind == "tri" || m.kind == "tet" ? compression::MomentSpace::total
                                                                                 : compression::MomentSpace::tensor;
        quadrature::QuadratureRules<double> compressed;
        compression::CompressionStats stats;
        const auto t1 = std::chrono::steady_clock::now();
        compression::compress_rules(rules, cfg.compress, space, compressed, stats);
        s.compress_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
        s.compressed_points = static_cast<long>(compressed._weights.size());
        s.residual = stats.max_residual;
    }
    return s;
}
} // namespace

int main(int argc, char** argv)
{
    Config cfg;
    try
    {
        cfg = parse_args(argc, argv);
    }
    catch (const std::exception& e)
    {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }
    SweepMesh m;
    build_mesh(cfg, m);
    const double r = cfg.radius;
    const V3 dir = {1.0, 0.5 * (std::sqrt(5.0) - 1.0), std::sqrt(2.0) - 1.0};
    const char* shape = m.tdim == 2 ? "circle" : "sphere";
    std::printf("P%d level set%s. ", cfg.degree, cfg.tilt != 0 ? (" (tilt " + std::to_string(cfg.tilt) + ")").c_str() : "");
    std::printf("Sweeping %s of radius %g on %s, n = %d (h = %g, %d cells): %d centres c = c0 + s h (1, 0.618, "
                "0.414), s = k / %d, c0 = (%g, %g, %g)\n",
                shape, r, cfg.mesh.c_str(), cfg.n, m.h, m.view.num_cells(), cfg.steps, cfg.steps, cfg.centre[0],
                cfg.centre[1], m.tdim == 3 ? cfg.centre[2] : 0.0);
    std::printf("Errors: total and per-cell L1 relative to the %s, worst cell relative to h^%s; max over the sweep.\n",
                m.tdim == 2 ? "disk area (circle length)" : "ball volume (sphere area)",
                m.tdim == 2 ? "2 (h)" : "3 (h^2)");

    std::ofstream csv;
    if (!cfg.csv.empty())
    {
        csv.open(cfg.csv);
        csv << "mesh,degree,n,q,generator,part,step,s,cut_cells,points,total_error,l1,worst,negative,nonfinite,"
               "outside_cell,outside_part,failed,seconds,compressed_points,compress_seconds\n";
    }

    const std::vector<std::string> parts = {"phi < 0", "phi = 0"};
    // [q][generator][part]
    std::vector<std::vector<std::vector<Metrics>>> sums(
        cfg.qs.size(), std::vector<std::vector<Metrics>>(cfg.generators.size(), std::vector<Metrics>(parts.size())));
    const LevelSetMeshData<double, int> layout = create_level_set_mesh_data<double, int>(m.view, cfg.degree);
    // k = -1: a warm-up run at the first centre, not counted (first-call setup)
    for (int k = -1; k < cfg.steps; ++k)
    {
        const double s = static_cast<double>(std::max(k, 0)) / cfg.steps;
        V3 c = cfg.centre;
        for (int j = 0; j < m.tdim; ++j)
            c[j] += s * m.h * dir[j];
        if (m.tdim == 2)
            c[2] = 0;
        LevelSetMeshData<double, int> data = layout;
        std::vector<double> values(static_cast<std::size_t>(data.num_dofs()));
        for (int dof = 0; dof < data.num_dofs(); ++dof)
        {
            const double* x = data.dof_coordinate(dof);
            double v = -r * r;
            for (int j = 0; j < m.tdim; ++j)
            {
                v += (x[j] - c[j]) * (x[j] - c[j]);
                for (int i = 0; i < j; ++i)
                    v += cfg.tilt * (x[i] - c[i]) * (x[j] - c[j]);
            }
            values[static_cast<std::size_t>(dof)] = v;
        }
        const std::vector<LevelSetFunction<double, int>> ls
            = {create_level_set_function<double, int>(std::move(data), values, "phi")};
        const part::CutResult<double, int> result = part::cut<double, int>(m.view, ls);
        CutCells cut;
        cut.where.assign(static_cast<std::size_t>(result.num_cells), -1);
        for (const int cell : result.cut_cells)
        {
            cut.where[static_cast<std::size_t>(cell)] = static_cast<int>(cut.infos.size());
            cut.infos.push_back(cell_info(m, cell, c, r));
        }
        if (cfg.degree == 1)
            p1_reference(m, result, cfg.ref_q, cut);
        else
        {
            cut.whole_volume = m.tdim == 2 ? M_PI * r * r : 4.0 / 3.0 * M_PI * r * r * r;
            cut.whole_surface = m.tdim == 2 ? 2 * M_PI * r : 4 * M_PI * r * r;
        }
        for (std::size_t iq = 0; iq < cfg.qs.size(); ++iq)
            for (std::size_t g = 0; g < cfg.generators.size(); ++g)
                for (std::size_t p = 0; p < parts.size(); ++p)
                {
                    const Metrics st = run_step(cfg, m, result, cut, parts[p], cfg.generators[g], cfg.qs[iq], c);
                    if (k < 0)
                        continue;
                    merge(sums[iq][g][p], st);
                    if (csv)
                        csv << cfg.mesh << ',' << cfg.degree << ',' << cfg.n << ',' << cfg.qs[iq] << ',' << cfg.generators[g] << ','
                            << parts[p] << ',' << k << ',' << s << ',' << st.cut << ',' << st.points << ','
                            << st.total_error << ',' << st.l1 << ',' << st.worst << ',' << st.negative << ','
                            << st.nonfinite << ',' << st.outside_cell << ',' << st.outside_part << ',' << st.failed
                            << ',' << st.seconds << ',' << st.compressed_points << ',' << st.compress_seconds
                            << '\n';
                }
    }

    std::printf("%-4s %2s %-9s %-8s | %6s | %7s | %8s %8s %8s | %5s %5s %5s %5s %4s | %8s", "mesh", "q", "generator",
                "part", "cut", "pts/cut", "total", "L1", "worst", "neg", "nan", "cell", "part", "fail", "us/cut");
    if (cfg.compress >= 0)
        std::printf(" | p = %d: %7s %8s %8s", cfg.compress, "pts/cut", "us/cut", "residual");
    std::printf("\n");
    for (std::size_t iq = 0; iq < cfg.qs.size(); ++iq)
        for (std::size_t g = 0; g < cfg.generators.size(); ++g)
            for (std::size_t p = 0; p < parts.size(); ++p)
            {
                const Metrics& t = sums[iq][g][p];
                const double cut = static_cast<double>(std::max<long>(t.cut, 1));
                std::printf("%-4s %2d %-9s %-8s | %6.0f | %7.1f | %8.1e %8.1e %8.1e | %5ld %5ld %5ld %5ld %4ld | %8.1f",
                            cfg.mesh.c_str(), cfg.qs[iq], cfg.generators[g].c_str(), parts[p].c_str(),
                            static_cast<double>(t.cut) / std::max<long>(t.steps, 1), t.points / cut, t.total_error,
                            t.l1, t.worst, t.negative, t.nonfinite, t.outside_cell, t.outside_part, t.failed,
                            1e6 * t.seconds / cut);
                if (cfg.compress >= 0)
                    std::printf(" | %13.1f %8.1f %8.1e", t.compressed_points / cut, 1e6 * t.compress_seconds / cut,
                                t.residual);
                std::printf("\n");
            }
    return 0;
}
