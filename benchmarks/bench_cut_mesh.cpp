#include <benchmark/benchmark.h>

#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

#include "cell_types.h"
#include "cut_hexahedron.h"
#include "cut_mesh.h"

namespace
{
constexpr int kRepeats = 1'000'000;

struct RandomHexCellData
{
    std::array<double, 24> points;
    std::array<int, 8> connectivity;
    std::array<int, 2> offset;
    std::array<int, 1> vtk_type;
    std::array<double, 8> ls;
};

RandomHexCellData make_random_hex_cell_data(std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> shift(-10.0, 10.0);
    std::uniform_real_distribution<double> len(0.1, 2.0);
    std::uniform_real_distribution<double> w(-1.0, 1.0);
    std::uniform_real_distribution<double> d(-1.0, 1.0);

    RandomHexCellData data{
        .points = {},
        .connectivity = {0, 1, 2, 3, 4, 5, 6, 7},
        .offset = {0, 8},
        .vtk_type = {static_cast<int>(cutcells::cell::vtk_types::VTK_HEXAHEDRON)},
        .ls = {},
    };

    const double x0 = shift(rng);
    const double y0 = shift(rng);
    const double z0 = shift(rng);
    const double lx = len(rng);
    const double ly = len(rng);
    const double lz = len(rng);

    data.points = {
        x0, y0, z0,
        x0 + lx, y0, z0,
        x0 + lx, y0 + ly, z0,
        x0, y0 + ly, z0,
        x0, y0, z0 + lz,
        x0 + lx, y0, z0 + lz,
        x0 + lx, y0 + ly, z0 + lz,
        x0, y0 + ly, z0 + lz,
    };

    for (;;)
    {
        const double nx = w(rng);
        const double ny = w(rng);
        const double nz = w(rng);
        const double shift_plane = d(rng);

        bool has_pos = false;
        bool has_neg = false;
        bool near_zero = false;

        for (int i = 0; i < 8; ++i)
        {
            const double x = data.points[i * 3 + 0];
            const double y = data.points[i * 3 + 1];
            const double z = data.points[i * 3 + 2];
            data.ls[i] = nx * x + ny * y + nz * z + shift_plane;

            has_pos = has_pos || (data.ls[i] > 0.0);
            has_neg = has_neg || (data.ls[i] < 0.0);
            near_zero = near_zero || (std::abs(data.ls[i]) < std::numeric_limits<double>::epsilon());
        }

        if (has_pos && has_neg && !near_zero)
            return data;
    }
}

cutcells::mesh::CutCells<double> make_random_cut_cells(std::mt19937_64& rng)
{
    const auto data = make_random_hex_cell_data(rng);

    cutcells::cell::CutCell<double> cut_cell;
    cutcells::cell::hexahedron::cut<double>(data.points, 3, data.ls, "phi=0", cut_cell, false);

    cut_cell._parent_cell_type = cutcells::cell::type::hexahedron;
    cut_cell._parent_vertex_ids = std::vector<int>(data.connectivity.begin(), data.connectivity.end());

    cutcells::mesh::CutCells<double> cells;
    cells._cut_cells = {std::move(cut_cell)};
    cells._parent_map = {0};
    cells._types = {cutcells::cell::type::hexahedron};
    return cells;
}
} // namespace

static void BM_CutVtkMesh(benchmark::State& state)
{
    std::mt19937_64 rng(20260307);
    const auto data = make_random_hex_cell_data(rng);

    for (auto _ : state)
    {
        cutcells::mesh::CutMesh<double> out;
        for (int i = 0; i < kRepeats; ++i)
        {
            out = cutcells::mesh::cut_vtk_mesh<double>(data.ls, data.points, data.connectivity, data.offset,
                                                       data.vtk_type, "phi=0", false);
            benchmark::DoNotOptimize(out);
        }
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) * kRepeats);
}

static void BM_CreateCutMesh(benchmark::State& state)
{
    std::mt19937_64 rng(20260307);

    for (auto _ : state)
    {
        for (int i = 0; i < kRepeats; ++i)
        {
            auto cells = make_random_cut_cells(rng);
            auto out = cutcells::mesh::create_cut_mesh<double>(cells);
            benchmark::DoNotOptimize(out);
        }
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) * kRepeats);
}

BENCHMARK(BM_CutVtkMesh);
BENCHMARK(BM_CreateCutMesh);
BENCHMARK_MAIN();
