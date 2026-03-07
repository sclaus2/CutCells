#include <benchmark/benchmark.h>

#include <array>
#include <cmath>
#include <limits>
#include <random>

#include "cut_cell.h"
#include "cut_hexahedron.h"

namespace
{
constexpr int kRepeats = 1'000'000;

std::array<double, 24> make_random_hexahedron(std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> shift(-10.0, 10.0);
    std::uniform_real_distribution<double> len(0.1, 2.0);

    const double x0 = shift(rng);
    const double y0 = shift(rng);
    const double z0 = shift(rng);
    const double lx = len(rng);
    const double ly = len(rng);
    const double lz = len(rng);

    return {
        x0, y0, z0,
        x0 + lx, y0, z0,
        x0 + lx, y0 + ly, z0,
        x0, y0 + ly, z0,
        x0, y0, z0 + lz,
        x0 + lx, y0, z0 + lz,
        x0 + lx, y0 + ly, z0 + lz,
        x0, y0 + ly, z0 + lz,
    };
}

std::array<double, 8> make_intersecting_ls(const std::array<double, 24>& coords, std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> w(-1.0, 1.0);
    std::uniform_real_distribution<double> d(-1.0, 1.0);

    std::array<double, 8> ls{};
    for (;;)
    {
        const double nx = w(rng);
        const double ny = w(rng);
        const double nz = w(rng);
        const double shift = d(rng);

        bool has_pos = false;
        bool has_neg = false;
        bool near_zero = false;

        for (int i = 0; i < 8; ++i)
        {
            const double x = coords[i * 3 + 0];
            const double y = coords[i * 3 + 1];
            const double z = coords[i * 3 + 2];
            ls[i] = nx * x + ny * y + nz * z + shift;

            has_pos = has_pos || (ls[i] > 0.0);
            has_neg = has_neg || (ls[i] < 0.0);
            near_zero = near_zero || (std::abs(ls[i]) < std::numeric_limits<double>::epsilon());
        }

        if (has_pos && has_neg && !near_zero)
            return ls;
    }
}
} // namespace

static void BM_CutHexahedron(benchmark::State& state)
{
    std::mt19937_64 rng(20260307);

    const auto coords = make_random_hexahedron(rng);
    const auto ls = make_intersecting_ls(coords, rng);

    for (auto _ : state)
    {
        cutcells::cell::CutCell<double> out;
        for (int i = 0; i < kRepeats; ++i)
        {
            cutcells::cell::hexahedron::cut<double>(coords, 3, ls, "phi=0", out, false);
            benchmark::DoNotOptimize(out);
        }
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) * kRepeats);
}

BENCHMARK(BM_CutHexahedron);
BENCHMARK_MAIN();
