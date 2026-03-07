# CutCells Benchmarks

This folder contains Google Benchmark executables for:

- `cut_triangle`
- `cut_tetrahedron`
- `cut_hexahedron`
- `cut_vtk_mesh`
- `create_cut_mesh`

Each benchmark uses `benchmark::State` and performs `1e6` cut operations per state iteration.

## Build (Release only)

Benchmarks are intentionally enabled only in Release mode.

```bash
cd /Users/sclaus/Workspace/FEniCSx0.10/CutCells/cpp
cmake -S . -B build-bench-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-bench-release --target bench_cut_triangle bench_cut_tetrahedron bench_cut_hex bench_cut_mesh -j
```

## Run individual benchmarks

```bash
./build-bench-release/benchmarks/bench_cut_triangle
./build-bench-release/benchmarks/bench_cut_tetrahedron
./build-bench-release/benchmarks/bench_cut_hex
./build-bench-release/benchmarks/bench_cut_mesh
```

## Run all benchmarks with one target

Use the aggregate target:

```bash
cmake --build build-bench-release --target run_benchmarks -j
```

## Useful Google Benchmark flags

Examples:

```bash
./build-bench-release/benchmarks/bench_cut_hex --benchmark_filter=CutHexahedron
./build-bench-release/benchmarks/bench_cut_mesh --benchmark_repetitions=5 --benchmark_report_aggregates_only=true
```
