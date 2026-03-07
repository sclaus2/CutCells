cmake_minimum_required(VERSION 3.21)

if(CMAKE_CONFIGURATION_TYPES)
  message(FATAL_ERROR "Benchmarks are supported only for single-config Release builds. Configure with -DCMAKE_BUILD_TYPE=Release.")
endif()

if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
  message(FATAL_ERROR "Benchmarks compile only in Release mode. Set -DCMAKE_BUILD_TYPE=Release.")
endif()

include(FetchContent)
find_package(benchmark QUIET)

if(NOT benchmark_FOUND)
  FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.9.1
  )
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(benchmark)
endif()

function(add_cutcells_benchmark target source_file)
  add_executable(${target} ${source_file})
  target_link_libraries(${target} PRIVATE cutcells benchmark::benchmark)
  target_compile_features(${target} PRIVATE cxx_std_20)
endfunction()

add_cutcells_benchmark(bench_cut_triangle bench_cut_triangle.cpp)
add_cutcells_benchmark(bench_cut_tetrahedron bench_cut_tetrahedron.cpp)
add_cutcells_benchmark(bench_cut_hex bench_cut_hex.cpp)
add_cutcells_benchmark(bench_cut_mesh bench_cut_mesh.cpp)

add_custom_target(run_benchmarks
  COMMAND $<TARGET_FILE:bench_cut_triangle>
  COMMAND $<TARGET_FILE:bench_cut_tetrahedron>
  COMMAND $<TARGET_FILE:bench_cut_hex>
  COMMAND $<TARGET_FILE:bench_cut_mesh>
  DEPENDS bench_cut_triangle bench_cut_tetrahedron bench_cut_hex bench_cut_mesh
  USES_TERMINAL
  COMMENT "Build and run all CutCells benchmarks"
)
