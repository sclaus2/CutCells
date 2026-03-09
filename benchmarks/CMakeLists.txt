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
  # On macOS, CMake embeds the Conda lib dir (from the transitive OpenMP
  # dependency) in the binary RPATH before the build-tree src path, causing
  # the OS to load the *installed* libcutcells instead of the freshly-built one.
  #
  # Post-build: run a cmake -P script that uses otool to find every
  # Conda/miniforge RPATH entry and strips it with install_name_tool.
  # The build-tree path (build-bench-release/src) that CMake already embeds
  # then becomes the first — and only — place the OS looks for libcutcells.
  if(APPLE)
    add_custom_command(TARGET ${target} POST_BUILD
      COMMAND ${CMAKE_COMMAND}
        -D "TARGET_FILE=$<TARGET_FILE:${target}>"
        -P "${CMAKE_CURRENT_LIST_DIR}/fix_benchmark_rpath.cmake"
      COMMENT "Stripping conda rpaths from ${target}"
    )
  endif()
endfunction()

add_cutcells_benchmark(bench_cut_triangle bench_cut_triangle.cpp)
add_cutcells_benchmark(bench_cut_tetrahedron bench_cut_tetrahedron.cpp)
add_cutcells_benchmark(bench_cut_hex bench_cut_hex.cpp)
add_cutcells_benchmark(bench_cut_mesh bench_cut_mesh.cpp)

# run_benchmarks: also set DYLD_LIBRARY_PATH as an extra safety net in case the
# rpath ordering is still overridden by the system linker.
add_custom_target(run_benchmarks
  COMMAND ${CMAKE_COMMAND} -E env "DYLD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/src"
    $<TARGET_FILE:bench_cut_triangle>
  COMMAND ${CMAKE_COMMAND} -E env "DYLD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/src"
    $<TARGET_FILE:bench_cut_tetrahedron>
  COMMAND ${CMAKE_COMMAND} -E env "DYLD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/src"
    $<TARGET_FILE:bench_cut_hex>
  COMMAND ${CMAKE_COMMAND} -E env "DYLD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/src"
    $<TARGET_FILE:bench_cut_mesh>
  DEPENDS bench_cut_triangle bench_cut_tetrahedron bench_cut_hex bench_cut_mesh
  USES_TERMINAL
  COMMENT "Build and run all CutCells benchmarks"
)
