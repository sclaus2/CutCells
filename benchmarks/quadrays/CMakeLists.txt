# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT

# Comparison drivers for the quadrays backend. They share the exact
# references and test meshes in cpp/tests/quadrays/support/ with the tests.

if(CUTCELLS_WITH_ALGOIM)
  # ==========================================================================
  # algoim's multi-polynomial engine with three local changes, written to the
  # build tree so third_party/ stays untouched:
  #  1. the two-polynomial constructor with user masks does not compile
  #     (maskEmpty needs detail::); fixed here
  #  2. a constructor taking any number of polynomials with masks
  #  3. force_k[N], a per-dimension override of the height direction
  # ==========================================================================
  file(READ "${CUTCELLS_ALGOIM_DIR}/quadrature_multipoly.hpp" _mp)
  string(REPLACE "\r" "" _mp "${_mp}")

  function(_algoim_patch old new)
    string(FIND "${_mp}" "${old}" _pos)
    if(_pos EQUAL -1)
      message(FATAL_ERROR "algoim patch: text not found:\n${old}")
    endif()
    string(REPLACE "${old}" "${new}" _out "${_mp}")
    set(_mp "${_out}" PARENT_SCOPE)
  endfunction()

  _algoim_patch([==[                if (!maskEmpty(mask))]==]
                [==[                if (!detail::maskEmpty(mask))]==])

  _algoim_patch([==[    enum QuadStrategy { AlwaysGL, AlwaysTS, AutoMixed };]==]
                [==[    // [CutCells benchmarks] force the height direction chosen at dimension N (-1 = algoim's choice)
    inline thread_local int force_k[4] = {-1, -1, -1, -1};

    enum QuadStrategy { AlwaysGL, AlwaysTS, AutoMixed };]==])

  _algoim_patch([==[                k = argmax(score);
                detail::eliminate_axis(phi, k, base.phi);]==]
                [==[                k = argmax(score);
                if (force_k[N] >= 0)
                    k = force_k[N];
                detail::eliminate_axis(phi, k, base.phi);]==])

  _algoim_patch([==[        // Assuming phi has been instantiated, determine elimination axis and build base]==]
                [==[        // [CutCells benchmarks] any number of polynomials with user-defined masks
        ImplicitPolyQuadrature(int count, const xarray<real,N>* const* ps, const booluarray<N,ALGOIM_M>* masks)
        {
            for (int i = 0; i < count; ++i)
            {
                auto mask = detail::nonzeroMask(*ps[i], masks[i]);
                if (!detail::maskEmpty(mask))
                    phi.push_back(*ps[i], mask);
            }
            build(true, false);
        }

        // Assuming phi has been instantiated, determine elimination axis and build base]==])

  set(QUADRAYS_ALGOIM_EXT_DIR "${CMAKE_CURRENT_BINARY_DIR}/algoim_ext")
  file(WRITE "${QUADRAYS_ALGOIM_EXT_DIR}/quadrature_multipoly.hpp" "${_mp}")

  # The patched header must be found before the vendored one.
  add_library(quadrays_algoim INTERFACE)
  target_include_directories(quadrays_algoim BEFORE INTERFACE
    "${QUADRAYS_ALGOIM_EXT_DIR}")
  target_link_libraries(quadrays_algoim INTERFACE benchmarks_algoim)

  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_INCLUDES
      "${QUADRAYS_ALGOIM_EXT_DIR};${CUTCELLS_ALGOIM_DIR};${CUTCELLS_SOURCE_DIR}/third_party/lapacke_compat")
  set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})
  check_cxx_source_compiles([=[
    #include <unordered_map>
    #include "quadrature_multipoly.hpp"

    int main()
    {
      double c[8] = {-1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
      algoim::xarray<algoim::real, 3> phi(c, algoim::uvector<int, 3>(2, 2, 2));
      const algoim::xarray<algoim::real, 3>* polys[1] = {&phi};
      const algoim::booluarray<3, ALGOIM_M> masks[1]
          = {algoim::booluarray<3, ALGOIM_M>(true)};
      algoim::force_k[3] = 0;
      algoim::ImplicitPolyQuadrature<3> q(1, polys, masks);
      algoim::force_k[3] = -1;
      double sum = 0.0;
      q.integrate(algoim::AlwaysGL, 2,
                  [&](const auto&, algoim::real w) { sum += w; });
      return sum > 0.0 ? 0 : 1;
    }
  ]=] QUADRAYS_PATCHED_ALGOIM_COMPILES)
  unset(CMAKE_REQUIRED_INCLUDES)
  unset(CMAKE_REQUIRED_LIBRARIES)
  if(NOT QUADRAYS_PATCHED_ALGOIM_COMPILES)
    message(FATAL_ERROR
      "CUTCELLS_WITH_ALGOIM=ON but the patched quadrature_multipoly.hpp did "
      "not compile or link (see CMakeFiles/CMakeConfigureLog.yaml)")
  endif()
endif()

# ============================================================================
# Drivers
# ============================================================================
set(QUADRAYS_SUPPORT_DIR "${CUTCELLS_SOURCE_DIR}/cpp/tests/quadrays/support")

# quadrays_benchmark(<name> <source>): a driver linked to the installed CutCells
# and, with CUTCELLS_WITH_ALGOIM, to algoim's generators.
function(quadrays_benchmark name source)
  add_executable(${name} ${source})
  target_include_directories(${name} PRIVATE "${QUADRAYS_SUPPORT_DIR}")
  target_link_libraries(${name} PRIVATE CUTCELLS::cutcells)
  if(CUTCELLS_WITH_ALGOIM)
    target_sources(${name} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/generators.cpp")
    target_link_libraries(${name} PRIVATE quadrays_algoim)
  endif()
endfunction()

quadrays_benchmark(quadrays_study study.cpp)
quadrays_benchmark(quadrays_robustness_report robustness_report.cpp)

# a sweeping circle and sphere through the front end: quadrays, and algoim on
# the same parts with CUTCELLS_WITH_ALGOIM
add_executable(quadrays_sweep sweep.cpp)
target_include_directories(quadrays_sweep PRIVATE "${QUADRAYS_SUPPORT_DIR}")
target_link_libraries(quadrays_sweep PRIVATE CUTCELLS::cutcells)
if(CUTCELLS_WITH_ALGOIM)
  target_sources(quadrays_sweep PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/mesh_part_algoim.cpp")
  target_link_libraries(quadrays_sweep PRIVATE quadrays_algoim)
endif()

# algoim's quadrature on the parts of the front end (formerly the library's
# backend='algoim'), and its checks
if(CUTCELLS_WITH_ALGOIM)
  add_executable(mesh_part_algoim_check mesh_part_algoim_check.cpp mesh_part_algoim.cpp)
  target_link_libraries(mesh_part_algoim_check PRIVATE CUTCELLS::cutcells quadrays_algoim)
  add_test(NAME mesh_part_algoim_check COMMAND mesh_part_algoim_check)
endif()
