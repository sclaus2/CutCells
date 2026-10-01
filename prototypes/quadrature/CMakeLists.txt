# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.21)
project(cutcells_quadrature_prototype LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CUTCELLS_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../..)
set(ALGOIM_DIR ${CUTCELLS_ROOT}/third_party/algoim/algoim)

# ============================================================================
# algoim's multi-polynomial engine with three local changes, written to the
# build tree so third_party/ stays untouched:
#  1. the two-polynomial constructor with user masks does not compile
#     (maskEmpty needs detail::); fixed here
#  2. a constructor taking any number of polynomials with masks
#  3. force_k[N], a per-dimension override of the height direction
# ============================================================================
file(READ ${ALGOIM_DIR}/quadrature_multipoly.hpp _mp)
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
              [==[    // [CutCells prototype] force the height direction chosen at dimension N (-1 = algoim's choice)
    inline thread_local int force_k[4] = {-1, -1, -1, -1};

    enum QuadStrategy { AlwaysGL, AlwaysTS, AutoMixed };]==])

_algoim_patch([==[                k = argmax(score);
                detail::eliminate_axis(phi, k, base.phi);]==]
              [==[                k = argmax(score);
                if (force_k[N] >= 0)
                    k = force_k[N];
                detail::eliminate_axis(phi, k, base.phi);]==])

_algoim_patch([==[        // Assuming phi has been instantiated, determine elimination axis and build base]==]
              [==[        // [CutCells prototype] any number of polynomials with user-defined masks
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

file(WRITE ${CMAKE_BINARY_DIR}/algoim_ext/quadrature_multipoly.hpp "${_mp}")

# Prefer the active conda environment's liblapack/libblas over the system's Accelerate.
if(NOT DEFINED BLA_VENDOR AND DEFINED ENV{CONDA_PREFIX})
  set(BLA_VENDOR Generic)
endif()
find_package(LAPACK REQUIRED)

add_executable(quadrature_study
  src/study.cpp
  src/clipped_box.cpp
  src/generators.cpp
  ${CUTCELLS_ROOT}/cpp/src/selection_expr.cpp)

# The patched header must be found before the vendored one.
target_include_directories(quadrature_study PRIVATE
  ${CMAKE_BINARY_DIR}/algoim_ext
  ${ALGOIM_DIR}
  ${CUTCELLS_ROOT}/third_party/lapacke_compat
  ${CUTCELLS_ROOT}/cpp/src
  ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(quadrature_study PRIVATE ${LAPACK_LIBRARIES})

add_executable(test_exact_reference src/test_exact_reference.cpp)
target_include_directories(test_exact_reference PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(test_leaf_ordering src/test_leaf_ordering.cpp)
target_include_directories(test_leaf_ordering PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
