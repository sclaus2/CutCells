cmake_minimum_required(VERSION 3.21)

project(CutCells VERSION "0.4.0" LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_SHARED_LIBS "Build CutCells with shared libraries." ON)
option(CUTCELLS_WITH_ALGOIM
  "Enable the optional Algoim-backed quadrature implementation."
  OFF)
option(CUTCELLS_ENABLE_IPO
  "Enable interprocedural optimization in release builds when supported."
  OFF)
option(BUILD_TESTING "Build the C++ tests in cpp/tests (run with ctest)." OFF)

include(GNUInstallDirs)
include(CheckIPOSupported)

add_subdirectory(src)

if(CUTCELLS_WITH_ALGOIM)
  set(CUTCELLS_ALGOIM_DIR
      "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/algoim")
  set(CUTCELLS_LAPACKE_COMPAT_DIR
      "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/lapacke_compat")
  if(NOT EXISTS "${CUTCELLS_ALGOIM_DIR}/algoim/quadrature_general.hpp")
    message(FATAL_ERROR
      "CUTCELLS_WITH_ALGOIM=ON requires the vendored Algoim headers at "
      "${CUTCELLS_ALGOIM_DIR}/algoim")
  endif()
  if(NOT EXISTS "${CUTCELLS_ALGOIM_DIR}/algoim/quadrature_multipoly.hpp")
    message(FATAL_ERROR
      "CUTCELLS_WITH_ALGOIM=ON requires Algoim quadrature_multipoly.hpp")
  endif()
  if(NOT EXISTS "${CUTCELLS_LAPACKE_COMPAT_DIR}/lapacke.h")
    message(FATAL_ERROR
      "CUTCELLS_WITH_ALGOIM=ON requires the local LAPACKE compatibility header")
  endif()

  set(CUTCELLS_ALGOIM_LAPACK_LIBRARIES "" CACHE STRING
    "BLAS/LAPACK libraries for the Algoim-backed quadrature backend")
  if(NOT CUTCELLS_ALGOIM_LAPACK_LIBRARIES)
    set(_cutcells_lapack_hints)
    foreach(_cutcells_prefix IN LISTS CMAKE_PREFIX_PATH)
      list(APPEND _cutcells_lapack_hints "${_cutcells_prefix}")
    endforeach()

    if(_cutcells_lapack_hints)
      find_library(CUTCELLS_BLAS_LIBRARY
        NAMES blas openblas
        PATHS ${_cutcells_lapack_hints}
        PATH_SUFFIXES lib
        NO_DEFAULT_PATH)
      find_library(CUTCELLS_LAPACK_LIBRARY
        NAMES lapack
        PATHS ${_cutcells_lapack_hints}
        PATH_SUFFIXES lib
        NO_DEFAULT_PATH)
      find_library(CUTCELLS_LAPACKE_LIBRARY
        NAMES lapacke
        PATHS ${_cutcells_lapack_hints}
        PATH_SUFFIXES lib
        NO_DEFAULT_PATH)
    endif()

    if(CUTCELLS_BLAS_LIBRARY)
      list(APPEND CUTCELLS_ALGOIM_LAPACK_LIBRARIES
        "${CUTCELLS_BLAS_LIBRARY}")
    endif()
    if(CUTCELLS_LAPACKE_LIBRARY)
      list(APPEND CUTCELLS_ALGOIM_LAPACK_LIBRARIES
        "${CUTCELLS_LAPACKE_LIBRARY}")
    endif()
    if(CUTCELLS_LAPACK_LIBRARY)
      list(APPEND CUTCELLS_ALGOIM_LAPACK_LIBRARIES
        "${CUTCELLS_LAPACK_LIBRARY}")
    endif()
  endif()

  if(NOT CUTCELLS_ALGOIM_LAPACK_LIBRARIES)
    find_package(LAPACK REQUIRED)
    list(APPEND CUTCELLS_ALGOIM_LAPACK_LIBRARIES ${LAPACK_LIBRARIES})
  endif()
  if(NOT CUTCELLS_ALGOIM_LAPACK_LIBRARIES)
    message(FATAL_ERROR
      "CUTCELLS_WITH_ALGOIM=ON requires LAPACK/LAPACKE libraries")
  endif()

  include(CheckCXXSourceCompiles)
  set(_cutcells_saved_required_includes "${CMAKE_REQUIRED_INCLUDES}")
  set(_cutcells_saved_required_libraries "${CMAKE_REQUIRED_LIBRARIES}")
  set(CMAKE_REQUIRED_INCLUDES
      "${CUTCELLS_LAPACKE_COMPAT_DIR};${CUTCELLS_ALGOIM_DIR}")
  set(CMAKE_REQUIRED_LIBRARIES ${CUTCELLS_ALGOIM_LAPACK_LIBRARIES})
  check_cxx_source_compiles([=[
    #include <algoim/quadrature_general.hpp>

    struct Phi
    {
      template <typename T>
      T operator()(const algoim::uvector<T, 2>& x) const
      {
        return x(0) + x(1) - T(0.5);
      }

      template <typename T>
      algoim::uvector<T, 2> grad(const algoim::uvector<T, 2>&) const
      {
        return algoim::uvector<T, 2>(T(1), T(1));
      }
    };

    int main()
    {
      auto rule = algoim::quadGen<2>(
          Phi{}, algoim::HyperRectangle<double, 2>(0.0, 1.0),
          -1, -1, 2);
      return rule.nodes.empty() ? 1 : 0;
    }
  ]=] CUTCELLS_ALGOIM_GENERAL_COMPILES)
  if(NOT CUTCELLS_ALGOIM_GENERAL_COMPILES)
    message(FATAL_ERROR
      "CUTCELLS_WITH_ALGOIM=ON but Algoim quadrature_general.hpp did not compile")
  endif()
  check_cxx_source_compiles([=[
    #include <unordered_map>
    #include <algoim/quadrature_multipoly.hpp>

    int main()
    {
      double coeffs[4] = {-1.0, 1.0, 1.0, 1.0};
      algoim::xarray<algoim::real, 2> phi(
          coeffs, algoim::uvector<int, 2>(2, 2));
      algoim::ImplicitPolyQuadrature<2> q(phi);
      double sum = 0.0;
      q.integrate(algoim::AutoMixed, 2, [&](const auto&, algoim::real w) {
        sum += w;
      });
      return sum <= 0.0 ? 1 : 0;
    }
  ]=] CUTCELLS_ALGOIM_MULTIPOLY_COMPILES)
  set(CMAKE_REQUIRED_INCLUDES "${_cutcells_saved_required_includes}")
  set(CMAKE_REQUIRED_LIBRARIES "${_cutcells_saved_required_libraries}")
  if(NOT CUTCELLS_ALGOIM_MULTIPOLY_COMPILES)
    message(FATAL_ERROR
      "CUTCELLS_WITH_ALGOIM=ON but Algoim quadrature_multipoly.hpp did not compile/link")
  endif()

  target_include_directories(cutcells SYSTEM PRIVATE
    "${CUTCELLS_LAPACKE_COMPAT_DIR}"
    "${CUTCELLS_ALGOIM_DIR}")
  target_link_libraries(cutcells PRIVATE ${CUTCELLS_ALGOIM_LAPACK_LIBRARIES})
  target_compile_definitions(cutcells PRIVATE CUTCELLS_HAS_ALGOIM=1)
  target_compile_definitions(cutcells PRIVATE CUTCELLS_HAS_ALGOIM_MULTIPOLY=1)
endif()

set(CUTCELLS_BUILD_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
file(MAKE_DIRECTORY "${CUTCELLS_BUILD_INCLUDE_DIR}/cutcells")
# A header keeps its path below src/, in the build tree and when installed:
# src/quadrays/rules.h becomes <cutcells/quadrays/rules.h>.
foreach(_cutcells_header IN LISTS HEADERS)
  file(RELATIVE_PATH _cutcells_header_path
       "${CMAKE_CURRENT_SOURCE_DIR}/src" "${_cutcells_header}")
  get_filename_component(_cutcells_header_dir "${_cutcells_header_path}"
                         DIRECTORY)
  file(WRITE "${CUTCELLS_BUILD_INCLUDE_DIR}/cutcells/${_cutcells_header_path}"
       "#pragma once\n#include \"${_cutcells_header}\"\n")
  install(FILES "${_cutcells_header}"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cutcells/${_cutcells_header_dir}"
    COMPONENT Development)
endforeach()
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/src/generated"
  DESTINATION "${CUTCELLS_BUILD_INCLUDE_DIR}/cutcells"
  FILES_MATCHING PATTERN "*.h")

find_package(OpenMP QUIET COMPONENTS CXX)
if(OpenMP_CXX_FOUND AND TARGET OpenMP::OpenMP_CXX)
  target_link_libraries(cutcells PUBLIC OpenMP::OpenMP_CXX)
endif()

if(CUTCELLS_ENABLE_IPO)
  check_ipo_supported(RESULT CUTCELLS_IPO_SUPPORTED OUTPUT CUTCELLS_IPO_ERROR)
  if(CUTCELLS_IPO_SUPPORTED)
    set_property(TARGET cutcells PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
  else()
    message(STATUS "IPO requested but not supported: ${CUTCELLS_IPO_ERROR}")
  endif()
endif()

target_include_directories(cutcells PUBLIC
  $<BUILD_INTERFACE:${CUTCELLS_BUILD_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()

# Install generated headers produced at build time (tables, etc.)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/generated
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cutcells
  COMPONENT Development
  FILES_MATCHING PATTERN "*.h")

# Install the cutcells library
install(TARGETS cutcells
  EXPORT CutCellsTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT RuntimeExecutables
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

set(_cutcells_license_install_dir
  "${CMAKE_INSTALL_DATADIR}/licenses/cutcells")
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE"
  DESTINATION ${_cutcells_license_install_dir}
  RENAME CutCells-LICENSE.txt)
install(FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/THIRD_PARTY_NOTICES.md"
  "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/VTK-Copyright.txt"
  "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/Basix-LICENSE.txt"
  "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/algoim/UPSTREAM_COMMIT"
  "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/algoim/PATCHES.md"
  DESTINATION ${_cutcells_license_install_dir})
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/algoim/LICENSE"
  DESTINATION ${_cutcells_license_install_dir}
  RENAME Algoim-LICENSE.txt)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(CutCellsConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)
configure_package_config_file(CutCellsConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/CutCellsConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cutcells)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/CutCellsConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/CutCellsConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cutcells COMPONENT Development)
install(EXPORT CutCellsTargets
  FILE CutCellsTargets.cmake
  NAMESPACE CUTCELLS::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cutcells)
