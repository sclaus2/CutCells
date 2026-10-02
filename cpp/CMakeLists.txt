cmake_minimum_required(VERSION 3.21)

project(CutCells VERSION "0.4.0" LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_SHARED_LIBS "Build CutCells with shared libraries." ON)
option(CUTCELLS_ENABLE_IPO
  "Enable interprocedural optimization in release builds when supported."
  OFF)
option(BUILD_TESTING "Build the C++ tests in cpp/tests (run with ctest)." OFF)

include(GNUInstallDirs)
include(CheckIPOSupported)

add_subdirectory(src)

set(CUTCELLS_BUILD_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
# Rewritten on every configure, so headers that moved or went leave nothing behind
file(REMOVE_RECURSE "${CUTCELLS_BUILD_INCLUDE_DIR}/cutcells")
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
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/src/lut/generated"
  DESTINATION "${CUTCELLS_BUILD_INCLUDE_DIR}/cutcells/lut"
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

# Install the generated tables: quadrature rules, and the lookup tables
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/generated
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cutcells
  COMPONENT Development
  FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/lut/generated
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cutcells/lut
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
  DESTINATION ${_cutcells_license_install_dir})

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
