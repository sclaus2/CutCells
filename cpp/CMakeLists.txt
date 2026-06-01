cmake_minimum_required(VERSION 3.21)

project(CutCells VERSION "0.4.0" LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_SHARED_LIBS "Build CutCells with shared libraries." ON)

include(GNUInstallDirs)
include(CheckIPOSupported)

add_subdirectory(src)

set(CUTCELLS_BUILD_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
file(MAKE_DIRECTORY "${CUTCELLS_BUILD_INCLUDE_DIR}/cutcells")
foreach(_cutcells_header IN LISTS HEADERS)
  get_filename_component(_cutcells_header_name "${_cutcells_header}" NAME)
  file(WRITE "${CUTCELLS_BUILD_INCLUDE_DIR}/cutcells/${_cutcells_header_name}"
       "#pragma once\n#include \"${_cutcells_header}\"\n")
endforeach()
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/src/generated"
  DESTINATION "${CUTCELLS_BUILD_INCLUDE_DIR}/cutcells"
  FILES_MATCHING PATTERN "*.h")

find_package(OpenMP QUIET COMPONENTS CXX)
if(OpenMP_CXX_FOUND AND TARGET OpenMP::OpenMP_CXX)
  target_link_libraries(cutcells PUBLIC OpenMP::OpenMP_CXX)
endif()

check_ipo_supported(RESULT CUTCELLS_IPO_SUPPORTED OUTPUT CUTCELLS_IPO_ERROR)
if(CUTCELLS_IPO_SUPPORTED)
  set_property(TARGET cutcells PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

target_include_directories(cutcells PUBLIC
  $<BUILD_INTERFACE:${CUTCELLS_BUILD_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)


install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cutcells
  COMPONENT Development)

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
