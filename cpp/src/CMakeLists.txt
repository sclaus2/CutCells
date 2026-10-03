include(GNUInstallDirs)

add_library(cutcells)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.h.in version.h)

target_include_directories(cutcells PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

# Sub-folders add their sources to the cutcells target and return their
# headers (installed under include/cutcells/<folder>/).
add_subdirectory(quadrays)
add_subdirectory(lut)
add_subdirectory(part)
add_subdirectory(compression)

set(HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_types.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_flags.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_topology.h
  ${CMAKE_CURRENT_SOURCE_DIR}/reference_cell.h
  ${CMAKE_CURRENT_SOURCE_DIR}/mapping.h
  ${CMAKE_CURRENT_SOURCE_DIR}/bernstein.h
  ${CMAKE_CURRENT_SOURCE_DIR}/level_set.h
  ${CMAKE_CURRENT_SOURCE_DIR}/level_set_cell.h
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh_view.h
  ${CMAKE_CURRENT_SOURCE_DIR}/selection_expr.h
  ${CMAKE_CURRENT_SOURCE_DIR}/quadrature.h
  ${CMAKE_CURRENT_SOURCE_DIR}/quadrature_tables.h
  ${CMAKE_CURRENT_SOURCE_DIR}/row_view.h
  ${CMAKE_CURRENT_SOURCE_DIR}/span_math.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utils.h
  ${CMAKE_CURRENT_SOURCE_DIR}/write_tikz.h
  ${CMAKE_CURRENT_SOURCE_DIR}/write_vtk.h
  ${QUADRAYS_HEADERS}
  ${LUT_HEADERS}
  ${PART_HEADERS}
  ${COMPRESSION_HEADERS}
  PARENT_SCOPE)

target_sources(cutcells PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/mapping.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bernstein.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/level_set.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/level_set_cell.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/selection_expr.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/quadrature.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/write_vtk.cpp
)

target_compile_options(cutcells PRIVATE
  $<$<CONFIG:Release>:-O3>)
