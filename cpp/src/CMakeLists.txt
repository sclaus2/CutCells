# ------------------------------------------------------------------------------
include(GNUInstallDirs)

# ------------------------------------------------------------------------------
# Declare the library (target)
add_library(cutcells)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.h.in version.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_types.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_flags.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_triangle.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_interval.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_tetrahedron.h
  ${CMAKE_CURRENT_SOURCE_DIR}/triangulation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_cell.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_mesh.h
  ${CMAKE_CURRENT_SOURCE_DIR}/write_tikz.h
  ${CMAKE_CURRENT_SOURCE_DIR}/write_vtk.h
  PARENT_SCOPE)

target_sources(cutcells PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_interval.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_triangle.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_tetrahedron.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_cell.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_mesh.cpp
)



