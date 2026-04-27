include(GNUInstallDirs)

add_library(cutcells)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.h.in version.h)

# Make headers in this directory discoverable for generated table headers
# (e.g. generated/*_tables.h includes "cell_types.h").
target_include_directories(cutcells PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

set(HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_types.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_flags.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_subdivision.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_topology.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_triangle.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_interval.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_quadrilateral.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_hexahedron.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_prism.h
  ${CMAKE_CURRENT_SOURCE_DIR}/quad_midpoint_split.h
  ${CMAKE_CURRENT_SOURCE_DIR}/prism_midpoint_split.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_pyramid.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_tetrahedron.h
  ${CMAKE_CURRENT_SOURCE_DIR}/triangulation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_cell.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_mesh.h
  ${CMAKE_CURRENT_SOURCE_DIR}/mapping.h
  ${CMAKE_CURRENT_SOURCE_DIR}/bernstein.h
  ${CMAKE_CURRENT_SOURCE_DIR}/level_set.h
  ${CMAKE_CURRENT_SOURCE_DIR}/level_set_cell.h
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh_view.h
  ${CMAKE_CURRENT_SOURCE_DIR}/adapt_cell.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ho_cut_mesh.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ho_mesh_part_output.h
  ${CMAKE_CURRENT_SOURCE_DIR}/curving.h
  ${CMAKE_CURRENT_SOURCE_DIR}/selection_expr.h
  ${CMAKE_CURRENT_SOURCE_DIR}/edge_certification.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_certification.h
  ${CMAKE_CURRENT_SOURCE_DIR}/refine_cell.h
  ${CMAKE_CURRENT_SOURCE_DIR}/quadrature.h
  ${CMAKE_CURRENT_SOURCE_DIR}/quadrature_tables.h
  ${CMAKE_CURRENT_SOURCE_DIR}/span_math.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utils.h
  ${CMAKE_CURRENT_SOURCE_DIR}/write_tikz.h
  ${CMAKE_CURRENT_SOURCE_DIR}/write_vtk.h
  PARENT_SCOPE)

target_sources(cutcells PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_interval.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_triangle.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_quadrilateral.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_hexahedron.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_prism.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/quad_midpoint_split.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/prism_midpoint_split.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_pyramid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_tetrahedron.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_cell.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cut_mesh.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mapping.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bernstein.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/level_set.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/level_set_cell.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/adapt_cell.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ho_cut_mesh.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ho_mesh_part_output.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/curving.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/selection_expr.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/edge_certification.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_certification.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/refine_cell.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/quadrature.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/write_vtk.cpp
)

target_compile_options(cutcells PRIVATE
  $<$<CONFIG:Release>:-O3>)
