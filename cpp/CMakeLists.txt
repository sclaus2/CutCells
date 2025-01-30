cmake_minimum_required(VERSION 3.21)

# ------------------------------------------------------------------------------
# Set project name and version number
project(CutCells VERSION "0.2.0" LANGUAGES CXX C)

# Use C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options
option(BUILD_SHARED_LIBS "Build CutCells with shared libraries." ON)
#add_feature_info(BUILD_SHARED_LIBS BUILD_SHARED_LIBS "Build cutfemx with shared libraries.")

include(FeatureSummary)

# Source files
add_subdirectory(src)

set_target_properties(cutcells PROPERTIES PRIVATE_HEADER "${HEADERS}")
target_include_directories(cutcells PUBLIC
  $<INSTALL_INTERFACE:include>
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_SOURCE_DIR}>")


install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cutcells/${DIR}
COMPONENT Development)

# Install the cutcells library
install(TARGETS cutcells
  EXPORT CutCellsTargets
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/src
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT RuntimeExecutables
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

# Configure CMake helpers
include(CMakePackageConfigHelpers)
write_basic_package_version_file(CutCellsConfigVersion.cmake VERSION ${PACKAGE_VERSION}
  COMPATIBILITY AnyNewerVersion)
configure_package_config_file(CutCellsConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/CutCellsConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cutcells)

# Install CMake files
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CutCellsConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/CutCellsConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cutcells COMPONENT Development)
install(EXPORT CutCellsTargets FILE CutCellsTargets.cmake NAMESPACE CUTCELLS::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cutcells)


