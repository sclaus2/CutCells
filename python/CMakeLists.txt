cmake_minimum_required(VERSION 3.21)

project(cutcells_nanobind VERSION "0.5.0" LANGUAGES CXX)
include(GNUInstallDirs)

# nanobind uses aligned deallocators only present on macOS > 10.14
if(APPLE AND NOT DEFINED CMAKE_OSX_DEPLOYMENT_TARGET)
  set(CMAKE_OSX_DEPLOYMENT_TARGET "10.14")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Imported CutCells targets may depend on OpenMP::OpenMP_CXX.
# Keep this optional and quiet so builds without OpenMP still proceed.
find_package(OpenMP QUIET COMPONENTS CXX)

# Wheels and source distributions (the pyproject.toml in the repository root)
# build the C++ core into the extension module. The default is to use an
# installed CutCells CMake package, as downstream packages do.
option(CUTCELLS_PYTHON_BUNDLE_CORE
  "Build the C++ core from ../cpp into the extension instead of using an installed CutCells."
  OFF)

if(CUTCELLS_PYTHON_BUNDLE_CORE AND NOT TARGET CUTCELLS::cutcells)
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  # EXCLUDE_FROM_ALL keeps the core's install rules (headers, CMake package)
  # out of the wheel.
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../cpp" cutcells_core EXCLUDE_FROM_ALL)
  add_library(CUTCELLS::cutcells ALIAS cutcells)
endif()

if(NOT TARGET CUTCELLS::cutcells)
  find_package(CutCells 0.5.0 REQUIRED CONFIG)
endif()

if(CutCells_FOUND)
  message(STATUS "Found CutCells in: ${CutCells_DIR}")
endif()

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
message(STATUS "Using Python executable: ${Python_EXECUTABLE}")

# Detect the installed nanobind package and import it into CMake
execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import nanobind; print(nanobind.cmake_dir())"
  OUTPUT_VARIABLE NB_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_cutcellscpp cutcells/wrapper.cpp)
target_link_libraries(_cutcellscpp PRIVATE CUTCELLS::cutcells)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wall -Wextra -Wno-comment -pedantic" HAVE_PEDANTIC)
if(HAVE_PEDANTIC)
  target_compile_options(_cutcellscpp PRIVATE -Wall -Wextra -Wno-comment)
endif()

if(CUTCELLS_PYTHON_BUNDLE_CORE)
  # The core is linked statically into the module: nothing to find at run time.
  install(TARGETS _cutcellscpp LIBRARY DESTINATION cutcells)
else()
  set(_build_rpath "$<TARGET_FILE_DIR:CUTCELLS::cutcells>")
  if(APPLE)
    set(_install_rpath "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
  elseif(UNIX)
    set(_install_rpath "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
  else()
    set(_install_rpath "")
  endif()

  set_target_properties(_cutcellscpp PROPERTIES
    BUILD_RPATH "${_build_rpath}"
    INSTALL_RPATH "${_install_rpath}"
    BUILD_WITH_INSTALL_RPATH FALSE
    INSTALL_RPATH_USE_LINK_PATH FALSE
  )

  install(TARGETS _cutcellscpp LIBRARY DESTINATION cutcells)
  install(IMPORTED_RUNTIME_ARTIFACTS CUTCELLS::cutcells
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
