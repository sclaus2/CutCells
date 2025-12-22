cmake_minimum_required(VERSION 3.19)

# nanobind uses aligned deallocators only present on macOS > 10.14
if(APPLE)
  set(CMAKE_OSX_DEPLOYMENT_TARGET "10.14")
endif()

project(cutcells_nanobind VERSION "0.2.0" LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (TARGET cutcells)
  add_library(CUTCELLS::cutcells ALIAS cutcells)
else()
  # Find CutCells (C++)
  # Prefer the in-repo build if present, so `pip install -e python` uses local C++ changes.
  if (NOT DEFINED CutCells_DIR)
    set(_local_cutcells_prefix "${CMAKE_CURRENT_LIST_DIR}/../cpp/build-dir/prefix")
    if (EXISTS "${_local_cutcells_prefix}/lib/cmake/cutcells/CutCellsConfig.cmake")
      list(PREPEND CMAKE_PREFIX_PATH "${_local_cutcells_prefix}")
    endif()
  endif()
  find_package(CutCells REQUIRED)
endif()

if(CUTCELLS_FOUND)
  message(STATUS "Found CutCells at ${CUTCELLS_DIR}")
endif()

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# Detect the installed nanobind package and import it into CMake
execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import nanobind; print(nanobind.cmake_dir())"
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

# Create the binding library
nanobind_add_module(_cutcellscpp cutcells/wrapper.cpp)
target_link_libraries(_cutcellscpp PRIVATE CUTCELLS::cutcells)

# Add strict compiler flags
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-Wall -Werror -Wextra -Wno-comment -pedantic" HAVE_PEDANTIC)

if(HAVE_PEDANTIC)
  target_compile_options(_cutcellscpp PRIVATE -Wall;-Wextra;-Werror;-Wno-comment)
endif()

get_target_property(_location CUTCELLS::cutcells LOCATION)
get_filename_component(_cutcells_dir ${_location} DIRECTORY)

# Ensure we load the CutCells library we just found/linked against.
# In conda environments it's common to also have a (possibly older) libcutcells.dylib
# available on the default RPATH; if that appears before our in-repo prefix, dyld
# will pick the wrong one at runtime.
set_target_properties(_cutcellscpp PROPERTIES
  INSTALL_RPATH "${_cutcells_dir}"
  BUILD_RPATH "${_cutcells_dir}"
  BUILD_WITH_INSTALL_RPATH TRUE
  INSTALL_RPATH_USE_LINK_PATH FALSE)

if(APPLE AND DEFINED ENV{CONDA_PREFIX})
  # scikit-build/nanobind often injects ${CONDA_PREFIX}/lib as an LC_RPATH ahead of ours.
  # If that environment contains another libcutcells.dylib, dyld will load the wrong one.
  # Reorder RPATHs post-link so the in-repo CutCells lib directory is searched first.
  set(_conda_lib "$ENV{CONDA_PREFIX}/lib")
  add_custom_command(TARGET _cutcellscpp POST_BUILD
    COMMAND install_name_tool -delete_rpath "${_conda_lib}" $<TARGET_FILE:_cutcellscpp> || true
    COMMAND install_name_tool -delete_rpath "${_cutcells_dir}" $<TARGET_FILE:_cutcellscpp> || true
    COMMAND install_name_tool -add_rpath "${_cutcells_dir}" $<TARGET_FILE:_cutcellscpp>
    COMMAND install_name_tool -add_rpath "${_conda_lib}" $<TARGET_FILE:_cutcellscpp>
    VERBATIM)
endif()

install(TARGETS _cutcellscpp LIBRARY DESTINATION cutcells)
