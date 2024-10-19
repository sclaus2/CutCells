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
  find_package(CutCells REQUIRED)
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
set_target_properties(_cutcellscpp PROPERTIES INSTALL_RPATH ${_cutcells_dir})

install(TARGETS _cutcellscpp LIBRARY DESTINATION cutcells)
