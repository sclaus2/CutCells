cmake_minimum_required(VERSION 3.19)

project(cutcells_nanobind VERSION "0.2.0" LANGUAGES CXX)

if(DEFINED ENV{CONDA_PREFIX} AND NOT "$ENV{CONDA_PREFIX}" STREQUAL "")
  list(PREPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(Python_FIND_VIRTUALENV ONLY)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import sys; print(sys.prefix)"
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE PY_PREFIX
)
if(PY_PREFIX)
  list(PREPEND CMAKE_PREFIX_PATH "${PY_PREFIX}")
endif()

find_package(OpenMP REQUIRED)
find_package(
  CutCells CONFIG REQUIRED
  HINTS "${PY_PREFIX}" "$ENV{CONDA_PREFIX}"
  PATH_SUFFIXES "lib/cmake/cutcells"
)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import nanobind; print(nanobind.cmake_dir())"
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE NB_DIR
)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_cutcellscpp cutcells/wrapper.cpp)
target_link_libraries(_cutcellscpp PRIVATE CUTCELLS::cutcells)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wall -Wextra -Wno-comment" HAVE_WARN_FLAGS)
if(HAVE_WARN_FLAGS)
  target_compile_options(_cutcellscpp PRIVATE -Wall -Wextra -Wno-comment)
endif()

install(TARGETS _cutcellscpp LIBRARY DESTINATION cutcells)
