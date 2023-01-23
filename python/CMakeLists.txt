cmake_minimum_required(VERSION 3.16)

PROJECT(cutcells_pybind11 VERSION "0.1"  LANGUAGES CXX)

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

# Correct COMPONENT specification Development.Module added only in CMake 3.18
# and above.
if(NOT SKBUILD)
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.18.0")
    find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
  else()
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
  endif()
endif()

execute_process(
  COMMAND
    "python3" -c
    "import pybind11; print(pybind11.get_cmake_dir())"
  OUTPUT_VARIABLE _tmp_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")

find_package(pybind11 REQUIRED CONFIG HINTS ${PYBIND11_DIR} ${PYBIND11_ROOT}
  $ENV{PYBIND11_DIR} $ENV{PYBIND11_ROOT})

# Create the binding library
pybind11_add_module(_cutcellscpp MODULE cutcells/wrapper.cpp)
target_link_libraries(_cutcellscpp PRIVATE CUTCELLS::cutcells)

# scikit-build specific steps
if (SKBUILD)
  # Add relative rpath so _basixcpp can find the Basix::basix library
  # when the build is relocated
  if (APPLE)
    set_target_properties(_cutcellscpp PROPERTIES INSTALL_RPATH "@loader_path/lib")
  else()
    set_target_properties(_cutcellscpp PROPERTIES INSTALL_RPATH "$ORIGIN/lib")
  endif()

  install(TARGETS _cutcellscpp DESTINATION .)
endif()
