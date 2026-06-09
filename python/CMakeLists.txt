cmake_minimum_required(VERSION 3.19)

project(cutcells_nanobind VERSION "0.4.0" LANGUAGES CXX)

# nanobind uses aligned deallocators only present on macOS > 10.14
if(APPLE AND NOT DEFINED CMAKE_OSX_DEPLOYMENT_TARGET)
  set(CMAKE_OSX_DEPLOYMENT_TARGET "10.14")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(_conda_prefix "$ENV{CONDA_PREFIX}")
if(_conda_prefix)
  # Prefer the active conda environment over the base installation.
  if(NOT DEFINED Python_ROOT_DIR)
    set(Python_ROOT_DIR "${_conda_prefix}")
  endif()
  if(NOT DEFINED Python_FIND_VIRTUALENV)
    set(Python_FIND_VIRTUALENV ONLY)
  endif()
  if(NOT DEFINED Python_FIND_STRATEGY)
    set(Python_FIND_STRATEGY LOCATION)
  endif()
  if(NOT DEFINED Python_EXECUTABLE)
    if(WIN32)
      set(_conda_python "${_conda_prefix}/python.exe")
    else()
      set(_conda_python "${_conda_prefix}/bin/python")
    endif()
    if(EXISTS "${_conda_python}")
      set(Python_EXECUTABLE "${_conda_python}")
    endif()
  endif()

  list(PREPEND CMAKE_PREFIX_PATH
    "${_conda_prefix}"
    "${_conda_prefix}/lib/cmake"
    "${_conda_prefix}/Library"
    "${_conda_prefix}/Library/lib/cmake")
endif()

# Imported CutCells targets may depend on OpenMP::OpenMP_CXX.
# Keep this optional and quiet so builds without OpenMP still proceed.
find_package(OpenMP QUIET COMPONENTS CXX)

option(CUTCELLS_PYTHON_USE_LOCAL_CPP
  "Build the bundled CutCells C++ library from source with the same compiler as the Python extension."
  ON)

if(CUTCELLS_PYTHON_USE_LOCAL_CPP
   AND EXISTS "${CMAKE_CURRENT_LIST_DIR}/../cpp/CMakeLists.txt"
   AND NOT DEFINED CutCells_DIR
   AND NOT TARGET CUTCELLS::cutcells
   AND NOT TARGET cutcells)
  add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../cpp"
                   "${CMAKE_CURRENT_BINARY_DIR}/cutcells_cpp")
endif()

if(TARGET cutcells)
  add_library(CUTCELLS::cutcells ALIAS cutcells)
  if(APPLE)
    set_target_properties(cutcells PROPERTIES
      INSTALL_NAME_DIR "@loader_path/../lib")
  endif()
elseif(NOT TARGET CUTCELLS::cutcells)
  find_package(CutCells 0.4.0 REQUIRED CONFIG)
endif()

if(CutCells_FOUND)
  message(STATUS "Found CutCells in: ${CutCells_DIR}")
endif()

set(_cutcells_lib "")
if(TARGET cutcells)
  set(_cutcells_lib "$<TARGET_FILE:cutcells>")
elseif(TARGET CUTCELLS::cutcells)
  get_target_property(_cutcells_lib_rel CUTCELLS::cutcells IMPORTED_LOCATION_RELEASE)
  get_target_property(_cutcells_lib_dbg CUTCELLS::cutcells IMPORTED_LOCATION_DEBUG)
  get_target_property(_cutcells_lib_any CUTCELLS::cutcells IMPORTED_LOCATION)
  if(_cutcells_lib_rel)
    set(_cutcells_lib "${_cutcells_lib_rel}")
  elseif(_cutcells_lib_dbg)
    set(_cutcells_lib "${_cutcells_lib_dbg}")
  elseif(_cutcells_lib_any)
    set(_cutcells_lib "${_cutcells_lib_any}")
  endif()
endif()

if(_cutcells_lib AND NOT TARGET cutcells)
  message(STATUS "Using CutCells library: ${_cutcells_lib}")
  find_program(NM_EXECUTABLE nm)
  if(NM_EXECUTABLE)
    execute_process(
      COMMAND "${NM_EXECUTABLE}" -gU "${_cutcells_lib}"
      OUTPUT_VARIABLE _nm_out
      ERROR_QUIET
    )
    string(FIND "${_nm_out}" "create_level_set_function" _has_create_level_set_function)
    string(FIND "${_nm_out}" "write_level_set_vtu" _has_write_level_set_vtu)
    if(_has_create_level_set_function EQUAL -1 OR _has_write_level_set_vtu EQUAL -1)
      message(FATAL_ERROR
        "The selected libcutcells is stale and missing required symbols "
        "(create_level_set_function / write_level_set_vtu).\n"
        "Selected library: ${_cutcells_lib}\n"
        "Rebuild and reinstall the C++ core, then rerun the Python build.")
    endif()
  endif()
endif()

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
message(STATUS "Using Python executable: ${Python_EXECUTABLE}")
if(_conda_prefix)
  message(STATUS "Using active conda prefix: ${_conda_prefix}")
endif()

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

set(_build_rpath "$<TARGET_FILE_DIR:CUTCELLS::cutcells>")
if(TARGET cutcells)
  if(APPLE)
    set(_install_rpath "@loader_path/../lib")
  else()
    set(_install_rpath "$ORIGIN/../lib")
  endif()
else()
  set(_install_rpath "$<TARGET_FILE_DIR:CUTCELLS::cutcells>")
endif()

set_target_properties(_cutcellscpp PROPERTIES
  BUILD_RPATH "${_build_rpath}"
  INSTALL_RPATH "${_install_rpath}"
  BUILD_WITH_INSTALL_RPATH FALSE
  INSTALL_RPATH_USE_LINK_PATH FALSE
)

install(TARGETS _cutcellscpp LIBRARY DESTINATION cutcells)
