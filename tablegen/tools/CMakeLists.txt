cmake_minimum_required(VERSION 3.21)
project(cutcells_tablegen_tools LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(VTK REQUIRED COMPONENTS
  CommonCore
  CommonDataModel)

add_executable(vtk_iso_refine_dump
  ${CMAKE_CURRENT_SOURCE_DIR}/vtk_iso_refine_dump.cpp)

target_link_libraries(vtk_iso_refine_dump PRIVATE ${VTK_LIBRARIES})

