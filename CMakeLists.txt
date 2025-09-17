cmake_minimum_required(VERSION 3.10)

project(FizzyCompiledSource)

set(CMAKE_CXX_STANDARD 17)
# Find HDF5
set(HDF5_PREFER_PARALLEL true)
find_package(HDF5 REQUIRED)

# Find MPI
find_package(MPI REQUIRED)

if(NOT DEFINED OPENMC_LIB_DIR)
  message(FATAL_ERROR "NO OPENMC_LIB_DIR SET")
endif(NOT DEFINED OPENMC_LIB_DIR)

if(NOT DEFINED OPENMC_INCLUDE_DIR)
  message(FATAL_ERROR "NO OPENMC_INCLUDE_DIR SET")
endif(NOT DEFINED OPENMC_INCLUDE_DIR)

if(NOT DEFINED LIBMESH_INCLUDE)
  message(FATAL_ERROR "NO LIBMESH_INCLUDE SET")
endif(NOT DEFINED LIBMESH_INCLUDE)

add_compile_definitions(OPENMC_LIBMESH_ENABLED)
add_compile_definitions(OPENMC_MPI)

add_library(CompiledSource SHARED ${CMAKE_SOURCE_DIR}/src/FispactCompiledSource.cpp)
target_link_directories(CompiledSource PUBLIC ${OPENMC_LIB_DIR})
target_link_libraries(CompiledSource PUBLIC stdc++ openmc)
target_include_directories(CompiledSource PUBLIC ${CMAKE_SOURCE_DIR}/include ${OPENMC_INCLUDE_DIR} ${HDF5_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_DIRS} ${LIBMESH_INCLUDE} /usr/include/eigen3/ /home/bill/Projects/moose/petsc/arch-moose/include/ /home/bill/Projects/moose/petsc/include/ )
