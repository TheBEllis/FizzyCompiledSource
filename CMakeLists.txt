cmake_minimum_required(VERSION 3.10)

project(FizzyCompiledSource)

set(CMAKE_CXX_STANDARD 17)

if(NOT DEFINED CARDINAL_DIR)
  message(WARNING "NO CARDINAL_DIR SET, MIGHT NOT FIND OPENMC")
endif(NOT DEFINED CARDINAL_DIR)

find_package(PkgConfig REQUIRED)

find_package(OpenMC REQUIRED
             HINTS ${CARDINAL_DIR}/install/cmake/OpenMC)

pkg_check_modules(LIBMESH REQUIRED libmesh-opt)

# Need compile definitions to open up correct code paths
add_compile_definitions(OPENMC_LIBMESH_ENABLED)
add_compile_definitions(OPENMC_MPI)

add_library(CompiledSource SHARED ${CMAKE_SOURCE_DIR}/src/FispactCompiledSource.cpp)
target_link_libraries(CompiledSource PUBLIC stdc++ OpenMC::libopenmc)
target_include_directories(CompiledSource PUBLIC ${CMAKE_SOURCE_DIR}/include ${OPENMC_INCLUDE_DIR}  ${LIBMESH_INCLUDE_DIRS})
