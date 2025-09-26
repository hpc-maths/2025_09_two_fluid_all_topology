cmake_minimum_required(VERSION 3.20)

project(relaxation)

find_package(samurai 0.26.0)
find_package(CLI11)
find_package(nlohmann_json REQUIRED)
find_package(MPI)

#include(FindPkgConfig)
#pkg_check_modules(PETSC PETSc)

set(CMAKE_CXX_STANDARD 23)

add_executable(All_topology_test All_topology_solver.cpp)
target_link_libraries(All_topology_test samurai::samurai CLI11::CLI11 nlohmann_json::nlohmann_json)
target_include_directories(All_topology_test PRIVATE ../include)
