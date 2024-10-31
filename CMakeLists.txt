cmake_minimum_required(VERSION 3.21.3 FATAL_ERROR)

project(ThreeBody VERSION 1.0)

set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)

add_executable(ThreeBody
               src/main.cpp
               src/BasicFuncs.cpp
               src/BasicFuncs.h
               src/free_2part_spec/FreeSpec.h
               src/free_2part_spec/FreeSpec.cpp)

set(CMAKE_CXX_FLAGS "-O2")
set(CMAKE_C_FLAGS "-O2")
set(CMAKE_EXE_LINKER_FLAGS "-O2")

find_package(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
target_link_libraries(ThreeBody MPI::MPI_CXX)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(ThreeBody Eigen3::Eigen)