set(EXES
    free_2part
    free_coupled22
)

foreach(EXE ${EXES})
    add_executable(${EXE} ${EXE}.cpp)
    target_link_libraries(${EXE} MPI::MPI_CXX)
    target_link_libraries(${EXE} Eigen3::Eigen)
    target_link_libraries(${EXE} 3Body)
endforeach()