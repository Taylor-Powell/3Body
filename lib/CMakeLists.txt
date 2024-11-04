add_library(3Body STATIC)

target_include_directories(3Body PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(3Body MPI::MPI_CXX)
target_link_libraries(3Body Eigen3::Eigen)

set(3BODY_HDRS
    BasicFuncs.h
    FreeSpec.h
    Free2Particle.h
)

install(FILES ${3BODY_HDRS} DESTINATION include/3Body)

target_sources(3Body PUBLIC
    BasicFuncs.cpp
    Free2Particle.cpp
)

install(TARGETS 3Body)