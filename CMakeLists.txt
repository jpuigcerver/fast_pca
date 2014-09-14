cmake_minimum_required (VERSION 2.8)
project (fast_pca)

list(APPEND CMAKE_CXX_FLAGS "-std=c++0x -Wall -pedantic")

find_package(LAPACK REQUIRED)
message(STATUS "LAPACK_LINKER_FLAGS: ${LAPACK_LINKER_FLAGS}")
message(STATUS "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")

include_directories(${PROJECT_SOURCE_DIR})
add_subdirectory(fast_pca)
