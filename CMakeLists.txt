cmake_minimum_required (VERSION 2.8)
project (fast_pca)

# Include custom CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

list(APPEND CMAKE_CXX_FLAGS "-std=c++0x -Wall -pedantic")

find_package(LAPACK REQUIRED)
message(STATUS "LAPACK_LINKER_FLAGS: ${LAPACK_LINKER_FLAGS}")
message(STATUS "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")

include_directories(${PROJECT_SOURCE_DIR})
add_subdirectory(fast_pca)

enable_testing()
add_subdirectory(tests)