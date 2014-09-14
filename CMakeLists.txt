cmake_minimum_required (VERSION 2.8)
project (fast_pca)

list(APPEND CMAKE_CXX_FLAGS "-std=c++0x -Wall -pedantic")

find_package(LAPACK REQUIRED)

include_directories(${PROJECT_SOURCE_DIR})
add_subdirectory(fast_pca)
