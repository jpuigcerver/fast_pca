cmake_minimum_required (VERSION 2.8.11)
project (fast_pca)

find_package(BLAS REQUIRED)

add_library(math OBJECT math.h math.cc)
add_library(file OBJECT file.h file.cc)


add_executable(fast_pca_partial
  fast_pca_partial.cc partial.h
  $<TARGET_OBJECTS:math>
  $<TARGET_OBJECTS:file>)
target_link_libraries(fast_pca_partial ${BLAS_LIBRARIES})


add_executable(fast_pca_reduce
  fast_pca_reduce.cc partial.h
  $<TARGET_OBJECTS:math>
  $<TARGET_OBJECTS:file>)
target_link_libraries(fast_pca_reduce ${BLAS_LIBRARIES})