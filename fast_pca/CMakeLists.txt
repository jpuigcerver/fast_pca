add_library(math OBJECT math.h math.cc)
add_library(file OBJECT
  file.h file.cc
  file_ascii.cc
  file_binary.cc
  file_vbosch.cc
  file_octave.h file_octave.cc
  file_htk.h file_htk.cc
  file_mat4.h file_mat4.cc
  )

add_executable(fast_pca
  fast_pca.cc
  $<TARGET_OBJECTS:math>
  $<TARGET_OBJECTS:file>)
target_link_libraries(fast_pca ${LAPACK_LIBRARIES})


add_executable(fast_pca_map
  fast_pca_map.cc
  $<TARGET_OBJECTS:math>
  $<TARGET_OBJECTS:file>)
target_link_libraries(fast_pca_map ${LAPACK_LIBRARIES})


add_executable(fast_pca_reduce
  fast_pca_reduce.cc
  $<TARGET_OBJECTS:math>
  $<TARGET_OBJECTS:file>)
target_link_libraries(fast_pca_reduce ${LAPACK_LIBRARIES})

install_targets(/bin fast_pca fast_pca_map fast_pca_reduce)
