get_property(fast_pca_path TARGET fast_pca PROPERTY LOCATION)

add_test(test_gauss2d_ascii "${CMAKE_CURRENT_SOURCE_DIR}/test_ascii.sh" "${fast_pca_path}" )
add_test(test_gauss2d_binary "${CMAKE_CURRENT_SOURCE_DIR}/test_binary.sh" "${fast_pca_path}" )
add_test(test_gauss2d_vbosch "${CMAKE_CURRENT_SOURCE_DIR}/test_vbosch.sh" "${fast_pca_path}" )
add_test(test_gauss2d_octave "${CMAKE_CURRENT_SOURCE_DIR}/test_octave.sh" "${fast_pca_path}" )
add_test(test_gauss2d_htk "${CMAKE_CURRENT_SOURCE_DIR}/test_htk.sh" "${fast_pca_path}" )
add_test(test_gauss2d_mat4 "${CMAKE_CURRENT_SOURCE_DIR}/test_mat4.sh" "${fast_pca_path}" )