cmake_minimum_required (VERSION 2.6)
project (core_tests_index)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_index_creation
                          test_index_creation.cpp
                          test_index_creation.h)

seqan_add_test_executable(test_index_crosscompare_char
                          test_index_crosscompare_char.cpp
                          test_crosscompare.h)

seqan_add_test_executable(test_index_crosscompare_dna
                          test_index_crosscompare_dna.cpp
                          test_crosscompare.h)

seqan_add_test_executable(test_index_qgram
                          test_index_qgram.cpp
                          test_qgram_index.h)

seqan_add_test_executable(test_index_sa_bwtwalk
                          test_index_sa_bwtwalk.cpp
                          test_sa_bwtwalk.h)

seqan_add_test_executable(test_index_shapes
                          test_index_shapes.cpp
                          test_shapes.h)

seqan_add_test_executable(test_index_drawing
                          test_index_drawing.cpp)

seqan_add_test_executable(test_index_fm
                          test_index_fm.cpp
			  test_index_fm.h
			  test_index_fm_compressed_sa_iterator.h
			  test_index_fm_compressed_sa.h
			  test_index_fm_prefix_sum_table.h
			  test_index_fm_rank_support_bit_string_iterator.h
			  test_index_fm_rank_support_bit_string.h
			  test_index_fm_right_array_binary_tree_iterator.h
			  test_index_fm_right_array_binary_tree.h
			  test_index_fm_sparse_string.h
			  test_index_fm_stree.h
			  test_index_fm_wavelet_tree.h)

# GCC < 4.5 segfaults when compiling the test_index_stree_iterators test.
# Enable if the compiler is not GCC or the version is >= 4.5.
if (NOT CMAKE_COMPILER_IS_GNUCXX OR (450 LESS _GCC_VERSION))
  seqan_add_test_executable(
    test_index_stree_iterators
    test_index_stree_iterators.cpp
    test_stree_iterators.h)
endif ()

if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)

seqan_add_test_executable(test_index_repeats
                          test_index_repeats.cpp
                          test_index_repeats.h)
