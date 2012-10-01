cmake_minimum_required (VERSION 2.6)
project (extras_tests_index)

# Fixing MSVC error C1128: number of sections exceeded object file format limit.
if (MSVC)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif (MSVC)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_index_fm test_index_fm.cpp test_index_fm.h test_index_fm_compressed_sa.h test_index_fm_prefix_sum_table.h test_index_fm_right_array_binary_tree.h test_index_fm_right_array_binary_tree_iterator.h test_index_fm_sparse_string.h test_index_fm_stree.h test_index_fm_wavelet_tree.h test_index_fm_rank_support_bit_string.h test_index_fm_rank_support_bit_string_iterator.h)
