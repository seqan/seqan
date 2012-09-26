cmake_minimum_required (VERSION 2.6)
project (sandbox_singer_tests_fm_index_beta)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_fm_index_beta test_fm_index_beta.cpp test_fm_index_beta.h test_prefix_sum_table_beta.h test_wavelet_tree_structure_beta.h)
