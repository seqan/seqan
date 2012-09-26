cmake_minimum_required (VERSION 2.6)
project (core_tests_seeds2)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_seeds2 test_seeds2.cpp
					test_align_chain_banded.h
					test_align_dynprog_affine.h
					test_align_dynprog_banded_affine.h
					test_align_dynprog_banded_linear.h
					test_align_dynprog_linear.h
					test_basic_iter_indirect.h
					test_seeds_combination.h
					test_seeds_extension.h
					test_seeds_global_chaining.h
					test_seeds_seed_base.h
					test_seeds_seed_chained.h
					test_seeds_seed_diagonal.h
					test_seeds_seed_set_base.h
					test_seeds_seed_set_unordered.h
					test_seeds_seed_simple.h)
