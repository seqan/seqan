cmake_minimum_required (VERSION 2.6)
project (core_tests_seeds)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_seeds test_seeds.cpp
				     test_seeds_banded_align.h
				     test_seeds_global_seed_chain.h
				     test_seeds.h
				     test_seeds_memory_manager.h
				     test_seeds_seed_set.h)
