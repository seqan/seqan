cmake_minimum_required (VERSION 2.6)
project (core_tests_random)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_random test_random.cpp
				      test_random_dists.h
				      test_random_rng.h
				      test_random_shuffle.h)
