cmake_minimum_required (VERSION 2.6)
project (core_tests_parallel)

# Enable OpenMP if found.
if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_parallel test_parallel.cpp
                          test_parallel_atomic_misc.h
					      test_parallel_atomic_primitives.h
					      test_parallel_splitting.h)
