# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/tests/seeds/CMakeLists.txt
#
# CMakeLists.txt file for the seeds module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_core_tests_seeds)
message (STATUS "Configuring core/tests/seeds")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your test.
add_executable (test_seeds
               test_seeds.cpp
               test_basic_iter_indirect.h
               test_seeds_combination.h
               test_seeds_extension.h
               test_seeds_global_chaining.h
               test_seeds_seed_base.h
               test_seeds_seed_chained.h
               test_seeds_seed_diagonal.h
               test_seeds_seed_set_base.h
               test_seeds_seed_set_unordered.h
               test_seeds_seed_simple.h
               test_align_banded_chain_impl.h
               test_banded_chain_alignment_interface.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_seeds ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_seeds COMMAND $<TARGET_FILE:test_seeds>)
