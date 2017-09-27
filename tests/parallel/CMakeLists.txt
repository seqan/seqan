# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /tests/parallel/CMakeLists.txt
#
# CMakeLists.txt file for the parallel module tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_tests_parallel CXX)
message (STATUS "Configuring tests/parallel")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (OpenMP)
    find_package (SeqAn REQUIRED)
endif ()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Update the list of file names below if you add source files to your test.
add_executable (test_parallel
               test_parallel.cpp
               test_parallel_atomic_misc.h
               test_parallel_atomic_primitives.h
               test_parallel_splitting.h
               test_parallel_queue.h
               test_parallel_thread_pool.h
               test_parallel_enumerable_thread_local.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_parallel ${SEQAN_LIBRARIES})

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_parallel COMMAND $<TARGET_FILE:test_parallel>)
