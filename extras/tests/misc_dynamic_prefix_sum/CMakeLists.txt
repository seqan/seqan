# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/tests/misc_dynamic_prefix_sum/CMakeLists.txt
#
# CMakeLists.txt file for the misc_dynamic_prefix_sum module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras_tests_misc_dynamic_prefix_sum)
message (STATUS "Configuring extras/tests/misc_dynamic_prefix_sum")

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
add_executable (test_misc_dynamic_prefix_sum
                test_misc_dynamic_prefix_sum.cpp
                test_misc_dynamic_prefix_sum.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_misc_dynamic_prefix_sum ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_misc_dynamic_prefix_sum COMMAND $<TARGET_FILE:test_misc_dynamic_prefix_sum>)
