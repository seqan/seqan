# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/tests/journaled_string_tree/CMakeLists.txt
#
# CMakeLists.txt file for the journaled_string_tree module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras_tests_journaled_string_tree)
message (STATUS "Configuring extras/tests/journaled_string_tree")

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
add_executable (test_journaled_string_tree
                test_journaled_string_tree.cpp
                test_delta_map.h
                test_journaled_string_tree_traverse.h
                test_journaled_string_tree_mock_generator.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_journaled_string_tree ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Enable OpenMP.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_journaled_string_tree COMMAND $<TARGET_FILE:test_journaled_string_tree>)
