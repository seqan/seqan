# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: tests/journaled_string_tree/CMakeLists.txt
#
# CMakeLists.txt file for the journaled_string_tree module tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_tests_journaled_string_tree CXX)
message (STATUS "Configuring tests/journaled_string_tree")

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

add_executable (test_basic_observable
                test_basic_observable.cpp
                test_basic_observable.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_basic_observable ${SEQAN_LIBRARIES})

# Update the list of file names below if you add source files to your test.
add_executable (test_journaled_string_tree
                test_journaled_string_tree.cpp
                test_journaled_string_tree.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_journaled_string_tree ${SEQAN_LIBRARIES})

add_executable (test_delta_map
                test_delta_map.cpp
                test_delta_store.h
                test_delta_map.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_delta_map ${SEQAN_LIBRARIES})

add_executable (test_journaled_string_tree_traverser
                test_journaled_string_tree_traverser.cpp
                test_journaled_string_tree_traverser.h
                test_journaled_string_tree_mock.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_journaled_string_tree_traverser ${SEQAN_LIBRARIES})

add_executable (test_journaled_string_tree_find
                test_journaled_string_tree_find.cpp
                test_journaled_string_tree_find.h
                test_journaled_string_tree_mock.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_journaled_string_tree_find ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_journaled_string_tree COMMAND $<TARGET_FILE:test_journaled_string_tree>)
add_test (NAME test_test_delta_map COMMAND $<TARGET_FILE:test_delta_map>)
add_test (NAME test_journaled_string_tree_traverser COMMAND $<TARGET_FILE:test_journaled_string_tree_traverser>)
add_test (NAME test_journaled_string_tree_find COMMAND $<TARGET_FILE:test_journaled_string_tree_find>)
