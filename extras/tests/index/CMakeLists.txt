# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/tests/index/CMakeLists.txt
#
# CMakeLists.txt file for the index module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras_tests_index)
message (STATUS "Configuring extras/tests/index")

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
add_executable (test_index_view
                test_index_view.cpp
                test_index_view.h)
add_executable (test_find_backtracking
                test_find_backtracking.cpp
                test_find_backtracking.h
                test_find_base.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_index_view ${SEQAN_LIBRARIES})
target_link_libraries (test_find_backtracking ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_index_view COMMAND $<TARGET_FILE:test_index_view>)
add_test (NAME test_test_find_backtracking COMMAND $<TARGET_FILE:test_find_backtracking>)
