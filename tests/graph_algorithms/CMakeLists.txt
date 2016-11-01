# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /tests/graph_algorithms/CMakeLists.txt
#
# CMakeLists.txt file for the graph_algorithms module tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_tests_graph_algorithms CXX)
message (STATUS "Configuring tests/graph_algorithms")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (SeqAn REQUIRED)
endif ()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your test.
add_executable (test_graph_algorithms
               test_graph_algorithms.cpp
               test_graph_algorithms.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_graph_algorithms ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_graph_algorithms COMMAND $<TARGET_FILE:test_graph_algorithms>)
