# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: tests/journaled_string_tree/CMakeLists.txt
#
# CMakeLists.txt file for the journaled_string_tree module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_tests_journaled_string_tree)
message (STATUS "Configuring tests/journaled_string_tree")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

find_package(CXX11 REQUIRED)

if (WIN32)
    message (STATUS "   Not building journaled string tree test under Windows. Incomplete C++11 support.")
    return()
endif (WIN32)

# Enable module only with CXX11 support.
if (NOT CXX11_FOUND)
    message (STATUS "  Cannot build journaled string tree module tests.")
    return ()
endif (NOT CXX11_FOUND)

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES OpenMP)
find_package (SeqAn REQUIRED)

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
                test_journaled_string_tree_traverser.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_journaled_string_tree_traverser ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS}")

# Enable OpenMP.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_journaled_string_tree COMMAND $<TARGET_FILE:test_journaled_string_tree>)
add_test (NAME test_test_delta_map COMMAND $<TARGET_FILE:test_delta_map>)
add_test (NAME test_journaled_string_tree_traverser COMMAND $<TARGET_FILE:test_journaled_string_tree_traverser>)
