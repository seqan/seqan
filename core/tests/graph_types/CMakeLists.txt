# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/tests/graph_types/CMakeLists.txt
#
# CMakeLists.txt file for the graph_types module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_core_tests_graph_types)

message (STATUS "Configuring core/tests/graph_types")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# Clang 3.2.0 has a bug that prevents test_graph_types from compiling.  See
# http://llvm.org/bugs/show_bug.cgi?id=14538
#
# Still present in Clang 3.3.0 but apparently fixed for 3.4.0.
if (COMPILER_IS_CLANG AND ("${_GCC_VERSION}" STREQUAL "320" OR
                           "${_GCC_VERSION}" STREQUAL "330"))
  message (STATUS "  Skipping test_graph_types in Clang 3.2.0/3.3.0.")
  return ()
endif ()


# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your test.
add_executable (test_graph_types
               test_graph_types.cpp
               test_graph_derived.h
               test_graph_iterators.h
               test_graph_properties.h
               test_graph_types.h
               test_graph_utils.h
               test_graph_basic.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_graph_types ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_graph_types COMMAND $<TARGET_FILE:test_graph_types>)
