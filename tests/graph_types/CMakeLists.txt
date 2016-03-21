# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /tests/graph_types/CMakeLists.txt
#
# CMakeLists.txt file for the graph_types module tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_tests_graph_types CXX)

message (STATUS "Configuring tests/graph_types")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (SeqAn REQUIRED)
endif ()

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

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add targets and register tests.
set (SOURCES property_map derived iterators types utils)
foreach (NAME ${SOURCES})
  add_executable (test_graph_types_${NAME} test_graph_types_${NAME}.cpp)
  target_link_libraries (test_graph_types_${NAME} ${SEQAN_LIBRARIES})
  # Register with CTest
  add_test (NAME test_test_graph_types_${NAME} COMMAND $<TARGET_FILE:test_graph_types_${NAME}>)
endforeach ()

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

