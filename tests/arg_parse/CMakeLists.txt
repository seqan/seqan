# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /tests/arg_parse/CMakeLists.txt
#
# CMakeLists.txt file for the arg_parse module tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_tests_arg_parse CXX)
message (STATUS "Configuring tests/arg_parse")

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
add_executable (test_arg_parse
                test_arg_parse.cpp
                test_extensions.h
                test_arg_parse.h
                test_arg_parse_argument.h
                test_arg_parse_option.h
                test_argument_parser.h
                test_arg_parse_ctd_support.h)

add_executable (test_arg_parse_version_check
               test_arg_parse_version_check.cpp
               test_arg_parse_version_check.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_arg_parse ${SEQAN_LIBRARIES})
target_link_libraries (test_arg_parse_version_check ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_arg_parse COMMAND $<TARGET_FILE:test_arg_parse>)
add_test (NAME test_test_arg_parse_version_check COMMAND $<TARGET_FILE:test_arg_parse_version_check>)
