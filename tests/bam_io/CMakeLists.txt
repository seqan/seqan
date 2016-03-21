# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /tests/bam_io/CMakeLists.txt
#
# CMakeLists.txt file for the bam_io module tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_tests_bam_io CXX)
message (STATUS "Configuring tests/bam_io")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (ZLIB)
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
add_executable (test_bam_io
               test_bam_io.cpp
               test_bam_alignment_record.h
               test_bam_header_record.h
               test_bam_index.h
               test_bam_io_context.h
               test_bam_sam_conversion.h
               test_bam_tags_dict.h
               test_read_sam.h
               test_write_bam.h
               test_write_sam.h
               test_bam_file.h
               test_bam_index.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_bam_io ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_bam_io COMMAND $<TARGET_FILE:test_bam_io>)
