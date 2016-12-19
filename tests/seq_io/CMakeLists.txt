# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /tests/seq_io/CMakeLists.txt
#
# CMakeLists.txt file for the seq_io module tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_tests_seq_io CXX)
message (STATUS "Configuring tests/seq_io")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (ZLIB)
    find_package (BZip2)
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
add_executable (test_seq_io
                test_seq_io.cpp
                test_fai_index.h
                test_genomic_region.h
                test_read_bam.h
                test_write_bam.h
                test_sequence_file.h
                test_stream_read_embl.h
                test_stream_read_genbank.h
                test_stream_write_fasta.h
                test_tag_select_intersect.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_seq_io ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_seq_io COMMAND $<TARGET_FILE:test_seq_io>)
