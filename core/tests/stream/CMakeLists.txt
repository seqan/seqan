# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/tests/stream/CMakeLists.txt
#
# CMakeLists.txt file for the stream module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_core_tests_stream)
message (STATUS "Configuring core/tests/stream")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES ZLIB BZip2)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your test.
add_executable (test_stream
               test_stream.cpp
               test_stream_adapt_cstdio.h
               test_stream_adapt_fstream.h
               test_stream_adapt_sstream.h
               test_stream_bz2_file.h
               test_stream_char_array.h
               test_stream_gz_file.h
               test_stream_file_stream.h
               test_stream_tokenize.h
               test_stream_generic.h
               test_stream_lexical_cast.h
               test_stream_record_reader.h
               test_stream_bgzf.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_stream ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_stream COMMAND $<TARGET_FILE:test_stream>)
