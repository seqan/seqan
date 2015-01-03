# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /manual/CMakeLists.txt
#
# CMakeLists.txt file for manual.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_manual)

# Look for Sphinx and balk out if we cannot find it.

find_package (Sphinx)

if (NOT SPHINX_FOUND)
    message (STATUS "  Sphinx not found, not building sphinx manual as a test.")
    return ()
endif (NOT SPHINX_FOUND)

# Add building Sphinx documentation as a test.

set (SEQAN_MANUAL_PATH "${CMAKE_SOURCE_DIR}/manual")

add_test (NAME build_manual
          COMMAND ${SPHINX_EXECUTABLE} -W -n
          -b html -d ${CMAKE_BINARY_DIR}/manual/doctrees
          ${CMAKE_CURRENT_SOURCE_DIR}/source
          ${CMAKE_BINARY_DIR}/manual/html
          WORKING_DIRECTORY ${SEQAN_MANUAL_PATH})

      set_tests_properties (build_manual PROPERTIES depends build_dox)
