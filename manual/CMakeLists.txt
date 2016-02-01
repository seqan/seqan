# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /manual/CMakeLists.txt
#
# CMakeLists.txt file for manual.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_manual CXX)

# Look for Sphinx and balk out if we cannot find it.

find_package (Sphinx)

if (NOT SPHINX_FOUND)
    message (STATUS "  Sphinx not found, not building sphinx manual as a test.")
    return ()
endif (NOT SPHINX_FOUND)

# Add building Sphinx documentation as a test.

add_test (build_manual
          ${SPHINX_EXECUTABLE} -W -n
          -b html -d ${CMAKE_BINARY_DIR}/manual/doctrees
          ${CMAKE_CURRENT_SOURCE_DIR}/source
          ${CMAKE_BINARY_DIR}/manual/html)

# Building the manual depends on the JSON file from dox.
set_tests_properties(build_manual PROPERTIES DEPENDS build_dox)

