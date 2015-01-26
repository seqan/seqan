# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /dox/CMakeLists.txt
#
# CMakeLists.txt file for dox.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_dox)

find_package (PythonInterp)
if (NOT PYTHONINTERP_FOUND)
    message (STATUS "  Python not found, not building dox as a tests.")
    return ()
endif (NOT PYTHONINTERP_FOUND)

# Add building the documentation as a test.
add_test (build_dox
  ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../util/bin/dox.py
  -b ${CMAKE_CURRENT_SOURCE_DIR}/..
  -i ${CMAKE_CURRENT_SOURCE_DIR}/../include/seqan
  -i ${CMAKE_CURRENT_SOURCE_DIR}/pages
  --image-dir ${CMAKE_CURRENT_SOURCE_DIR}/images)
