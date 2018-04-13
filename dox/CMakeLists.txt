# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /dox/CMakeLists.txt
#
# CMakeLists.txt file for dox.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_dox CXX)

# require python 2.7, not python3
set(PythonInterp_FIND_VERSION 2.7)
set(PythonInterp_FIND_VERSION_MAJOR 2)
set(PythonInterp_FIND_VERSION_MINOR 7)
set(PythonInterp_FIND_VERSION_COUNT 2)

find_package (PythonInterp)
if (NOT PYTHONINTERP_FOUND OR NOT(${PYTHON_VERSION_MAJOR} MATCHES "2"))
    message (STATUS "  You need Python 2.x for building dox. (skip the tests)")
    return ()
endif ()

# Add building the documentation as a test.
add_test (build_dox
  ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../util/bin/dox.py
  -b ${CMAKE_CURRENT_SOURCE_DIR}/..
  -i ${CMAKE_CURRENT_SOURCE_DIR}/../include/seqan
  -i ${CMAKE_CURRENT_SOURCE_DIR}/pages
  --image-dir ${CMAKE_CURRENT_SOURCE_DIR}/images)

if (${SEQAN_BUILD_SYSTEM} MATCHES "SEQAN_RELEASE_LIBRARY") # includes SEQAN_RELEASE_LIBRARY
    add_custom_target (dox ALL
        cd ${CMAKE_CURRENT_SOURCE_DIR} && ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../util/bin/dox.py
        -b ${CMAKE_CURRENT_SOURCE_DIR}/..
        -i ${CMAKE_CURRENT_SOURCE_DIR}/../include/seqan
        -i ${CMAKE_CURRENT_SOURCE_DIR}/pages
        --image-dir ${CMAKE_CURRENT_SOURCE_DIR}/images
        --out-dir ${CMAKE_BINARY_DIR}/dox/html)

    install (DIRECTORY ${CMAKE_BINARY_DIR}/dox/html
             DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif ()
