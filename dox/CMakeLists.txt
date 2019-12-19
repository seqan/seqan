# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /dox/CMakeLists.txt
#
# CMakeLists.txt file for dox.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_dox CXX)

find_package (Python3)
if (NOT Python3_Interpreter_FOUND)
    message (STATUS "  You need Python for building dox. (skip the tests)")
    return ()
endif ()

# require jinja2 python package
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import jinja2"
                RESULT_VARIABLE PYTHON_JINJA2_NOT_FOUND
                OUTPUT_VARIABLE _IGNORED
                ERROR_VARIABLE _IGNORED)

if (PYTHON_JINJA2_NOT_FOUND)
  message (STATUS "  The Python package jinja2 is required for building dox ('import jinja2' failed). (skip the tests)")
  return ()
endif (PYTHON_JINJA2_NOT_FOUND)

# Add building the documentation as a test.
add_test (build_dox
  ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../util/bin/dox.py
  -b ${CMAKE_CURRENT_SOURCE_DIR}/..
  -i ${CMAKE_CURRENT_SOURCE_DIR}/../include/seqan
  -i ${CMAKE_CURRENT_SOURCE_DIR}/pages
  --image-dir ${CMAKE_CURRENT_SOURCE_DIR}/images)

if (${SEQAN_BUILD_SYSTEM} MATCHES "SEQAN_RELEASE_LIBRARY") # includes SEQAN_RELEASE_LIBRARY
    add_custom_target (dox ALL
        cd ${CMAKE_CURRENT_SOURCE_DIR} && ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../util/bin/dox.py
        -b ${CMAKE_CURRENT_SOURCE_DIR}/..
        -i ${CMAKE_CURRENT_SOURCE_DIR}/../include/seqan
        -i ${CMAKE_CURRENT_SOURCE_DIR}/pages
        --image-dir ${CMAKE_CURRENT_SOURCE_DIR}/images
        --out-dir ${CMAKE_BINARY_DIR}/dox/html)

    install (DIRECTORY ${CMAKE_BINARY_DIR}/dox/html
             DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif ()
