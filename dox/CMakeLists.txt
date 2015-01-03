# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /dox/CMakeLists.txt
#
# CMakeLists.txt file for dox.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_dox)

# Skip building dox if building app.
if (NOT ("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE") AND
    NOT ("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE_LIBRARY") AND
    NOT ("${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP"))
    return ()
endif ()

find_package (PythonInterp)

set (SEQAN_DOX_PATH "${CMAKE_SOURCE_DIR}/dox")

add_test (NAME build_dox
          COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/util/bin/dox.py
                  -b ${CMAKE_SOURCE_DIR}/core -b ${CMAKE_SOURCE_DIR}/extras
                  -i ${CMAKE_SOURCE_DIR}/core/include/seqan
                  -i ${CMAKE_SOURCE_DIR}/extras/include/seqan
                  -i ${SEQAN_DOX_PATH}/pages
                  --image-dir ${SEQAN_DOX_PATH}/../docs2/images
          WORKING_DIRECTORY ${SEQAN_DOX_PATH})
