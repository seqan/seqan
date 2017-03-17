# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /dox/CMakeLists.txt
#
# CMakeLists.txt file for dox.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_dox CXX)

find_package (PythonInterp)
if (NOT PYTHONINTERP_FOUND)
    message (STATUS "  Python not found, not building dox as a tests.")
    return ()
endif ()

# find python 2.x
if(NOT(${PYTHON_VERSION_MAJOR} MATCHES "2"))
    set(_PYTHON2_VERSIONS 2.7 2.6 2.5 2.4 2.3 2.2 2.1 2.0)
    unset(_Python_NAMES)
    foreach(_CURRENT_VERSION IN LISTS _PYTHON2_VERSIONS)
        list(APPEND _Python_NAMES "python${_CURRENT_VERSION}")
    endforeach()
    unset(_PYTHON2_VERSIONS)

    find_program(PYTHON2_EXECUTABLE NAMES ${_Python_NAMES})
    if(PYTHON2_EXECUTABLE MATCHES "PYTHON2_EXECUTABLE-NOTFOUND")
        message (STATUS "You need Python 2.x for building dox. (skip the tests)")
        return()
    endif()
    set(PYTHON_EXECUTABLE "${PYTHON2_EXECUTABLE}")
    message (STATUS "Found Python 2.x: ${PYTHON_EXECUTABLE}") 
endif()

# find python 2.x
if(NOT(${PYTHON_VERSION_MAJOR} MATCHES "2"))
    set(_PYTHON2_VERSIONS 2.7 2.6 2.5 2.4 2.3 2.2 2.1 2.0)
    unset(_Python_NAMES)
    foreach(_CURRENT_VERSION IN LISTS _PYTHON2_VERSIONS)
        list(APPEND _Python_NAMES "python${_CURRENT_VERSION}")
    endforeach()
    unset(_PYTHON2_VERSIONS)

    find_program(PYTHON2_EXECUTABLE NAMES ${_Python_NAMES})
    if(PYTHON2_EXECUTABLE MATCHES "PYTHON2_EXECUTABLE-NOTFOUND")
        message (STATUS "You need Python 2.x for building dox. (skip the tests)")
        return()
    endif()
    set(PYTHON_EXECUTABLE "${PYTHON2_EXECUTABLE}")
    message (STATUS "Found Python 2.x: ${PYTHON_EXECUTABLE}")

    # change shell scripts for building dox
    file(READ "${CMAKE_CURRENT_SOURCE_DIR}/../dox/run.sh.in" run_content)
    file(READ "${CMAKE_CURRENT_SOURCE_DIR}/../dox/dox_only.sh.in" dox_only_content)
    string(REPLACE "python2 ../" "${PYTHON_EXECUTABLE} ../" run_content "${run_content}")
    string(REPLACE "python2 ../" "${PYTHON_EXECUTABLE} ../" dox_only_content "${dox_only_content}")
    file(WRITE "${CMAKE_CURRENT_SOURCE_DIR}/../dox/run.sh" "${run_content}")
    file(WRITE "${CMAKE_CURRENT_SOURCE_DIR}/../dox/dox_only.sh" "${dox_only_content}")
endif()

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
             DESTINATION share/doc/seqan)
endif ()
