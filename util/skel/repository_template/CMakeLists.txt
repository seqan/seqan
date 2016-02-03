# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /%(CMAKE_PROJECT_PATH)s/CMakeLists.txt
#
# CMakeLists.txt file for %(CMAKE_PROJECT_PATH)s.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_%(CMAKE_PROJECT_NAME)s CXX)

# Add the paths %(CMAKE_PROJECT_PATH)s/include and %(CMAKE_PROJECT_PATH)s/include to the paths that CMake
# searches for libraries.
set (CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/include
                                              ${CMAKE_CURRENT_SOURCE_DIR}/../../include
                                              ${CMAKE_CURRENT_SOURCE_DIR}/../../include)

# Setup the library modules from core and extras.
seqan_setup_library (%(REPOSITORY_PSEUDO_TARGET_NAME)s core extras)

# Add subdirectory for apps.
add_subdirectory (apps)

# Demos are required when doing a Whole SeqAn Release (copy demos) or when
# developing (build demos).
if (("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE") OR
    ("${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP"))
    add_subdirectory (demos)
endif ()

# Tests are only built when building in DEVLOP mode.
if ("${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    add_subdirectory (tests)
endif ()
