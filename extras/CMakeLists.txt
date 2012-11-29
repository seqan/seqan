# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/CMakeLists.txt
#
# CMakeLists.txt file for extras.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras)

# Add the paths core/include and extras/include to the paths that CMake
# searches for libraries.
set (CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/include
                                              ${CMAKE_CURRENT_SOURCE_DIR}/../core/include)

# Setup the library modules from extras.
seqan_setup_library (extras core)

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
