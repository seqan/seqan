# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/breakpoint_calculator/CMakeLists.txt
#
# CMakeLists.txt file for Breakpoint Calculator.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (extras_apps_breakpoint_calculator)
message (STATUS "Configuring extras/apps/breakpoint_calculator")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# The FindLemon.cmake file is in the current directory.  Thus, add the latter to
# the CMAKE_MODULE_PATH.
set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}")
find_package (LEMON QUIET)

# If the LEMON library could not be found then do not configure Breakpoint
# Calculator.
if (NOT LEMON_FOUND)
  message (STATUS "  LEMON library not found: Not building breakpoint_calculator.")
  return ()
endif (NOT LEMON_FOUND)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS}
                     ${LEMON_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (breakpoint_calculator breakpoint_calculator.cpp
                                     breakpoint_calculator.h
                                     parse_alignment.h
                                     breakpoint_counts.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (breakpoint_calculator ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (breakpoint_calculator)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install breakpoint_calculator in ${PREFIX}/bin directory
install (TARGETS breakpoint_calculator
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/breakpoint_calculator for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES tests/alignment.maf
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (breakpoint_calculator)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:breakpoint_calculator")
    set (CPACK_PACKAGE_NAME "breakpoint_calculator")
    set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "breakpoint_calculator - Alignment Free Sequence Comparison")
    set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
    set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

    seqan_configure_cpack_app (breakpoint_calculator "breakpoint_calculator")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:breakpoint_calculator")

