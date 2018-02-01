# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/param_chooser/CMakeLists.txt
#
# CMakeLists.txt file for param_chooser.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_param_chooser CXX)
message (STATUS "Configuring apps/param_chooser")

set (SEQAN_APP_VERSION "0.0.9")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (SeqAn CONFIG REQUIRED)
endif ()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Update the list of file names below if you add source files to your application.
add_executable (param_chooser param_chooser.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (param_chooser ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add binary for tests.
add_executable (test_funcs_param_chooser test_param_chooser.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_funcs_param_chooser ${SEQAN_LIBRARIES})

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (param_chooser)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install param_chooser in ${PREFIX}/bin directory
install (TARGETS param_chooser
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/param_chooser for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

# Unit tests for the alignment free tests.
add_test (NAME test_test_funcs_param_chooser COMMAND $<TARGET_FILE:test_funcs_param_chooser>)

#seqan_add_app_test (param_chooser)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:param_chooser")
  set (CPACK_PACKAGE_NAME "param_chooser")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "param_chooser - Parameter Choice for RazerS")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (param_chooser "param_chooser")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:param_chooser")

