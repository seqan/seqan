# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/apps/seqcons/CMakeLists.txt
#
# CMakeLists.txt file for SeqCons.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (core_apps_seqcons)
message (STATUS "Configuring core/apps/seqcons")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (seqcons seqcons.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (seqcons ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (seqcons)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install seqcons in ${PREFIX}/bin directory
install (TARGETS seqcons
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/seqcons for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

# SeqCons does not have any app tests right now.
#seqan_add_app_test (seqcons)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:seqcons")
  set (CPACK_PACKAGE_NAME "seqcons")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "SeqCons - Read alignment via realignment or MSA")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "David Weese <david.weese@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (seqcons "seqcons")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:seqcons")

