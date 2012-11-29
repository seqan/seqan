# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/apps/tree_recon/CMakeLists.txt
#
# CMakeLists.txt file for tree_recon.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (core_apps_tree_recon)
message (STATUS "Configuring core/apps/tree_recon")

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
add_executable (tree_recon tree_recon.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (tree_recon ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (tree_recon)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install tree_recon in ${PREFIX}/bin directory
install (TARGETS tree_recon
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/tree_recon for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES example/example.dist
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (tree_recon)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:tree_recon")
  set (CPACK_PACKAGE_NAME "tree_recon")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "tree_recon - Tree Reconstruction")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "David Weese <david.weese@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (tree_recon "tree_recon")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:tree_recon")

