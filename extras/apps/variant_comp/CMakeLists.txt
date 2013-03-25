# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/variant_comp/CMakeLists.txt
#
# CMakeLists.txt file for variant_comp.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (extras_apps_variant_comp)
message (STATUS "Configuring extras/apps/variant_comp")

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
add_executable (variant_comp variant_comp.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (variant_comp ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (variant_comp)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install variant_comp in ${PREFIX}/bin directory
install (TARGETS variant_comp
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/variant_comp for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (variant_comp)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:variant_comp")
  set (CPACK_PACKAGE_NAME "variant_comp")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "variant_comp - Alignment Free Sequence Comparison")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (variant_comp "variant_comp")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:variant_comp")

