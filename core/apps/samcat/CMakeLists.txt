# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/apps/samcat/CMakeLists.txt
#
# CMakeLists.txt file for samcat.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_core_apps_samcat)
message (STATUS "Configuring core/apps/samcat")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
#set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (samcat samcat.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (samcat ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (samcat)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install samcat in ${PREFIX}/bin directory
install (TARGETS samcat
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/samcat for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (samcat)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:samcat")
  set (CPACK_PACKAGE_NAME "samcat")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Concatenate SAM or BAM files")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "David Weese <david.weese@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "David Weese <david.weese@fu-berlin.de>")

  seqan_configure_cpack_app (samcat "samcat")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:samcat")

