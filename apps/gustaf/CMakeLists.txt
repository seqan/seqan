# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/gustaf/CMakeLists.txt
#
# CMakeLists.txt file for gustaf.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (extras__apps_gustaf)
message (STATUS "Configuring extras/apps/gustaf")

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
add_executable (gustaf gustaf.cpp)

add_executable (gustaf_mate_joining join_mates.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (gustaf ${SEQAN_LIBRARIES})
target_link_libraries (gustaf_mate_joining ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (gustaf)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install gustaf in ${PREFIX}/bin directory
install (TARGETS gustaf gustaf_mate_joining
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/gustaf for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES example/adeno.fa
               example/adeno_modified.fa
               example/adeno_modified_reads.fa
               example/stellar.gff
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (gustaf)


# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:gustaf")
  set (CPACK_PACKAGE_NAME "gustaf")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Gustaf - Generic multi-split alignment finder")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Kathrin Trappe <kathrin.trappe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Kathrin Trappe <kathrin.trappe@fu-berlin.de>")

  seqan_configure_cpack_app (gustaf "gustaf")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:gustaf")

