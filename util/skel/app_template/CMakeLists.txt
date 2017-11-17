# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /%(CMAKE_PROJECT_PATH)s/CMakeLists.txt
#
# CMakeLists.txt file for %(NAME)s.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_%(CMAKE_PROJECT_NAME)s CXX)
message (STATUS "Configuring %(CMAKE_PROJECT_PATH)s")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS_MAIN})
include_directories (SYSTEM ${SEQAN_INCLUDE_DIRS_DEPS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (%(NAME)s %(NAME)s.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (%(NAME)s ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (%(NAME)s)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install %(NAME)s in ${PREFIX}/bin directory
install (TARGETS %(NAME)s
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/%(NAME)s for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (%(NAME)s)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:%(NAME)s")
  set (CPACK_PACKAGE_NAME "%(NAME)s")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "%(NAME)s")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "%(AUTHOR)s")
  set (CPACK_PACKAGE_VENDOR "%(AUTHOR)s")

  seqan_configure_cpack_app (%(NAME)s "%(NAME)s")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:%(NAME)s")

